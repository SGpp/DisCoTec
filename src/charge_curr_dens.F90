#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the charge and current densities from g_1
module charge_curr_dens
  use par_mod
  use discretization
  use discretization_adptv_module
  use vel_space
  use gyro_average_ff_mod, only: jfac, I1_factor, get_jfac
  use gyro_average_df_mod, only: gyro_average_df,&
       &get_nDeriv_base,basex, estimate_gyromatrix_maxbw, xgrid
  use communications
  use geometry
  use aux_func
  use hybrid, only: trap  
  USE MatrixModule
!use adiabatic_response
  implicit none
  public:: initialize_charge_curr_dens, calc_charge_curr_dens, finalize_charge_curr_dens, mem_est_charge_curr_dens
  public:: calculate_density, calculate_current, calculate_density_all, calculate_density_trap
  ! routines for adaptive grids
  public:: calc_charge_curr_dens_adptv
  
  private

  Integer :: init_status = 0

  Real, Dimension(:,:,:,:,:,:,:), Allocatable:: mmat,mmat_perf
  Complex, Dimension(:,:,:,:,:), Allocatable :: vmoments

  ! for hybrid model
  complex, dimension(:,:,:), allocatable :: p_phi_int_trap_pass, p_phi_int2_trap_pass
  Real, Dimension(:,:,:), Allocatable:: fac_inv2_trap_pass


#ifdef G_DAGGER
  character, parameter :: G_op = 'C' !take conjugate transpose 
  !e.g., for first matrix in dot_multiply or in gyro averages
  !if G_op = 'C'
#else
  character, parameter :: G_op = 'N' !take original matrix
#endif


  !TYPE(Matrix), SAVE :: inv_p_doubleprime_ions

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_charge_curr_dens(mem_req_in)
    real:: mem_req_in
    real:: mem_loc

    !vmoments
    mem_loc=SIZE_OF_COMPLEX_MB*lijk0*n_fields*ln0

    if (xy_local) then
       !mmat
       mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lklmn0*n_fields
       !mmat_perf
       if (perf_vec(9).eq.2) then
          mem_loc=mem_loc+mem_est_calc_moments_2()
          !mmat_perf
          mem_loc=mem_loc+SIZE_OF_REAL_MB*lijklmn0*n_fields
       endif
    else

       if (trap_pass) then
          !avg_gamma0_ion
          mem_loc=mem_loc+SIZE_OF_REAL_MB*n_spec*li0*lj0
          !p_phi_int_trap_pass+p_phi_int2
          mem_loc=mem_loc+2*SIZE_OF_COMPLEX_MB*lijk0
          !fac_inv
          mem_loc=mem_loc+SIZE_OF_REAL_MB*lijk0
       endif
       
       !densmu (also denscurr, but not at the same time)
       mem_loc = mem_loc+SIZE_OF_COMPLEX_MB*lijk0*lm0*ln0
       !bardensmu (...)
       mem_loc = mem_loc+SIZE_OF_COMPLEX_MB*lijk0
       
    endif

    mem_est_charge_curr_dens=mem_req_in+mem_loc
  End Function mem_est_charge_curr_dens

  subroutine initialize_charge_curr_dens
    integer:: o,n,i,j,k,l,m
    
    if ((init_status.ne.perf_vec(9)).and.(init_status.gt.0)) &
         & call finalize_charge_curr_dens
    
    if (init_status.eq.0) then
       allocate(vmoments(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2))

       if(xy_local) then
          allocate(mmat(pi1:pi2,pj1:pj2,lk1:lk2,1:n_fields,ll1:ll2,lm1:lm2,ln1:ln2))       
          do o=1,n_fields
             do n=ln1,ln2     
                if (spec(n)%passive) then
                   mmat(:,:,:,o,:,:,n) =0.
                else
                   do k=lk1,lk2
                      do l=ll1,ll2
                         do m=lm1,lm2
                            select case(o)
                            case(1) 
                               if(trap_pass.and.spec(n)%charge.eq.-1) then
                                  mmat(:,:,k,o,l,m,n) = spec(n)%charge*spec(n)%dens*mat_00_trap(:,:,k,l,m)
                               else
                                  mmat(:,:,k,o,l,m,n) = spec(n)%charge*spec(n)%dens*mat_00(:,:,k,l,m)
                               end if
                            case(2)
                               mmat(:,:,k,o,l,m,n) = beta/2*spec(n)%charge*spec(n)%dens*mat_10(:,:,k,l,m,n)
                            CASE(3)
                               mmat(:,:,k,o,l,m,n) = spec(n)%dens * spec(n)%temp * mat_001(:,:,k,l,m)
                            end select
                         enddo
                      enddo
                   enddo
                endif
             enddo
          enddo
          
          if(perf_vec(9).eq.2) then
             allocate(mmat_perf(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,1:n_fields,ln1:ln2))
             do n=ln1,ln2
                do o=1,n_fields
                   do m=lm1,lm2
                      do l=ll1,ll2
                         do k=lk1,lk2
                            do j=lj1,lj2
                               do i=li1,li2
                                  IF (o .EQ. 3) THEN
                                     mmat_perf(i,j,k,l,m,o,n) = mmat(pi1,pj1,k,o,l,m,n) * I1_factor(i,j,k,m,n)
                                  ELSE
                                     mmat_perf(i,j,k,l,m,o,n) = mmat(pi1,pj1,k,o,l,m,n) * jfac(i,j,k,m,n)
                                  END IF
                               enddo
                            enddo
                         enddo
                      enddo
                   enddo
                enddo
             enddo
             deallocate(mmat)
          
             call initialize_calc_moments_2
             
          endif
          
       end if
       
       init_status = perf_vec(9)
    endif
    
  end subroutine initialize_charge_curr_dens

#define NEW_DENS_CALC
  subroutine calc_charge_curr_dens(p_g_1,moments)
    !>modified distribution function of type g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields), intent(out) :: moments
    !local var
    integer, dimension(1:n_fields) :: gy_op
    integer:: n, o

    PERFON('ccdens')
    if (xy_local) then
       if (perf_vec(9).eq.1) then    
          gy_op = 0
          if (n_fields.eq.3) gy_op(n_fields) = 1
          call calc_moments(n_fields,.false.,p_g_1,mmat,vmoments,p_gy_op=gy_op)
       else
          call calc_moments_2(lijk0,llm0,n_fields,p_g_1,mmat_perf,vmoments)
       end if
       if (only_passing) then
         ! in density electron contribution has to be neglected, but it has already been computed, so
         moments=(0.0,0.0)
         do  n=ln1,ln2
             if(spec(n)%charge.ne.-1) then
               moments=moments+vmoments(:,:,:,:,n)
             end if
         end do
       else
         moments=sum(vmoments,5)
       end if
       call my_complex_sum_vwspec(moments,n_fields*lijk0)
    else
       if (trap_pass) then
          do n=ln1,ln2    
             if(spec(n)%charge.eq.-1) then
                call calculate_density_trap(p_g_1,vmoments,n)
             else
                call calculate_density_all(p_g_1,vmoments,n)
             end if
          end do
       else
          call calculate_density(p_g_1, vmoments)
       end if

       if (n_fields.gt.1) call calculate_current(p_g_1, vmoments)

       moments = (0.0,0.0)
       do n=ln1,ln2
          do o=1,n_fields
             if (.not.spec(n)%passive) then
                moments(:,:,:,o) = moments(:,:,:,o)  + spec(n)%charge*spec(n)%dens*vmoments(:,:,:,o,n)
             end if
          end do
       end do       

#ifdef NEW_DENS_CALC
       call my_complex_sum_wspec(moments,n_fields*lijk0)
#else
       call my_complex_sum_vwspec(moments,n_fields*lijk0)
#endif
    end if

   
    PERFOFF

  end subroutine calc_charge_curr_dens


  subroutine finalize_charge_curr_dens
    
    if (xy_local) then
       if(init_status.eq.1) then
          deallocate(mmat)
       else
          deallocate(mmat_perf)
          call finalize_calc_moments_2
       endif       
    end if

    deallocate(vmoments)
    init_status = 0

    !if(allocated(dens_out)) deallocate(dens_out)    
  end subroutine finalize_charge_curr_dens


  !>Calculate gyrocenter density
  !!\param p_g_1 distribution function
  !!\param dens gyrocenter density
  SUBROUTINE calculate_density(p_g_1, vmoms)
    COMPLEX, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(IN) :: p_g_1
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), INTENT(OUT) :: vmoms
    
    COMPLEX, DIMENSION(:,:,:,:,:),allocatable :: densmu
    COMPLEX, DIMENSION(:,:,:),allocatable :: bardensmu
    INTEGER:: j,k,l,m,n,ierror
    logical,parameter :: OUTPUT=.false.
    real :: local_sum,global_sum
    !complex,dimension(li1:li2,lj1:lj2,lk1:lk2) :: recvbuf
    
    ALLOCATE(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
    ALLOCATE(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    IF (OUTPUT) THEN
       CALL calculate_test_sum(p_g_1,local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "calculate_density: p_g_1 = ",global_sum
    END IF

    ! first we do the v_par integration of the distribution function
    densmu = (0.0,0.0)
    DO n=ln1,ln2
       DO m=lm1,lm2
          DO l=ll1,ll2
             Do k=lk1,lk2
                Do j=lj1,lj2
                   densmu(li1:li2,j,k,m,n) = densmu(li1:li2,j,k,m,n) + &
                        &mat_00(:,pj1,k,l,m) * p_g_1(:,j,k,l,m,n)
                END DO
             END DO
          END DO
#ifdef NEW_DENS_CALC
          call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m,n),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
!          call MPI_Allreduce(densmu(li1,lj1,lk1,m,n),recvbuf,lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
!          densmu(:,:,:,m,n)=recvbuf
#endif
       END DO
    END DO

    IF (OUTPUT) THEN
       CALL calculate_test_sum(densmu(:,:,:,:,ln1),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "densmu(ln1) = ",global_sum

       CALL calculate_test_sum(densmu(:,:,:,:,ln2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "densmu(ln2) = ",global_sum
    END IF

    ! now densmu is the same on all vpar processes, so that we can reuse the vpar processes to
    ! speed up the gyro_average routine
    DO n=ln1,ln2
       vmoms(:,:,:,1,n) = (0.0,0.0)
       DO m=lm1,lm2
          CALL gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m,n),bardensmu,m,n,.false.,G_op)
          !write(*,"(3I3,A,ES17.10)") mype,n,m,": bardensmu = ",real(sum(conjg(bardensmu)*bardensmu))
          vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
       END DO
    END DO
    
    DEALLOCATE(densmu,bardensmu)
  END SUBROUTINE calculate_density
  
  
  !> Computes trapped fraction of the density of species n
  !! it is used in the trapped/passing hybrid model for the 
  !! the computation of qn_rhs
  SUBROUTINE calculate_density_trap(p_g_1,vmoms,n)
    COMPLEX, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(IN) :: p_g_1
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), INTENT(OUT) :: vmoms
    
    COMPLEX, DIMENSION(:,:,:,:),allocatable :: densmu
    COMPLEX, DIMENSION(:,:,:),allocatable :: bardensmu
    INTEGER :: j,k,l,m,n, ierror
    
    ALLOCATE(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2))
    ALLOCATE(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    ! first we do the v_par integration of the distribution function
    densmu = CMPLX(0.0,0.0)
    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                densmu(li1:li2,j,k,m) = densmu(li1:li2,j,k,m)+&
                     &mat_00_trap(:,pj1,k,l,m)*p_g_1(:,j,k,l,m,n)
             END DO
          END DO
       END DO
#ifdef NEW_DENS_CALC
       call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
    END DO
    
    vmoms(:,:,:,1,n)= 0.0
    DO m=lm1,lm2
       CALL gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m),bardensmu,m,n,.true.,G_op)
       vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
    END DO
    
    DEALLOCATE(densmu,bardensmu)
    
  END SUBROUTINE calculate_density_trap
  

  !> Compute full (i.e., trapped + passing) density for species n
  !! (should be called for all ions in the trapped/passing hybrid model)
  SUBROUTINE calculate_density_all(p_g_1,vmoms,n)
    COMPLEX, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(IN) :: p_g_1
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), INTENT(OUT) :: vmoms
    
    COMPLEX, DIMENSION(:,:,:,:),allocatable :: densmu
    COMPLEX, DIMENSION(:,:,:),allocatable :: bardensmu
    INTEGER:: j,k,l,m,n,ierror
    
    ALLOCATE(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2))
    ALLOCATE(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    ! first we do the v_par integration of the distribution function
    densmu = CMPLX(0.0,0.0)
    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                densmu(li1:li2,j,k,m) = densmu(li1:li2,j,k,m) + &
                    &mat_00(:,pj1,k,l,m) * p_g_1(:,j,k,l,m,n)
             END DO
          END DO
       END DO
#ifdef NEW_DENS_CALC
       call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
    END DO
    
    vmoms(:,:,:,1,n) = 0.0
    DO m=lm1,lm2
       CALL gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m),bardensmu,m,n,.true.,G_op)
       vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
    END DO
    
    DEALLOCATE(densmu,bardensmu)
   
  END SUBROUTINE calculate_density_all
  

  !> calculate the gyrocenter current along the field lines
  !! \param p_g_1 perturbed distribution function g
  !! \param curr gyrocenter parallel current
  SUBROUTINE calculate_current(p_g_1, vmoms)
    COMPLEX, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(IN) :: p_g_1
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), INTENT(OUT) :: vmoms

    ! Local variables
    COMPLEX, DIMENSION(:,:,:,:,:),allocatable :: currmu
    COMPLEX, DIMENSION(:,:,:),allocatable :: barcurrmu
    INTEGER :: j,k,l,m,n,ierror

    ALLOCATE(currmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
    ALLOCATE(barcurrmu(li1:li2,lj1:lj2,lk1:lk2))

    ! first we do the v_par integration of the distribution function
    !densmu = CMPLX(0.0,0.0)
    currmu = CMPLX(0.0,0.0)
    DO n=ln1,ln2
       DO m=lm1,lm2
          DO l=ll1,ll2
             Do k=lk1,lk2
                Do j=lj1,lj2
                   currmu(li1:li2,j,k,m,n) = currmu(li1:li2,j,k,m,n) + &
                        &mat_10(:,pj1,k,l,m,n) * p_g_1(:,j,k,l,m,n)
                End do
             END DO
          END DO
#ifdef NEW_DENS_CALC
          call MPI_Allreduce(MPI_IN_PLACE,currmu(li1,lj1,lk1,m,n),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
       END DO
    END DO

    DO n=ln1,ln2
       vmoms(:,:,:,2,n)=(0.,0.)
       DO m=lm1,lm2
          CALL gyro_average_df(currmu(li1:li2,lj1:lj2,lk1:lk2,m,n),barcurrmu,m,n,.false.,G_op)
          vmoms(:,:,:,2,n) = vmoms(:,:,:,2,n) + 0.5*beta*barcurrmu
       END DO
    END DO

    DEALLOCATE(currmu, barcurrmu)

  END SUBROUTINE calculate_current

  ! ---- modified routines for adaptive grids ----

  subroutine calc_charge_curr_dens_adptv(p_g_1,moments)
    !>modified distribution function of type g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields), intent(out) :: moments
    !local var
    integer:: n, o

    PERFON('ccdens')
    if (xy_local) then
       stop "no adaptive grids for xy local version!"
    else
       if (trap_pass) then
          do n=ln1,ln2    
             if(spec(n)%charge.eq.-1) then
                call calculate_density_trap_adptv(p_g_1,vmoments,n)
             else
                call calculate_density_all_adptv(p_g_1,vmoments,n)
             end if
          end do
       else
          call calculate_density_adptv(p_g_1, vmoments)
       end if

       if (n_fields.gt.1) call calculate_current_adptv(p_g_1, vmoments)

       moments = (0.0,0.0)
       do n=ln1,ln2
          do o=1,n_fields
             if (.not.spec(n)%passive) then
                moments(:,:,:,o) = moments(:,:,:,o)  + spec(n)%charge*spec(n)%dens*vmoments(:,:,:,o,n)
             end if
          end do
       end do       

#ifdef NEW_DENS_CALC
       call my_complex_sum_wspec(moments,n_fields*lijk0)
#else
       call my_complex_sum_vwspec(moments,n_fields*lijk0)
#endif
    end if
   
    PERFOFF

  end subroutine calc_charge_curr_dens_adptv

  !> Computes trapped fraction of the density of species n
  !! it is used in the trapped/passing hybrid model for the 
  !! the computation of qn_rhs
  subroutine calculate_density_trap_adptv(p_g_1,vmoms,n)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), intent(out) :: vmoms
    
    complex, dimension(:,:,:,:),allocatable :: densmu
    complex, dimension(:,:,:),allocatable :: bardensmu
    integer :: j,k,l,m,n, ierror, i1, i2, p1, p2
    
    allocate(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2))
    allocate(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    ! first we do the v_par integration of the distribution function
    densmu = cmplx(0.0,0.0)
    do m=lm1,lm2
       do l=ll1,ll2
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p1 = pi1_vwadp(l,m)
          p2 = pi2_vwadp(l,m)
          do k=lk1,lk2
             do j=lj1,lj2
                densmu(i1:i2,j,k,m) = densmu(i1:i2,j,k,m)+&
                     &mat_00_trap(p1:p2,pj1,k,l,m)*p_g_1(i1:i2,j,k,l,m,n)
             end do
          end do
       end do
#ifdef NEW_DENS_CALC
       call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
    end do
    
    vmoms(:,:,:,1,n)= 0.0
    ! TODO: return to block-structured grid and gyro-average computation after consultation!
    do m=lm1,lm2
       call gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m),bardensmu,m,n,.true.,g_op)
       vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
    end do
    
    deallocate(densmu,bardensmu)
    
  end subroutine calculate_density_trap_adptv

  !> Compute full (i.e., trapped + passing) density for species n
  !! (should be called for all ions in the trapped/passing hybrid model)
  subroutine calculate_density_all_adptv(p_g_1,vmoms,n)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), intent(out) :: vmoms
    
    complex, dimension(:,:,:,:),allocatable :: densmu
    complex, dimension(:,:,:),allocatable :: bardensmu
    integer:: j,k,l,m,n,ierror, i1, i2, p1, p2
    
    allocate(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2))
    allocate(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    ! first we do the v_par integration of the distribution function
    densmu = cmplx(0.0,0.0)
    do m=lm1,lm2
       do l=ll1,ll2
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p1 = pi1_vwadp(l,m)
          p2 = pi2_vwadp(l,m)
          do k=lk1,lk2
             do j=lj1,lj2
                densmu(i1:i2,j,k,m) = densmu(i1:i2,j,k,m) + &
                    &mat_00(p1:p2,pj1,k,l,m) * p_g_1(i1:i2,j,k,l,m,n)
             end do
          end do
       end do
#ifdef NEW_DENS_CALC
       call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
    end do
    
    vmoms(:,:,:,1,n) = 0.0
    ! TODO: return to block-structured grid and gyro-average computation after consultation!
    do m=lm1,lm2
       call gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m),bardensmu,m,n,.true.,g_op)
       vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
    end do
    
    deallocate(densmu,bardensmu)
   
  end subroutine calculate_density_all_adptv

  !>Calculate gyrocenter density
  !!\param p_g_1 distribution function
  !!\param dens gyrocenter density
  subroutine calculate_density_adptv(p_g_1, vmoms)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), intent(out) :: vmoms
    
    complex, dimension(:,:,:,:,:),allocatable :: densmu
    complex, dimension(:,:,:),allocatable :: bardensmu
    integer:: j,k,l,m,n,ierror, i1, i2, p1, p2
    logical,parameter :: output=.false.
    real :: local_sum,global_sum
      
    allocate(densmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
    allocate(bardensmu(li1:li2,lj1:lj2,lk1:lk2))
    
    if (output) then
       call calculate_test_sum(p_g_1,local_sum, global_sum)
       if (mype.eq.0) write(*,"(a,es20.12)") "calculate_density: p_g_1 = ",global_sum
    end if

    ! first we do the v_par integration of the distribution function
    densmu = (0.0,0.0)
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p1 = pi1_vwadp(l,m)
             p2 = pi2_vwadp(l,m)
             do k=lk1,lk2
                do j=lj1,lj2
                   densmu(i1:i2,j,k,m,n) = densmu(i1:i2,j,k,m,n) + &
                        &mat_00(p1:p2,pj1,k,l,m) * p_g_1(i1:i2,j,k,l,m,n)
                end do
             end do
          end do
#ifdef NEW_DENS_CALC
          call MPI_Allreduce(MPI_IN_PLACE,densmu(li1,lj1,lk1,m,n),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
       end do
    end do

    if (output) then
       call calculate_test_sum(densmu(:,:,:,:,ln1),local_sum, global_sum)
       if (mype.eq.0) write(*,"(a,es20.12)") "densmu(ln1) = ",global_sum

       call calculate_test_sum(densmu(:,:,:,:,ln2),local_sum, global_sum)
       if (mype.eq.0) write(*,"(a,es20.12)") "densmu(ln2) = ",global_sum
    end if

    ! now densmu is the same on all vpar processes, so that we can reuse the vpar processes to
    ! speed up the gyro_average routine
    do n=ln1,ln2
       vmoms(:,:,:,1,n) = (0.0,0.0)
       ! TODO: return to block-structured grid and gyro-average computation after consultation!
       do m=lm1,lm2
          call gyro_average_df(densmu(li1:li2,lj1:lj2,lk1:lk2,m,n),bardensmu,m,n,.false.,g_op)
          vmoms(:,:,:,1,n) = vmoms(:,:,:,1,n) + bardensmu
       end do
    end do
    
    deallocate(densmu,bardensmu)
  end subroutine calculate_density_adptv

  !> calculate the gyrocenter current along the field lines
  !! \param p_g_1 perturbed distribution function g
  !! \param curr gyrocenter parallel current
  subroutine calculate_current_adptv(p_g_1, vmoms)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in) :: p_g_1
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,1:n_fields,ln1:ln2), intent(out) :: vmoms

    ! local variables
    complex, dimension(:,:,:,:,:),allocatable :: currmu
    complex, dimension(:,:,:),allocatable :: barcurrmu
    integer :: j,k,l,m,n,ierror, i1, i2, p1, p2

    allocate(currmu(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
    allocate(barcurrmu(li1:li2,lj1:lj2,lk1:lk2))

    ! first we do the v_par integration of the distribution function
    !densmu = cmplx(0.0,0.0)
    currmu = cmplx(0.0,0.0)
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p1 = pi1_vwadp(l,m)
             p2 = pi2_vwadp(l,m)
             do k=lk1,lk2
                do j=lj1,lj2
                   currmu(i1:i2,j,k,m,n) = currmu(i1:i2,j,k,m,n) + &
                        &mat_10(p1:p2,pj1,k,l,m,n) * p_g_1(i1:i2,j,k,l,m,n)
                end do
             end do
          end do
#ifdef NEW_DENS_CALC
          call MPI_Allreduce(MPI_IN_PLACE,currmu(li1,lj1,lk1,m,n),lijk0,MPI_COMPLEX_TYPE,MPI_SUM,mpi_comm_v,ierror)
#endif
       end do
    end do

    do n=ln1,ln2
       vmoms(:,:,:,2,n)=(0.,0.)
       ! TODO: return to block-structured grid and gyro-average computation after consultation!
       do m=lm1,lm2
          call gyro_average_df(currmu(li1:li2,lj1:lj2,lk1:lk2,m,n),barcurrmu,m,n,.false.,g_op)
          vmoms(:,:,:,2,n) = vmoms(:,:,:,2,n) + 0.5*beta*barcurrmu
       end do
    end do

    deallocate(currmu, barcurrmu)

  end subroutine calculate_current_adptv
  
end module charge_curr_dens
