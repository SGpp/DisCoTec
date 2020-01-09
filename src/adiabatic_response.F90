#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"

module adiabatic_response
  use discretization
  use par_in
  use par_other
  use coordinates
  use communications
  use geometry
  use aux_func
  use gyro_average_ff_mod
  use MatrixModule
  use spatial_averages
  use hybrid, only: trap

  implicit none
  public:: mem_est_adiabatic_response, check_for_adiabatic_species, adiabatic_electrons
  public:: initialize_adiabatic_response, finalize_adiabatic_response
  public:: dens_electrons, temp_electrons, inv_p_doubleprime_ions_m, inv_p_doubleprime_ions_v
  public:: add_adiabatic
  public:: add_adiabatic_passing_electrons
  private

  real, dimension(:), allocatable:: dens_electrons, temp_electrons
  logical:: adiabatic_electrons=.true., adiabatic_ions=.true.

  type(matrix), save:: rhsvec, solvec
  type(matrix), save:: inv_p_doubleprime_ions_m
  complex, dimension(:), allocatable:: inv_p_doubleprime_ions_v
  real, dimension(:), allocatable:: myden 
  real, dimension(:,:,:), allocatable:: sum_1
  real, dimension(:,:,:), allocatable:: tmp_sum1
  
  contains
  
  !>Give an estimate of the memory requirements of this module
  real function mem_est_adiabatic_response(mem_req_in)
    real:: mem_req_in
    real:: mem_loc
   
    mem_loc=0.
    if(.not.x_local) then
       ! inv_p_doubleprime_ions_m
       mem_loc = mem_loc + static_size_of_matrix()/1024./1024.&
            & + SIZE_OF_COMPLEX_MB*ni0*li0/n_procs_x
       ! rhsvec, solvec
       mem_loc = mem_loc + 2.*static_size_of_matrix()/1024./1024.    
    end if

    mem_loc=mem_loc+SIZE_OF_REAL_MB*n_spec*li0*lj0   

    mem_est_adiabatic_response=mem_req_in + mem_loc

  end function mem_est_adiabatic_response

  !>Checks quasineutrality of the equilibrium profiles both for
  !!local and global simulations.
  subroutine check_for_adiabatic_species
    real,dimension(px0):: qn, qnomn
    logical:: write_pe
    integer:: n

    write_pe = ((mype.le.0).AND.(print_ini_msg))

    !qn = sum_j q_j n_j
    !qnomn = sum_j q_j n_j omn_j = sum_j q_j dr n_j
    qn = 0.0; qnomn = 0.0
    do n=0,n_spec-1
       if (.not.spec(n)%passive) then
          if (spec(n)%charge.eq.-1) adiabatic_electrons=.false.
          if (spec(n)%charge.gt.0) adiabatic_ions=.false.
          qn=qn+ spec(n)%charge*spec(n)%dens*spec(n)%dens_prof
          if (x_local) then
             qnomn = qnomn + spec(n)%charge*&
                  &spec(n)%dens*spec(n)%omn
          else
             qnomn = qnomn+spec(n)%charge*&
                  &spec(n)%dens*spec(n)%dens_prof*&
                  &spec(n)%omn_prof
          endif
       endif
    enddo

    if (adiabatic_electrons.or.no_electron_response) then
       qn = 0.0 !qn - spec(0)%dens*spec(0)%dens_prof !\TODO multiple ion species -> qn = 0.0?
       qnomn = 0.0 !\TODO add qnomn check for ad. electrons
    endif
    if (adiabatic_ions) then
       qn = qn + spec(0)%dens*spec(0)%dens_prof
       qnomn = 0.0 !\TODO add qnomn check for ad. ions
    endif
       
    IF(adiabatic_electrons.and.trap_pass) THEN
       STOP "Hybrid model should be run with electron species" 
    END IF

    if (any(abs(qn).gt.1e-6)) then
       print *, 'Maximum violation of quasineutrality: ', maxval(abs(qn))
       stop 'Violation of quasineutrality, please change charge and/or dens'       
    endif
    if (any(abs(qnomn).gt.1e-6)) then
       print *, 'Maximum violation of quasineutrality in gradients: ', &
            &maxval(abs(qnomn))
       stop 'Violation of quasineutrality in density gradients'
    endif    
    
    if (adiabatic_electrons) then
       if (write_pe) write(*,"(A)") "Electrons are adiabatic, beta is set to zero"
       if (approx_avg) then
          if (write_pe) write(*,"(A)") "Approximating <AB>=<A><B> for  <phi>"
       endif
       beta=0.0
    elseif (adiabatic_ions) then
       if (write_pe) write(*,"(A)") "Ions are adiabatic"
    elseif (trap_pass) then
       if (write_pe) write(*,"(A)") "Trapped kinetic & Passing adiabatic, beta is set to zero"
       beta=0.0
       if (approx_avg) then
         if (write_pe) write(*,"(A)") "Approximating <AB>=<A><B> for  <phi>"
       endif
       if (only_passing) then
         if (write_pe) write(*,"(A)") "Neglecting trapped electron contribution"
       endif
    elseif (no_electron_response) then
       if (write_pe) write(*,"(A)") "No-response electrons, beta is set to zero"
       beta=0.0
    else
       if (write_pe) write(*,"(A)") "Kinetic electrons and ions"
    endif

  end subroutine check_for_adiabatic_species


  subroutine initialize_adiabatic_response
    integer:: n

    allocate(temp_electrons(pi1gl:pi2gl), dens_electrons(pi1gl:pi2gl))
    if (adiabatic_electrons) then
       !for the moment the electron temp. is equal to the first ion species
       temp_electrons = spec(0)%temp*spec(0)%temp_prof
       dens_electrons = 0.0
       do n=ln1,ln2
          if (.not.spec(n)%passive) then
             dens_electrons = dens_electrons+&
                  &spec(n)%charge*spec(n)%dens*spec(n)%dens_prof
          end if
       end do
    else
       do n=0,n_spec-1
          if (spec(n)%charge.eq.-1) then 
             temp_electrons = spec(n)%temp*spec(n)%temp_prof
             dens_electrons = spec(n)%dens*spec(n)%dens_prof
          end if
       enddo
    endif

    if (.not.x_local) then
       call initialize(inv_p_doubleprime_ions_m,ni0,ni0)
       call allocate(inv_p_doubleprime_ions_m)
       call set_zero(inv_p_doubleprime_ions_m)
       call initialize(rhsvec,ni0)
       call initialize(solvec,ni0)
       allocate(inv_p_doubleprime_ions_v(lg1:lg2))
    else
       if (.not.y_local) then
          allocate(inv_p_doubleprime_ions_v(lg1:lg2))
          inv_p_doubleprime_ions_v=0.
       else
          allocate(myden(lg1:lg2))
          allocate(sum_1(li1:li2,lj1:lj2,lk1:lk2))
          myden=0.
          !       allocate(elec_cont(li1:li2,lj1:lj2,lk1:lk2))
          if (xy_local) call initialize_add_adiabatic_ff
       endif
    end if

  end subroutine initialize_adiabatic_response


  !>Adds the contribution of the adiabatic species (if any) to the charge density 
  subroutine add_adiabatic(dens)
    !>charge density to be modified
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(inout) :: dens

    ! Local variables
    complex, dimension(lg1:lg2):: response
    integer:: j,k

    PERFON('addadi')


    call compute_adiabatic(dens,response)
        
    if (yx_order) then
       if (y_local) then
          if (p_has_0_mode) then 
             do k=lk1,lk2
                do j=lj1,lj2
                   dens(li1,j,k) = dens(li1,j,k) + response(j)
                end do
             end do
          end if
       else
          do k=lk1,lk2
             do j=lj1,lj2
                dens(:,j,k) = dens(:,j,k) + response(j)
             end do
          end do
       end if
    else
       if (p_has_0_mode) then 
          if (x_local) then
             do k=lk1,lk2
                dens(:,lj1,k) = dens(:,lj1,k) + response    
             end do
          else
             do k=lk1,lk2
                dens(:,lj1,k) = dens(:,lj1,k) + response * &
                     dens_electrons(li1:li2)/temp_electrons(li1:li2)
             end do
          end if
       end if
    end if
    PERFOFF

  end subroutine add_adiabatic


  subroutine add_adiabatic_passing_electrons(dens)
    !>charge density to be modified
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(inout) :: dens

    ! Local variables
    complex, dimension(lg1:lg2):: response
    integer:: j, k, n
    
    call compute_adiabatic(dens,response)
 
    do n = 0, n_spec-1
       if(spec(n)%charge.eq.-1) then
          if (yx_order) then            
             if (y_local) then
                if (p_has_0_mode) then
                   do k = lk1, lk2
                      do j = lj1, lj2
                         dens(li1, j, k) = dens(li1, j, k) + (1.0-trap(pi1,k)) * &
                              &dens_electrons(pi1)/temp_electrons(pi1) * response(j)
                      end do
                   end do
                end if
             else
                do k = lk1, lk2
                   do j = lj1, lj2
                      dens(:, j, k) = dens(:, j, k) + (1.0-trap(pi1,k)) * &
                           &dens_electrons(pi1)/temp_electrons(pi1) * response(j)
                   end do
                end do
             end if
          else
             if (p_has_0_mode) then
                if (x_local) then
                   do k = lk1, lk2
                      dens(:, lj1, k) = dens(:, lj1, k) + (1.0-trap(pi1,k)) * &
                           &dens_electrons(pi1)/temp_electrons(pi1)* response
                   end do
                else
                   do k = lk1, lk2
                      dens(:, lj1, k) = dens(:, lj1, k) + (1.0-trap(li1:li2,k)) * &
                           &dens_electrons(li1:li2)/temp_electrons(li1:li2)* response
                   end do
                end if
             end if
          end if
       end if
    end do
    
  end subroutine add_adiabatic_passing_electrons


  subroutine compute_adiabatic(dens,response) 
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(in):: dens
    complex, dimension(lg1:lg2),intent(OUT):: response
    ! Local variables
    complex, dimension(lg1:lg2):: dens_avg

    PERFON_I('compadi')
    
    if (.not.x_local) then
       call flux_surface_average(dens,.false.,dens_avg)
       call solve_for_fs_avg(inv_p_doubleprime_ions_m, dens_avg, response)
       ! in the adiabatic electrons case, the temperatures are normalized to the 
       ! electron temperature. dens_electrons is calculated in the initializiation to
       ! be such that quasineutrality is fulfilled for the equilibrium.
    else
       if (.not.y_local) then
          call flux_surface_average(dens,.false.,dens_avg)
          response=dens_avg*inv_p_doubleprime_ions_v
       else
          !solve in local case
          if (approx_avg) then
             call flux_surface_average(dens, .false., dens_avg)
          else
             call flux_surface_average(dens/sum_1(:,:,lk1:lk2), .false., dens_avg)
          end if
          response=dens_avg*myden
       endif
    end if

    PERFOFF_I

  end subroutine compute_adiabatic

  subroutine finalize_adiabatic_response
    deallocate(dens_electrons, temp_electrons)
    
    if (.not.x_local) then
       call finalize(inv_p_doubleprime_ions_m)
       call finalize(rhsvec)
       call finalize(solvec)
       deallocate(inv_p_doubleprime_ions_v)
    else
       if (.not.y_local) then 
          deallocate(inv_p_doubleprime_ions_v)
       else
          deallocate(myden,sum_1)
       endif
    endif
   
  end subroutine finalize_adiabatic_response



  !some auxiliary routines------------------------------------------------

  !>Solves for flux surface average (required for adiabatic electrons)
  !!\param mat inverse left hand side operator
  !!\param rhs the right hand side of the field equations
  !!\param res the resulting potential/field
  subroutine solve_for_fs_avg(mat, rhs, res)
    type(matrix) :: mat
    complex, dimension(li1:li2), intent(inout) :: rhs
    complex, dimension(li1:li2), intent(inout) :: res

    integer :: ierr !, irow,icol
    complex :: mean_rhs_loc, mean_rhs

    if ((rad_bc_type.eq.0).and.(p_has_0_mode)) then
       ! delete the mean of the rhs
       mean_rhs = sum(rhs)
       call mpi_allreduce(mean_rhs_loc, mean_rhs, 1, MPI_COMPLEX_TYPE,&
            &mpi_sum, mpi_comm_x, ierr)
       rhs = rhs - mean_rhs/ni0
       if (li1 == 0) rhs(li1) = (0.0,0.0)
    end if
    call attach(rhsvec,rhs)
    call attach(solvec,res)
    call dot_multiply(mat,rhsvec,solvec)

  end subroutine solve_for_fs_avg



  !initialization for xy_local
  subroutine initialize_add_adiabatic_ff
    real:: b2spec
    real, dimension(lg1:lg2,0:n_spec-1) :: avg_gamma0
    real, dimension(li1:li2,lj1:lj2,lk1:lk2):: gamma
    real, dimension(li1:li2,lj1:lj2,lk1:lk2,0:n_spec-1):: my_gamma
    real, dimension(li1:li2,lj1:lj2,lk1:lk2)::elec_cont
    real :: k_perp2
    integer:: i,j,k,m,n
    real :: zint_mu, zarg
#ifdef EXTERNAL_ERF
  EXTERNAL erf
#endif
    real :: erf

    ! Calculate the flux surface average of Gamma_0(b_i) which is
    ! necessary for calculation of flux surface average of phi

    avg_gamma0 = 0.0d0
    my_gamma=0.0d0
    do n=0,n_spec-1
       if (spec(n)%charge .ne. -1) then
          do k=lk1,lk2
             do j=lj1,lj2
                do i=li1,li2
                   k_perp2 = geom%gii(pi1,pj1,k)*ki(i)**2+&
                        &2*geom%gij(pi1,pj1,k)*ki(i)*kj(j)+&
                        &geom%gjj(pi1,pj1,k)*kj(j)**2
                   b2spec = k_perp2*spec(n)%mass*spec(n)%temp/&
                        &(spec(n)%charge**2*geom%Bfield(pi1,pj1,k)**2)
                   gamma(i,j,k) = gamma0(b2spec)
                   enddo
             enddo
          end do
          my_gamma(:,:,:,n)=gamma 
       end if
    end do

    if (trap_pass) then
       ! Now we have to compute the flux surface average of the electrons contribution
       ! given in subroutine initialize_field_solve_ff of field_solve_ff.F90 
       ! term 1 averaged. Here all the term built in initialize_fieldsolve_ff 
       ! should be built again. 
       elec_cont=0.0d0
       do n=ln1,ln2
          IF (spec(n)%charge .eq. -1) THEN
             do i=li1,li2
                do j=lj1,lj2
                   do k=lk1,lk2
                      zint_mu = 0.0
                      do m = lm1, lm2
                         zarg = sqrt(mu(m)*geom%Bfield(pi1,pj1,k)) * &
                              trap(pi1,k) / sqrt(1.0 - trap(pi1,k)**2)
                         zint_mu = zint_mu + mu_weight(m) * &
                              get_jfac(i,j,k,m,n)**2 * erf(zarg) * &
                              EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)))
                      end do
                      
                      ! ... now the global sum over mu-domains
                      CALL my_real_sum_w_0d(zint_mu)
                      
                      zint_mu = geom%Bfield(pi1,pj1,k) * zint_mu

                      elec_cont(i,j,k)= spec(n)%charge**2 * spec(n)%dens / spec(n)%temp *zint_mu
                   end do
                end do
             end do
          ENDIF
       end do
       call my_sum_to_all_real(elec_cont, SIZE(ELEC_CONT), MPI_COMM_SPEC)
    end if
    ! build sum1 is used only in the local code; its value depends on the 
    ! model as well as the way the fsavg is computed, making the initialization
    ! a bit messy, only_passing means hybrid but neglecting kinetic contribution

    if (trap_pass) then
      if (approx_avg) then      
        if (only_passing) then
          sum_1=0.0d0
        else
          do k=lk1,lk2
            sum_1(:,:,k)=trap(pi1,k)*dens_electrons(pi1)/temp_electrons(pi1)-elec_cont(:,:,k)
          end do
        end if
      else
        if (only_passing) then
          do k=lk1,lk2
            sum_1(:,:,k)=(1-trap(pi1,k))*dens_electrons(pi1)/temp_electrons(pi1)
          end do
        else
          do k=lk1,lk2
            sum_1(:,:,k)=dens_electrons(pi1)/temp_electrons(pi1)-elec_cont(:,:,k)
          end do
        end if
       end if
    else
       !full adiabatic case, in this case electron is reference
       if (approx_avg) then
         sum_1=0.0d0
       else
         sum_1=1.0d0
       end if
    end if 

    do n=0,n_spec-1
       if (.not.spec(n)%passive) then
          if (spec(n)%charge .ne. -1) then
             sum_1=sum_1 + spec(n)%charge**2 * spec(n)%dens/&
                  &spec(n)%temp * (1.0d0 - my_gamma(:,:,:,n))
          end if
        end if
    end do
 
    !build the multiplying factor of gyrodensity for fsavg phi
    !this has to be done accordingly to way of average and model

    if (approx_avg) then
       !this builds the fsavg as done in past, i.e. <AB>=<A><B>
       if (trap_pass) then
         call flux_surface_average(sum_1, .false., myden)
       else
         call flux_surface_average(sum_1, .false., myden)
       endif
       myden=1.0d0/myden
     else
      ! exact evaluation of fsavg phi
      if (trap_pass) then
         allocate(tmp_sum1(li1:li2,lj1:lj2,lk1:lk2))
         tmp_sum1=0 
         do k=lk1,lk2 
           tmp_sum1(:,:,k)=dens_electrons(pi1)/temp_electrons(pi1)*(1-trap(pi1,k))/sum_1(:,:,k)
         end do
         call flux_surface_average(tmp_sum1, .false., myden)
         deallocate(tmp_sum1)
      else
         call flux_surface_average(1/sum_1, .false., myden)
      end if  
      if (p_has_00_mode) then
         myden(lg1)=0.
         myden(lg1+1:)=1.0d0/(1.0d0-myden(lg1+1:))
      else
         myden=1.0d0/(1.0d0-myden)
      endif
    end if
   

  end subroutine initialize_add_adiabatic_ff


end module adiabatic_response
