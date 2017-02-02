#include "redef.h"
MODULE boundary_exchange_z
  use BoundaryDescriptionModule
  use boundary_exchange_general
  use box_data_module
  use discretization
  use par_in, only: shifted_metric
  use par_other, only: pi, imag, print_ini_msg
  use communications
  use fourier

  implicit none

  PRIVATE
  PUBLIC :: exchange_z_df, exchange_z_nopb, initialize_boundary_exchange_z, finalize_boundary_exchange_z
  PUBLIC :: exchange_z_df_noinner
  PUBLIC :: exchange_z_3d_equil
  
  INTERFACE exchange_z_df
     MODULE PROCEDURE exchange_z_3D, exchange_z_4D,exchange_z_5D, exchange_z_6D
  END INTERFACE

  INTERFACE exchange_z_nopb
     MODULE PROCEDURE exchange_z_3D_nopb
  END INTERFACE
  
  INTERFACE exchange_z_df_noinner
     MODULE PROCEDURE exchange_z_3D_noinner
  END INTERFACE

  INTERFACE initialize_boundary_exchange_z
     module procedure bez_initialize_boundary_exchange_z
  END INTERFACE

  !INTEGER :: ubexc

  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: pb_phase_fac

CONTAINS
  
  !\todo add qprof to simulation_box type?
  SUBROUTINE bez_initialize_boundary_exchange_z(ai_comm_z, ai_simulation_box,&
       &q_prof)
    integer :: ai_comm_z
    type(box_data_type) :: ai_simulation_box
    real, dimension(pi1gl:pi2gl) :: q_prof

    ! Local variables
    INTEGER:: i, j
    REAL:: q_value, xposition
    Complex, Dimension(:,:), Allocatable :: temp_fac, inv_temp_fac
    integer :: n0_global
    REAL :: kymin,lx,shat,q0,rhostar
    logical :: lilo


    ! global simulation
    lilo      = get_lilo(ai_simulation_box)
    n0_global = get_n0_global(ai_simulation_box)
    q0        = get_q0(ai_simulation_box)
    rhostar   = get_rhostar(ai_simulation_box)

    IF (.NOT.HasMagneticProfile(ai_simulation_box)) THEN
       kymin     = get_kymin(ai_simulation_box)
       lx        = get_lx(ai_simulation_box)
       shat      = get_shat(ai_simulation_box)
    END IF

    allocate(temp_fac(li1:li2,lj1:lj2),inv_temp_fac(li1:li2,lj1:lj2))
    allocate(pb_phase_fac(li1:li2,lj1:lj2,2))
    pb_phase_fac = CMPLX(0.,0.)

    Do j=lj1,lj2
       Do i=li1,li2
          if (HasMagneticProfile(ai_simulation_box)) then
             q_value = q_prof(i)
          else
             IF (lilo) THEN !use "local: parallel boundary condition
                xposition = get_xposition(ai_simulation_box,i-pi1gl+1)
                q_value = kymin*shat*xposition
                !for shifted metric, we need periodic boundaries; 
                !setting q_value=1 eliminates the shifts (n0*q0 is kept!)
                if (shifted_metric) q_value=1.
                if (n0_global.ne.0.and.n0_global.ne.-1111) q_value = q_value + n0_global*q0
                temp_fac(i,j) = exp(imag*(j+ky0_ind)*2.0*pi*q_value)
                inv_temp_fac(i,j) = exp(-imag*(j+ky0_ind)*2.0*pi*q_value)
             else
                ! local code pbc
                !xposition = (-0.5 + (i-pi1gl)/REAL(nx0-1))*lxarr(j)
                xposition = get_xposition(ai_simulation_box,i-pi1gl+1)
                !xposition = (-0.5 + (i-pi1gl)/REAL(nx0-1))*lx
                q_value = q0*(1.0 + shat*xposition*rhostar/ai_simulation_box%x0)
             END IF
          end if
          IF (.NOT.lilo) THEN
             temp_fac(i,j) = EXP(imag*(j+ky0_ind)*2.0D0*pi*n0_global*q_value)
             !for shifted metric, we only need n0*q0 here, shear is taken account
             !in the calculation of the phases applied in the parallel derivatives
             if (shifted_metric) temp_fac(i,j)=exp(imag*(j+ky0_ind)*2.0D0*pi*n0_global*q0)
             inv_temp_fac(i,j) = CONJG(temp_fac(i,j)) 
                        !EXP(-imag*j*2.0D0*pi*n0_global*q_value)
          END IF
       End Do
    end Do

    pb_phase_fac(li1:li2,lj1:lj2,2) = temp_fac(li1:li2,lj1:lj2)
    pb_phase_fac(li1:li2,lj1:lj2,1) = inv_temp_fac(li1:li2,lj1:lj2)

    deallocate(temp_fac, inv_temp_fac)

  END SUBROUTINE bez_initialize_boundary_exchange_z

  SUBROUTINE finalize_boundary_exchange_z
    IF (ALLOCATED(pb_phase_fac)) DEALLOCATE(pb_phase_fac)
  END SUBROUTINE finalize_boundary_exchange_z

  SUBROUTINE exchange_z_3D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u))&
         &.OR.bdesc%count.ne.1) THEN
       STOP "Wrong BoundaryDescription in call to exchange_z_3D!"
    END IF

    CALL exchange_z_general(bdesc, u,1)

  END SUBROUTINE exchange_z_3D

  SUBROUTINE exchange_z_3D_nopb(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u))&
         &.OR.bdesc%count.ne.1) THEN
       STOP "Wrong BoundaryDescription in call to exchange_z_3D!"
    END IF

    CALL exchange_z_general_nopb(bdesc,u,1)

  END SUBROUTINE exchange_z_3D_nopb


  SUBROUTINE exchange_z_4D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3))&
         &.OR.bdesc%count.gt.SIZE(u,4)) THEN
       PRINT*,"Wrong BoundaryDescription in call to exchange_z_4D!"
       PRINT*,"bdesc%n_points = ",bdesc%n_points,", bdesc%count = ", bdesc%count
       STOP
    END IF
    
    CALL exchange_z_general(bdesc, u,SIZE(u,4))
  END SUBROUTINE exchange_z_4D

  SUBROUTINE exchange_z_5D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3))&
         &.OR.bdesc%count.GT.SIZE(u,4)*SIZE(u,5)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_z_5D!"
    END IF
    
    CALL exchange_z_general(bdesc, u,SIZE(u,4)*SIZE(u,5))
  END SUBROUTINE exchange_z_5D

  SUBROUTINE exchange_z_6D(bdesc,u)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(:,:,:,:,:,:) :: u

    IF ((bdesc%n_points.NE.SIZE(u,1)*SIZE(u,2)*SIZE(u,3))&
         &.OR.bdesc%count.GT.SIZE(u,4)*SIZE(u,5)*SIZE(u,6)) THEN
       STOP "Wrong BoundaryDescription in call to exchange_z_6D!"
    END IF
    
    CALL exchange_z_general(bdesc, u,SIZE(u,4)*SIZE(u,5)*SIZE(u,6))
  END SUBROUTINE exchange_z_6D

  !> General z exchange routine without application of boundary conditions (e.g. for exchange of metric)
  SUBROUTINE exchange_z_general_nopb(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    COMPLEX, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange

    PERFON('exz_nopb')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_z_general! We cannot exchange ",&
            & n_dim_2," times in z direction in blocks of ",bdesc%count,". Aborting!"
       PRINT*, "bdesc = ",bdesc
       STOP
    END IF

    bdesc%exchange_direction = 3 ! set to z direction
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       CALL exchange_general(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO
    PERFOFF
  END SUBROUTINE exchange_z_general_nopb


  !this routine exchanges z boundary points for the case of
  !each z process having a full aray (only external boundary, not x_local)
  !as it is necessary for read-in of checkpoints with changing z resolution
  !note that the mpi communicator in z direction has to be ignored,
  !which makes this task incompatible with the exchange_general subroutine
  SUBROUTINE exchange_z_3d_noinner(u)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz) :: u
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,0:nzb-1) :: tmp
    integer :: k

    !lower bc
    tmp = u(li1:li2,lj1:lj2,lk2-nzb+1:lk2)
    do k=0,nzb-1
      tmp(:,:,k)=tmp(:,:,k)*pb_phase_fac(li1:li2,lj1:lj2,2)
    enddo
    u(li1:li2,lj1:lj2,lbz:lk1-1)=tmp(:,:,:)
    !upper bc
    tmp = u(li1:li2,lj1:lj2,lk1:lk1+nzb-1)
    do k=0,nzb-1
      tmp(:,:,k)=tmp(:,:,k)*pb_phase_fac(li1:li2,lj1:lj2,1)
    enddo
    u(li1:li2,lj1:lj2,lk2+1:ubz)=tmp(:,:,:)

  END SUBROUTINE exchange_z_3d_noinner

  ! =======================================================

  SUBROUTINE exchange_z_general(bdesc,u, n_dim_2)
    TYPE(BoundaryDescription) :: bdesc
    INTEGER :: n_dim_2
    COMPLEX, DIMENSION(bdesc%n_points,1:n_dim_2) :: u

    ! Local variables
    INTEGER :: n_exchanges, i_exchange

    PERFON('exz')
    IF (MOD(n_dim_2,bdesc%count).NE.0) THEN
       WRITE(*,"(2(A,I6),A)") "Wrong BoundaryDescription in call to exchange_z_general! We cannot exchange ",&
            & n_dim_2," times in z direction in blocks of ",bdesc%count,". Aborting!"
       PRINT*, "bdesc = ",bdesc
       STOP
    END IF

    ! Here we apply the physical boundary conditions as pre-process of the exchange
    ! End of application of physical boundary conditions

    !WRITE(*,"(A,ES20.10)") "-deb- in exchange_z_general, u = ", DBLE(SUM(u*CONJG(u)))
    !PRINT*, "size = ", SIZE(u,1),SIZE(u,2)
    bdesc%exchange_direction = 3 ! set to z direction
    !bdesc%communicator = my_comm_z
    n_exchanges = n_dim_2/bdesc%count
    DO i_exchange=0,n_exchanges-1
       !PRINT*,i_exchange,1+i_exchange*bdesc%count,(i_exchange+1)*bdesc%count
       CALL exchange_general(bdesc,u(:,1+i_exchange*bdesc%count:(i_exchange+1)*bdesc%count))
    END DO
    !WRITE(*,"(A,ES20.10)") "-deb- after periodic exchange, u = ", DBLE(SUM(u*CONJG(u)))

    ! Here we apply the physical boundary conditions as post-process of the exchange.
    ! It is only applied to the outer processes.
    IF (my_pez.EQ.0) THEN
       CALL post_apply_boundary_condition_z(bdesc,u,2)
    END IF
    IF (my_pez.EQ.n_procs_z-1) THEN
       CALL post_apply_boundary_condition_z(bdesc,u,1)
    END IF
    ! End of application of physical boundary conditions
    PERFOFF
  END SUBROUTINE exchange_z_general

  ! =======================================================
  !  now we add the physical boundary conditions
  ! =======================================================
  SUBROUTINE post_apply_boundary_condition_z(bdesc,u,side)
    TYPE(BoundaryDescription) :: bdesc
    COMPLEX, DIMENSION(0:,0:) :: u
    INTEGER :: side
    
    ! Local variables
    INTEGER :: i_index, i_subarray, start_index, end_index

    IF (side.EQ.1) THEN
       ! modify the upper boundary
       !PRINT*,"Modifying the upper boundary."
       DO i_index=0,SIZE(u,2)-1
          DO i_subarray = 1,bdesc%n_upper_subarrays
             start_index = bdesc%innerlast+(i_subarray-1)*bdesc%subarray_size+1
             end_index   = bdesc%innerlast+i_subarray*bdesc%subarray_size
             CALL multiply_array(&
                  & u(start_index:end_index,i_index),&
                  & pb_phase_fac(:,:,side),bdesc%subarray_size)
          END DO
       END DO
    ELSE
       !PRINT*,"Modifying the lower boundary."
       ! modify the lower boundary
       DO i_index=0,SIZE(u,2)-1
          DO i_subarray = 1,bdesc%n_lower_subarrays
             start_index = (i_subarray-1)*bdesc%subarray_size
             end_index   = i_subarray*bdesc%subarray_size-1
             CALL multiply_array(&
                  & u(start_index:end_index,i_index),&
                  & pb_phase_fac(:,:,side), bdesc%subarray_size)
          END DO
       END DO
    END IF
  END SUBROUTINE post_apply_boundary_condition_z


  SUBROUTINE multiply_array(u,phase_fac,end_index)
    INTEGER,intent(IN) :: end_index
    COMPLEX, DIMENSION(0:end_index-1),intent(INOUT) :: u
    COMPLEX, DIMENSION(0:end_index-1),intent(IN) :: phase_fac
    
    !n_multiplications = SIZE(u)/SIZE(phase_fac)
    !DO i_multiplication = 1,n_multiplications
    u = u * phase_fac
    !u = u((i_multiplication-1)*SIZE(phase_fac):i_multiplication*SIZE(phase_fac))*phase_fac
    !END DO
  END SUBROUTINE multiply_array


!--------------------- Exchange for equilibrium quantities --------------------
  Subroutine exchange_z_3d_equil(u,n0q0)
    real, dimension(pi1:pi2,pj1:pj2,lbz:ubz),intent(inout) :: u
    real, intent(in) :: n0q0

    Integer:: i,k
    Integer:: tag, dest_pes, recv_pes, nelem
    Integer:: stat(MPI_STATUS_SIZE), ierr
    real, dimension(pi1:pi2,pj1:pj2, 1:nzb):: tmps_r
    complex,dimension(:),allocatable:: u_c, u_fft

    If (n_procs_z.Eq.1) Then       
       u(:,:,lbz:lk1-1)=u(:,:,lk2-nzb+1:lk2)
       u(:,:,lk2+1:ubz)=u(:,:,lk1:lk1+nzb-1)
    Else
       nelem = pmi0*pj0*nzb
       dest_pes = mod(my_pez+1,n_procs_z)
       recv_pes = mod(my_pez-1+n_procs_z,n_procs_z)
       tag = 333
       Call mpi_sendrecv(u(pi1,pj1,lk2-nzb+1), nelem, &
            &MPI_REAL_TYPE, dest_pes, tag,&
            &tmps_r(pi1,pj1,1), nelem, MPI_REAL_TYPE, recv_pes, tag,&
            &mpi_comm_z, stat, ierr)
       u(:,:,lbz:lk1-1) = tmps_r
       dest_pes = mod(my_pez-1+n_procs_z,n_procs_z)
       recv_pes = mod(my_pez+1,n_procs_z)
       tag = 340
       Call mpi_sendrecv(&
            u(pi1,pj1,lk1), nelem, MPI_REAL_TYPE, dest_pes, tag,&
            tmps_r(pi1,pj1,1), nelem, MPI_REAL_TYPE, recv_pes, tag,&
            mpi_comm_z, stat, ierr)
       u(:,:,lk2+1:ubz) = tmps_r
    end if

    !include toroidal shift for nonaxisymmetric equilibria
    if (.not.y_local.and.(my_pez==0.or.my_pez==n_procs_z-1)) then
       allocate(u_c(li1:li2),u_fft(li1:li2))
       call initialize_fourier_x_1d
       !left z-boundary
       if (my_pez==0) then
          do k=lbz,lk1-1
             u_c=u(:,pj1,k)
             call to_fourier_x_1d(u_c,u_fft)
             do i=li1,hky
                u_fft(i)=u_fft(i)*exp(imag*2*pi*i*n0q0)
             enddo
             do i=hky+1,li2
                u_fft(i)=u_fft(i)*exp(imag*2*pi*(i-nky0)*n0q0)
             enddo
             call to_real_x_1d(u_fft,u_c)
             u(:,pj1,k)=u_c
          enddo
       endif
       !right z-boundary
       if (my_pez==n_procs_z-1) then
          do k=lk2+1,ubz
             u_c=u(:,pj1,k)
             call to_fourier_x_1d(u_c,u_fft)
             do i=li1,hky
                u_fft(i)=u_fft(i)*exp(-imag*2*pi*i*n0q0)
             enddo
             do i=hky+1,li2
                u_fft(i)=u_fft(i)*exp(-imag*2*pi*(i-nky0)*n0q0)
             enddo
             call to_real_x_1d(u_fft,u_c)
             u(:,pj1,k)=u_c
          enddo
       endif
       deallocate(u_c,u_fft)
       call finalize_fourier_x_1d
    endif

  end subroutine exchange_z_3d_equil

END MODULE boundary_exchange_z
