#include "switches.h"
!>Solver for the (nonlocal) field equations
!!
!! We solve the two gyrokinetic field equations (quasi-neutrality
!! and Ampere's law) in this module.
#include "intrinsic_sizes.h"
#include "redef.h"
#undef FLMSLV

MODULE fieldsolver_df
  USE spectype_mod
  USE discretization
  USE coordinates
  use par_in, only: trap_pass, debye2, beta, tau, rad_bc_type, &
       &spec, del_phi, only_Er, delzonal, add_zonal_phi
  USE par_other, ONLY: n_fields, pi, p_has_0_mode
  use geometry, only: geom
  USE vel_space
  Use mpi
  USE communications
  USE gyro_average_df_mod, ONLY: get_nDeriv_base,basex, estimate_gyromatrix_maxbw, xgrid !, get_gyromatrix_dagger
  USE debug_output!, ONLY: output_data
  USE MatrixModule
  USE BandedMatrixModule
  use Grid1DModule
  use DerivativeMatrixModule
  use flr_corr
  use hybrid, only: trap
  use adiabatic_response
  use spatial_averages

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: field_solve_df, &
       &initialize_field_solve_df, finalize_field_solve_df
  public :: estimate_memory_fieldsolver_df

  INTERFACE output_local_sum_jkn
     MODULE PROCEDURE output_local_sum_jkn_banded, output_local_sum_jkn_full
  END INTERFACE

  ! constants
  COMPLEX, PARAMETER :: imag=(0.0,1.0)
  LOGICAL,PARAMETER :: transposed = .TRUE.
#ifdef G_DAGGER
  CHARACTER, PARAMETER :: G_op = 'C' !take conjugate transpose 
  !e.g., for first matrix in dot_multiply or in gyro averages
  !if G_op = 'C'
#else
  CHARACTER, PARAMETER :: G_op = 'N' !take original matrix
#endif

  INTEGER :: derivative_order=4

  ! private variables II: calculated from the given variables
  INTEGER :: nCoeffs !< number of independent base functions in x direction
  LOGICAL :: periodic_boundary
  ! private variables III: real proprietary variables of the module

  TYPE(matrix), DIMENSION(:,:), ALLOCATABLE :: inv_qn_matrix, inv_amp_matrix
  TYPE(Matrix), SAVE :: rhsvec, solvec
  TYPE(DerivativeMatrix), DIMENSION(:), allocatable :: derivmat, derivmat_vN
  integer :: nDeriv

CONTAINS

  !> This is the main routine of the module. 
  !! It is called to calculate the two fields out of the distribution function. 
  !! \param a_g_1 The distribution function g
  !! \param a_emfields An array for the electromagnetic fields, which will be calculated out of the given distribution function.
  SUBROUTINE field_solve_df(moments,a_emfields)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,1:n_fields), INTENT(INOUT) :: moments
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), INTENT(OUT) :: a_emfields
    complex, dimension(lg1:lg2):: phi_fs_avg
    integer :: i,j,k

    PERFON('fldsolve')

    CALL solve_fieldequations(inv_qn_matrix, moments(:,:,:,1), a_emfields(:,:,:,1))  
    IF (n_fields.gt.1) CALL solve_fieldequations(inv_amp_matrix, moments(:,:,:,2), a_emfields(:,:,:,2)) 


    !below: modifications for neoclassical and zonal flow studies
    if (del_phi) a_emfields(:,:,:,1)=0.0
    if (delzonal) then
       call flux_surface_average(a_emfields(:,:,:,1),.true.,phi_fs_avg)
       if (.not.x_local) then
          if (p_has_0_mode) then
             do k=lbz,ubz
                a_emfields(:,lj1,k,1)=a_emfields(:,lj1,k,1)-phi_fs_avg(:)
             enddo
          endif
       else if (.not.y_local) then
          do k=lbz,ubz
             do j=lj1,lj2
                do i=li1,li2
                   a_emfields(i,j,k,1)=a_emfields(i,j,k,1)-phi_fs_avg(j)
                enddo
             enddo
          enddo
       endif
    endif       
    if (only_Er) then
       if (y_local) then
          if (p_has_0_mode) then
             call flux_surface_average(a_emfields(:,:,:,1),.true.,phi_fs_avg)
             do k=lbz,ubz
                a_emfields(:,lj1,k,1)=phi_fs_avg(:)
             enddo
             a_emfields(:,lj1+1:lj2,:,1)=0.0
          else
             a_emfields(:,:,:,1) = 0.0
          endif
       else
          call flux_surface_average(a_emfields(:,:,:,1),.true.,phi_fs_avg)
          do k=lbz,ubz
             do i=pi1,pi2
                a_emfields(i,:,k,1)=phi_fs_avg(:)
             enddo
          enddo
       endif
    endif
    if (abs(add_zonal_phi).gt.epsilon(0.0)) then
       if (y_local) then
          if (p_has_0_mode) then
             a_emfields(:,lj1,:,1) = a_emfields(:,lj1,:,1) + add_zonal_phi
          endif
       else
          a_emfields(:,lj1+1:lj2,:,1) = a_emfields(:,lj1+1:lj2,:,1) + &
               &add_zonal_phi
       endif
    endif

    PERFOFF

  END SUBROUTINE field_solve_df


  !>Solves the field equations
  !!
  !!\param mat inverse left hand side operator
  !!\param rhs the right hand side of the field equations
  !!\param sol the resulting potential/field
  SUBROUTINE solve_fieldequations(mat, rhs, sol)
    TYPE(Matrix), DIMENSION(lj1:lj2,lk1:lk2) :: mat
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2), INTENT(INOUT),TARGET :: rhs
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(INOUT),TARGET :: sol

    INTEGER :: j,k,ierr
    complex :: mean_rhs_loc, mean_rhs

    IF ((periodic_boundary).AND.(p_has_0_mode)) THEN
       ! delete the mean of the rhs
       DO k=lk1,lk2
          mean_rhs_loc = SUM(rhs(:,lj1,k))
          Call MPI_Allreduce(mean_rhs_loc, mean_rhs, 1, MPI_COMPLEX_TYPE,&
               MPI_SUM, mpi_comm_x, ierr)
          rhs(:,lj1,k) = rhs(:,lj1,k)-mean_rhs/ni0
          If (li1 == 0) rhs(li1,lj1,k) = (0.0,0.0)
       END DO
    END IF

    ! Solve the LU-factored equation systems.
    ! This is done in 3 distinct phases in order to run in parallel.
    ! See the source code of LU_solve for details!

    DO k=lk1,lk2
       DO j=lj1,lj2
          CALL attach(rhsvec,rhs(:,j,k))
          CALL LU_solve(mat(j,k),rhsvec,solvec,1)
       END DO
    END DO
    DO k=lk1,lk2
       DO j=lj1,lj2
          CALL LU_solve(mat(j,k),rhsvec,solvec,2)
       END DO
    END DO
    DO k=lk1,lk2
       DO j=lj1,lj2
          CALL attach(solvec,sol(li1:li2,j,k))
          CALL LU_solve(mat(j,k),rhsvec,solvec,3)
       END DO
    END DO

  END SUBROUTINE solve_fieldequations

  !> Function to estimate the memory consumption of the fieldsolver_df module
  !! It uses the function estimate_gyromatrix_maxbw from the gyro_average_df module to
  !! estimate the bandwidth of the gyromatrix. The bandwidth of the field matrices is then estimated
  !! as twice the bandwidths of the gyromatrix.
  !! \param xmatrix_dimension Number of points in x direction globally.
  FUNCTION estimate_memory_fieldsolver_df(xmatrix_dimension) RESULT(memory_need)
    integer :: xmatrix_dimension
    REAL :: memory_need
    
    ! Local variables
    REAL :: mem_xmatrix
    REAL :: mem_banded_xmatrix
    INTEGER :: gm_maxbw, gamma_maxbw

    mem_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*xmatrix_dimension*xmatrix_dimension/n_procs_x
    
    gm_maxbw = estimate_gyromatrix_maxbw()
    ! gamma0 and gamma_mu are squared gyromatrices integrated over mu. The bandwidth
    ! is therefore the double of the gyromatrix bandwidth
    gamma_maxbw = 2*(gm_maxbw-1)+1

    mem_banded_xmatrix = static_size_of_BandedMatrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*xmatrix_dimension*gamma_maxbw/n_procs_x

    memory_need = 0.0
    ! modular variables

    ! inv_qn_matrix, inv_amp_matrix
    memory_need = memory_need + n_fields*lj0*lk0*mem_xmatrix

    ! rhsvec, solvec
    memory_need = memory_need + 2.*static_size_of_matrix()/1024./1024.

    ! local variables in field_solve_df
    ! dens, curr
    memory_need = memory_need + 2.*li0*lj0*lk0*ln0*SIZE_OF_COMPLEX_MB
    ! qn_rhs, amp_rhs
    memory_need = memory_need + 2*lijk0*SIZE_OF_COMPLEX_MB
    ! qn_rhs_avg, phi_fs_avg
    memory_need = memory_need + 2*li0*SIZE_OF_COMPLEX_MB

    ! local variables in calculate_density, calculate_current
    ! densmu, currmu
    memory_need = memory_need + 2.*lijk0*lm0*ln0*SIZE_OF_COMPLEX_MB
    ! bardensmu, barcurrmu
    memory_need = memory_need + 2.*lijk0*SIZE_OF_COMPLEX_MB

    ! local variables in calculate_premats
    ! p_doubleprime
    memory_need = memory_need + mem_xmatrix
    ! inv_p_doubleprime, temparr
    memory_need = memory_need + mem_xmatrix
    ! xgrid
    memory_need = memory_need + size_of_grid1d(ni0)/1024./1024.
    ! derivmat
    memory_need = memory_need + 2*mem_xmatrix
    ! derivmat_vN
    IF (rad_bc_type.eq.2) memory_need = memory_need + 2*mem_xmatrix
    ! laplacian_perp
    memory_need = memory_need + mem_xmatrix
    IF (n_fields.GT.1) THEN
       memory_need = memory_need + ln0*mem_banded_xmatrix
    END IF

  END FUNCTION estimate_memory_fieldsolver_df

  SUBROUTINE initialize_field_solve_df
    INTEGER :: j, k, ideriv
    
    PERFON('fld_init')

    periodic_boundary = (rad_bc_type.EQ.0)

    nCoeffs = ni0

    CALL my_barrier()
    
    CALL initialize(rhsvec,ni0)
    CALL initialize(solvec,ni0)

    ALLOCATE(inv_qn_matrix(lj1:lj2,lk1:lk2))
    IF (n_fields.GT.1) ALLOCATE(inv_amp_matrix(lj1:lj2,lk1:lk2))

    ! Please note: The data of inv_qn_matrix and inv_amp_matrix is not allocated here.
    ! It will be allocated just before usage since these matrices are LU factored
    ! immediatly and thus the regular matrix data can be deallocated directly
    ! after allocation (in favour of the smaller LU factorization).
    DO k=lk1,lk2
       DO j=lj1,lj2
          CALL initialize(inv_qn_matrix(j,k),nCoeffs,nCoeffs)
          IF (n_fields.GT.1) &
             & CALL initialize(inv_amp_matrix(j,k),nCoeffs,nCoeffs)
       END DO
    END DO

    call initialize_flr_corr_moms

    
    !IF ((n_fields.GT.1).OR.(debye2.GT.0)) THEN
    ! allocate the derivative matrices for all derivatives and also the 
    ! integral matrix contains the theta integral
    nDeriv = get_nDeriv_base()
    
    ALLOCATE(derivmat(1:nDeriv))
    IF (rad_bc_type.EQ.2) ALLOCATE(derivmat_vN(1:nDeriv))
    
    DO iDeriv=1,nDeriv
       CALL initialize(derivmat(iDeriv),xgrid,derivative_order,transposed)
       IF (rad_bc_type.EQ.2) THEN
          CALL calculate(derivmat(iDeriv),iDeriv,1)
          CALL initialize(derivmat_vN(iDeriv),xgrid,derivative_order,transposed)
          CALL calculate(derivmat_vN(iDeriv),iDeriv,rad_bc_type)
       ELSE
          CALL calculate(derivmat(iDeriv),iDeriv,rad_bc_type)
       ENDIF
    END DO

    CALL calculate_premats

    PERFOFF
  END SUBROUTINE initialize_field_solve_df


  !> Removes and deallocates all data structures
  !! used in fieldsolver_df
  SUBROUTINE finalize_field_solve_df
    INTEGER :: j,k,ideriv

    call finalize(rhsvec)
    call finalize(solvec)

    call finalize_flr_corr_moms

    DO k=lk1,lk2
       DO j=lj1,lj2
          CALL finalize(inv_qn_matrix(j,k))
          IF (n_fields.GT.1) CALL finalize(inv_amp_matrix(j,k))
       END DO
    END DO

    DEALLOCATE(inv_qn_matrix)
    IF (n_fields.gt.1) DEALLOCATE(inv_amp_matrix)

    DO iDeriv=1,nDeriv
       CALL finalize(derivmat(iDeriv))
       IF (rad_bc_type.EQ.2) CALL finalize(derivmat_vN(iDeriv))
    END DO
    DEALLOCATE(derivmat)
    IF (rad_bc_type.EQ.2) DEALLOCATE(derivmat_vN)

  END SUBROUTINE finalize_field_solve_df


  !========================================
  ! internal subroutines follow
  !========================================

  !> Calculation of the matrices needed for solution of 
  !! field equations. 
  SUBROUTINE calculate_premats

    ! Local variables
    TYPE(BandedMatrix), DIMENSION(:,:,:), POINTER :: mom00_flrcorr, &
         &gamma0_num_matrix
    TYPE(Matrix) :: p_doubleprime

    TYPE(Matrix) :: temparr   !,fs_avg_p_doubleprime
    TYPE(Matrix) :: laplacian_perp
    INTEGER :: i,ibase,j,k,n !, info, np_bw
    INTEGER :: irow,icol
!    Integer:: ierr
#ifdef APPROX_METRIC
    real :: gii,gij
#endif

    PERFON('cal_prem')

    mom00_flrcorr=>get_mom00_flrcorr()
    IF (n_fields.gt.1) gamma0_num_matrix=>get_gamma0_num_matrix()

    CALL initialize(temparr,nCoeffs,nCoeffs)
    call allocate(temparr)


    IF ((n_fields.GT.1).OR.(debye2.GT.0)) THEN  !build laplacian_perp matrix
       CALL initialize(laplacian_perp,nCoeffs,nCoeffs)
       CALL ALLOCATE(laplacian_perp)
    END IF


    CALL initialize(p_doubleprime,nCoeffs,nCoeffs)
    CALL ALLOCATE(p_doubleprime)

    PERFON('cp_p1')
    ! calculate the square of this matrix
    DO k=lk1,lk2
       DO j=lj1,lj2

          CALL set_zero(p_doubleprime)

          !build laplacian_perp matrix, if necessary
          IF ((n_fields.GT.1).OR.(debye2.GT.0)) THEN
             CALL set_zero(laplacian_perp)
             DO i=my_pex*li0+1,(my_pex+1)*li0
                DO ibase=1,nCoeffs
#ifdef APPROX_METRIC
                   !here we take only the hermitian part of the Laplace operator
                   !by computing A=1/2*(A+A^H)
                   gii=(geom%gii(i-1+pi1gl,pj1,k)+geom%gii(ibase-1+pi1gl,pj1,k))/2.
                   gij=(geom%gij(i-1+pi1gl,pj1,k)+geom%gij(ibase-1+pi1gl,pj1,k))/2.
                   IF ((rad_bc_type.EQ.2).AND.p_has_0_mode.AND.(j==lj1)) THEN
                      CALL add_value(laplacian_perp,i,ibase,&
                           &gii*mat_get_value(derivmat_vN(2),i,ibase) &
                           & + 2.0*gij&
                           &   * mat_get_value(derivmat_vN(1),i,ibase)*imag*kj(j))
                   ELSE
                      CALL add_value(laplacian_perp,i,ibase,&
                           &gii*mat_get_value(derivmat(2),i,ibase) &
                           & + 2.0*gij &
                           &   * mat_get_value(derivmat(1),i,ibase)*imag*kj(j))
                   ENDIF
#else
                   ! this loop is the analogue of multiplying diag(gii),diag(gij)
                   ! the following line does imply to have the derivmat stored in not transposed form
                   ! this is because it is not implemented if the mat_get_value function references an entry
                   ! which is stored on another processor
                   IF ((rad_bc_type.EQ.2).AND.p_has_0_mode.AND.(j==lj1)) THEN
                      CALL add_value(laplacian_perp,i,ibase,&
                           &geom%gii(i-1+pi1gl,pj1,k)*mat_get_value(derivmat_vN(2),i,ibase) &
                           & + 2.0*geom%gij(i-1+pi1gl,pj1,k) &
                           &   * mat_get_value(derivmat_vN(1),i,ibase)*imag*kj(j))
                   ELSE
                      CALL add_value(laplacian_perp,i,ibase,&
                           &geom%gii(i-1+pi1gl,pj1,k)*mat_get_value(derivmat(2),i,ibase) &
                           & + 2.0*geom%gij(i-1+pi1gl,pj1,k) &
                           &   * mat_get_value(derivmat(1),i,ibase)*imag*kj(j))
                   ENDIF
#endif

                END DO
                call commit_values(laplacian_perp)
                CALL add_value(laplacian_perp,i,i,&
                     &-geom%gjj(i-1+pi1gl,pj1,k)*kj(j)**2)
             END DO
             call commit_values(laplacian_perp)
          ENDIF

          DO n=ln1,ln2
             !build p_doubleprime matrix
             DO i=my_pex*li0+1,(my_pex+1)*li0
                DO ibase=1,nCoeffs
                   CALL add_value(p_doubleprime,i,ibase,&
                        & -spec(n)%charge * spec(n)%dens * &
                        & mat_get_value(mom00_flrcorr(j,k,n),i,ibase))
                END DO
             ENDDO
             call commit_values(p_doubleprime)
          END DO
          ! parallel sum over the species
          CALL my_sum(p_doubleprime,mpi_comm_spec)

          IF (debye2.GT.0) THEN
             CALL multiply_matrix_with_scalar(laplacian_perp,debye2,temparr)
             CALL subtract_matrix(p_doubleprime,temparr)
          END IF

         
          IF (adiabatic_electrons.or.trap_pass) THEN
             !flux surface average of p_doubleprime, this is needed for the 
             !computation of the flux surface averaged phi
             IF (y_local) THEN
                !global in x
                IF (p_has_0_mode.AND.(j.EQ.lj1)) THEN
                   ! calculate the flux surface average of p_doubleprime
                   DO icol=1,nCoeffs
                      DO irow=my_pex*li0+1,(my_pex+1)*li0
                         !! old field solver (assuming vector instead of operator type flux surface
                         !! average)
#ifdef FLMSLV
                         CALL add_value(inv_p_doubleprime_ions_m,irow,icol,&
                              &mat_get_value(p_doubleprime,irow,icol) &
                              &*geom%jacobian(irow-1+pi1gl,pj1,k)/geom%avg_jaco_yz(irow-1+pi1gl)&
                              &*geom%avg_jaco_yz(icol-1+pi1gl)/geom%jacobian(icol-1+pi1gl,pj1,k)/nz0)
#else
                         CALL add_value(inv_p_doubleprime_ions_m,irow,icol,&
                              &mat_get_value(p_doubleprime,irow,icol)*geom%jacobian(irow-1+pi1gl,pj1,k)&
                              /(nz0*geom%avg_jaco_yz(irow-1+pi1gl)))
#endif
                      END DO
                   END DO
                   call commit_values(inv_p_doubleprime_ions_m)
                END IF
             ELSE
                DO icol=1,nCoeffs
                   DO irow=my_pex*li0+1,(my_pex+1)*li0
#ifdef FLMSLV
                      inv_p_doubleprime_ions_v(j) = inv_p_doubleprime_ions_v(j) &
                           +mat_get_value(p_doubleprime,irow,icol)*geom%jacobian(irow-1+pi1gl,pj1,k)&
                           /geom%jacobian(icol-1+pi1gl,pj1,k)/nz0/nky0
#else
                      inv_p_doubleprime_ions_v(j) = inv_p_doubleprime_ions_v(j) &
                           +mat_get_value(p_doubleprime,irow,icol)*geom%jacobian(irow-1+pi1gl,pj1,k)&
                           /(nz0*nky0*geom%avg_jaco_yz(pj1))
#endif
                   END DO
                END DO
             END IF
          END IF

          ! p_doubleprime now contains the value for P'' from my script. It has now 
          ! to be inverted and multiplied with the gyromatrix to give the final
          ! premat_dens. For the inversion, we first have to do a LU factorization
          ! with a LAPACK routine. For the matrices with periodic boundaries, we
          ! have to use a general matrix solver, but later, we can hopefully use
          ! a banded matrix solver which is much faster.
          
          CALL allocate(inv_qn_matrix(j,k))
          inv_qn_matrix(j,k) = p_doubleprime

          IF (adiabatic_electrons) THEN
             DO i=1,nCoeffs
                CALL add_value(inv_qn_matrix(j,k),i,i,&
                     &dens_electrons(i-1+pi1gl)/temp_electrons(i-1+pi1gl))
             END DO
          ELSE IF (trap_pass) THEN
             DO i=1,nCoeffs
                CALL add_value(inv_qn_matrix(j,k),i,i,&
                     &(1-trap(i-1+pi1gl,k))*dens_electrons(i-1+pi1gl)/temp_electrons(i-1+pi1gl))
             END DO
          ELSE
             IF (n_spec.EQ.1) THEN !adiabatic ions
                DO i=1,nCoeffs
                   CALL add_value(inv_qn_matrix(j,k),i,i,&
                        &spec(0)%dens_prof(i-1+pi1gl)/spec(0)%temp_prof(i-1+pi1gl)*tau)
                END DO
             ENDIF
          END IF
          call commit_values(inv_qn_matrix(j,k))

          ! before we can continue, we have to modify the matrix in the periodic
          ! boundary case as there the matrix does not have full rank. This comes
          ! from the fact, that the solution is determined only up to a constant.
          ! We now modify the matrix so that the average of the solution is zero. 
          ! This is done by setting the first row to ones and the corresponding
          ! rhs entry to zero.
          IF (p_has_0_mode.AND.periodic_boundary.AND.(j.EQ.lj1)) THEN
             DO ibase=1,nCoeffs
                CALL set_value(inv_qn_matrix(lj1,k),1,ibase,1.0)
             END DO
             CALL commit_values(inv_qn_matrix(j,k))
          END IF

          call LU_factor(inv_qn_matrix(j,k))

          ! the same procedure for the matrix for Ampere's law
          IF (n_fields.GT.1) THEN
             CALL allocate(inv_amp_matrix(j,k))
             CALL set_zero(inv_amp_matrix(j,k))
             
             DO n=ln1,ln2
                DO i=my_pex*li0+1,(my_pex+1)*li0
                   DO ibase=1,nCoeffs
                      CALL add_value(inv_amp_matrix(j,k),i,ibase,&
                           &mat_get_value(gamma0_num_matrix(j,k,n),i,ibase))
                   END DO
                END DO
                CALL commit_values(inv_amp_matrix(j,k))
             END DO
             
             CALL multiply_matrix_with_scalar(inv_amp_matrix(j,k),beta)
             CALL my_sum(inv_amp_matrix(j,k),mpi_comm_spec)
             
             CALL subtract_matrix(inv_amp_matrix(j,k),laplacian_perp)

             IF ((periodic_boundary).AND.(p_has_0_mode).AND.(j.EQ.lj1)) THEN
                DO ibase=1,nCoeffs
                   CALL set_value(inv_amp_matrix(lj1,k),1,ibase,1.0)
                END DO
             END IF

             CALL LU_factor(inv_amp_matrix(j,k))
             
          ENDIF !n_fields>1
          
       END DO ! j
    END DO ! k

    IF ((n_fields.GT.1).OR.(debye2.GT.0)) THEN
       CALL finalize(laplacian_perp)
    END IF

    call finalize(temparr)
    call finalize(p_doubleprime)


    PERFOFF

    PERFON('cp_p2')
    IF (adiabatic_electrons.or.trap_pass) THEN
       ! In the adiabatic electron case, the matrix p_doubleprime now contains
       ! only the ion contributions, as the electrons are not a separate species. 
       ! For the solution of the adiabatic electron case, we need the flux surface
       ! average of the p_doubleprime matrix.
       IF (.NOT.x_local) THEN
          IF (p_has_0_mode) THEN
             CALL my_sum(inv_p_doubleprime_ions_m, mpi_comm_z)
             IF (periodic_boundary) THEN
                DO ibase=1,nCoeffs
                   CALL set_value(inv_p_doubleprime_ions_m,1,ibase,1.0)
                END DO
             END IF
             call invert(inv_p_doubleprime_ions_m)
          END IF
       ELSE IF (.not.y_local) THEN
          !sum over parallel y nodes
          CALL sum_to_all_complex(inv_p_doubleprime_ions_v, lg0, mpi_comm_x)
          CALL sum_to_all_complex(inv_p_doubleprime_ions_v, lg0, mpi_comm_z)
          !invert
          inv_p_doubleprime_ions_v=1./inv_p_doubleprime_ions_v
       END IF
    END IF !adiabatic_electrons or trap_pass
    PERFOFF

    PERFOFF
  END SUBROUTINE calculate_premats



!------------------------------------------------------------------------------------
!---------------------------- DEBUGGING & TESTING -----------------------------------
!------------------------------------------------------------------------------------

  SUBROUTINE output_local_sum_jkn_banded(bmat,label)
    TYPE(BandedMatrix),DIMENSION(lj1:lj2,lk1:lk2,ln1:ln2) :: bmat
    CHARACTER(len=*) :: label

    ! Local variables
    real :: local_sum
    INTEGER :: j,k,n

    local_sum = 0.0D0
    DO n=ln1,ln2
       DO k=lk1,lk2
          DO j=lj1,lj2
             local_sum = local_sum + get_local_abs_square_sum(bmat(j,k,n))
          END DO
       END DO
    END DO
    WRITE(*,"(I3,3A,ES20.8)") mype,": ",label," = ",local_sum
  END SUBROUTINE output_local_sum_jkn_banded

  SUBROUTINE output_local_sum_jkn_full(mat,label)
    TYPE(Matrix),DIMENSION(lj1:lj2,lk1:lk2,ln1:ln2) :: mat
    CHARACTER(len=*) :: label

    ! Local variables
    real :: local_sum
    INTEGER :: j,k,n

    local_sum = 0.0D0
    DO n=ln1,ln2
       DO k=lk1,lk2
          DO j=lj1,lj2
             local_sum = local_sum + get_local_abs_square_sum(mat(j,k,n))
          END DO
       END DO
    END DO
    WRITE(*,"(I3,3A,ES20.8)") mype,": ",label," = ",local_sum
  END SUBROUTINE output_local_sum_jkn_full

  SUBROUTINE output_local_sum_jk(mat,label)
    TYPE(Matrix),DIMENSION(lj1:lj2,lk1:lk2) :: mat
    CHARACTER(len=*) :: label

    ! Local variables
    real :: local_sum
    INTEGER :: j,k

    local_sum = 0.0D0
    DO k=lk1,lk2
       DO j=lj1,lj2
          local_sum = local_sum + get_local_abs_square_sum(mat(j,k))
       END DO
    END DO
    WRITE(*,"(I3,3A,ES20.8)") mype,": ",label," = ",local_sum
  END SUBROUTINE output_local_sum_jk

#ifdef WITH_TESTROUTINES
  !> To test, if the fieldsolver works correctly, we can use the calculated
  !! fields and put them into the matrix-vector equation and get out
  !! the original rhs. This is done in this subroutine. We assume that the 
  !! fields (phi and psi) are already calculated.
  SUBROUTINE test_fieldsolve(p_g_1,p_phi,p_psi)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),INTENT(IN) :: p_g_1
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(IN) :: p_phi,p_psi

    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: locresd,locresc
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2) :: dens,curr
    LOGICAL :: OUTPUT=.TRUE., op
    INTEGER :: thisunit,i,j,k,n

    DO k=lk1,lk2
       DO j=lj1,lj2
          locresd(:,j,k) = MATMUL(qn_matrix(:,:,j,k),p_phi(li1:li2,j,k))
          locresc(:,j,k) = MATMUL(amp_matrix(:,:,j,k),p_psi(li1:li2,j,k))
       END DO
    END DO

    CALL calculate_density(p_g_1, dens)
    CALL calculate_current(p_g_1, curr)
    
    IF (OUTPUT) THEN
       thisunit=30
       DO 
          INQUIRE(thisunit,opened=op)
          IF (op) THEN
             thisunit = thisunit+1
          ELSE 
             EXIT
          END IF
       END DO
       
       OPEN(thisunit,file='./diff_fieldsolve.dat')
       WRITE(thisunit,"(5I3)") li2-li1+1,lj2-lj1+1,lk2-lk1+1,ln2-ln1+1
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li1,li2
                WRITE(thisunit,"(2ES20.10)") locresd(i,j,k)
             END DO
          END DO
       END DO
       DO n=ln1,ln2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=li1,li2
                   WRITE(thisunit,"(2ES20.10)") dens(i,j,k,n)
                END DO
             END DO
          END DO
       END DO
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li1,li2
                WRITE(thisunit,"(2ES20.10)") locresc(i,j,k)
             END DO
          END DO
       END DO
       DO n=ln1,ln2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=li1,li2
                   WRITE(thisunit,"(2ES20.10)") curr(i,j,k,n)
                END DO
             END DO
          END DO
       END DO
       CLOSE(thisunit)
    END IF

  END SUBROUTINE test_fieldsolve
#endif

END MODULE fieldsolver_df
