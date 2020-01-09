#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"
Module flr_corr
  Use par_mod
  use communications
  USE gyro_average_df_mod
  use MatrixModule
  USE BandedMatrixModule
  USE VectorModule, ONLY: Vector, static_size_of_vector
  USE geometry, ONLY: geom
  use vel_space
  use aux_func, only: gamma0, gamma1
  use hybrid, only: Bfieldmax, trap
  use equilibrium_fields, only: dens_co

  IMPLICIT NONE
  PRIVATE

  PUBLIC:: initialize_flr_corr_moms, get_mom00_flrcorr, &
       &get_gamma0_num_matrix, finalize_flr_corr_moms

  PUBLIC:: mem_est_flr_corrs, initialize_flr_corrs, &
       &finalize_flr_corrs, correct_for_FLR, correct_for_FLR_bpar

  PUBLIC :: n_flr_corrs

  TYPE(BandedMatrix), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: mom00_flrcorr, &
       &mom01_flrcorr, mom20_flrcorr, gamma0_num_matrix
  logical:: exist_gammas=.false.

#ifdef GDAGGER_G
  CHARACTER, PARAMETER :: G_op = 'C' !take conjugate transpose 
  !e.g., for first matrix in dot_multiply or in gyro averages
  !if G_op = 'C'
#else
  CHARACTER, PARAMETER :: G_op = 'N'
#endif

  COMPLEX, DIMENSION(:,:,:,:,:), ALLOCATABLE:: corr_mom_00_df,&
       &corr_mom_20_df,corr_mom_01_df
  TYPE(Vector),SAVE :: resvec, tmpvec

  REAL, DIMENSION(:,:,:,:), ALLOCATABLE:: corr_mom_00_ff, &
       & corr_mom_01_ff,corr_mom_20_ff
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: corr_mom_00_ff_bpar, &
    corr_mom_01_ff_bpar, corr_mom_20_ff_bpar
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: corr_momI_00_ff, &
    corr_momI_01_ff, corr_momI_00_ff_bpar, corr_momI_01_ff_bpar, &
    corr_momI_20_ff, corr_momI_20_ff_bpar
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: corr_mom_10_ff, &
    corr_mom_11_ff, corr_mom_30_ff
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: corr_mom_10_ff_bpar, &
    corr_mom_11_ff_bpar, corr_mom_30_ff_bpar
  REAL, DIMENSION(:,:), ALLOCATABLE :: aux_j0_flrcorr

  INTEGER :: n_flr_corrs = 3

CONTAINS
!!!******************************************************************!!!
!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of FLR correction terms
  Real Function mem_est_flr_corr_df(mem_req_in)
    real:: mem_req_in
    real:: mem_loc
    integer :: gm_maxbw, gamma_maxbw
    real :: mem_xmatrix
    real :: mem_banded_xmatrix

    mem_loc = mem_req_in
    ! get max bandwidth of gyromatrix
    gm_maxbw = estimate_gyromatrix_maxbw()
    ! max bandwidth of gamma matrices is double of gyromatrix
    gamma_maxbw = 2*(gm_maxbw-1)+1
    mem_banded_xmatrix = static_size_of_BandedMatrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*ni0*gamma_maxbw/n_procs_x
    mem_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*ni0*ni0/n_procs_x

    !mom_flrcorr and gamma_0_num matrices
    mem_loc = mem_loc + 4.*mem_banded_xmatrix*lj0*lk0*ln0
    !tmpvec, resvec
    mem_loc = mem_loc + 2.*static_size_of_vector()/(1024.)**2

    !local variables in initialize_FLR_corrs
    !single_gyromatrix, full_gyromatrix2, profilematrix, profdotgyro
    mem_loc = mem_loc + 4.*mem_xmatrix    
    !gyromatrix2
    mem_loc = mem_loc + mem_banded_xmatrix

    !in correct_for_FLR
    !loc_phi
    mem_loc = mem_loc + lij0*lz0*SIZE_OF_COMPLEX_MB
    !corrections
    mem_loc = mem_loc + lijk0*ln0*3*SIZE_OF_COMPLEX_MB

    mem_est_flr_corr_df=mem_loc
  End Function mem_est_flr_corr_df


  !> A getter function for the gamma0 matrix, needed for diagnostics and field solver.
  !! This function is necessary to use the modular gamma matrices outside the module.
  !! It returns a pointer to the gamma0 matrix as banded matrix.
  FUNCTION get_mom00_flrcorr() RESULT(res)
    TYPE(BandedMatrix),DIMENSION(:,:,:),POINTER :: res

    IF (exist_gammas) THEN
       res => mom00_flrcorr
    ELSE
       res => NULL()
    END IF
  END FUNCTION get_mom00_flrcorr

  !> A getter function for the gamma_0_num matrix, needed for the electromagnetic field solver.
  !! This function is necessary to use the modular gamma matrices outside the module.
  !! It returns a pointer to the gamma_mu matrix  as banded matrix.
  FUNCTION get_gamma0_num_matrix() RESULT(res)
    TYPE(BandedMatrix),DIMENSION(:,:,:),POINTER :: res

    IF (exist_gammas) THEN
       res => gamma0_num_matrix
    ELSE
       res => NULL()
    END IF
  END FUNCTION get_gamma0_num_matrix


!>Intialization of the FLR correction operators being applied on phi for different
!!velocity space moments: mom(ab)_flrcorr (first index is vpar, second index is mu)
!!The normalization is species dependent (!):
!!n_{0j} (x_0) * v_{Tj}^{a+2b}(x_0) * \rho_{ref} / L_{ref}
!!(compare sec. 2.3 in tbg's thesis)
!!
!!Note: 
!!(1) The mom11 and mom31 corrections are zero as long as we do not consider j_0\|
!!(2) We additionally compute a gamma0-like operator which appears in Ampere's law
!!    when transforming from f_1 to g_1. A vpar integration is performed numerically
!!    to avoid the Ampere cancellation problem
  SUBROUTINE initialize_flr_corr_moms
    ! Local variables
    TYPE(BandedMatrix), DIMENSION(:,:,:,:),POINTER :: gyromatrix, gyromatrix_dagger
    !TYPE(Matrix) :: single_gyromatrix, full_gyromatrix2,  profdotgyro
    TYPE(BandedMatrix) :: profdotgyro
    type(BandedMatrix) :: profilematrix
    TYPE(BandedMatrix) :: gyromatrix2
    INTEGER :: i,j,k,l,m,n,index
#ifndef GDAGGER_G
    INTEGER :: ibase
#endif
    INTEGER :: nCoeffs
    LOGICAL :: transposed=.true.
    COMPLEX :: tmpvar
    REAL :: tmpfac
#ifdef EXTERNAL_ERF
  EXTERNAL erf
#endif
    real :: erf

    nCoeffs=ni0
    ALLOCATE(mom00_flrcorr(lj1:lj2,lk1:lk2,ln1:ln2))
    ALLOCATE(mom01_flrcorr(lj1:lj2,lk1:lk2,ln1:ln2))
    ALLOCATE(mom20_flrcorr(lj1:lj2,lk1:lk2,ln1:ln2))
    IF (n_fields.gt.1) ALLOCATE(gamma0_num_matrix(lj1:lj2,lk1:lk2,ln1:ln2))

    DO n=ln1,ln2
       DO k=lk1,lk2
          DO j=lj1,lj2
             CALL initialize(mom00_flrcorr(j,k,n),nCoeffs,nCoeffs,transposed)
             CALL initialize(mom01_flrcorr(j,k,n),nCoeffs,nCoeffs,transposed)
             IF (n_fields.gt.1) CALL initialize(gamma0_num_matrix(j,k,n),nCoeffs,nCoeffs,transposed)
             CALL initialize(mom20_flrcorr(j,k,n),nCoeffs,nCoeffs,transposed)

             CALL ALLOCATE(mom00_flrcorr(j,k,n))
             CALL ALLOCATE(mom01_flrcorr(j,k,n))
             IF (n_fields.gt.1) CALL ALLOCATE(gamma0_num_matrix(j,k,n))
             CALL ALLOCATE(mom20_flrcorr(j,k,n))
          END DO
       END DO
    END DO

    !-------------------------- INITIALIZATION -------------------------------!
    gyromatrix => get_gyromatrix()
    IF (.NOT.ASSOCIATED(gyromatrix)) THEN
       PRINT*,"You have to run gyro_average_initialize before running any field solve routine."
       STOP
    END IF

#ifdef GDAGGER_G
    gyromatrix_dagger => get_gyromatrix_dagger()
    IF (.NOT.ASSOCIATED(gyromatrix_dagger)) THEN
       PRINT*,"You have to run gyro_average_initialize before running any field solve routine."
       STOP
    END IF
#endif

    CALL initialize(profilematrix,nCoeffs,nCoeffs,transposed)
    call allocate(profilematrix)

    !CALL initialize(profdotgyro,nCoeffs,nCoeffs)
    !call allocate(profdotgyro)
    CALL initialize(profdotgyro,nCoeffs,nCoeffs,transposed)
    call allocate(profdotgyro)

    !CALL initialize(single_gyromatrix,nCoeffs,nCoeffs)
    !call allocate(single_gyromatrix)

    CALL initialize(gyromatrix2,nCoeffs,nCoeffs,transposed)
    CALL allocate(gyromatrix2)
    !CALL initialize(full_gyromatrix2,nCoeffs,nCoeffs)
    !call allocate(full_gyromatrix2)

    !---------------------- EXECUTION ----------------------------!    
    do k=lk1,lk2
       do j=lj1,lj2
          do n=ln1,ln2
             CALL set_zero(mom00_flrcorr(j,k,n))
             CALL set_zero(mom01_flrcorr(j,k,n))
             IF (n_fields.gt.1) CALL set_zero(gamma0_num_matrix(j,k,n))
             !---------------------- new field solver matrices ----------------------------! 

#ifdef GDAGGER_G
             do m=lm1,lm2
                !CALL convert_Banded_to_Full(gyromatrix(j,k,m,n),single_gyromatrix)

                do index=1,4
                   !index 3 (numerical vpar integration) needs an extra loop and is therefore excluded here
                   if (index.ne.3) then
                      !---- Build diagonal profilematrix according to currently selected index ------!
                      CALL set_zero(profilematrix)
                      DO i=my_pex*li0+1,(my_pex+1)*li0
                         select case (index)
                         case(1) !prefactors * (Gamma0 - 1)
                            tmpvar = mu_weight(m)*geom%Bfield(i-1+pi1gl,pj1,k)*&
                                 &EXP(-mu(m)*geom%Bfield(i-1+pi1gl,pj1,k)/spec(n)%temp_prof(i-1+pi1gl))&
                                 &* spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl) &
                                 &/ (spec(n)%temp*(spec(n)%temp_prof(i-1+pi1gl))**2) 
                         case(2) !prefactors * (Gamma_mu - 1)
                            tmpvar = mu_weight(m)*mu(m)*(geom%Bfield(i-1+pi1gl,pj1,k))**2*&
                                 &EXP(-mu(m)*geom%Bfield(i-1+pi1gl,pj1,k)/spec(n)%temp_prof(i-1+pi1gl))&
                                 &* spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl)&
                                 &/ (spec(n)%temp*(spec(n)%temp_prof(i-1+pi1gl))**2)
                         case(4) !prefactors * (Gamma_vp2 - 1)
                            tmpvar = 0.5*mu_weight(m)*geom%Bfield(i-1+pi1gl,pj1,k)*&
                                 &EXP(-mu(m)*geom%Bfield(i-1+pi1gl,pj1,k)/spec(n)%temp_prof(i-1+pi1gl))&
                                 &* spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl) &
                                 &/ (spec(n)%temp*spec(n)%temp_prof(i-1+pi1gl)) 
                         end select

                         !\todo check whether the trap_pass modification should really
                         !! be applied to all moments (will affect the diagnostics)
                         IF(trap_pass .and. (spec(n)%charge .eq. -1)) tmpvar = tmpvar * &
                              &erf(sqrt((Bfieldmax(i-1+pi1gl)-geom%Bfield(i-1+pi1gl,pj1,k))* &
                              &mu(m)/(spec(n)%temp_prof(i-1+pi1gl))))

                         CALL set_value(profilematrix,i,i,tmpvar)
                         !diag_profilematrix(i) = tmpvar
                      ENDDO
                      call commit_values(profilematrix)
                      !--------------------- End profile matrix part --------------------------------!
                      !--------------------- Generate squared gyromatrix ----------------------------!
                      !dot-multiply gyromatrix(dagger)*profilematrix*gyromatrix
                      !CALL dot_multiply(profilematrix,single_gyromatrix,profdotgyro,'N')
                      CALL dot_multiply(profilematrix,gyromatrix(j,k,m,n),profdotgyro,'N')
                      !CALL dot_multiply(single_gyromatrix,profdotgyro,full_gyromatrix2,G_op)
                      !CALL dot_multiply(gyromatrix(j,k,m,n),profdotgyro,gyromatrix2,G_op)
                      CALL dot_multiply(gyromatrix_dagger(j,k,m,n),profdotgyro,gyromatrix2,'N')
                      !CALL convert_Full_to_Banded(full_gyromatrix2,gyromatrix2)

                      !-------------- Generate flr corrections matrices containing profiles ---------!
                      select case(index)
                      case (1)
                         DO i=my_pex*li0+1,(my_pex+1)*li0
                            CALL row_axpy(mom00_flrcorr(j,k,n),i,gyromatrix2,cmplx(1.0,0.0))
                         ENDDO
                         CALL commit_values(mom00_flrcorr(j,k,n))
                      case (2)
                         DO i=my_pex*li0+1,(my_pex+1)*li0
                            CALL row_axpy(mom01_flrcorr(j,k,n),i,gyromatrix2,cmplx(1.0,0.0))
                         ENDDO
                         CALL commit_values(mom01_flrcorr(j,k,n))
                      case (4)
                         DO i=my_pex*li0+1,(my_pex+1)*li0
                            CALL row_axpy(mom20_flrcorr(j,k,n),i,gyromatrix2,cmplx(1.0,0.0))
                         ENDDO
                         CALL commit_values(mom20_flrcorr(j,k,n))
                      end select
                   elseif (n_fields.gt.1) then
                      do l=ll1,ll2
                         !---------- Build profilematrix for Gamma0_num_matrix ----------------------!
                         CALL set_zero(profilematrix)
                         DO i=my_pex*li0+1,(my_pex+1)*li0
                            tmpvar = spec(n)%charge**2*spec(n)%dens/&
                                 &(spec(n)%mass*spec(n)%temp_prof(i-1+pi1gl))*&
                                 & fm(i-1+pi1gl,pj1,k,l,m,n)*&
                                 & vp(l)**2 * mat_00(i-1+pi1gl,pj1,k,l,m)
                            CALL set_value(profilematrix,i,i,tmpvar)
                         ENDDO
                         call commit_values(profilematrix)
                         !-------------------------- End profilematrix part -------------------------!


                         !----------------------- Generate squared gyromatrix -----------------------!
                         !CALL dot_multiply(profilematrix,single_gyromatrix,profdotgyro,'N')
                         !CALL dot_multiply(single_gyromatrix,profdotgyro,full_gyromatrix2,G_op)
                         !CALL convert_Full_to_Banded(full_gyromatrix2,gyromatrix2)
                         
                         CALL dot_multiply(profilematrix,gyromatrix(j,k,m,n),profdotgyro,'N')
                         !CALL dot_multiply(gyromatrix(j,k,m,n),profdotgyro,gyromatrix2,G_op)
                         CALL dot_multiply(gyromatrix_dagger(j,k,m,n),profdotgyro,gyromatrix2,'N')
                         
                         !-------------------------- Build gamma0_num_matrix ------------------------!
                         do i=my_pex*li0+1,(my_pex+1)*li0
                            CALL row_axpy(gamma0_num_matrix(j,k,n),i,gyromatrix2,cmplx(1.0,0.0))
                         enddo
                         CALL commit_values(gamma0_num_matrix(j,k,n))
                      enddo
                   endif
                enddo
             enddo
             !mu integration
             CALL my_sum(mom20_flrcorr(j,k,n), mpi_comm_w)

#else
             !---------------------- old field solver matrices ----------------------------! 
             DO m=lm1,lm2
                !CALL convert_Banded_to_Full(gyromatrix(j,k,m,n),single_gyromatrix)
                !CALL dot_multiply(single_gyromatrix,single_gyromatrix,full_gyromatrix2,G_op)
                !CALL convert_Full_to_Banded(full_gyromatrix2,gyromatrix2)
                CALL dot_multiply(gyromatrix(j,k,m,n),gyromatrix(j,k,m,n),gyromatrix2,G_op)

                DO i=my_pex*li0+1,(my_pex+1)*li0
                   tmpvar = mu_weight(m)*geom%Bfield(i-1+pi1gl,pj1,k)*&
                        &EXP(-mu(m)*geom%Bfield(i-1+pi1gl,pj1,k)/&
                        &spec(n)%temp_prof(i-1+pi1gl)) * & !new:
                        &spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl)/&
                        & (spec(n)%temp*spec(n)%temp_prof(i-1+pi1gl)**2)

                   !\todo check whether the trap_pass modification should really
                   !! be applied to all moments (will affect the diagnostics)
                   IF(trap_pass .and. (spec(n)%charge .eq. -1)) tmpvar = tmpvar * &
                        &erf(sqrt((Bfieldmax(i-1+pi1gl)-geom%Bfield(i-1+pi1gl,pj1,k))* &
                        &mu(m)/(spec(n)%temp_prof(i-1+pi1gl))))
                   
                   CALL row_axpy(mom00_flrcorr(j,k,n),i,gyromatrix2,tmpvar)
                   CALL row_axpy(mom01_flrcorr(j,k,n),i,gyromatrix2,&
                        & tmpvar*mu(m)*geom%Bfield(i-1+pi1gl,pj1,k))
                ENDDO
                
                PERFON('cmt_val')
                CALL commit_values(mom00_flrcorr(j,k,n))
                CALL commit_values(mom01_flrcorr(j,k,n))
                PERFOFF

                IF (n_fields.gt.1) THEN
                   !numerical vpar integration (to avoid Ampere cancellation problem; analytical result can
                   !be switched on below)
                   DO l=ll1,ll2
                      DO i=pi1,pi2
                         tmpvar = fm(i,pj1,k,l,m,n)&
                              & *vp(l)**2 * mat_00(i,pj1,k,l,m) *&
                              &spec(n)%charge**2*spec(n)%dens/&
                              &(spec(n)%mass*spec(n)%temp_prof(i))
                         CALL row_axpy(gamma0_num_matrix(j,k,n),i-pi1gl+1,gyromatrix2,&
                              & tmpvar )
                      END DO
                   END DO ! l
                   
                END IF !n_fields > 1
                
             END DO ! m
#endif
             CALL my_sum(mom00_flrcorr(j,k,n), mpi_comm_w)
             CALL my_sum(mom01_flrcorr(j,k,n), mpi_comm_w)
             
             !vparallel and mu integration
             IF (n_fields.GT.1) THEN
                !switch on for analytic result (& comment out sum_vw!):
                !       VIntCorr=0.5*mom00_flrcorr
                !CALL my_complex_sum_vw(VIntCorr,SIZE(VIntCorr))
                CALL my_sum(gamma0_num_matrix(j,k,n), mpi_comm_vw)
             ENDIF

             !add the profile matrix
             DO i=my_pex*li0+1,(my_pex+1)*li0
                !\todo check whether the trap_pass modification should really
                !! be applied to all moments (will affect the diagnostics)
                IF (trap_pass .and. (spec(n)%charge .eq. -1)) THEN
                   tmpfac = trap(i-1+pi1gl,k)
                ELSE
                   tmpfac = 1.0
                ENDIF

                CALL add_value(mom00_flrcorr(j,k,n),i,i,&
                     & -spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl)/&
                     & (spec(n)%temp*spec(n)%temp_prof(i-1+pi1gl)) * tmpfac)
                CALL add_value(mom01_flrcorr(j,k,n),i,i,&
                     & -spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl)/&
                     & (spec(n)%temp) * tmpfac)
#ifdef GDAGGER_G
                Call add_value(mom20_flrcorr(j,k,n),i,i,&
                     & -0.5*spec(n)%charge * spec(n)%dens_prof(i-1+pi1gl)/&
                     & (spec(n)%temp) * tmpfac)
#else
                !in this case, mom20_flrcorr is simply Tp(x)/2*mom00_flrcorr
                DO ibase=1,nCoeffs
                   CALL set_value(mom20_flrcorr(j,k,n),i,ibase,&
                        & 0.5*spec(n)%temp_prof(i-1) * &
                        & mat_get_value(mom00_flrcorr(j,k,n),i,ibase) * tmpfac)
                END DO
#endif
             END DO

             Call commit_values(mom00_flrcorr(j,k,n))
             Call commit_values(mom01_flrcorr(j,k,n))
             Call commit_values(mom20_flrcorr(j,k,n))

          END DO ! n   
       END DO !j
    END DO !k

    exist_gammas=.TRUE.    


    !-------------------------- FINALIZATION ------------------------------- 
    CALL finalize(profilematrix)
    CALL finalize(profdotgyro)

    !CALL finalize(single_gyromatrix)
    !CALL finalize(full_gyromatrix2)
    CALL finalize(gyromatrix2)

  END SUBROUTINE initialize_flr_corr_moms
     
  SUBROUTINE finalize_flr_corr_moms
    integer:: j,k,n

    DO n=ln1,ln2
       DO k=lk1,lk2
          DO j=lj1,lj2
             CALL finalize(mom00_flrcorr(j,k,n))
             CALL finalize(mom01_flrcorr(j,k,n))
             CALL finalize(mom20_flrcorr(j,k,n))
             IF (n_fields.gt.1) CALL finalize(gamma0_num_matrix(j,k,n))
          END DO
       END DO
    END DO
    DEALLOCATE(mom00_flrcorr, mom01_flrcorr, mom20_flrcorr)
    IF (n_fields.gt.1) DEALLOCATE(gamma0_num_matrix)

    exist_gammas = .FALSE.

  end SUBROUTINE finalize_flr_corr_moms

!!!**************************************************************************************!!!
!!!**************************************************************************************!!!
!!!**************************************************************************************!!!

  Real Function mem_est_flr_corrs(mem_req_in)
    Real :: mem_req_in
    IF (xy_local) THEN
       mem_est_flr_corrs = mem_est_flr_corr_ff(mem_req_in)
    ELSE
       mem_est_flr_corrs = mem_est_flr_corr_df(mem_req_in)
    ENDIF 
  End function mem_est_flr_corrs

  Subroutine initialize_FLR_corrs
    IF (xy_local) THEN
       call initialize_flr_corr_ff
    ELSE
       call initialize_flr_corr_df
    ENDIF
  end Subroutine initialize_FLR_corrs

  Subroutine finalize_FLR_corrs
    IF (xy_local) THEN
       call finalize_flr_corr_ff
    ELSE
       call finalize_flr_corr_df
    ENDIF
  end Subroutine Finalize_FLR_corrs
     


!!!******************************************************************!!!

  !> Calculates FLR correction terms used e.g. in NRG and MOM diagnostics
  SUBROUTINE initialize_flr_corr_df
    IMPLICIT NONE

    IF (mype.EQ.0) THEN
       write(*,"(A)") "Initializing the FLR Correction prefactors for global calculation."
    END IF

    ! initialize the two vectors needed in the real calculation of 
    ! the FLR correction terms.
    CALL initialize(tmpvec,ni0)
    CALL initialize(resvec,ni0)

  END SUBROUTINE initialize_flr_corr_df

  SUBROUTINE finalize_flr_corr_df
    CALL finalize(tmpvec)
    CALL finalize(resvec)
  END SUBROUTINE finalize_flr_corr_df
  

!!!**************************************************************************************!!!
!!!**************************************************************************************!!!
  !>Give an estimate of the memory requirements of FLR correction terms
  Real Function mem_est_flr_corr_ff(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !corr_mom_aa_ff
    mem_loc = 3.*lijk0*ln0*SIZE_OF_COMPLEX_MB
    IF (equil_par_curr) mem_loc = mem_loc + 3.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB
    !corr_mom_aa_ff_bpar
    IF (n_fields .GT. 2) THEN
      mem_loc = mem_loc + 3.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB
      IF (equil_par_curr) mem_loc = mem_loc + 3.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB
    END IF
    !corr_mom_bb_ff and _bpar for Bpar (each have 2)
    IF (n_fields .GT. 2) mem_loc = mem_loc + 4.0 * lijk0 * ln0 * SIZE_OF_COMPLEX_MB

    !in correct_for_FLR
    !loc_phi
    mem_loc = mem_loc + li0*lj0*lz0*SIZE_OF_COMPLEX_MB
    !corrections
    mem_loc = mem_loc + lijk0*ln0*3*SIZE_OF_COMPLEX_MB

    mem_est_flr_corr_ff = mem_loc+mem_req_in

  End Function mem_est_flr_corr_ff

!!!******************************************************************!!!

  !> Calculates (local) FLR correction terms used e.g. in NRG and MOM diagnostics
  Subroutine initialize_flr_corr_ff
    Implicit None
    
    Real :: b2, b2spec(0:n_spec-1), delta(0:n_spec-1)
    Real :: k_perp2
    Integer::   i, j, k, n

    IF (n_fields.GT.2) n_flr_corrs = 6
    !parallel current, thus B_{0\parallel}^\ast
    IF (equil_par_curr) THEN
      ALLOCATE(aux_j0_flrcorr(lk1:lk2,ln1:ln2))

      n_flr_corrs = n_flr_corrs + 3

      DO n = ln1, ln2
        DO k = lk1, lk2
          aux_j0_flrcorr(k,n) = beta * currdens_par(pi1) / &
            (geom%Bfield(pi1,pj1,k) * geom%Bfield(pi1,pi2,k) * spec(n)%charge) * &
            SQRT(0.5*spec(n)%temp*spec(n)%mass)
        END DO
      END DO
    END IF

    Allocate(corr_mom_00_ff(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
    Allocate(corr_mom_01_ff(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
    Allocate(corr_mom_20_ff(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
    IF (equil_par_curr) THEN
      ALLOCATE(corr_mom_10_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_mom_30_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_mom_11_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
    END IF
    
    IF (n_fields .GT. 2) THEN
      ALLOCATE(corr_mom_00_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_mom_20_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_mom_01_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))

      ALLOCATE(corr_momI_00_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_momI_20_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_momI_01_ff(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_momI_00_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_momI_20_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      ALLOCATE(corr_momI_01_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))

      IF (equil_par_curr) THEN
        ALLOCATE(corr_mom_10_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
        ALLOCATE(corr_mom_30_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
        ALLOCATE(corr_mom_11_ff_bpar(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2))
      END IF
    END IF



    Do k=lk1,lk2
       Do j = lj1,lj2
          Do i = li1,li2
             k_perp2 = geom%gii(pi1,pj1,k)*ki(i)**2+2.0*geom%gij(pi1,pj1,k)*&
                  & ki(i)*kj(j)+geom%gjj(pi1,pj1,k)*kj(j)**2
             ! switch off Bfield for benchmark with Dimits
             b2 = k_perp2/(geom%Bfield(pi1,pj1,k)**2)
             do n=0,n_spec-1
                ! b2spec is identical to b_j in script/thesis
                b2spec(n) = b2*spec(n)%mass*spec(n)%temp/(spec(n)%charge**2)
                delta(n) = gamma0(b2spec(n)) - gamma1(b2spec(n))
             enddo

             If ((ki(i).Ne.0.).Or.(kj(j).Ne.0.)) Then
                do n=ln1,ln2
                   ! prefactor for polarisation density, used in diag_nrg and diag_mom
                   ! (corrects 00-moment for FLR terms)
                   corr_mom_00_ff(i,j,k,n) = -spec(n)%charge/spec(n)%temp*(1.0-gamma0(b2spec(n)))*dens_co(pi1,k,n)
                   ! prefactor for correction of par. temperature, par velocity, 
                   ! par heat transport (corrects 20-moment)
                   corr_mom_20_ff(i,j,k,n) = 0.5*corr_mom_00_ff(i,j,k,n)
                   ! prefactor for correction of perp. temperature, 
                   corr_mom_01_ff(i,j,k,n)   = -spec(n)%charge/spec(n)%temp &
                        & *(1.0 - gamma0(b2spec(n)) + b2spec(n)*delta(n))*dens_co(pi1,k,n)
                   ! prefactors for correction of q_par & q_perp are zero
                enddo
             ELSE
                corr_mom_00_ff(i,j,k,:) = 0.0
                corr_mom_01_ff(i,j,k,:) = 0.0
                corr_mom_20_ff(i,j,k,:) = 0.0
             End If

             IF (n_fields .GT. 2) THEN
               IF ((ki(i) .NE. 0.0) .OR. (kj(j) .NE. 0.0)) THEN
                 DO n = ln1, ln2
                   corr_mom_00_ff_bpar(i,j,k,n) = delta(n) / geom%Bfield(pi1,pj1,k) * dens_co(pi1,k,n)
                   corr_mom_20_ff_bpar(i,j,k,n) = 0.5 * corr_mom_00_ff_bpar(i,j,k,n)
                   corr_mom_01_ff_bpar(i,j,k,n) = (delta(n) * (1.0 - 2.0 * b2spec(n)) + &
                        &gamma0(b2spec(n))) / geom%Bfield(pi1,pj1,k) * dens_co(pi1,k,n)

                   corr_momI_00_ff(i,j,k,n) = spec(n)%charge/spec(n)%temp * delta(n) * dens_co(pi1,k,n)
                   corr_momI_00_ff_bpar(i,j,k,n) = 2.0 / geom%Bfield(pi1,pj1,k) * delta(n) * dens_co(pi1,k,n)
                   corr_momI_20_ff(i,j,k,n) = 0.5 * corr_momI_00_ff(i,j,k,n)
                   corr_momI_20_ff_bpar(i,j,k,n) = 0.5 * corr_momI_00_ff_bpar(i,j,k,n)
                   corr_momI_01_ff(i,j,k,n) = spec(n)%charge/spec(n)%temp * &
                        (gamma0(b2spec(n))+(1-b2spec(n))*delta(n)) * dens_co(pi1,k,n)
                   corr_momI_01_ff_bpar(i,j,k,n) = 2.0/geom%Bfield(pi1,pj1,k) * &
                        (gamma0(b2spec(n))+2*(1-b2spec(n))*delta(n)) * dens_co(pi1,k,n)
                 END DO
               ELSE
                 corr_mom_00_ff_bpar(i,j,k,:) = 0.0
                 corr_mom_01_ff_bpar(i,j,k,:) = 0.0
                 corr_mom_20_ff_bpar(i,j,k,:) = 0.0

                 corr_momI_00_ff(i,j,k,:) = 0.0
                 corr_momI_00_ff_bpar(i,j,k,:) = 0.0
                 corr_momI_20_ff(i,j,k,:) = 0.0
                 corr_momI_20_ff_bpar(i,j,k,:) = 0.0
                 corr_momI_01_ff(i,j,k,:) = 0.0
                 corr_momI_01_ff_bpar(i,j,k,:) = 0.0
               END IF
             END IF

             IF (equil_par_curr) THEN !need corrections for Apar terms: 10,30,11
               IF ((ki(i) .NE. 0.0) .OR. (kj(j) .NE. 0.0)) THEN
                 DO n = ln1, ln2
                   corr_mom_10_ff(i,j,k,n) = 0.5 * aux_j0_flrcorr(k,n) * &
                     spec(n)%charge / spec(n)%temp * (gamma0(b2spec(n)) - 1) * dens_co(pi1,k,n)
                   corr_mom_30_ff(i,j,k,n) = 0.5 * aux_j0_flrcorr(k,n) * &
                     spec(n)%charge / spec(n)%temp * &
                     (gamma0(b2spec(n)) - 1 - b2spec(n) * delta(n))
                   corr_mom_11_ff(i,j,k,n) = 0.75 * aux_j0_flrcorr(k,n) * &
                     spec(n)%charge / spec(n)%temp * (gamma0(b2spec(n)) - 1) * dens_co(pi1,k,n)
                 END DO
               ELSE
                 corr_mom_10_ff(i,j,k,:) = 0.0
                 corr_mom_30_ff(i,j,k,:) = 0.0
                 corr_mom_11_ff(i,j,k,:) = 0.0
               END IF
               IF (n_fields .GT. 2) THEN
                 IF ((ki(i) .NE. 0.0) .OR. (kj(j) .NE. 0.0)) THEN
                   DO n = ln1, ln2
                     corr_mom_10_ff_bpar(i,j,k,n) = 0.5 * aux_j0_flrcorr(k,n) / geom%Bfield(pi1,pj1,k) * &
                          delta(n) * dens_co(pi1,k,n)
                     corr_mom_30_ff_bpar(i,j,k,n) = 0.5 * aux_j0_flrcorr(k,n) / geom%Bfield(pi1,pj1,k) * &
                          ((1 - b2spec(n)) * delta(n) + gamma0(b2spec(n))) * dens_co(pi1,k,n)
                     corr_mom_11_ff_bpar(i,j,k,n) = 0.75 * aux_j0_flrcorr(k,n) / geom%Bfield(pi1,pj1,k) * &
                          delta(n) * dens_co(pi1,k,n)
                   END DO
                 ELSE
                   corr_mom_10_ff_bpar(i,j,k,:) = 0.0
                   corr_mom_30_ff_bpar(i,j,k,:) = 0.0
                   corr_mom_11_ff_bpar(i,j,k,:) = 0.0
                 END IF
               END IF
             END IF
          End Do
       End Do
    End Do
  End Subroutine initialize_flr_corr_ff


  !> Deallocates FLR correction terms
  Subroutine finalize_flr_corr_ff
    Implicit None
    
    Deallocate(corr_mom_00_ff,corr_mom_01_ff,corr_mom_20_ff)
    IF (equil_par_curr) &
      DEALLOCATE(corr_mom_10_ff,corr_mom_30_ff,corr_mom_11_ff)
    IF (n_fields .GT. 2) THEN
       DEALLOCATE(corr_mom_00_ff_bpar,corr_mom_20_ff_bpar,corr_mom_01_ff_bpar)
       DEALLOCATE(corr_momI_00_ff,corr_momI_20_ff,corr_momI_01_ff,&
            corr_momI_00_ff_bpar,corr_momI_20_ff_bpar,corr_momI_01_ff_bpar)
       
       IF (equil_par_curr) THEN
          DEALLOCATE(aux_j0_flrcorr)
          DEALLOCATE(corr_mom_10_ff_bpar,corr_mom_30_ff_bpar,corr_mom_11_ff_bpar)
       ENDIF
    END IF

  End Subroutine finalize_flr_corr_ff


!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Compute FLR corrections as functions of phi
  SUBROUTINE correct_for_FLR(locphi, corrections)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz),TARGET :: locphi
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,1:n_flr_corrs),TARGET :: corrections
    INTEGER :: j,k,n
    !****

    IF (xy_local) THEN
       DO n=ln1,ln2
          Do k=lk1,lk2
             corrections(:,:,k,n,1) = locphi(li1:li2,:,k) * corr_mom_00_ff(:,:,k,n)
             ! Add FLR Correction terms to the temperatures and heat fluxes
             corrections(:,:,k,n,2) = locphi(li1:li2,:,k) * corr_mom_20_ff(:,:,k,n)
             corrections(:,:,k,n,3) = locphi(li1:li2,:,k) * corr_mom_01_ff(:,:,k,n)
          END DO
       END DO

       IF (equil_par_curr) THEN
         DO n = ln1, ln2
           DO k = lk1, lk2
             corrections(:,:,k,n,4) = locphi(li1:li2,:,k) * corr_mom_10_ff(:,:,k,n)
             corrections(:,:,k,n,5) = locphi(li1:li2,:,k) * corr_mom_30_ff(:,:,k,n)
             corrections(:,:,k,n,6) = locphi(li1:li2,:,k) * corr_mom_11_ff(:,:,k,n)
           END DO
         END DO
       END IF
    ELSE
       DO n=ln1,ln2
          DO k=lk1,lk2
             DO j=lj1,lj2
                CALL attach(tmpvec,locphi(li1:li2,j,k))

                CALL attach(resvec,corrections(li1:li2,j,k,n,1)) 
                CALL dot_multiply(mom00_flrcorr(j,k,n),tmpvec,resvec)

                CALL attach(resvec,corrections(li1:li2,j,k,n,2))
                CALL dot_multiply(mom20_flrcorr(j,k,n),tmpvec,resvec)

                CALL attach(resvec,corrections(li1:li2,j,k,n,3))
                CALL dot_multiply(mom01_flrcorr(j,k,n),tmpvec,resvec)

             END DO
          END DO
       END DO
    ENDIF
    
  END SUBROUTINE correct_for_FLR

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Compute FLR corrections as functions of phi and bpar
  SUBROUTINE correct_for_FLR_bpar(locphi,locbpar,corrections)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), TARGET :: locphi, locbpar
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,1:n_flr_corrs), TARGET :: corrections
    INTEGER :: k, n
    !****

    IF (xy_local) THEN
      DO n = ln1, ln2
        Do k = lk1, lk2
          corrections(:,:,k,n,1) = locphi(li1:li2,:,k) * corr_mom_00_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_mom_00_ff_bpar(:,:,k,n)
          corrections(:,:,k,n,2) = locphi(li1:li2,:,k) * corr_mom_20_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_mom_20_ff_bpar(:,:,k,n)
          corrections(:,:,k,n,3) = locphi(li1:li2,:,k) * corr_mom_01_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_mom_01_ff_bpar(:,:,k,n)
          !correction of the momI_00 moment (called N00 in the doc/observables.tex script)
          corrections(:,:,k,n,4) = locphi(li1:li2,:,k) * corr_momI_00_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_momI_00_ff_bpar(:,:,k,n)
          !correction of the momI_20 moment (called N20 in the doc/observables.tex script)
          corrections(:,:,k,n,5) = locphi(li1:li2,:,k) * corr_momI_20_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_momI_20_ff_bpar(:,:,k,n)
          !correction of the momI_01 moment (called N02 in the doc/observables.tex script)
          corrections(:,:,k,n,6) = locphi(li1:li2,:,k) * corr_momI_01_ff(:,:,k,n) + &
                                   locbpar(li1:li2,:,k) * corr_momI_01_ff_bpar(:,:,k,n)
          IF (equil_par_curr) THEN
            corrections(:,:,k,n,7) = locphi(li1:li2,:,k) * corr_mom_10_ff(:,:,k,n) + &
                                     locbpar(li1:li2,:,k) * corr_mom_10_ff_bpar(:,:,k,n)
            corrections(:,:,k,n,8) = locphi(li1:li2,:,k) * corr_mom_30_ff(:,:,k,n) + &
                                     locbpar(li1:li2,:,k) * corr_mom_30_ff_bpar(:,:,k,n)
            corrections(:,:,k,n,9) = locphi(li1:li2,:,k) * corr_mom_11_ff(:,:,k,n) + &
                                     locbpar(li1:li2,:,k) * corr_mom_11_ff_bpar(:,:,k,n)
          END IF
        END DO
      END DO
    ELSE
      IF (mype .LE. 0) WRITE(*,"(A)") "FLRcorr: no Bpar for global code"
      STOP
    ENDIF
    
  END SUBROUTINE correct_for_FLR_bpar



End Module flr_corr 
