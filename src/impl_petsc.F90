#include "redef.h"
!>Derived rhs computation used for IE1p
!!
!!Computes one timestep for the implicit Euler scheme, i.e. \f$ g^{n+1}=(1-dt*L)^{-1}g^n\f$. 
!!The inversion is done iteratively with matrix-free PETSc, at the moment, no preconditioner is used
module impl_petsc
  use par_mod
  use petsc_aux
  use calc_rhs
  use petsc_precond
  use convergence_monitoring, only: omega_prec

  implicit none  
!#include "finclude/petsc.h"
  public:: initialize_CalFullRhs_iPETSc,CalFullRhs_iPETSc,finalize_CalFullRhs_iPETSc
  private

contains
  !>Initializes the necessary PETSc routines
  subroutine initialize_CalFullRhs_iPETSc

    PETSC_COMM_WORLD=MY_MPI_COMM_WORLD

    call PetscInitialize(PETSC_NULL_CHARACTER,globerr)
    call initialize_petsc_mat(.true.)
    call initialize_petsc_vec

    call KSPCreate(MY_MPI_COMM_WORLD,ksp_obj,globerr)
    call KSPSetTolerances(ksp_obj,omega_prec*0.1,1.e-35,1.e5,300000,globerr)

!    call KSPSetType(ksp_obj,KSPGMRES,globerr)
!    call KSPSetType(ksp_obj,KSPBCGS,globerr)
!    call KSPSetInitialGuessNonzero(ksp_obj,PETSC_TRUE,globerr)

    if(pc_type.eq.'none') then
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>4)  
       call KSPSetOperators(ksp_obj,shellmat,shellmat,globerr)
#else
       call KSPSetOperators(ksp_obj,shellmat,shellmat,DIFFERENT_NONZERO_PATTERN,globerr)
#endif
    else
       call initialize_L_g
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>4)  
       call KSPSetOperators(ksp_obj,shellmat,L_g_mat,globerr)
#else
       call KSPSetOperators(ksp_obj,shellmat,L_g_mat,DIFFERENT_NONZERO_PATTERN,globerr)
#endif
       call KSPGetPC(ksp_obj,pc_obj,globerr)
        !call KSPSetFromOptions(ksp_obj,globerr)
        !call PCSetFromOptions(pc_obj,globerr)
        !call PCSetUp(pc_obj,globerr)
        !call KSPSetUp(ksp_obj,globerr)
        !call finalize_L_g
       call set_pc_obj
    end if

  end subroutine initialize_CalFullRhs_iPETSc

  !>Computes loc_k=(1-dt*L)^(-1)*loc_g iteratively
  subroutine CalFullRhs_iPETSc(loc_g,loc_rhs)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in):: loc_g !<modified distribution function
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: loc_rhs !<result of the multiplication

    call fo2pe(loc_g,glob_g)
    call KSPSolve(ksp_obj,glob_g,glob_rhs,globerr)
    call pe2fo(glob_rhs,loc_rhs)

  end subroutine CalFullRhs_iPETSc     

  !>Deletes the PETSc objects
  subroutine finalize_CalFullRhs_iPETSc

    if(pc_type.ne.'none') call finalize_L_g
    call KSPDestroy(ksp_obj,globerr)
    call finalize_petsc_mat
    call finalize_petsc_vec

    call PetscFinalize(globerr)

  end subroutine finalize_CalFullRhs_iPETSc

end module impl_petsc
