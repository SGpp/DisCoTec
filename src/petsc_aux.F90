#include "redef.h"
#include "petscversion.h"
!>Routines for handling the KSP and PC objects used in PETSc (and SLEPc) 
!!and routines for the conversion between fortran and PETSc vector formats
!!
!!It provides a module variable ksp_obj which should be used by all routines based on PETSc 
!!and a variable shellmat which defines the operator for the matrix-free PETSc methods used 
!!in GENE 
module petsc_aux
  use par_other, only: dt, prec, print_ini_msg
  use par_in
  use eigen_parameters
  use discretization
  use calc_rhs
  use communications, only: MY_MPI_COMM_WORLD
  use petsc
  use petscvec
  use petscmat
  use petscksp
  use petscpc
  implicit none

  public:: pe2fo, fo2pe, initialize_petsc_mat, finalize_petsc_mat, initialize_petsc_vec, finalize_petsc_vec,&
       initialize_petsc_flags
  public:: glob_g, glob_rhs, v0, ksp_obj, shellmat, globerr,a_mat,b_mat, impl_shift
#ifdef COMBI
  public:: v0_out, petsc_shift,ksp_test_obj
#endif  
  public:: matop_mult, different_nonzero_pattern, same_nonzero_pattern, petsc_null_object,&
       petsc_null_character, petsc_null_scalar, petsc_decide, petsc_null_integer, petsc_comm_world,&
       kspreason,petsc_t,petsc_f
  public:: my_mpi_comm_world

  private
#if (PETSC_VERSION_MINOR<1)&&(PETSC_VERSION_MAJOR<4)
#include "finclude/petscdef.h"
#else
#include "finclude/petscsysdef.h"
#endif
#include "finclude/petscvecdef.h"   
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"

  Mat shellmat, L_g_mat,a_mat,b_mat
  Vec glob_g, glob_rhs, v0(5)
  KSP ksp_obj
  KSPConvergedReason kspreason
  PetscErrorCode globerr
#ifdef COMBI
  Vec v0_out(5)
  KSP ksp_test_obj
  PetscScalar petsc_shift
#endif

!petsc<3.2 does not know PetscBool but PetscTruth
#ifndef PetscBool
#define PetscBool PetscTruth
#endif

  PetscBool petsc_t, petsc_f

  logical:: impl_shift=.false.

contains
  !>Converts an array of type g_1 from PETSc to Fortran format
  subroutine pe2fo(petg,locg)
    Vec petg !<petsc representation of a g_1-type array
    !>fortran array of type g_1
    complex, dimension(lijklmn0),intent(out):: locg 
    PetscScalar,pointer:: arr_point(:)
    integer:: ierr
   
    PERFON('pefo')
    call VecGetArrayF90(petg,arr_point,ierr)   

    if ((evenx.eq.1).and.(xy_local)) then
       call insert_zero_kx(arr_point,locg)
    else
       locg=arr_point
    endif
    call VecRestoreArrayF90(petg,arr_point,ierr)
    
    PERFOFF
  end subroutine pe2fo

  !>Converts a vector of type g_1 from Fortran to PETSc format
  subroutine fo2pe(locg,petg) !convert g from fortran to petsc format
    !>fortran array of type g_1
    complex, dimension(lijklmn0),intent(in):: locg 
    Vec petg !<petsc representation of a g_1-type array
    PetscScalar,pointer:: arr_point(:)
    integer:: ierr

    PERFON('fope')
    call VecGetArrayF90(petg,arr_point,ierr)   

    if ((evenx.eq.1).and.(xy_local)) then
       call remove_zero_kx(locg,arr_point)
    else
       arr_point=locg
    endif
    call VecRestoreArrayF90(petg,arr_point,ierr)
    PERFOFF

  end subroutine fo2pe

  subroutine remove_zero_kx(locg,locsmallg)
    complex, dimension(lg1:lg2,ljklmn0), intent(in):: locg 
    complex, dimension(lg1:lg2-1,ljklmn0), intent(out):: locsmallg 
    integer:: ind

    do ind=1,ljklmn0
       locsmallg(lg1:hkx,ind)=locg(lg1:hkx,ind)
       locsmallg(lkx-1:lg2-1,ind)=locg(lkx:lg2,ind)
    enddo

  end subroutine remove_zero_kx

  subroutine insert_zero_kx(locsmallg, locg)
    complex, dimension(lg1:lg2-1,ljklmn0), intent(in):: locsmallg 
    complex, dimension(lg1:lg2,ljklmn0), intent(out):: locg 
    integer:: ind

    do ind=1,ljklmn0
       locg(lg1:hkx,ind)=locsmallg(lg1:hkx,ind)
       locg(hkx+1,ind)=0
       locg(lkx:lg2,ind)=locsmallg(lkx-1:lg2-1,ind)
    enddo

  end subroutine insert_zero_kx

  !>Initializes flags for calling petsc routines (only logicals right now..)
  subroutine initialize_petsc_flags
    petsc_t = PETSC_TRUE
    petsc_f = PETSC_FALSE
  end subroutine initialize_petsc_flags
  
  !>Initializes shellmat for time step computation (comp_timestep.F90) or the IE1p 
  !!time stepping scheme (impl_petsc.F90) and creates the necessary PETSc vectors
  subroutine initialize_petsc_mat(mat_invert,for_nlev)
    !>switch between initialization for matrix inversion and time step computation
    logical,intent(in):: mat_invert
    !>switch for diag_nlev : nonlinear eigenvalues
    logical,intent(in),optional:: for_nlev
    integer :: petsc_scalar=0, petsc_complex=0

#ifdef PETSC_SCALAR
    petsc_scalar = PETSC_SCALAR
    petsc_complex = PETSC_COMPLEX
#endif
    if ((mype==0).and.print_ini_msg) &
         & write(*,'(A,I1,A,I1,A,I1,A,I2.2)') "Using PETSc version ",PETSC_VERSION_MAJOR,&
         &'.',PETSC_VERSION_MINOR,'.',PETSC_VERSION_SUBMINOR,'-p',PETSC_VERSION_PATCH
    if (petsc_scalar.ne.petsc_complex) &
       & stop 'Sorry, you need to configure PETSc using -with-scalar-type=complex'

    if(present(for_nlev))then
       call initialize_CalFullRhs(for_petsc=.true.,for_nlev=for_nlev)
    else
       call initialize_CalFullRhs(for_petsc=.true.)
    endif

    call MatCreateShell(MY_MPI_COMM_WORLD,vlen,vlen,vlen*n_procs_sim,&
         vlen*n_procs_sim,PETSC_NULL_INTEGER,shellmat,globerr)
    call MatShellSetOperation(shellmat,MATOP_MULT,imp_matmul,globerr)

    call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,glob_g,globerr)

    impl_shift=mat_invert

  end subroutine initialize_petsc_mat

  !>Deletes shellmat and the PETSc vectors
  !!also finalizes CalFullRhs (unless the optional argument suppresses this)
  subroutine finalize_petsc_mat(final_rhs)
    implicit none
    !> switch for finalizing rhs computations, default is finalize 
    logical,intent(in),optional::final_rhs
    logical::finalize_rhs

    if (present(final_rhs)) then
       finalize_rhs = final_rhs
    else
       finalize_rhs = .true.
    endif

    call MatDestroy(shellmat,globerr)
    call VecDestroy(glob_g,globerr)
    if (finalize_rhs) call finalize_CalFullRhs

  end subroutine finalize_petsc_mat

  subroutine initialize_petsc_vec
    call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,glob_rhs,globerr)
  end subroutine initialize_petsc_vec

  subroutine finalize_petsc_vec
    call VecDestroy(glob_rhs,globerr)
  end subroutine finalize_petsc_vec

  !>Defines the matrix-free matrix vector multiplication by directly using CalFullRhs
  subroutine imp_matmul(locmat,vec_in,vec_out,ierr)
    
    Mat locmat !<PETSc matrix object which is defined in this routine
    Vec vec_in !<PETSc input vector
    Vec vec_out !<PETSc vector containing the result of locmat*vec_in
    PetscErrorCode ierr !<PETSC error code
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2):: loc_g,loc_res

    PERFON('rhsmult')

    call pe2fo(vec_in,loc_g)
    Call CalFullRhs(loc_g,loc_res,0)

    if (impl_shift) then
       !compute (1-dt*L)g
       loc_res=loc_g-dt*loc_res
    endif

    call fo2pe(loc_res,vec_out)
    ierr=0
    PERFOFF
  end subroutine imp_matmul

end module petsc_aux
