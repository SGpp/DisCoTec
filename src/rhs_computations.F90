#include "redef.h"

!>Wrapper for all available methods for the computation of the right hand side of the 
!!gyrokinetic equation
!!
!!Can use the primary implementations in calc_rhs, but also derived methods (for linear computations) 
!!like matrix-free PETSc, or methods based on an explicit representation of the linear operator
module rhs_computations
  use par_mod
  use discretization
  use communications  
  use calc_rhs
  use full_rhs_mat
#ifdef WITHSCAL
  use impl_scalapack
#endif
#ifdef WITHSLEPC
  use impl_petsc
#endif
  implicit none
  public:: initialize_calc_k, calc_k, finalize_calc_k, calc_k_coll
  ! routines for adaptive grids
  public:: calc_k_adptv
  private
  
  Integer :: init_status = 0

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_rhs_comput(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_loc=0
    mem_loc=mem_est_calc_rhs(mem_loc)

    mem_est_rhs_comput=mem_req_in+mem_loc
  End Function mem_est_rhs_comput 


  !>Calls the initialization routines of the various methods
  !!\todo implement init_status to _fullmat, _LU, _iPETSc
  subroutine initialize_calc_k

    if (direct_rhs) then
       call initialize_CalFullRhs
    else
       select case(derived_method)
       case(1) 
          call initialize_CalFullRhs_fullmat 
#ifdef WITHSCAL
       case(2) 
          call initialize_CalFullRhs_LU
#endif
#ifdef WITHSLEPC
       case(3)
          call initialize_CalFullRhs_iPETSc
#endif
       case default
          If(mype==0) Write(*,"(A)") "Necessary libraries (SCALAPACK/PETSc) are not linked!"
          Stop
       end select
    endif
    
    init_status = 1

  end subroutine initialize_calc_k

  !Computes the right hand side with the chosen method
  subroutine calc_k(p_g_, rhs, stage)
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(OUT) :: rhs
    Integer:: stage

    if (direct_rhs) then
       call CalFullRhs (p_g_,rhs,stage)
    else 
       select case(derived_method)
       case(1) 
          call CalFullRhs_fullmat(p_g_,rhs)
#ifdef WITHSCAL
       case(2) 
          call CalFullRhs_LU(p_g_,rhs)
#endif
#ifdef WITHSLEPC
       case(3) 
          call CalFullRhs_iPETSc(p_g_,rhs)
#endif
       end select
    endif

  end subroutine calc_k

  !Computes the right hand side with the chosen method
  subroutine calc_k_coll(p_g_, rhs, stage)
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(OUT) :: rhs
    Integer:: stage
    
    rhs_only_coll = .true.
    call CalFullRhs (p_g_,rhs,stage)
    rhs_only_coll = .false.
    !derived methods (IE1 for example) are not implemented at the moment

  end subroutine calc_k_coll
  
  !>Calls the cleanup routine for the chosen method 
  subroutine finalize_calc_k
    if (direct_rhs) then
       call finalize_CalFullRhs
    else
       select case(derived_method)
       case(1) 
          call finalize_CalFullRhs_fullmat
#ifdef WITHSCAL
       case(2) 
          call finalize_CalFullRhs_LU
#endif
#ifdef WITHSLEPC
       case(3) 
          call finalize_CalFullRhs_iPETSc
#endif
       end select
    endif

    init_status = 0
  end subroutine finalize_calc_k

  ! --- modified routines for adaptive grids
  
  !Computes the right hand side with the chosen method
  subroutine calc_k_adptv(p_g_, rhs, stage)
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(INOUT):: p_g_
    
    ! rhs is the right hand side, which is to be calculated from
    ! the given quantities
    Complex, Dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),Intent(OUT) :: rhs
    Integer:: stage

    if (direct_rhs) then
       call CalFullRhs_adptv (p_g_,rhs,stage)
    else
       stop "only direct_rhs is supporterd for adpative grids. Others: _fullmat, _LU, _iPETSc"
    endif
       
  end subroutine calc_k_adptv

end module rhs_computations
