#include "redef.h"
!>Wrapper for eigenvalue computations
module eigenvalue_comp

  use parameters_IO, only: write_parameters,final_init
  use eigen_parameters
  use phys_ini
  use aux_fields

#ifdef WITHSLEPC 
  use   eigen_iterative
#endif
#ifdef WITHSCAL
  use   eigen_direct
#endif
  implicit none
  public:: comp_eigenvalues
  private

contains

  !>Wrapper for the different eigenvalue solvers
  subroutine comp_eigenvalues
    if (which_ev.ne.'none') then
       
       final_init=.true.
       call initialize_current_parall
       call initialize_calc_aux_fields

       if (mype.eq.0) call write_parameters

       if (which_ev.eq.'all_mpl') then
#ifdef WITHSCAL
          call ev_direct
#else
          if (mype.eq.0) print*, 'all_mpl not possible without SCALAPACK! '
#endif
       else
#ifdef WITHSLEPC
#ifdef COMBI
        select case(which_ev)
        case('rhs_only')
                call ev_operator_only
        case('solve_system')
                call solve_for_residual
        case default
                call ev_iterative
        end select
#else
          call ev_iterative
#endif
#else
          if (mype.eq.0) print*, 'eigenvalue computation not possible without PETSc/SLEPc! '
#endif
       end if
       call finalize_calc_aux_fields
       call finalize_current_parall
    endif
  end subroutine comp_eigenvalues
  
end module eigenvalue_comp
