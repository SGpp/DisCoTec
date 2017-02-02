#include "redef.h"

module rhs_term_mod
  implicit none
  
  type, abstract, public :: rhs_term_t
     integer :: init_status=0
     integer,public :: lbg0
     contains
       procedure (mem_est_if),deferred :: mem_est
       procedure (sub_empty_if), deferred :: initialize
       procedure (sub_empty_if), deferred :: finalize
       procedure (fun_empty_if), deferred :: getType
#ifdef GENERAL_STATE
       procedure (add_term_if), deferred :: add
#endif
       procedure :: set_blocksize => set_standard_blocksize
  end type rhs_term_t

  abstract interface
     function mem_est_if(this,mem_req_in) 
       import
       class(rhs_term_t) :: this
       real,intent(IN) :: mem_req_in
       real :: mem_est_if
     end function mem_est_if

     function fun_empty_if(this)
       import
       class(rhs_term_t) :: this
       character(len=MAX_TYPENAME_LENGTH) :: fun_empty_if
     end function fun_empty_if
     subroutine sub_empty_if(this)
       import
       class(rhs_term_t) :: this
     end subroutine sub_empty_if

#ifdef GENERAL_STATE
     subroutine add_term_if(this,p_in_state,p_out_state)
       import
       class(rhs_term_t),intent(INOUT) :: this
       class(state_t),intent(INOUT) :: p_in_state
       class(state_t),intent(INOUT) :: p_out_state
     end subroutine add_term_if
#endif
  end interface

  private
contains
  subroutine set_standard_blocksize(this,blocksize)
    class(rhs_term_t) :: this
    integer :: blocksize

    this%lbg0 = blocksize
  end subroutine set_standard_blocksize
end module rhs_term_mod

module advective_term_mod
  use rhs_term_mod
  implicit none

  type, abstract, public, extends(rhs_term_t) :: advective_term_t
  end type advective_term_t

end module advective_term_mod

module diffusive_term_mod
  use rhs_term_mod
  implicit none

  type, abstract, public, extends(rhs_term_t) :: diffusive_term_t
  end type diffusive_term_t

end module diffusive_term_mod


