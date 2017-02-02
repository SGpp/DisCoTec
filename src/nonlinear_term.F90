module nonlinear_term_mod
  use discretization, only: li1,li2,lj1,lj2,lbz,ubz,lbi,ubi, ly0da, li0da
  use par_other,only: n_fields
  use advective_term_mod
  implicit none

  type, abstract, public, extends(advective_term_t) :: nonlinear_term_t
   contains
     procedure (add_nonlin_if),deferred,pass :: add
     procedure (construct_if),deferred :: construct
     procedure (destruct_if),deferred :: destruct
     procedure :: allocate_arrays => allocate_arrays_standard
     procedure :: free_arrays => free_arrays_standard
  end type nonlinear_term_t

  abstract interface
     subroutine add_nonlin_if(this,g_block,p_dgdxy,p_emfields,p_barchi,p_dbarchidxy,p_rhs,lb1,lb2,stage)
       import 
       class(nonlinear_term_t) :: this
       integer,intent(in) :: lb1,lb2
       complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout), target:: p_rhs    
       integer, intent(in):: stage
       complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
       Complex, Dimension(:,:,:), pointer :: p_barchi
       Complex, Dimension(:,:,:,:), pointer :: p_dbarchidxy, p_dgdxy
       complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in), target:: g_block
     end subroutine add_nonlin_if

     subroutine construct_if(this,equil_par_curr)
       import
       class(nonlinear_term_t) :: this
       logical :: equil_par_curr
     end subroutine construct_if

     ! if supported by all compilers, this will be a final subroutine
     subroutine destruct_if(this)
       import
       class(nonlinear_term_t) :: this
     end subroutine destruct_if

  end interface

  private

contains
  subroutine allocate_arrays_standard(this,vexy,dgdxy)
    class(nonlinear_term_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy,dgdxy

    allocate(vexy(li1:li2,lj1:lj2,2,this%lbg0))
    allocate(dgdxy(li1:li2,lj1:lj2,2,this%lbg0))
  end subroutine allocate_arrays_standard

  subroutine free_arrays_standard(this,vexy,dgdxy)
    class(nonlinear_term_t) :: this
    complex, dimension(:,:,:,:),allocatable :: vexy, dgdxy

    deallocate(vexy)
    deallocate(dgdxy)
  end subroutine free_arrays_standard
end module nonlinear_term_mod
