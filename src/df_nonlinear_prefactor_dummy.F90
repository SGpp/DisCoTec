#include "intrinsic_sizes.h"
module df_nonlinear_prefactor_dummy_mod
  use discretization, only: pi0, pi1,pi2, li1, li2,lj1,lj2
  use df_nonlinear_prefactor_mod
  implicit none

  type, public, extends(df_nonlinear_prefactor_t) :: df_nonlinear_prefactor_dummy_t
   contains
     procedure :: initialize => initialize_dummy_prefactor
     procedure :: finalize => finalize_dummy_prefactor
     procedure :: multiply_with => multiply_with_dummy_prefactor
     procedure :: mem_est => mem_est_dummy_prefactor
  end type df_nonlinear_prefactor_dummy_t

contains
  function mem_est_dummy_prefactor(this,mem_req_in) 
    class(df_nonlinear_prefactor_dummy_t) :: this
    real :: mem_req_in
    real :: mem_loc
    real :: mem_est_dummy_prefactor

    mem_loc = SIZE_OF_REAL_MB*pi0

    mem_est_dummy_prefactor = mem_req_in + mem_loc
  end function mem_est_dummy_prefactor

  subroutine initialize_dummy_prefactor(this)
    class(df_nonlinear_prefactor_dummy_t) :: this

    Allocate(this%pnl_1d(pi1:pi2))
    this%pnl_1d(pi1:pi2) = 1.0

  end subroutine initialize_dummy_prefactor

  subroutine finalize_dummy_prefactor(this)
    class(df_nonlinear_prefactor_dummy_t) :: this

    deallocate(this%pnl_1d)
  end subroutine finalize_dummy_prefactor

  subroutine multiply_with_dummy_prefactor(this,nonlin,localrhs,howmany)
    CLASS(df_nonlinear_prefactor_dummy_t),intent(IN) :: this
    INTEGER,intent(IN) :: howmany
    complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(IN) :: nonlin
    complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(INOUT) :: localrhs
    real :: testval
    integer :: klmn,j, omp_get_thread_num
#if 0
    !$OMP CRITICAL
    write(*,"(I2,2(A,ES20.10),A,I6)") omp_get_thread_num(),&
         &": in multiply_with_dummy_prefactor: nonlin = ",&
         &sum(abs(nonlin)),&
         &", localrhs = ",sum(abs(localrhs)),", howmany = ", howmany
    !write(*,"(I2,2ES20.10)") omp_get_thread_num(),sum(abs(localrhs(:,:,1))),&
    !     &sum(abs(localrhs(:,:,howmany)))
    do klmn=1,howmany
       testval = sum(abs(localrhs(:,:,klmn)))
       if (testval.gt.1e-5) then
          write(*,"(I2,I6,ES20.10)") omp_get_thread_num(),klmn,testval
       end if
    end do
    !$OMP END CRITICAL
#endif

    !$OMP BARRIER
    if (yx_order) then
       !$OMP DO
       do klmn=1,howmany
          do j=lj1,lj2
             localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) &
                   & -this%pnl_1d(li1:li2)*nonlin(li1:li2,j,klmn)
          end do
       end do
       !$OMP END DO
    else
       !$OMP DO
       do klmn=1,howmany
          do j=lj1,lj2
             localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) &
                  & + this%pnl_1d(li1:li2)*nonlin(li1:li2,j,klmn)
          end do
       end do
       !$OMP END DO
    endif
  end subroutine multiply_with_dummy_prefactor

end module df_nonlinear_prefactor_dummy_mod
