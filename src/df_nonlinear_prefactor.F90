#include "intrinsic_sizes.h"
#include "redef.h"
module df_nonlinear_prefactor_mod
  use par_mod, only: imag,pi
  use discretization, only: pi0,pi1,pi2, li0,li1,li2, pj0,pj1,pj2, lj0,lj1,lj2, lk0,lk1,lk2, &
       & ll0,ll1,ll2, ln0,ln1,ln2, nj0, yx_order, my_pey, n_procs_y, ly0da, mype
  use par_in, only: turbdeal
  use geometry, only: C_xy
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif

  implicit none

  type, public :: df_nonlinear_prefactor_t
     Real, Dimension(:), Allocatable,public :: pnl_1d !<nonlinearity prefactor if B_{0\}_^*=B_0
     !complex, Dimension(:), Allocatable,public :: pnl_1d !<nonlinearity prefactor if B_{0\}_^*=B_0
     ! for turbdeal=.true. the following two arrays are needed
     complex,dimension(:,:), allocatable :: shiftmat_df
     COMPLEX,DIMENSION(:),pointer:: ranshift_to_df, ranshift_back_df

   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: multiply_with
     procedure :: multiply_max
     procedure :: mem_est
  end type df_nonlinear_prefactor_t

  interface
     subroutine zavdxpy(n, a,x,y) bind(c)
       use, intrinsic :: iso_c_binding
       integer(C_INT),value,intent(IN) :: n
       real(C_DOUBLE),dimension(*),intent(IN) :: a
       complex(C_DOUBLE_COMPLEX), dimension(*),intent(IN) :: x
       complex(C_DOUBLE_COMPLEX), dimension(*),intent(INOUT) :: y
     end subroutine zavdxpy
  end interface

contains
  function mem_est(this,mem_req_in) 
    class(df_nonlinear_prefactor_t) :: this
    real :: mem_req_in
    real :: mem_loc
    real :: mem_est

    mem_loc = SIZE_OF_REAL_MB*pi0

    if (turbdeal) then
       ! Two ranshift arrays
       mem_loc = mem_loc + 2*SIZE_OF_COMPLEX_MB*lj0
       ! shiftmat
       mem_loc = mem_loc + 2*SIZE_OF_COMPLEX_MB*lj0
    end if

    mem_est = mem_req_in + mem_loc
  end function mem_est

  subroutine initialize(this)
    class(df_nonlinear_prefactor_t) :: this

    integer :: j
    complex :: prefay

    Allocate(this%pnl_1d(pi1:pi2))
    this%pnl_1d(pi1:pi2) = 1.0/C_xy(pi1:pi2)

    if (turbdeal) then
       !turbo dealiasing
       allocate(this%ranshift_to_df(lj1:lj2), this%ranshift_back_df(lj1:lj2))
       allocate(this%shiftmat_df(lj1:lj2,0:1))
       prefay=imag*pi/2/(2*nj0)
       do j=lj1,lj2
          this%shiftmat_df(j,0)=exp(+j*prefay)
          this%shiftmat_df(j,1)=exp(-j*prefay)
       enddo
     endif

  end subroutine initialize

  subroutine finalize(this)
    class(df_nonlinear_prefactor_t) :: this

    if (turbdeal) then
       deallocate(this%shiftmat_df)
       deallocate(this%ranshift_to_df, this%ranshift_back_df)
    endif

    deallocate(this%pnl_1d)
  end subroutine finalize

  subroutine multiply_with(this,nonlin,localrhs,howmany,lb1,lb2)
    class(df_nonlinear_prefactor_t),intent(IN) :: this
    integer,intent(IN) :: howmany, lb1, lb2
    !complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(IN) :: nonlin
    !complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(INOUT) :: localrhs
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2),intent(IN) :: nonlin
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2),intent(INOUT) :: localrhs

    integer :: klmn,j,i

    LIKWID_ON('mult_pre')
    PERFON('mult_pre')
    if (yx_order) then
       if (turbdeal) then
          !$OMP DO
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) - &
                     this%pnl_1d(li1:li2)*nonlin(li1:li2,j,klmn)*this%ranshift_back_df(j)
             end do
          end do
          !$OMP END DO
       else
          !$OMP DO
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) - &
                     this%pnl_1d(li1:li2)*nonlin(li1:li2,j,klmn)
             end do
          end do
          !$OMP END DO
       endif
    else
       if (turbdeal) then
          !$OMP DO
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) + &
                     this%pnl_1d(li1:li2)*nonlin(li1:li2,j,klmn)*this%ranshift_back_df(j)
             end do
          end do
          !$OMP END DO
       else
          !$OMP DO collapse(2)
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
#if 0
                call zavdxpy(li0,this%pnl_1d,nonlin(:,j,klmn),localrhs(:,j,klmn))
#else
                !DIR$ simd VECTORLENGTH(2) 
                do i=li1,li2
                   localrhs(i,j,klmn) = localrhs(i,j,klmn) + &
                        this%pnl_1d(i)*nonlin(i,j,klmn)
                end do
#endif
             end do
          end do
          !$OMP END DO
       endif
    endif
    PERFOFF
    LIKWID_OFF('mult_pre')
  end subroutine multiply_with

  !>multiplies vexy and nonlinear prefactor in real space for the nonlocal version
  !!and computes the max. abs. value used for the timestep estimate
  subroutine multiply_max(this,vexy_max,vexy_re,lb1,lb2)
    class(df_nonlinear_prefactor_t),intent(IN) :: this
    integer,intent(IN) :: lb1, lb2
    real, dimension(0:ly0da-1,0:li0/n_procs_y-1,1:2,lb1:lb2),intent(IN) :: vexy_re
    real, dimension(1:2),intent(out) :: vexy_max
    real, dimension(0:li0/n_procs_y-1,1:2) :: vexy_max_tmp
    integer :: i, li

    li = pi1+my_pey*(pi0/n_procs_y) !lower index for relevant pnl_1d part

    LIKWID_ON('multmax')
    PERFON('multmax')
    !the real-space direction is distributed among y procs
    do i=0,li0/n_procs_y-1
       vexy_max_tmp(i,1)=maxval(abs(vexy_re(:,i,1,lb1:lb2)))*abs(this%pnl_1d(i+li))
       vexy_max_tmp(i,2)=maxval(abs(vexy_re(:,i,2,lb1:lb2)))*abs(this%pnl_1d(i+li))
    enddo
    vexy_max(1) = maxval(vexy_max_tmp(:,2))  !ve_x
    vexy_max(2) = maxval(vexy_max_tmp(:,1))  !ve_y
    PERFOFF
    LIKWID_OFF('multmax')

  end subroutine multiply_max

end module df_nonlinear_prefactor_mod
