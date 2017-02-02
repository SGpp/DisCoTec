#include "intrinsic_sizes.h"
#include "redef.h"
module df_epc_nl_prefactor_mod
  use df_nonlinear_prefactor_mod
  use par_in, only: beta, turbdeal
  use blockindex, only: sk,sl,sn
  use par_other, only: currdens_par
  use geometry, only: geom, C_xy
  use par_mod, only: spec,imag,pi
  use coordinates, only: vp
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
  implicit none

  type, public,extends(df_nonlinear_prefactor_t) :: df_epc_nl_prefactor_t
     Real, Dimension(:,:,:,:,:), Allocatable :: pnl !<nonlinearity prefactor

   contains
     procedure :: initialize => initialize_epc_nl
     procedure :: finalize => finalize_epc_nl
     procedure :: multiply_with => multiply_with_epc_nl
     procedure :: multiply_max => multiply_max_epc_nl
     procedure :: mem_est => mem_est_epc_nl
  end type df_epc_nl_prefactor_t

contains
  function mem_est_epc_nl(this,mem_req_in) 
    class(df_epc_nl_prefactor_t) :: this
    real :: mem_req_in
    real :: mem_loc
    real :: mem_est_epc_nl

    mem_loc = SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*ln0

    if (turbdeal) then
       ! Two ranshift arrays
       mem_loc = mem_loc + 2*SIZE_OF_COMPLEX_MB*lj0
       ! shiftmat
       mem_loc = mem_loc + 2*SIZE_OF_COMPLEX_MB*lj0
    end if

    mem_est_epc_nl = mem_req_in + mem_loc
  end function mem_est_epc_nl

  subroutine initialize_epc_nl(this)
    class(df_epc_nl_prefactor_t) :: this

    integer :: i,j,k,l,n
    REAL :: B_Bstar
    complex:: prefay


    Allocate(this%pnl(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,ln1:ln2))
    do n=ln1,ln2
       do l=ll1,ll2
          do k=lk1,lk2
             do j=pj1,pj2
                do i=pi1,pi2
                   B_Bstar = 1.0D0/(1.0D0+beta*&
                        &sqrt(0.5*spec(n)%mass*spec(n)%temp)*vp(l)*&
                        &currdens_par(i)/(spec(n)%charge*geom%Bfield(i,j,k)**2))
                   this%pnl(i,j,k,l,n)  = B_Bstar/C_xy(i)
                enddo
             enddo
          enddo
       enddo
    enddo

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

   end subroutine initialize_epc_nl

  subroutine finalize_epc_nl(this)
    class(df_epc_nl_prefactor_t) :: this

    if (turbdeal) then
       deallocate(this%shiftmat_df)
       deallocate(this%ranshift_to_df, this%ranshift_back_df)
    endif

    deallocate(this%pnl)
  end subroutine finalize_epc_nl

  !multilpy prefactor and add to rhs
  !
  subroutine multiply_with_epc_nl(this,nonlin,localrhs,howmany,lb1,lb2)
    CLASS(df_epc_nl_prefactor_t),INTENT(IN) :: this
    INTEGER,intent(IN) :: howmany, lb1, lb2
    !complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(IN) :: nonlin
    !complex, dimension(li1:li2,lj1:lj2,1:howmany),intent(INOUT) :: localrhs
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2),intent(IN) :: nonlin
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2),intent(INOUT) :: localrhs

    integer :: klmn,j

    if (yx_order) then
       if (turbdeal) then
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) - &
                     this%pnl(li1:li2,pj1,sk(klmn),sl(klmn),sn(klmn))*nonlin(li1:li2,j,klmn)*this%ranshift_back_df(j)
             end do
          end do
       else
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) - &
                     this%pnl(li1:li2,pj1,sk(klmn),sl(klmn),sn(klmn))*nonlin(li1:li2,j,klmn)
             end do
          end do
       endif
    else
       if (turbdeal) then
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) + &
                     this%pnl(li1:li2,pj1,sk(klmn),sl(klmn),sn(klmn))*nonlin(li1:li2,j,klmn)*this%ranshift_back_df(j)
             end do
          end do
       else
          !do klmn=1,howmany
          do klmn=lb1,lb2
             do j=lj1,lj2
                localrhs(li1:li2,j,klmn) = localrhs(li1:li2,j,klmn) + &
                     this%pnl(li1:li2,pj1,sk(klmn),sl(klmn),sn(klmn))*nonlin(li1:li2,j,klmn)
             end do
          end do
       endif
    endif
  end subroutine multiply_with_epc_nl

  !>multiplies vexy and nonlinear prefactor in real space for the nonlocal version
  !!and computes the max. abs. value used for the timestep estimate
  subroutine multiply_max_epc_nl(this,vexy_max,vexy_re,lb1,lb2)
    class(df_epc_nl_prefactor_t),intent(in) :: this
    integer,intent(IN) :: lb1,lb2
    real, dimension(0:ly0da-1,0:li0/n_procs_y-1,1:2,lb1:lb2),intent(IN) :: vexy_re
    real, dimension(1:2),intent(out) :: vexy_max
    real:: tmp  !temporary ve
    !real, dimension(0:li0/n_procs_y-1,1:2) :: vexy_max_tmp
    integer :: klmn, j, i, li

    li = pi1+my_pey*(pi0/n_procs_y) !lower bound of i index in prefactor

    LIKWID_ON('multmax')
    PERFON('multmax')
    !the real-space direction is distributed among y procs
    !it is located at the second index in vexy_re
    !klmn blocks locally run 1:howmany, but pnl is defined with indices k l n
    vexy_max = 0.0
    do klmn=lb1,lb2  !need to consider start of vexy_re
       do i=0,li0/n_procs_y-1
          do j=0,ly0da-1
             !ve_x  (vexy_re index 2)
             tmp=abs(vexy_re(j,i,2,klmn-lb1)*this%pnl(i+li,pj1,sk(klmn),sl(klmn),sn(klmn)))
             if (tmp>vexy_max(1)) vexy_max(1)=tmp
             !ve_y  (vexy_re index 1)
             tmp=abs(vexy_re(j,i,1,klmn-lb1)*this%pnl(i+li,pj1,sk(klmn),sl(klmn),sn(klmn)))
             if (tmp>vexy_max(2)) vexy_max(2)=tmp
          enddo
       enddo
    enddo
    !vexy_max(1) = maxval(vexy_max_tmp(:,2))  !ve_x
    !vexy_max(2) = maxval(vexy_max_tmp(:,1))  !ve_y
    PERFOFF
    LIKWID_OFF('multmax')

  end subroutine multiply_max_epc_nl

end module df_epc_nl_prefactor_mod
