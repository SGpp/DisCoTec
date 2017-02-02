#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the parallel nonlinearity
!!This module represents a first effort to consider this additional
!!O(epsilon) term.
!!\todo In a second step, it might be integrated in existing modules to 
!!avoid some overhead (e.g., computation of derivatives; only one 
!!backtransform combined with the ExB nonlinearity)
module parallel_nonlin
  use par_mod, only: imag, pi, ve_pnl_max, phase_1, phase_2, spec,&
       &par_sten_bound, par_sten, vderivative
  use par_other, only: equil_par_curr, n_fields, currdens_par, xy_local
  use par_in, only: fourier2d
  use coordinates, only: kx, ky, kymin, vp
  use communications, only: mpi_comm_y,&
       &my_real_max_to_all, my_barrier
  use discretization
  use blockindex
  use fourier
  use prefactors
  use gyro_average_ff_mod, only: jfac, I1_factor
  !USE x_derivatives, only: x_deriv_exc
  use geometry, only: geom, C_xy, rhostar, minor_r
  !USE nonlinearity
  use axpy
  use dchidz_term, only: z_deriv
  use mpi
  USE all_rhs_terms,only:this_nonlinear_term
  USE ff_nonlinear_term_mod,only: ff_nonlinear_term_t
  USE df_nonlinear_term_mod,only: df_nonlinear_term_t

  implicit none
  public:: initialize_parallel_nonlin, add_parallel_nonlin, finalize_parallel_nonlin,&
       &mem_est_parallel_nonlin, nl_to_fourier_xy

  private
  
  complex,dimension(:),allocatable:: zeros,val

  real,dimension(:,:,:,:),allocatable:: pdphidx,pdphidy
  complex,dimension(:,:,:,:),allocatable:: pdphidxy
  Real, Dimension(:), Allocatable:: pdphidz
  complex, dimension(:,:,:), allocatable:: dbarphidz

  integer:: lzero1, lzero2, init_status = 0
  COMPLEX,DIMENSION(:,:),POINTER :: ranshift_to,ranshift_back
  COMPLEX,DIMENSION(:),POINTER :: ranshift_to_df,ranshift_back_df

  REAL :: rho_Lref
contains
  
  Real function mem_est_parallel_nonlin(mem_req_in)
    real:: mem_req_in
    
    if (xy_local) then
       mem_est_parallel_nonlin = mem_est_parallel_nonlin_ff(mem_req_in)
    else
       mem_est_parallel_nonlin = mem_est_parallel_nonlin_df(mem_req_in)
    endif

    !pdphidz
    mem_est_parallel_nonlin = mem_est_parallel_nonlin + &
         & ln0*SIZE_OF_REAL_MB
  end function mem_est_parallel_nonlin

  subroutine initialize_parallel_nonlin
    Integer :: n

    if (init_status==1) return

    rho_lref = rhostar*minor_r

    if (equil_par_curr) then 
       stop 'equil_par_curr currently not implemented in parallel_nonlin'
    endif

    if (xy_local) then
       call initialize_parallel_nonlin_ff
    else
       call initialize_parallel_nonlin_df
    endif

    IF (turbdeal) THEN
       SELECT TYPE (this_nonlinear_term)
       CLASS is (ff_nonlinear_term_t)
          ranshift_to   => this_nonlinear_term%ranshift_to
          ranshift_back => this_nonlinear_term%ranshift_back
       CLASS is (df_nonlinear_term_t)
          ranshift_to_df   => this_nonlinear_term%prefactor%ranshift_to_df
          ranshift_back_df => this_nonlinear_term%prefactor%ranshift_back_df
       END SELECT
    END IF


    allocate(pdphidz(ln1:ln2))
    do n=ln1,ln2
       pdphidz(n) = -rho_lref * spec(n)%charge/spec(n)%mass/sqrt(2*spec(n)%temp/spec(n)%mass)
    enddo

    init_status = 1
  end subroutine initialize_parallel_nonlin


  subroutine add_parallel_nonlin(p_f,p_emfields,p_bar_emfields,p_dbarchidxy,p_rhs,lb1,lb2,stage)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    Complex, Dimension(:,:,:,:), pointer :: p_dbarchidxy
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout), target:: p_rhs    
    integer, intent(in):: stage

    !local variables
    real :: ve_pnl_max_sav
    logical:: first

    PERFON_I('add_pnl')

    if (stage.eq.1) then
       first=.true.
       ve_pnl_max_sav = ve_pnl_max
    else
       first=.false.
    end if

    if (xy_local) then
       ! add the ExB parallel_nonlin, calculation is performed in real space
       if(yx_order) then
          stop 'yx_order not implemented'
       else
          call calc_parallel_nonlin_ff(p_f,p_emfields,p_rhs,stage,lb1,lb2)
       endif
    else
       ! add the ExB parallel_nonlin, calculation is performed in real space
       call calc_parallel_nonlin_df(p_f,p_bar_emfields,p_dbarchidxy,p_rhs,lb1,lb2,first)
    endif

    if (first) then
       ! now in ve_pnl_max we find the maximal ExB velocity in the whole 
       ! volume. The reduction on the phase space has already been 
       ! performed, but we have to do it still over the species.
       ve_pnl_max=max(ve_pnl_max,ve_pnl_max_sav)
       !exchange only at the end of stage 1
       if(lb2.eq.lklmn0) then
          Call my_real_max_to_all(ve_pnl_max)
       endif
    endif
    
    PERFOFF_I

  end subroutine add_parallel_nonlin


  subroutine finalize_parallel_nonlin
    if (xy_local) then
       call finalize_parallel_nonlin_ff
    else
       call finalize_parallel_nonlin_df
    endif

    deallocate(pdphidz)

    init_status = 0
  end subroutine finalize_parallel_nonlin

!--------------------------------------------------------------------------------------

  Real function mem_est_parallel_nonlin_ff(mem_req_in)
    real:: mem_req_in

    !val, dphidxy
    mem_est_parallel_nonlin_ff = mem_req_in + &
         &li0da*(1+lj0*lk0*ll0*ln0)*SIZE_OF_COMPLEX_MB

    !in calc_parallel_nonlin_ff
    mem_est_parallel_nonlin_ff = mem_est_parallel_nonlin_ff + &
         & (2.*li0da*lj0+0.5*ly0da*li0da/n_procs_y+2.0*lij0)*SIZE_OF_COMPLEX_MB

    !in par_nl_to_direct_xy
    mem_est_parallel_nonlin_ff = mem_est_parallel_nonlin_ff + &
         &ly0da*li0da/n_procs_y*SIZE_OF_COMPLEX_MB

  end function mem_est_parallel_nonlin_ff
  
  subroutine initialize_parallel_nonlin_ff
    integer:: i,j,k,l

    lzero1=hkx+1
    lzero2=lkx+kx_offset-1
    allocate(zeros(lzero1:lzero2))
    zeros=(0.,0.)
    allocate(val(li1da:li2da))
    val=(0.,0.)

    allocate(pdphidxy(li1da:li2da,lj1:lj2,lk1:lk2,ll1:ll2))
    pdphidxy=0.

    if (mype==0) write(*,'(A,F13.6)') 'Initializing parallel nonlinearity with rho/Lref = ', rho_lref

    do l=ll1,ll2
       do k=lk1,lk2
          do j=lj1,lj2
             do i=li1,hkx
                pdphidxy(i,j,k,l) = -rho_lref * vp(l)/geom%Bfield(pi1,pj1,k) * &
                     &imag*(geom%K_i(pi1,pj1,k)*kx(i)+geom%K_j(pi1,pj1,k)*ky(j))       
             enddo
             do i=lkx,li2
                pdphidxy(i+kx_offset,j,k,l) = -rho_lref * vp(l)/geom%Bfield(pi1,pj1,k) * &
                     &imag*(geom%K_i(pi1,pj1,k)*kx(i)+geom%K_j(pi1,pj1,k)*ky(j))   
             enddo
          enddo
       enddo
    enddo

  end subroutine initialize_parallel_nonlin_ff


  subroutine finalize_parallel_nonlin_ff

    deallocate(zeros,val,pdphidxy)
    
  end subroutine finalize_parallel_nonlin_ff


  subroutine calc_parallel_nonlin_ff(p_f,p_emfields,p_rhs,stage,lb1,lb2)
    Integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    Complex,Dimension(li1:li2,lj1:lj2,lb1:lb2),intent(inout) :: p_rhs
    integer,intent(in):: stage

    complex, Dimension(li1da:li2da,lj1:lj2):: dphibar_terms, dfdv_ext
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: nonlin_real
    Complex,Dimension(li1:li2,lj1:lj2) :: nonlin, dfdv
    integer:: j,k,l,m,n,klmn,lb,ub,sten
    logical:: first

    PERFON_I('calcpnl')

    if (n_fields.gt.1) stop 'Parallel nonlinearity only available for electrostatic runs'

    do klmn=lb1,lb2

       k=sk(klmn)
       l=sl(klmn)
       m=sm(klmn)
       n=sn(klmn)
       
       !compute dphidx, dphidy, dphidz terms
       !\todo: skip l loop as barphi is independent of vpar?
       Do j=lj1,lj2

          val(li1:hkx)= p_emfields(li1:hkx,j,k,1)*jfac(li1:hkx,j,k,m,n)
          val(lkx+kx_offset:li2da)= p_emfields(lkx:li2,j,k,1)*jfac(lkx:li2,j,k,m,n)
          
          !dphidx/y terms
          dphibar_terms(:,j) = pdphidxy(:,j,k,l)*val

          !dphidz term
          do sten=-par_sten_bound,par_sten_bound
             val(li1:hkx)= p_emfields(li1:hkx,j,k+sten,1)*jfac(li1:hkx,j,k+sten,m,n)
             val(lkx+kx_offset:li2da)= p_emfields(lkx:li2,j,k+sten,1)*jfac(lkx:li2,j,k+sten,m,n)
             call axpy_ij(li0da,pdphidz(n)*par_sten(sten),val,dphibar_terms(:,j))
          enddo

       enddo
       
       lb=max(-l,-2)
       ub=min(nv0-1-l,2)
       sten=lb
       dfdv = vderivative(sten,pi1,k,l) * p_f(:,:,k,l+sten,m,n)
       do sten=lb+1,ub
          call axpy_ij(lij0,vderivative(sten,pi1,k,l),p_f(:,:,k,l+sten,m,n),dfdv)
       enddo

       do j=lj1,lj2
          dfdv_ext(li1:hkx,j) = dfdv(li1:hkx,j)
          dfdv_ext(lzero1:lzero2,j)=zeros
          dfdv_ext(lkx+kx_offset:li2da,j) = dfdv(lkx:li2,j) 
       end do
       
 !      PERFON_I('p_shift')
       if (turbdeal) then
          dphibar_terms = dphibar_terms * ranshift_to
          dfdv_ext = dfdv_ext * ranshift_to
       endif


 !      PERFON_I('rest')
       first=(stage.eq.1)
       Call par_nl_to_direct_xy(dphibar_terms,dfdv_ext,nonlin_real,first)
       Call nl_to_fourier_xy(nonlin_real, nonlin)


       if (turbdeal) then
          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + nonlin*ranshift_back
       else
          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + nonlin
       endif

    enddo
    PERFOFF_I

  end subroutine calc_parallel_nonlin_ff

  Subroutine par_nl_to_direct_xy(inarr1,inarr2,outarr,first)
    Complex,dimension(li1da:li2da, lj1:lj2), Intent(In):: inarr1, inarr2
    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1),Intent(inout):: outarr
    Logical,intent(in):: first

    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: rearr1,rearr2
    Real:: ve_pnl_max_loc

    ! Do inverse fourier transform to get real values
    if(fourier2d) then
       call nl_to_direct_xy_2d(inarr1,rearr1) ! dphiterms
       call nl_to_direct_xy_2d(inarr2,rearr2) ! dfdv
    else
       call nl_to_direct_xy_1d(inarr1,rearr1) ! dphiterms
       call nl_to_direct_xy_1d(inarr2,rearr2) ! dfdv
    end if

    ! Find maximum for ve_pnl
    if (first) then
       ve_pnl_max_loc=maxval(rearr1)
       ve_pnl_max=max(ve_pnl_max,ve_pnl_max_loc)
    endif

    !multiply the real arrays
    outarr=rearr1*rearr2

  end subroutine par_nl_to_direct_xy


!--------------------------------------------------------------------------------------
  subroutine initialize_parallel_nonlin_df
    integer :: l
    
    allocate(pdphidx(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2))
    allocate(pdphidy(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2))
    allocate(dbarphidz(li1:li2,lj1:lj2,lk1:lk2))

    do l=ll1,ll2
       pdphidx(:,:,:,l) = -rho_lref * vp(l)/geom%Bfield(pi1:pi2,pj1:pj1,lk1:lk2) * &
            &geom%K_i(pi1:pi2,pj1:pj1,lk1:lk2)
       pdphidy(:,:,:,l) = -rho_lref * vp(l)/geom%Bfield(pi1:pi2,pj1:pj1,lk1:lk2) * &
            &geom%K_j(pi1:pi2,pj1:pj1,lk1:lk2)
    enddo

  end subroutine initialize_parallel_nonlin_df


  Real function mem_est_parallel_nonlin_df(mem_req_in)
    real:: mem_loc, mem_req_in
    !local variables
    mem_loc = 3.*9.*nj0*li0da/n_procs_y*SIZE_OF_REAL_MB

    mem_loc = mem_loc + (lij0+2.*nj0*li0da/n_procs_y+&
         &(li0da/n_procs_y+2.*li0da+2*nib)*3.*nj0)*SIZE_OF_COMPLEX_MB

    mem_est_parallel_nonlin_df = mem_req_in + mem_loc

  end function mem_est_parallel_nonlin_df

  !this routine is only called if turbdeal=.t.
  subroutine finalize_parallel_nonlin_df 
    deallocate(pdphidx,pdphidy,dbarphidz)
  end subroutine finalize_parallel_nonlin_df


  !>Computes the parallel nonlinearity for the (x) nonlocal version
  Subroutine calc_parallel_nonlin_df(p_f,p_bar_emfields,p_dbarchidxy,p_rhs,lb1,lb2,first)
    integer,intent(in) :: lb1,lb2
    Complex,Dimension(li1:li2,lj1:lj2,lb1:lb2),intent(inout) :: p_rhs
    logical, intent(in):: first
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    Complex, Dimension(:,:,:,:), pointer :: p_dbarchidxy
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f

    !Local variables
    complex, Dimension(li1:li2,lj1:lj2) :: pnonlin, dphibar_terms, dfdv
    complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr
    real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: dphibar_terms_re, dfdv_re
    real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: pnonlin_re
    !Complex, Dimension(lbida2:ubida2,lj1:lj2) :: nl_tmp

    Real :: ve_pnl_max_loc
    Integer :: j,k,l,m,n,klmn,lb,ub,sten

    PERFON_I('calcpnl_df')

    do klmn=lb1,lb2

       k=sk(klmn)
       l=sl(klmn)
       m=sm(klmn)
       n=sn(klmn)

       !dphidx/y terms
       do j=lj1,lj2
          dphibar_terms(:,j) = pdphidx(:,pj1,k,l)*p_dbarchidxy(:,j,1,klmn-lb1+1)+&
               &pdphidy(:,pj1,k,l)*p_dbarchidxy(:,j,2,klmn-lb1+1)
       enddo

       !dphidz term
       !recompute dbardphidz only if the m or n index has changed
       if(compute_gy_av(klmn)) &
            &call z_deriv(p_bar_emfields(:,:,:,m,n,1),dbarphidz,l)
       dphibar_terms = dphibar_terms + pdphidz(n)*dbarphidz(:,:,k)

       !dfdv term
       lb=max(-l,-2)
       ub=min(nv0-1-l,2)
       do j=lj1,lj2
          sten=lb
          dfdv(:,j) = vderivative(sten,:,k,l) * p_f(:,j,k,l+sten,m,n)
          do sten=lb+1,ub
             dfdv(:,j) = dfdv(:,j)+vderivative(sten,:,k,l)*p_f(:,j,k,l+sten,m,n)
          enddo
       enddo

       if (turbdeal) then
          do j=lj1,lj2
             dphibar_terms(:,j)=dphibar_terms(:,j)*ranshift_to_df(j)
             dfdv(:,j)=dfdv(:,j)* ranshift_to_df(j)
          enddo
       endif

       ! anti-aliasing and y fourier transform
       PERFON_I('pnl_deal_FFT')
       Call gto_real(dphibar_terms,dphibar_terms_re,1)
       Call gto_real(dfdv,dfdv_re,1)
       PERFOFF_I

       ! get max ExB velocity
       if (first) then
          ve_pnl_max_loc=maxval(dphibar_terms_re)
          ve_pnl_max=max(ve_pnl_max,ve_pnl_max_loc)
       endif

       ve_pnl_max=maxval(dphibar_terms_re)

       PERFON_I('comp_pnl')
       pnonlin_re = dphibar_terms_re*dfdv_re

       ! transform back to fourier space
       Call to_fourier_y(pnonlin_re,tmp_arr)
       
       ! Transpose and remove zeros in y for dealiasing
       Call transpose_cmplx(nj0, li0da, tmp_arr, 0, pnonlin, 0)
       PERFOFF_I
       

       PERFON_I('plus')
       if (turbdeal) then
          do j=lj1,lj2
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + pnonlin(:,j)*ranshift_back_df(j)
          enddo
       else
          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + pnonlin
       endif
       PERFOFF_I

    enddo

  End Subroutine calc_parallel_nonlin_df

  !
  ! compute NL contributions in Fourier space using 1d fourier transforms.
  !
  Subroutine nl_to_direct_xy_1d(inarr,rearr)

    Complex,dimension(li1da:li2da, lj1:lj2), Intent(In):: inarr
    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1), Intent(InOut):: rearr

    Complex,dimension(0:ly0da/2, 0:li0da/n_procs_y-1):: temparr
    Complex:: rbuf(0:lj0-1,0:li0da/n_procs_y-1,0:n_procs_y-1)
    Complex:: tempypar(lj1:lj2, li1da:li2da)
    integer:: i, pe, ierr

    !Fourier transform array (first x direction, then transposition, then y)
    if (n_procs_y.gt.1) then
       Call fft_kx_to_x(inarr,tempypar)
#ifdef WITHOMP_BLOCKLOOP
       call mpi_alltoall(tempypar(lj1,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,&
            rbuf(0,0,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,threadlocal_mpi_comm_y,ierr)
#else
       call mpi_alltoall(tempypar(lj1,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,&
            rbuf(0,0,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,mpi_comm_y,ierr)
#endif
       do i=0,li0da/n_procs_y-1
          do pe=0,n_procs_y-1
             temparr(pe*lj0:(pe+1)*lj0-1,i)=rbuf(:,i,pe)
          enddo
       enddo
    else
       Call fft_kx_to_x(inarr,temparr)
    endif

    temparr(nky0:, :) = 0.

    Call fft_ky_to_y(temparr,rearr)
    
  End Subroutine nl_to_direct_xy_1d

  Subroutine nl_to_direct_xy_2d(inarr,rearr)

    Complex,dimension(li1da:li2da,lj1:lj2), Intent(In):: inarr
    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1), Intent(InOut):: rearr

    Real    :: realtmp(0:ly0da-1,0:li0da-1)
    Complex :: temparr(0:ly0da/2, 0:li0da-1)
    Complex :: rbuf(0:li0da-1,0:ly0da/2-1)
    integer :: ierr

    ! if n_procs_y > 1 get fourier components from other processes
    if (n_procs_y.gt.1) then
       call mpi_allgather(inarr,li0da*lj0,MPI_COMPLEX_TYPE, &
            rbuf,li0da*lj0,MPI_COMPLEX_TYPE,mpi_comm_y,ierr)       
       ! Transpose array
       temparr(0:ly0da/2-1,0:li0da-1) = transpose(rbuf)
    else       
       ! Transpose array if only 1 y proc
       temparr(0:ly0da/2-1,0:li0da-1) = transpose(inarr)
    endif

    ! For 3/2 dealiasing temparray holds 50% extra modes, 
    ! for turbo just one extra in the first dimension
    temparr(nky0:,:) = 0.0 ! Zero the extra fourier modes

    ! Use 2d inverse fourier transform
    if(n_procs_y.gt.1) then
       call fft_ff_to_xy(temparr,realtmp)
       ! Have to use the process id to chop up the output array
       rearr = realtmp(:,my_pey*li0da/n_procs_y:(my_pey+1)*li0da/n_procs_y-1)
    else
       call fft_ff_to_xy(temparr,rearr)
    endif
    
  End Subroutine nl_to_direct_xy_2d

  !
  ! Transform the nonlinearity back to Fourier space
  !
  Subroutine nl_to_fourier_xy(inarr, outarr)

    real, Intent(In) :: inarr(0:ly0da-1, 0:li0da/n_procs_y-1)
    complex, Intent(Out) :: outarr(li1:li2, lj1:lj2)

    if(fourier2d) then ! Choose 1D or 2D FFTs
       call nl_to_fourier_xy_2d(inarr, outarr)
    else
       call nl_to_fourier_xy_1d(inarr, outarr)
    end if

  End Subroutine nl_to_fourier_xy

  !
  ! Transform the nonlinearity back to Fourier space using 1d fourier transforms
  !
  Subroutine nl_to_fourier_xy_1d(inarr, outarr)

    real, Intent(In):: inarr(0:ly0da-1, 0:li0da/n_procs_y-1)
    complex, Intent(Out):: outarr(li1:li2, lj1:lj2)

    Complex:: temp_1(0:ly0da/2, 0:li0da/n_procs_y-1)
    Complex:: temp_3(li1da:li2da,lj1:lj2)
    Complex:: sendbuf(0:lj0-1,0:li0da/n_procs_y-1,0:n_procs_y-1)
    Complex:: rbuf(lj1:lj2,li1da:li2da)
    integer:: j,spo,pe,i,ierr

    spo=li0da-(nx0-1)/2

    !Fourier transform (first y direction, then transposition, then x)
    Call fft_y_to_ky(inarr,temp_1)

    if (n_procs_y.ne.1) then
       do i=0,li0da/n_procs_y-1
          do pe=0,n_procs_y-1
             sendbuf(:,i,pe)=temp_1(pe*lj0:(pe+1)*lj0-1,i)
          enddo
       enddo

#ifdef WITHOMP_BLOCKLOOP
       call mpi_alltoall(sendbuf(0,0,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,&
            rbuf(lj1,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,threadlocal_mpi_comm_y,ierr)
#else
       call mpi_alltoall(sendbuf(0,0,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,&
            rbuf(lj1,0),lj0*li0da/n_procs_y,MPI_COMPLEX_TYPE,mpi_comm_y,ierr)
#endif

       call fft_x_to_kx(rbuf,temp_3)
    else
       call fft_x_to_kx(temp_1,temp_3)
    endif

    !extract relevant modes
    if (evenx==1) then
       do j=lj1,lj2
          call ccopy((nx0-1)/2+1,temp_3(0,j),1,outarr(li1,j),1)
          outarr((nx0-1)/2+1,j)=0.
          call ccopy((nx0-1)/2,temp_3(spo,j),1,outarr(lkx,j),1)
       enddo
    else 
       do j=lj1,lj2
          call ccopy((nx0-1)/2+1,temp_3(0,j),1,outarr(li1,j),1)
          call ccopy((nx0-1)/2,temp_3(spo,j),1,outarr(lkx,j),1)
       enddo
    endif

  End Subroutine nl_to_fourier_xy_1d


 !
 ! Transform the nonlinearity back to Fourier space using 2d fourier transforms
 !
   Subroutine nl_to_fourier_xy_2d(inarr, outarr)

    real, Intent(In) :: inarr(0:ly0da-1, 0:li0da/n_procs_y-1)
    complex, Intent(Out) :: outarr(li1:li2, lj1:lj2)

    complex:: temp_1(0:ly0da/2, 0:li0da-1)  ! holds the transformed complex data before the transpose
    real   :: rbuf(0:ly0da-1, 0:li0da-1)     ! receives real data from other processes in the y direction
    integer:: i,j,spo,ierr, spodiff

    spo=li0da-(nx0-1)/2 
    spodiff=spo-nx0/2-1 ! The number of columns to miss

    ! Collect real data from other y comm processes
    if (n_procs_y.gt.1) then
       call mpi_allgather(inarr,ly0da*li0da/n_procs_y,MPI_REAL_TYPE, &
            rbuf,ly0da*li0da/n_procs_y,MPI_REAL_TYPE,mpi_comm_y,ierr)
       call fft_xy_to_ff(rbuf,temp_1)
    else
       call fft_xy_to_ff(inarr,temp_1)
    end if

    ! Transpose array back to y/2 format
    if(turbdeal) then
       outarr = transpose(temp_1(lj1:lj2,:))
    else
       do j=lj1,lj2                 ! For 3/2 scheme you ignore the central columns
          do i=0,nx0/2
             outarr(i,j) = temp_1(j,i)
          enddo
          do i=spo,li0da-1
             outarr(i-spodiff,j) = temp_1(j,i)
          enddo
       enddo
    end if

  End Subroutine nl_to_fourier_xy_2d

  !> Go to real space
  subroutine gto_real(inarr,outarr,howmany)
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
    real, dimension(:,:), pointer :: p_out
    Integer:: klmn, i_block,number_of_lbg0_blocks
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num

    my_thread = omp_get_thread_num()
#endif
    number_of_lbg0_blocks = howmany/lbg0
    if (number_of_lbg0_blocks.eq.0) then
       ! Transpose x-y
       Call transpose_cmplx(li0da, nj0, inarr(:,:,1), 0,tmp_arr2 , 0)
       ! Fourier transfrom in y (include dealiasing step)
       Call to_real_y(tmp_arr2,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks
          ! Transpose x-y
          do klmn=1,lbg0
             Call transpose_cmplx(li0da, nj0, inarr(:,:,(i_block-1)*lbg0+klmn), 0,tmp_arr2 , 0)

             ! Fourier transfrom in y (include dealiasing step)
             p_out => outarr(:,:,(i_block-1)*lbg0+klmn)
             Call to_real_y(tmp_arr2,p_out) !outarr(:,:,(i_block-1)*lbg0+klmn))
          end do
       end do
    end if

  end subroutine gto_real

  !----------------------------------------------------------------------
  !> Transpose the Matrix in_matrix and leave Result in transposed_matrix
  !! Use mpi_alltoall, pack data
  !----------------------------------------------------------------------
  Subroutine transpose_cmplx(n1, n2, in_matrix, ssoff, transposed_matrix, ddoff)

    Integer, Intent(in):: n1, n2, ssoff, ddoff
    Complex, Intent(in):: in_matrix(0:, 0:)
    Complex, Intent(out):: transposed_matrix(0:, 0:)

    Complex, Dimension(n2/n_procs_y, n1/n_procs_y, 0:n_procs_y-1) :: sbuf, rbufc
    Integer:: n1l, n2l, i1
    Integer:: pp, ierr

    n1l = n1/n_procs_y
    n2l = n2/n_procs_y

    PERFON_I("TransC")
    if (n_procs_y.eq.1) then
       transposed_matrix(ddoff:ddoff+n2-1,:)=Transpose(in_matrix(ssoff:ssoff+n1-1,:))
    else
       Do pp = 0, n_procs_y-1
          i1 = ssoff+pp*n1l
          sbuf(:,:,pp) = Transpose(in_matrix(i1:i1+n1l-1, :))
       Enddo
#ifdef WITHOMP_BLOCKLOOP
       Call mpi_alltoall(&
            sbuf, n1l*n2l, MPI_COMPLEX_TYPE,&
            rbufc, n1l*n2l, MPI_COMPLEX_TYPE,&
            threadlocal_mpi_comm_y, ierr)
#else
       Call mpi_alltoall(&
            sbuf, n1l*n2l, MPI_COMPLEX_TYPE,&
            rbufc, n1l*n2l, MPI_COMPLEX_TYPE,&
            mpi_comm_y, ierr)
#endif
       Do pp = 0, n_procs_y-1
          i1 = ddoff + pp*n2l
          transposed_matrix(i1:i1+n2l-1,:) = rbufc(:,:,pp)
       Enddo
    endif
    PERFOFF_I

  End Subroutine transpose_cmplx

end module parallel_nonlin
