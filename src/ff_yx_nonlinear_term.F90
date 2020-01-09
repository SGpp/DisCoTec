#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

module ff_yx_nonlinear_term_mod
  use par_mod, only: imag, pi, ve_max, phase_1, phase_2, spec
  use par_other, only: equil_par_curr, n_fields, currdens_par, xy_local
  use par_in, only: fourier2d, beta
  use coordinates, only: kx, ky, kjmin, vp
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use blockindex
  use fourier
  use prefactors
  use gyro_average_ff_mod, only: jfac, I1_factor
  USE x_derivatives, only: x_deriv_exc
  use geometry, only: geom, C_xy
  use mpi
  use ff_nonlinear_term_mod
  implicit none


  type, public, extends(ff_nonlinear_term_t) :: ff_yx_nonlinear_term_t
   contains
     procedure :: initialize => initialize_nonlinearity
     procedure :: finalize => finalize_nonlinearity
     procedure, private :: calc => calc_nonlinearity_ff_t
     procedure :: getType => getThisType
     procedure :: to_direct_xy_t => nl_to_direct_xy_t
  end type ff_yx_nonlinear_term_t

  private

contains
  
  function getThisType(this)
    class(ff_yx_nonlinear_term_t) :: this
    character(MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "ff_yx_nonlinear_term_t"
  end function getThisType


  subroutine initialize_nonlinearity(this)
    class(ff_yx_nonlinear_term_t) :: this

    REAL :: B_Bstar
    Integer :: i,j,k,l,n
    complex:: prefax, prefay
    integer:: modenumx

    if (this%init_status==1) return
    !call initialize_da_bounds
    call initialize_fourier(li0da,ly0da)

    if (equil_par_curr) then 
       Allocate(this%pnl(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,ln1:ln2))
       do n=ln1,ln2
          do l=ll1,ll2
             do k=lk1,lk2
                do j=pj1,pj2
                   do i=pi1,pi2
                      B_Bstar = 1.0D0/(1.0D0+beta*&
                           &sqrt(0.5*spec(n)%mass*spec(n)%temp)*vp(l)*&
                           &currdens_par(i)/(spec(n)%charge*geom%Bfield(i,j,k)**2))
                      this%pnl(i,j,k,l,n) = B_Bstar/C_xy(i)
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       Allocate(this%pnl_1d(pi1:pi2))
       this%pnl_1d(pi1:pi2) = 1.0/C_xy(pi1:pi2)
    endif

    !----------------------------------------------
    !dealiasing
    
    allocate(this%ranshift_to(li1:li2,lj1:lj2), this%ranshift_back(li1:li2,lj1:lj2))
    allocate(this%shiftmat(li1:li2,lj1:lj2,0:3))
    prefax=imag*pi/2/nx0
    prefay=imag*pi/2/(2*nky0)
    do j=lj1,lj2
       do i=li1,li2
          modenumx=mod(j+hkx,nx0)-hkx
          this%shiftmat(i,j,0)=exp(modenumx*prefax+i*prefay)
          this%shiftmat(i,j,1)=exp(modenumx*prefax-i*prefay)
          this%shiftmat(i,j,2)=exp(-modenumx*prefax-i*prefay)
          this%shiftmat(i,j,3)=exp(-modenumx*prefax+i*prefay)
       enddo
    enddo
    if (evenx.eq.1) then
       this%shiftmat(:,hkx+1,:)=0
    endif
    !---------------------------------------------------

    this%init_status = 1
  end subroutine initialize_nonlinearity

  subroutine finalize_nonlinearity(this)
    class(ff_yx_nonlinear_term_t) :: this

    deallocate(this%shiftmat)
    deallocate(this%ranshift_to, this%ranshift_back)

    if (equil_par_curr) then 
       deallocate(this%pnl)
    else
       deallocate(this%pnl_1d)
    endif
    call finalize_fourier

    this%init_status = 0
  end subroutine finalize_nonlinearity


  subroutine calc_nonlinearity_ff_t(this,p_g,p_emfields,p_barchi,p_rhs,stage,lb1,lb2)
    class(ff_yx_nonlinear_term_t) :: this
    Integer,intent(in) :: lb1,lb2
    complex,Dimension(li1:li2,lj1:lj2,1:lklmn0),intent(in) :: p_g
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    Complex, Dimension(:,:,:), pointer :: p_barchi
    Complex,Dimension(li1:li2,lj1:lj2,lb1:lb2),intent(inout) :: p_rhs
    integer,intent(in):: stage

    complex,Dimension(li1:li2):: val
    complex, Dimension(li1:li2+1,lj1:lj2):: ve_x, ve_y, dG_x, dG_y
    Real, Dimension(0:ly0da-1,0:lj0-1) :: nonlin_real
    Complex,Dimension(li1:li2+1,lj1:lj2) :: nonlin
    integer:: j,k,l,m,n,klmn
    logical:: first

    PERFON_I('calcnl')
    do klmn=lb1,lb2

!       PERFON_I('recve')

       ! calculates nonlin_real = -ve_y_real*dG_y_real + ve_x_real*dG_x_real
       !compute ve_x,ve_y
       if (associated(p_barchi).and.(stage.ne.0)) then
          Do j=lj1,lj2
             ve_x(li1:li2,j) = imag*ky(li1:li2)*p_barchi(:,j,klmn-lb1+1)* this%ranshift_to(:,j)
             ve_y(li1:li2,j) = imag*kx(j)*p_barchi(:,j,klmn-lb1+1)* this%ranshift_to(:,j)
          enddo          
       else
          k=sk(klmn)
          l=sl(klmn)
          m=sm(klmn)
          n=sn(klmn)
          if (n_fields.eq.1) then
             Do j=lj1,lj2
                val= p_emfields(:,j,k,1)*jfac(:,j,k,m,n)
                ve_x(li1:li2,j) = imag*ky(li1:li2)*val* this%ranshift_to(:,j)
                ve_y(li1:li2,j) = imag*kx(j)*val* this%ranshift_to(:,j)
             enddo
          else
             Do j=lj1,lj2
                val = (p_emfields(:,j,k,1) - &
                     vTvpar(l,n) * p_emfields(:,j,k,2))*jfac(:,j,k,m,n)
                if (n_fields.gt.2) then
                   val = val + mu_tjqj(m,n) * p_emfields(:,j,k,3)*&
                        I1_factor(:,j,k,m,n)
                endif
                ve_x(li1:li2,j) = imag*ky(li1:li2)*val* this%ranshift_to(:,j)
                ve_y(li1:li2,j) = imag*kx(j)*val* this%ranshift_to(:,j)
             enddo
          endif
       endif
       !      PERFOFF_I
 !      PERFON_I('recdg')
       do j=lj1,lj2
          dG_x(li1:li2,j) = imag*kx(j)*p_g(:,j,klmn)* this%ranshift_to(:,j) 
          dG_y(li1:li2,j) = imag*ky(li1:li2)*p_g(:,j,klmn)* this%ranshift_to(:,j)
       end do
 !      PERFOFF_I

 !      PERFON_I('rest')
       first=(stage.eq.1)
       Call this%to_direct_xy_t(ve_x,dG_x,0,nonlin_real,first)
       Call this%to_direct_xy_t(ve_y,dG_y,-1,nonlin_real,first)
       Call fft_xy_to_ff_t(nonlin_real,nonlin)

       !parallel current, thus B_{0\parallel}^\ast
       IF (equil_par_curr) nonlin(li1:li2,:) = nonlin(li1:li2,:) * this%pnl(pi1,pj2,k,l,n)

       p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + nonlin(li1:li2,:)*this%ranshift_back

    enddo
    PERFOFF_I
  end subroutine calc_nonlinearity_ff_t

  !
  ! compute NL contributions in Fourier space.
  !
  Subroutine nl_to_direct_xy_t(this,inarr1,inarr2,was, outarr,first)
    class(ff_yx_nonlinear_term_t) :: this
    Complex,dimension(li1:li2+1, lj1:lj2), Intent(InOut):: inarr1, inarr2
    Real, dimension(0:ly0da-1, 0:lj0-1),Intent(inout):: outarr
    Integer,intent(in):: was
    Logical,intent(in):: first

    Real, dimension(0:ly0da-1,0:lj0-1) :: rearr1,rearr2
    Real:: ve_x_max_loc, ve_y_max_loc

    ! Zero the extra fourier mode
    inarr1(ly0da/2,:)=0.0
    inarr2(ly0da/2,:)=0.0

    ! Do inverse fourier transform to get real values
     call fft_ff_to_xy_t(inarr1,rearr1)
     call fft_ff_to_xy_t(inarr2,rearr2)

    ! Find maximum for ve_x, ve_y
    if (first) then
       if (was.eq.0) then
          ve_x_max_loc=maxval(rearr1(0:ly0da-1,:))
          ve_max(1)=max(ve_max(1),ve_x_max_loc)
       else
          ve_y_max_loc=maxval(rearr1(0:ly0da-1,:))
          ve_max(2)=max(ve_max(2),ve_y_max_loc)
       endif
    endif

    !multiply the real arrays
    if (was.eq.0) then
       outarr=rearr1*rearr2
    elseif (was.eq.-1) then
       outarr=outarr-rearr1*rearr2
    endif

  end subroutine nl_to_direct_xy_t

  !
  ! compute NL contributions in Fourier space using 1d fourier transforms.
  !
  Subroutine nl_to_direct_xy_1d(this,inarr,rearr)

    class(ff_yx_nonlinear_term_t) :: this
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

  !
  ! compute NL contributions in Fourier space using 2d fourier transforms.
  !
  Subroutine nl_to_direct_xy_2d(this,inarr,rearr)

    class(ff_yx_nonlinear_term_t) :: this
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
  Subroutine nl_to_fourier_xy(this,inarr, outarr)

    class(ff_yx_nonlinear_term_t) :: this
    real, Intent(In) :: inarr(0:ly0da-1, 0:li0da/n_procs_y-1)
    complex, Intent(Out) :: outarr(li1:li2, lj1:lj2)

    if(fourier2d) then ! Choose 1D or 2D FFTs
       call this%to_fourier_xy_2d(inarr, outarr)
    else
       call this%to_fourier_xy_1d(inarr, outarr)
    end if

  End Subroutine nl_to_fourier_xy

  !
  ! Transform the nonlinearity back to Fourier space using 1d fourier transforms
  !
  Subroutine nl_to_fourier_xy_1d(this,inarr, outarr)

    class(ff_yx_nonlinear_term_t) :: this
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
   Subroutine nl_to_fourier_xy_2d(this,inarr, outarr)

    class(ff_yx_nonlinear_term_t) :: this
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

end module ff_yx_nonlinear_term_mod
