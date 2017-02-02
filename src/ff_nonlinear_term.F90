#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

module ff_nonlinear_term_mod
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
  use nonlinear_term_mod
  implicit none


  type, public, extends(nonlinear_term_t) :: ff_nonlinear_term_t
     complex,dimension(:,:,:), allocatable :: shiftmat
     complex,dimension(:,:),allocatable:: ranshift_to, ranshift_back
     complex,dimension(:),allocatable:: zeros
     !complex,dimension(:),allocatable:: val
     
     complex,dimension(:),allocatable:: ikxda  
     Real, Dimension(:,:,:,:,:), Allocatable :: pnl !nonlinearity prefactor
     Real, Dimension(:), Allocatable :: pnl_1d !nonlinearity prefactor if B_{0\}_^*=B_0
     integer:: lzero1, lzero2
     real :: ve_max_omp_x, ve_max_omp_y

   contains
     procedure :: initialize => initialize_nonlinearity
     procedure :: finalize => finalize_nonlinearity
     procedure :: add => add_nonlinearity
     procedure,private :: calc => calc_nonlinearity_ff
     procedure :: mem_est => mem_est_nonlinearity
     procedure :: construct => construct_ff_nonlinear_term
     procedure :: destruct => destruct_ff_nonlinear_term
     procedure :: getType => getThisType
     procedure:: to_direct_xy => nl_to_direct_xy 
     procedure:: to_direct_xy_1d => nl_to_direct_xy_1d
     procedure:: to_direct_xy_2d => nl_to_direct_xy_2d
     procedure:: to_fourier_xy => nl_to_fourier_xy 
     procedure:: to_fourier_xy_1d => nl_to_fourier_xy_1d
     procedure:: to_fourier_xy_2d => nl_to_fourier_xy_2d

  end type ff_nonlinear_term_t

  private

contains
  
  function getThisType(this)
    class(ff_nonlinear_term_t) :: this
    character(MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "ff_nonlinear_term_t"
  end function getThisType

  Real function mem_est_nonlinearity(this,mem_req_in)
    class(ff_nonlinear_term_t) :: this
    real,intent(IN):: mem_req_in
    
    mem_est_nonlinearity = mem_est_nonlinearity_ff(mem_req_in)

    !pnl
    if (equil_par_curr) mem_est_nonlinearity=mem_est_nonlinearity+&
         &SIZE_OF_REAL_MB*pi0*pj0*lk0*ll0*ln0
    
  end function mem_est_nonlinearity

  ! In the constructor we do work to instantiate the whole
  ! object, mainly for being able to call later on the
  ! mem_est routines.
  ! construct is not equal to initialize. construct just sets up
  ! the object, but does not allocate any data arrays.
  subroutine construct_ff_nonlinear_term(this,equil_par_curr)
    class(ff_nonlinear_term_t) :: this
    logical :: equil_par_curr
  end subroutine construct_ff_nonlinear_term

  subroutine destruct_ff_nonlinear_term(this)
    class(ff_nonlinear_term_t) :: this
    !if (allocated(this%prefactor)) deallocate(this%prefactor)
  end subroutine destruct_ff_nonlinear_term


  subroutine initialize_nonlinearity(this)
    class(ff_nonlinear_term_t) :: this

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

    !-------------------------------------
    !dealiasing
    if (turbdeal) then
       allocate(this%ranshift_to(li1:li2,lj1:lj2), this%ranshift_back(li1:li2,lj1:lj2))
       allocate(this%shiftmat(li1:li2,lj1:lj2,0:3))
       prefax=imag*pi/2/nx0
       prefay=imag*pi/2/(2*nky0)
       do j=lj1,lj2
          do i=li1,li2
             modenumx=mod(i+hkx,nx0)-hkx
             this%shiftmat(i,j,0)=exp(modenumx*prefax+j*prefay)
             this%shiftmat(i,j,1)=exp(modenumx*prefax-j*prefay)
             this%shiftmat(i,j,2)=exp(-modenumx*prefax-j*prefay)
             this%shiftmat(i,j,3)=exp(-modenumx*prefax+j*prefay)
          enddo
       enddo
       if (evenx.eq.1) then
          this%shiftmat(hkx+1,:,:)=0
       endif
    endif

    this%lzero1=hkx+1
    this%lzero2=lkx+kx_offset-1
    allocate(this%zeros(this%lzero1:this%lzero2))
    this%zeros=(0.,0.)
    !allocate(this%val(li1da:li2da))
    !this%val=(0.,0.)
    allocate(this%ikxda(li1da:li2da))
    this%ikxda=0.
    do i=li1,hkx
       this%ikxda(i)=imag*kx(i)
    enddo
    do i=lkx,li2
       this%ikxda(i+kx_offset)=imag*kx(i)
    enddo
    
    !------------------------
    this%init_status = 1
  end subroutine initialize_nonlinearity

  subroutine add_nonlinearity(this,g_block,p_dgdxy,p_emfields,p_barchi,p_dbarchidxy,p_rhs,lb1,lb2,stage)
    class(ff_nonlinear_term_t) :: this
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout), target:: p_rhs    
    integer, intent(in):: stage
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    Complex, Dimension(:,:,:), pointer :: p_barchi
    Complex, Dimension(:,:,:,:), pointer :: p_dbarchidxy, p_dgdxy
    complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in), target:: g_block

    real, dimension(2) :: ve_max_sav
    integer:: modenumx
    integer:: i,j
    logical:: first

    if (stage.eq.1) then
       first=.true.
       ve_max_sav = ve_max
    else
       first=.false.
    end if

    if (turbdeal) then
       !the shift matrix does not depend on the klmn index and has to be computed only for the first klmn block
#ifndef WITHOMP_BLOCKLOOP
       if(lb1.eq.1) then
#endif
          if(yx_order) then
             do j=lj1,lj2
                do i=li1,li2
                   modenumx=mod(j+hkx,nx0)-hkx
                   this%ranshift_to(i,j)=exp(imag*2.*pi*(modenumx*phase_1+i*phase_2))*this%shiftmat(i,j,stage-1)
                enddo
             enddo
          else
             do j=lj1,lj2
                do i=li1,li2
                   modenumx=mod(i+hkx,nx0)-hkx
                   this%ranshift_to(i,j)=exp(imag*2.*pi*(modenumx*phase_1+j*phase_2))*this%shiftmat(i,j,stage-1)
                enddo
             enddo
          end if
          this%ranshift_back=conjg(this%ranshift_to) 
#ifndef WITHOMP_BLOCKLOOP
       endif
#endif
    endif

    ! add the ExB nonlinearity, calculation is performed in real space
    call this%calc(g_block,p_emfields,p_barchi,p_rhs,stage,lb1,lb2)
    !if(yx_order) then
    !   call calc_nonlinearity_ff_t(g_block,p_emfields,p_barchi,p_rhs,stage,lb1,lb2)
    !else
    !   call calc_nonlinearity_ff(g_block,p_emfields,p_barchi,p_rhs,stage,lb1,lb2)
    !endif
    
    if (first) then
       ! now in ve_max we find the maximal ExB velocity in the whole 
       ! volume. The reduction on the phase space has already been 
       ! performed, but we have to do it still over the species.
       ve_max(1)=max(ve_max(1),ve_max_sav(1))
       ve_max(2)=max(ve_max(2),ve_max_sav(2))
       !exchange only at the end of stage 1
       if(lb2.eq.lklmn0) then
          Call my_2reals_max_to_all(ve_max)
       endif
    endif
  end subroutine add_nonlinearity

  subroutine finalize_nonlinearity(this)
    class(ff_nonlinear_term_t) :: this

    if (turbdeal) then
       deallocate(this%shiftmat)
       deallocate(this%ranshift_to, this%ranshift_back)
    end if
    deallocate(this%zeros,this%ikxda)
    !deallocate(this%val)

    if (equil_par_curr) then 
       deallocate(this%pnl)
    else
       deallocate(this%pnl_1d)
    endif
    call finalize_fourier

    this%init_status = 0
  end subroutine finalize_nonlinearity

!--------------------------------------------------------------------------------------

  Real function mem_est_nonlinearity_ff(mem_req_in)
    real:: mem_req_in

    !dealiasing
    if (turbdeal) then
       mem_est_nonlinearity_ff = mem_req_in + &
            &6.*lij0*SIZE_OF_COMPLEX_MB
    else
       mem_est_nonlinearity_ff = mem_req_in + &
            &(li0da*(lj0+1)+li0)*SIZE_OF_COMPLEX_MB
    end if

    !in calc_nonlinearity
    mem_est_nonlinearity_ff = mem_est_nonlinearity_ff + &
         & (4.*li0da*lj0+0.5*ly0da*li0da/n_procs_y+lij0)*SIZE_OF_COMPLEX_MB

    !in nl_to_direct_xy / nl_to_fourier_xy
    mem_est_nonlinearity_ff = mem_est_nonlinearity_ff + ((ly0da/2+1)*(li0da/n_procs_y)+&
         & lj0*li0da+li0da*lj0+ly0da*li0da/n_procs_y)*SIZE_OF_COMPLEX_MB

  end function mem_est_nonlinearity_ff
  

  subroutine calc_nonlinearity_ff(this,p_g,p_emfields,p_barchi,p_rhs,stage,lb1,lb2)
    class(ff_nonlinear_term_t) :: this
    Integer,intent(in) :: lb1,lb2
    complex,Dimension(li1:li2,lj1:lj2,1:lklmn0),intent(in) :: p_g
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    Complex, Dimension(:,:,:), pointer :: p_barchi
    Complex,Dimension(li1:li2,lj1:lj2,lb1:lb2),intent(inout) :: p_rhs
    integer,intent(in):: stage

    complex, Dimension(li1da:li2da,lj1:lj2):: ve_x, ve_y, dG_x, dG_y
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: nonlin_real
    Complex,Dimension(li1:li2,lj1:lj2) :: nonlin
    complex,dimension(li1da:li2da) :: val
    integer:: j,k,l,m,n,klmn,my_thread
#ifdef WITHOMP_NONLIN
    integer :: omp_get_thread_num
#endif
    logical:: first

    PERFON_I('calcnl')
    first=(stage.eq.1)
    this%ve_max_omp_x = 0.0
    this%ve_max_omp_y = 0.0
#ifdef WITHOMP_NONLIN
    !$OMP PARALLEL default(none) &
    !$OMP shared(li1,li2,hkx,kx_offset,lkx,li2da,this,n_fields,jfac,vtvpar,mu_tjqj)&
    !$OMP shared(I1_Factor,kx,p_g,turbdeal,equil_par_curr,pi1,pj2,p_rhs,mype)&
    !$OMP shared(lb1,lb2,lj1,lj2,sk,sl,sm,sn,stage,p_barchi,ky,p_emfields,first)&
    !$OMP private(klmn,j,k,l,m,n,ve_x,ve_y,dG_x,dG_y,nonlin_real,nonlin,val,my_thread) &
    !$OMP reduction(MAX:this%ve_max_omp_x) &
    !$OMP reduction(MAX:this%ve_max_omp_y)
    if (first) val=cmplx(0.0,0.0)
    my_thread = omp_get_thread_num()
    !$OMP DO
#else
    my_thread = 0
#endif
    do klmn=lb1,lb2

       PERFON_I('recve')
       k=sk(klmn)
       l=sl(klmn)
       m=sm(klmn)
       n=sn(klmn)
       ! calculates nonlin_real = -ve_y_real*dG_y_real + ve_x_real*dG_x_real
       if (associated(p_barchi).and.(stage.ne.0)) then
          !write(*,"(2I3,A,ES17.10)") mype,my_thread,", p_barchi = ",real(sum(conjg(p_barchi)*p_barchi))
          Do j=lj1,lj2
             ve_x(li1:hkx,j) = imag*ky(j) * p_barchi(li1:hkx,j,klmn-lb1+1)
             ve_x(hkx+1:lkx+kx_offset-1,j) = (0.0,0.0)
             ve_x(lkx+kx_offset:li2da,j) = imag*ky(j) * p_barchi(lkx:li2,j,klmn-lb1+1)
             ve_y(li1:hkx,j) = this%ikxda(li1:hkx) * p_barchi(li1:hkx,j,klmn-lb1+1)
             ve_y(hkx+1:lkx+kx_offset-1,j) = (0.0,0.0)
             ve_y(lkx+kx_offset:li2da,j) = this%ikxda(lkx+kx_offset:li2da) * p_barchi(lkx:li2,j,klmn-lb1+1)
          enddo
       else
          !compute ve_x,ve_y
          !write(*,"(2I3,A,ES17.10,A,6I3)") mype,my_thread,", p_emfields 1 = ",&
          !     &real(sum(conjg(p_emfields(:,:,:,1))*p_emfields(:,:,:,1))),&
          !     &", li1,hkx,lkx,kx_offset,li2,li2da, = ",li1,hkx,lkx,kx_offset,li2,li2da
          
          if (n_fields.eq.1) then
             Do j=lj1,lj2
                val(li1:hkx)= p_emfields(li1:hkx,j,k,1)*jfac(li1:hkx,j,k,m,n)
                val(hkx+1:lkx+kx_offset-1) = (0.0,0.0)
                val(lkx+kx_offset:li2da)= p_emfields(lkx:li2,j,k,1)*jfac(lkx:li2,j,k,m,n)
                ve_x(:,j) = imag*ky(j)*val
                ve_y(:,j) = this%ikxda(:)*val
             enddo
          else
             Do j=lj1,lj2
                val(li1:hkx)= (p_emfields(li1:hkx,j,k,1) - &
                     vTvpar(l,n) * p_emfields(li1:hkx,j,k,2) )*jfac(li1:hkx,j,k,m,n)
                val(hkx+1:lkx+kx_offset-1) = (0.0,0.0)
                val(lkx+kx_offset:li2da)=( p_emfields(lkx:li2,j,k,1) - &
                     vTvpar(l,n) * p_emfields(lkx:li2,j,k,2) )*jfac(lkx:li2,j,k,m,n)
                if (n_fields.gt.2) then
                   val(li1:hkx) = val(li1:hkx) + mu_tjqj(m,n) * p_emfields(li1:hkx,j,k,3)*&
                        I1_factor(li1:hkx,j,k,m,n)
                   val(lkx+kx_offset:li2da) = val(lkx+kx_offset:li2da) + mu_tjqj(m,n) * p_emfields(lkx:li2,j,k,3)*&
                        I1_factor(lkx:li2,j,k,m,n)
                endif
                ve_x(:,j) = imag*ky(j)*val
                ve_y(:,j) = this%ikxda(:)*val
             enddo
          endif
       endif
       
 !      PERFOFF_I
       
 !      PERFON_I('recdg')
       do j=lj1,lj2
          dG_x(li1:hkx,j) = imag*kx(li1:hkx)*p_g(li1:hkx,j,klmn)
          dG_x(this%lzero1:this%lzero2,j)=this%zeros
          dG_x(lkx+kx_offset:li2da,j) = imag*kx(lkx:li2)*p_g(lkx:li2,j,klmn) 
          
          dG_y(li1:hkx,j) = imag*ky(j)*p_g(li1:hkx,j,klmn) 
          dG_y(this%lzero1:this%lzero2,j)=this%zeros
          dG_y(lkx+kx_offset:li2da,j) = imag*ky(j)*p_g(lkx:li2,j,klmn)
       end do
 !      PERFOFF_I
       
 !      PERFON_I('shift')
       if (turbdeal) then
          ve_x=ve_x* this%ranshift_to
          ve_y=ve_y* this%ranshift_to
          dG_x=dG_x* this%ranshift_to
          dG_y=dG_y* this%ranshift_to 
       endif
       PERFOFF_I

       !write(*,"(2I3,2(A,2ES17.10))") mype,my_thread,&
       !     &", ve_x,ve_y = ",sqrt(real(sum(conjg(ve_x)*ve_x))),sqrt(real(sum(conjg(ve_y)*ve_y))),&
       !     &", dG_x,dG_y = ",sqrt(real(sum(conjg(dG_x)*dG_x))),sqrt(real(sum(conjg(dG_y)*dG_y)))
       Call this%to_direct_xy(ve_x,dG_x,0,nonlin_real,first)
       !write(*,"(2I3,A,ES17.10)") mype,my_thread,", nonlin_real 1 = ",sqrt(sum(nonlin_real**2))
       Call this%to_direct_xy(ve_y,dG_y,-1,nonlin_real,first) 
       !write(*,"(2I3,A,ES17.10)") mype,my_thread,", nonlin_real 2 = ",sqrt(sum(nonlin_real**2))
       Call this%to_fourier_xy(nonlin_real, nonlin)

       !parallel current, thus B_{0\parallel}^\ast
       IF (equil_par_curr) nonlin = nonlin * this%pnl(pi1,pj2,k,l,n)

       if (turbdeal) then
          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + nonlin*this%ranshift_back
       else
          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + nonlin
       endif
    enddo
#ifdef WITHOMP_NONLIN
    !$OMP END DO
    !$OMP END PARALLEL
#endif
    if (first) then
       ve_max(1) = max(ve_max(1),this%ve_max_omp_x)
       ve_max(2) = max(ve_max(2),this%ve_max_omp_y)
       !write(*,"(I3,A,2ES10.3)") mype,"ve_max = ",ve_max(1),ve_max(2)
    end if
    PERFOFF_I

  end subroutine calc_nonlinearity_ff

  Subroutine nl_to_direct_xy(this,inarr1,inarr2,was, outarr,first)
    class(ff_nonlinear_term_t) :: this
    Complex,dimension(li1da:li2da, lj1:lj2), Intent(In):: inarr1, inarr2
    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1),Intent(inout):: outarr
    integer,intent(in):: was
    Logical,intent(in):: first

    Real, dimension(0:ly0da-1, 0:li0da/n_procs_y-1) :: rearr1,rearr2
    Real:: ve_x_max_loc, ve_y_max_loc

    ! Do inverse fourier transform to get real values
    if(fourier2d) then
       call this%to_direct_xy_2d(inarr1,rearr1) ! ve_x/y
       call this%to_direct_xy_2d(inarr2,rearr2) ! dG_x/y
    else
       call this%to_direct_xy_1d(inarr1,rearr1) ! ve_x/y
       call this%to_direct_xy_1d(inarr2,rearr2) ! dG_x/y
    end if

    ! Find maximum for ve_x, ve_y
    if (first) then
       if (was.eq.0) then
          ve_x_max_loc=maxval(rearr1)
          !ve_x_max_loc=maxval(abs(rearr1))
          !ve_max(1)=max(ve_max(1),ve_x_max_loc)
          this%ve_max_omp_x=max(this%ve_max_omp_x,ve_x_max_loc)
          !write(*,"(2I3,A,ES10.3)") mype,omp_get_thread_num(),"ve_max_omp_x = ",ve_max_omp_x
       else
          ve_y_max_loc=maxval(rearr1)
          !ve_y_max_loc=maxval(abs(rearr1))
          !ve_max(2)=max(ve_max(2),ve_y_max_loc)
          this%ve_max_omp_y=max(this%ve_max_omp_y,ve_y_max_loc)
          !write(*,"(2I3,A,ES10.3)") mype,omp_get_thread_num(),"ve_max_omp_y = ",ve_max_omp_y
       endif
    endif

    !multiply the real arrays
    if (was.eq.0) then
       outarr=rearr1*rearr2
    elseif (was.eq.-1) then
       outarr=outarr-rearr1*rearr2
    endif

  end subroutine nl_to_direct_xy


  !
  ! compute NL contributions in Fourier space using 1d fourier transforms.
  !
  Subroutine nl_to_direct_xy_1d(this,inarr,rearr)
    class(ff_nonlinear_term_t) :: this

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
    class(ff_nonlinear_term_t) :: this

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
    class(ff_nonlinear_term_t) :: this

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
    class(ff_nonlinear_term_t) :: this

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
    class(ff_nonlinear_term_t) :: this

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

end module ff_nonlinear_term_mod
