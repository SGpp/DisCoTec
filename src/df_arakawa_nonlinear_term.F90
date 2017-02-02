#include "redef.h"
#include "intrinsic_sizes.h"

module df_arakawa_nonlinear_term_mod
  use df_nonlinear_term_mod

  use par_mod, only: imag, pi, ve_max
  use coordinates, only: kjmin
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use blockindex, only: sk,sl,sn
  use fourier
  use prefactors
  USE x_derivatives, only: x_deriv_exc
  use geometry, only: geom, C_xy
  use mpi
  implicit none

  type, public, extends(df_nonlinear_term_t) :: df_arakawa_nonlinear_term_t
   contains
     procedure :: calc => calc_arakawa_nonlinearity_df
     procedure :: getType => getThisType
  end type df_arakawa_nonlinear_term_t

private
  complex,dimension(:),allocatable:: ranshift_to_df, ranshift_back_df

contains

  function getThisType(this)
    class(df_arakawa_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_arakawa_nonlinear_term_t"
  end function getThisType

  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_arakawa_nonlinearity_df(this,gy_chi,g_block,vexy,dgdxy,localrhs,first,lb1,lb2)
    class(df_arakawa_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in):: first
    Integer, Intent(in):: lb1,lb2

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    !Complex, Dimension(0:nj0-1, 0:li0/n_procs_y-1) :: tmp_arr
    Complex, Dimension(0:nj0-1, 0:li0/n_procs_y-1,1:this%lbg0) :: nonlin3,tmp_arr
    complex, dimension(li1:li2, lj1:lj2,1:this%lbg0):: g2d
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,2,1:this%lbg0) ::  vexy_re, dgdxy_re
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,1:this%lbg0) ::  gy_chi_re, g2d_re
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,1:this%lbg0) ::  nonlin1_re, nonlin2_re, nonlin3_re
    Real, Dimension(0:li0/n_procs_y-1,0:ly0da-1):: tmp
    Real, Dimension(-nib:li0-1+nib,0:ly0da-1):: nonlin2_exc
    Real, Dimension(0:li0-1,0:ly0da-1):: nonlin2_dx
    Real, Dimension(1:2)  ::  vexy_max_loc
    Integer :: i,j, ierr,klmn
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num
    my_thread = omp_get_thread_num()
#endif

    if (turbdeal) then
       do j=lj1,lj2
          vexy(:,j,:,1:this%lbg0)=vexy(:,j,:,1:this%lbg0)* ranshift_to_df(j)
          dgdxy(:,j,:,1:this%lbg0)=dgdxy(:,j,:,1:this%lbg0)* ranshift_to_df(j)
       enddo
    endif

    ! anti-aliasing and y fourier transform
    PERFON_I('deal_FFT')
    Call gto_real(this,vexy,vexy_re,2*this%lbg0)
    !Call gto_real(ve_y,ve_y_re)
    Call gto_real(this,dgdxy,dgdxy_re,2*this%lbg0)
    !Call gto_real(dg_dy,dg_dy_re)
    PERFOFF_I

    ! get max ExB velocity for timestep estimate
    if (first) then
       call this%prefactor%multiply_max(vexy_max_loc,vexy_re,lb1,lb2)
       ve_max(1)=max(ve_max(1),vexy_max_loc(1)) !ve_x
       ve_max(2)=max(ve_max(2),vexy_max_loc(2)) !ve_y
    end if

    ! compute the 'standard' nonlinear term, which we can also use here 
    ! in the Arakawa representation
    PERFON_I('comp_nl')
    nonlin1_re = -vexy_re(:,:,1,:)*dgdxy_re(:,:,2,:) + vexy_re(:,:,2,:)*dgdxy_re(:,:,1,:)

    PERFON_I('add_Arakawa')
    if (turbdeal) then
       do j=lj1,lj2
          g2d(:,j,1:this%lbg0)=g_block(li1:li2,j,1:this%lbg0)* ranshift_to_df(j)
          gy_chi(:,j,1:this%lbg0)=gy_chi(:,j,1:this%lbg0)* ranshift_to_df(j)
       enddo
    else
       ! the following statement has been moved from df_nonlinear_term%add, as
       ! it is needed only for arakawa scheme
       g2d = g_block(li1:li2,:,1:this%lbg0)
    endif
    !we need the distribution function and the potential in real space       
    Call gto_real(this,g2d,g2d_re,this%lbg0)
    Call gto_real(this,gy_chi,gy_chi_re,this%lbg0)
    
    nonlin2_re = g2d_re*vexy_re(:,:,2,:) - gy_chi_re*dgdxy_re(:,:,2,:)
    nonlin3_re = gy_chi_re*dgdxy_re(:,:,1,:) - g2d_re*vexy_re(:,:,1,:)
    
    !compute x derivative of nonlin_re2 with 2D-arrays for faster x parallelization
    !PERFON_I('dx_Arakawa')
    do klmn=1,this%lbg0
       !transpose arrays for mpi exchange 
       do i=0,ly0da-1
          do j=0,li0/n_procs_y-1
             tmp(j,i)=nonlin2_re(i,j,klmn)
          enddo
       enddo
       if (n_procs_y.eq.1) then
          do i=0,ly0da-1
             nonlin2_exc(0:li0-1,i) = tmp(0:li0-1,i)
          enddo
       else
          !gather from all y processes
          do i=0,ly0da-1
#ifdef WITHOMP_BLOCKLOOP
             call mpi_allgather(tmp(0,i),li0/n_procs_y, MPI_REAL_TYPE,&
                  nonlin2_exc(0,i),li0/n_procs_y,MPI_REAL_TYPE, threadlocal_mpi_comm_y, ierr) 
#else
             call mpi_allgather(tmp(0,i),li0/n_procs_y, MPI_REAL_TYPE,&
                  nonlin2_exc(0,i),li0/n_procs_y,MPI_REAL_TYPE, mpi_comm_y, ierr) 
#endif
          enddo
       end if
       
       call x_deriv_exc(nonlin2_exc,nonlin2_dx,li0)
          
       !another transposition
       do i=0,ly0da-1
          do j=0,li0/n_procs_y-1
             nonlin2_re(i,j,klmn)=nonlin2_dx(j+li0/n_procs_y*my_pey,i)
          enddo
       enddo
    end do
    !PERFOFF_I
    
    do klmn=1,this%lbg0
       ! fourier transform nonlin_re3 back to ky space and compute its y derivative
       Call to_fourier_y(nonlin3_re(:,:,klmn),nonlin3(:,:,klmn))
    end do
    Do j=0,nj0-1
       nonlin3(j,:,:)=imag*kjmin*j*nonlin3(j,:,:)
    Enddo
    
    !add the two parts which can be fourier transformed in one step
    nonlin1_re=nonlin1_re+nonlin2_re
    PERFOFF_I
    
    do klmn=1,this%lbg0
       ! transform back to fourier space
       Call to_fourier_y(nonlin1_re(:,:,klmn),tmp_arr(:,:,klmn))
    end do
    ! add all terms and divide by three to calculate the average
    tmp_arr=1/3.*tmp_arr+1/3.*nonlin3
    
    ! Transpose and remove zeros in y for dealiasing
    do klmn=1,this%lbg0
       Call transpose_cmplx(nj0, li0, tmp_arr(:,:,klmn), 0, nonlin(:,:,klmn), 0)
    end do
    
    PERFOFF_I

    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0,lb1,lb2)

  End Subroutine calc_arakawa_nonlinearity_df

  !> Go to real space
  !! This routine is equal to the one of df_nonlinear_term, with the
  !! small exception, that here it is clear, that we do not have
  !! explicit x-dealiasing. This saves us one copy operation.
  subroutine gto_real(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1:li2,lj1:lj2,1:howmany) :: tmp_arr1
    Complex, Dimension(0:nj0-1, 0:li0/n_procs_y-1) :: tmp_arr2
    real, dimension(:,:), pointer :: p_out
    Integer:: klmn, i_block,number_of_lbg0_blocks
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num

    my_thread = omp_get_thread_num()
#endif
    number_of_lbg0_blocks = howmany/this%lbg0
    if (number_of_lbg0_blocks.eq.0) then
       ! Transpose x-y
       Call transpose_cmplx(li0, nj0, inarr(:,:,1), 0,tmp_arr2 , 0)
       ! Fourier transfrom in y (include dealiasing step)
       Call to_real_y(tmp_arr2,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks

          ! Transpose x-y
          do klmn=1,this%lbg0
             Call transpose_cmplx(li0, nj0, inarr(:,:,(i_block-1)*this%lbg0+klmn), 0,tmp_arr2 , 0)

             ! Fourier transfrom in y (include dealiasing step)
             p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
             Call to_real_y(tmp_arr2,p_out) !outarr(:,:,(i_block-1)*this%lbg0+klmn))
          end do
       end do
    end if

  end subroutine gto_real

end module df_arakawa_nonlinear_term_mod
