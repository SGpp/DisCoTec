#include "redef.h"
#include "intrinsic_sizes.h"

module df_acc_nonlinear_term_mod
  use df_nonlinear_term_mod
  use cufft

  use par_mod, only: imag, pi, ve_max, phase_1, phase_2, spec
  use par_other, only: n_fields
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use fourier
  use mpi
  implicit none

  type, public, extends(df_nonlinear_term_t) :: df_acc_nonlinear_term_t
   contains
     procedure :: initialize => initialize_acc_nonlinearity
     procedure :: finalize => finalize_acc_nonlinearity
     procedure :: calc => calc_acc_nonlinearity_df
     procedure :: getType => getThisType
  end type df_acc_nonlinear_term_t
  private

  !Complex, Dimension(:,:,:) :: dev_full_tmp2

contains

  function getThisType(this)
    class(df_acc_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_acc_nonlinear_term_t"
  end function getThisType
    
  subroutine initialize_acc_nonlinearity(this)
    class(df_acc_nonlinear_term_t) :: this

    if (this%init_status==1) return
    !write(*,"(3A,I5)") "Initializing ",this%getType()," with blocksize = ",this%lbg0
    call initialize_fourier(li0da,ly0da)

    call this%prefactor%initialize()

    this%init_status = 1
  end subroutine initialize_acc_nonlinearity

  subroutine finalize_acc_nonlinearity(this)
    class(df_acc_nonlinear_term_t) :: this

    call this%prefactor%finalize()

    call finalize_fourier

    this%init_status = 0
  end subroutine finalize_acc_nonlinearity

#undef DEBUGGING
  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_acc_nonlinearity_df(this,gy_chi,g2d,vexy,dgdxy,localrhs,first)
    class(df_acc_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi, g2d
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    !Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1,1:this%lbg0) :: tmp_arr
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:this%lbg0) ::  vexy_re, dgdxy_re
    !Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:this%lbg0) ::  nonlin1_re
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) ::  nonlin1_re
    Real:: ve_x_max_loc, ve_y_max_loc
    Integer :: i,j, klmn

    print*,"vexy = ",sum(abs(vexy(:,:,:,1:this%lbg0))),", dgdxy = ",sum(abs(dgdxy(:,:,:,1:this%lbg0)))

    ! anti-aliasing and y fourier transform
    Call gto_real(this,vexy(:,:,:,1:this%lbg0),vexy_re,2*this%lbg0)
    Call gto_real(this,dgdxy(:,:,:,1:this%lbg0),dgdxy_re,2*this%lbg0)


    print*,"vexy_re = ",sum(abs(vexy_re)),", dgdxy_re = ",sum(abs(dgdxy_re))
    ! get max ExB velocity
    if (first) then
       PERFON_I("ve_max")
       ve_x_max_loc=maxval(vexy_re(:,:,2,:))
       ve_y_max_loc=maxval(vexy_re(:,:,1,:))
       !write(*,"(A,2ES17.10)") "ve_max_loc = ",ve_x_max_loc, ve_y_max_loc
       ve_max(1)=max(ve_max(1),ve_x_max_loc)        
       ve_max(2)=max(ve_max(2),ve_y_max_loc)
       PERFOFF_I
    end if

    !compute the 'standard' nonlinear term
    do klmn=1,this%lbg0
       nonlin1_re = -vexy_re(:,:,1,klmn)*dgdxy_re(:,:,2,klmn) &
            & + vexy_re(:,:,2,klmn)*dgdxy_re(:,:,1,klmn)
       ! transform back to fourier space
       Call to_fourier_y(nonlin1_re,tmp_arr)
       nonlin(:,:,klmn) = transpose(tmp_arr)
    end do

    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0)

  End Subroutine calc_acc_nonlinearity_df

  !> Go to real space
  subroutine gto_real(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
    real, dimension(:,:), pointer :: p_out
    !real :: local_sum,thread_global_sum(0:7,0:n_procs_x*n_procs_y)
    Integer:: klmn, i_block,number_of_lbg0_blocks
    integer :: c2r_plan_y, err
    !Complex, Dimension(0:nj0-1,0:li0da/n_procs_y-1,1:howmany) :: dev_full_tmp2
    Complex, device, Dimension(:,:,:) :: dev_full_tmp2

    !$acc data copyin(inarr) copyout(outarr) create(dev_full_tmp2)
    !$acc kernels 
    !$acc loop independent
    do klmn=1,howmany
       dev_full_tmp2(:,:,klmn)=transpose(inarr(:,:,klmn))
    end do
    !$acc end kernels
    !$acc end data
    err = cufftPlan1d(c2r_plan_y,ly0da,CUFFT_Z2D,li0da*howmany);

    err = cufftExecZ2D(c2r_plan_y, dev_full_tmp2, outarr);

    !number_of_lbg0_blocks = howmany/this%lbg0
    !do i_block=1,number_of_lbg0_blocks
    !   do klmn=1,this%lbg0
          ! Transpose x-y
          !tmp_arr2 = transpose(inarr(:,:,(i_block-1)*this%lbg0+klmn))
          ! Fourier transfrom in y (include dealiasing step)
          !p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
    !      Call to_real_y(dev_full_tmp2(:,:,(i_block-1)*this%lbg0+klmn),outarr(:,:,(i_block-1)*this%lbg0+klmn))
    !   end do
    !end do
  end subroutine gto_real

end module df_acc_nonlinear_term_mod
