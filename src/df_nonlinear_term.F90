#include "redef.h"
#include "intrinsic_sizes.h"

module df_nonlinear_term_mod
  use nonlinear_term_mod
  use df_nonlinear_prefactor_mod
  use df_epc_nl_prefactor_mod

  use par_mod, only: imag, pi, ve_max, phase_1, phase_2, spec
  use par_other, only: n_fields
  use communications, only: MPI_COMPLEX_TYPE, MPI_REAL_TYPE, mpi_comm_y,&
       &my_2reals_max_to_all, my_barrier, mpi_comm_xy, reduce_per_thread,&
       &threadlocal_mpi_comm_y
  use discretization
  use fourier
  !use prefactors
  use mpi
#ifdef WITH_LIKWID
#include "likwid_f90.h"
#endif
  implicit none

  type, public, extends(nonlinear_term_t) :: df_nonlinear_term_t
     class(df_nonlinear_prefactor_t),allocatable,public :: prefactor
   contains
     procedure :: initialize => initialize_nonlinearity
     procedure :: finalize => finalize_nonlinearity
     procedure :: add => add_nonlinearity
     procedure :: calc => calc_nonlinearity_df
     procedure :: gto_real => gto_real_nonlinearity_df
     procedure :: mem_est => mem_est_nonlinearity
     procedure :: construct => construct_df_nonlinear_term
     !final :: destruct
     procedure :: destruct => destruct_df_nonlinear_term
     procedure :: getType => getThisType
     PROCEDURE :: get_ptr_to_ranshift_to => get_ptr_to_ranshift_to_df
     procedure :: get_ptr_to_ranshift_back => get_ptr_to_ranshift_back_df
  end type df_nonlinear_term_t

  public :: transpose_cmplx, transpose_and_extend_cmplx

  private

contains

  function getThisType(this)
    class(df_nonlinear_term_t) :: this
    character(len=MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_nonlinear_term_t"
  end function getThisType
    
  ! In the constructor we do work to instantiate the whole
  ! object, mainly for being able to call later on the
  ! mem_est routines.
  ! construct is not equal to initialize. construct just sets up
  ! the object, but does not allocate any data arrays.
  subroutine construct_df_nonlinear_term(this,equil_par_curr)
    class(df_nonlinear_term_t) :: this
    logical :: equil_par_curr

    !print*,"Constructing ",this%getType()

    if (equil_par_curr) then
       allocate(df_epc_nl_prefactor_t::this%prefactor)
    else
       allocate(df_nonlinear_prefactor_t::this%prefactor)
    end if
  end subroutine construct_df_nonlinear_term

  subroutine destruct_df_nonlinear_term(this)
    class(df_nonlinear_term_t) :: this
    !write(*,"(3A)",advance="no") "Destructing ",trim(this%getType()),", prefactor is "
    if (allocated(this%prefactor)) then
       !write(*,"(A)") "allocated. Now deallocating."
       deallocate(this%prefactor)
    !else
    !   write(*,"(A)") "not allocated."
    end if
  end subroutine destruct_df_nonlinear_term

  Real function mem_est_nonlinearity(this,mem_req_in)
    class(df_nonlinear_term_t) :: this
    real, intent(IN) :: mem_req_in
    real :: mem_loc
    
    mem_loc = 3.*9.*nj0*li0da/n_procs_y*SIZE_OF_REAL_MB

    mem_loc = mem_loc + (lij0+2.*nj0*li0da/n_procs_y+&
         &(li0da/n_procs_y+2.*li0da+2*nib)*3.*nj0)*SIZE_OF_COMPLEX_MB

    mem_est_nonlinearity = mem_req_in + this%prefactor%mem_est(mem_loc)
    
  end function mem_est_nonlinearity

  subroutine initialize_nonlinearity(this)
    class(df_nonlinear_term_t) :: this

    if (this%init_status==1) return
    !write(*,"(3A,I5)") "Initializing ",this%getType()," with blocksize = ",this%lbg0
    call initialize_fourier(li0da,ly0da)

    call this%prefactor%initialize()

    this%init_status = 1
  end subroutine initialize_nonlinearity

  subroutine finalize_nonlinearity(this)
    class(df_nonlinear_term_t) :: this

    call this%prefactor%finalize()

    call finalize_fourier

    this%init_status = 0
  end subroutine finalize_nonlinearity

  FUNCTION get_ptr_to_ranshift_to_df(this)
    class(df_nonlinear_term_t) :: this
    COMPLEX,DIMENSION(:),POINTER :: get_ptr_to_ranshift_to_df

    get_ptr_to_ranshift_to_df => this%prefactor%ranshift_to_df
  END FUNCTION get_ptr_to_ranshift_to_df

  FUNCTION get_ptr_to_ranshift_back_df(this)
    class(df_nonlinear_term_t) :: this
    COMPLEX,DIMENSION(:),pointer :: get_ptr_to_ranshift_back_df

    get_ptr_to_ranshift_back_df => this%prefactor%ranshift_back_df
  END FUNCTION get_ptr_to_ranshift_back_df
  
  subroutine add_nonlinearity(this,g_block,p_dgdxy,p_emfields,p_barchi,p_dbarchidxy,p_rhs,lb1,lb2,stage)
    class(df_nonlinear_term_t) :: this
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout), target:: p_rhs    
    integer, intent(in):: stage
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    Complex, Dimension(:,:,:), pointer :: p_barchi
    Complex, Dimension(:,:,:,:), pointer :: p_dbarchidxy, p_dgdxy
    complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in), target:: g_block

    ! Local variables
    real, dimension(2) :: ve_max_sav
    !complex, dimension(li1:li2, lj1:lj2,lb1:lb2):: p_g2d
    integer:: j
    logical:: first

    LIKWID_ON('add_nl')
    if (stage.eq.1) then
       first=.true.
       ve_max_sav = ve_max
    else
       first=.false.
    end if

    if (turbdeal) then
       if(lb1.eq.1) then
          if (stage.lt.3) then
             do j=lj1,lj2
                this%prefactor%ranshift_to_df(j)=exp(imag*2.*pi*j*phase_1)*this%prefactor%shiftmat_df(j,stage-1)
             enddo
          else
             do j=lj1,lj2
                this%prefactor%ranshift_to_df(j)=exp(imag*2.*pi*j*phase_1)*this%prefactor%shiftmat_df(j,stage-3)
             enddo
          endif
          this%prefactor%ranshift_back_df=conjg(this%prefactor%ranshift_to_df)
       endif
    endif

    ! add the ExB nonlinearity, calculation is performed in real space
    !p_g2d=g_block(li1:li2,:,:)

    PERFON('calc_nl')
    !print*,"size(p_dbarchidxy) = ",size(p_dbarchidxy,1),size(p_dbarchidxy,2),&
    !     &size(p_dbarchidxy,3),size(p_dbarchidxy,4)
    !write(*,"(2(A,ES17.10))") "vexy = ",real(sum(conjg(p_dbarchidxy)*p_dbarchidxy)),&
    !     &", dgdxy = ",real(sum(conjg(p_dgdxy)*p_dgdxy))
    !CALL this%calc(p_barchi,p_g2d,p_dbarchidxy,p_dgdxy,p_rhs,first)
    CALL this%calc(p_barchi,g_block,p_dbarchidxy,p_dgdxy,p_rhs,first,lb1,lb2)
    PERFOFF

    if (first) then
       ! now in ve_max we find the maximal ExB velocity in the whole 
       ! volume. The reduction on the phase space has already been 
       ! performed, but we have to do it still over the species.
       !print*,"ve_max = ",ve_max(1),ve_max(2)
       ve_max(1)=max(ve_max(1),ve_max_sav(1))
       ve_max(2)=max(ve_max(2),ve_max_sav(2))
       !exchange only at the end of stage 1
       if(lb2.eq.lklmn0) then
          Call my_2reals_max_to_all(ve_max)
       endif
    endif
    LIKWID_OFF('add_nl')
  end subroutine add_nonlinearity

#undef DEBUGGING
  !>Computes the nonlinearity for the (x) nonlocal version
  !!\todo Get rid of transposition using striding? make some arrays for arakawa allocatable
  !!\todo Check whether the nonlinearity prefactor is important for the CFL criterion
  Subroutine calc_nonlinearity_df(this,gy_chi,g_block,vexy,dgdxy,localrhs,first,lb1,lb2)
    class(df_nonlinear_term_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: gy_chi
    complex, dimension(lbi:ubi,lj1:lj2,1:*),intent(in) :: g_block
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first
    Integer, Intent(in) :: lb1,lb2

    ! Local variables
    Complex, Dimension(li1:li2,lj1:lj2,1:this%lbg0) :: nonlin    
    !Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1,1:this%lbg0) :: tmp_arr
    Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,2,1:this%lbg0) ::  vexy_re, dgdxy_re
    !Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:this%lbg0) ::  nonlin1_re
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1) ::  nonlin1_re
    Real, Dimension(1:2)  ::  vexy_max_loc
    Integer :: i,j, klmn
#if DEBUGGING
    real :: sum_of_part
    integer :: nXYPlanesPerPart, nParts, iPart, start_index, end_index
#endif
#ifdef WITHOMP
    integer :: my_thread, omp_get_thread_num
    my_thread = omp_get_thread_num()
#endif

    if (turbdeal) then
       do j=lj1,lj2
          vexy(:,j,:,1:this%lbg0)  = vexy(:,j,:,1:this%lbg0)* this%prefactor%ranshift_to_df(j)
          dgdxy(:,j,:,1:this%lbg0) = dgdxy(:,j,:,1:this%lbg0)* this%prefactor%ranshift_to_df(j)
       enddo
    endif

    ! anti-aliasing and y fourier transform
    PERFON('deal_FFT')
    Call this%gto_real(vexy,vexy_re,2*this%lbg0)
    !Call gto_real(ve_y,ve_y_re)
    Call this%gto_real(dgdxy,dgdxy_re,2*this%lbg0)
    !Call gto_real(dg_dy,dg_dy_re)
    PERFOFF

    ! get max ExB velocity
    if (first) then
       PERFON("ve_max")
#if 0
       print*,"First part:"
       do klmn=1,this%lbg0
          write(*,"(F10.6)",advance="no") maxval(vexy_re(:,:,1,klmn))
       end do
       write(*,"(A)") ""
       
       print*,"Second part:"
       do klmn=1,this%lbg0
          write(*,"(F10.6)",advance="NO") maxval(vexy_re(:,:,2,klmn))
       end do
       write(*,"(A)") ""
#endif
       call this%prefactor%multiply_max(vexy_max_loc,vexy_re,lb1,lb2)
       ve_max(1)=max(ve_max(1),vexy_max_loc(1)) !ve_x
       ve_max(2)=max(ve_max(2),vexy_max_loc(2)) !ve_y
       PERFOFF
    end if

    !compute the 'standard' nonlinear term
    PERFON('nonlin')
    !nonlin1_re = -vexy_re(:,:,1,:)*dgdxy_re(:,:,2,:) + vexy_re(:,:,2,:)*dgdxy_re(:,:,1,:)
#if DEBUGGING
       nParts = 4
       nXYPlanesPerPart = this%lbg0/nParts
       sum_of_part = 0.0D0
#endif
       do klmn=1,this%lbg0
          PERFON('nl_kernl')
          nonlin1_re = -vexy_re(:,:,1,klmn)*dgdxy_re(:,:,2,klmn) &
               & + vexy_re(:,:,2,klmn)*dgdxy_re(:,:,1,klmn)
          PERFOFF
          ! transform back to fourier space
          Call to_fourier_y(nonlin1_re,tmp_arr)
#if DEBUGGING
          sum_of_part = sum_of_part + sum(abs(tmp_arr))
          if (modulo(klmn,nXYPlanesPerPart).eq.0) then
             write(*,"(A,ES17.10)") "tmp_arr after to_fourier_y is ",sum_of_part
             sum_of_part = 0.0D0
          end if
#endif
          Call transpose_cmplx(nj0, li0da, tmp_arr, 0, nonlin(:,:,klmn), 0)
       end do
    PERFOFF

#if DEBUGGING
    do iPart=0,3
       start_index = iPart*this%lbg0/nParts+1
       end_index = (iPart+1)*this%lbg0/nParts
       !print*,start_index, end_index
       write(*,"(A,ES17.10)") "after backtranspose: nonlin = ",&
            & sum(abs(nonlin(:,:,start_index:end_index)))
    end do
#endif
    PERFON('pref_mul')
    ! multiply with the prefactor and write to localrhs
    call this%prefactor%multiply_with(nonlin,localrhs,this%lbg0,lb1,lb2)
    PERFOFF
#if DEBUGGING
    write(*,"(A,E17.10)") "localrhs after pref_mult is ",sum(abs(localrhs(:,:,1:this%lbg0)))
#endif
  End Subroutine calc_nonlinearity_df


  !> Go to real space
  subroutine gto_real_nonlinearity_df(this,inarr,outarr,howmany)
    class(df_nonlinear_term_t) :: this
    integer, intent(IN) :: howmany
    Complex, Dimension(li1:li2,lj1:lj2,1:howmany),Intent(in) :: inarr
    Real, Dimension(0:ly0da-1, 0:li0da/n_procs_y-1,1:howmany), target, Intent(out) :: outarr

    ! Local variables
    !Complex, Dimension(li1da:li2da,lj1:lj2,1:howmany) :: tmp_arr1
    !Complex, Dimension(0:nj0-1, 0:li0da/n_procs_y-1) :: tmp_arr2
    complex, dimension(0:ly0da/2,0:li0da/n_procs_y-1) :: extended_tarray

    real, dimension(:,:), pointer :: p_out
    !real :: local_sum,thread_global_sum(0:7,0:n_procs_x*n_procs_y)
    Integer:: klmn, i_block,number_of_lbg0_blocks

    number_of_lbg0_blocks = howmany/this%lbg0
    if (number_of_lbg0_blocks.eq.0) then
       ! Transpose x-y
       Call transpose_and_extend_cmplx(li0da, nj0, inarr(:,:,1),&
            &extended_tarray , ly0da/2+1)
       ! Fourier transfrom in y (include dealiasing step)
       !Call to_real_y(tmp_arr2,outarr(:,:,1))
       Call to_real_y_no_extend(extended_tarray,outarr(:,:,1))
    else
       do i_block=1,number_of_lbg0_blocks
          do klmn=1,this%lbg0
             !tmp_arr1(:,:,(i_block-1)*this%lbg0+klmn)=inarr(:,:,(i_block-1)*this%lbg0+klmn)
             ! Transpose x-y
             !Call transpose_cmplx(li0da, nj0, inarr(:,:,(i_block-1)*this%lbg0+klmn), 0,tmp_arr2 , 0)
             Call transpose_and_extend_cmplx(li0da, nj0, inarr(:,:,(i_block-1)*this%lbg0+klmn),&
                  &extended_tarray , ly0da/2+1)
             
             ! Fourier transfrom in y (include dealiasing step)
             !PERFON('to_re_y1')    
             !extended_tarray(0:nj0-1,:)=tmp_arr2
             ! Dealiasing in y direction
             !extended_tarray(nj0:,:)=cmplx(0.0,0.0)
             !PERFOFF
             p_out => outarr(:,:,(i_block-1)*this%lbg0+klmn)
             Call to_real_y_no_extend(extended_tarray,p_out) !outarr(:,:,(i_block-1)*this%lbg0+klmn))
          end do
       end do
    end if

  end subroutine gto_real_nonlinearity_df

  subroutine transpose_and_extend_cmplx(n1,n2,in_matrix,transposed_matrix,ldb)
    Integer, Intent(in) :: n1, n2, ldb
    Complex, Intent(in) :: in_matrix(0:, 0:)
    Complex, Intent(out):: transposed_matrix(0:, 0:)

    Complex, Dimension(n2/n_procs_y, n1/n_procs_y, 0:n_procs_y-1) :: sbuf, rbufc
    Integer:: n1l, n2l, i1
    Integer:: pp, ierr

    PERFON("TandExt")
    if (n_procs_y.eq.1) then
       transposed_matrix(0:n2-1,:) = transpose(in_matrix)
       transposed_matrix(n2:ldb-1,:)=cmplx(0.0,0.0)
    else
       n1l = n1/n_procs_y
       n2l = n2/n_procs_y
       Do pp = 0, n_procs_y-1
          i1 = pp*n1l
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
          i1 = pp*n2l
          transposed_matrix(i1:i1+n2l-1,:) = rbufc(:,:,pp)
       Enddo
       transposed_matrix(n2:ldb-1,:) = cmplx(0.0,0.0)
    end if
    PERFOFF
  end subroutine transpose_and_extend_cmplx

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

    PERFON("TransC")
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
    PERFOFF

  End Subroutine transpose_cmplx

end module df_nonlinear_term_mod
