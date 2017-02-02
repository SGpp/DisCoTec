#include "redef.h"
#ifdef WITHPAPI
#include "fpapi.h"
#endif

PROGRAM perf_MatrixMatrixMultiplication
  use MatrixModule
  USE BandedMatrixModule
  use VectorModule
  USE mpi
  IMPLICIT NONE

  !TYPE(Vector) :: vec, res
  INTEGER :: n_procs, rank, ierr, comm_cart
  !COMPLEX, DIMENSION(:,:), ALLOCATABLE :: localfullmat
  !COMPLEX, DIMENSION(:), allocatable :: localfullvec
#ifdef WITHPAPI
  integer :: papi_ncounters, check
  integer,parameter :: array_len=6
  !PAPI_TOT_CYC PAPI_L1_TCA PAPI_L1_TCM PAPI_TOT_INS
  !PAPI_TOT_CYC PAPI_L2_TCA PAPI_L2_TCM PAPI_TLB_TL PAPI_TLB_DM
  !integer,parameter,dimension(array_len) :: events = (/PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_L2_TCA, PAPI_L2_TCM/)
  !PAPI_L3_TCA PAPI_L3_TCM PAPI_L3_TCR
  !integer,parameter,dimension(array_len) :: events = (/PAPI_TOT_CYC, PAPI_L3_TCA, PAPI_L3_TCM, PAPI_L3_TCR/)
  !PAPI_TOT_CYC PAPI_TOT_INS PAPI_L3_TCA PAPI_L3_LDM PAPI_L3_DCW PAPI_L3_TCM
  !integer,parameter,dimension(array_len) :: events = (/PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_L3_TCA, PAPI_L3_LDM, PAPI_L3_DCW, PAPI_L3_TCM/)
  !PAPI_TOT_CYC PAPI_TOT_INS PAPI_L2_TCA PAPI_L2_LDM PAPI_L2_STM PAPI_L2_ICA
  integer,parameter,dimension(array_len) :: events = (/PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_L1_TCA, PAPI_L1_LDM, PAPI_L1_STM,PAPI_L1_ICA/)
  !PAPI_TOT_CYC PAPI_TOT_INS PAPI_VEC_DP PAPI_FP_OPS PAPI_RES_STL
  !integer,parameter,dimension(array_len) :: events = (/PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_VEC_DP, PAPI_FP_OPS, PAPI_RES_STL/)
  integer(kind=8),dimension(array_len) :: values

  call PAPIF_num_counters(papi_ncounters)
  print*,"We have ",papi_ncounters," PAPI counters."
#endif
  PERFINIT
  PERFON('all')

  call mpi_init(ierr)  
  CALL mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
  CALL mpi_cart_create(MPI_COMM_WORLD, 1, (/n_procs/), (/.FALSE./), .FALSE., comm_cart,ierr)

  call initialize_matrix_module(comm_cart)
  CALL initialize_BandedMatrix_module
  call initialize_Vector_module

  CALL mpi_comm_size(comm_cart,n_procs, ierr)
  CALL mpi_comm_rank(comm_cart,rank, ierr)
  write(*,"(A,I4,A)") "Using ",n_procs," processors."

  ! first test is FullMatrix x FullMatrix => FullMatrix
  !call test_FullFullFull
  ! second test is BandedMatrix(rmo) x FullMatrix => FullMatrix
  !call test_Banded_rmo_FullFull
  ! third test is TNT BandedMatrix(rmo) x BandedMatrix(cmo) =>BandedMatrix(rmo)
  !call test_Banded_rmo_Banded_cmo_Banded_rmo
  ! fourth test is TTT BandedMatrix(rmo) x BandedMatrix(rmo) => BandedMatrix(rmo)
  !call test_Banded_rmo_Banded_rmo_Banded_rmo
  !call test_NTT

  ! test the multiplication of a row-major matrix with a vector
  ! this is the operation in GENE used for gyro-averaging
  call test_TMatVec
  call mpi_finalize(ierr)
  PERFOFF
  PERFOUT('rmofufu')

contains
  subroutine test_FullFullFull
    Type(Matrix) :: Amat, Bmat,Cmat

    integer,parameter :: Ndims(8)=(/16,128,512,1024,2048,3072,4096,5120/)
    integer :: Ndim, irow,icol,i_loop
    real(8) :: time_start, time_end
    character(len=30) :: filename

    if (rank.eq.0) then
       write(filename,"(A,I3.3,A)") "fullfullfull_",n_procs,".dat"
       open(30,file=trim(filename))
    end if
    do i_loop=1,8
       Ndim=Ndims(i_loop)

       call initialize(Amat,Ndim,Ndim)
       call allocate(Amat)
       call initialize(Bmat,Ndim,Ndim)
       call allocate(Bmat)
       call initialize(Cmat,Ndim,Ndim)
       call allocate(Cmat)

       do irow=1,Ndim
          do icol=1,Ndim
             call set_value(Amat,irow,icol,cmplx(irow,icol))
             call set_value(Bmat,irow,icol,cmplx(irow,0.1*icol))
          end do
       end do
       call commit_values(Amat)
       call commit_values(Bmat)

       call mpi_barrier(comm_cart,ierr)
       time_start=MPI_Wtime()
       call dot_multiply(Amat,Bmat,Cmat,'T')
       call mpi_barrier(comm_cart,ierr)
       time_end  =MPI_Wtime()
       if (rank.eq.0) then
          write(*,"(A,I6,A,ES11.4,A)") "FullFullFull(",Ndim,") needed ",time_end-time_start," seconds."
          write(30,"(I6,2X,ES11.4)") Ndim,time_end-time_start
       end if
       call finalize(Amat)
       call finalize(Bmat)
       call finalize(Cmat)
    end do
    if (rank.eq.0) close(30)
  end subroutine test_FullFullFull

  subroutine test_Banded_rmo_FullFull
    type(BandedMatrix) :: Amat
    type(Matrix) :: Bmat,Cmat

    integer,parameter :: Ndims(8)=(/16,128,512,1024,2048,3072,4096,5120/)
    real, parameter :: bandratio(7)=(/0.001,0.01,0.02,0.05,0.1,0.2,0.5/)
    integer :: Ndim, irow,icol,i_loop, i_band,bandwidth
    real(8) :: time_start, time_end
    logical :: transpose=.true.
    character(len=30) :: filename
    EPIK_USER_REG(dotmul,"TFF_dot_mult")

    PERFON('TFF')

    if (rank.eq.0) then
       write(filename,"(A,I3.3,A)") "rmofullfull_",n_procs,".dat"
       open(30,file=trim(filename))
    end if
    do i_band=5,5
       do i_loop=7,7
          Ndim=Ndims(i_loop)

          call initialize(Amat,Ndim,Ndim,transpose)
          call allocate(Amat)
          call initialize(Bmat,Ndim,Ndim)
          call allocate(Bmat)
          call initialize(Cmat,Ndim,Ndim)
          call allocate(Cmat)

          bandwidth = Ndim*bandratio(i_band)/2
          do irow=1,Ndim
             do icol=1,Ndim
                if ((icol.gt.irow-bandwidth).and.(icol.lt.irow+bandwidth)) then
                   call set_value(Amat,irow,icol,cmplx(irow,icol))
                end if
                call set_value(Bmat,irow,icol,cmplx(irow,0.1*icol))
             end do
          end do
          call commit_values(Amat)
          call commit_values(Bmat)

          call mpi_barrier(comm_cart,ierr)
          PERFON('TFF_dm')
          time_start=MPI_Wtime()
          EPIK_USER_START(dotmul)
          call dot_multiply(Amat,Bmat,Cmat)
          EPIK_USER_END(dotmul)
          call mpi_barrier(comm_cart,ierr)
          time_end  =MPI_Wtime()
          PERFOFF
          if (rank.eq.0) then
             write(*,"(A,2I6,F6.3,A,ES11.4,A)") "rmoFullFull(",Ndim,bandwidth,bandratio(i_band),") needed ",time_end-time_start," seconds."
             write(30,"(2I6,2X,ES11.4)") Ndim,bandwidth,time_end-time_start
          end if
          call finalize(Amat)
          call finalize(Bmat)
          call finalize(Cmat)
       end do
       write(30,"(A)") ""
    end do
    if (rank.eq.0) close(30)
    PERFOFF
  end subroutine test_Banded_rmo_FullFull

  subroutine test_Banded_rmo_Banded_cmo_Banded_rmo
  end subroutine test_Banded_rmo_Banded_cmo_Banded_rmo

  subroutine test_Banded_rmo_Banded_rmo_Banded_rmo
    type(BandedMatrix) :: Amat,Bmat,Cmat

    integer,parameter :: Ndims(10)=(/16,128,512,1024,2048,3072,4096,5120,8192,16384/)
    real, parameter :: bandratio(7)=(/0.001,0.01,0.02,0.05,0.1,0.2,0.5/)
    integer :: Ndim, irow,icol,i_loop, i_band,bandwidth
    real(8) :: time_start, time_end
    logical :: transpose=.true.
    character(len=30) :: filename
    EPIK_USER_REG(dotmul,"TTT_dot_mult")

    PERFON('TTT')

    if (rank.eq.0) then
       write(filename,"(A,I3.3,A)") "rmormormo_",n_procs,".dat"
       open(30,file=trim(filename))
    end if
    do i_band=5,5
       do i_loop=7,7
          Ndim=Ndims(i_loop)

          call initialize(Amat,Ndim,Ndim,transpose)
          call allocate(Amat)
          call initialize(Bmat,Ndim,Ndim,transpose)
          call allocate(Bmat)
          call initialize(Cmat,Ndim,Ndim,transpose)
          call allocate(Cmat)

          bandwidth = Ndim*bandratio(i_band)/2
          do irow=1,Ndim
             do icol=1,Ndim
                if ((icol.gt.irow-bandwidth).and.(icol.lt.irow+bandwidth)) then
                   call set_value(Amat,irow,icol,cmplx(irow,icol))
                   call set_value(Bmat,irow,icol,cmplx(irow,0.1*icol))
                end if
             end do
          end do
          call commit_values(Amat)
          call commit_values(Bmat)

          call mpi_barrier(comm_cart,ierr)
          PERFON('TTT_dm')
          time_start=MPI_Wtime()
          EPIK_USER_START(dotmul)
          call dot_multiply(Amat,Bmat,Cmat)
          EPIK_USER_END(dotmul)
          call mpi_barrier(comm_cart,ierr)
          time_end  =MPI_Wtime()
          PERFOFF
          if (rank.eq.0) then
             write(*,"(A,2I6,F6.3,A,ES11.4,A)") "rmormormo(",Ndim,bandwidth,bandratio(i_band),") needed ",time_end-time_start," seconds."
             write(30,"(2I6,F6.3,2X,ES11.4)") Ndim,bandwidth,bandratio(i_band),time_end-time_start
          end if
          call finalize(Amat)
          call finalize(Bmat)
          call finalize(Cmat)
       end do
       write(30,"(A)") ""
    end do
    if (rank.eq.0) close(30)
    PERFOFF
  end subroutine test_Banded_rmo_Banded_rmo_Banded_rmo

  subroutine test_NTT
    type(BandedMatrix) :: Amat,Bmat,Cmat

    !integer,parameter :: Ndims(13)=(/16,128,512,1024,2048,3072,4096,5120,8192,16384,20000,32768,40000/)
    integer,allocatable, dimension(:) :: Ndims
    real, parameter :: bandratio(6)=(/0.001,0.01,0.02,0.05,0.1,0.2/)
    integer :: Ndim, irow,icol,i_loop, i_band,bandwidth, max_Ndim, min_Ndim, steps_Ndim,ival
    real(8) :: time_start, time_end
    real :: step_Ndim
#ifdef WITHPERF
    real :: ti,mf
#endif
    logical :: transpose=.true.
    character(len=30) :: filename
    EPIK_USER_REG(dotmul,"NTT_dot_mult")

    PERFON('NTT')

    if (rank.eq.0) then
       write(filename,"(A,I3.3,A)") "cmormormo_",n_procs,".dat"
       open(unit=30,file=trim(filename))
    end if
    max_Ndim = 2000
    min_Ndim = 100
    steps_Ndim = 100
    if (steps_Ndim.gt.1) then
       step_Ndim = (log(real(max_Ndim))-log(real(min_Ndim)))/real(steps_Ndim-1)
       !print*,(log(real(max_Ndim))-log(real(min_Ndim)))/(steps_Ndim-1)
       print*,"step_Ndim = ",step_Ndim
    else
       step_Ndim = 0.0
    end if
    allocate(Ndims(steps_Ndim))
    do i_loop=1,steps_Ndim
       Ndims(i_loop) = ceiling(exp(log(real(min_Ndim))+(i_loop-1)*step_Ndim))
    end do
    print*,"Ndims = ",Ndims
    do i_band=5,5
       do i_loop=1,steps_Ndim
          Ndim=Ndims(i_loop)

          call initialize(Amat,Ndim,Ndim)
          call allocate(Amat)
          call initialize(Bmat,Ndim,Ndim,transpose)
          call allocate(Bmat)
          call initialize(Cmat,Ndim,Ndim,transpose)
          call allocate(Cmat)

          bandwidth = Ndim*bandratio(i_band)/2
          do irow=1,Ndim
             do icol=1,Ndim
                if ((icol.gt.irow-bandwidth).and.(icol.lt.irow+bandwidth)) then
                   call set_value(Amat,irow,icol,cmplx(irow,icol))
                   call set_value(Bmat,irow,icol,cmplx(irow,0.1*icol))
                end if
             end do
          end do
          call commit_values(Amat)
          call commit_values(Bmat)

          call mpi_barrier(comm_cart,ierr)
          PERFON('NTT_dm')
          time_start=MPI_Wtime()
          EPIK_USER_START(dotmul)
#ifdef WITHPAPI
          call PAPIF_start_counters( events,  array_len,  check )
#endif
          call dot_multiply(Amat,Bmat,Cmat,'N','N')
#ifdef WITHPAPI
          call PAPIF_stop_counters( values,  array_len,  check )
#endif
          EPIK_USER_END(dotmul)
          call mpi_barrier(comm_cart,ierr)
          time_end  =MPI_Wtime()
          PERFOFF
          PERF_GET('NTT_dm',ti,mf)
!          PERF_RESET('NTT_dm')
          if (rank.eq.0) then
             write(*,"(A,2I6,F6.3,A,ES11.4,A)",advance='no') &
                  &"cmormormo(",Ndim,bandwidth,bandratio(i_band),") needed ",&
                  &time_end-time_start," seconds, "
             write(30,"(2I6,F6.3,2X,ES11.4)",advance='no') Ndim,bandwidth,bandratio(i_band),time_end-time_start
#ifdef WITHPERF
             write(*,"(A,F8.4,A,F8.1)") "perftime = ",ti,", MFlop/s = ",mf
             write(30,"(F8.4,F8.1)") ti, mf
#endif
#ifdef WITHPAPI
             do ival=1,array_len
                write(*,"(I20)",advance='no') values(ival)
                write(30,"(I20)",advance='no') values(ival)
             end do
             write(*,"(A)") ""
             write(30,"(A)") ""
#endif
          end if
          call finalize(Amat)
          call finalize(Bmat)
          call finalize(Cmat)
       end do
       write(30,"(A)") ""
    end do
    if (rank.eq.0) close(30)
    PERFOFF
  end subroutine test_NTT

  subroutine test_TMatVec
    type(BandedMatrix) :: Amat
    type(Vector) :: Vec, Res

    !integer,parameter :: Ndims(13)=(/16,128,512,1024,2048,3072,4096,5120,8192,16384,20000,32768,40000/)
    integer,allocatable, dimension(:) :: Ndims
    real, parameter :: bandratio(6)=(/0.001,0.01,0.02,0.05,0.1,0.2/)
    integer :: Ndim, irow,icol,i_loop, i_band,bandwidth, max_Ndim, min_Ndim, steps_Ndim,ival
    real(8) :: time_start, time_end
    real :: step_Ndim
#ifdef WITHPERF
    real(8) :: ti,mf
#endif
    logical :: transpose=.true.
    character(len=30) :: filename
    EPIK_USER_REG(dotmul,"NTT_dot_mult")

    PERFON('NTT')

    if (rank.eq.0) then
       write(filename,"(A,I3.3,A)") "TMatVec_",n_procs,".dat"
       open(unit=30,file=trim(filename))
    end if
    min_Ndim = 50
    max_Ndim = 10000
    steps_Ndim = 23
    if (steps_Ndim.gt.1) then
       step_Ndim = (log(real(max_Ndim))-log(real(min_Ndim)))/real(steps_Ndim-1)
       !print*,(log(real(max_Ndim))-log(real(min_Ndim)))/(steps_Ndim-1)
       print*,"step_Ndim = ",step_Ndim
    else
       step_Ndim = 0.0
    end if
    allocate(Ndims(steps_Ndim))
    do i_loop=1,steps_Ndim
       Ndims(i_loop) = ceiling(exp(log(real(min_Ndim))+(i_loop-1)*step_Ndim))
    end do
    print*,"Ndims = ",Ndims
    do i_band=5,5
       do i_loop=1,steps_Ndim
          Ndim=Ndims(i_loop)

          call initialize(Amat,Ndim,Ndim,transpose)
          call allocate(Amat)
          call initialize(Vec,Ndim)
          call allocate(Vec)
          call initialize(Res,Ndim)
          call allocate(Res)

          bandwidth = Ndim*bandratio(i_band)/2
          do irow=1,Ndim
             do icol=irow-bandwidth,irow+bandwidth
                !if ((icol.gt.irow-bandwidth).and.(icol.lt.irow+bandwidth)) then
                call set_value(Amat,irow,icol,cmplx(irow,icol))
                !end if
             end do
             call commit_values(Amat)
          end do
          call commit_values(Amat)

          !call autotune(Amat,Vec,Res)
          
          call mpi_barrier(comm_cart,ierr)
          PERFON('TMV_dm')
          time_start=MPI_Wtime()
          EPIK_USER_START(dotmul)
#ifdef WITHPAPI
          call PAPIF_start_counters( events,  array_len,  check )
#endif
          call dot_multiply(Amat,Vec,Res)
#ifdef WITHPAPI
          call PAPIF_stop_counters( values,  array_len,  check )
#endif
          EPIK_USER_END(dotmul)
          call mpi_barrier(comm_cart,ierr)
          time_end  =MPI_Wtime()
          PERFOFF
          !PERF_GET('TMV_dm',ti,mf)
          !PERF_RESET('TMV_dm')
          PERF_GET('sbm_dm',ti,mf)
          PERF_RESET('sbm_dm')
          !PERF_GET('subtask',ti,mf)
          !PERF_RESET('subtask')
          if (rank.eq.0) then
             write(*,"(A,2I6,F6.3,A,ES11.4,A)",advance='no') &
                  &"TMatVec(",Ndim,bandwidth,bandratio(i_band),") needed ",&
                  &time_end-time_start," seconds"
             write(30,"(2I6,F6.3,2X,ES11.4)",advance='no') Ndim,bandwidth,bandratio(i_band),time_end-time_start
#ifdef WITHPERF
             write(*,"(A,F8.4,A,F8.1)",advance='no') ", perftime = ",ti,", MFlop/s = ",mf
             write(30,"(F8.4,F8.1)",advance='no') ti, mf
#endif
#ifdef WITHPAPI
             do ival=1,array_len
                write(*,"(I20)",advance='no') values(ival)
                write(30,"(I20)",advance='no') values(ival)
             end do
#endif
             write(*,"(A)") ""
             write(30,"(A)") ""
          end if
          call finalize(Amat)
          call finalize(Vec)
          call finalize(Res)
       end do
       write(30,"(A)") ""
    end do
    if (rank.eq.0) close(30)
    PERFOFF
  end subroutine test_TMatVec
end PROGRAM perf_MatrixMatrixMultiplication
