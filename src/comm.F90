#include "redef.h"
!>Module containing MPI communicators and several wrappers
Module communications
  Use par_mod
  Use mpi
  USE MatrixModule
  use VectorModule
  USE, INTRINSIC :: iso_c_binding

  Implicit None
  
  public:: initialize_comm_scan, initialize_comm_sim, split_comm, free_comm, finalize_comm_sim, &
       finalize_comm_scan, set_my_mpi_comm_world, omp_num_threads, omp_level
  public:: my_barrier, calculate_test_sum, get_systime, my_secnds, my_real_get, pexyzvwspec, &
       my_complex_sum_w, my_real_sum_w, my_complex_sum_vw, my_complex_sum_vwspec, my_complex_sum_spec,&
       sum_vw, sum_to_all_complex,  my_real_max_to_all, my_sum_to_all_real, my_2reals_max_to_all, &
       my_real_min_to_all, my_sum_to_0_real, my_sum_to_0_complex, my_real_sum_vw, my_real_sum_spec,&
       &my_complex_sum_wspec, my_complex_sum_v, my_sum_to_all_int
  public::my_real_gather_to_0,my_complex_gather_to_0,t_n_procs
  public::my_real_sum_w_0d, my_sum_to_all_real_0d, my_real_bcast
  public:: mpi_comm_x, mpi_comm_y, mpi_comm_z,  mpi_comm_yz, mpi_comm_xyz, mpi_comm_v,&
       mpi_comm_w, mpi_comm_vw, mpi_comm_spec, mpi_comm_vwspec, mpi_comm_xyzvw, mpi_comm_xy_for_x, &
       mpi_comm_zvw, mpi_comm_xyzvwspec, mpi_comm_world, mpi_comm_xyzw, mpi_comm_xy,mpi_comm_zvwspec
  public:: threadlocal_mpi_comm_y
  public:: communicators, comm_x, comm_z, comm_v, comm_w, comm_cart, n_pe_dims
  public:: MPI_REAL_TYPE, MPI_COMPLEX_TYPE, mpi_integer, mpi_logical, mpi_sum, &
       & mpi_min, mpi_max, mpi_land, mpi_lor, mpi_status_size, mpi_status_ignore, &
       & MY_MPI_COMM_WORLD,MPI_ANY_SOURCE, MPI_PROC_NULL, MPI_IN_PLACE
  public:: read_parall_nml !public for interfaces
  public :: reduce_per_thread

#if defined(AIX) || defined(BGP)
  public:: flush
#endif

  public :: sync_time

  private

  INTEGER, PARAMETER :: N_COMM=12 !< number of communicators, used in the code
  INTEGER, PARAMETER :: N_PE_DIMS=6 !< number of dimensions for the Cartesian topology

#ifdef FORTRAN2003_ENUM
  enum,bind(c)
  enumerator :: COMM_SPEC=1, COMM_W, COMM_V, COMM_Z, COMM_Y, COMM_X
  enumerator :: COMM_WSPEC, COMM_VW, COMM_YZ
  enumerator :: COMM_VWSPEC, COMM_XYZ, COMM_XYZVW
  endenum
#else
  INTEGER, parameter :: COMM_SPEC=1, COMM_W=2, COMM_V=3, COMM_Z=4, COMM_Y=5, COMM_X=6
  INTEGER, PARAMETER :: COMM_WSPEC=7, COMM_VW=8, COMM_YZ=9, COMM_VWSPEC=10, COMM_XYZ=11, COMM_XYZVW=12
#endif

  Integer:: mpi_comm_x, mpi_comm_y, mpi_comm_z,&
       &mpi_comm_v, mpi_comm_w, mpi_comm_w_spec, mpi_comm_zvw
  integer :: threadlocal_mpi_comm_y
#ifdef WITHOMP_BLOCKLOOP
  !$OMP THREADPRIVATE(threadlocal_mpi_comm_y)
#endif
  INTEGER:: mpi_comm_vw, mpi_comm_vwspec, mpi_comm_xy,mpi_comm_xy_for_x
  INTEGER:: mpi_comm_xyz, mpi_comm_yz, mpi_comm_spec, mpi_comm_xyzvw !,mpi_comm_yzs
  INTEGER:: mpi_comm_xyzvwspec, mpi_comm_xyzw, mpi_comm_zvwspec
  INTEGER, DIMENSION(N_COMM) :: communicators

  !INTEGER :: my_cart_rank
  INTEGER:: MY_MPI_COMM_WORLD, comm_cart
  Logical:: wtime_is_global
  Integer:: t_n_procs
  integer:: omp_num_threads, omp_level


  INTERFACE my_sum_w
     module procedure my_complex_sum_w, my_real_sum_w, my_real_sum_w_0d
  END INTERFACE

  INTERFACE calculate_test_sum
     MODULE PROCEDURE calculate_test_sum_4D,calculate_test_sum_6D,calculate_test_sum_2D,calculate_test_sum_3D
     module procedure calculate_test_sum_5D, calculate_test_sum_real_5D
  END INTERFACE

Contains
#if defined(AIX) || defined(BGP)
  Subroutine flush (unit)
    Integer:: unit
    Call flush_(unit)
  End Subroutine flush
#endif


  !>splits MPI_COMM_WORLD into n_parallel_sims and 
  !!n_procs_sim communicators
  !!
  subroutine initialize_comm_scan(split_comm,parall_comm)
    integer:: split_comm, parall_comm
    integer:: n_procs, rank
    integer:: ierr
    
    call mpi_comm_size (mpi_comm_world, n_procs, ierr)
    if (ierr /= 0) Stop 'mpi_comm_size failed!'
    call mpi_comm_rank (mpi_comm_world, mype_gl, ierr)

    if(mype_gl.eq.0) then
       select case(omp_level)
       case (MPI_THREAD_SINGLE)
          print*,"MPI_THREAD_SINGLE initialized."
       case (MPI_THREAD_FUNNELED)
          print*,"MPI_THREAD_FUNNELED initialized."
       case (MPI_THREAD_SERIALIZED)
          print*,"MPI_THREAD_SERIALIZED initialized."
       case (MPI_THREAD_MULTIPLE)
          print*,"MPI_THREAD_MULTIPLE initialized."
       end select
    end if

    my_sim=mype_gl/n_procs_sim
    rank=mype_gl-my_sim*n_procs_sim

    call mpi_comm_split(mpi_comm_world,my_sim,rank,&
         split_comm,ierr)
    call mpi_comm_rank (split_comm, mype, ierr)

    call mpi_comm_split(mpi_comm_world,rank,my_sim,&
         parall_comm,ierr)
    
  end subroutine initialize_comm_scan

  subroutine finalize_comm_scan(split_comm,parall_comm)
    integer:: split_comm, parall_comm, ierr

    call mpi_comm_free(split_comm,ierr)
    if (ierr /= 0) then
       print*,mype,": error with comm_free(split_comm)"
    end if

    call mpi_comm_free(parall_comm,ierr)
    if (ierr /= 0) then
       print*,mype,": error with comm_free(parall_comm)"
    end if
    
  end subroutine finalize_comm_scan


!>Initializes MY_MPI_COMM_WORLD which is the (sub) MPI_COMM_WORLD
!!which is employed by the current GENE instance.
!!
!! \param MPI_COMM_WORLD_in MPI_COMM_WORLD for the current GENE instance. 
!! Set to -1 for taking the whole environment supplied by mpi_init(_thread)
!!\todo y parallelization and OMP should work if MPI_THREAD_MULTIPLE is provided
!!\todo y parallelization and OMP do work for linear runs (move the check and parameters read to parameters_IO)
  Subroutine initialize_comm_sim(MPI_COMM_SIM)
    Integer,intent(in) :: MPI_COMM_SIM
    Integer:: ierr
    INTEGER:: n_procs_already_set
    Integer(KIND=MPI_ADDRESS_KIND):: attrval
    Logical::   attrflag

    ! (re-)Set mype as rank of the process (not yet set for trinity or itm)
    ! and t_n_procs as the total number of processes used

    Call mpi_comm_rank (MPI_COMM_SIM, mype, ierr)
    If (ierr /= 0) Stop 'mpi_comm_rank failed!'
    Call mpi_comm_size (MPI_COMM_SIM, t_n_procs, ierr)
    If (ierr /= 0) Stop 'mpi_comm_size failed!'
    omp_num_threads = check_omp()
    IF ((mype.EQ.0).and.(print_ini_msg)) then
#ifdef WITHOMP
       write(*,"(A,I7,A,I2,A,I7,A)") "We have ",t_n_procs," MPI tasks, each with ",&
            &omp_num_threads," OpenMP threads, thus using ",t_n_procs*omp_num_threads," cores."
#else
       write(*,"(A,I7,A)") "We have ",t_n_procs," MPI tasks."
#endif
    end IF

!    if ((omp_level.ne.MPI_THREAD_MULTIPLE).and.(omp_num_threads.gt.1)) then
!    if(omp_num_threads.gt.1) then
!       n_procs_y=1
!       if (mype.eq.0) PRINT*,"Set n_procs_y=1 because MPI_THREAD_MULTIPLE is not supported"
!       if (mype.eq.0) print*, 'OMP and y parallelization cannnot be combined, set n_procs_y to 1'
!    endif

    call read_parall_nml(par_in_dir)

    IF (n_procs_sim.LE.0) n_procs_sim=t_n_procs

    If (t_n_procs.lt.n_procs_sim)&
         stop 'The number of processes specified in the parameters file is greater than the number of mpi processes available!'
    
    If ((n_procs_sim.gt.0).and.(t_n_procs.gt.n_procs_sim).and.(mype.eq.0)) then
       WRITE(*,*)
       WRITE(*,"(3A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       WRITE(*,"(A,I3,A)") '!!! too many CPUs requested,',t_n_procs-n_procs_sim,' processors will idle !!!'  
       WRITE(*,"(3A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    endif

    if(n_procs_sim.eq.-1) n_procs_sim = t_n_procs

    ! Determine the number of processes already set by the parameters file
    n_procs_already_set=1
    call multiply_if_not_zero(n_procs_already_set,n_procs_s)
    call multiply_if_not_zero(n_procs_already_set,n_procs_w)
    call multiply_if_not_zero(n_procs_already_set,n_procs_v)
    call multiply_if_not_zero(n_procs_already_set,n_procs_x)
    call multiply_if_not_zero(n_procs_already_set,n_procs_y)
    call multiply_if_not_zero(n_procs_already_set,n_procs_z)
    
    ! check if any of the processnumbers are not set to a finite value. 
    IF (  (n_procs_s.LT.1).OR. &
         & (n_procs_w.LT.1).OR. &
         & (n_procs_v.LT.1).OR. &
         & (n_procs_x.LT.1).OR. &
         & (n_procs_y.LT.1).OR. &
         & (n_procs_z.LT.1) ) THEN

        ! if there are processes left, we set the autoparallelization flag
       IF (n_procs_already_set.LT.n_procs_sim) THEN
          auto_parall=.TRUE.
       elseif (n_procs_already_set.EQ.n_procs_sim) THEN
          ! if no processes left, we set the remaining number of processes to 1
          auto_parall=.FALSE.
          if (n_procs_s.lt.1) n_procs_s=1
          if (n_procs_w.lt.1) n_procs_w=1
          if (n_procs_v.lt.1) n_procs_v=1
          if (n_procs_x.lt.1) n_procs_x=1
          if (n_procs_y.lt.1) n_procs_y=1
          if (n_procs_z.lt.1) n_procs_z=1
       else
          stop 'The number of processes specified in the parameters file exceeds the number of available mpi processes!'
       endif
    else
       if (n_procs_already_set.ne.n_procs_sim) then
          stop "chosen parallelization is not consistent with the total number of processes"
       else
          auto_parall=.FALSE.
       end if
    endif
    
    call set_my_mpi_comm_world(MPI_COMM_SIM)
    
    ! is timer globally synchronized?
    ! The attribute MPI_WTIME_IS_GLOBAL is attached to MPI_COMM_WORLD.
    ! See http://www.mpi-forum.org/docs/mpi21-report-bw/node176.htm#Node176
    call mpi_comm_get_attr(MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL,&
         attrval,attrflag,ierr)
    
    wtime_is_global = (attrval.EQ.1)
  End Subroutine initialize_comm_sim

  subroutine set_my_mpi_comm_world(temporary_mpi_comm_world)
    INTEGER, intent(IN) :: temporary_mpi_comm_world

    ! Local variables
    integer:: active_cpu, ierr

    if (mype.lt.n_procs_sim) then
       my_mpi_comm_world_member=.true.
    else
       my_mpi_comm_world_member=.false.
    endif

    active_cpu=1
    if (.not.my_mpi_comm_world_member) active_cpu=MPI_UNDEFINED

    call mpi_comm_split(temporary_mpi_comm_world,active_cpu,mype,&
         &MY_MPI_COMM_WORLD,ierr)

    If (ierr /= 0) Stop 'my_mpi_comm_world split failed!'
  end subroutine set_my_mpi_comm_world

  subroutine multiply_if_not_zero(npset,nprocs)
    integer,intent(inout):: npset
    integer,intent(in):: nprocs 
    
    if(nprocs.gt.0) then 
       npset=npset*nprocs
    endif
  end subroutine multiply_if_not_zero
  
  Subroutine split_comm

    ! Local variables
    Integer:: ierr
    INTEGER :: my_pe_coords(N_PE_DIMS), pe_dims(N_PE_DIMS)
    LOGICAL :: periods(N_PE_DIMS), remains_dim(N_PE_DIMS)
    integer :: i_pe_dim

    if (my_mpi_comm_world_member) then

       ! This is the attempt to setup a Cartesian topology of the processes.
       pe_dims = (/n_procs_s,n_procs_w,n_procs_v,n_procs_z,n_procs_y,n_procs_x/)
       IF (rad_bc_type.NE.0) THEN
          ! x direction is not periodic
          periods = (/.FALSE.,.FALSE.,.FALSE.,.TRUE.,.FALSE.,.FALSE./)
       ELSE
          ! x direction is periodic
          periods = (/.FALSE.,.FALSE.,.FALSE.,.TRUE.,.FALSE.,.TRUE./)
       END IF

       CALL mpi_cart_create(MY_MPI_COMM_WORLD, N_PE_DIMS, pe_dims, periods, .false., comm_cart, ierr)
       if (ierr /= 0) stop 'error during mpi_cart_create'

       ! Create the communicators of the 1D subgrids (the directions)
       DO i_pe_dim=1,N_PE_DIMS
          remains_dim=.FALSE.
          remains_dim(i_pe_dim)=.TRUE.
          CALL mpi_cart_sub(comm_cart,remains_dim,communicators(i_pe_dim),ierr)
       END DO

       mpi_comm_spec = communicators(COMM_SPEC)
       mpi_comm_w    = communicators(COMM_W)
       mpi_comm_v    = communicators(COMM_V)
       mpi_comm_z    = communicators(COMM_Z)
       mpi_comm_y    = communicators(COMM_Y)
       mpi_comm_x    = communicators(COMM_X)

       ! 2D subgrids
       remains_dim=.FALSE.
       remains_dim(COMM_SPEC)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_w_spec,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_V)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_vw,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_Y)=.TRUE.
       remains_dim(COMM_Z)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_yz,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_X)=.TRUE.
       remains_dim(COMM_Y)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_xy,ierr)

       ! 3D subgrids
       remains_dim=.FALSE.
       remains_dim(COMM_Z)=.TRUE.
       remains_dim(COMM_V)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_zvw,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_SPEC)=.TRUE.
       remains_dim(COMM_V)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_vwspec,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_X)=.TRUE.
       remains_dim(COMM_Y)=.TRUE.
       remains_dim(COMM_Z)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_xyz,ierr)
       if (ierr.ne.0) stop 'error in mpi_cart_sub_xyz'

       ! 4D subgrid
       remains_dim=.FALSE.
       remains_dim(COMM_X)=.TRUE.
       remains_dim(COMM_Y)=.TRUE.
       remains_dim(COMM_Z)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_xyzw,ierr)
       remains_dim=.FALSE.
       remains_dim(COMM_Z)=.TRUE.
       remains_dim(COMM_V)=.TRUE.
       remains_dim(COMM_W)=.TRUE.
       remains_dim(COMM_SPEC)=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_zvwspec,ierr)

       ! 5D subgrid
       remains_dim=.TRUE.
       remains_dim(COMM_SPEC)=.FALSE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_xyzvw,ierr)
       if (ierr.ne.0) stop 'error in mpi_cart_sub_xyzvw'

       CALL mpi_comm_rank(comm_cart,mype,ierr)
       CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,my_pe_coords,ierr)

       ! 6D subgrids (for checkpoint HDF5 output file)
       remains_dim=.TRUE.
       CALL mpi_cart_sub(comm_cart,remains_dim,mpi_comm_xyzvwspec,ierr)
       if (ierr.ne.0) stop 'error in mpi_cart_sub_xyzvwspec'

       my_pespec = my_pe_coords(1)
       my_pew    = my_pe_coords(2)
       my_pev    = my_pe_coords(3)
       my_pez    = my_pe_coords(4)
       my_pey    = my_pe_coords(5)
       my_pex    = my_pe_coords(6)

       ! setup a communicator in the x-y plane with reordered ranks, to be used in the 
       ! Arakawa scheme, where we have the whole x dimension parallelized over
       ! n_procs_x*n_procs_y ranks in a special order.
       CALL mpi_comm_split(mpi_comm_xy,1,my_pex*n_procs_y+my_pey,mpi_comm_xy_for_x,ierr)

       ! initialize the process grid for the BLACS and ScaLapack routines
       IF (.NOT.xy_local) THEN
          CALL initialize_matrix_module(mpi_comm_x)
          call initialize_vector_module()
       END IF

    endif
  End Subroutine split_comm


  Subroutine free_comm
    Integer:: ierr

    if (my_mpi_comm_world_member) then
       if (.not.xy_local) THEN
          CALL finalize_matrix_module()
          call finalize_vector_module()
       END IF
       Call mpi_comm_free (mpi_comm_w,ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_w failed!'

       Call mpi_comm_free (mpi_comm_vw,ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_vw failed!'

       Call mpi_comm_free (mpi_comm_w_spec, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_w_spec failed!'

       Call mpi_comm_free (mpi_comm_zvw, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_zvw failed!'

       Call mpi_comm_free (mpi_comm_vwspec, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_vwspec failed!'

       Call mpi_comm_free (mpi_comm_v, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_v failed!'

       Call mpi_comm_free (mpi_comm_z, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_z failed!'

       Call mpi_comm_free (mpi_comm_y, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_y failed!'

       Call mpi_comm_free (mpi_comm_x, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_x failed!'

       Call mpi_comm_free (mpi_comm_yz, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_yz failed!'

       Call mpi_comm_free (mpi_comm_xy, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_xy failed!'

       Call mpi_comm_free (mpi_comm_xyz, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_xyz failed!'

       Call mpi_comm_free (mpi_comm_xyzw, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_xyzw failed!'

       Call mpi_comm_free (mpi_comm_zvwspec, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_zvwspec failed!'

       Call mpi_comm_free (mpi_comm_xyzvw, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_xyzvw failed!'

       Call mpi_comm_free (mpi_comm_xy_for_x, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_xy_for_x failed!'

       !CALL mpi_comm_free (mpi_comm_yzs, ierr)
       !If (ierr /= 0) Stop 'mpi_comm_free_yzs failed!'

       Call mpi_comm_free (mpi_comm_spec, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_spec failed!'

       Call mpi_comm_free (comm_cart, ierr)
       If (ierr /= 0) Stop 'mpi_comm_free_comm_cart failed!'

       Call mpi_comm_free (mpi_comm_xyzvwspec, ierr) 
       If (ierr /= 0) Stop 'mpi_comm_free_xyzvwspec failed!'
    endif
  End Subroutine free_comm


  ! ===================================================
  ! == MPI WRAPPER  ===================================
  ! ===================================================

  Subroutine my_barrier
    integer :: ierr
    Call mpi_barrier (MY_MPI_COMM_WORLD, ierr)
  End Subroutine my_barrier
  !
  !------------------ my_real_max_to_all --------------------
  !
  Subroutine my_real_max_to_all(areal)
    Real, Intent(Inout) :: areal
    Real :: areal2
    integer :: ierr
    Call mpi_allreduce (areal, areal2, 1, MPI_REAL_TYPE, MPI_MAX,&
         MY_MPI_COMM_WORLD, ierr)
    areal = areal2
  End Subroutine my_real_max_to_all

  SUBROUTINE my_2reals_max_to_all(areal)
    REAL(C_REAL_TYPE), DIMENSION(2), INTENT(Inout) :: areal
    REAL(C_REAL_TYPE), DIMENSION(2) :: areal2
    integer :: ierr

    CALL mpi_allreduce (areal, areal2, 2, MPI_REAL_TYPE, MPI_MAX,&
         MY_MPI_COMM_WORLD, ierr)
    areal = areal2
  END SUBROUTINE my_2reals_max_to_all

  Subroutine my_real_min_to_all(areal)
    Real, Intent(Inout) :: areal
    Real :: areal2
    integer :: ierr
    Call mpi_allreduce (areal, areal2, 1, MPI_REAL_TYPE, MPI_MIN,&
         MY_MPI_COMM_WORLD, ierr)
    areal = areal2
  End Subroutine my_real_min_to_all

  subroutine my_real_bcast(areal,pe)
    real, intent(inout):: areal
    integer, intent(in):: pe
    integer:: ierr
    call mpi_bcast(areal,1,MPI_DOUBLE_PRECISION,&
         &pe,MY_MPI_COMM_WORLD,ierr)

  end subroutine my_real_bcast

  !********************************************************!
  !************ Subroutine for Real Data Transfer *********!
  !********************************************************!
  Subroutine my_real_get(dest, src, len, pedest, pesrc)
    Real::    dest(*), src(*)
    Integer:: len, pedest, pesrc, ierr
    Integer:: stat(MPI_STATUS_SIZE)

    If (pedest == pesrc) Then
       If (mype == pedest) dest(1:len) = src(1:len)
       Return
    Endif
    If (mype == pesrc) Then
       !       write(*,"(I2,A)") mype,": senden"
       Call mpi_send (src, len, MPI_REAL_TYPE, pedest,&
            888, MY_MPI_COMM_WORLD, ierr)
       !       write(*,"(I2,A)") mype,": senden ENDE"
    Endif
    If (mype == pedest) Then
       !       write(*,"(I2,A)") mype,": empfangen"
       Call mpi_recv (dest, len, MPI_REAL_TYPE, pesrc,&
            888, MY_MPI_COMM_WORLD, stat, ierr)
       !       write(*,"(I2,A)") mype,": empfangen ENDE"
    Endif
  End Subroutine my_real_get

  !********************************************************!
  !********* Subroutine for Complex Data Transfer *********!
  !********************************************************!
  Subroutine my_complex_get(dest, src, len, pedest, pesrc)
    complex:: dest(*), src(*)
    Integer:: len, pedest, pesrc, ierr
    Integer:: stat(MPI_STATUS_SIZE)

    If (pedest == pesrc) Then
       If (mype == pedest) dest(1:len) = src(1:len)
       Return
    Endif
    If (mype == pesrc) Then
       Call mpi_send(src, len, MPI_COMPLEX_TYPE, pedest,&
            888, MY_MPI_COMM_WORLD, ierr)
    Endif
    If (mype == pedest) Then
       Call mpi_recv(dest, len, MPI_COMPLEX_TYPE, pesrc,&
            888, MY_MPI_COMM_WORLD, stat, ierr)
    Endif
  End Subroutine my_complex_get

  !**************************************************************!
  ! Subroutine for Complex Sum over n_procs_v, n_procs_w and specs
  !**************************************************************!
  Subroutine my_complex_sum_vwspec(localsum,len)
    Integer, Intent(In):: len
    Complex, Intent(Inout):: localsum(len)
    Complex               :: totalsum(len)
    Integer:: ierr

    PERFON("sum_vwsp")
    Call mpi_allreduce( &
         localsum(1), totalsum(1), len, &
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_vwspec, ierr)
    localsum = totalsum
    PERFOFF
  End Subroutine my_complex_sum_vwspec

  Subroutine my_complex_sum_wspec(localsum,len)
    Integer, Intent(In)    :: len
    Complex, Intent(Inout) :: localsum(len)
    Complex                :: totalsum(len)
    Integer :: ierr
    PERFON("sum_wsp")
    Call mpi_allreduce( &
         localsum(1),totalsum(1), len, &
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_w_spec, ierr)
    localsum = totalsum
    PERFOFF
  End Subroutine my_complex_sum_wspec


  SUBROUTINE my_complex_sum_spec(localsum,len)
    Integer, Intent(In):: len
    Complex, Intent(Inout):: localsum(len)
    Complex               :: totalsum(len)
    Integer:: ierr

    PERFON("sum_sp")
    Call mpi_allreduce(&
         localsum(1), totalsum(1), len,&
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_spec, ierr)
    localsum = totalsum
    PERFOFF
  END SUBROUTINE my_complex_sum_spec

  Subroutine my_complex_sum_vw(localsum, len)
    Integer, Intent(In)::    len
    complex, Intent(Inout):: localsum(len)
    complex::  totalsum(len)
    integer :: ierr

    ! Sum over all local sums of a certain species
    Call mpi_allreduce (localsum(1), totalsum(1), len,&
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_vw, ierr)
    localsum = totalsum
  End Subroutine my_complex_sum_vw

  Subroutine my_complex_sum_v(localsum,len)
    Integer, Intent(In)    :: len
    Complex, Intent(InOut) :: localsum(len)
    Complex                :: totalsum(len)
    Integer :: ierr
    Call mpi_allreduce (localsum(1), totalsum(1), len,&
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_v, ierr)
    localsum = totalsum
  End Subroutine my_complex_sum_v

  Subroutine my_complex_sum_w(localsum,len)
    Integer, Intent(In)    :: len
    Complex, Intent(InOut) :: localsum(len)
    Complex                :: totalsum(len)
    Integer :: ierr
    Call mpi_allreduce (localsum, totalsum, len,&
         MPI_COMPLEX_TYPE, MPI_SUM, mpi_comm_w, ierr)
    localsum = totalsum
  End Subroutine my_complex_sum_w

  Subroutine my_real_sum_w(localsum,len)
    Integer, Intent(In)    :: len
    Real, Intent(InOut) :: localsum(len)
    Real                :: totalsum(len)
    Integer :: ierr
    Call mpi_allreduce (localsum, totalsum, len,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_w, ierr)
    localsum = totalsum
  End Subroutine my_real_sum_w

  Subroutine my_real_sum_w_0d(localsum)
    Real, Intent(InOut) :: localsum
    Real :: totalsum
    Integer :: ierr
    Call mpi_allreduce (localsum, totalsum, 1,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_w, ierr)
    localsum = totalsum
  End Subroutine my_real_sum_w_0d

  Subroutine my_real_sum_vw(localsum,len)
    Integer, Intent(In):: len
    real, Intent(Inout):: localsum(len)
    real::  totalsum(len)
    integer :: ierr

    Call mpi_allreduce (localsum(1), totalsum(1), len,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_vw, ierr)
    localsum = totalsum
  End Subroutine my_real_sum_vw

  SUBROUTINE my_real_sum_spec(localsum,len)
    Integer, Intent(In):: len
    Real, Intent(Inout):: localsum(len)
    Real               :: totalsum(len)
    Integer:: ierr

    PERFON("sum_sp")
    Call mpi_allreduce(&
         localsum(1), totalsum(1), len,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_spec, ierr)
    localsum = totalsum
    PERFOFF
  END SUBROUTINE my_real_sum_spec



  Real Function sum_vw (rnum)
    Real, Intent(In):: rnum
    Real:: rtmp
    integer :: ierr

    Call mpi_allreduce (rnum, rtmp, 1,&
         MPI_REAL_TYPE, MPI_SUM, mpi_comm_vw, ierr)
    sum_vw = rtmp
  End Function sum_vw

  subroutine my_sum_to_all_real(localsum,len,communicator)
    Integer, Intent(In):: len, communicator
    Real, Intent(Inout):: localsum(len)
    Real               :: totalsum(len)
    integer:: ierr
    Call mpi_allreduce(localsum(1), totalsum(1), len,&
         MPI_REAL_TYPE, MPI_SUM, communicator, ierr)
    localsum=totalsum
  end subroutine my_sum_to_all_real

  subroutine my_sum_to_all_int(localsum,len,communicator)
    Integer, Intent(In):: len, communicator
    integer, Intent(Inout):: localsum(len)
    integer               :: totalsum(len)
    integer:: ierr
    Call mpi_allreduce(localsum(1), totalsum(1), len,&
         MPI_INTEGER, MPI_SUM, communicator, ierr)
    localsum=totalsum
  end subroutine my_sum_to_all_int

  Subroutine my_sum_to_all_real_0d(localsum,communicator)
    Integer, Intent(In):: communicator
    Real, Intent(InOut) :: localsum
    Real :: totalsum
    Integer :: ierr
    Call mpi_allreduce (localsum, totalsum, 1,&
         MPI_REAL_TYPE, MPI_SUM, communicator, ierr)
    localsum = totalsum
  End Subroutine my_sum_to_all_real_0d

  
  subroutine sum_to_all_complex(arr,size_arr,communicator)
    integer:: size_arr, communicator, ierr 
    complex, Intent(Inout) :: arr(size_arr)
    complex:: temparr(size_arr)
    
    Call mpi_allreduce(arr, temparr, size_arr,&
         MPI_COMPLEX_TYPE, MPI_SUM, communicator, ierr)
    arr=temparr
  end subroutine sum_to_all_complex
  
  subroutine my_sum_to_0_complex(arr,size_arr,communicator)
    integer:: size_arr, communicator, ierr, rank 
    complex,dimension(size_arr):: arr, temparr
    
    Call mpi_reduce(arr(1), temparr(1), size_arr,&
         MPI_COMPLEX_TYPE, MPI_SUM, 0, communicator, ierr)

    call mpi_comm_rank(communicator,rank,ierr)

    if (rank.eq.0) then
       arr=temparr
    else
       arr=0.
    endif
  end subroutine my_sum_to_0_complex
  
  subroutine my_sum_to_0_real(arr,size_arr,communicator)
    integer:: size_arr, communicator, ierr, rank
    real,dimension(size_arr):: arr, temparr
    
    Call mpi_reduce(arr(1), temparr(1), size_arr,&
         MPI_REAL_TYPE, MPI_SUM, 0, communicator, ierr)

    call mpi_comm_rank(communicator,rank,ierr)

    if (rank.eq.0) then
       arr=temparr
    else
       arr=0.
    endif
  end subroutine my_sum_to_0_real

  subroutine my_real_gather_to_0(arr,start,loc_len,tot_len,communicator)
    Integer, Intent(In):: start, loc_len, tot_len, communicator
    Real, Intent(Inout):: arr(tot_len)
    Real               :: loc_arr(loc_len)
    integer:: ierr
    loc_arr = arr(start:start+loc_len-1)

    Call mpi_gather(loc_arr(1), loc_len, MPI_REAL_TYPE, &
         &arr(1), loc_len, MPI_REAL_TYPE, 0, communicator, ierr)

  end subroutine my_real_gather_to_0

  subroutine my_complex_gather_to_0(arr,start,loc_len,tot_len,communicator)
    Integer, Intent(In):: start, loc_len, tot_len, communicator
    Complex, Intent(Inout):: arr(tot_len)
    Complex               :: loc_arr(loc_len)
    integer:: ierr
    loc_arr = arr(start:start+loc_len-1)

    Call mpi_gather(loc_arr(1), loc_len, MPI_COMPLEX_TYPE, &
         &arr(1), loc_len, MPI_COMPLEX_TYPE, 0, communicator, ierr)

  end subroutine my_complex_gather_to_0

  SUBROUTINE calculate_test_sum_6D(arr6d,local_sum,global_sum)
    COMPLEX, DIMENSION(:,:,:,:,:,:),INTENT(IN) :: arr6d
    REAL,INTENT(OUT) :: local_sum, global_sum

    !Local variables
    integer :: ierr
    
    local_sum = REAL(SUM(arr6d*CONJG(arr6d)))
    CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,MY_MPI_COMM_WORLD,ierr)
  END SUBROUTINE calculate_test_sum_6D

  SUBROUTINE calculate_test_sum_4D(arr,local_sum,global_sum)
    COMPLEX, DIMENSION(:,:,:,:),INTENT(IN) :: arr
    REAL,INTENT(OUT) :: local_sum, global_sum

    !Local variables
    integer :: ierr

    local_sum = REAL(SUM(arr*CONJG(arr)))
    CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,mpi_comm_xyzw,ierr)
  END SUBROUTINE calculate_test_sum_4D

  SUBROUTINE calculate_test_sum_5D(arr,local_sum,global_sum)
    COMPLEX, DIMENSION(:,:,:,:,:),INTENT(IN) :: arr
    REAL,INTENT(OUT) :: local_sum, global_sum

    !Local variables
    integer :: ierr

    local_sum = REAL(SUM(arr*CONJG(arr)))
    CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,mpi_comm_xyzvw,ierr)
  END SUBROUTINE calculate_test_sum_5D

  SUBROUTINE calculate_test_sum_real_5D(arr,local_sum,global_sum)
    real, DIMENSION(:,:,:,:,:),INTENT(IN) :: arr
    REAL,INTENT(OUT) :: local_sum, global_sum

    !Local variables
    integer :: ierr

    local_sum = SUM(arr*arr)
    CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,mpi_comm_xyzvw,ierr)
  END SUBROUTINE calculate_test_sum_real_5D

  SUBROUTINE calculate_test_sum_2D(arr,local_sum,global_sum)
    COMPLEX, DIMENSION(:,:),INTENT(IN) :: arr
    REAL,INTENT(OUT) :: local_sum, global_sum

    !Local variables
    integer :: ierr
    
    local_sum = REAL(SUM(arr*CONJG(arr)))
    !CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,MY_MPI_COMM_WORLD,ierr)
    CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,communicators(COMM_X),ierr)
  END SUBROUTINE calculate_test_sum_2D

  SUBROUTINE calculate_test_sum_3D(arr,local_sum,global_sum)
    COMPLEX, DIMENSION(:,:,:),INTENT(IN) :: arr
    REAL,INTENT(OUT) :: local_sum, global_sum
    
    !Local variables
    integer :: ierr
    REAL :: global_sumx, global_sumxy
    
    local_sum = REAL(SUM(arr*CONJG(arr)))
    !CALL mpi_reduce(local_sum,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,MY_MPI_COMM_WORLD,ierr)
    CALL mpi_reduce(local_sum,global_sumx,1,MPI_REAL_TYPE,MPI_SUM,0,communicators(COMM_X),ierr)
    IF (my_pex.EQ.0) THEN
       CALL mpi_reduce(global_sumx,global_sumxy,1,MPI_REAL_TYPE,MPI_SUM,0,communicators(COMM_Y),ierr)
       IF (my_pey.EQ.0) THEN
          CALL mpi_reduce(global_sumxy,global_sum,1,MPI_REAL_TYPE,MPI_SUM,0,communicators(COMM_Z),ierr)
       END IF
    END IF
  END SUBROUTINE calculate_test_sum_3D

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Computes process number depending on y,z,v,w,spec process coordinates
  !\param pex process coordinate in x direction
  !\param pey process coordinate in y direction
  !\param pez process coordinate in z direction
  !\param pev process coordinate in v direction
  !\param process coordinate in w direction
  !\param pespec process coordinate in s direction
  !\param pexyzvwspec process number corresponding to specified coordinates
  INTEGER FUNCTION pexyzvwspec(pex, pey, pez, pev, pew, pespec)
    INTEGER::   pex, pey, pez, pev, pew, pespec

    pexyzvwspec = &
         & pex + n_procs_x*(pey + n_procs_y*&
         & (pez + n_procs_z*(pev + n_procs_v*(pew + n_procs_w*pespec))))
  End Function pexyzvwspec


  Subroutine get_systime(time)
    double precision, Intent(Out):: time
    Integer :: ierror

    time = MPI_Wtime()  
    !print*,"time = ",time
!Commented the following line as the intel compiler
!gave the wrong value for wtime_is_global on Hydra
!    IF (.NOT.wtime_is_global) &
         CALL mpi_bcast(time,1,MPI_DOUBLE_PRECISION,&
         &0,MY_MPI_COMM_WORLD,ierror)
    !print*,"time = ",time
  End Subroutine get_systime

  !>Mimics the 'secnds' intrinsic function which is not available on all compilers
  !\param start_time reference start time
  REAL FUNCTION my_secnds(start_time)
    real, intent(in) :: start_time
    double precision :: current_time
    
    call get_systime(current_time)

    my_secnds = current_time - start_time
  END FUNCTION my_secnds

  Subroutine finalize_comm_sim
    Integer :: ierror

    !PRINT*,"Calling free_comm on PE ", mype
    ! free_comm is called in finalize_current_parall
    !CALL free_comm

!In the global code, blacs_exit should be called at this point
!       if (.not.xy_local) call blacs_exit(1)
!However, in this case grid_init crashes with strange error messages
!Let's keep fingers crossed that skipping blacs_exit does not produce
!large memory leaks ...
    CALL mpi_comm_free(MY_MPI_COMM_WORLD,ierror)
    IF (ierror /= 0) STOP 'mpi_comm_free MY_MPI_COMM_WORLD failed!'

  End Subroutine finalize_comm_sim

  !>Determines number of OMP threads
integer function check_omp()
#ifdef WITHOMP
    integer OMP_GET_MAX_THREADS
    external OMP_GET_MAX_THREADS
    
    omp_num_threads=OMP_GET_MAX_THREADS()
    !if ((mype==0).AND.(print_ini_msg)) write(*,"(a,i3,a)") 'running on ',omp_num_threads,' omp threads'
#else
    omp_num_threads=1
#endif
    check_omp=omp_num_threads
  end function check_omp


  subroutine read_parall_nml(par_dir)
    character(Len=*):: par_dir
    integer:: PARFILE=30
    integer:: ioerr
    

    namelist /parallelization/ &
         n_procs_s, n_procs_w, n_procs_v, n_procs_z, n_procs_y, n_procs_x, &
         & n_procs_sim, n_parallel_sims,&
         & max_nps, max_npw, max_npv, max_npz, max_npy, max_npx,&
         & min_nps, min_npw, min_npv, min_npz, min_npy, min_npx

    if (TRIM(par_dir).ne.'skip_parfile') then
       !use autoparallelization as default
       !(important to set it here if gene is called multiple times)
       n_procs_s = -1
       n_procs_w = -1
       n_procs_v = -1
       n_procs_y = -1
       n_procs_z = -1
       n_procs_x = -1
       
       if (par_dir.ne.'') then
          open(PARFILE, file=trim(par_dir)//'/parameters'//&
               &trim(file_extension), form='formatted',&
               &iostat=ioerr, status='old')
       else
          open(PARFILE,file='parameters', &
               &action='read', form='formatted',&
               &iostat=ioerr, status='old')
       end if
       
       If (ioerr /= 0) Then
          Write(*,"(A)") "could neither read parameters from stdin "//&
               &"nor from ./parameters file"
          STOP
       endif
       
       read(PARFILE, nml=parallelization,iostat=ioerr)
       if (ioerr.ne.0) stop &
            'on i/o error: failed to read parallelization namelist'
       close(PARFILE)

       !parallelization parameters are already set by interface
    endif

  end subroutine read_parall_nml

  SUBROUTINE sync_time(time)
    real, intent(inout) :: time
    real :: tmp
    integer :: ierr
    
    ! On input, time is a value which may be slightly different on every CPU.
    ! On output, time is equal on each CPU so that it can be used for
    ! switches etc without the risk of running in different branches
    ! on different CPUs.

    ! Here we use the maximum value for output, one could also use the mean or the minimum

    tmp = time
    CALL MPI_Allreduce(tmp, time, 1, MPI_REAL_TYPE, MPI_MAX, MY_MPI_COMM_WORLD, ierr)
  END SUBROUTINE sync_time

  subroutine reduce_per_thread(local_sum,subcomm,outstring)
    real,intent(IN) :: local_sum
    integer, intent(IN) :: subcomm
    character(len=*) :: outstring

    integer :: tag, my_thread, p_recv, rank,ierr,status(MPI_STATUS_SIZE)
    integer :: n_procs
    real,allocatable,dimension(:,:) :: thread_global_sum
#ifdef WITHOMP
    integer :: omp_get_thread_num

    my_thread = omp_get_thread_num()
#else
    my_thread = 0
#endif
    !write(*,"(2I4,A,ES17.10)") mype,my_thread, ": local_sum = ", local_sum
    call mpi_comm_rank(subcomm,rank,ierr)
    call mpi_comm_size(subcomm,n_procs,ierr)

    allocate(thread_global_sum(0:7,0:n_procs-1))
    tag = 100+my_thread
    if (rank.ne.0) then
       call mpi_send(local_sum,1,MPI_REAL,0,tag,subcomm,ierr)
    else
       thread_global_sum(my_thread,:) = 0.0
       thread_global_sum(my_thread,0) = local_sum
       do p_recv=1,n_procs-1
          call mpi_recv(thread_global_sum(my_thread,p_recv),1,&
               & MPI_REAL,p_recv,tag,subcomm,status,ierr)
       end do
       write(*,"(2I4,A,A15,A,ES17.10)") rank,my_thread," : ",outstring,&
            & " = ",sum(thread_global_sum(my_thread,:))
    end if
    deallocate(thread_global_sum)
  end subroutine reduce_per_thread

End Module communications
