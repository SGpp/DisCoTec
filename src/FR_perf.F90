!!>Performance measurement routines based on Fortran MPI_WTIME
!! Code added or modified under the EUFORIA project, 
!! FP7-INFRASTRUCTURES-2007-1, Grant 211804
!!
!! Code added or modified by Fiona Reid from EPCC, The University of Edinburgh,
!! on 01/06/2010
!
!! Code modified (introduced new Module FR_perf_mod) by Tobias Goerler
!! on 01/07/2010, 05/07/2010
!!
Module FR_perf_mod
  Implicit None
  
  
  ! kind specific variables added by Fiona Reid to avoid using
  ! non standard kind = 8, real(4) etc. 
  integer, parameter :: my_kind_sp = kind(1.0)
  integer, parameter :: my_kind_dp = kind(1.0d0)
  integer, parameter :: my_mpi_offset_kind = kind(1.0d0)
  integer, parameter :: my_mpi_address_kind = kind(1.0d0)
  
  ! FR: Variables for use with FR's timing calls April 2010
  !   
  integer, parameter :: MAX_TIMERS=200                ! Max no. of timers
  integer :: timer_id(1:MAX_TIMERS)                   ! Array to store timer handles
  character (LEN=50) :: timer_name(1:MAX_TIMERS)      ! Array to store timer names
  real(kind=my_kind_dp) :: timer_value(1:MAX_TIMERS)  ! Arrays to store the times
  real(kind=my_kind_dp) :: start_time(1:MAX_TIMERS)
  real(kind=my_kind_dp) :: end_time(1:MAX_TIMERS)
  integer :: call_count(1:MAX_TIMERS)
  integer :: current_timer_id = 0
  integer :: max_timer_id = 0
  
  !
  ! Variables inserted so that manual timers can be inserted
  !
  real(kind=my_kind_dp) :: fr_start_time(1:10)
  real(kind=my_kind_dp) :: fr_end_time(1:10)
  real(kind=my_kind_dp) :: fr_time_taken(1:10)
  
End Module FR_perf_mod


SUBROUTINE perfinit
  use FR_perf_mod
  implicit none
 
  ! Zero all index array, timer arrays and call count array before use
  timer_id(:) = 0 
  timer_value(:) = 0.0d0
  start_time(:) = 0.0d0
  end_time(:) = 0.0d0
  call_count(:) = 0

  return

END SUBROUTINE perfinit

SUBROUTINE perfinit_manual
  use FR_perf_mod
  implicit none
 
  ! Zero the timer arrays
  fr_start_time(:) = 0.0d0
  fr_end_time(:) = 0.0d0
  fr_time_taken(:) = 0.0d0

  return

END SUBROUTINE perfinit_manual

SUBROUTINE perfon(str)
  use FR_perf_mod
  use mpi
  implicit none 
  integer :: i
  logical :: new_timer
  CHARACTER(LEN=*):: str

  ! Current_timer_id tells us how many timers are currently active
  ! or in terms of stacks, its the stack pointer  
  current_timer_id = current_timer_id + 1

  ! Test whether this timer has been used before
  new_timer = .true.
  do i = 1, MAX_TIMERS
    if (str==timer_name(i)) then 
      timer_id(current_timer_id) = i
      new_timer = .false.
    end if
  end do
  ! For new timer need to increment max_timer_id and update timer_id() array
  if (new_timer) then ! First time timer has been used
      max_timer_id = max_timer_id + 1
      timer_id(current_timer_id) = max_timer_id
  end if 
  
  ! Set name of timer, actually only needs to be done for new timers, thus can 
  ! go in if statment above
  timer_name(timer_id(current_timer_id)) = str
  ! Increment the call_count for the current timer
  call_count(timer_id(current_timer_id)) =  call_count(timer_id(current_timer_id)) + 1
 
#ifdef DEBUG_FR
  write(*,*)"Current timer has name= ",str," and id = ",&
    timer_id(current_timer_id)
#endif

  ! Start timer
  start_time(timer_id(current_timer_id)) = MPI_Wtime()

  return

END SUBROUTINE perfon

SUBROUTINE perfoff
  use FR_perf_mod
  use mpi
  implicit none
  real(kind=my_kind_dp) :: timediff
  integer :: running_timer
  timediff = 0.0d0 ! initialise to zero

  ! Running_timer = id of the currently running timer, i.e. the one of the
  ! top of the stack
  running_timer = timer_id(current_timer_id)
#ifdef DEBUG_FR
  write(*,*)"Stopping the current timer with id = ",running_timer !current_timer_id
#endif

  ! Stop the timer
  end_time(running_timer) = MPI_Wtime()
  ! Compute time difference for this instance
  timediff = end_time(running_timer) - start_time(running_timer)
  ! Compute the running total 
  timer_value(running_timer) = timer_value(running_timer) + timediff

  ! Decrement current_timer_id now that we're done with this timer
  current_timer_id = current_timer_id - 1 

  return

END SUBROUTINE perfoff

SUBROUTINE perfout(str)
  use FR_perf_mod
  implicit none
  integer :: i
  CHARACTER(LEN=*):: str

  if (trim(str).eq.'GENE') then
     call perfout_all
  else
     ! Print out the timing results to file or screen for a single timer
     ! Find correct data by searching timer_name() array
     do i = 1, MAX_TIMERS
        if(str==trim(timer_name(i))) then
           write(*,99)i,timer_name(i),timer_value(i),call_count(i)
        end if
     end do
  endif

  99 format(i4,1x,a20,1x,f16.5,1x,i10)

  return

END SUBROUTINE perfout

SUBROUTINE perfout_all
  use FR_perf_mod
  implicit none
  integer :: i

  ! Print out all the timing results to file or screen
  write (*,'(A)') ''
  write (*,'(A)') '*********** performance results ***********'
  write (*,'(A)') '                                   Incl.   '
  write (*,'(A)') 'Timer name        #calls  Time(s)    %     '
  write (*,'(A)') '-------------------------------------------'
  do i = 1, max_timer_id
      write(*,99) timer_name(i),call_count(i),timer_value(i),100.*timer_value(i)/timer_value(1)
  end do 
  write (*,'(A)') ''

  99 format(a16,1x,i7,1x,f8.3,1x,f5.1)
  
  return

END SUBROUTINE perfout_all


SUBROUTINE perfout_all_manual
  use FR_perf_mod
  implicit none
  integer :: i
   
  ! Print out all the timing results to file or screen
  do i = 1, 10
      write(*,*)i," ",fr_time_taken(i)
  end do 
  
  return

END SUBROUTINE perfout_all_manual



