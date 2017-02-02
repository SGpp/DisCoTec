program get_mpi_info
  use mpi
  implicit none

  integer :: ierr,provided

  call mpi_init_thread(MPI_THREAD_MULTIPLE,provided,ierr)

  select case (provided)
  case(MPI_THREAD_MULTIPLE)
     print*,"MPI supports MPI_THREAD_MULTIPLE."
  case(MPI_THREAD_SERIALIZED)
     print*,"MPI supports MPI_THREAD_SERIALIZED."
  case(MPI_THREAD_FUNNELED)
     print*,"MPI supports MPI_THREAD_FUNNELED."
  case(MPI_THREAD_SINGLE)
     print*,"MPI supports MPI_THREAD_SINGLE."
  end select
  call mpi_finalize(ierr)
end program get_mpi_info
