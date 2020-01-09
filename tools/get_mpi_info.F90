program get_mpi_info
  implicit none

  integer :: info,ierr,nkeys,ikey
  character(len=255) :: key

  call mpi_init(ierr)

  call mpi_info_create(info,ierr)
  call mpi_info_set(info,"romio_ds_write", "disable",ierr)
  call mpi_info_set(info,"romio_ds_read", "disable",ierr)
  call mpi_info_get_nkeys(info,nkeys,ierr)
  print*,"We have ",nkeys," keys in info."
  do ikey=0,nkeys-1
     call mpi_info_get_nthkey(info,ikey,key,ierr)
     print*,"Key no. ",ikey," is ",key
  end do
  call mpi_info_free(info,ierr)

  call mpi_finalize(ierr)
end program get_mpi_info
