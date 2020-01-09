#include "redef.h"
#include "intrinsic_sizes.h"
!>Contains routines for file i/o
!!
!!In the long term, this module should provide wrappers
!!which either call binary output, MPI-IO, or 
!!FUTILS (HDF5) routines
Module file_io
  Use discretization
  Use communications
  USE,intrinsic :: iso_fortran_env
  Implicit None

  Public :: get_unit_nr, open_file, write3d, stop_file_found, &
       & erase_stop_file, avgflux_stime_file_found, ExB_stime_file_found,&
       & read_in_1d, read_in_2d, write_out_1d, write_out_2d, &
       & get_nr_of_lines, geT_nr_of_substr, create_finished_file
  
  Private

Contains

  !>Returns a file unit number which has not been used before
  !!Warning: subsequent calls will thus return the same number
  !!if no file has been opened in between
  subroutine get_unit_nr(number)
    integer,intent(out):: number

    ! Local variables
    logical :: already_used

    number=31
    DO
       IF ((number.EQ.ERROR_UNIT).OR.(number.EQ.INPUT_UNIT).OR.(number.EQ.OUTPUT_UNIT)) THEN
          number = number + 1
          CYCLE
       END IF
       INQUIRE(unit=number,opened=already_used)
       IF (already_used) THEN
          number = number + 1
       ELSE
          ! found a free file unit
          EXIT
       END IF
    END DO
    
  end subroutine get_unit_nr

  subroutine open_file
    
  end subroutine open_file


  !> Returns the number of lines until the 
  !! first appearance of an empty line or 
  !! eof in a given file
  integer FUNCTION get_nr_of_lines(fileunit) 
    Integer, intent(in) :: fileunit
    Character(Len=128) :: line
    integer :: i, iostat

    i = 0

    read(fileunit,"(A)",iostat=iostat) line
    do while ((iostat.eq.0).and.(TRIM(line).ne.''))
       i = i + 1
       read(fileunit,"(A)",iostat=iostat) line
    enddo
    
    get_nr_of_lines = i

  end FUNCTION get_nr_of_lines


  !>Returns the number of substrings being seperated by
  !!predefined separator in a given string
  !!(used to identify columns in formatted files)
  integer FUNCTION get_nr_of_substr(str) 
    Character(Len=*),intent(in) :: str
    Character :: separator = ' ', tab, ch
    logical :: in_substr
    integer :: i, nr

    tab = ACHAR(9)

    in_substr = .false.
    nr = 0
    do i=1, LEN(str)
       ch = str(i:i)
       if ((str(i:i).eq.separator).or.(str(i:i).eq.tab)) then
          in_substr = .false.
       elseif (.not.in_substr) then 
          nr = nr+1
          in_substr = .true.
       endif
    enddo

    get_nr_of_substr = nr

  end FUNCTION get_nr_of_substr

  !!!******************************************************************!!!
  !> Accepts local 3D space array and either writes with MPI-IO or copies all data to a global array before writing to process 0
  subroutine write3d(filehandle,arr3d,write_pe)
  !****
    Complex, Dimension(li1:li2,lj1:lj2,lk1:lk2):: arr3d
    Complex, Dimension(0:ni0-1,lj1:lj2,lk1:lk2):: arr3dx
    Complex, Dimension(0:ni0-1, 0:nj0-1, lk1:lk2):: arr3df
    Complex, Dimension(0:ni0-1, 0:nj0-1, 0:nz0-1)::  dest3d
    Integer :: j,k,filehandle,write_pe,ierr

    DEBUG(1,"========== BEGIN write3d =========")

    ! Gather arr3d over x-Distribution on pex==0.
    Do k = lk1, lk2
       Do j = lj1, lj2
          Call mpi_gather(arr3d(li1,j,k), li0, MPI_COMPLEX_TYPE,&
               arr3dx(0,j,k), li0, MPI_COMPLEX_TYPE,&
               0, mpi_comm_x, ierr)
       Enddo
    Enddo
    
    ! Gather arr3d over y-Distribution on pey==0.
    Do k = lk1, lk2
       Call mpi_gather(arr3dx(0,lj1,k), Size(arr3dx(:,:,k)), MPI_COMPLEX_TYPE,&
            arr3df(0,0,k), Size(arr3dx(:,:,k)), MPI_COMPLEX_TYPE,&
            0, mpi_comm_y, ierr)
    Enddo
    
    ! Gather arr3df over z-Distribution on pez=0
    CALL mpi_gather(arr3df(0,0,lk1),SIZE(arr3df(:,:,:)), MPI_COMPLEX_TYPE,&
         dest3d(0,0,0),SIZE(arr3df),MPI_COMPLEX_TYPE,&
         0, mpi_comm_z, ierr)

    IF (mype.eq.write_pe) WRITE(filehandle) dest3d

    DEBUG(1,"========== END write3d =========")
  End Subroutine write3d

  
  !>Check for STOP files
  FUNCTION stop_file_found(file_extension,simtimelim) result(found)
    character(len=*), intent(in) :: file_extension
    Integer :: STOPFILE, ierr
    Logical :: found
    real :: simtimelim, new_simtimelim
    found = .false.
    new_simtimelim = -1

    if (mype==0) then
       Inquire(file='GENE.stop',exist=found)

       If ((.not.found).and.(trim(file_extension).ne.'.dat')) Then
          Inquire(file='GENE.stop'//trim(file_extension),exist=found)
          !Those files can be erased immediately while GENE.stop
          !has to be kept until all instances (scangene or GENE library
          !mode) are returning
          If ((found).and.(mype==0)) then
             CALL get_unit_nr(STOPFILE)
             OPEN(STOPFILE,FILE='GENE.stop'//trim(file_extension))
             READ(STOPFILE, "(F13.6)",iostat=ierr) new_simtimelim
             IF (ierr.ne.0) new_simtimelim = -1
             CLOSE(STOPFILE,STATUS='delete')   
          Endif
       Endif
    endif
    
    !broadcast found value
    call mpi_bcast(found,1,MPI_LOGICAL,0,MY_MPI_COMM_WORLD,ierr)

    if (found) then
       call mpi_bcast(new_simtimelim,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)
       
       if (new_simtimelim.ge.0.0) then
          found=.false.
          simtimelim=new_simtimelim
       endif
    endif
    
  End Function stop_file_found
  
  !>Create GENE.finished file, if not found
  Subroutine create_finished_file
    integer :: FINFILE
    logical :: found

    if (mype==0) then
       Inquire(file='GENE.finished',exist=found)
       IF (.not.found) THEN
          CALL get_unit_nr(FINFILE)
          OPEN(FINFILE,FILE='GENE.finished')
          CLOSE(FINFILE)
       ENDIF
    endif
  End Subroutine create_finished_file


  !>Erase GENE.stop file if found
  Subroutine erase_stop_file
    integer :: STOPFILE
    logical :: found

    if (mype_gl==0) then
       Inquire(file='GENE.stop',exist=found)
       IF (found) THEN
          CALL get_unit_nr(STOPFILE)
          OPEN(STOPFILE,FILE='GENE.stop')
          CLOSE(STOPFILE,STATUS='delete')
          !    OPEN(STOPFILE,FILE='GENE.stop'//trim(file_extension))
          !    CLOSE(STOPFILE,STATUS='delete')    
       ENDIF
    endif
  End Subroutine erase_stop_file

  !>Check for time averages start files
  FUNCTION avgflux_stime_file_found(file_extension,time,avgflux_stime) result(found)
    character(len=*), intent(in) :: file_extension
    real, intent(in) :: time
    real, intent(inout) :: avgflux_stime
    real :: new_avgflux_stime
    Integer :: STOPFILE, ierr
    Logical :: found
    found = .false.
    
    if (mype==0) then
       Inquire(file='avgflux_stime'//trim(file_extension),exist=found)
       IF (found) THEN
          CALL get_unit_nr(STOPFILE)
          OPEN(STOPFILE,FILE='avgflux_stime'//trim(file_extension))

          READ(STOPFILE, "(F13.6)",iostat=ierr) new_avgflux_stime
          IF (ierr.eq.0) THEN
             avgflux_stime=new_avgflux_stime
          ELSE !if file empty use current time
             avgflux_stime=time
          ENDIF
          CLOSE(STOPFILE,STATUS='delete')

       ENDIF
    endif
    
    !broadcast found value
    call mpi_bcast(found,1,MPI_LOGICAL,0,MY_MPI_COMM_WORLD,ierr)
    !broadcast avgflux_stime value
    if (found) then
       call mpi_bcast(avgflux_stime,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)    
    endif

  End Function avgflux_stime_file_found

  
  !>Check for time averages start files
  FUNCTION ExB_stime_file_found(file_extension,time,ExB_stime) result(found)
    character(len=*), intent(in) :: file_extension
    real, intent(in) :: time
    real, intent(inout) :: ExB_stime
    real :: new_ExB_stime
    Integer :: STOPFILE, ierr
    Logical :: found
    found = .false.
    
    if (mype==0) then
       Inquire(file='ExB_stime'//trim(file_extension),exist=found)
       IF (found) THEN
          CALL get_unit_nr(STOPFILE)
          OPEN(STOPFILE,FILE='ExB_stime'//trim(file_extension))          
          READ(STOPFILE, "(F13.6)",iostat=ierr) new_ExB_stime
          IF (ierr.eq.0) THEN
             ExB_stime=new_ExB_stime
          ELSE !if file empty use current time
             ExB_stime=time
          ENDIF
          CLOSE(STOPFILE,STATUS='delete')
       ENDIF
    endif
    
    !broadcast found value
    call mpi_bcast(found,1,MPI_LOGICAL,0,MY_MPI_COMM_WORLD,ierr)
    !broadcast ExB_stime value
    if (found) then
       call mpi_bcast(ExB_stime,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr) 
    endif
    
  End Function ExB_stime_file_found


  !>Reads a 1D field from a given file
  SUBROUTINE READ_IN_1D(data_name,data,N,fileunit)
    IMPLICIT NONE
    INTEGER :: N, i, fileunit, iostat
    REAL, DIMENSION(1:N) :: data
    CHARACTER(*) :: data_name
    CHARACTER(len=128) :: line
    
    line = ''
    iostat = 0
    rewind(fileunit)
    do while ((iostat.eq.0).and.(trim(line).ne.trim(data_name)))
       read(fileunit,"(A)",iostat=iostat) line
    end do

    IF (iostat.eq.0) THEN
       READ(fileunit,*) (data(i),i=1,N)
    ELSE
       WRITE(*,"((3A))") "No data found for ",data_name," in read_in_1d"
       STOP
    ENDIF
       
  END SUBROUTINE READ_IN_1D


  !>Reads a 2D field from a given file.
  SUBROUTINE READ_IN_2D(data_name,data,N1,N2,fileunit)
    IMPLICIT NONE
    INTEGER :: N1, N2,i ,j, fileunit, iostat
    REAL, DIMENSION(1:N1,1:N2) :: data
    CHARACTER(*) :: data_name
    CHARACTER(len=128) :: line

    line = ''
    iostat = 0
    rewind(fileunit)
    do while ((iostat.eq.0).and.(trim(line).ne.trim(data_name)))
       read(fileunit,"(A)",iostat=iostat) line
    end do

    IF (iostat.eq.0) THEN
       READ(fileunit,"(16ES20.10)") ((data(i,j),i=1,N1),j=1,N2)
    ELSE
       WRITE(*,"((3A))") "No data found for ",data_name," in read_in_2d"
       STOP
    ENDIF
  END SUBROUTINE READ_IN_2D


  !>Writes a 1D field to the given file.  
  SUBROUTINE WRITE_OUT_1D(data_name,data,N,fileunit)
    IMPLICIT NONE
    INTEGER :: N, i, fileunit
    REAL, DIMENSION(1:N) :: data
    CHARACTER(*) :: data_name
    ! WRITE ASCII
    WRITE(fileunit,"(A)") data_name
    WRITE(fileunit,*) (data(i),i=1,N)
  END SUBROUTINE WRITE_OUT_1D


  !>Writes a 2D field to the given file.
  SUBROUTINE WRITE_OUT_2D(data_name,data,N1,N2,fileunit)
    IMPLICIT NONE
    INTEGER :: N1, N2,i ,j, fileunit
    REAL, DIMENSION(1:N1,1:N2) :: data
    CHARACTER(*) :: data_name
    ! WRITE ASCII
    WRITE(fileunit,"(A)") data_name
    WRITE(fileunit,"(16ES20.10)") ((data(i,j),i=1,N1),j=1,N2)
  END SUBROUTINE WRITE_OUT_2D
  
End Module file_io
