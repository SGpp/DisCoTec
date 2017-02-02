! vtrace_perf.F90
! Uses perfon and perfoff entries to set user defined
! interfaces to the Intel Trace Collector (installed on HLRB)
! Usage: switch on PERFLIB in linux.mk
!        'module load mpi_tracing' before compilation
! Add. info: some global variables are stored in par_mod

SUBROUTINE perfinit

  USE par_mod

!  include 'VT.inc' 
  call vtclassdef('mylib',mylib_handle,ierr)
  handle_index=0

END SUBROUTINE perfinit


SUBROUTINE perfon(str)

  USE par_mod

  CHARACTER(LEN=*):: str
  INTEGER :: name_handle

  handle_index=handle_index+1  
  IF (handle_index.le.100) THEN
    call vtfuncdef(str,mylib_handle,name_handle,ierr)
    call vtbegin(name_handle,ierr)
    name_handle_arr(handle_index)=name_handle  
  END IF

END SUBROUTINE perfon


SUBROUTINE perfoff
  
  USE par_mod

  IF (handle_index.le.100) THEN
    call vtend(name_handle_arr(handle_index), ierr)
  END IF
  handle_index=handle_index-1

END SUBROUTINE perfoff


SUBROUTINE perfout
END SUBROUTINE perfout
