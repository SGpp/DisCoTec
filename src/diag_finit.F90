#include "redef.h"
#include "intrinsic_sizes.h"
#undef FLR_NUMERIC
!> This module is responsible for the fielddiff.dat output
!! which computes the fields from g1(t) - g1(0)
MODULE diag_finit
  USE par_mod
  USE file_io
  USE vel_space
  USE communications
  USE aux_fields

  IMPLICIT NONE

  PUBLIC :: initialize_all_diags_finit, exec_all_diags_finit, &
    finalize_all_diags_finit, mem_est_diag_finit

  PRIVATE 

  INTEGER :: istep_field = -1
  CHARACTER(LEN=8) :: filestat='replace', filepos='rewind'
  INTEGER :: FIELDDIFFFILE
  COMPLEX, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: g_1_init, g_1_mod
  LOGICAL :: first_fielddiff = .true.

  !administration
  REAL :: last_exec_diag_time = -3.14159265

CONTAINS

!!!******************************************************************!!!

  !>Give an estimate of the memory requirements of this module
  REAL FUNCTION mem_est_diag_finit(mem_req_in)
    REAL :: mem_req_in
    REAL :: mem_loc = 0.0

    IF (istep_field .GT. 0) mem_loc = mem_est_diag_fielddiff(mem_loc)

    mem_est_diag_finit = mem_req_in + mem_loc
  END FUNCTION mem_est_diag_finit

!!!******************************************************************!!!

  !!Each individual GENE diagnostic initialization should be called 
  !!by this routine
  SUBROUTINE initialize_all_diags_finit(p_istep_field)
    INTEGER, INTENT(IN) :: p_istep_field

    istep_field = p_istep_field
    IF (istep_field .GT. 0) CALL initialize_diag_fielddiff

  END SUBROUTINE initialize_all_diags_finit

!!!******************************************************************!!!
  SUBROUTINE exec_all_diags_finit(itime,time)
    INTEGER, INTENT(IN) :: itime
    REAL, INTENT(IN) :: time

    !avoid double entries if exec_all_diags is, e.g., called during
    !convergence exit
    IF (time .EQ. last_exec_diag_time) THEN
      RETURN
    ELSE
      last_exec_diag_time = time
    ENDIF

    IF (istep_field .GT. 0) THEN
      IF (MODULO(itime,istep_field) .EQ. 0) THEN
        g_1_mod = g_1 - g_1_init
        CALL calc_aux_fields(g_1_mod,emfields,f_,.false.)
        CALL diag_fielddiff
        CALL calc_aux_fields(g_1,emfields,f_,.false.)
      END IF
    END IF

  END SUBROUTINE exec_all_diags_finit

!!!******************************************************************!!!

  !>Finalizes all GENE internal diagnostics
  SUBROUTINE finalize_all_diags_finit
    IF (istep_field .GT. 0) CALL finalize_diag_fielddiff
  END SUBROUTINE finalize_all_diags_finit

!!!******************************************************************!!!
  
  !>Give an estimate of the memory requirements of diag_field
  REAL FUNCTION mem_est_diag_fielddiff(mem_req_in)
    REAL :: mem_req_in
    REAL :: mem_loc = 0.0

    !g_1_init,_mod
    mem_loc = 2 * SIZE_OF_COMPLEX_MB * lijklmn0

    mem_est_diag_fielddiff = mem_loc + mem_req_in

  End Function mem_est_diag_fielddiff

!!!******************************************************************!!!

  SUBROUTINE initialize_diag_fielddiff
    IMPLICIT NONE

    call get_unit_nr(FIELDDIFFFILE)

    IF (mype .EQ. 0) OPEN(FIELDDIFFFILE,FILE=TRIM(diagdir)//'/fielddiff'//&
      TRIM(file_extension),FORM='unformatted',STATUS=filestat,POSITION=filepos)

    ALLOCATE(g_1_init(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    ALLOCATE(g_1_mod(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))

  END SUBROUTINE initialize_diag_fielddiff

  !*****************************************************************
  !>FIELDDIFF diagnostic
  !!
  !!Writes 3D data (kx,ky,z) of recalibrated electrostatic field
  !!phi and electromagnetic field Apar (if beta > 0)
  !!and Bpar (if beta > 0 and bpar = .t.)
  !!to fielddiff.dat file
  SUBROUTINE diag_fielddiff
    IMPLICIT NONE
    INTEGER :: o

    PERFON('d_fielddiff')

    IF ((my_pev+my_pew+my_pespec).eq.0) THEN
      IF (first_fielddiff) THEN
        first_fielddiff = .false.
        g_1_init = g_1
      END IF
      IF (mype.eq.0) WRITE(FIELDDIFFFILE) time

      DO o = 1, n_fields
        CALL write3d(FIELDDIFFFILE,emfields(li1:li2,lj1:lj2,lk1:lk2,o),0)
      END DO
      IF (mype.eq.0) call flush(FIELDDIFFFILE)
    END IF

    PERFOFF

  END SUBROUTINE diag_fielddiff

  SUBROUTINE finalize_diag_fielddiff
    IF (mype .EQ. 0) CLOSE(FIELDDIFFFILE)
    DEALLOCATE(g_1_init,g_1_mod)
  END SUBROUTINE finalize_diag_fielddiff

END MODULE diag_finit
