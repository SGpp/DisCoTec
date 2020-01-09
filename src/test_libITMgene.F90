#include "redef.h"
!>Test program for the ITM interface library libITMgene
!!Compile with: gmake -f ../makefile lib LIBN=libITMgene
PROGRAM test_libITMgene

  USE ITM_Types
  USE ITM_Constants
  USE Euitm_schemas
  use xml_file_reader
  use read_structures
  use write_structures

  IMPLICIT NONE

  TYPE (type_equilibrium), pointer ::  eq(:)
  TYPE (type_coreprof), pointer :: coreprof(:)
  TYPE (type_coretransp), pointer :: coretransp(:)
  type (type_param) :: code_parameters

  interface
    subroutine ITMgene(eq, coreprof, coretransp, code_parameters)
      use euitm_schemas
      type (type_equilibrium), pointer :: eq(:)
      type (type_coreprof), pointer :: coreprof(:)
      type (type_coretransp), pointer :: coretransp(:)
      type (type_param) :: code_parameters
    end subroutine ITMgene
  end interface

  INTEGER(ITM_I4) :: i
  INTEGER(ITM_I4) :: nrho_prof = 20, nion = 1, npsi_eq = 32, &
       &nangle_eq = 24, neta_eq = 64, nrho_transp

  REAL(r8) :: pi=itm_pi
  REAL(r8) :: tpi=2.0_r8*itm_pi

  REAL(r8) :: mu_0 = itm_mu0
  REAL(r8) :: kb=itm_ev
  REAL(r8) :: ee=itm_qe
  REAL(r8) :: cc=1.0_r8, lcoul=14.0_r8

  REAL(r8) :: me=itm_me
  REAL(R8), DIMENSION(:), POINTER :: eta,ra,psi,qq,jphi

  INTEGER :: ierr

  call mpi_init(ierr)
  
!!!4 lines by anders nilsson
  
  allocate (coreprof(1))
  
  allocate (eq(1))
  
  call open_read_file         (10,'IMP4Init')
  
  call read_cpo               (coreprof(1),'IMP4InitCoreprof')
  
  call read_cpo               (eq(1),'IMP4InitEquil')
  
  call close_read_file
  
  
!...  get code specific parameters into xml files
!...  the argument left blank is the defaults here set by the compiler
!...  for Kepler they want this put in manually as well for reproducibility
!...  the last argument is the XML schema file

!!  CALL Get_Code_Parms(code_parameters, 'etaigb.xml', '', 'etaigb.xsd')
  CALL Fill_Param(code_parameters, 'gene_in.xml', '', 'gene_in.xsd')

!...  call transport routine


  CALL ITMgene(eq, coreprof, coretransp, code_parameters)


  nrho_transp = SIZE(coretransp(1)%rho_tor)



  !...  write profiles
!IF (mype == 0) THEN
#if 0
  call open_write_file(10,'IMP4geneResult')
  call write_cpo(coretransp(1),'IMP4geneCoretransp')
  call close_write_file
#endif
!END IF

#if 0
  OPEN (10, file = 'fluxes.dat', form = 'formatted')
  WRITE(10,*) ' i   rho_tor       Gamma        Q_e      Q_i'
  DO i=1,nrho_transp
     WRITE (10,"(I3,4G11.3)") i, &
          coretransp(1)%rho_tor(i), &
          coretransp(1)%ne_transp%flux(i), &
          coretransp(1)%te_transp%flux(i), &
          coretransp(1)%ti_transp%flux(i,1)
  END DO
  CLOSE (10)
#endif

#if 0
  OPEN (10, file = 'prof.dat', form = 'formatted')
  WRITE(10,*) ' i   rho_tor       ne         Te         Ti       D        Q_e      Q_i'
  DO i=1,nrho_prof
!print*, i, (i-1)*nrho_transp/nrho_prof+1
     WRITE (10,"(I3,7G11.3)") i, &
!.     WRITE (10,*) &
          coreprof(1)%rho_tor(i), &
          coreprof(1)%ne%value(i), &
          coreprof(1)%te%value(i), &
          corepro(1)%ti%value(i,1), &
!For testing, coretransp is evaluated on a different grid compared to coreprof
          coretransp(1)%ne_transp%flux((i-1)*nrho_transp/nrho_prof+1), &
          coretransp(1)%te_transp%flux((i-1)*nrho_transp/nrho_prof+1), &
          coretransp(1)%ti_transp%flux((i-1)*nrho_transp/nrho_prof+1,1)
  END DO
  CLOSE (10)

#endif

!print*, '... done'
  call mpi_finalize(ierr)

END PROGRAM test_libITMgene
