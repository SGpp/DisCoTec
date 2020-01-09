!
#ifndef HLST_HAC_USE_REAL
!#define HLST_HAC_USE_REAL 0
#error
#endif
#ifndef HLST_HAC_USE_CMSK
!#define HLST_HAC_USE_CMSK 1
#error
#endif
#ifndef HLST_HAC_USE_ADIM
!#define HLST_HAC_USE_ADIM 6
#error
#endif
!
! ##################### BEGIN USER DOCUMENTATION #######################
!>  Fortran module for ADIOS based checkpoint of MPI applications; 
!! developed by the EFDA High Level Support Team within the HLST-ITM-ADIOS project.
!!
!! Checkpoint / restart functionality is provided via the following routines, each documented individually.
!! \li #hac_init
!! \li #hac_write
!! \li #hac_info
!! \li #hac_read 
!! \li #hac_exit
!!
!! This module is specialized in high throughput checkpointing of a single \e global array on a parallel filesystem.
!!
!!
!! The user is expected to specify the supported numerical type and number of 
!! dimensions at module compile time via three preprocessor symbols, here described:
!! \li If \c HLST_HAC_USE_REAL is 1, will use REAL (otherwise COMPLEX) as type.
!! \li If \c HLST_HAC_USE_CMSK is 1, will not specify the type's KIND, leaving it compiler set; otherwise a value of 8 will be used.
!! \li Value of \c HLST_HAC_USE_ADIM is the dimension of the supported array; either of 2,4,6.
!!
!! Before usage, the #hac_init initialization subroutine shall be called.
!! 
!! The module supports a single global array per checkpoint file.
!! So to checkpoint two global arrays, two separate checkpoints are necessary.
!! To checkpoint a local array, this module is not appropriate.
!! 
!! The checkpoint file format is specific to the ADIOS library.
!! It is possible to use tools distributed with ADIOS to manipulate that file; 
!! however this should not be necessary, because this module contains all that is needed to read the data back, e.g. in a converter.
!! 
!! Written metadata are kept at a minimum: program's metadata are expected to be saved by the user separately.
!! Only information about the global array is written to the checkpoint file: nothing about the originating application instance.
!! A given checkpoint file can be read back from a number of tasks independent from the number at write time.
!! 
!! A multidimensional global array shape is defined by its lower and upper bounds indices.
!! Size on each dimension is equal to the difference of upper and lower bounds plus one.
!! No dimension can be zero; there is no hardcoded upper size limit.
!! A global array is represented by \e slice arrays allocated on each participating MPI task.
!! Each slice array can have arbitrary lower and upper bounds and size (accessible via the LBOUND, UBOUND, SHAPE intrinsics).
!! 
!! When writing (with #hac_write), lower bounds must specified explicitly by the user; 
!! the local slice size will be assumed to be the array size unless explicitly provided.
!! It is responsibility of the user to provide non-overlapping array slices.
!!
!! The global array's upper and lower bounds will be computed from the array slices bounds.
!! The upper bounds will be assumed to be the maximal upper bounds across the array slices; 
!! the lower bounds will be assumed to be the minimal lower bounds across the array slices.
!!
!! Two conventions are supported for the global array's lower bounds: 1-based or 0-based.
!! Either of the two conventions can be used, but shall be consistent on all of the dimensions.
!! When writing, the convention will be determined according to the minimal lower bound encountered.
!! When reading, 0 will be assumed to be the base unless otherwise specified.
!! 
!! When writing or reading, it is not necessary for all of the array slices to be of the same size, 
!! but in no case can a dimension be zero.
!! So, reads (with #hac_read) can be performed on any portion of the global array, 
!! as long as indices within its lower and upper bounds are specified.
!! Before reading a checkpoint file, a user may inquire (with #hac_info) into the dimensions and type of the global array.
!!
!! After usage, the #hac_exit initialization subroutine shall be called.
!!
!! No environment variable is accessed directly by this module, exception made for the subroutine #hac_test.
!!
!! Official documentation consists of a man page (this document) and the following example files:
!!
!! \li \e hac_gene.F90
!! \li \e hac_gemza.F90
!! \li \e hac_nemorb.F90
!! \li \e hac_dump.F90
!! \li \e Makefile
!! 
!!
!! \warning ADIOS-1.5 needs to allocate a buffer at initialization time or
!!     before the first I/O; as a consequence, one cannot allocate a
!!     larger buffer at the time of e.g. the second I/O.
!!     <em> So in case the second I/O is expected to be larger than the first, 
!!     one should allocate enough at the beginning explicitly</em>.
!!     Otherwise, expect failures.
!! 
!! \remark Only \c hac_* prefixed Fortran identifiers are made \c PUBLIC by this module.
!! \remark This module uses ADIOS, so many \c adios_ prefixed C symbols may appear in your linked application.
!! \remark As a clean programming practice (to avoid any potential name clash), 
!! \remark please don't name any of your identifiers as prefixed with \c 'hac_' or \c 'adios_'.
!!
!! \note
!! \li This module requires ADIOS-1.5.
!! \li In the current setup, a file (with name specified by the user) and a 
!! \li likely named (with ".dir" suffixed) directory will be created. 
!! \li They shall always be used together, and addressed by the file path.
!! \li Currently using hardcoded Lustre filesystem related values in A_NO, A_SC, A_NA.
!! \li Local I/O times are measured without \c FSYNC() semantics.
!! \li Handling only of certain simple errors like e.g. 'no file found' is supported; 
!! \li failure to comply with the specifications may go undetected and behaviour undefined.
!! \li One precision \c KIND and either \c REAL or \c COMPLEX is supported.
!! \li Assuming \c INTEGER is 4 bytes. Will break otherwise.
!! \li If the \c OPTIONAL error status variables are omitted, errors are fatal and handled via the Fortran \c STOP statement.
!! \li Reading a checkpoint file written with a different type or KIND may go undetected; you may end up reading garbage.
!! 
!! \version 20140310
!! 
!! \author Michele Martone, EFDA High Level Support Team
!! 
! ###################### END USER DOCUMENTATION ########################
! ############## WHAT FOLLOWS IS INTERNAL DOCUMENTATION ################
!
! Throughout this module, UPPERCASE text has been used for: PARAMETER 
! identifiers, language tokens, intrinsic functions; an underscore is 
! used trailing optional arguments identifiers.
! When USE'ing modules, relevant module identifiers follow in an ONLY
! list.
!
MODULE hlst_adios_checkpoint
!
! What follows in this module is all internals, except for the PUBLIC
! interface mentioned before.
!
#if (HLST_HAC_USE_REAL==1)
#define HAC_FIELD REAL
#else
#define HAC_FIELD COMPLEX
#endif
#if (HLST_HAC_USE_CMSK==1)
#define HAC_KSPEC
#else
#define HAC_KSPEC (8)
#endif
#define HAC_NTYPE HAC_FIELD HAC_KSPEC
#define HLST_HAC_COMPILING_WITH_GNU_STANDARD 1
#define HLST_HAC_GOTO_END_ON_ERROR(LABEL) CALL hac_herr; IF (ierr.NE.OKV_.OR.aerr.NE.OKV_) GOTO LABEL
#define HLST_HAC_HANDLE_IERR(IERR_,IERR) IF (PRESENT(IERR_)) IERR_ = IERR
#define HAC_SET_IF_PRESENT(VAR_,VAR) IF(PRESENT(VAR_))VAR=VAR_
!
#if   (HLST_HAC_USE_ADIM==6)
#define HAC_DSPEC DIMENSION(:,:,:,:,:,:)
#elif (HLST_HAC_USE_ADIM==4)
#define HAC_DSPEC DIMENSION(:,:,:,:)
#elif (HLST_HAC_USE_ADIM==2)
#define HAC_DSPEC DIMENSION(:,:)
#else
#error "Unsupported HLST_HAC_USE_ADIM case."
#endif
#define HLST_HAC_EXTRA_METADATA 1
!
  USE adios_write_mod, ONLY: adios_write, adios_close, adios_open,&
          & adios_write, adios_group_size, adios_define_var, &
          & adios_finalize , adios_select_method, adios_declare_group,&
          & adios_allocate_buffer
  ! Note: in ADIOS-1.5 adios_init_noxml is defined only in C sources.
  USE adios_defs_mod, ONLY: adios_double, adios_long, adios_complex, &
        & adios_unknown, adios_real 
  USE adios_read_mod, ONLY: adios_read_init_method,&
        & adios_read_open_file, adios_get_scalar, ADIOS_READ_METHOD_BP,&
        & adios_get_scalar, adios_read_close, adios_read_open_file,&
        & adios_selection_boundingbox, adios_schedule_read, &
        & adios_perform_reads, adios_read_close, &
        & adios_selection_delete, adios_selection_delete,&
        & adios_read_finalize_method
  USE mpi, ONLY: MPI_INTEGER8, MPI_INTEGER, MPI_COMM_NULL, MPI_SUM, &
          &MPI_WTIME, MPI_COMM_WORLD, MPI_MAX ! MPI_Allreduce is not here !?
  USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT, OUTPUT_UNIT
  IMPLICIT NONE
  PRIVATE
  SAVE
  REAL HAC_KSPEC,PARAMETER,PRIVATE :: SAMPLEVAR = 0 
! ############## BEGIN OF PUBLIC DECLARATIONS ##########################
  PUBLIC :: hac_init,hac_exit,hac_read,hac_write !< Init/exit/read/write
  PUBLIC :: hac_info, hac_test !
  INTEGER,PARAMETER,PUBLIC :: HAC_KIND = KIND(SAMPLEVAR) !< Supported KIND of the module's numerical type.
  INTEGER,PARAMETER,PUBLIC :: HAC_QUIET = 0   !< Verbosity level for hac_init: quiet value.
  INTEGER,PARAMETER,PUBLIC :: HAC_VERBOSE = 1 !< Verbosity level for hac_init: verbose value.
  INTEGER,PARAMETER,PUBLIC :: HAC_VERBOSE_ALL = 1000000000 !< Verbosity level for hac_init: verbose value, output from all tasks.
  INTEGER,PARAMETER,PUBLIC :: HAC_DEBUG = -1 !< Verbosity level for hac_init: debug value.
  INTEGER,PARAMETER,PUBLIC :: HAC_DEBUG_ALL = -1000000000 !< Verbosity level for hac_init: debug value, output from all tasks.
! ############## END OF PUBLIC DECLARATIONS ############################
  PRIVATE hac_adios_lib_exit !< Here to force Doxygen-1.8.4 to document correctly.
  INTEGER,PARAMETER :: AMFNL_ = 128  !< Maximal file name length
  INTEGER,PARAMETER :: MEVVL_ = 128  !< Maximal environment variable value length
  ! Note: adios_double_complex(=HAC_ADIOS_DOUBLE_COMPLEX) is missing or not PUBLIC !
  INTEGER,PARAMETER :: HAC_ADIOS_DOUBLE_COMPLEX = 11 !< (16) ! type code for double complex
  INTEGER,PARAMETER :: OFF_K = 8 !< KIND for long integers
  INTEGER,PARAMETER :: HAC_A_INT_T = adios_long !< ADIOS type for integers
  INTEGER,PARAMETER :: HAC_M_INT_T = MPI_INTEGER8 !< MPI type for integers
  INTEGER,PARAMETER :: TIME_K = 8 !< Time REAL kind
  INTEGER,PARAMETER :: AI_K=8 !< INTEGER kind for ADIOS interface
  INTEGER,PARAMETER :: K_ = 1024, M_ = K_*K_, MB_ = M_ !< KB, MB
  INTEGER :: A_SS = 4*MB_, A_BS = 4*MB_ !< ADIOS  Stripe Size, Block Size
  INTEGER,PARAMETER :: OKV_ = 0 !< No error value
  INTEGER,PARAMETER :: ERRV_G = -1 !< General error value
  INTEGER,PARAMETER :: ERRV_M = -2 !< Dimensions mismatch value
  INTEGER,PARAMETER :: ERRV_F = -3 !< Invoking after post-finalization error value
  INTEGER,PARAMETER :: ERRV_I = -4 !< Multiple finalization error value
  INTEGER,PARAMETER :: ERRV_D = -5 !< Dimensions Limits related error (MAXLD trespassed)
  INTEGER,PARAMETER :: MAXLD = K_*K_*K_ !< Maximum Local Dimension
  INTEGER,PARAMETER :: RR = 0 !< Root Rank
  INTEGER :: A_NA = 1024 !< ADIOS aggregators (output subfiles)
  INTEGER :: A_NO = 128, A_SC=16 !< Lustre OST's and stripe count
  INTEGER :: eanf = 0 !< Errors Are Non Fatal (if 1)
  CHARACTER(LEN=*),PARAMETER :: A_GN = "all" !< ADIOS group name
  CHARACTER(LEN=*),PARAMETER :: DAC = "dac"  !< dimensions array count
  INTEGER,PARAMETER :: EU = ERROR_UNIT, OU=OUTPUT_UNIT !< output units
  INTEGER,PARAMETER :: MDIMS_ = 7, VIL = 3, ADT_ = 3 !< Maximum Dimensions, Variable Identifier Length, Array Dimension Types
  INTEGER,PARAMETER :: OFF_ = 1, GLO_ = 2, LOC_ = 3 !< Offset, Global, Local
  INTEGER,PARAMETER :: MDS_ = 100 !< Max Dumpable local Size (for debug)
  INTEGER :: DODEBUG_ = 0 !< Quiet mode if DODEBUG_.LT.1
  CHARACTER(LEN=VIL),PARAMETER :: AID = "GAS" !< Global Array Slice
  INTEGER(KIND=OFF_K),PARAMETER :: NDIMS = HLST_HAC_USE_ADIM !< dimensions of the I/O array
  CHARACTER(LEN=*),PARAMETER :: BNR = "hlst_adios_checkpoint: " !< Banner
  CHARACTER(LEN=*),PARAMETER :: ARMP = "verbose=3"!read init method parms
  INTEGER :: verbose = 0 !< Quiet mode if .LT.1
  INTEGER :: atss, itss ! Array/Index type storage size
  REAL(KIND=TIME_K) :: wt, rt !< Write time, read time
  INTEGER :: a_buf_mb = 0 ! ADIOS buffer size in MB
  INTEGER(KIND=AI_K) :: agid, avid !< ADIOS Group ID, ADIOS Variable ID
  CHARACTER(LEN=AMFNL_) :: fhn, aps, awm!< First host name, ADIOS write parameter string, ADIOS write method
  INTEGER :: acomm = MPI_COMM_NULL !< Communicator for ADIOS operations
  INTEGER :: asize = -1 !< ADIOS communicator size
  INTEGER :: arank !< Rank within the ADIOS communicator
  INTEGER :: aerr = 0 !< ADIOS routines error variable
  INTEGER :: ierr = 0 !< Generic error variable
  INTEGER :: hac_a_num_t = adios_unknown !< ADIOS array numerical type
  LOGICAL :: initialized = .FALSE., finalized = .FALSE. !< helper vars
  LOGICAL :: initialized_vars = .FALSE. !< helper vars
  LOGICAL :: abao = .FALSE. !< Allocate Buffer At Once
  LOGICAL :: amroot = .FALSE. !< Is This Rank Root ?
 CONTAINS
!
  SUBROUTINE hac_herr
!> Error handler. It performs collective operations.
   IMPLICIT NONE
   ! CHARACTER(LEN=AMFNL_) :: aes = ''
   CHARACTER(LEN=*),PARAMETER :: BBC = " *** hlst_adios_checkpoint *** "
   CHARACTER(LEN=*),PARAMETER :: ENR = " *** "
   INTEGER :: gierr = 0, ferr = 0
!
   CALL MPI_Allreduce(ierr, gierr, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
   ierr = gierr 
   CALL MPI_Allreduce(aerr, gierr, 1, MPI_INTEGER, MPI_SUM, acomm, ierr)
   aerr = gierr 
!
   IF ( amroot ) THEN
    SELECT CASE (ierr)
     CASE (ERRV_M)
     WRITE(EU,*)BBC,'Read Data Dimensions Mismatch !',ENR
     CASE (ERRV_D)
     WRITE(EU,*)BBC,'Read Data Dimensions Off-Limits !',ENR
     CASE (ERRV_F)
     WRITE(EU,*)BBC,'More than one finalization is not allowed !',ENR
     CASE (ERRV_I)
     WRITE(EU,*)BBC,'Initializing after finalization not allowed !',ENR
     !CASE (errv_a)
     !WRITE(EU,*)'Generic Error!'
     ! MPI Error ?!
    END SELECT
   END IF
   IF (ierr .NE. OKV_) THEN
    IF ( amroot ) THEN
     WRITE(EU,*)BBC,'Detected Error Code: ',ierr,ENR
    END IF
    ferr = ierr
   END IF
!
   IF (aerr .LT. OKV_) THEN
    IF ( amroot ) THEN
     ! CALL adios_errmsg(aes, aesl) ! This is not part of the interface.
     WRITE(EU,*) BBC,'Detected an ADIOS Error: ',aerr,ENR
     ! The following are not PUBLIC members of adios_defs_mod, with no adios_*
     ! prefixed symbol (as of ADIOS-1.5).
     SELECT CASE (aerr)
      CASE (-1) 
       WRITE(EU,*) BBC,'Presumably a Memory Allocation Error!',ENR
      CASE (-2)
       WRITE(EU,*) BBC,'Presumably a File Open Error!',ENR
      CASE (-3)
       WRITE(EU,*) BBC,'Presumably a File Not Found Error!',ENR
      CASE (-24)
       WRITE(EU,*) BBC,'Presumably a Too Many Files Error!',ENR
      CASE DEFAULT
     END SELECT
    END IF
    ferr = aerr
   END IF
!
   IF ( eanf .EQ. 0 .AND. ferr .NE. 0 ) THEN
    IF ( amroot ) WRITE(EU,*)eanf
    IF ( amroot ) WRITE(EU,*)BBC,'Will abort job with code',ferr,'.',ENR
    CALL MPI_Abort(acomm, ferr, ierr)
   END IF
  END SUBROUTINE hac_herr
!
  SUBROUTINE hac_rerr
   aerr = OKV_
   ierr = OKV_
  END SUBROUTINE hac_rerr
!
  PURE SUBROUTINE hac_gadn ( vid, adi, adj)
!> ADIOS (array) Dimension(s) Name (string)
   IMPLICIT NONE
   CHARACTER(LEN=*),INTENT(INOUT) :: vid
   CHARACTER(LEN=VIL+1) :: vidt
   INTEGER,INTENT(IN) :: adi, adj
   INTEGER :: i,u,l
!
   IF (adi.LT.0) THEN
    l = 1
    u = - adi
   ELSE
    l = adi
    u = adi
   END IF
   WRITE (vid,'(a)') TRIM("")
   DO i = l, u
     IF (adj==OFF_) WRITE(vidt,'(a,i0)') TRIM("od"),i
     IF (adj==GLO_) WRITE(vidt,'(a,i0)') TRIM("gd"),i
     IF (adj==LOC_) WRITE(vidt,'(a,i0)') TRIM("ld"),i
     vid = TRIM(vid) // TRIM(vidt)
     IF (i.LT.u) vid = TRIM(vid) // TRIM(',')
   END DO
  END SUBROUTINE hac_gadn
!
#if 0
  SUBROUTINE hac_defv_ra (lds,gds,ods)
   IMPLICIT NONE
   CHARACTER(LEN=MDIMS_*(VIL+1)-1),INTENT(IN) :: lds,gds,ods
!
   CALL adios_define_var (agid, AID,"",hac_a_num_t,lds,gds,ods,avid)
  END SUBROUTINE hac_defv_ra
#endif
!
!< Define the necessary ADIOS variables.
!! Cannot be called before hac_init.
  SUBROUTINE hac_defv( A_NA_,A_NO_,A_SC_,A_SS_,A_BS_,ierr_ )
   IMPLICIT NONE
!
   INTEGER,INTENT(IN),OPTIONAL :: A_NA_,A_NO_,A_SC_,A_SS_,A_BS_,ierr_
!
   CHARACTER(LEN=VIL) :: vid
   INTEGER :: i,j
   CHARACTER(LEN=MDIMS_*(VIL+1)-1) :: lds,gds,ods
!
   IF ( initialized_vars .EQV. .FALSE. ) THEN
!
    HAC_SET_IF_PRESENT(A_NA_,A_NA)
    HAC_SET_IF_PRESENT(A_NO_,A_NO)
    HAC_SET_IF_PRESENT(A_SS_,A_SS)
    HAC_SET_IF_PRESENT(A_SC_,A_SC)
    HAC_SET_IF_PRESENT(A_BS_,A_BS)
!
    CALL adios_declare_group (agid, A_GN,"time0",1, aerr)
    aerr = 0 ! FIXME
    CALL hac_herr
    IF (.FALSE.) THEN
     A_NA = 16
     A_NO = 64
     A_SC = 16
     A_SS = 4*MB_
     A_BS = 4*MB_
     WRITE(aps,'(a,i0,a,i0,a,i0,a,i0,a,i0)') &
     &TRIM("num_aggregators="),A_NA,&
     &TRIM(";num_ost="),A_NO,&
     &TRIM(",stripe_count="),A_SC,&
     &TRIM(",stripe_size="),A_SS,&
     &TRIM(",block_size="),A_BS
     WRITE(awm,'(a)') "MPI_AGGREGATE"
    ELSE
     WRITE(aps,'(a,i0,a,i0)') &
     &TRIM("num_aggregators="),A_NA,&
     &TRIM(";num_ost="),A_NO
     WRITE(awm,'(a)') "MPI_AMR"
    END IF
     IF (DODEBUG_.NE.0.AND.amroot) THEN
      WRITE(OU,'(a,a,a)') BNR,'Using ADIOS method: ',TRIM(awm)
      WRITE(OU,'(a,a,a)') BNR,'Using ADIOS params: ',TRIM(aps)
     END IF
     CALL adios_select_method (agid,TRIM(awm),TRIM(aps),"",aerr)
     aerr = 0 ! FIXME
     CALL adios_define_var (agid,DAC,"",HAC_A_INT_T,"","","",avid)
     DO j = 1, ADT_
      DO i = 1, NDIMS
       CALL hac_gadn(vid,i,j)
       CALL adios_define_var (agid, vid,"",HAC_A_INT_T,"","","",avid)
      END DO
     END DO
     CALL hac_herr
     CALL hac_gadn(ods,INT(-NDIMS),OFF_)
     CALL hac_gadn(gds,INT(-NDIMS),GLO_)
     CALL hac_gadn(lds,INT(-NDIMS),LOC_)
#if HLST_HAC_USE_REAL
     IF (KIND(SAMPLEVAR).EQ.4) THEN
      hac_a_num_t = adios_real
     ELSE
      hac_a_num_t = adios_double
     ENDIF
#else
     IF (KIND(SAMPLEVAR).EQ.4) THEN
      hac_a_num_t = adios_complex
     ELSE
      hac_a_num_t = HAC_ADIOS_DOUBLE_COMPLEX
     ENDIF
#endif
     !
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,'(a,a,i0,a,a,a,a,a,a,a,a,a)') BNR,('on rank '),&
      & arank," aid:'",aid,"' lds:'",TRIM(lds),"' gds:'",TRIM(gds)&
      &,"' ods:'",TRIM(ods),"'"
     ENDIF
     !
     CALL adios_define_var (agid, AID,"",hac_a_num_t,lds,gds,ods,avid)
#if (HLST_HAC_EXTRA_METADATA==1)
     CALL adios_define_var (agid, "st","",adios_double,"","","",avid)
     CALL adios_define_var (agid, "ts","",adios_double,"","","",avid)
#endif
     ! call hac_defv_ra (lds,gds,ods)
     CALL hac_herr
     initialized_vars = .TRUE.
     HAC_SET_IF_PRESENT(ierr_,ierr)
    END IF
!
  END SUBROUTINE hac_defv
!
!> Initializes module hlst_adios_checkpoint.
!! It must be called once, after MPI_init.
!! 
!! In verbose mode, supported arrays type/dimensions information will be printed out.
!! 
  SUBROUTINE hac_init ( acomm_, verbose_, ierr_ )
   USE mpi, ONLY: MPI_Comm_dup 
   IMPLICIT NONE
!
INTEGER,OPTIONAL,INTENT(IN) :: acomm_ !< Checkpoint specific MPI communicator. If not provided, MPI_COMM_WORLD will be used.
   INTEGER,OPTIONAL,INTENT(IN) :: verbose_ !< Requested verbosity level. Default is not verbose (0).
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
!
   IF ( finalized .EQV. .TRUE. ) THEN
    ierr = ERRV_I
   END IF
!
   IF ( initialized .EQV. .FALSE. ) THEN
    IF ( PRESENT ( verbose_ ) ) THEN
     verbose = MAX(0, ABS(verbose_))
     IF ( verbose_ .LT. 0 ) DODEBUG_ = 1
    END IF
!
    IF ( PRESENT ( acomm_ ) ) THEN
     ! CALL MPI_Comm_dup (acomm_, acomm, ierr)
     acomm = acomm_
    ELSE
     CALL MPI_Comm_dup (MPI_COMM_WORLD, acomm, ierr)
     ! CALL hac_herr ! not allowed yet
    END IF
    CALL MPI_Comm_rank (acomm, arank, ierr)
    IF ( arank .EQ. RR ) amroot = .TRUE.
    CALL hac_herr
    CALL MPI_Comm_size (acomm, asize, ierr)
    CALL hac_herr
    CALL adios_init_noxml (acomm, aerr)
    CALL hac_herr
    a_buf_mb = 0 ! Will be allocated at the next call of adios_allocate_buffer
    CALL adios_allocate_buffer (a_buf_mb, aerr)
    CALL hac_rerr ! FIXME: resetting aerr. This is tolerated by ADIOS-1.5.
    CALL hac_herr
    ! CALL hac_defv
    itss = hac_isz()
    atss = hac_asz()
    eanf = 0 ! errors are fatal
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD 
!#ifdef __GFORTRAN__
     ! It seems like ifort stands this. Compatibility mode by default !?
     IF ( verbose .GT. 0 ) THEN
      CALL HOSTNM(fhn)
      fhn = TRIM (fhn)
#if (HLST_HAC_USE_REAL==1)
      IF (amroot) WRITE(OU,'(a,a,i0,a,i0,a)')BNR,'Using REAL*',atss ,&
#else
      IF (amroot) WRITE(OU,'(a,a,i0,a,i0,a)')BNR,'Using COMPLEX*',atss,&
#endif
              &' with ', HLST_HAC_USE_ADIM,' dimensions as array type.'
      IF (amroot) WRITE(OU,'(a,a,a,a)') BNR,'First  task on host ',&
              &TRIM(fhn),'.'
      IF (arank==1) WRITE(OU,'(a,a,a,a)') BNR,'Second task on host ',&
              &TRIM(fhn),'.'
      IF (asize.GT.1.AND.arank.EQ.asize-1) WRITE(OU,'(a,a,i0,a,a,a)')&
              & BNR,'Last task of ',asize,' on host  ',TRIM(fhn),'.'
      IF ( amroot .AND. DODEBUG_ .NE. 0 ) WRITE(OU,'(a,a,i0,a,i0)')BNR,&
              & 'verbosity level: ',verbose_,', debug level:  ',dodebug_
     END IF
!#endif
#endif
    initialized = .TRUE.
    finalized = .FALSE.
    ierr = OKV_
    aerr = OKV_
   END IF
   HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_init
!
!!> \internal
  SUBROUTINE hac_adios_lib_exit
   IF ( finalized .EQV. .FALSE. ) THEN
    CALL adios_finalize (arank, ierr)
    CALL hac_herr
    finalized = .TRUE.
   ELSE
    ierr = ERRV_F
   END IF
   CALL hac_herr
  END SUBROUTINE hac_adios_lib_exit
!
  INTEGER FUNCTION check_dda(LBA)
   IMPLICIT NONE
   INTEGER :: ierr = OKV_
   INTEGER,INTENT(IN),DIMENSION(:) :: LBA
!
   IF ( PRODUCT (INT8(LBA(1:NDIMS)) ) .EQ. 0) THEN
    ierr = ERRV_M
   END IF
!
   IF ( MINVAL (LBA(1:NDIMS) ) .GT. 1) THEN
    ierr = ERRV_M
   END IF
!
   IF ( MAXVAL (LBA(1:NDIMS) ) .GT. MAXLD ) THEN
    ierr = ERRV_D
   END IF
!
   check_dda = ierr
  END FUNCTION
!
  SUBROUTINE check_ada(GASS, ada)
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: GASS(:) !< Global Array Slice Shape
   INTEGER(KIND=OFF_K),INTENT(IN) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
   INTEGER :: lerr = OKV_
!
   lerr = check_dda(INT(ada(1:NDIMS, OFF_)))
   IF ( PRODUCT (ada(1:NDIMS, LOC_) ) .EQ. 0) THEN
    lerr = ERRV_M
   END IF
   IF ( PRODUCT (ada(1:NDIMS, GLO_) ) .EQ. 0) THEN
    lerr = ERRV_M
   END IF
   IF (PRODUCT (ada(1:NDIMS,LOC_)) .GT. PRODUCT(ada(1:NDIMS,GLO_))) THEN
    lerr = ERRV_M
   END IF
   IF ( SUM (GASS(1:NDIMS) - ada(1:NDIMS, LOC_) ) .NE. 0) THEN
    ierr = ERRV_M
   END IF
   ierr = lerr
  END SUBROUTINE check_ada
!
INTEGER FUNCTION hac_isz()
  IMPLICIT NONE
  INTEGER :: sz
!
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD
     sz = INT(SIZEOF(NDIMS))
#else
     ! This Requires Fortran-2008
     sz = STORAGE_SIZE(NDIMS)/ 8
#endif
   hac_isz = sz
  END FUNCTION hac_isz
!
INTEGER FUNCTION hac_asz()
  IMPLICIT NONE
  INTEGER :: sz
  HAC_NTYPE,PARAMETER :: XYE = 0.0 !< A sample value.
!
#if HLST_HAC_COMPILING_WITH_GNU_STANDARD
     sz = INT(SIZEOF(XYE))
#else
     ! This Requires Fortran-2008
     sz = STORAGE_SIZE(XYE)  / 8
#endif
   hac_asz = sz
  END FUNCTION hac_asz
!
!> Reads header information from the checkpoint file.
!! Can be called repeatedly to get information out of a checkpoint file.
  SUBROUTINE hac_info(fname, ierr_, GADA_, ntc_, ts_, st_)
   IMPLICIT NONE
!
   INTEGER,INTENT(OUT),OPTIONAL :: GADA_(:) !< Global Array's Dimensions Array. 
!! The user can use a subset of these dimensions to ALLOCATE the right sized (or, SHAPE'd) local array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   CHARACTER,OPTIONAL,INTENT(OUT) :: ntc_ !< Numerical type character. Will be either S,D,C,Z.
#if (HLST_HAC_EXTRA_METADATA==1)
   REAL(KIND=8),OPTIONAL,INTENT(OUT) :: ts_ !< Meta-data scalar: Time Step. Can be useful to certain users.
   REAL(KIND=8),OPTIONAL,INTENT(OUT) :: st_ !< Meta-data scalar: Simulation time. Can be useful to certain users.
#endif
!
   CHARACTER :: ntc !< Numerical type
   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
   INTEGER :: i,j
   INTEGER(KIND=AI_K) :: fh ! 
   CHARACTER(LEN=VIL) :: vid
   INTEGER(KIND=OFF_K) :: ndims_v
#if (HLST_HAC_EXTRA_METADATA==1)
   REAL(KIND=8) :: ts = 0.0, st = 0.0 ! DOUBLE PRECISION is not supported by ADIOS
#endif
!
   IF (finalized) GOTO 9999
!
   CALL hac_defv
   rt = - MPI_Wtime()
   CALL adios_read_init_method (ADIOS_READ_METHOD_BP, acomm, ARMP, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, DAC, ndims_v, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( ndims_v .GT. NDIMS ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( DODEBUG_ .GT. arank ) THEN
    WRITE (OU,*) "READ ",DAC,": ",NDIMS," on ", NDIMS
   END IF
   !IF ( PRODUCT (SHAPE (GAS) ) .EQ. 0 ) THEN
   ! ierr = ERRV_M
   !END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = 1 ! LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(INT(ada(1:NDIMS,OFF_))) ! FIXME: overflow danger
   !ierr = check_dda(LBA)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   !CALL check_ada(SHAPE(GAS), ada)
   DO j = 1, ADT_
    DO i = 1, NDIMS
     IF (j.EQ.OFF_) CYCLE ! this variable should be locally set: read values won't match
     CALL hac_gadn (vid,i,j)
     CALL adios_get_scalar (fh, vid, ada(i,j), ierr)
     HLST_HAC_GOTO_END_ON_ERROR(9999)
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,*) "READ ",vid,": ",ada(i,j)," on ", arank
     END IF
    END DO
   END  DO
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = 1 ! LBA(1:NDIMS) + 1 - ozb
#if (HLST_HAC_EXTRA_METADATA==1)
   CALL adios_get_scalar (fh, "ts", ts, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, "st", st, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
#endif
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF (.NOT.PRESENT(GADA_)) THEN
    ! WRITE (OU,*) "GADA ",ada(1:NDIMS,GLO_)
   ELSE
    GADA_(1:NDIMS) = INT(ada(1:NDIMS,GLO_)) ! FIXME: conversion danger + overflow (no check on GADA_)
   END IF
!
   ntc = '?'
   SELECT CASE (hac_asz())
#if (HLST_HAC_USE_REAL==1)
    CASE (4)
     ntc = 'S'
    CASE (8)
     ntc = 'D'
#else
    CASE (8)
     ntc = 'C'
    CASE (16)
     ntc = 'Z'
#endif
   END SELECT
!
   IF (PRESENT(ntc_)) THEN
    ntc_ = ntc
   ELSE
    ! WRITE (OU,*) "TYPE: ",ntc, hac_asz()
   END IF
!
   IF (PRESENT(ierr_)) THEN
    ierr_ = ierr
   END IF
#if (HLST_HAC_EXTRA_METADATA==1)
   IF (PRESENT(ts_)) THEN
    ts_ = ts
   END IF
   IF (PRESENT(st_)) THEN
    st_ = st
   END IF
#endif
!
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_info
!
!> Prints Read/Write Statistics
  SUBROUTINE hac_prtsts(lcib,row,et,fname)
   INTEGER(KIND=OFF_K),INTENT(IN) :: lcib ! local contribution in bytes
   LOGICAL, INTENT(IN) :: row ! read (.TRUE.) or write (.FALSE.)
   ! CHARACTER(LEN=*),INTENT(IN) :: ms ! message string
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   REAL(KIND=TIME_K),INTENT(IN) :: et !< Elapsed time
!
   INTEGER(KIND=OFF_K) :: tcib,tcim ! local contribution in bytes/megabytes
!
   tcib = 0
   CALL MPI_Allreduce( lcib, tcib, 1, HAC_M_INT_T, MPI_SUM, acomm, ierr)
   CALL hac_herr
!
   tcim = tcib / MB_
   IF ( verbose.GT.0 .AND. amroot ) THEN
    IF(row)&
    &WRITE(OU,'(a,a,i0,a,f10.3,a,f10.2,a,a)') BNR,"Read    ",&
             &tcim," MB in ",et," s, at ",(1.0/et)*tcim," MB/s from ",&
             &TRIM(fname)
!
    IF(.NOT.row)&
    &WRITE(OU,'(a,a,i0,a,f10.3,a,f10.2,a,a)') BNR,&
     &"Written ",tcim," MB in ",wt," s, at ",(1.0/et)*tcim,&
     &" MB/s  to  ",TRIM(fname)
   END IF
!
  END SUBROUTINE hac_prtsts
!
 SUBROUTINE hac_wrt_sda(TB,TDA,TA)
  INTEGER,INTENT(IN) :: TDA(:) !
  CHARACTER(LEN=*),INTENT(IN) :: TB,TA ! before, after
!
   SELECT CASE (SIZE(TDA))
    CASE (2) 
     WRITE(*,'(a,"[",&
      &i0,",",i0,"]",a)')TB,&
     &TDA(1),TDA(2),TA
    CASE (4) 
     WRITE(*,'(a,"[",&
      &i0,",",i0,",",i0,",",i0,"]",a)')TB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TA
    CASE (6) 
     WRITE(*,'(a,"[",&
      &i0,",",i0,",",i0,",",i0,",",i0,",",i0,"]",a)')TB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),TA
    CASE (8) 
     WRITE(*,'(a,"[",&
      &i0,":",i0,",",i0,":",i0,",",i0,":",i0,",",i0,":",i0,"]",a)')TB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),TDA(7),TDA(8),TA
    CASE (12) 
     WRITE(*,'(a,"[",&
      &i0,":",i0,",",i0,":",i0,",",i0,":",i0,",",i0,":",i0,","&
     &,i0,":",i0,",",i0,":",i0,"]",a)')TB,&
     &TDA(1),TDA(2),TDA(3),TDA(4),TDA(5),TDA(6),TDA(7),&
     &TDA(8),TDA(9),TDA(10),TDA(11),TDA(12),TA
   END SELECT
 END SUBROUTINE hac_wrt_sda
!
!
!> Reads the distributed array from the checkpoint file.
!! An arbitrary array slice can be read from whatever rank.
!! For this reason unlike hac_write, here the index base cannot be unambiguously 
!! inferred from the input, and it is advised to specify it.
!!
  SUBROUTINE hac_read(GAS, LBA, fname, LDA_, ozb_, ierr_)
   IMPLICIT NONE
!
   HAC_NTYPE,HAC_DSPEC,INTENT(INOUT) :: GAS !< Global Array Slice. 
!! \c SHAPE(GAS) will be assumed to be local slice dimensions; LBOUND(GAS) its global lower bounds, 
!! unless \a LBA is provided. Expected to be already allocated.
   INTEGER,INTENT(IN),DIMENSION(:) :: LBA !< Lower Bounds Array. Either 1- or 0-based, on *all* the dimensions. 
!! Offset of the local slice in GAS with respect to the global array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(IN),DIMENSION(:) :: LDA_ !< Local Dimensions Array. 
!! Dimensions of the local slice in GAS. If omitted, SHAPE(GAS) will be used instead.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   INTEGER,OPTIONAL,INTENT(IN) :: ozb_ !< One or Zero index Base for the global array: 
!! can be either 1 or 0; if not specified, this will be taken to be the minimum value across LBA instances.
!
   CHARACTER(LEN=VIL) :: vid
   INTEGER :: i,j
   INTEGER(KIND=AI_K),DIMENSION(MDIMS_) :: off = -1, rds = -1!offset, read size
   INTEGER(KIND=AI_K) :: sel  ! ADIOS selection object
   INTEGER(KIND=AI_K) :: lcib  ! local contribution in bytes
   INTEGER(KIND=AI_K) :: fh ! 
   INTEGER(KIND=OFF_K) :: ndims_v
   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_) !< Array's Dimensions Array
   INTEGER(KIND=OFF_K) :: lmo, ozb = 0 ! Local Minimum Offset, One or Zero Base ?
!
   IF (finalized) GOTO 9999
!
   IF (PRESENT(OZB_)) THEN
    lmo = ozb_ ! Either 1 or 0; FIXME: no check here.
    ELSE 
    lmo = MINVAL(LBA(1:NDIMS) )
   END IF
   lmo = MIN(1_OFF_K,MAX(lmo,0_OFF_K))
   CALL MPI_Allreduce( lmo, ozb, 1, HAC_M_INT_T, MPI_MAX, acomm, ierr )
   CALL hac_defv
   rt = - MPI_Wtime()
   CALL adios_read_init_method (ADIOS_READ_METHOD_BP, acomm, ARMP, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_get_scalar (fh, DAC, ndims_v, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( ndims_v .GT. NDIMS ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF ( DODEBUG_ .GT. arank ) THEN
    WRITE (OU,*) "READ ",DAC,": ",NDIMS," on ", NDIMS
   END IF
   IF ( PRODUCT (SHAPE (GAS) ) .EQ. 0 ) THEN
    ierr = ERRV_M
   END IF
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(INT(ada(1:NDIMS,OFF_))) ! FIXME: overflow danger
   !ierr = check_dda(LBA)
   CALL hac_herr
   CALL check_ada(SHAPE(GAS), ada)
   ada(1:NDIMS,LOC_) = SHAPE(GAS) ! The local SHAPE will determine how much we read
   DO j = 1, ADT_
    DO i = 1, NDIMS
     IF (j.EQ.OFF_) CYCLE ! this variable should be locally set: read values won't match
     IF (j.EQ.LOC_) CYCLE ! this variable should be locally set: depends on SHAPE(GAS)
     CALL hac_gadn (vid,i,j)
     CALL adios_get_scalar (fh, vid, ada(i,j), ierr)
     HLST_HAC_GOTO_END_ON_ERROR(9999)
     IF ( DODEBUG_ .GT. arank ) THEN
      WRITE (OU,*) "READ ",vid,": ",ada(i,j)," on ", arank
     END IF
    END DO
   END  DO
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
!
   CALL check_ada(SHAPE(GAS), ada)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_open_file (fh,fname,ADIOS_READ_METHOD_BP,acomm,ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   off(1:NDIMS) = ada(1:NDIMS, OFF_) - 1
   IF (PRESENT(LDA_)) THEN
    IF ( dodebug_.GT.0 .AND. amroot ) &
     &WRITE(OU,*)'Will read custom shape:', LDA_
    ada(1:NDIMS, LOC_) = LDA_(1:NDIMS) ! Overriding SHAPE(GAS)
   END IF
   rds(1:NDIMS) = ada(1:NDIMS, LOC_)
!
   CALL adios_selection_boundingbox (sel, INT(NDIMS), off, rds)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_schedule_read (fh, sel, TRIM("/" // AID), 0, 1, GAS, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_perform_reads (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_close (fh, ierr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_selection_delete (sel)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_read_finalize_method (ADIOS_READ_METHOD_BP, ierr)
!
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_ ) THEN
    WRITE (OU,*) "READ ",AID,": ",GAS," on ", arank
   END IF
!
   rt = rt + MPI_Wtime()
   lcib = itss * (1 + NDIMS * ADT_) + atss * PRODUCT(ada(1:NDIMS,LOC_))
!     
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_ ) THEN
    IF ( amroot ) WRITE (OU,*)"SIZEOF(NDIMS):", itss
    IF ( amroot ) WRITE (OU,*)"SIZEOF(TYPE):", atss
    ! IF ( amroot ) WRITE (OU,*)"SIZEOF(ada(1,1)):", itss
   END IF
!
   CALL hac_prtsts(lcib,.TRUE.,rt,TRIM(fname))
!
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_read
!
!> Writes the global array to the checkpoint file.
  SUBROUTINE hac_write (GAS, LBA, fname, ts_, st_, a_buf_mb_, ierr_)
   IMPLICIT NONE
!
   HAC_NTYPE,HAC_DSPEC,INTENT(IN) :: GAS !< Global Array Slice. 
!! \c SHAPE(GAS) will be assumed to be the local slice dimensions.
   INTEGER,INTENT(IN),DIMENSION(:) :: LBA !< Lower Bounds Array. 
!! Either 1- or 0-based, on *all* the dimensions. Offset of the local slice in GAS with respect to the global array.
   CHARACTER(LEN=*),INTENT(IN) :: fname !< Path to the checkpoint file. Any length is allowed.
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
   INTEGER,OPTIONAL,INTENT(IN) :: a_buf_mb_ !< ADIOS buffer size in MB. 
!! If omitted, it will be automatically set at the first write. If you foresee a subsequent write to be larger, 
!! submit an upper limit here at the first call.
#ifdef HLST_HAC_EXTRA_METADATA
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: ts_ !< Meta-data scalar: Time Step. Can be useful to certain users.
   REAL(KIND=8),OPTIONAL,INTENT(IN) :: st_ !< Meta-data scalar: Simulation time. Can be useful to certain users.
!
   REAL(KIND=8) :: ts = 0.0, st = 0.0 ! DOUBLE PRECISION is not supported by ADIOS
#endif
!
   INTEGER(KIND=AI_K) :: ah,ats,ags=0 !ADIOS handle, total size, group size
   INTEGER :: i,j
   CHARACTER(LEN=VIL) :: vid
   INTEGER(KIND=OFF_K) :: ada(MDIMS_,ADT_), adas(MDIMS_) !< Array's Dimensions Array (and send buffer)
   INTEGER(KIND=OFF_K) :: lmo, ozb = 0 ! Local Minimum Offset, One or Zero Base ?
!
   IF(finalized) GOTO 9999
   wt = - MPI_Wtime()
   lmo = MINVAL(LBA(1:NDIMS) )
   CALL MPI_Allreduce( lmo, ozb, 1, HAC_M_INT_T, MPI_MAX, acomm, ierr )
!
   ada(1:NDIMS,OFF_) = LBA(1:NDIMS) + 1 - ozb
   ierr = check_dda(INT(ada(1:NDIMS,OFF_))) ! FIXME: overflow danger
   ! ierr = check_dda(LBA)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL hac_defv
   CALL adios_open (ah, A_GN, TRIM(fname), "w", acomm, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ada(:,:) = 1
   ada(1:NDIMS, OFF_) = LBA(1:NDIMS) + 1 - ozb
   ada(1:NDIMS, LOC_) = SHAPE(GAS)
   ada(1:NDIMS, GLO_) = ada(1:NDIMS,OFF_) + ada(1:NDIMS,LOC_) - 1 
   adas(1:NDIMS) = ada(1:NDIMS,GLO_)
   CALL MPI_Allreduce( ada(1:NDIMS,GLO_), adas(1:NDIMS), NDIMS,&
            & HAC_M_INT_T, MPI_MAX, acomm, ierr )
   ada(1:NDIMS,GLO_) = adas(1:NDIMS)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   ags = itss * (1 + NDIMS * ADT_) + atss * PRODUCT( ada(1:NDIMS, LOC_))
!
   IF ( DODEBUG_ .GT. arank ) THEN
    IF ( amroot )WRITE (OU,*)"WRITE NDIMS: ",NDIMS," ags: ",&
            & ags, " PROD(dims):", PRODUCT(ada(1:NDIMS,LOC_)),&
            & ags, " offsets:", (ada(1:NDIMS,OFF_)),&
            & ags, " gl.dims:", (ada(1:NDIMS,GLO_))
   END IF
   IF (abao.EQV..FALSE.) THEN
    IF (PRESENT(a_buf_mb_)) THEN
      a_buf_mb = INT( 1 + ( a_buf_mb_ ) ) ! The extra 1 is to round up.
    ELSE
      a_buf_mb = INT( 1 + ( ags / MB_ ) )
    END IF
    CALL adios_allocate_buffer (a_buf_mb, aerr)
    ! ADIOS-1.5 allows adios_allocate_buffer only once.
    abao = .TRUE.
   END IF
   CALL adios_group_size (ah, ags, ats, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_write (ah, DAC, NDIMS, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL check_ada(SHAPE(GAS), ada)
   ada(:,OFF_) = ada(:,OFF_) - 1 ! for ADIOS
   DO j = 1, ADT_
    DO i = 1, NDIMS
     CALL hac_gadn(vid,i,j)
     CALL adios_write (ah, vid, ada(i,j), aerr)
     HLST_HAC_GOTO_END_ON_ERROR(9999)
    END DO
   END DO
#if (HLST_HAC_EXTRA_METADATA==1)
   IF (PRESENT(ts_)) ts = ts_
   IF (PRESENT(st_)) st = st_
   CALL adios_write (ah, "ts", ts, aerr)
   CALL adios_write (ah, "st", st, aerr)
#endif
   ada(:,OFF_) = ada(:,OFF_) + 1 ! for ADIOS
!
   CALL adios_write (ah, AID, GAS, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   CALL adios_close (ah, aerr)
   HLST_HAC_GOTO_END_ON_ERROR(9999)
   IF (DODEBUG_.GT.arank.AND.PRODUCT(ada(:,LOC_))*asize.LT.MDS_) THEN
    WRITE (OU,*) "WRITE ",AID,": ",GAS," on ", arank
   END IF
!
   wt = wt + MPI_Wtime()
   CALL hac_prtsts(ags,.FALSE.,wt,TRIM(fname))
!
   IF ( verbose.GT.0 .AND. amroot ) THEN
   END IF
9999 HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_write
!
!> Finalizes module hlst_adios_checkpoint.
!! It must be called once, before MPI_Finalize.
!! \note #hac_exit reports successful operation as default.
  SUBROUTINE hac_exit(ierr_)
    IMPLICIT NONE
!
   INTEGER,OPTIONAL,INTENT(OUT) :: ierr_ !< Error status variable. If omitted, errors will be fatal.
!
    CALL hac_adios_lib_exit
    IF ( PRESENT(ierr_) ) THEN
     ierr_ = OKV_
    END IF
    HLST_HAC_HANDLE_IERR(ierr_,ierr)
  END SUBROUTINE hac_exit
!
!> Gets an INTEGER environment variable value.
INTEGER FUNCTION get_env_opt_i(name,dval_)
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: name
  CHARACTER(LEN=MEVVL_) :: value
  INTEGER :: i, status = 0
  INTEGER,INTENT(IN),OPTIONAL :: dval_
!
  ! CALL GETENV(name,value) 
  CALL GET_ENVIRONMENT_VARIABLE(name,value,status=status) 
  READ (value,'(i10)') i
  IF ( PRESENT(dval_) ) THEN
    !IF (status.EQ.1 .AND. i .EQ. 0 ) i = dval_
    IF (status.EQ.1 ) i = dval_ ! if var unset, override
    ! IF (i==0) i = dval_
  END IF
  get_env_opt_i = i
END FUNCTION get_env_opt_i
!
! Single write and read.
SUBROUTINE hac_laa(AS,LB,UB,erank)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: LB(NDIMS) !< Lower Bounds
  INTEGER,INTENT(IN) :: UB(NDIMS) !< Upper Bounds
  INTEGER,INTENT(IN) :: erank !
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE,INTENT(INOUT) :: AS !< Array Slice
!
#if   (HLST_HAC_USE_ADIM==6)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)&
           &,LB(4):UB(4),LB(5):UB(5),LB(6):UB(6)))
#elif   (HLST_HAC_USE_ADIM==4)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)&
           &,LB(4):UB(4)))
#elif   (HLST_HAC_USE_ADIM==2)
 ALLOCATE(AS(LB(1):UB(1),LB(2):UB(2)))
#else
#error !
#endif
#if HLST_HAC_USE_REAL
   AS =        1000*(erank + 1)
#else
   AS = CMPLX( 1000*(erank + 1), erank + 1 )
#endif
!
END SUBROUTINE hac_laa
!
! Single write and read.
 SUBROUTINE hac_swar(ecomm,LB,UB)
  IMPLICIT NONE
!
  INTEGER,INTENT(INOUT) :: ecomm
  INTEGER,INTENT(IN) :: LB(HLST_HAC_USE_ADIM)
  INTEGER,INTENT(IN) :: UB(HLST_HAC_USE_ADIM)
!
  INTEGER :: GADA(NDIMS) !< Global Array's Dimensions Array
  CHARACTER :: ntc
  CHARACTER(LEN=AMFNL_) :: fname = 'checkpoint.dat'
  INTEGER :: erank
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE :: GAS,GCP !< Global Array Slice
!
  CALL MPI_Comm_rank (ecomm, erank, ierr)
  fname = 'checkpoint.dat'
!
  CALL hac_laa(GAS,LB,UB,erank)
  CALL hac_laa(GCP,LB,UB,erank)
  ! CALL hac_laa(GAS,(/lb(1),(ib,i=2,NDIMS)/),ub-lb+1,erank)
  ! CALL hac_laa(GCP,(/lb(1),(ib,i=2,NDIMS)/),ub-lb+1,erank)
!
  GCP = GAS
  IF ( verbose.GT.erank ) THEN
   CALL hac_wrt_sda(BNR//'Local Array Dimensions    :',SHAPE(GAS),'')
  END IF
  CALL hac_write (GAS, LBOUND(GAS), TRIM(fname)) ! FIXME: shall test with writing st, ts
  GAS = 0
  CALL hac_info (fname, ierr, GADA_=GADA )
  IF ( verbose.GT.erank ) THEN
   CALL hac_wrt_sda(BNR//'Global Array Dimensions: ',GADA,'')
   CALL hac_wrt_sda(BNR//'Global Array Local Bounds: ',LBOUND(GAS),'')
  END IF
  CALL hac_info (fname, ierr, ntc_=ntc)
  IF ( verbose.GT.erank .AND. amroot ) THEN
   WRITE(OU,'(a,a)') BNR//'Numerical Type Character: ',NTC
  END IF
  CALL hac_info (fname, ierr)
  CALL hac_read (GAS, LBOUND(GAS), TRIM(fname))
  IF (SUM(GAS-GCP).NE.0.0) STOP "Critical error in restart !"
  DEALLOCATE(GAS)
  DEALLOCATE(GCP)
 END SUBROUTINE hac_swar
!
!> Performs basic error testing.
 SUBROUTINE hac_ex_ckp_err_test
  IMPLICIT NONE
!
  eanf = 1 ! 
!
  ! CALL hac_info ('/dev/zero', ierr) ! This crashes ADIOS. No workaround.
  ! CALL hac_info ('/dev/non-existent-file', ierr) ! This does not crash ADIOS but leads it to an inconsistent state.
  ierr = OKV_
  CALL hac_rerr
  eanf = 1 - eanf
!
 END SUBROUTINE hac_ex_ckp_err_test
!
!> Performs a basic correctness test run.
!! This is a complete MPI program and can only be called once, e.g.: in an empty program body.
 SUBROUTINE hac_test
  IMPLICIT NONE
  INTEGER :: loff = 0, losz = 0, hac_mb, rep, reps=10
  INTEGER :: ib ! Index Base
  INTEGER :: miub = 0, maub = 2, acub ! minimal/maximal/actual unbalance
  INTEGER :: ti ! test iterations
  INTEGER :: erank, ecomm
  HAC_NTYPE,HAC_DSPEC,ALLOCATABLE :: GAS,GCP !< Global Array Slice
#if HLST_HAC_USE_REAL
  INTEGER,PARAMETER :: HAC_1MB = 1024*1024/(1*HAC_KIND)
#else
  INTEGER,PARAMETER :: HAC_1MB = 1024*1024/(2*HAC_KIND)
#endif
  INTEGER :: tdim, ld(HLST_HAC_USE_ADIM)
  CHARACTER(AMFNL_) :: fname = 'checkpoint.dat'
  CHARACTER :: ntc
  INTEGER :: GADA(NDIMS) !< Global Array's Dimensions Array
  INTEGER :: i
  INTEGER :: LB(NDIMS), UB(NDIMS)
!
  ecomm = MPI_COMM_WORLD
  CALL MPI_Init (ierr)
  CALL MPI_Comm_rank (ecomm, erank, ierr)
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a)')'+--------------------------------------------------+'
   WRITE(OU,'(a)')'| MODULE hlst_adios_checkpoint test program        |'
   WRITE(OU,'(a)')'+--------------------------------------------------+'
   WRITE(OU,'(a)')'| # To influence it, set the HAC_MB env. var.:     |'
   WRITE(OU,'(a)')'| export HAC_MB=1  # if>0,  size of I/O tests      |'
   WRITE(OU,'(a)')'| export HAC_MB=-1 # if<=1, auto choice            |'
   WRITE(OU,'(a)')'| export HAC_ST=0 # stress test; default is 1      |'
   WRITE(OU,'(a)')'| export HAC_ET=0 # error test; default is 1       |'
   WRITE(OU,'(a)')'| export HAC_VL=0 # verbosity level; default is 1  |'
   WRITE(OU,'(a)')'| # you can increase it, negative values will      |'
   WRITE(OU,'(a)')'| # request extra (debug) output.                  |'
   WRITE(OU,'(a)')'| # And additionally (ADIOS parameters--internals):|'
   WRITE(OU,'(a)')'| export HAC_A_NA=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_NO=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_SC=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_SS=... #                            |'
   WRITE(OU,'(a)')'| export HAC_A_BS=... #                            |'
   WRITE(OU,'(a)')'+--------------------------------------------------+'
  END IF
!
  CALL hac_init (ecomm, get_env_opt_i("HAC_VL",HAC_VERBOSE), ierr)
!
  CALL hac_defv( &
         & A_NA_ = get_env_opt_i("HAC_A_NA",A_NA),&
         & A_NO_ = get_env_opt_i("HAC_A_NO",A_NO),&
         & A_SC_ = get_env_opt_i("HAC_A_SC",A_SC),&
         & A_SS_ = get_env_opt_i("HAC_A_SS",A_SS),&
         & A_BS_ = get_env_opt_i("HAC_A_BS",A_BS),&
         & ierr_=ierr)
!
  CALL hac_herr
  ! HLST_HAC_GOTO_END_ON_ERROR(9999)
!
  hac_mb = get_env_opt_i("HAC_MB",1)
  IF (hac_mb .LE. 0) THEN
   hac_mb = 1 ! Start here to improve the 'auto choice' option.
  END IF
  losz = hac_mb * HAC_1MB
  fname = 'checkpoint.dat'
  IF ( erank .EQ. 0 ) &
   & WRITE(OU,'(a,i0)') BNR//'BEGIN test with HAC_MB = ',hac_mb
  DO ib = 0, 1
  DO tdim = 2, HLST_HAC_USE_ADIM
   ld = 1 + ib - 1
   loff = erank + ib
   ld(tdim) = losz + ib - 1 ! 
   ld(1) = loff             !
!
   CALL hac_laa(GAS,(/ld(1),(ib,i=2,NDIMS)/),ld,erank)
   CALL hac_laa(GCP,(/ld(1),(ib,i=2,NDIMS)/),ld,erank)
!
   GCP = GAS
   CALL hac_write (GAS, LBOUND(GAS), TRIM(fname), st_=1.0_8, ts_=1.0_8) !
   GAS = 0
   CALL hac_info (fname, ierr, GADA_=GADA )
   ! IF (erank .EQ. 0) WRITE(OU,*)'GADA ',GADA
   IF(erank.EQ.0) THEN
     CALL hac_wrt_sda(BNR//'Global Array Dimensions:',GADA,'')
     CALL hac_wrt_sda(BNR//'Global Array Local Bounds: ',LBOUND(GAS),'')
   END IF
   CALL hac_info (fname, ierr, ntc_=ntc)
   IF (erank.EQ.0)WRITE(OU,'(a,a)')BNR//'Numerical Type Character: ',NTC
   CALL hac_info (fname, ierr)
   CALL hac_read (GAS, LBOUND(GAS), TRIM(fname))
   CALL hac_info (fname, ierr)
   CALL hac_read (GAS, LBOUND(GAS), TRIM(fname), ozb_=ib)
   IF (SUM(GAS-GCP).NE.0.0) STOP "Critical error in restart !"
   DEALLOCATE(GAS)
   DEALLOCATE(GCP)
  END DO
  END DO
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a)')BNR//'Stress testing...'
   WRITE(OU,'(a)')BNR//'(may be damaging for a mechanical hard drive).'
  END IF
  ! verbose = HAC_QUIET
  ! verbose = HAC_VERBOSE
  ! verbose = 3
!
  IF ( get_env_opt_i("HAC_ST",1) .NE. 0 ) THEN
!
  ti = 0
  DO rep = 1, reps
  DO tdim = 1, HLST_HAC_USE_ADIM
  DO losz = 1, 4 !
  DO acub = miub, maub !
  DO ib = 0, 1
   ! All dimensions are replicated and fixed.
   LB = (/(ib,     i=1,NDIMS)/)
   UB = (/(ib+losz,i=1,NDIMS)/)
   ! ... except this one which is distributed and if acub.GT.0, skewed.
   LB(tdim) = ib + (erank+0) * losz * (1+acub * (erank+0))
   UB(tdim) = ib + (erank+1) * losz * (1+acub * (erank+1)) - 1
   CALL hac_swar(ecomm,LB,UB)
   IF (verbose.GT.0) THEN
    ! IF (erank.EQ.0) WRITE(OU,*)' ',losz,'/',100,' ib=',ib
   END IF
   ti = ti + 1
  END DO
  END DO
  END DO
  END DO
!  IF(erank.EQ.0.AND.MOD(rep,reps/10).EQ.0)WRITE(OU,'(a,i0,a,i0,a,i0)')&
!  'at ',rep,' of ',reps,', sz=',losz
  END DO
  END IF
!
  IF (erank .EQ. 0) THEN
   WRITE(OU,'(a,i0,a)')&
    &BNR//'Stress testing done: ',ti,' tries successful.'
   WRITE(OU,'(a)')BNR//'END Benchmarking'
  END IF
!
  IF ( get_env_opt_i("HAC_ET",1) .NE. 0 ) THEN
   WRITE(OU,'(a)')BNR//'BEGIN Error test.'
   CALL hac_ex_ckp_err_test
   WRITE(OU,'(a)')BNR//'END Error test.'
  END IF
!
  CALL hac_exit(ierr)
  CALL MPI_Finalize (ierr)
 END SUBROUTINE hac_test
!
#undef HAC_FIELD
#undef HAC_KSPEC
#undef HAC_NTYPE
#undef HLST_HAC_USE_ADIM 
#undef HLST_HAC_USE_CMSK 
!
END MODULE hlst_adios_checkpoint
!
