#define COMMENT(comm) 
COMMENT(/* The last definition is necessary to include this file in C and in Fortran files. */)
#ifdef EPIK
#include "epik_user.inc"
#else
#define EPIK_FUNC_REG(x) 
#define EPIK_FUNC_START()
#define EPIK_FUNC_END()
#define EPIK_USER_REG(x,y)
#define EPIK_USER_START(x)
#define EPIK_USER_END(x)
#endif

#ifdef VTRACE
#include "VT.inc"
#else
#define VTTRACEON()
#define VTTRACEOFF()
#endif


#ifdef OUTPUT
#define DEBUG(level,str) if (debuglevel>=level) write(*,"(A)") str
#define DEBUGASTART(level) if (debuglevel>=level) then
#define DEBUGAEND end if
#else
#define DEBUG(level,str)
#define DEBUGASTART(level)
#define DEBUGAEND
#endif

#ifdef DOUBLE_PREC
COMMENT(/* FFTW routines */)
#define sfftw_plan_many_dft dfftw_plan_many_dft
#define sfftw_plan_dft_1d dfftw_plan_dft_1d
#define sfftw_plan_many_dft_r2c dfftw_plan_many_dft_r2c 
#define sfftw_plan_many_dft_c2r dfftw_plan_many_dft_c2r 
#define sfftw_plan_dft_c2r_2d dfftw_plan_dft_c2r_2d 
#define sfftw_plan_dft_r2c_2d dfftw_plan_dft_r2c_2d 
#define sfftw_execute_dft dfftw_execute_dft
#define sfftw_execute dfftw_execute
#define sfftw_execute_dft_c2r dfftw_execute_dft_c2r
#define sfftw_execute_dft_r2c dfftw_execute_dft_r2c
#define sfftw_destroy_plan dfftw_destroy_plan
#define sfftw_init_threads dfftw_init_threads
#define sfftw_plan_with_nthreads dfftw_plan_with_nthreads
#define sfftw_cleanup_threads dfftw_cleanup_threads
COMMENT(/* BLAS routines */)
#define saxpy	daxpy
#define caxpy	zaxpy
#define ccopy	zcopy
#define szaxpy	dzaxpy
#define cdotu   zdotu
#define sdot    ddot
#define cgemv   zgemv
#define cvem    zvem
COMMENT(/* Fourier DFT (MKL) routines */)
#define scft	dcft
#define srcft	drcft
#define scrft	dcrft
COMMENT(/* LAPACK routines */)
#define cgetrf  zgetrf
#define cgetri  zgetri
#define sgetrf  dgetrf
#define sgetri  dgetri
#define cgetmo zgetmo
#define sgetmo dgetmo
#define cgbtrf zgbtrf
#define sgbtrf dgbtrf
#define cgbtrs zgbtrs
#define sgbtrs dgbtrs
COMMENT(/* SCALAPACK routines */)
#define PCLAHQR PZLAHQR
#define PCGEHRD PZGEHRD
#define PCGEMR2D PZGEMR2D
#define PCUNMHR PZUNMHR
#define PCSCAL PZSCAL
#define PCTRSV PZTRSV
#define PCGEMM PZGEMM
#define PSGEMM PDGEMM
#define pcgemr2d pzgemr2d
#define pcgetrf pzgetrf
#define psgetrf pdgetrf
#define pcgetri pzgetri
#define psgetri pdgetri
#define pcgetrs pzgetrs
#define pcgemm pzgemm
#define pcgemv pzgemv
#define psgemm pdgemm
#define psgemv pdgemv
#define pcgbtrf pzgbtrf
#define pcgbtrs pzgbtrs
#define psgbtrf pdgbtrf
#define psgbtrs pdgbtrs
COMMENT(/* MPI types */)
#define MPI_REAL_TYPE	MPI_DOUBLE_PRECISION
#define MPI_COMPLEX_TYPE	MPI_DOUBLE_COMPLEX
#define DFTI_SINGLE	DFTI_DOUBLE
COMMENT(/* C interoperable types */)
#define C_REAL_TYPE C_DOUBLE
#else
#define MPI_REAL_TYPE	MPI_REAL
#define MPI_COMPLEX_TYPE	MPI_COMPLEX
#define C_REAL_TYPE C_FLOAT
#endif

#if defined(WITHPERF)
#define PERFINIT call perfinit
#define PERFON(str) call perfon(str)
#define PERFOFF call perfoff
#define C_PERFON(str,strlen) perfon(str,strlen)
#define C_PERFOFF() perfoff()
#define PERFOUT(str) call perfout(str)
#define PERF_GET(str,x,y) call perf_get(str,x,y)
#define PERF_RESET(str) call perf_reset(str)
#define PERF_CONTEXT_START(str) call perf_context_start(str)
#define PERF_CONTEXT_END call perf_context_end
#if (WITHPERF==2)
COMMENT( /*all routines are analyzed*/)
#define PERFON_I(str) call perfon(str)
#define PERFOFF_I call perfoff
#else
COMMENT(/* inner routines are not monitored to keep impact on runtime low */)
#define PERFON_I(str)
#define PERFOFF_I
#endif
#else
#define PERFINIT
#define PERFON(str)
#define PERFOFF
#define C_PERFON(str,strlen)
#define C_PERFOFF()
#define PERFON_I(str)
#define PERFOFF_I
#define PERFOUT(str)
#define PERF_GET(str,x,y)
#define PERF_RESET(str)
#define PERF_CONTEXT_START(str)
#define PERF_CONTEXT_END
#endif

#if defined(WITH_LIKWID)
#define LIKWID_INIT call likwid_markerInit()
#define LIKWID_ON(x) call likwid_markerStart(x)
#define LIKWID_OFF(x) call likwid_markerStop(x)
#define LIKWID_CLOSE call likwid_markerClose()
#else
#define LIKWID_HEADER
#define LIKWID_INIT
#define LIKWID_ON(x)
#define LIKWID_OFF(x)
#define LIKWID_CLOSE
#endif

COMMENT("!Preliminary implementation of the new boundary conditions:")
COMMENT("!Replace by define by undef to try new versions. The new parbc")
COMMENT(" seems to work, the new vpbc creates instabilities in some test")
COMMENT(" cases with two species.")
#define oldvpbc
#define oldparbc
#define MAX_SEND_REQUESTS 64
#define MAX_TYPENAME_LENGTH 40

#define WHICHMIC 0
#undef EN_BLOC
