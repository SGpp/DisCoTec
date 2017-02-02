#include "switches.h"
!> Routines for the gyroaveraging operation in the global code
!!
!! Gyro-averaging in real space. Both coordinates, x and y 
!!  are in real space. The only functions which can be called from 
!!  outside this module are the initialization function and the 
!!  gyro-average functions.
!!  This module uses LocalPolynomBase as base for the
!!  x direction.
!!
!\TODO (Re-) Implement a more efficient way for gyro-averages
!!  with periodic b.c. in the radial direction
!!  This functionality was present in the previous full matrix 
!!  implementation (see rev < 2205) but has not been added adapted
!!  to the sparse matrix implementation yet (it requires some 
!!  restructuring and is not top priority since it does only affect
!!  LILO simulations)
#include "redef.h"
#include "intrinsic_sizes.h"
MODULE gyro_average_dd_mod
  USE localpolynombase2d_mod
  USE mpi
  USE communications, ONLY: MY_MPI_COMM_WORLD,COMM_X, communicators,&
       &mpi_comm_x, mpi_comm_spec,my_barrier, calculate_test_sum,&
       &mpi_comm_v
  USE discretization
  USE coordinates, ONLY: lp1,lx,ly,lw,zval,mu,x0,kj,deli
  USE par_other, ONLY: imag,pi, print_ini_msg, p_has_0_mode
  USE par_in, ONLY: rad_bc_type, spec
  USE geometry, ONLY: geom, rhostar, magn_geometry, minor_r,major_r,C_j, q_prof,dqdx_prof
  use BoundaryDescriptionModule
  use BandedMatrixModule
  use MatrixModule
  use VectorModule
  use DerivativeMatrixModule
  use Grid1DModule
  USE file_io, only: get_unit_nr
  use fourier

  IMPLICIT NONE

  !Include "fftw3.f"

  ! variables also used in other modules, basex in the fieldsolver_df module
  TYPE(LocalPolynomBase2d), SAVE :: basex !xlf95 called for SAVE attribute

  ! private variables II: modular variables, set and calculated in this module
  ! set here the derivative_order used for gyro_averaging and field_solver
  REAL, DIMENSION(:,:,:),ALLOCATABLE, PRIVATE :: sqrtgdet
  INTEGER, PRIVATE :: mybw
  INTEGER,PRIVATE :: derivative_order = 4
  INTEGER,PRIVATE :: polynom_degree = 5, nDeriv
  LOGICAL,PRIVATE :: exist_gyromatrix=.FALSE., first=.true.
  !gyro_averaging : using xi,eta coordinate
  LOGICAL:: gyroav_in_xi_eta=.FALSE. ! only false should be use for the moment

  !LOGICAL, PRIVATE :: periodic_boundary=.TRUE.
  REAL,DIMENSION(:),ALLOCATABLE,private :: my_xval, my_yval

  ! gyromatrix contains the matrix for simple gyro-averaging (J0)
  TYPE(BandedMatrixReal), DIMENSION(:,:,:),ALLOCATABLE,TARGET,PRIVATE :: sparse_gyromatrix
#ifdef GDAGGER_G  
  TYPE(BandedMatrixReal), DIMENSION(:,:,:),ALLOCATABLE,TARGET,PRIVATE :: sparse_gyromatrix_dagger
#endif
  TYPE(VectorReal), SAVE :: vec_v,res_v

  INTEGER, DIMENSION(:,:,:),ALLOCATABLE,private :: bandwidth !maximum over x
  TYPE(BoundaryDescription), DIMENSION(:,:,:,:), ALLOCATABLE, PRIVATE :: bandwidth_boundary
  TYPE(Grid1D),save :: xgrid,ygrid

  ! private functions
  !PRIVATE :: new_calculate_J0matrix, old_calculate_J0matrix
  PRIVATE :: calculate_J0matrix
  PRIVATE :: gyro_average_2D, gyro_average_3D
  PRIVATE :: gyro_average_2D_wb,gyro_average_3D_wb
  private :: SparseMatMul


  !> This interface function is accessible from outside the module.
  !! It is the redirected to the gyro_average_nD (n=1-3) functions,
  !! which are private in the module. Another direction is to 
  !! gyro_average_fields.
  INTERFACE gyro_average_dd
     MODULE PROCEDURE gyro_average_2D, gyro_average_2D_real, gyro_average_3D
  END INTERFACE gyro_average_dd

  INTERFACE gyro_average_dd_wb
     MODULE PROCEDURE gyro_average_2D_wb, gyro_average_2D_wb_real, gyro_average_3D_wb
  END INTERFACE gyro_average_dd_wb
CONTAINS

  FUNCTION estimate_gyromatrix_maxbw() RESULT(total_maxbw)
    integer :: total_maxbw

    ! Local variables
    REAL :: upper_limit_gxx, lower_limit_Bfield, lower_limit_gxx, upper_limit_gxz, deli, delj
    INTEGER :: ierr,nj0

    nj0 = nky0

    if(y_local) then
       deli = lx/nx0
    elseif(x_local) then
       deli = ly/nky0
    endif

    delj = ly/nj0


    SELECT CASE (magn_geometry)
    CASE('circular') 
       lower_limit_gxx = 1.0
       upper_limit_gxx = 1.0
       upper_limit_gxz = 1.0
       lower_limit_Bfield = 0.5
    case default
       lower_limit_gxx = 0.2
       upper_limit_gxx = 20.0
       upper_limit_gxz = 1.0
       lower_limit_Bfield = 0.5
    END SELECT

    total_maxbw = 2*MAXVAL( CEILING( &
         & SQRT(2.0*spec(:)%mass*spec(:)%temp*lw*upper_limit_gxx&
         & /spec(:)%charge**2/lower_limit_Bfield)/deli  &
         & + 0.5*ni0*(upper_limit_gxz*rhostar*minor_r)**2/lower_limit_gxx ) ) + derivative_order + 1
    ! and then also globally
    CALL mpi_allreduce(MPI_IN_PLACE,total_maxbw,1,MPI_INTEGER,MPI_MAX,mpi_comm_spec,ierr)
  END FUNCTION estimate_gyromatrix_maxbw

  !>Give an estimate of the memory requirements of this module
  !!\todo add a better estimate for the sparse gyromatrix
  !!\todo upper limit for gxx has to be improved for strange geometries, same for lower limit of Bfield
  FUNCTION mem_est_gyro_average_dd(mem_req_in)  RESULT(memory_need)
    REAL:: mem_req_in !, rhomax

    ! Local variables
    REAL :: memory_need, mem_xmatrix, mem_xgrid1d
    real :: mem_banded_xmatrix
    INTEGER :: total_maxbw,nj0

    nj0 = nky0

    total_maxbw = estimate_gyromatrix_maxbw()

    mem_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*li0*ni0
    mem_banded_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*li0*total_maxbw

    mem_xgrid1d = size_of_grid1d(ni0)/1024./1024.

    ! modular variable
    ! basex
    memory_need = mem_req_in + size_of_localpolynombase(ni0,nj0)/1024./1024.
    ! sqrtgdet
    memory_need = memory_need + (pi2gl-pi1gl+1)*(pj2-pj1+1)*lk0*SIZE_OF_REAL_MB
    ! some scalars
    memory_need = memory_need + (5*SIZE_OF_INTEGER_MB + 3*SIZE_OF_LOGICAL_MB)
    ! my_xval
    memory_need = memory_need + SIZE_OF_REAL_MB*ni0

    ! sparse_gyromatrix 
    memory_need = memory_need + mem_banded_xmatrix*lj0*lk0*lm0*ln0
#ifdef GDAGGER_G
    ! sparse_gyromatrix_dagger
    memory_need = memory_need + mem_banded_xmatrix*lj0*lk0*lm0*ln0
#endif
    ! vec_m, res_m
    memory_need = memory_need + 2*static_size_of_vector()/1024./1024.

    ! bandwidth
    memory_need = memory_need + SIZE_OF_INTEGER_MB*lk0*lm0*ln0
    ! bandwidth_boundary
    memory_need = memory_need + static_size_of_BoundaryDescription()/1024./1024.*lk0*lm0*ln0

    ! calculate_J0matrix
    ! xgrid or tempmat, because not both at the same time allocated
    memory_need = memory_need + MAX(mem_xgrid1d,mem_xmatrix)
    ! sqrtgxx, prho_gdy_1/2, prho_gdz_1/2, prho_gdphi
    memory_need = memory_need + SIZE_OF_REAL_MB*(pi2gl-pi1gl)*lk0
    ! rho_gdx, rho_gdy, rho_gdz, rho_gdphi, ntheta=100 assumed
    memory_need = memory_need + SIZE_OF_REAL_MB*100
    ! sum_over_theta
    memory_need = memory_need + SIZE_OF_COMPLEX_MB*(nDeriv+1)
    ! integral
    memory_need = memory_need + mem_xmatrix*nDeriv
    ! derivmat
    memory_need = memory_need + (static_size_of_DerivativeMatrix()/1024./1024.&
         & + mem_banded_xmatrix + mem_xgrid1d)*nDeriv

    ! rho
    memory_need = memory_need + SIZE_OF_COMPLEX_MB*(pi2gl-pi1gl)
  END FUNCTION mem_est_gyro_average_dd


  SUBROUTINE SparseMatMul(vec,j,k,m,n,res, transA_in)
    INTEGER,INTENT(IN) :: j,k,m,n
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0)), INTENT(IN),target :: vec
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0)),INTENT(INOUT),target :: res
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER :: transA
    INTEGER :: nj0

    nj0 = nky0

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    IF (.NOT.isInitialized(vec_v)) THEN
       ! initialize the vec_v and res_v vectors
       CALL initialize(vec_v,ni0*nj0)
       CALL initialize(res_v,ni0*nj0)
    END IF

    CALL attach(vec_v,vec)

    CALL attach(res_v,res)

    if (transA.eq.'N') then
       CALL dot_multiply(sparse_gyromatrix(k,m,n),vec_v,res_v) !,transA)
#ifdef G_DAGGER
    elseif (transA.eq.'C') then
       CALL dot_multiply(sparse_gyromatrix_dagger(k,m,n),vec_v,res_v)
#endif
    else
       stop 'current choice of transA not supported in gyro_average_dd.F90/SparseMatMul'
    endif

  END SUBROUTINE SparseMatMul

#undef OLD_VPAR_FOR_GA

  SUBROUTINE gyro_average_2D(func, barfunc,k,m,n,v_dependence,transA_in)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
    REAL, DIMENSION(:,:), ALLOCATABLE :: func_real,barfunc_real
    REAL, DIMENSION(:), ALLOCATABLE :: func_temp,barfunc_temp
    INTEGER,INTENT(IN) :: k,m,n
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in
    ! Local variables
    INTEGER :: j,nj0,nx,ix,iy,ixy !,ierror, local_number_of_y_points
    logical :: v_dep
    CHARACTER :: transA

    nj0 = nky0
    nx = li2-li1+1

    ALLOCATE(func_real(li1:li2,1:nj0),func_temp(1:nx*nj0),&
         &barfunc_real(li1:li2,1:nj0),barfunc_temp(1:nx*nj0))

    DO ix=li1,li2
       CALL fft_ky_to_y(func(ix,:),func_real(ix,:))
    END DO

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          func_temp(ixy) = func_real(ix,iy)
       END DO
    END DO

    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga2D')
#ifdef OLD_VPAR_FOR_GA
    j=1
    CALL SparseMatMul(func_temp,j,k,m,n,barfunc_temp,transA)

#else
    !	TMB - not worrying about vpar parallelization for now.
    !
    !  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
    !     local_number_of_y_points = lj0/n_procs_v
    !     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
    !        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
    !     END DO
    ! now the vpar processes have calculated a y slice of the gyro-averaged
    ! vector barfunc. We now have to gather these slices and broadcast them
    ! to all other processes. We therefore need an allgather operation.
    !     call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
    !          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
    !  else
    j=1
    CALL SparseMatMul(func_temp,j,k,m,n,barfunc_temp,transA)

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          barfunc_real(ix,iy) = barfunc_temp(ixy)
       END DO
    END DO

    DO ix=li1,li2
       CALL fft_y_to_ky(barfunc_real(ix,:),barfunc(ix,:))
    END DO

    !  end if
#endif
    PERFOFF_I

  END SUBROUTINE gyro_average_2D

  SUBROUTINE gyro_average_2D_real(func_real, barfunc_real,k,m,n,v_dependence,transA_in)
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0),lk1:lk2), INTENT(IN) :: func_real
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0),lk1:lk2), INTENT(OUT) :: barfunc_real
    INTEGER,INTENT(IN) :: k,m,n
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in
    ! Local variables
    INTEGER :: j,nj0,nx !,ierror, local_number_of_y_points
    logical :: v_dep
    CHARACTER :: transA

    nj0 = nky0
    nx = li2-li1+1

    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga2D')
#ifdef OLD_VPAR_FOR_GA
    j=1
    CALL SparseMatMul(func_real,j,k,m,n,barfunc_real,transA)
#else
    !       TMB - not worrying about vpar parallelization for now.
    !
    !  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
    !     local_number_of_y_points = lj0/n_procs_v
    !     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
    !        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
    !     END DO
    ! now the vpar processes have calculated a y slice of the gyro-averaged
    ! vector barfunc. We now have to gather these slices and broadcast them
    ! to all other processes. We therefore need an allgather operation.
    !     call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
    !          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
    !  else
    j=1
    CALL SparseMatMul(func_real,j,k,m,n,barfunc_real,transA)

    !  end if
#endif
    PERFOFF_I

  END SUBROUTINE gyro_average_2D_real


  SUBROUTINE gyro_average_2D_wb(func, barfunc,k,m,n,v_dependence,transA_in)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
    REAL, DIMENSION(:,:), ALLOCATABLE :: func_real, barfunc_real
    REAL, DIMENSION(:), ALLOCATABLE :: func_temp, barfunc_temp
    INTEGER,intent(IN) :: k,m,n
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in

    ! Local variables
    INTEGER :: j,nj0,nx,ix,iy,ixy !,ierror, local_number_of_y_points
    CHARACTER :: transA
    logical :: v_dep

    nj0 = nky0
    nx = li2-li1+1

    ALLOCATE(func_real(li1:li2,1:nj0),func_temp(1:nx*nj0),&
         &barfunc_real(li1:li2,1:nj0),barfunc_temp(1:nx*nj0))

    DO ix=li1,li2
       CALL fft_ky_to_y(func(ix,:),func_real(ix,:))
    END DO

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          func_temp(ixy) = func_real(ix,iy)
       END DO
    END DO


    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga2Dwb')

    ! Distribute the following loop over all vpar processes as they
    ! all calculate redundantly the same, due to the fact that the 
    ! gyromatrix does not depend on the vpar direction.
    ! the old loop was
#ifdef OLD_VPAR_FOR_GA
    j=1
    CALL SparseMatMul(func_temp,j,k,m,n,barfunc_temp,transA)
#else
    !  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
    !     local_number_of_y_points = lj0/n_procs_v
    !     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
    !        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
    !     END DO
    !     ! now the vpar processes have calculated a y slice of the gyro-averaged
    !     ! vector barfunc. We now have to gather these slices and broadcast them
    !     ! to all other processes. We therefore need an allgather operation.
    !     call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
    !          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
    !  else
    j=1
    CALL SparseMatMul(func_temp,j,k,m,n,barfunc_temp,transA)

    !  end if
#endif
    PERFOFF_I

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          barfunc_real(ix,iy) = barfunc_temp(ixy)
       END DO
    END DO

    DO ix=li1,li2
       CALL fft_y_to_ky(barfunc_real(ix,:),barfunc(ix,:))
    END DO

  END SUBROUTINE gyro_average_2D_wb

  SUBROUTINE gyro_average_2D_wb_real(func_real,barfunc_real,k,m,n,v_dependence,transA_in)
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0),lk1:lk2), INTENT(IN) :: func_real
    REAL, DIMENSION(1:((li2-li1+1)*2*nky0),lk1:lk2), INTENT(OUT) :: barfunc_real
    INTEGER,intent(IN) :: k,m,n
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in

    ! Local variables
    INTEGER :: j,nj0,nx !,ierror, local_number_of_y_points
    CHARACTER :: transA
    logical :: v_dep

    nj0 = nky0
    nx = li2-li1+1

    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga2Dwb')
    ! Distribute the following loop over all vpar processes as they
    ! all calculate redundantly the same, due to the fact that the
    ! gyromatrix does not depend on the vpar direction.
    ! the old loop was
#ifdef OLD_VPAR_FOR_GA
    j=1
    CALL SparseMatMul(func_real,j,k,m,n,barfunc_real,transA)
#else
    !  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
    !     local_number_of_y_points = lj0/n_procs_v
    !     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
    !        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
    !     END DO
    !     ! now the vpar processes have calculated a y slice of the gyro-averaged
    !     ! vector barfunc. We now have to gather these slices and broadcast them
    !     ! to all other processes. We therefore need an allgather operation.
    !     call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
    !          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
    !  else
    j=1
    CALL SparseMatMul(func_real,j,k,m,n,barfunc_real,transA)

    !  end if
#endif
    PERFOFF_I

  END SUBROUTINE gyro_average_2D_wb_real

  SUBROUTINE gyro_average_3D(func, barfunc, m, n, v_dependence,transA_in)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2), INTENT(OUT) :: barfunc
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: func_real, barfunc_real
    REAL, DIMENSION(:,:), ALLOCATABLE :: func_temp, barfunc_temp
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in
    !local variables
    CHARACTER :: transA
    logical :: v_dep
    INTEGER :: k,m,n,nj0,nx,ix,iy,ixy,j

    nj0 = nky0
    nx = li2-li1+1

    ALLOCATE(func_real(li1:li2,1:nj0,lk1:lk2),func_temp(1:nx*nj0,lk1:lk2),&
         &barfunc_real(li1:li2,1:nj0,lk1:lk2),barfunc_temp(1:nx*nj0,lk1:lk2))

    DO k=lk1,lk2
       DO ix=li1,li2
          CALL fft_ky_to_y(func(ix,:,k),func_real(ix,:,k))
       END DO
    END DO

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          func_temp(ixy,lk1:lk2) = func_real(ix,iy,lk1:lk2)
       END DO
    END DO

    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga3D')

    j=1
    DO k=lk1,lk2
       CALL gyro_average_dd(func_temp,barfunc_temp,k,m,n,v_dep,transA)
    END DO


    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          barfunc_real(ix,iy,lk1:lk2) = barfunc_temp(ixy,lk1:lk2)
       END DO
    END DO

    DO k=lk1,lk2
       DO ix=li1,li2
          CALL fft_y_to_ky(barfunc_real(ix,:,k),barfunc(ix,:,k))
       END DO
    END DO

    PERFOFF_I
  END SUBROUTINE gyro_average_3D

  SUBROUTINE gyro_average_3D_wb(func, barfunc, m, n, v_dependence, transA_in)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(OUT) :: barfunc
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: func_real,barfunc_real
    REAL, DIMENSION(:,:), ALLOCATABLE :: func_temp,barfunc_temp
    logical, intent(IN), optional :: v_dependence
    CHARACTER, intent(IN), optional :: transA_in
    !local variables
    INTEGER :: k,m,n,nj0,nx,ix,iy,ixy
    CHARACTER :: transA
    logical :: v_dep

    nj0 = nky0
    nx = li2-li1+1

    ALLOCATE(func_real(li1:li2,1:nj0,lk1:lk2),func_temp(1:nx*nj0,lk1:lk2),&
         &barfunc_real(li1:li2,1:nj0,lk1:lk2),barfunc_temp(1:nx*nj0,lk1:lk2))

    DO k=lk1,lk2
       DO ix=li1,li2
          CALL fft_ky_to_y(func(ix,:,k),func_real(ix,:,k))
       END DO
    END DO

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          func_temp(ixy,lk1:lk2) = func_real(ix,iy,lk1:lk2)
       END DO
    END DO

    if (.not.present(v_dependence)) then
       v_dep = .true.
    else
       v_dep = v_dependence
    end if

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif

    PERFON_I('ga3Dwb')

    DO k=lk1,lk2
       CALL gyro_average_dd_wb(func_temp,barfunc_temp,k,m,n,v_dep,transA)
    END DO

    DO ix=li1,li2
       DO iy=1,nj0
          ixy = (ix-1)*nj0 + iy
          barfunc_real(ix,iy,lk1:lk2) = barfunc_temp(ixy,lk1:lk2)
       END DO
    END DO

    DO k=lk1,lk2
       DO ix=li1,li2
          CALL fft_y_to_ky(barfunc_real(ix,:,k),barfunc(ix,:,k))
       END DO
    END DO

    PERFOFF_I

  END SUBROUTINE gyro_average_3D_wb


  !> The gyromatrix is private in this module. If one wants to get
  !! it from outside, this function here has to be used. It is mainly
  !! called in the field solver and flr corr modules.
  FUNCTION get_gyromatrix() RESULT(ptr_gm)
    TYPE(BandedMatrixReal), DIMENSION(:,:,:), POINTER :: ptr_gm

    IF (exist_gyromatrix) THEN
       ptr_gm => sparse_gyromatrix
    ELSE
       ptr_gm => NULL()
    END IF
  END FUNCTION get_gyromatrix

  FUNCTION get_gyromatrix_dagger() RESULT(ptr_gm)
    TYPE(BandedMatrixReal), DIMENSION(:,:,:), POINTER :: ptr_gm
    ptr_gm => NULL()
#ifdef GDAGGER_G
    IF (exist_gyromatrix) &
         ptr_gm => sparse_gyromatrix_dagger
#else
    stop 'gyromatrix_dagger does not exist'
#endif
  END FUNCTION get_gyromatrix_dagger

  FUNCTION get_nDeriv_base()
    integer :: get_nDeriv_base

    get_nDeriv_base = get_nDeriv(basex)
  END FUNCTION get_nDeriv_base

  SUBROUTINE allocate_gyromatrix
    INTEGER :: n,m,k,nj0

    LOGICAL :: transposed = .TRUE.

    ALLOCATE(sparse_gyromatrix(lk1:lk2,lm1:lm2,ln1:ln2))
#ifdef GDAGGER_G
    ALLOCATE(sparse_gyromatrix_dagger(lk1:lk2,lm1:lm2,ln1:ln2))
#endif

    nj0 = nky0

    DO n=ln1,ln2
       DO m=lm1,lm2
          DO k=lk1,lk2

             CALL initialize(sparse_gyromatrix(k,m,n),ni0*nj0,ni0*nj0,transposed)
             CALL ALLOCATE(sparse_gyromatrix(k,m,n))
#ifdef GDAGGER_G
             CALL initialize(sparse_gyromatrix_dagger(k,m,n),ni0*nj0,ni0*nj0,transposed) !!,(.not.transposed)) !revert?
             CALL ALLOCATE(sparse_gyromatrix_dagger(k,m,n))
#endif

          END DO
       END DO
    END DO
  END SUBROUTINE allocate_gyromatrix


  SUBROUTINE initialize_gyro_average_dd
    INTEGER :: k,m,n,ix,iy,nj0

    nj0 = nky0

    !PERFON('ga_init')

    IF ((mype.EQ.0).and.(print_ini_msg)) &
         &PRINT*,"Initializing gyro_average_dd using BANDED matrices."
    ALLOCATE(sqrtgdet(pi1gl:pi2gl,pj1:pj2,lk1:lk2))

    sqrtgdet = sqrt(geom%gii*geom%gjj-geom%gij**2)

    ! from the main program we get the indices of the first and last
    ! point, li1 and li2, resp. 
    ! In the periodic case, li1 is the left boundary point, while
    ! li2 is NOT the right boundary point, but the butlast point. As we need
    ! the right point to lie on the right boundary, we have to add one point
    ! internally in this module.
    ! It is always assumed that the first and last point lie on the 
    ! boundaries. 
    ! However, in the periodic case, the first and last point are
    ! therefore identical and we have one independent coefficient less
    ! than we have points.

    ! initialize the underlying grid which is called xgrid
    CALL initialize(xgrid,ni0)
    CALL initialize(ygrid,nj0)
    ! set boundary of xgrid to 0 and lp1 and non-periodic boundary

    IF (rad_bc_type.EQ.0) THEN
       CALL set_boundaries(xgrid,x0/rhostar-lp1/2,x0/rhostar+lp1/2,.TRUE.)
    ELSEIF ((rad_bc_type.EQ.1).or.(rad_bc_type.EQ.2)) THEN
       CALL set_boundaries(xgrid,x0/rhostar-lp1/2,x0/rhostar+lp1/2,.FALSE.)
       !PRINT*,"lp1 = ", lp1
    END IF

    CALL set_boundaries(ygrid,-ly/2,ly/2,.TRUE.)

    ! set up the basis
    call initialize_BandedMatrix_module
    !CALL initialize_polybase2d(basex,polynom_degree,derivative_order,ni0,rad_bc_type)
    CALL initialize_polybase2d(basex,polynom_degree,derivative_order,xgrid,ygrid)
    nDeriv = get_nDeriv(basex)

    ALLOCATE(my_xval(1:ni0)) !my_xval is NOT identical to xval for LILO cases
    ALLOCATE(my_yval(1:nj0))
    CALL set_startindex(xgrid,1)
    CALL set_startindex(ygrid,1)
    DO ix=1,ni0
       my_xval(ix) = get_node(xgrid,ix)
       !    WRITE(*,*) my_xval(ix)
    END DO
    DO iy=1,nky0
       my_yval(iy) = get_node(ygrid,iy)
       !    WRITE(*,*) my_yval(iy)
    END DO

    ! the matrix gyromatrix can be smaller, the second dimension must only
    ! contain bandwidth elements. The first two dimensions span the matrix,
    ! then we need a matrix for each derivative, for all kj, all z,all mu and
    ! all species

    ! using the MatrixModule
    call allocate_gyromatrix
    ALLOCATE(bandwidth(lk1:lk2,lm1:lm2,ln1:ln2))
    ! Temporay variable for testing phi boundary

    !initialize bandwidth
    bandwidth=0

    call calculate_J0matrix

    ALLOCATE(bandwidth_boundary(lk1:lk2,lm1:lm2,ln1:ln2,2))

    DO n=ln1,ln2
       DO m= lm1,lm2
          DO k=lk1,lk2
             CALL initialize_type(bandwidth_boundary(k,m,n,1),1,li0+2*bandwidth(k,m,n),&
                  & bandwidth(k,m,n),bandwidth(k,m,n),lj0)
             CALL initialize_type(bandwidth_boundary(k,m,n,2),1,li0+2*bandwidth(k,m,n),&
                  & bandwidth(k,m,n),bandwidth(k,m,n),1)
          END DO
       END DO
    END DO

    IF (n_procs_x.GT.1) THEN
       ! bandwidth is now set and we can define the mpi datatypes for bandwidth-boundary exchange
       DO n=ln1,ln2
          DO m=lm1,lm2
             DO k=lk1,lk2
                CALL set_mpi_type(bandwidth_boundary(k,m,n,1))
                CALL set_mpi_type(bandwidth_boundary(k,m,n,2))
             END DO
          END DO
       END DO
    END IF
    !PERFOFF
  END SUBROUTINE initialize_gyro_average_dd

#ifdef GA_WITHTEST
  SUBROUTINE test_gyro_average
    COMPLEX, DIMENSION(:,:,:),ALLOCATABLE :: field
    COMPLEX, DIMENSION(:,:,:,:,:), ALLOCATABLE :: fieldbar
    REAL :: kx(0:ni0-1), kxmin, phase(0:ni0-1)
    INTEGER :: thisunit,nz_ind=0, ikx, k,j,i,m,mykx,icount,n
    LOGICAL :: op

    !ALLOCATE(field(li1:li2,lj1:lj2,lk1:lk2),field_deriv(li1:li2,lj1:lj2,lk1:lk2))
    ALLOCATE(field(li1:li2,lj1:lj2,lk1:lk2))
    ALLOCATE(fieldbar(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))

    kxmin = 2.0*pi/lp1
    DO ikx=0,ni0-1
       IF (ikx.LE.(ni0)/2) THEN
          kx(ikx) = ikx*kxmin
       ELSE
          kx(ikx) = (ikx-ni0-2)*kxmin
       END IF
    END DO

    !PRINT*,"kx = ",kx
    call random_number(phase)
    phase = phase*2*pi
    mykx = 1
    phase(0)=0.0

    ! set the field to a superposition of modes
    DO k=lk1,lk2
       DO j=lj1,lj2
          DO i=1,ni0
             field(li1+i-1,j,k) = SUM(CMPLX(COS(kx*my_xval(i)+phase),SIN(kx*,my_xval(i)+phase)))
             !field_deriv(i,j,k) = SUM(CMPLX(-kx*SIN(kx*my_xval(i)+phase),kx*COS(kx*my_xval(i)+phase)))
          END DO
       END DO
    END DO

    thisunit=30
    DO 
       INQUIRE(thisunit,opened=op)
       IF (op) THEN
          thisunit = thisunit+1
       ELSE 
          EXIT
       END IF
    END DO

    OPEN(thisunit,file='./ga_test_before.dat')
    WRITE(thisunit,"(4I4)") ni0,lj2-lj1+1,lk2-lk1+1,lm2-lm1+1
    DO i=1,ni0
       WRITE(thisunit,"(ES20.10)") my_xval(i)
       DO j=lj1,lj2
          WRITE(thisunit,"(ES20.10)") kj(j)
       END DO
    END DO
    DO k=lk1,lk2
       DO j=lj1,lj2
          DO i=li1,li2
             WRITE(thisunit,"(2ES20.10)") field(i,j,k)
          END DO
       END DO
    END DO
    CLOSE(thisunit)

    PERFON('ga_loop')
    DO icount=1,1
       !CALL gyro_average(field,fieldbar)
       DO n=ln1,ln2
          DO m=lm1,lm2
             CALL gyro_average_dd(field,fieldbar(:,:,:,m,n),m,n)
          END DO
       END DO
    END DO
    PERFOFF

    DO m=lm1,lm2
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li1,li2
                WRITE(thisunit,"(2ES20.10E3)") fieldbar(i,j,k,m,ln1)
             END DO
          END DO
       END DO
    END DO

    CLOSE(thisunit)
  END SUBROUTINE test_gyro_average
#endif

#ifdef GA_OUTPUT
  SUBROUTINE print_gyro_matrix_dd
    REAL, DIMENSION(:,:),ALLOCATABLE :: gyro_matrix
    INTEGER :: thisunit,nz_ind=0, ikx, k,j,i,m,mykx,icount,n
    INTEGER :: ix1,ix2,ik,im,in,ij 
    LOGICAL :: op

    ALLOCATE(gyro_matrix(1:ni0*nky0,1:ni0*nky0))

    ik=0
    im=lm2
    in=0
    ij=0

    DO ix1=1,ni0*nky0
       DO ix2=1,ni0*nky0
          gyro_matrix(ix1,ix2) = mat_get_value(sparse_gyromatrix(ik,im,in),ix1,ix2)
       END DO
    END DO

    thisunit=30
    DO
       INQUIRE(thisunit,opened=op)
       IF (op) THEN
          thisunit = thisunit+1
       ELSE
          EXIT
       END IF
    END DO

    OPEN(thisunit,file='./gyro_matrix.dat')

    DO ix1=1,ni0*nky0
       DO ix2=1,ni0*nky0
          WRITE(thisunit,"(ES20.10)") gyro_matrix(ix1,ix2)
       END DO
    END DO

    CLOSE(thisunit)
  END SUBROUTINE print_gyro_matrix_dd
#endif

  !> DESCRIPTION
  !! In this subroutine, the gyromatrix is calculated, according
  !! to the method described in the documentation.
  SUBROUTINE calculate_J0matrix
    ! setup of the gyroaveraging matrix, which is described in my
    ! script for all x positions
    LOGICAL :: transposed=.true.
    REAL :: xprime,yprime,delj
    INTEGER :: startbase_x,startbase_y,endbase_x,endbase_y
    INTEGER :: ispec,j,k,m,ix,itheta,ibase_x,ibase_y, loc_bandwidth_x,iy,&
         &loc_bandwidth_y,ixy,ibase_xy
    INTEGER :: ntheta,ideriv,pii, ndtheta, nj0
    REAL :: dtheta, theta, rhocostheta, rhosintheta !,invN 
    REAL :: rho_Lref !,local_sum
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: sqrtgxx, prho_gdy_1,&
         &prho_gdy_2,prho_gdz_1,prho_gdz_2,prho_gdphi
    REAL,DIMENSION(:), ALLOCATABLE ::  rho_gdx,rho_gdy,rho_gdz,rho_gdphi
    REAL,DIMENSION(:,:), ALLOCATABLE :: C_j_temp
    REAL, DIMENSION(:),ALLOCATABLE :: sum_over_theta
    TYPE(BandedMatrixReal), DIMENSION(:), ALLOCATABLE :: integral
    type(BandedMatrixReal) :: tempmat
    ! derivmat is only a real matrix, but at the moment, we have only a complex
    ! matrix type, so we use it as a complex matrix
    TYPE(DerivativeMatrixReal), DIMENSION(:), allocatable :: derivmat, derivmat_vN
    !  TYPE(DerivativeMatrix), DIMENSION(:), allocatable :: derivmat, derivmat_vN

    REAL,dimension(0:nDeriv) :: res
    REAL, DIMENSION(:), ALLOCATABLE :: rho

    nj0 = nky0

    !PERFON('ga_mat')

    ! allocate the derivative matrices for all derivatives and also the 
    ! integral matrix contains the theta integral
    ALLOCATE(derivmat(1:nDeriv),integral(1:nDeriv))
    IF (rad_bc_type.eq.2) ALLOCATE(derivmat_vN(1:nDeriv))
    ALLOCATE(sum_over_theta(0:nDeriv))
    ALLOCATE(rho(pi1gl:pi2gl))

    DO iDeriv=1,nDeriv
       CALL initialize(derivmat(iDeriv),xgrid,derivative_order)
       IF (rad_bc_type.eq.2) THEN
          CALL initialize(derivmat_vN(iDeriv),xgrid,derivative_order)
          CALL calculate_real(derivmat(iDeriv),iDeriv,1) 
          CALL calculate_real(derivmat_vN(iDeriv),iDeriv,rad_bc_type)
       ELSE
          CALL calculate_real(derivmat(iDeriv),iDeriv,rad_bc_type)
       ENDIF
       CALL initialize(integral(iDeriv),ni0,ni0,transposed)
       CALL ALLOCATE(integral(iDeriv))
    END DO

    ! initialize a temporary matrix
    CALL initialize(tempmat,ni0, ni0, transposed)
    CALL allocate(tempmat)
    pii = pi1

    !TMB - for now, don't allocate y from 1:nj0.  Needs to be fixed though. 
    !some temporary variables and prefactors
    ALLOCATE(sqrtgxx(pi1gl:pi2gl,1,lk1:lk2),prho_gdy_1(pi1gl:pi2gl,1,lk1:lk2),&
         &prho_gdy_2(pi1gl:pi2gl,1,lk1:lk2),prho_gdz_1(pi1gl:pi2gl,1,lk1:lk2),&
         &prho_gdz_2(pi1gl:pi2gl,1,lk1:lk2),prho_gdphi(pi1gl:pi2gl,1,lk1:lk2),&
         &C_j_temp(pi1gl:pi2gl,1))

    !TMB - same here
    !  DO iy=1
    C_j_temp(:,1) = C_j(:)
    !  ENDDO

    rho_Lref = rhostar*minor_r

    Do k=lk1,lk2
       sqrtgxx(:,:,k) = SQRT(geom%gii(:,:,k))
       prho_gdy_1(:,:,k) = geom%gij(:,:,k)/sqrtgxx(:,:,k)
       prho_gdy_2(:,:,k) = sqrtgdet(:,:,k)/sqrtgxx(:,:,k)
       ! gxz and gyz from the geometry module are normalized to L_ref. 
       ! They need to be normalized here to rho_ref, since the larmor 
       ! radius is normalized to rho_ref
       prho_gdz_1(:,:,k) = geom%giz(:,:,k)*rho_Lref/sqrtgxx(:,:,k)
       prho_gdz_2(:,:,k) = (geom%gii(:,:,k)*geom%gjz(:,:,k)-geom%gij(:,:,k)*geom%giz(:,:,k)) &
            &*rho_Lref/(sqrtgdet(:,:,k)*sqrtgxx(:,:,k))
       prho_gdphi(:,:,k) = - sqrtgxx(:,:,k)/sqrtgdet(:,:,k)*C_j_temp(:,:)/&
            &major_R**2*rho_Lref
    Enddo

    !local_sum = 0.0D0
    DO ispec=ln1,ln2
       DO m=lm1,lm2
          DO k=lk1,lk2
             rho(pi1gl:pi2gl) = sqrt(2.0*spec(ispec)%mass*spec(ispec)%temp*mu(m)&
                  /spec(ispec)%charge**2/geom%Bfield(pi1gl:pi2gl,0,k))
             DO j=lj1,lj2
                CALL set_zero(sparse_gyromatrix(k,m,ispec))
#ifdef GDAGGER_G
                CALL set_zero(sparse_gyromatrix_dagger(k,m,ispec))
#endif
                DO iDeriv=1,nDeriv
                   CALL set_zero(integral(iDeriv))
                END DO
                ! Attention, the indices for the integral dimension are running from 1 to ni0
                ! not from li1 to li2 !!!
                DO ix=1,ni0
                   DO iy=1,nj0
                      pii=ix-1+pi1gl !li1
                      ixy = (ix-1)*ni0 + iy
                      !ntheta = CEILING(sqrtgxx(pii,k)*rho(pii)/deli) * 2 + 40

                      ntheta = MAX(40,CEILING(sqrtgxx(pii,1,k)*rho(pii)/deli * 40))
                      dtheta = 2*pi/ntheta
                      ALLOCATE(rho_gdx(0:ntheta-1),rho_gdy(0:ntheta-1),rho_gdz(0:ntheta-1)&
                           &,rho_gdphi(0:ntheta-1))
                      DO itheta=0,ntheta-1
                         theta = itheta*dtheta
                         rhocostheta = rho(pii)*COS(theta)              
                         rhosintheta = rho(pii)*SIN(theta)              

                         rho_gdx(itheta) = sqrtgxx(pii,1,k) * rhocostheta         ! rho.grad(x)
                         rho_gdy(itheta) = prho_gdy_1(pii,1,k) * rhocostheta + &
                              & prho_gdy_2(pii,1,k) * rhosintheta                 ! rho.grad(y)
                         rho_gdz(itheta) = prho_gdz_1(pii,1,k) * rhocostheta + &
                              & prho_gdz_2(pii,1,k) * rhosintheta                 ! rho.grad(z)
                         rho_gdphi(itheta) = prho_gdphi(pii,1,k) * rhosintheta
                      END DO
                      !invN = 1.0D0/ntheta
                      !    delj = ly/nj0
                      delj = ly/nky0

                      loc_bandwidth_x=CEILING(rho_gdx(0)/deli)
                      loc_bandwidth_y=CEILING(rho_gdx(0)/delj)

                      !WRITE(*,"(6I3,I4)") mype,ispec,m,k,j,ix,loc_bandwidth

                      startbase_x = ix-loc_bandwidth_x
                      endbase_x   = ix+loc_bandwidth_x

                      startbase_y = iy-loc_bandwidth_y
                      endbase_y = iy+loc_bandwidth_y

                      IF (startbase_x.LT.1) THEN
                         startbase_x = 1
                         IF (rad_bc_type.eq.0) endbase_x=ni0 !periodic boundaries
                      ENDIF
                      IF (endbase_x.GT.ni0) THEN
                         endbase_x=ni0
                         IF (rad_bc_type.eq.0) startbase_x=1 !periodic boundaries
                      ENDIF

                      IF (startbase_y.LT.1) THEN
                         startbase_y = 1
                         IF (rad_bc_type.eq.0) endbase_y=ni0 !periodic boundaries
                      ENDIF
                      IF ( endbase_y.GT.nj0) THEN
                         endbase_y=nj0
                         IF (rad_bc_type.eq.0) startbase_y=1 !periodic boundaries
                      ENDIF


                      ! Due to the multiplication with the derivative matrix,
                      ! which is also banded with a derivative_order/2 upper
                      ! diagonals and the same number of lower diagonals, we
                      ! have to extend the bandwidth for the multiplication:
                      bandwidth(k,m,ispec)=max(bandwidth(k,m,ispec),loc_bandwidth_x+&
                           &derivative_order/2)

                      DO ibase_x=startbase_x,endbase_x
                         DO ibase_y=startbase_y,endbase_y 
                            ibase_xy = (ibase_x-1)*ni0 + ibase_y
                            sum_over_theta = CMPLX(0.0,0.0,KIND(sum_over_theta))
                            ndtheta = 0
                            DO itheta = 0,ntheta-1

                               ! WRITE(*,*) my_xval(ix),my_yval(iy),zval(k)

                               ! compute particle position
                               CALL compute_part_position(xprime,yprime,my_xval(ix),my_yval(iy),zval(k),&
                                    &rho_gdx(itheta),rho_gdy(itheta),rho_gdz(itheta),&
                                    &rho_gdphi(itheta),q_prof(pii), &
                                    &dqdx_prof(pii),C_j(pii))

                               !"mirror" function values at (inner) boundary for von Neumann bc
                               !on kj=0 mode
                               if ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(j==lj1)) then
                                  if (xprime.lt.my_xval(1)) then
                                     xprime = my_xval(1)+(my_xval(1)-xprime)
                                  endif
                               endif

                               !WRITE(*,*) ix,iy,k

                               !WRITE(*,*) xprime,yprime,ibase_x,ibase_y

                               CALL get_value2d(basex,ibase_x,xprime,ibase_y,yprime,res)

                               DO ideriv=0,nderiv
                                  sum_over_theta(ideriv) = sum_over_theta(ideriv) &
                                       & + res(ideriv) 
                               END DO
                            END DO !itheta



                            sum_over_theta = sum_over_theta/ntheta
                            ! WRITE(*,*) sum_over_theta
                            Call set_value(sparse_gyromatrix(k,m,ispec),ixy,ibase_xy,sum_over_theta(0))
                            DO ideriv=1,nderiv
                               CALL set_value(integral(ideriv),ixy,ibase_xy,sum_over_theta(ideriv))
                            ENDDO

                         END DO !ibase_y 
                      END DO !ibase_x

                      DEALLOCATE(rho_gdx,rho_gdy,rho_gdz,rho_gdphi)
                   END DO !ix
                END DO !iy

                CALL commit_values(sparse_gyromatrix(k,m,ispec))

                DO ideriv=1,nderiv
                   CALL commit_values(integral(ideriv))
                END DO

                ! construct the final gyro_average matrix
                DO ideriv=1,nderiv !2,nderiv
                   !CALL set_zero(tempmat)  -> this has been put in the dot_multiply routine
                   IF ((rad_bc_type.eq.2).and.(p_has_0_mode.and.j==lj1)) THEN
                      CALL dot_multiply(integral(ideriv),derivmat_vN(ideriv)%Data,tempmat)
                   ELSE
                      CALL dot_multiply(integral(ideriv),derivmat(ideriv)%Data,tempmat)
                   ENDIF

                   CALL add_matrix(sparse_gyromatrix(k,m,ispec), tempmat)

                END DO
#ifdef GDAGGER_G
                !              call set_dagger_matrix(sparse_gyromatrix(k,m,ispec),&
                !                   &sparse_gyromatrix_dagger(k,m,ispec),bandwidth(k,m,ispec))
#endif

             END DO !j
          END DO ! k
       END DO !m
    END DO ! ispec

    deallocate(sqrtgxx, prho_gdy_1,prho_gdy_2,prho_gdz_1,prho_gdz_2,prho_gdphi)
    deallocate(sum_over_theta)
    call finalize(tempmat)

    DO iDeriv=1,nDeriv
       call finalize(derivmat(iDeriv))
       IF (rad_bc_type.EQ.2) call finalize(derivmat_vN(iDeriv))
       call finalize(integral(iDeriv))
    END DO
    DEALLOCATE(integral,derivmat,rho)
    IF (rad_bc_type.eq.2) DEALLOCATE(derivmat_vN)

    exist_gyromatrix = .TRUE.

  END SUBROUTINE calculate_J0matrix

#ifdef GDAGGER_G
  SUBROUTINE set_dagger_matrix(mat_in,mat_dagger,bandwidth)
    TYPE(BandedMatrixReal),INTENT(IN) :: mat_in
    TYPE(BandedMatrixReal),INTENT(INOUT) :: mat_dagger

    INTEGER, INTENT(IN) :: bandwidth
    !local variables
    INTEGER :: ix, ibase, startbase, endbase
    REAL, DIMENSION(1:ni0*nky0,1:ni0*nky0) :: temp_fullmat

    call get_global_matrix_locally(mat_in,temp_fullmat)
    !!\todo check whether get_global_matrix_locally works on IBM
    !!compilers (because something is strange if GDAGGER is set
    !!on those machines)

    DO ix=1,ni0
       startbase = ix-bandwidth
       endbase   = ix+bandwidth

       IF (startbase.LT.1) THEN
          startbase = 1
          IF (rad_bc_type.eq.0) endbase=ni0 !periodic boundaries
       ENDIF
       IF (endbase.GT.ni0) THEN
          endbase=ni0
          IF (rad_bc_type.eq.0) startbase=1 !periodic boundaries
       ENDIF
       DO ibase=startbase,endbase
          Call set_value(mat_dagger,ibase,ix,&
               &temp_fullmat(ix,ibase))
       ENDDO
    ENDDO

    CALL commit_values(mat_dagger)

  END SUBROUTINE set_dagger_matrix
#endif

  !>Clean up all arrays etc. employed by this module
  SUBROUTINE finalize_gyro_average_dd
    INTEGER :: k,m,n

    ! finalize the xgrid
    call finalize(xgrid)

    DEALLOCATE(my_xval,my_yval,sqrtgdet)

    CALL finalize(vec_v)
    call finalize(res_v)
    CALL finalize_localpolynombase(basex)

    IF (ALLOCATED(sparse_gyromatrix)) THEN
       DO n=ln1,ln2
          DO m=lm1,lm2
             DO k=lk1,lk2

                CALL finalize(sparse_gyromatrix(k,m,n))
#ifdef GDAGGER_G
                CALL finalize(sparse_gyromatrix_dagger(k,m,n))
#endif

             END DO
          END DO
       END DO
       DEALLOCATE(sparse_gyromatrix)
#ifdef GDAGGER_G
       DEALLOCATE(sparse_gyromatrix_dagger)
#endif
    END IF
    DEALLOCATE(bandwidth)


    !periodic_boundary=.TRUE.
    exist_gyromatrix=.FALSE.
    IF (n_procs_x.GT.1) THEN
       DO n=ln1,ln2
          DO m=lm1,lm2
             DO k=lk1,lk2
                CALL finalize_type(bandwidth_boundary(k,m,n,1))
                CALL finalize_type(bandwidth_boundary(k,m,n,2))
             END DO
          END DO
       END DO
    END IF
    DEALLOCATE(bandwidth_boundary)
    call finalize(xgrid)

  END SUBROUTINE finalize_gyro_average_dd


  SUBROUTINE compute_part_position(ao_xprime,ao_yprime_y, ai_x, ai_y, ai_z, ai_rhogdx, ai_rhogdy, ai_rhogdz,&
       &ai_rhogdphi,ai_q,ai_dqdx,ai_Cy)
    REAL, INTENT(out) :: ao_xprime,ao_yprime_y
    REAL, INTENT(in)  :: ai_x, ai_z,ai_rhogdx,ai_rhogdy, ai_rhogdz, &
         &ai_rhogdphi,ai_q,ai_dqdx,ai_Cy,ai_y
    REAL              :: zprime,qxprime

    IF (.NOT.gyroav_in_xi_eta) THEN
       ! using linearisation of the x,y,z metric
       ao_xprime = ai_x + ai_rhogdx
       ao_yprime_y = ai_y + ai_rhogdy
    ELSE

       ! x'=x(X+rho),y'=y(X+rho), z'=z(X+rho) are evaluated using an intermediate 
       ! transformation to pseudo cartesian coord. Xi,eta (see notes on gyroaveraging).
       ! q(x')=q(x)+(x'-x)*dqdx

       ao_xprime = SQRT(ai_x**2+2.0D0*ai_x*ai_rhogdx + &                     
            ai_rhogdx**2 + ai_x**2*ai_rhogdz**2)
       zprime = ATAN2(SIN(ai_z)*(ai_x+ai_rhogdx)+ai_x*COS(ai_z)*ai_rhogdz,&  
            & COS(ai_z)*(ai_x+ai_rhogdx) -ai_x*SIN(ai_z)*ai_rhogdz)
       ! avoid discontinuties around z=pi or -pi :
       IF ((ABS(zprime)>1.5).AND.(ai_z*zprime).LT.0) THEN
          zprime=zprime+2*pi*ai_z/ABS(ai_z)
       END IF
       qxprime = ai_q+(ao_xprime-ai_x)*rhostar*minor_r*ai_dqdx          
       ao_yprime_y = 1.0/(rhostar*minor_r)*ai_Cy*(qxprime*zprime - &                  
            ai_q*ai_z-ai_rhogdphi)
    END IF

  END SUBROUTINE compute_part_position

END MODULE gyro_average_dd_mod
