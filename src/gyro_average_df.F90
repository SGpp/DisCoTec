#include "switches.h"
!> Routines for the gyroaveraging operation in the global code
!!
!! Gyro-averaging in mixed space. The radial coordinate is 
!!  in real space, while the y coordinate is still in Fourier
!!  space. The only functions which can be called from outside
!!  this module are the initialization function and the 
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
MODULE gyro_average_df_mod
  USE localpolynombase_mod
  USE mpi
  USE communications, ONLY: MY_MPI_COMM_WORLD,COMM_X, communicators,&
       &mpi_comm_x, mpi_comm_spec,my_barrier, calculate_test_sum,&
       &mpi_comm_v
  USE discretization
  use coordinates, only: lp1,lx,lx_a,ly,lw,zval,mu,x0,kj,deli
  USE par_other, ONLY: imag,pi, print_ini_msg, p_has_0_mode
  use par_in, only: rad_bc_type, spec, ga_spatial_var
  USE geometry, ONLY: geom, rhostar, magn_geometry, minor_r,major_r,&
       C_y, q_prof, dqdx_prof
  use BoundaryDescriptionModule
  use BandedMatrixModule
  use MatrixModule
  use VectorModule
  use DerivativeMatrixModule
  use Grid1DModule
  USE file_io, only: get_unit_nr
  use lagrange_interpolation
  IMPLICIT NONE

  ! variables also used in other modules, basex in the fieldsolver_df module
  TYPE(LocalPolynomBase), SAVE :: basex !xlf95 called for SAVE attribute

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
  REAL,DIMENSION(:),ALLOCATABLE,private :: my_xval

  ! gyromatrix contains the matrix for simple gyro-averaging (J0)
  TYPE(BandedMatrix), DIMENSION(:,:,:,:),ALLOCATABLE,TARGET,PRIVATE :: sparse_gyromatrix
#ifdef GDAGGER_G  
  TYPE(BandedMatrix), DIMENSION(:,:,:,:),ALLOCATABLE,TARGET,PRIVATE :: sparse_gyromatrix_dagger
#endif
  TYPE(Vector), SAVE :: vec_v, res_v

  INTEGER, DIMENSION(:,:,:),ALLOCATABLE,private :: bandwidth !maximum over x
  TYPE(BoundaryDescription), DIMENSION(:,:,:,:), ALLOCATABLE, PRIVATE :: bandwidth_boundary
  TYPE(Grid1D),save :: xgrid

  ! private functions
  !PRIVATE :: new_calculate_J0matrix, old_calculate_J0matrix
  PRIVATE :: calculate_J0matrix
  PRIVATE :: gyro_average_1D, gyro_average_2D, gyro_average_3D
  PRIVATE :: gyro_average_1D_wb,gyro_average_2D_wb,gyro_average_3D_wb
  private :: SparseMatMul


  !> This interface function is accessible from outside the module.
  !! It is the redirected to the gyro_average_nD (n=1-3) functions,
  !! which are private in the module. Another direction is to 
  !! gyro_average_fields.
  INTERFACE gyro_average_df
     MODULE PROCEDURE gyro_average_1D, gyro_average_2D, gyro_average_3D
  END INTERFACE

  INTERFACE gyro_average_df_wb
     MODULE PROCEDURE gyro_average_1D_wb, gyro_average_2D_wb,gyro_average_3D_wb
  END INTERFACE
  
CONTAINS

  FUNCTION estimate_gyromatrix_maxbw() RESULT(total_maxbw)
    integer :: total_maxbw

    ! Local variables
    REAL :: upper_limit_gxx, lower_limit_Bfield, lower_limit_gxx, upper_limit_gxz, deli
    INTEGER :: ierr

    if(y_local) then
       deli = lx/nx0
    elseif(x_local) then
       deli = ly/nky0
    endif

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

    !if we specified lx_a instead of lx, we don't know the actual 
    !lx yet in the first call to this routine,
    !so we have to guess a hopefully safe value for deli
    if (deli==0.) deli=lx_a*100/nx0

    total_maxbw = 2*MAXVAL( CEILING( &
         & SQRT(2.0*spec(:)%mass*spec(:)%temp*lw*upper_limit_gxx&
         & /spec(:)%charge**2/lower_limit_Bfield)/deli  &
         & + 0.5*ni0*(upper_limit_gxz*rhostar*minor_r)**2/lower_limit_gxx ) ) + derivative_order + 1
    total_maxbw = MIN(total_maxbw, ni0)
    ! and then also globally
    CALL mpi_allreduce(MPI_IN_PLACE,total_maxbw,1,MPI_INTEGER,MPI_MAX,mpi_comm_spec,ierr)
  END FUNCTION estimate_gyromatrix_maxbw

  !>Give an estimate of the memory requirements of this module
  !!\todo add a better estimate for the sparse gyromatrix
  !!\todo upper limit for gxx has to be improved for strange geometries, same for lower limit of Bfield
  FUNCTION mem_est_gyro_average_df(mem_req_in)  RESULT(memory_need)
    REAL:: mem_req_in !, rhomax

    ! Local variables
    REAL :: memory_need, mem_xmatrix, mem_xgrid1d
    real :: mem_banded_xmatrix
    INTEGER :: total_maxbw

    total_maxbw = estimate_gyromatrix_maxbw()

    mem_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*li0*ni0
    mem_banded_xmatrix = static_size_of_matrix()/1024./1024.&
         & + SIZE_OF_COMPLEX_MB*li0*total_maxbw

    mem_xgrid1d = size_of_grid1d(ni0)/1024./1024.

    ! modular variable
    ! basex
    memory_need = mem_req_in + size_of_localpolynombase(ni0)/1024./1024.
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
    ! temp_dagger_matrix
    memory_need = memory_need + mem_banded_xmatrix
    ! for the transpose_storage routine in storebandedmatrix
    ! we must account for two non-distributed local banded matrices
    memory_need = memory_need + 2*SIZE_OF_COMPLEX_MB*nx0*total_maxbw
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
  END FUNCTION mem_est_gyro_average_df


  SUBROUTINE SparseMATMUL(vec,j,k,m,n,res, transA_in)
    INTEGER,INTENT(IN) :: j,k,m,n
    COMPLEX, DIMENSION(li1:li2), INTENT(IN),target :: vec
    COMPLEX, DIMENSION(li1:li2),INTENT(INOUT),target :: res
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER :: transA
    
    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    
    IF (.NOT.isInitialized(vec_v)) THEN
       ! initialize the vec_v and res_v vectors
       CALL initialize(vec_v,ni0)
       CALL initialize(res_v,ni0)
    END IF
    
    CALL attach(vec_v,vec)
    
    CALL attach(res_v,res)
    
    if (transA.eq.'N') then
       CALL dot_multiply(sparse_gyromatrix(j,k,m,n),vec_v,res_v) !,transA)
#ifdef G_DAGGER
    elseif (transA.eq.'C') then
       CALL dot_multiply(sparse_gyromatrix_dagger(j,k,m,n),vec_v,res_v)
#endif
    else
       stop 'current choice of transA not supported in gyro_average_df.F90/SparseMatMul'
    endif
    
  END SUBROUTINE SparseMATMUL

  SUBROUTINE gyro_average_1D(func, barfunc,j,k,m,n, transA_in)
    COMPLEX, DIMENSION(li1:li2), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2), INTENT(INOUT) :: barfunc
    INTEGER,intent(IN) :: j,k,m,n
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER :: transA

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    
    PERFON('ga1D')
    
    CALL SparseMatMul(func(li1:li2),j,k,m,n,barfunc(li1:li2),transA)

    PERFOFF
    
  END SUBROUTINE gyro_average_1D

  SUBROUTINE gyro_average_1D_wb(func, barfunc,j,k,m,n, transA_in)
    COMPLEX, DIMENSION(li1:li2), INTENT(IN) :: func
    COMPLEX, DIMENSION(li1:li2), INTENT(OUT) :: barfunc
    INTEGER,intent(IN) :: j,k,m,n
    CHARACTER, intent(IN), optional :: transA_in
    CHARACTER :: transA

    if (present(transA_in)) then
       transA = transA_in
    else
       transA = 'N'
    endif
    
    PERFON('ga1Dwb')
    
    CALL SparseMatMul(func(li1:li2),j,k,m,n,barfunc(li1:li2),transA)

    PERFOFF
    
  END SUBROUTINE gyro_average_1D_wb
#undef OLD_VPAR_FOR_GA

SUBROUTINE gyro_average_2D(func, barfunc,k,m,n,v_dependence,transA_in)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
  INTEGER,INTENT(IN) :: k,m,n
  logical, intent(IN), optional :: v_dependence
  CHARACTER, intent(IN), optional :: transA_in
  ! Local variables
  INTEGER :: j, ierror, local_number_of_y_points
  logical :: v_dep
  CHARACTER :: transA
  complex, dimension(li1:li2,lj1:lj2) :: tempbarfunc
  
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
  DO j=lj1,lj2
     CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
  END DO
#else
  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
     local_number_of_y_points = lj0/n_procs_v
     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
        !CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,tempbarfunc(li1:li2,j),transA)
     END DO
     ! now the vpar processes have calculated a y slice of the gyro-averaged
     ! vector barfunc. We now have to gather these slices and broadcast them
     ! to all other processes. We therefore need an allgather operation.
     !call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
     !     & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
     call MPI_Allgather(tempbarfunc(li1,lj1+my_pev*local_number_of_y_points),&
          &local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
  else
     DO j=lj1,lj2
        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
     END DO
  end if
#endif
  PERFOFF_I
  
END SUBROUTINE gyro_average_2D

SUBROUTINE gyro_average_2D_wb(func, barfunc,k,m,n,v_dependence,transA_in)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2), INTENT(OUT) :: barfunc
  INTEGER,intent(IN) :: k,m,n
  logical, intent(IN), optional :: v_dependence
  CHARACTER, intent(IN), optional :: transA_in

  ! Local variables
  INTEGER :: j, local_number_of_y_points,ierror
  CHARACTER :: transA
  logical :: v_dep
  complex,dimension(li1:li2,lj1:lj2) :: tempbarfunc

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
  DO j=lj1,lj2
     CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
  END DO
#else
  if ((.not.v_dep).and.(n_procs_v.gt.1).and.(mod(lj0,n_procs_v).eq.0)) then
     local_number_of_y_points = lj0/n_procs_v
     DO j=lj1+my_pev*local_number_of_y_points,lj1+(my_pev+1)*local_number_of_y_points-1
        !CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,tempbarfunc(li1:li2,j),transA)
     END DO
     ! now the vpar processes have calculated a y slice of the gyro-averaged
     ! vector barfunc. We now have to gather these slices and broadcast them
     ! to all other processes. We therefore need an allgather operation.
     !call MPI_Allgather(MPI_IN_PLACE,local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
     !     & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)
     call MPI_Allgather(tempbarfunc(li1,lj1+my_pev*local_number_of_y_points),&
          &local_number_of_y_points*li0,MPI_COMPLEX_TYPE,&
          & barfunc(li1,lj1),local_number_of_y_points*li0,MPI_COMPLEX_TYPE,mpi_comm_v,ierror)

  else
     DO j=lj1,lj2
        CALL SparseMatMul(func(li1:li2,j),j,k,m,n,barfunc(li1:li2,j),transA)
     END DO
  end if
#endif
  PERFOFF_I

END SUBROUTINE gyro_average_2D_wb


SUBROUTINE gyro_average_3D(func, barfunc, m, n, v_dependence,transA_in)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2), INTENT(OUT) :: barfunc
  logical, intent(IN), optional :: v_dependence
  CHARACTER, intent(IN), optional :: transA_in
  !local variables
  CHARACTER :: transA
  logical :: v_dep
  INTEGER :: k,m,n

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

  !PERFON_I('ga3D')

  DO k=lk1,lk2
     CALL gyro_average_df(func(:,:,k),barfunc(:,:,k),k,m,n,v_dep,transA)
  END DO

  !PERFOFF_I
END SUBROUTINE gyro_average_3D

SUBROUTINE gyro_average_3D_wb(func, barfunc, m, n, v_dependence, transA_in)
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(IN) :: func
  COMPLEX, DIMENSION(li1:li2,lj1:lj2,lbz:ubz), INTENT(OUT) :: barfunc
  logical, intent(IN), optional :: v_dependence
  CHARACTER, intent(IN), optional :: transA_in
  !local variables
  INTEGER :: k,m,n
  CHARACTER :: transA
  logical :: v_dep

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

  !PERFON_I('ga3Dwb')

  DO k=lk1,lk2
     CALL gyro_average_df_wb(func(:,:,k),barfunc(:,:,k),k,m,n,v_dep,transA)
  END DO
  !PERFOFF_I

END SUBROUTINE gyro_average_3D_wb


!> The gyromatrix is private in this module. If one wants to get
!! it from outside, this function here has to be used. It is mainly
!! called in the field solver and flr corr modules.
FUNCTION get_gyromatrix() RESULT(ptr_gm)
  TYPE(BandedMatrix), DIMENSION(:,:,:,:), POINTER :: ptr_gm

  IF (exist_gyromatrix) THEN
     ptr_gm => sparse_gyromatrix
  ELSE
     ptr_gm => NULL()
  END IF
END FUNCTION get_gyromatrix

FUNCTION get_gyromatrix_dagger() RESULT(ptr_gm)
  TYPE(BandedMatrix), DIMENSION(:,:,:,:), POINTER :: ptr_gm
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
  INTEGER :: n,m,k,j

  LOGICAL :: transposed = .TRUE.

  ALLOCATE(sparse_gyromatrix(lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
#ifdef GDAGGER_G
  ALLOCATE(sparse_gyromatrix_dagger(lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2))
#endif
  DO n=ln1,ln2
     DO m=lm1,lm2
        DO k=lk1,lk2
           DO j=lj1,lj2
              CALL initialize(sparse_gyromatrix(j,k,m,n),ni0,ni0,transposed)
              CALL ALLOCATE(sparse_gyromatrix(j,k,m,n))
#ifdef GDAGGER_G
              CALL initialize(sparse_gyromatrix_dagger(j,k,m,n),ni0,ni0,transposed)
              CALL ALLOCATE(sparse_gyromatrix_dagger(j,k,m,n))
#else
              ! without GDAGGER_G, the sparse_gyromatrix is the matrix which
              ! is used for multiplication with the vectors, hence it should
              ! be in transposed format
              !CALL initialize(sparse_gyromatrix(j,k,m,n),ni0,ni0,transposed)
              !CALL ALLOCATE(sparse_gyromatrix(j,k,m,n))
#endif
           END DO
        END DO
     END DO
  END DO
END SUBROUTINE allocate_gyromatrix


SUBROUTINE initialize_gyro_average_df
  INTEGER :: k,m,n,ix

  !PERFON('ga_init')

  IF ((mype.EQ.0).and.(print_ini_msg)) &
       & write (*,"(A)") "Initializing gyro_average_df using BANDED matrices."
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
  ! set boundary of xgrid to 0 and lp1 and non-periodic boundary

  IF (rad_bc_type.EQ.0) THEN
     if (y_local) then
        call set_boundaries(xgrid,x0/rhostar-lp1/2,x0/rhostar+lp1/2,.true.)
     else
        call set_boundaries(xgrid,-lp1/2,+lp1/2,.true.)
     endif
  ELSEIF ((rad_bc_type.EQ.1).or.(rad_bc_type.EQ.2)) THEN
     CALL set_boundaries(xgrid,x0/rhostar-lp1/2,x0/rhostar+lp1/2,.FALSE.)
     !PRINT*,"lp1 = ", lp1
  END IF

  ! set up the basis
  call initialize_BandedMatrix_module
  !CALL initialize_polybase(basex,polynom_degree,derivative_order,ni0,rad_bc_type)
  CALL initialize_polybase(basex,polynom_degree,derivative_order,xgrid)
  nDeriv = get_nDeriv(basex)

  ALLOCATE(my_xval(1:ni0)) !my_xval is NOT identical to xval for LILO cases
  CALL set_startindex(xgrid,1)
  DO ix=1,ni0
     my_xval(ix) = get_node(xgrid,ix)
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
END SUBROUTINE initialize_gyro_average_df

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
           CALL gyro_average_df(field,fieldbar(:,:,:,m,n),m,n)
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
SUBROUTINE print_gyro_matrix_df
  COMPLEX, DIMENSION(:,:),ALLOCATABLE :: gyro_matrix
  INTEGER :: thisunit,nz_ind=0, ikx, k,j,i,m,mykx,icount,n
  INTEGER :: ix1,ix2,ik,im,in,ij
  LOGICAL :: op

  ALLOCATE(gyro_matrix(1:ni0,1:ni0))

  ik=0
  im=lm2
  in=0
  ij=0

  DO ix1=1,ni0
     DO ix2=1,ni0   
        gyro_matrix(ix1,ix2) = mat_get_value(sparse_gyromatrix(ij,ik,im,in),ix1,ix2)
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

  DO ix1=1,ni0
     DO ix2=1,ni0
        WRITE(thisunit,"(ES20.10)") gyro_matrix(ix1,ix2)
     END DO
  END DO

  CLOSE(thisunit)
END SUBROUTINE print_gyro_matrix_df
#endif




!> DESCRIPTION
!! In this subroutine, the gyromatrix is calculated, according
!! to the method described in the documentation.
SUBROUTINE calculate_J0matrix
  ! setup of the gyroaveraging matrix, which is described in my
  ! script for all x positions
  LOGICAL :: transposed=.true.
  REAL :: xprime,yprime_y
  integer :: startbase, endbase
  INTEGER :: ispec,j,k,m,ix,itheta,ibase, loc_bandwidth
  INTEGER :: ntheta,ideriv,pii, ndtheta
  REAL :: dtheta, theta, rhocostheta, rhosintheta !,invN 
  REAL :: rho_Lref!,local_sum
  REAL,DIMENSION(:,:),ALLOCATABLE :: sqrtgxx, prho_gdy_1,&
       &prho_gdy_2,prho_gdz_1,prho_gdz_2,prho_gdphi
  REAL,DIMENSION(:),ALLOCATABLE ::  rho_gdx,rho_gdy,rho_gdz,rho_gdphi
  COMPLEX, DIMENSION(:),ALLOCATABLE :: sum_over_theta
  TYPE(BandedMatrix), DIMENSION(:), ALLOCATABLE :: integral
  type(BandedMatrix) :: tempmat
  ! derivmat is only a real matrix, but at the moment, we have only a complex
  ! matrix type, so we use it as a complex matrix
  TYPE(DerivativeMatrix), DIMENSION(:), allocatable :: derivmat, derivmat_vN
!  TYPE(DerivativeMatrix), DIMENSION(:), allocatable :: derivmat, derivmat_vN
#ifdef GDAGGER_G
  type(BandedMatrix) :: temp_dagger_matrix
#endif
  REAL,dimension(0:nDeriv) :: res
  real, dimension(:), allocatable :: rho, rho_x, rho_y

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
        CALL calculate(derivmat(iDeriv),iDeriv,1) 
        CALL calculate(derivmat_vN(iDeriv),iDeriv,rad_bc_type)
     ELSE
        CALL calculate(derivmat(iDeriv),iDeriv,rad_bc_type)
     ENDIF
     CALL initialize(integral(iDeriv),ni0,ni0,transposed)
     CALL ALLOCATE(integral(iDeriv))
  END DO

  ! initialize a temporary matrix
  CALL initialize(tempmat,ni0, ni0, transposed)
  CALL allocate(tempmat)
#ifdef GDAGGER_G
  call initialize(temp_dagger_matrix,ni0,ni0)
  call allocate(temp_dagger_matrix)
#endif
  pii = pi1

  !some temporary variables and prefactors
  ALLOCATE(sqrtgxx(pi1gl:pi2gl,lk1:lk2),prho_gdy_1(pi1gl:pi2gl,lk1:lk2),&
       &prho_gdy_2(pi1gl:pi2gl,lk1:lk2),prho_gdz_1(pi1gl:pi2gl,lk1:lk2),&
       &prho_gdz_2(pi1gl:pi2gl,lk1:lk2),prho_gdphi(pi1gl:pi2gl,lk1:lk2))

  rho_Lref = rhostar*minor_r
  Do k=lk1,lk2
     ! gxz and gyz from the geometry module are normalized to L_ref. 
     ! They need to be normalized here to rho_ref, since the larmor 
     ! radius is normalized to rho_ref
     if (yx_order) then
        sqrtgxx(:,k) = SQRT(geom%gjj(:,pj1,k))
        prho_gdz_1(:,k) = geom%gjz(:,pj1,k)*rho_Lref/sqrtgxx(:,k)
        prho_gdz_2(:,k) = (geom%gjj(:,pj1,k)*geom%giz(:,pj1,k)-geom%gij(:,pj1,k)*geom%gjz(:,pj1,k)) &
             &*rho_Lref/(sqrtgdet(:,pj1,k)*sqrtgxx(:,k))
     else
        sqrtgxx(:,k) = SQRT(geom%gii(:,pj1,k))
        prho_gdz_1(:,k) = geom%giz(:,pj1,k)*rho_Lref/sqrtgxx(:,k)
        prho_gdz_2(:,k) = (geom%gii(:,pj1,k)*geom%gjz(:,pj1,k)-geom%gij(:,pj1,k)*geom%giz(:,pj1,k)) &
             &*rho_Lref/(sqrtgdet(:,pj1,k)*sqrtgxx(:,k))
     endif
     prho_gdphi(:,k) = - sqrtgxx(:,k)/sqrtgdet(:,pj1,k)*C_y(:)/&
          &major_R**2*rho_Lref
     prho_gdy_1(:,k) = geom%gij(:,pj1,k)/sqrtgxx(:,k)
     prho_gdy_2(:,k) = sqrtgdet(:,pj1,k)/sqrtgxx(:,k)
  Enddo

  !local_sum = 0.0D0
  DO ispec=ln1,ln2
     DO m=lm1,lm2
        DO k=lk1,lk2
           rho(pi1gl:pi2gl) = sqrt(2.0*spec(ispec)%mass*spec(ispec)%temp*mu(m)&
                /spec(ispec)%charge**2/geom%Bfield(pi1gl:pi2gl,pj1,k))
           DO j=lj1,lj2
              CALL set_zero(sparse_gyromatrix(j,k,m,ispec))
#ifdef GDAGGER_G
              CALL set_zero(sparse_gyromatrix_dagger(j,k,m,ispec))
#endif
              DO iDeriv=1,nDeriv
                 CALL set_zero(integral(iDeriv))
              END DO
              ! Attention, the indices for the integral dimension are running from 1 to ni0
              ! not from li1 to li2 !!!
              DO ix=1,ni0
                 pii=ix-1+pi1gl !li1
                 IF (y_local) THEN
                    ntheta = MAX(40,CEILING(sqrtgxx(pii,k)*rho(pii)/deli * 40))
                 ELSE
                    ntheta = MAX(40,CEILING(MAX(ABS(prho_gdy_1(pii,k)),&
                         ABS(prho_gdy_2(pii,k)))*rho(pii)/deli*40))
                 ENDIF
                 dtheta = 2*pi/ntheta

                 
                 allocate(rho_gdx(0:ntheta-1),rho_gdy(0:ntheta-1),rho_gdz(0:ntheta-1),rho_gdphi(0:ntheta-1))
                 if (ga_spatial_var) then 
                    allocate(rho_x(0:ntheta-1), rho_y(0:ntheta-1))
                 endif

                 if (ga_spatial_var) loc_bandwidth=0.
                 DO itheta=0,ntheta-1
                    theta = itheta*dtheta
                    if (ga_spatial_var) then
                       call compute_part_position_varmetric(theta,pii,k,rho(pii),sqrtgxx&
                            ,rho_x(itheta),rho_y(itheta),loc_bandwidth)
                    else
                       rhocostheta = rho(pii)*cos(theta)              
                       rhosintheta = rho(pii)*sin(theta)              
                       
                       rho_gdx(itheta) = sqrtgxx(pii,k) * rhocostheta         ! rho.grad(x)
                       rho_gdy(itheta) = prho_gdy_1(pii,k) * rhocostheta - &
                            & prho_gdy_2(pii,k) * rhosintheta                 ! rho.grad(y)
                       rho_gdz(itheta) = prho_gdz_1(pii,k) * rhocostheta - &
                            & prho_gdz_2(pii,k) * rhosintheta                 ! rho.grad(z)
                       rho_gdphi(itheta) = prho_gdphi(pii,k) * rhosintheta
                    endif
                 END DO
                 !invN = 1.0D0/ntheta

                 if (.not.ga_spatial_var) then
                    if (.not.gyroav_in_xi_eta) then
                       if (y_local) then
                          loc_bandwidth=ceiling(rho_gdx(0)/deli)
                       else
                          loc_bandwidth=ceiling(maxval(abs(rho_gdy))/deli)
                       endif
                    else
                       IF (.not.y_local) stop 'gyroav_in_xi_eta not implemented yet for y_global'
                       loc_bandwidth=ceiling(rho_gdx(0)/deli+0.5*ix*rho_gdz(0)**2)
                       !WRITE(*,"(6I3,I4)") mype,ispec,m,k,j,ix,loc_bandwidth
                    endif
                 endif
                 startbase = ix-loc_bandwidth
                 endbase   = ix+loc_bandwidth

                 IF (startbase.LT.1) THEN
                    startbase = 1
                    IF (rad_bc_type.eq.0) endbase=ni0 !periodic boundaries
                 ENDIF
                 IF (endbase.GT.ni0) THEN
                    endbase=ni0
                    IF (rad_bc_type.eq.0) startbase=1 !periodic boundaries
                 ENDIF

                 ! Due to the multiplication with the derivative matrix,
                 ! which is also banded with a derivative_order/2 upper
                 ! diagonals and the same number of lower diagonals, we
                 ! have to extend the bandwidth for the multiplication:
                 bandwidth(k,m,ispec)=max(bandwidth(k,m,ispec),loc_bandwidth+&
                      &derivative_order/2)

                 DO ibase=startbase,endbase 
                    sum_over_theta = CMPLX(0.0,0.0,KIND(sum_over_theta))
                    ndtheta = 0
                    DO itheta = 0,ntheta-1
                       if (.not.ga_spatial_var) then
                          ! compute particle position
                          call compute_part_position(xprime,yprime_y,my_xval(ix),zval(k),&
                               &rho_gdx(itheta),rho_gdy(itheta),rho_gdz(itheta),&
                               &rho_gdphi(itheta),q_prof(pii), &
                               &dqdx_prof(pii),C_y(pii))
                          
                          !"mirror" function values at (inner) boundary for von Neumann bc
                          !on kj=0 mode
                          if ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(j==lj1)) then
                             if (xprime.lt.my_xval(1)) then
                                xprime = my_xval(1)+(my_xval(1)-xprime)
                             endif
                          endif
                       else
                          if (y_local) then
                             xprime=my_xval(ix)+rho_x(itheta)
                             yprime_y=rho_y(itheta)
                          else
                             xprime=my_xval(ix)+rho_y(itheta)
                             yprime_y=rho_x(itheta)
                          endif
                       endif
                       CALL get_value(basex,ibase,xprime, res)

                       DO ideriv=0,nderiv
                          sum_over_theta(ideriv) = sum_over_theta(ideriv) &
                               & + res(ideriv) * EXP(imag*kj(j)*yprime_y)
                       END DO
                    END DO !itheta

                    sum_over_theta = sum_over_theta/ntheta
                    Call set_value(sparse_gyromatrix(j,k,m,ispec),ix,ibase,sum_over_theta(0))
                    DO ideriv=1,nderiv
                       CALL set_value(integral(ideriv),ix,ibase,sum_over_theta(ideriv))
                    ENDDO
      
                 END DO !ibase 

                 deallocate(rho_gdx,rho_gdy,rho_gdz,rho_gdphi)
                 if (ga_spatial_var) deallocate(rho_x,rho_y)
              END DO !ix

              CALL commit_values(sparse_gyromatrix(j,k,m,ispec))

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

                 CALL add_matrix(sparse_gyromatrix(j,k,m,ispec), tempmat)
              END DO
                    
#ifdef GDAGGER_G
              call transpose_and_conjugate(sparse_gyromatrix(j,k,m,ispec),&
                   temp_dagger_matrix)
              call transpose_storage(temp_dagger_matrix,sparse_gyromatrix_dagger(j,k,m,ispec))
#endif

           END DO !j
        END DO ! k
     END DO !m
  END DO ! ispec
  deallocate(sqrtgxx, prho_gdy_1,prho_gdy_2,prho_gdz_1,prho_gdz_2,prho_gdphi)
  deallocate(sum_over_theta)
#ifdef GDAGGER_G
  call finalize(temp_dagger_matrix)
#endif
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
  TYPE(BandedMatrix),INTENT(IN) :: mat_in
  TYPE(BandedMatrix),INTENT(INOUT) :: mat_dagger

  INTEGER, INTENT(IN) :: bandwidth
  !local variables
  INTEGER :: ix, ibase, startbase, endbase
  COMPLEX, DIMENSION(1:ni0,1:ni0) :: temp_fullmat
  
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
             &CONJG(temp_fullmat(ix,ibase)))
     ENDDO
  ENDDO
  
  CALL commit_values(mat_dagger)

END SUBROUTINE set_dagger_matrix
#endif

!>Clean up all arrays etc. employed by this module
SUBROUTINE finalize_gyro_average_df
  INTEGER :: j,k,m,n

  ! finalize the xgrid
  call finalize(xgrid)

  DEALLOCATE(my_xval,sqrtgdet)

  CALL finalize(vec_v)
  call finalize(res_v)
  CALL finalize_localpolynombase(basex)

  IF (ALLOCATED(sparse_gyromatrix)) THEN
     DO n=ln1,ln2
        DO m=lm1,lm2
           DO k=lk1,lk2
              DO j=lj1,lj2
                 CALL finalize(sparse_gyromatrix(j,k,m,n))
#ifdef GDAGGER_G
                 CALL finalize(sparse_gyromatrix_dagger(j,k,m,n))
#endif
              END DO
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

END SUBROUTINE finalize_gyro_average_df


SUBROUTINE compute_part_position(ao_xprime,ao_yprime_y, ai_x, ai_z, ai_rhogdx, ai_rhogdy, ai_rhogdz,&
     &ai_rhogdphi,ai_q,ai_dqdx,ai_Cy)
  REAL, INTENT(out) :: ao_xprime,ao_yprime_y
  REAL, INTENT(in)  :: ai_x, ai_z,ai_rhogdx,ai_rhogdy, ai_rhogdz, &
       &ai_rhogdphi,ai_q,ai_dqdx,ai_Cy
  REAL              :: zprime,qxprime
  
  IF (.NOT.gyroav_in_xi_eta) THEN ! using linearisation of the x,y,z metric
     IF (y_local) THEN
        ao_xprime = ai_x + ai_rhogdx
        ao_yprime_y = ai_rhogdy
     ELSE
        ao_xprime = ai_x + ai_rhogdy
        ao_yprime_y = ai_rhogdx
     ENDIF
  ELSE
     IF (.not.y_local) stop 'gyroav_in_xi_eta not implemented yet for y_global'
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


!!\Todo Do we need 100 integration points? 
!!can we choose a fixed ntheta and pull this routine out of the j loop? 
subroutine compute_part_position_varmetric(theta,pii,k,rho,sqrtgxx,rho_x,rho_y,loc_bandwidth)
  integer, intent(in):: pii, k
  real, intent(in):: theta, rho
  real, dimension(pi1gl:pi2gl,lk1:lk2):: sqrtgxx
  real, intent(out):: rho_x, rho_y
  integer, intent(inout):: loc_bandwidth
  integer:: n_integral, interp=1, pi_low, pi_high, i_int, i
  real:: drho, pi_int, dist_left, dist_right, g_fac1, g_fac2, g_fac3, drho_x, drho_y
  real, dimension(1):: gij_int, sqrtgxx_int, sqrtgdet_int, pi_arr_out
  real, dimension(pi1gl:pi2gl):: pi_arr_in
  
  n_integral=100

!interp:  0=linear, 1=Lagrange

  drho=rho/(n_integral-1)
  pi_int=pii
  rho_x=0.
  rho_y=0.
  do i=pi1gl,pi2gl
     pi_arr_in(i)=real(i)
  enddo
  do i_int=1,n_integral-1
     if (interp==0) then
        pi_low=floor(pi_int)
        pi_high=ceiling(pi_int)
        if (rad_bc_type==0) then
           pi_low=modulo(floor(pi_int),ni0)
           pi_high=modulo(ceiling(pi_int),ni0)
        else
           pi_low=min(ni0-1,pi_low)
           pi_high=min(ni0-1,pi_high)
           pi_low=max(0,pi_low)
           pi_high=max(0,pi_high)
        endif
        dist_left=pi_int-pi_low
        dist_right=pi_high-pi_int
        if (dist_left.lt.0) then
           dist_left=dist_left+ni0
        endif
        if (dist_right.lt.0) then
           dist_right=dist_right+ni0
        endif
        if (pi_low==pi_high) then
           g_fac1=geom%gij(pi_low,pj1,k)/sqrtgxx(pi_low,k)
           g_fac2=sqrtgdet(pi_low,pj1,k)/sqrtgxx(pi_low,k)
           g_fac3=sqrtgxx(pi_low,k)
        else
           g_fac1=dist_left*geom%gij(pi_high,pj1,k)/sqrtgxx(pi_high,k)+dist_right*geom%gij(pi_low,pj1,k)/sqrtgxx(pi_low,k)
           g_fac2=dist_left*sqrtgdet(pi_high,pj1,k)/sqrtgxx(pi_high,k)+dist_right*sqrtgdet(pi_low,pj1,k)/sqrtgxx(pi_low,k)
           g_fac3=dist_left*sqrtgxx(pi_high,k)+dist_right*sqrtgxx(pi_low,k)
        endif
     elseif (interp==1) then
        pi_arr_out(1)=pi_int
        call lag3interp(geom%gij(:,pj1,k),pi_arr_in,ni0,gij_int,pi_arr_out,1)
        call lag3interp(sqrtgxx(:,k),pi_arr_in,ni0,sqrtgxx_int,pi_arr_out,1)
        call lag3interp(sqrtgdet(:,pj1,k),pi_arr_in,ni0,sqrtgdet_int,pi_arr_out,1)
        g_fac1=gij_int(1)/sqrtgxx_int(1)
        g_fac2=sqrtgdet_int(1)/sqrtgxx_int(1)
        g_fac3=sqrtgxx_int(1)
     else
        stop 'invalid interpolation type in gyroaverage'
     endif
     drho_y = drho * (cos(theta) * g_fac1 - sin(theta) * g_fac2)
     drho_x = drho * cos(theta) * g_fac3
     rho_y=rho_y+drho_y
     rho_x=rho_x+drho_x
     if (y_local) then
        pi_int=pi_int+drho_x/deli
     else
        pi_int=pi_int+drho_y/deli
     endif
     if (rad_bc_type==0) then
        if (pi_int.ge.ni0) pi_int=pi_int-ni0
        if (pi_int.lt.0) pi_int=pi_int+ni0
     endif
  enddo
  if (y_local) then
     loc_bandwidth=max(loc_bandwidth,ceiling(abs(rho_x)/deli))
  else
     loc_bandwidth=max(loc_bandwidth,ceiling(abs(rho_y)/deli))
  endif
  
end subroutine compute_part_position_varmetric

END MODULE gyro_average_df_mod
