PROGRAM test_boundary_exchange
  use mpi
  use mpi_environment
  USE boundary_exchange
  use BoundaryDescriptionModule
  IMPLICIT NONE

  LOGICAL :: passed_all_tests=.TRUE.
  integer :: ierr
  INTEGER :: MY_MPI_COMM_WORLD, comm_cart
  LOGICAL :: wtime_is_global,pe_periods(6)
  integer :: pe_dims(6)
  INTEGER :: communicators(N_COMM),mype, n_procs

  ! setup MPI
  pe_dims=(/1,1,1,1,1,8/)
  pe_periods=.FALSE.
  pe_periods(4)=.TRUE. ! set z direction to periodic
  !pe_periods(6)=.TRUE. ! set x direction to periodic
  CALL initialize_mpi_environment(pe_dims,MY_MPI_COMM_WORLD,wtime_is_global)
  IF (MY_MPI_COMM_WORLD.NE.MPI_COMM_NULL) THEN
     CALL split_comm(MY_MPI_COMM_WORLD,pe_dims,pe_periods,comm_cart,communicators)
     CALL mpi_comm_rank(comm_cart, mype, ierr)
     CALL mpi_comm_size(comm_cart, n_procs, ierr)
     PRINT*,"This is process ",mype," of ",n_procs," processes."

     !-----------------------------------------
     ! Testing interface exchange_x
     !-----------------------------------------
     if (test_exchange_x_1D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_1D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_1D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if

     if (test_exchange_x_2D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_2D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_2D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if

     if (test_exchange_x_3D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_3D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_3D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if

     if (test_exchange_x_4D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_4D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_4D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if

     if (test_exchange_x_5D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_5D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_5D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if

     if (test_exchange_x_6D()) then
        !write(*,"(A40,I3,A10)") "test_exchange_x_6D on process ",mype," passed"
     else
        write(*,"(A40,I3,A10)") "test_exchange_x_6D on process ",mype,"  FAILED"
        passed_all_tests=.false.
     end if


     IF (passed_all_tests) THEN
        print*,'All tests on process ',mype,' passed.'
     else
        print*,'Some tests on process ',mype,' FAILED.'
     END IF
  END IF
  CALL finalize_mpi_environment

contains

  FUNCTION test_exchange_x_1D() RESULT(passed)
    LOGICAL :: passed

    ! Local variables
    INTEGER :: nx0, nxb, lbx, ubx, li1, li2, li0, lx0
    INTEGER ::  ierr
    COMPLEX, ALLOCATABLE, DIMENSION(:) :: array
    type(BoundaryDescription) :: bdesc
    INTEGER :: pe_coords(N_PE_DIMS), i, i_proc, rad_bc_type
    logical :: output=.false.
    
    passed=.TRUE.

    CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,pe_coords,ierr)
    ! pe_coords contains the process coordinates in the Cartesion topology
    ! pe_coords(1) -> my_pespec, ..., pe_coords(6) -> my_pex

    ! initialize exchange
    !PRINT*,mype,"Calling initialize_boundary_exchange"
    CALL initialize_boundary_exchange(comm_cart, rad_bc_type, communicators)

    ! initialize array
    nx0 = 32

    ! specify boundary widths
    nxb = 2

    li0 = nx0 / pe_dims(6)
    li1 = pe_coords(6)*li0+1
    lbx = li1-nxb; li2 = li1+li0-1; ubx = li2+nxb
    lx0 = ubx-lbx+1

    rad_bc_type = 1 ! Dirichlet boundary condition


    ! 1d array
    allocate(array(lbx:ubx))
    array = CMPLX(0,0,KIND(array))
    DO i=li1,li2
       array(i) = CMPLX(i,0,kind(array))
    END DO

    CALL initialize_type(bdesc,lx0,nxb,nxb,1,5)
    call set_mpi_type(bdesc)

    CALL exchange_x(bdesc,array)

    ! check the result
    ! check lower boundary
    IF (mype.EQ.0) THEN
       DO i=lbx,li1-1
          IF (ABS(array(i)) .GT.1e-15) passed=.FALSE.
       END DO
    ELSE
       DO i=lbx,li1-1
          IF (ABS(array(i)-CMPLX(i,0,KIND(array))).GT.1e-15) passed=.FALSE.
       END DO
    END IF

    ! check upper boundary
    IF (mype.EQ.n_procs-1) THEN
       DO i=li2+1,ubx
          IF (ABS(array(i)).GT.1e-15) passed=.FALSE.
       END DO
    ELSE
       DO i=li2+1,ubx
          IF (ABS(array(i)-CMPLX(i,0,KIND(array))).GT.1e-15) passed=.FALSE.
       END DO
    END IF

    call finalize_type(bdesc)
    call finalize_boundary_exchange

end function test_exchange_x_1D

function test_exchange_x_2D() result(passed)
	logical :: passed

    ! Local variables
    integer :: nx0, nxb, lbx, ubx, li1, li2, li0, lx0
    INTEGER :: ny0, nyb, lby, uby, lj1, lj2, lj0, ly0
    INTEGER ::  ierr
    COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: array
    type(BoundaryDescription) :: bdesc
    INTEGER :: pe_coords(N_PE_DIMS), i,j, i_proc, rad_bc_type
    logical :: output=.false.

    passed = .TRUE.
    CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,pe_coords,ierr)
    ! pe_coords contains the process coordinates in the Cartesion topology
    ! pe_coords(1) -> my_pespec, ..., pe_coords(6) -> my_pex

    ! initialize exchange
    !PRINT*,mype,"Calling initialize_boundary_exchange"
    CALL initialize_boundary_exchange(comm_cart, rad_bc_type, communicators)

    ! initialize array
    nx0 = 32; ny0 = 4

    ! specify boundary widths
    nxb = 2; nyb = 0

    li0 = nx0 / pe_dims(6)
    li1 = pe_coords(6)*li0+1
    lbx = li1-nxb; li2 = li1+li0-1; ubx = li2+nxb
    lx0 = ubx-lbx+1

    lj0 = ny0 / pe_dims(5)
    lj1 = pe_coords(5)*lj0+1
    lby = lj1-nyb; lj2 = lj1+lj0-1; uby = lj2+nyb
    ly0 = uby-lby+1

    rad_bc_type = 1 ! Dirichlet boundary condition

    ! test 2D array
    ALLOCATE(array(lbx:ubx,lby:uby))
    array = CMPLX(0.0,0.0,KIND(array))

    DO j=lj1,lj2
       DO i=li1,li2
          array(i,j) = CMPLX(i,j,kind(array))
       END DO
    END DO

    !PRINT*,mype,"Initializing the BoundaryDescription"
    ! set the BoundaryDescription
    CALL initialize_type(bdesc,lx0,nxb,nxb,ly0,5)
    CALL set_mpi_type(bdesc)

    CALL exchange_x(bdesc,array)

    ! check the result
    ! check lower boundary
    IF (mype.EQ.0) THEN
       DO j=lj1,lj2
          DO i=lbx,li1-1
             !PRINT*,i,j,ABS(array(i,j))
             IF (ABS(array(i,j)) .GT.1e-15) passed=.FALSE.
          END DO
       END DO
    ELSE
       DO j=lj1,lj2
          DO i=lbx,li1-1
             !PRINT*,i,j,ABS(array(i,j)-CMPLX(i,j,KIND(array)))
             IF (ABS(array(i,j)-CMPLX(i,j,kind(array))).GT.1e-15) passed=.FALSE.
          END DO
       END DO
    END IF

    ! check upper boundary
    IF (mype.EQ.n_procs-1) THEN
       DO j=lj1,lj2
          DO i=li2+1,ubx
             IF (ABS(array(i,j)).GT.1e-15) passed=.false.
          END DO
       END DO
    ELSE
       DO j=lj1,lj2
          DO i=li2+1,ubx
             !PRINT*,i,j,ABS(array(i,j)-CMPLX(i,j,KIND(array)))
             IF (ABS(array(i,j)-CMPLX(i,j,kind(array))).GT.1e-15) passed=.FALSE.
          END DO
       END DO
    END IF

    call finalize_type(bdesc)
    call finalize_boundary_exchange

end function test_exchange_x_2D

function test_exchange_x_3D() result(passed)
  LOGICAL :: passed

    ! Local variables
    integer :: nx0, nxb, lbx, ubx, li1, li2, li0, lx0
    INTEGER :: ny0, nyb, lby, uby, lj1, lj2, lj0, ly0
    INTEGER :: nz0, nzb, lbz, ubz, lk1, lk2, lk0, lz0
    INTEGER ::  ierr
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: array
    type(BoundaryDescription) :: bdesc
    INTEGER :: pe_coords(N_PE_DIMS), i,j, k,i_proc, rad_bc_type
    logical :: output=.false.

	passed=.true.
    CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,pe_coords,ierr)
    ! pe_coords contains the process coordinates in the Cartesion topology
    ! pe_coords(1) -> my_pespec, ..., pe_coords(6) -> my_pex

    ! initialize exchange
    !PRINT*,mype,"Calling initialize_boundary_exchange"
    CALL initialize_boundary_exchange(comm_cart, rad_bc_type, communicators)

    ! initialize array
    nx0 = 32; ny0 = 4; nz0 = 8

    ! specify boundary widths
    nxb = 2; nyb = 0; nzb = 2

    li0 = nx0 / pe_dims(6)
    li1 = pe_coords(6)*li0+1
    lbx = li1-nxb; li2 = li1+li0-1; ubx = li2+nxb
    lx0 = ubx-lbx+1

    lj0 = ny0 / pe_dims(5)
    lj1 = pe_coords(5)*lj0+1
    lby = lj1-nyb; lj2 = lj1+lj0-1; uby = lj2+nyb
    ly0 = uby-lby+1

    lk0 = nz0 / pe_dims(4)
    lk1 = pe_coords(4)*lk0+1
    lbz = lk1-nzb; lk2 = lk1+lk0-1; ubz = lk2+nzb
    lz0 = ubz-lbz+1

    rad_bc_type = 1 ! Dirichlet boundary condition


    ! test 3D array
    ALLOCATE(array(lbx:ubx,lby:uby,lbz:ubz))
    array = CMPLX(0.0,0.0,KIND(array))

    DO k=lk1,lk2
       DO j=lj1,lj2
          DO i=li1,li2
             array(i,j,k) = CMPLX((k-lk1)*li0+i,j,kind(array))
          END DO
       END DO
    END DO

    !PRINT*,mype,"Initializing the BoundaryDescription"
    ! set the BoundaryDescription
    CALL initialize_type(bdesc,lx0,nxb,nxb,ly0,5)
    CALL set_mpi_type(bdesc)

    CALL exchange_x(bdesc,array)

    ! check the result
    ! check lower boundary
    IF (mype.EQ.0) THEN
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=lbx,li1-1
                IF (ABS(array(i,j,k)) .GT.1e-15) passed=.FALSE.
             END DO
          END DO
       END DO
    ELSE
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=lbx,li1-1
                IF (ABS(array(i,j,k)-CMPLX((k-lk1)*li0+i,j,KIND(array))).GT.1e-15) passed=.FALSE.
             END DO
          END DO
       END DO
    END IF

    ! check upper boundary
    IF (mype.EQ.n_procs-1) THEN
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li2+1,ubx
                IF (ABS(array(i,j,k)).GT.1e-15) passed=.FALSE.
             END DO
          END DO
       END DO
    ELSE
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li2+1,ubx
                IF (ABS(array(i,j,k)-CMPLX((k-lk1)*li0+i,j,KIND(array))).GT.1e-15) passed=.FALSE.
             END DO
          END DO
       END DO
    END IF

    CALL finalize_type(bdesc)
    CALL finalize_boundary_exchange
end function test_exchange_x_3D

function test_exchange_x_4D() result(passed)
	logical :: passed

    ! Local variables
    integer :: nx0, nxb, lbx, ubx, li1, li2, li0, lx0
    INTEGER :: ny0, nyb, lby, uby, lj1, lj2, lj0, ly0
    INTEGER :: nz0, nzb, lbz, ubz, lk1, lk2, lk0, lz0
    INTEGER :: nv0, nvb, lbv, ubv, ll1, ll2, ll0, lv0
    INTEGER :: nw0, nwb, lbw, ubw, lm1, lm2, lm0, lw0
    INTEGER :: nspec, ln1, ln2, ln0
    INTEGER ::  ierr
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: array
    type(BoundaryDescription) :: bdesc
    INTEGER :: pe_coords(N_PE_DIMS), i,j, k,l,i_proc, rad_bc_type
    logical :: output=.false.

	passed=.true.
    CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,pe_coords,ierr)
    ! pe_coords contains the process coordinates in the Cartesion topology
    ! pe_coords(1) -> my_pespec, ..., pe_coords(6) -> my_pex

    ! initialize exchange
    !PRINT*,mype,"Calling initialize_boundary_exchange"
    CALL initialize_boundary_exchange(comm_cart, rad_bc_type, communicators)

    ! initialize array
    nx0 = 32; ny0 = 4; nz0 = 8; nv0 = 4; nw0 = 4; nspec = 2

    ! specify boundary widths
    nxb = 2; nyb = 0; nzb = 2; nvb = 1; nwb = 1

    li0 = nx0 / pe_dims(6)
    li1 = pe_coords(6)*li0+1
    lbx = li1-nxb; li2 = li1+li0-1; ubx = li2+nxb
    lx0 = ubx-lbx+1
    !WRITE(*,"(I3,2I4,4(A,I3),A)") mype,lx0,li0,&
    !     &"outer = (",lbx,":",ubx,"), inner = (",li1,":",li2,")"

    lj0 = ny0 / pe_dims(5)
    lj1 = pe_coords(5)*lj0+1
    lby = lj1-nyb; lj2 = lj1+lj0-1; uby = lj2+nyb
    ly0 = uby-lby+1

    lk0 = nz0 / pe_dims(4)
    lk1 = pe_coords(4)*lk0+1
    lbz = lk1-nzb; lk2 = lk1+lk0-1; ubz = lk2+nzb
    lz0 = ubz-lbz+1

    ll0 = nv0 / pe_dims(3)
    ll1 = pe_coords(3)*ll0+1
    lbv = ll1-nvb; ll2 = ll1+ll0-1; ubv = ll2+nvb
    lv0 = ubv-lbv+1

    lm0 = nw0 / pe_dims(2)
    lm1 = pe_coords(2)*lm0+1
    lbw = lm1-nwb; lm2 = lm1+lm0-1; ubw = lm2+nwb
    lw0 = ubw-lbw+1

    ln0 = nspec / pe_dims(1)
    ln1 = pe_coords(1)*ln0+1
    ln2 = ln1+ln0-1

    rad_bc_type = 1 ! Dirichlet boundary condition


    ! test 3D array
    ALLOCATE(array(lbx:ubx,lby:uby,lbz:ubz,lbv:ubv))
    array = CMPLX(0.0,0.0,KIND(array))

    DO l=ll1,ll2
       DO k=lk1,lk2
          DO j=lj1,lj2
             DO i=li1,li2
                array(i,j,k,l) = CMPLX((k-lk1)*li0+i,(l-ll1)*lj0+j,kind(array))
             END DO
          END DO
       END DO
    END DO

    !PRINT*,mype,"Initializing the BoundaryDescription"
    ! set the BoundaryDescription
    CALL initialize_type(bdesc,lx0,nxb,nxb,ly0*lz0*lv0,5)
    CALL set_mpi_type(bdesc)

    CALL exchange_x(bdesc,array)

    ! check the result
    ! check lower boundary
    IF (mype.EQ.0) THEN
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=lbx,li1-1
                   IF (ABS(array(i,j,k,l)) .GT.1e-15) passed=.FALSE.
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=lbx,li1-1
                   IF (ABS(array(i,j,k,l)-CMPLX((k-lk1)*li0+i,(l-ll1)*lj0+j,KIND(array))).GT.1e-15) &
                        &passed=.FALSE.
                END DO
             END DO
          END DO
       END DO
    END IF

    ! check upper boundary
    IF (mype.EQ.n_procs-1) THEN
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=li2+1,ubx
                   IF (ABS(array(i,j,k,l)).GT.1e-15) passed=.FALSE.
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=li2+1,ubx
                   IF (ABS(array(i,j,k,l)-CMPLX((k-lk1)*li0+i,(l-ll1)*lj0+j,KIND(array))).GT.1e-15) &
                        &passed=.FALSE.
                END DO
             END DO
          END DO
       END DO
    END IF

    CALL finalize_type(bdesc)
    CALL finalize_boundary_exchange
end function test_exchange_x_4D

function test_exchange_x_5D() result(passed)
	logical :: passed

    ! Local variables
    integer :: nx0, nxb, lbx, ubx, li1, li2, li0, lx0
    INTEGER :: ny0, nyb, lby, uby, lj1, lj2, lj0, ly0
    INTEGER :: nz0, nzb, lbz, ubz, lk1, lk2, lk0, lz0
    INTEGER :: nv0, nvb, lbv, ubv, ll1, ll2, ll0, lv0
    INTEGER :: nw0, nwb, lbw, ubw, lm1, lm2, lm0, lw0
    INTEGER :: nspec, ln1, ln2, ln0
    INTEGER ::  ierr
    COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: array
    type(BoundaryDescription) :: bdesc
    INTEGER :: pe_coords(N_PE_DIMS), i,j, k,l,m,i_proc, rad_bc_type
    logical :: output=.false.

	passed=.true.
    CALL mpi_cart_coords(comm_cart,mype,N_PE_DIMS,pe_coords,ierr)
    ! pe_coords contains the process coordinates in the Cartesion topology
    ! pe_coords(1) -> my_pespec, ..., pe_coords(6) -> my_pex

    ! initialize exchange
    !PRINT*,mype,"Calling initialize_boundary_exchange"
    CALL initialize_boundary_exchange(comm_cart, rad_bc_type, communicators)

    ! initialize array
    nx0 = 32; ny0 = 4; nz0 = 8; nv0 = 4; nw0 = 4; nspec = 2

    ! specify boundary widths
    nxb = 2; nyb = 0; nzb = 2; nvb = 1; nwb = 1

    li0 = nx0 / pe_dims(6)
    li1 = pe_coords(6)*li0+1
    lbx = li1-nxb; li2 = li1+li0-1; ubx = li2+nxb
    lx0 = ubx-lbx+1
    !WRITE(*,"(I3,2I4,4(A,I3),A)") mype,lx0,li0,&
    !     &"outer = (",lbx,":",ubx,"), inner = (",li1,":",li2,")"

    lj0 = ny0 / pe_dims(5)
    lj1 = pe_coords(5)*lj0+1
    lby = lj1-nyb; lj2 = lj1+lj0-1; uby = lj2+nyb
    ly0 = uby-lby+1

    lk0 = nz0 / pe_dims(4)
    lk1 = pe_coords(4)*lk0+1
    lbz = lk1-nzb; lk2 = lk1+lk0-1; ubz = lk2+nzb
    lz0 = ubz-lbz+1

    ll0 = nv0 / pe_dims(3)
    ll1 = pe_coords(3)*ll0+1
    lbv = ll1-nvb; ll2 = ll1+ll0-1; ubv = ll2+nvb
    lv0 = ubv-lbv+1

    lm0 = nw0 / pe_dims(2)
    lm1 = pe_coords(2)*lm0+1
    lbw = lm1-nwb; lm2 = lm1+lm0-1; ubw = lm2+nwb
    lw0 = ubw-lbw+1

    ln0 = nspec / pe_dims(1)
    ln1 = pe_coords(1)*ln0+1
    ln2 = ln1+ln0-1

    rad_bc_type = 1 ! Dirichlet boundary condition


    ! test 5D array
    ALLOCATE(array(lbx:ubx,lby:uby,lbz:ubz,lbv:ubv,lbw:ubw))
    array = CMPLX(0.0,0.0,KIND(array))

    DO m=lm1,lm2
       DO l=ll1,ll2
          DO k=lk1,lk2
             DO j=lj1,lj2
                DO i=li1,li2
                   array(i,j,k,l,m) = CMPLX((m-lm1)*lk0*li0+(k-lk1)*li0+i,(l-ll1)*lj0+j,KIND(array))
                END DO
             END DO
          END DO
       END DO
    END DO

    !PRINT*,mype,"Initializing the BoundaryDescription"
    ! set the BoundaryDescription
    CALL initialize_type(bdesc,lx0,nxb,nxb,ly0*lz0*lv0*lw0,5)
    CALL set_mpi_type(bdesc)

    CALL exchange_x(bdesc,array)

    ! check the result
    ! check lower boundary
    IF (mype.EQ.0) THEN
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                DO j=lj1,lj2
                   DO i=lbx,li1-1
                      IF (ABS(array(i,j,k,l,m)) .GT.1e-15) passed=.FALSE.
                   END DO
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                DO j=lj1,lj2
                   DO i=lbx,li1-1
                      IF (ABS(array(i,j,k,l,m)-CMPLX((m-lm1)*lk0*li0+(k-lk1)*li0+i,(l-ll1)*lj0+j,KIND(array))).GT.1e-15) passed=.FALSE.
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF

    ! check upper boundary
    IF (mype.EQ.n_procs-1) THEN
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                DO j=lj1,lj2
                   DO i=li2+1,ubx
                      IF (ABS(array(i,j,k,l,m)).GT.1e-15) passed=.FALSE.
                   END DO
                END DO
             END DO
          END DO
       END DO
    ELSE
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                DO j=lj1,lj2
                   DO i=li2+1,ubx
                      IF (ABS(array(i,j,k,l,m)-CMPLX((m-lm1)*lk0*li0+(k-lk1)*li0+i,(l-ll1)*lj0+j,KIND(array))).GT.1e-15) passed=.FALSE.
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF

    CALL finalize_type(bdesc)
    CALL finalize_boundary_exchange
end function test_exchange_x_5D

function test_exchange_x_6D() result(passed)
	logical :: passed

 passed = .TRUE.
end function test_exchange_x_6D

function test_mp_initialize_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_initialize_type

function test_mp_finalize_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_finalize_type

function test_mp_set_mpi_type() result(passed)
	logical :: passed

	passed=.true.
end function test_mp_set_mpi_type
end program
