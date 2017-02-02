!>ITM interface library
!!Compile with gmake -f ../makefile lib LIBN=libITMgene
Subroutine ITMgene(eq, coreprof, coretransp, code_parameters)
  USE mpi
  USE ITM_Types
  USE ITM_Constants
  USE Euitm_schemas
  Use gene_subroutine
  Use libITMgene_aux
  Use par_geom
#ifdef WITHSLEPC
  use slepc_aux
#endif

  Implicit None

  TYPE (type_equilibrium), pointer :: eq(:)
  TYPE (type_coreprof), pointer :: coreprof(:)
  TYPE (type_coretransp), pointer :: coretransp(:)
  TYPE (type_param) :: code_parameters

  INTEGER(ITM_I4) :: i,jm,j0,j1,j2
  INTEGER(ITM_I4) :: nrho_prof,nrho_eq,nion_prof
  REAL(R8) :: x,aintm,aint0,aint1,aint2,xm,ITM_x0,x1,x2
  REAL(R8) :: a00,b00,r00,hra
  REAL(R8) :: rho_ref,c_ref,zeff,qq

  REAL(R8), DIMENSION(:), ALLOCATABLE :: flux_eq, angle_eq

  integer(ITM_I4) :: return_status

  character(len = 132), target :: codename(1) = 'GENE11'
  character(len = 132), target :: codeversion(1) = '1.4'
  character(len = 32) :: cpofile

  integer(ITM_I4) :: file_length, i_line, n_lines
  integer(ITM_I4) :: ios
  
  character(len=128):: ITM_par_in_dir='skip_parfile'
  character(len=128):: ITM_file_extension, &
       &ITM_checkpoint_in=''

  integer :: n, nky0_sav, ierr, n_parallel_fluxtubes
  real :: nexc_sav, kymin_sav, ql_diff_max

  real, dimension(:,:), allocatable :: gene_fluxes

  !geom_container for the full 2D geometry parameters
  TYPE(geomtype) :: geom_container

#if defined(SVN_REV)
!  write(codeversion(1),"(4A)") codeversion(1),'(',SVN_REV,')'
#endif

!  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,ierr)

#if defined(WITHSLEPC)
  if (.not.slepc_restartable()) &
       call my_SlepcInitialize(MPI_COMM_WORLD,.true.)
#endif

!...  find grid size for profiles
!...  default sets grid size to that of profiles
!...  find number of ion species also from coreprof

  nrho_prof=SIZE(coreprof(1)%rho_tor) !number of radial input grid points
  nion_prof=SIZE(coreprof(1)%ni%value)/nrho_prof  !number of ions
  nrho_eq=SIZE(eq(1)%profiles_1d%rho_vol)

!...  allocations

  IF (.NOT. ASSOCIATED(coretransp)) THEN
     allocate(coretransp(1))
     allocate(coretransp(1)%codeparam%codename(1))
     allocate(coretransp(1)%codeparam%codeversion(1))
     if (.not. associated(code_parameters%parameters)) then
        write(*,*) 'ERROR: GENE parameters not associated!'
        stop
     else
        allocate(coretransp(1)%codeparam%parameters(size( &
             code_parameters%parameters)))
     end if

     !-- add to coretransp
     coretransp(1)%codeparam%codename = codename
     coretransp(1)%codeparam%codeversion = codeversion
     coretransp(1)%codeparam%parameters = code_parameters%parameters

     IF (nrho_transp == 0) then
        nrho_transp = nrho_prof
        if (nrho_transp.ne.1) nrho_transp = 4
!        ! ... for testing, we use a reduced (1/5) grid
!        nrho_transp = nrho_transp/5
!        ! ... for benchmark we use 4 points: .4,.5,.6,.7 *minor_r
!        !or 1 point the .7
!        nrho_transp = 1
     ENDIF
     IF (nion == 0) nion = nion_prof

     n_spec = nion + 1
     allocate(spec(0:n_spec-1))
     allocate(gene_fluxes(0:n_spec-1,2))

     call set_par_in_defaults
     call set_cm_defaults
     call read_xml_parameters(code_parameters, return_status)
     if (return_status /= 0) then
        write(*,*) 'ERROR: Could not assign GENE parameters!'
        stop
     end if
     if (.not.nonlinear) quasilin_model = 1

     ! initializes what should later be turbulence cpo turbulence_cpo (?)
     CALL Turb_Constructor(coretransp(1), 1, nrho_transp, nion)

     Call initialize_geomtype(geom_container,1,nrho_transp,1,nrho_transp,&
          &1,nrho_transp,1,1,1,1,1,nz0,2)
      
     !...  done initialisation
  END IF

  !...  find grid sizes for transp
  
  nrho_transp=SIZE(coretransp(1)%rho_tor)
  nion=SIZE(coretransp(1)%ti_transp%flux)/nrho_transp

  !constant normalization
  Lref = eq(1)%eqgeometry%a_minor !coreprof(1)%toroid_field%r0
  Bref = coreprof(1)%toroid_field%b0

  a00=eq(1)%eqgeometry%a_minor

  !normalized (constant) geometry parameters
  major_R = coreprof(1)%toroid_field%r0/Lref
  minor_r = eq(1)%eqgeometry%a_minor/Lref

  !...  set transport grid
  if ((nrho_transp.eq.1).and.(flux_pos>0)) then
     coretransp(1)%rho_tor_norm=(/ flux_pos /)
  else
     hra=1.0_R8/nrho_transp
     coretransp(1)%rho_tor_norm=(/ (hra*(i-0.5_R8),i=1,nrho_transp) /)
     
! for the imp4 benchmark  we want four specific radial positions simulated:
! rho_tor = (/0.4*a00,0.5*a00,0.6*a00,0.7*a00/)
!coretransp(1)%rho_tor_norm=(/0.4,0.5,0.6,0.7/)
  endif

  coretransp(1)%rho_tor=coretransp(1)%rho_tor_norm*a00

!...  allocations with coretransp grid

  IF (.NOT. ALLOCATED(rho_tor)) THEN
     ALLOCATE(rho_tor(nrho_transp))
     ALLOCATE(rho_tor_norm(nrho_transp))
     ALLOCATE(nnex(nrho_transp))
     ALLOCATE(ttex(nrho_transp))
     ALLOCATE(nnix(nrho_transp))
     ALLOCATE(ttix(nrho_transp))
     ALLOCATE(zeffx(nrho_transp))
     ALLOCATE(qqx(nrho_transp))
     ALLOCATE(rlnex(nrho_transp))
     ALLOCATE(rltex(nrho_transp))
     ALLOCATE(rltix(nrho_transp))
     ALLOCATE(shatx(nrho_transp))
  END IF

!...  do the interpolation to get parameters and gradients

  rho_tor=coretransp(1)%rho_tor
  rho_tor_norm=coretransp(1)%rho_tor_norm

  CALL L3interp( coreprof(1)%ne%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       nnex, rho_tor_norm, nrho_transp)
  CALL L3interp( coreprof(1)%te%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       ttex, rho_tor_norm, nrho_transp)
  CALL L3interp( coreprof(1)%ni%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       nnix, rho_tor_norm, nrho_transp)
  CALL L3interp( coreprof(1)%ti%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       ttix, rho_tor_norm, nrho_transp)

!take q profile from eq CPO? (-> see changes by mjpuesch in test_libITMgene.F90 )
!must be taken from eq CPO whith rho_vol axis!
  CALL L3interp( eq(1)%profiles_1d%q, &
       eq(1)%profiles_1d%rho_vol, nrho_eq, qqx, rho_tor_norm, nrho_transp)
 
  CALL L3interp( coreprof(1)%profiles1d%zeff%value, &
       coreprof(1)%rho_tor_norm, nrho_prof, zeffx, rho_tor_norm, nrho_transp)

  CALL L3deriv( coreprof(1)%ne%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       rlnex, rho_tor_norm, nrho_transp)
  CALL L3deriv( coreprof(1)%te%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       rltex, rho_tor_norm, nrho_transp)
  CALL L3deriv( coreprof(1)%ti%value, coreprof(1)%rho_tor_norm, nrho_prof, &
       rltix, rho_tor_norm, nrho_transp)

!take q profile from eq CPO? (-> see changes by mjpuesch in test_libITMgene.F90 )
!hkd: yes, q and shat from eq CPO (agreed Sept 2010)
  CALL L3deriv( eq(1)%profiles_1d%q, &
       eq(1)%profiles_1d%rho_vol, nrho_eq, shatx,rho_tor_norm, nrho_transp)
  shatx=shatx*rho_tor_norm/qqx

  !R/L_n=-dndrho_tor_norm*(drho_tor_norm/drho_tor)/n*R
  rlnex=-rlnex/a00/nnex*Lref
  rltex=-rltex/a00/ttex*Lref
  rltix=-rltix/a00/ttix*Lref

  !... get metric coefficients if available
  Call set_metric_coeffs(eq(1), coretransp(1),coreprof(1), geom_container, Lref)
  IF (mpi_rank.eq.0) &
       & WRITE(*,'(A,I3,A)') "Running GENE at ", nrho_transp, &
       & " radial position(s)"

  DO i=1,nrho_transp
     IF (mpi_rank.eq.0) &
          & WRITE(*,'(A,F12.6)') "Running GENE at rho_tor_norm = ", &
          & rho_tor_norm(i)

     !local normalization
     nref = nnex(i)
     Tref = ttex(i)
     mref = mD !for now deuterium; for later it's maybe better to use the 
               !the mass of the main ion species from the summary CPO
     !ions
     do n=0,n_spec-2
        write(spec(n)%name,'(A,I1.1)') 'ions',(n+1)
        spec(n)%dens = nnix(i)/nref
        spec(n)%temp = ttix(i)/Tref
        spec(n)%mass = mD/mref
        spec(n)%charge = 1.0 !zeffx(i) ? but then, dens has to be adapted
        spec(n)%omt = rltix(i)
        spec(n)%omn = rlnex(i)
     enddo

     !electrons (as last species)
     if (n_spec.gt.1) then
        n=n_spec-1
        spec(n)%name = 'electrons'
        spec(n)%dens = nnex(i)/nref
        spec(n)%temp = ttex(i)/Tref
        spec(n)%charge = -1
        spec(n)%mass = me/mref
        spec(n)%omt = rltex(i)
        spec(n)%omn = rlnex(i)
     endif

     !geometry parameters
     if (magn_geometry.eq.'eq_cpo') then
        call assign_local_metric(geom_container,i)
     else
        q0 = qqx(i)
        shat=shatx(i)
        trpeps = rho_tor_norm(i)*minor_r/major_R !r/R0
     end if
     
     !other local parameters
     beta = 2*mu_0*nref*kb*Tref/(Bref*Bref)
     coll = 2.3031E-5*Lref*(nref*1E-19)/(Tref*1E-3)**2*lcoul  !lcoul is 14(23.11.10)
     debye2 = 5.2936E-4 * (Bref)**2/(nref*1E-19) * itm_mp/mref

     !other local normalization parameters for conversion to SI units
     rho_ref = cc*SQRT(mref*kb*Tref)/(ee*Bref)
     c_ref = SQRT(kb*Tref/mref)
     chigb = c_ref*rho_ref**2/Lref
     Gamma_gb = nref*chigb/Lref
     Qgb = Gamma_gb * kb * Tref

!!\todo remove bad style: redefine reference quantities for consistent output in units namelist
     mref = mref/mD*2.0   !from kg to m_proton    
     nref = nref*1E-19    !from m^-3 to 1e19 m^-3
     Tref = Tref*1E-3     !from eV to keV
     !Bref is in T already
     !Lref is in m already

     call check_params(ITM_par_in_dir)
     
     nexc_sav = nexc
     if (.not.nonlinear) then
        nky0_sav = nky0
        nky0 = 1
        kymin_sav = kymin
        ql_diff_max = 0.0
        do j=1,nky0_sav
           kymin = kymin_sav*j
           write(ITM_file_extension,'(A,I3.3,A,I3.3)') '_',i,'_',j
           call rungene(MPI_COMM_WORLD, ITM_par_in_dir, ITM_file_extension, &
                &ITM_checkpoint_in)
           if (ql_diffusivity.gt.ql_diff_max) &
                &gene_fluxes = quasilin_flux
           deallocate(quasilin_flux)           
        enddo
        nky0=nky0_sav
        kymin=kymin_sav
     else
        write(ITM_file_extension,'(A,I3.3)') '_',i
        call rungene(MPI_COMM_WORLD, ITM_par_in_dir, ITM_file_extension, &
             &ITM_checkpoint_in)
        gene_fluxes = avgfluxes
     endif
     nexc=nexc_sav

     call fill_coretransp(coretransp(1),i,gene_fluxes)
  END DO

  !...  stamp time
  coretransp(1)%time=coreprof(1)%time 


  !... cleanup
  call finalize_geomtype(geom_container)

#if defined(WITHSLEPC)
  if (.not.slepc_restartable()) &
       call my_SlepcFinalize(.true.)
#endif 

  IF (mpi_rank.eq.0) print*, '... done'

!  call mpi_finalize()


End Subroutine ITMgene


