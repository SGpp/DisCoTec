#include "redef.h"
#include "intrinsic_sizes.h"
!>Contains the routines to set the geometric coefficients
!!\todo This module is not structured and has no comments, plus it has a lot of parameters
module geometry
  use par_mod
  use par_geom
  use file_io
  use profile_io
  use communications
  use chease_mod
  use circular_mod
  use tracer_mod
!  use x_derivatives
  use lagrange_interpolation, only: lag3deriv, lag3interp
#ifdef WITHFUTILS
  USE futils 
  USE hashtable
#endif
  USE miller_mod

  implicit none
  public:: geomtype, q_prof, dqdx_prof, C_y, C_xy, C_j, geom, geom_in

  !input parameters
  public:: parscale, q0, major_R, minor_r, trpeps, shat, amhd, rhostar, mag_prof,&
       &q_coeffs, magn_geometry, geomdir, geomfile, flux_label, x_def, flux_pos,&
       &Lref, Bref, edge_opt, norm_flux_projection, dpdx_pm, dpdx_term, ialpha


  !output parameters
  public:: dVdx, sqrtgxx_FS, area_fs, flux_geomfac, dpdx_pm_arr

  public:: initialize_geometry, finalize_geometry, write_geometry, mem_est_geometry,&
       & read_tracer_namelist, set_geometry_defaults, check_geometry

  private

  !public geomtype variable 
  TYPE(geomtype) :: geom, geom_in

  !geometry
  Real::    parscale=1.0

  !C coefficients used to define x,y coordinates
  real, dimension(:), allocatable :: C_y, C_xy, C_j
  !Cyq0_x0 accounts for definitions which do not have Cy=x0/q0, only
  !necessary for correct kx-shifts in local simulations
  real:: Cyq0_x0

  !q profile
  Real, Dimension(:), Allocatable :: q_prof, dqdx_prof
  Real, DIMENSION(0:5) :: q_coeffs  ! used to determine the q profile

  !pressure gradient profile
  Real, DIMENSION(:), Allocatable :: dpdx_pm_arr
  real :: dpdx_pm
  !define how the grad p term enters the drift velocity
  Character(len=13) :: dpdx_term = '' !on (or full_drift),curv_eq_gradB,gradB_eq_curv 

  !misc
  logical:: mag_prof=.false.
  real:: edge_opt=0.
  integer :: ialpha=0 !relevant for GIST

  !local variables
  real, dimension(:),allocatable :: eps_p,dqdx
  real, dimension(:),allocatable :: dpdx_pm_geom !introduced for bookkeeping 
                                                 !grad p from geometry even if dpdx_pm<>-2

  !auxiliary quantities
  Real, dimension(:),allocatable :: dVdx, sqrtgxx_FS, area_fs
  Real, dimension(:,:,:),allocatable :: flux_geomfac !!pure geometry factor 
                                                     !!appearing in transport fluxes
  logical :: write_pe, norm_flux_projection=.false.

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_geometry(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !geometric quantities (>=3D)
    mem_loc=14*SIZE_OF_REAL_MB*px0*pj0*lk0
    !local arrays in write_geometry; set_curvature and shift_metric also have local arrays,
    !but this is the largest number allocated simultaneously
    mem_loc=mem_loc+SIZE_OF_REAL_MB*(11*px0*pj0*nz0+2*px0)
    !complex arrays in geom
    mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*px0*lj0*lk0*(2*nzb+1)+SIZE_OF_COMPLEX_MB*px0*pj0*lz0
    !add also initialization of circular and chease + curvature 
    mem_est_geometry=mem_req_in+mem_loc

  End Function mem_est_geometry

  !>Set default values for geometry parameters
  subroutine set_geometry_defaults
    dpdx_pm = -2 !-1111.0
    norm_flux_projection = .false.
    geomdir = './'
    tracer_initialized=.false.
    call set_miller_defaults

    Bprof_coeffs(:) = 0.0
    ialpha = 0
  end subroutine set_geometry_defaults

  !>Check geometry parameters after reading input file
  subroutine check_geometry
    integer :: fileunit
    logical :: file_exists

    !check for geometry file
    if (((INDEX(magn_geometry,'s_alpha').ne.1).and.&
         &INDEX(magn_geometry,'miller').ne.1).and.&
         &(magn_geometry.ne.'circular').and.(INDEX(magn_geometry,'slab').ne.1)) then
       inquire(file=trim(geomdir)//'/'//trim(geomfile),exist=file_exists)
       if (file_exists) then
          if ((magn_geometry.eq.'tracer').or.(magn_geometry.eq.'gene')) then
             fileunit = -1
             if (trim(par_in_dir).ne.'skip_parfile') &
                  & call read_tracer_namelist(fileunit,.true.,.false.)
          endif
       else
          write(*,"(4A)") trim(geomdir),'/',trim(geomfile),' does not exist'
          stop
       endif
    endif

    if ((magn_geometry.eq.'circular').and.(amhd.ne.0.0)) then
       if (mype.eq.0) write(*,'(A)') "WARNING: amhd<>0 not supported by 'circular' geometry module"
       amhd = 0.0
    endif

    !in the present code version, pressure_term just monitors whether press needs to
    !be computed & added in prefactors.F90
    !however, in older revisions it used to be an input parameter which is why we need
    !some backward compatibility stuff here:
    
    !very old pressure_off parameter:
    if (.not.pressure_off) pressure_term=.true.

    if (pressure_term) dpdx_term = 'full_drift'

    if (trim(dpdx_term).eq.'on') dpdx_term = 'full_drift'
    if (trim(dpdx_term).eq.'off') dpdx_term = 'gradB_eq_curv'

    if ((trim(dpdx_term).ne.'full_drift').and.(trim(dpdx_term).ne.'curv_eq_gradB').and.&
         &(trim(dpdx_term).ne.'gradB_eq_curv')) then
       if (mype.eq.0)    write(*,'(A)') "WARNING: no dpdx_term treatment defined!!"
       if ((dpdx_pm.eq.0).or.(dpdx_pm.eq.-1111.0)) then
          if (mype.eq.0) write(*,'(A)') "         Setting dpdx_term = 'curv_eq_gradB' ... "
          dpdx_term = 'curv_eq_gradB'
       else
          if (bpar.or.(.not.bpar_off)) then
             if (mype.eq.0) write(*,'(A)') "         Setting dpdx_term = 'full_drift' ... "
             dpdx_term = 'full_drift'
          else
             if (mype.eq.0) write(*,'(A)') "         Setting dpdx_term = 'gradB_eq_curv' ... "
             dpdx_term = 'gradB_eq_curv'   
          endif
       endif
    endif

    if (trim(dpdx_term).eq.'curv_eq_gradB') dpdx_pm=0.0

    !switch on pressure term if dpdx_pm<>0
    pressure_term = ((dpdx_pm.ne.-1111.0).and.(dpdx_pm.ne.0)).or.pressure_term

    if (pressure_term.and.dpdx_pm.eq.-1111) then !set dpdx_pm defaults
       !switch to profile based dpdx for simple analytical geometries
       !and to MHD based dpdx for more sophisticated models/interfaces
       if ((INDEX(magn_geometry,'s_alpha').eq.1).or.&
            &(magn_geometry.eq.'circular').or.(INDEX(magn_geometry,'slab').eq.1)) then
          if (write_pe) write(*,"(A)") "evaluating dpdx_pm from grad n,T"
          dpdx_pm = -1
       else
          if (write_pe) write(*,"(A)") "evaluating dpdx_pm from geometry"
          dpdx_pm = -2
       endif
    endif

  end subroutine check_geometry

  !>Allocates and sets all geometry-related fields according to the chosen geometry.
  subroutine initialize_geometry
    implicit none
    integer :: ispec,unit
    real, dimension(1) :: dummy
    real :: p_spec, Lref_sav, Bref_sav

    Lref_sav = Lref
    Bref_sav = Bref

    write_pe = ((mype==0).and.(print_ini_msg))

    allocate(C_y(pi1gl:pi2gl),C_xy(pi1gl:pi2gl),C_j(pi1gl:pi2gl))
    allocate(q_prof(pi1gl:pi2gl),dqdx_prof(pi1gl:pi2gl),dpdx_pm_geom(pi1gl:pi2gl))
    allocate(dpdx_pm_arr(pi1gl:pi2gl))
    dpdx_pm_arr = 0.0

    call initialize_geomtype(geom,pmi1gl,pmi2gl,pg1gl,pg2gl,pmi1,pmi2,&
       lj1,lj2,pj1,pj2,lk1,lk2,nzb)
    
    allocate(eps_p(pi1gl:pi2gl),dqdx(pi1gl:pi2gl))
    if (.not.allocated(dVdx)) allocate(dVdx(pg1gl:pg2gl))
    if (.not.allocated(sqrtgxx_FS)) allocate(sqrtgxx_FS(pg1gl:pg2gl))
    if (.not.allocated(area_FS)) allocate(area_FS(pg1gl:pg2gl))
    allocate(flux_geomfac(pi1:pi2,pj1:pj2,lk1:lk2))

    magn_geo_geofile = ""
    C_xy=1.0; C_y=1.0; C_j=1.0; Cyq0_x0=1.0
    
    !if rhostar is explicitly set negative, calculate it automatically if 
    !tracer_efit/chease along with numerical input profiles are used
    !(could include more cases in the future)
    if (rhostar.lt.0.0) then
       if (magn_geometry.eq.'tracer_efit') then
          call get_efit_refvals(Bref,Lref,major_R,minor_r,geomdir,geomfile)
       else if (magn_geometry.eq.'chease') then
#ifdef WITHFUTILS
          call get_chease_refvals(Bref,Lref,major_R,minor_r,geomfile,geomdir)
#else
          stop 'Need to compile with FUTILS to use Chease geometry!'
#endif
       elseif (magn_geometry.eq.'tracer'.or.magn_geometry.eq.'gist') then
          unit=-1
          call read_tracer_namelist(unit,.true.,.false.)
       else
          !otherwise, Bref and Lref need to be set explicitly
       endif
       call calc_rhostar(Bref,Lref,minor_r,rhostar)
    endif
    !map lx_a to lx
    if (.not.x_local.and.lx.eq.0.0.and.lx_a.ne.0.) lx=lx_a/rhostar

    !cannot be placed in check_geometry as before since beta might be -1 there
    if ((INDEX(magn_geometry,'miller').eq.1).or.(INDEX(magn_geometry,'s_alpha').eq.1)) then
       ! set alpha_mhd according to beta and gradients if not set by hand
       if (amhd.lt.0.0) then
          if (beta.lt.0.0) then
             if ((nref.lt.0).or.(Tref.lt.0).or.(Bref.lt.0)) stop &
                  & 'Cannot compute amhd and beta without given nref, Tref, and Bref'
             beta = 403.0E-5*nref*Tref/(Bref*Bref)
          endif
          amhd = q0**2*major_R*beta*sum((spec%omn+spec%omt) * spec%temp * spec%dens)
          if (n_spec.eq.1) then
             if (write_pe) WRITE(*,'(A)') &
                  &'assuming same pressure gradient for adiabatic species in amhd computation'
             amhd = amhd * 2.0
          endif
          if (write_pe) WRITE(*,'(A, F12.6)') 'used n/T gradients to compute amhd=',amhd
       endif
    endif

    if ((INDEX(magn_geometry,'s_alpha').eq.1).or.&
         &(INDEX(magn_geometry,'slab').eq.1).or.(magn_geometry.eq.'circular')) then 
       !Analytical magnetic geometries:
       !all necessary parameters (e.g. q0, shat, n_pol) are already specified;
       !the coordinates are thus initialized now (if possible) since some of 
       !them (zval, for instance) are explicitly used

       if (x_local.or.lilo) then
          !need to evaluate q_profile at x0 here for correct q0, shat
          dummy = x0
          call set_q_profile(dummy)
          if (mag_prof) then
             C_y = x0*minor_r/q0
          else
             if (trpeps.gt.0) C_y = trpeps*major_R/q0
          endif
          !So far, the local magnetic geometry models do not explicitly use
          !C_y which would affect the y coordinate --> not true anymore due to
          !adapt_ly = T option; 
          !hence, C_y has been replaced by x0*minor_r/q0 here
          
          call set_y_coordinate_vars(n0_global, q0, &
               &C_y, rhostar*minor_r, lilo, write_pe)
       endif
       !the x, z coordinates depend on geometry coefficients which are
       !already known and (only in the local code) possibly on ky
       call set_z_coordinate_vars(n_pol)
       call set_x_coordinate_vars(nexc, n_pol, rad_bc_type, &
            &shat, rhostar, lilo, adapt_lx, mag_prof, &
            &only_zonal, write_pe,Cyq0_x0)
       if (.not.(x_local.or.lilo)) call set_q_profile(xval_a)

       select case (magn_geometry)
       case('s_alpha_B')
          !s_alpha_B is the same as s_alpha but does not assume that a
          !Bfield factor is equal to 1 on the Bfield derivative, i.e.
          !there's no additional division by Bfield in set_curvature.
          !This model could be the one being implemented in GYRO
          !(important for beta scans) while 's_alpha' seems to be the 
          !traditional but less logic(?) choice
          if(.not.mag_prof) then
             call set_s_alpha(geom)
          else
             stop 'the s_alpha model does not include magnetic x profiles!'
          endif          
       case('s_alpha')
          if(.not.mag_prof) then
             call set_s_alpha(geom)
          else
             stop 'the s_alpha model does not include magnetic x profiles!'
          endif
       case('slab')
          call set_slab(geom,.false.)
       case('slab_curv')   
          call set_slab(geom,.true.)
       case('circular') !flux coordinate x=r
          call init_circ(n_pol,q_prof,eps_p,dqdx,major_R,minor_R)
          call get_circ(geom,C_y,C_xy,q0)
       case default
          stop 'no valid choice for magn_geometry'
       end select

       if (.not.(x_local.or.lilo).or.mag_prof) then
          !So far, the local magnetic geometry models do not explicitly use
          !C_y which would affect the y coordinate
          call set_y_coordinate_vars(n0_global, q0, &
               &C_y, rhostar*minor_r, lilo, write_pe)
       endif
       dpdx_pm_geom = amhd / (q0**2 * major_R)
    else
       !Numerical equilibria:
       !Typically, those interfaces set/redefine variables as, e.g., shat, q0, etc.
       !Hence, they must be called before the coordinate initialization
       
       select case (magn_geometry)
       case('chease')
          p_spec=0
          do ispec=0,n_spec-1
             p_spec=p_spec+spec(ispec)%temp*spec(ispec)%dens
          end do
          
          call init_chease(x_def,x0,flux_pos,p_spec,pressure_term,rhostar,&
               &print_ini_msg,geomfile,geomdir,force_x,force_z)
#ifdef WITHFUTILS
          call get_chease(geom,q0,trpeps,shat,beta,dpdx_pm_geom,lx,x0,&
               &C_y,C_xy,q_prof,dqdx_prof,Lref,Bref,major_R,minor_r,nx0)      
#endif
          !normalize Jacobian to L_perp
          !       jacobian=q0*jacobian
       case('tracer')
          call read_tracer
          call set_z_coordinate_vars(n_pol)
       case('gene') !read GENE geometry output format (similar to tracer interface)
          if (edge_opt.gt.0) stop 'edge_opt option not implemented yet for this interface'
          call read_tracer
          call set_z_coordinate_vars(n_pol)
          Cyq0_x0 = C_y(pi1gl)*q0/x0
       case('tracer_efit')
          if (.not.tracer_initialized) call init_tracer(geomdir, geomfile, x0, &
               &lx, rhostar, minor_r, major_R, edge_opt, trpeps)
          call get_tracer(geom,C_y,C_xy,q_prof,dqdx_prof,dpdx_pm_geom,&
               &q0,shat,Lref,Bref,edge_opt)
          if (final_init) call finalize_tracer
          call set_z_coordinate_vars(n_pol)
       case('gist')
          call read_gist
          call set_z_coordinate_vars(n_pol)
       case('eq_cpo')
          call get_from_geom_in
          call set_z_coordinate_vars(n_pol)
       case default
          if (index(magn_geometry,'miller').eq.1) then
             call get_miller(geom,C_y,C_xy,dpdx_pm_geom,q_prof,dqdx_prof,Cyq0_x0,edge_opt)
             call set_z_coordinate_vars(n_pol)
          else
             stop 'no valid choice for magn_geometry'
          endif
       end select
       call set_y_coordinate_vars(n0_global, q0, &
            &C_y, rhostar*minor_r, lilo, write_pe)
 
       call set_x_coordinate_vars(nexc, n_pol, rad_bc_type, &
            &shat, rhostar, lilo, adapt_lx, mag_prof, &
            &only_zonal, write_pe,Cyq0_x0)
       
       if (.not.(x_local.or.lilo)) then
          if (magn_geometry.eq.'gene'.or.magn_geometry.eq.'tracer') then
             call lag3deriv(q_prof,xval_a,nx0,dqdx_prof,xval_a,nx0)
             dqdx_prof = dqdx_prof/minor_r
          endif
       endif
    endif
    
    if (y_local) C_j = C_y

    if (shifted_metric) call shift_metric
    
    call set_curvature

    call compute_aux_geom_quantities

    if (par_in_dir.eq.'skip_parfile') then
       if (((Lref_sav.gt.0).and.(abs((Lref_sav-Lref)/Lref).gt.0.01)).or.&
            &((Bref_sav.gt.0).and.(abs((Bref_sav-Bref)/Bref).gt.0.01))) then
          if (mype.eq.0) then
             write(*,'(A,2(F12.6))') 'old/new Lref: ', Lref_sav, Lref
             write(*,'(A,2(F12.6))') 'old/new Bref: ', Bref_sav, Bref
             write(*,'(A,F12.6,A)') '(trinity: phia_in =',Lref**2*pi*Bref,')'
          endif
          call my_barrier()
          stop 'ERROR: geometry reference value mismatch in interface (see above)'
       endif
    endif

    if (pressure_term) then
       if (dpdx_pm.eq.-2) then
          if (all(abs(dpdx_pm_geom).lt.epsilon(dpdx_pm))) then
             dpdx_pm_arr = 0
             dpdx_pm = 0
             pressure_term = .false.
          else
             dpdx_pm_arr = dpdx_pm_geom
          endif
          if (write_pe.and.x_local) write(*,"(A,F12.6)") &
               &'dpdx_pm taken from geometry: ', dpdx_pm_arr(pi1gl)
       else
          if (write_pe.and.x_local) write(*,"(A,F12.6)") &
               &'info: dpdx_pm from geometry: ', dpdx_pm_geom(pi1gl)
       endif
    endif

    deallocate(dpdx_pm_geom)

  end subroutine Initialize_geometry


  !>Reads in geometric coefficients from a TRACER output file.
  subroutine read_tracer
    
    integer :: fileunit, i,k, line
    logical:: read_nml
    real :: dummy
    real, dimension(1:16) :: dummy_arr
    real, dimension(pmi1gl:pmi2gl,pj1:pj2,lk1:lk2):: gxx, gxz, gyy, gyz, dBdx, dBdy
    real, dimension(pmi1gl:pmi2gl,0:nz0-1)::read_arr

    ! read in new parameters
    fileunit = -1
    read_nml = (par_in_dir.ne.'skip_parfile')
    call read_tracer_namelist(fileunit,read_nml,.true.)

    if (xy_local) then
       do line=0,nz0-1
          if (my_pez.eq.(line/lk0)) then
             k=lk1+mod(line,lk0)
             read(fileunit,*) dummy_arr !avoid temporary variables by using this dummy array
             gxx(:,:,k) = dummy_arr(1)
             geom%gij(:,:,k) = dummy_arr(2)
             gxz(:,:,k) = dummy_arr(3)
             gyy(:,:,k) = dummy_arr(4)
             gyz(:,:,k) = dummy_arr(5)
             geom%gzz(:,:,k) = dummy_arr(6)
             geom%Bfield(:,:,k) = dummy_arr(7)
             dBdx(:,:,k) = dummy_arr(8)
             dBdy(:,:,k) = dummy_arr(9)
             geom%dBdz(:,:,k) = dummy_arr(10)
             geom%jacobian(:,:,k) = dummy_arr(11)
             geom%R(pmi1gl,k) = dummy_arr(12)
             geom%PHI(pmi1gl,k) = dummy_arr(13)
             geom%Z(pmi1gl,k) = dummy_arr(14)
             geom%dxdR(pmi1gl,k) = dummy_arr(15)
             geom%dxdZ(pmi1gl,k) = dummy_arr(16)
          else
             read(fileunit,*) dummy
          endif
       enddo

       if ((x0.ne.-1).and.(all(C_y.eq.1.0))) then
          C_y=x0/q0
          C_xy=1.
          minor_r=1.
       endif
    else
       call read_in_1d("q",q_prof(:),px0,fileunit)
       call read_in_2d("gxx",read_arr,pmx0,nz0,fileunit)
       gxx(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)
       call read_in_2d("gxy",read_arr,pmx0,nz0,fileunit)
       geom%gij(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)       
       call read_in_2d("gxz",read_arr,pmx0,nz0,fileunit)
       gxz(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)  
       call read_in_2d("gyy",read_arr,pmx0,nz0,fileunit)
       gyy(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)  
       call read_in_2d("gyz",read_arr,pmx0,nz0,fileunit)
       gyz(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)         
       call read_in_2d("gzz",read_arr,pmx0,nz0,fileunit)
       geom%gzz(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)
       call read_in_2d("Bfield",read_arr,pmx0,nz0,fileunit)
       geom%Bfield(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)
       call read_in_2d("dBdx",read_arr,pmx0,nz0,fileunit)
       dBdx(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)       
       call read_in_2d("dBdy",read_arr,pmx0,nz0,fileunit)
       dBdy(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)         
       call read_in_2d("dBdz",read_arr,pmx0,nz0,fileunit)
       geom%dBdz(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)         
       call read_in_2d("jacobian",read_arr,pmx0,nz0,fileunit)
       geom%jacobian(:,pj1,lk1:lk2) = read_arr(:,lk1:lk2)
       call read_in_1d("C_y",C_y,px0,fileunit)
       call read_in_1d("C_xy",C_xy,px0,fileunit)
       call read_in_2d("geo_R",read_arr,pmx0,nz0,fileunit)
       geom%R = read_arr(:,lk1:lk2)
       call read_in_2d("geo_Z",read_arr,pmx0,nz0,fileunit)
       geom%Z = read_arr(:,lk1:lk2)
       call read_in_2d("geo_c1",read_arr,pmx0,nz0,fileunit)
       geom%dxdR = read_arr(:,lk1:lk2)
       call read_in_2d("geo_c2",read_arr,pmx0,nz0,fileunit)
       geom%dxdZ = read_arr(:,lk1:lk2)
    endif

    close(fileunit)

    if (yx_order) then
       geom%gii=gyy
       geom%giz=gyz
       geom%gjj=gxx
       geom%gjz=gxz
       geom%dBdj=dBdx
       geom%dBdi=dBdy
    else
       geom%gii=gxx
       geom%giz=gxz
       geom%gjj=gyy
       geom%gjz=gyz
       geom%dBdi=dBdx
       geom%dBdj=dBdy
    end if

    if (x_local.or.lilo) then
       do i=pi1gl,pi2gl
          xval_a(i)=x0
          q_prof(i)=q0 
          dqdx(i)=1/C_y(i)*shat
       end do
    endif

    geom%R_hat=geom%R/Lref
    geom%Z_hat=geom%Z/Lref

    if (write_pe) write(*,"(a)") "Geometric coefficients have been read."

  end subroutine read_tracer

  !> Subroutine to just read the TRACER file namelist 
  !! Variables (q0, shat, etc.) will be taken by GENE
  !! if assign_vars is true
  !! This subroutine is required for some interfaces
  subroutine read_tracer_namelist(fileunit,assign_vars,keep_file_open)
    integer, intent(inout) :: fileunit
    logical, intent(in) :: assign_vars,keep_file_open
    integer :: gridpoints, nalpha0
    real :: my_dpdx=-100., q0_sav, shat_sav, n_pol_sav, beta_sav
    real :: trpeps_sav, minor_r_sav
    real:: Tref, nref, my_currdens_par=-100., Cy, Cxy !Lref, Bref
    logical:: file_exists
    character(len=FILENAME_MAX) :: magn_geometry_sav
    real :: s0

    !NOTE: my_dpdx is the negative pressure gradient normalized to 
    !      the reference magnetic pressure and called dpdx (w/o _pm)
    !      for historic reasons

    namelist /parameters/ &
         & q0, gridpoints, shat, n_pol, my_dpdx, my_currdens_par, &
         & beta, Tref, nref, Bref, Lref, Cy, Cxy, magn_geometry, n0_global, &
         & s0, nalpha0, trpeps, minor_r

    file_exists = .false.
    if (fileunit.lt.0) then
       Inquire(file=trim(geomdir)//'/'//trim(geomfile),exist=file_exists)
       if (.not.file_exists) then
          WRITE(*,"(4A)") trim(geomdir),'/',trim(geomfile),' does not exist'
          stop
       endif
       call get_unit_nr(fileunit)
       OPEN(fileunit,file=trim(geomdir)//'/'//trim(geomfile))       
    endif
    
    !save current GENE parameter settings
    q0_sav    = q0
    shat_sav  = shat
    minor_r_sav = minor_r
    trpeps_sav = trpeps
    n_pol_sav = n_pol
    beta_sav = beta
    magn_geometry_sav = magn_geometry

    ! set defaults
    Cy = -1.0
    Cxy = -1.0
    s0 = -1.0
    nalpha0 = -1
    my_dpdx = -100.0
    my_currdens_par = -100.0

    ! read in new parameters
    read(fileunit,nml=parameters)
    IF (gridpoints.NE.nz0) THEN
       WRITE(*,"(A,I5,A,I5,A)") "Your geometry file uses ",gridpoints,&
            & " gridpoints, while the par file contains nz0=",nz0,&
            & ". These two numbers must match."
       STOP
    END IF

    IF ( (.not. y_local) .and. (nalpha0.GT.0).and.(nalpha0.NE.nky0) ) THEN
       WRITE(*,"(A,I5,A,I5,A)") "Your geometry file uses ",nalpha0,&
            & " points in y, while the par file contains nky0=",nky0,&
            & ".  These two numbers must match.  Do not pass go, do not collect 100 dollars."
       STOP
    END IF

    if (.not.assign_vars) then
       q0    = q0_sav
       shat  = shat_sav
       n_pol = n_pol_sav
       beta = beta_sav
       trpeps = trpeps_sav
       minor_r = minor_r_sav
    else
       if (my_dpdx.ne.-100.0) dpdx_pm_geom = my_dpdx
       if (my_currdens_par.ne.-100.0) then
          currdens_par=my_currdens_par
          if(xy_local) currdens_par_scal=currdens_par(pi1)
       end if
       if ((Cy.ne.-1.0).and.(xy_local)) C_y = Cy
       if ((Cy.eq.-1.0).and.(s0.ne.-1.0)) then
          x0 = s0**(0.5)
          C_y = x0 / q0
       endif
       if ((Cxy.ne.-1.0).and.(xy_local)) C_xy = Cxy
       magn_geo_geofile = magn_geometry
    endif    
    magn_geometry = magn_geometry_sav

    trpeps = trpeps_sav

    if (file_exists.AND.(.NOT.keep_file_open)) CLOSE(fileunit)

  end subroutine read_tracer_namelist


!############################################################################################
!> Read magnetic geometry from HDF5 or ASCII file generated by GIST
!!note that the GIST interface does not provide the individual dBdx/dBdy derivatives
!!or the gxz, gyz, gzz metric elements but directly the curvature terms as evaluated
!!in [P. Xanthopolous et al., Physics of Plasmas 16, 082303 (2009)]
!!!which are here called dBdj and dBdi to avoid additional arrays
  subroutine read_gist
#ifdef WITHFUTILS
    Use futils
#endif
    real, dimension(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxx, gxy, gyy, dBdx, dBdy
!    real, dimension(pi1gl:pi2gl,pj1:pj2,8) :: dumarr
    real, allocatable, dimension(:,:,:) :: dumarr
!    real, dimension(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxz, gyz
    integer :: nalpha, fileunit, i, k, line, ealpha, xind
    real :: dummy
    
!#ifdef WITHFUTILS    
!    integer :: hdf5_iout, in_nx0, in_nz0
!    logical :: file_exists
!    real, dimension(1:nx0) :: s0, r0

!    Inquire(file=trim(geomdir)//'/'//trim(geomfile)//'0001',exist=file_exists)
!    if (.not.file_exists) then
!       WRITE(*,"(4A)") trim(geomdir)//'/'//trim(geomfile)//'0001', 'does not exist'
!       stop       
!    endif

!    call openf(trim(geomdir)//'/'//trim(geomfile),hdf5_iout)

!    call getatt(hdf5_iout,'/data/var0d/','nx0',in_nx0)    
    
!    if ((in_nx0.gt.1).and.(xy_local)) then
!       WRITE(*,"(3A)") 'nx0 in ',trim(geomfile),' should be 1 for local runs'
!       stop       
!    endif

!    if (nx0.ne.in_nx0) then
!       WRITE(*,"(A,I4,A)") 'inconsistent nx0 (',in_nx0,') in ',trim(geomfile)
!       stop
!    endif
!    call getatt(hdf5_iout,'/data/var0d/','nz0',in_nz0)    
!    if (nz0.ne.in_nz0) then
!       WRITE(*,"(A,I4,A)") 'inconsistent nz0 (',in_nz0,') in ',trim(geomfile)
!       stop
!    endif
    
!    call getatt(hdf5_iout,'/data/var0d/','npol',n_pol)      

!    if (in_nx0.eq.1) then
!       CALL getatt(hdf5_iout, "/data/var0d", "q0", q0) 
!       CALL getatt(hdf5_iout, "/data/var0d", "shat", shat) 
!       CALL getarr(hdf5_iout, "/data/var1d/g11", gxx(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/g12", gxy(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/g13", gxz(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/g22", gyy(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/g23", gyz(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/g33", geom%gzz(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/Bfield", geom%Bfield(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/dbdx", dBdx(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/dbdy", dBdy(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/dbdz", geom%dBdz(:,pj1,:) )
!       CALL getarr(hdf5_iout, "/data/var1d/jac", geom%jacobian(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var1d/r", geom%R ) 
!       CALL getarr(hdf5_iout, "/data/var1d/z", geom%Z )
!       CALL getarr(hdf5_iout, "/data/var1d/phi", geom%PHI )
!       CALL getarr(hdf5_iout, "/data/var1d/c1", geom%dxdR )
!       CALL getarr(hdf5_iout, "/data/var1d/c2", geom%dxdZ)        
!    else
!       CALL getarr(hdf5_iout, "/data/grid/s0", s0 )
!       CALL getarr(hdf5_iout, "/data/var1d/r0", r0 )
!       CALL getarr(hdf5_iout, "/data/var1d/q0", q_prof ) 
!       CALL getarr(hdf5_iout, "/data/var1d/qprime", dqdx )
!       CALL getarr(hdf5_iout, "/data/var2d/g11", gxx(:,pj1,:) )
!       CALL getarr(hdf5_iout, "/data/var2d/g12", gxy(:,pj1,:) ) 
!       IF (isdataset(hdf5_iout,"/data/var2d/g13")) &
!            &CALL getarr(hdf5_iout, "/data/var2d/g13", gxz(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var2d/g22", gyy(:,pj1,:) ) 
!       IF (isdataset(hdf5_iout,"/data/var2d/g23")) &
!            &CALL getarr(hdf5_iout, "/data/var2d/g23", gyz(:,pj1,:) )
!       IF (isdataset(hdf5_iout, "/data/var2d/g33")) &
!            &CALL getarr(hdf5_iout, "/data/var2d/g33", geom%gzz(:,pj1,:) )        
!       CALL getarr(hdf5_iout, "/data/var2d/Bfield", geom%Bfield(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var2d/dbdx", dBdx(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var2d/dbdy", dBdy(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var2d/dbdz", geom%dBdz(:,pj1,:) ) 
!       CALL getarr(hdf5_iout, "/data/var2d/jac", geom%jacobian(:,pj1,:) ) 
!!       CALL getarr(hdf5_iout, "/data/var2d/r", dum2 )
!!       CALL getarr(hdf5_iout, "/data/var2d/z", dum2 )
!!       CALL getarr(hdf5_iout, "/data/var2d/phi", dum2 )
!!       CALL getarr(hdf5_iout, "/data/var2d/c1", dum2 )
!!       CALL getarr(hdf5_iout, "/data/var2d/c2", dum2 )
!       
!       !set lx, calculate dqdr

!    endif

!    C_y = (s0**(0.5)) / q0

!    call closef(hdf5_iout)

!    if (write_pe) write(*,"(a)") "Geometric coefficients have been read."
! #else
    

    if(lilo) then
        allocate(dumarr(0:0,pj1:pj2,8))
    else
        allocate(dumarr(pi1gl:pi2gl,pj1:pj2,8))
    end if
!    real, dimension(pi1gl:pi2gl,pj1:pj2,8) :: dumarr

    ! read in new parameters first and open file
    fileunit=-1
    call read_tracer_namelist(fileunit,.true.,.true.)

    ! Read in the geometric coefficients
    if ( y_local ) then
       ealpha = ialpha
       if (ialpha.EQ.-1) then
          ealpha = ni0-1
          gxx = 0.0
          gxy = 0.0
          gyy = 0.0
          geom%Bfield = 0.0
          geom%jacobian = 0.0
          dBdx = 0.0
          dBdy = 0.0
          geom%dBdz = 0.0
       endif
       do nalpha=0, ealpha
          do line=0,nz0-1
             if (my_pez.eq.(line/lk0)) then
                k=lk1+mod(line,lk0)
                
                read(fileunit,*) dumarr
                if(lilo) then
                   if (ialpha.ge.0) then
                      do xind = pi1gl, pi2gl
                         gxx(xind,:,k) = dumarr(0,:,1)
                         gxy(xind,:,k) = dumarr(0,:,2)
                         gyy(xind,:,k) = dumarr(0,:,3)
                         geom%Bfield(xind,:,k) = dumarr(0,:,4)
                         geom%jacobian(xind,:,k) = dumarr(0,:,5)
                         dBdx(xind,:,k) = dumarr(0,:,6)
                         dBdy(xind,:,k) = dumarr(0,:,7)
                         geom%dBdz(xind,:,k) = dumarr(0,:,8)
                      end do
                   else
                      do xind = pi1gl, pi2gl
                         gxx(xind,:,k) = gxx(xind,:,k) + dumarr(0,:,1)/real(ni0)
                         gxy(xind,:,k) = gxy(xind,:,k) + dumarr(0,:,2)/real(ni0)
                         gyy(xind,:,k) = gyy(xind,:,k) + dumarr(0,:,3)/real(ni0)
                         geom%Bfield(xind,:,k) = geom%Bfield(xind,:,k)+dumarr(0,:,4)/real(ni0)
                         geom%jacobian(xind,:,k) = geom%jacobian(xind,:,k)+dumarr(0,:,5)/real(ni0)
                         dBdx(xind,:,k) = dBdx(xind,:,k) + dumarr(0,:,6)/real(ni0)
                         dBdy(xind,:,k) = dBdy(xind,:,k) + dumarr(0,:,7)/real(ni0)
                         geom%dBdz(xind,:,k) = geom%dBdz(xind,:,k) + dumarr(0,:,8)/real(ni0)
                      end do
                   end if
                else
                   if (ialpha.ge.0) then
                         gxx(:,:,k) = dumarr(:,:,1)
                      gxy(:,:,k) = dumarr(:,:,2)
                      gyy(:,:,k) = dumarr(:,:,3)
                      geom%Bfield(:,:,k) = dumarr(:,:,4)
                      geom%jacobian(:,:,k) = dumarr(:,:,5)
                      dBdx(:,:,k) = dumarr(:,:,6)
                      dBdy(:,:,k) = dumarr(:,:,7)
                      geom%dBdz(:,:,k) = dumarr(:,:,8)
                   else
                      gxx(:,:,k) = gxx(:,:,k) + dumarr(:,:,1)/real(ni0)
                      gxy(:,:,k) = gxy(:,:,k) + dumarr(:,:,2)/real(ni0)
                      gyy(:,:,k) = gyy(:,:,k) + dumarr(:,:,3)/real(ni0)
                      geom%Bfield(:,:,k) = geom%Bfield(:,:,k)+dumarr(:,:,4)/real(ni0)
                      geom%jacobian(:,:,k) = geom%jacobian(:,:,k)+dumarr(:,:,5)/real(ni0)
                      dBdx(:,:,k) = dBdx(:,:,k) + dumarr(:,:,6)/real(ni0)
                      dBdy(:,:,k) = dBdy(:,:,k) + dumarr(:,:,7)/real(ni0)
                      geom%dBdz(:,:,k) = geom%dBdz(:,:,k) + dumarr(:,:,8)/real(ni0)
                   endif
                end if
             else
                read(fileunit,*) dummy
             endif
          enddo
       enddo
    else
       if (mype==0.and.print_ini_msg) write(*,"(a)") "Loading y-global GIST geometry information"
       
       do nalpha=0,nky0-1
          do line=0,nz0-1
             if (my_pez.eq.(line/lk0)) then
                k=lk1+mod(line,lk0)
                read(fileunit,*) gxx(nalpha,pj1,k), gxy(nalpha,pj1,k), gyy(nalpha,pj1,k), &
                     geom%Bfield(nalpha,pj1,k), geom%jacobian(nalpha,pj1,k), dBdx(nalpha,pj1,k), &
                     dBdy(nalpha,pj1,k), geom%dBdz(nalpha,pj1,k)
             else
                read(fileunit,*) dummy
             endif
          enddo
       enddo
    endif
    close(fileunit)

    if (write_pe) write(*,"(a)") "Geometric coefficients have been read."
! #endif

    if (yx_order) then
       geom%gii=gyy
       geom%giz=0.
       geom%gjj=gxx
       geom%gjz=0.
       geom%dBdj=dBdx
       geom%dBdi=dBdy
    else
       geom%gii=gxx
       geom%giz=0.
       geom%gjj=gyy
       geom%gjz=0.
       geom%dBdi=dBdx
       geom%dBdj=dBdy
    end if
    geom%gij=gxy

    if (x_local.or.lilo) then
       do i=pi1gl,pi2gl
          xval_a(i)=x0
          q_prof(i)=q0 
          dqdx(i)=1/C_y(i)*shat
       end do
    endif

  end subroutine read_gist



  !>Sets geometric quantities to the values given by the analytic s-alpha model.  
  subroutine set_s_alpha(geom)
    TYPE(geomtype) :: geom
    real:: arg
    integer:: k
    real, dimension(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxx, gxy, gxz, gyy, gyz, dBdx, dBdy

    ! initialise geometric coefficients etc.
    do k=lk1,lk2
       arg    = zval(k) * shat - amhd * sin(zval(k))
       gxx(:,:,k) = 1.0
       gxy(:,:,k) = arg
       gyy(:,:,k) = 1.0+arg**2                    !dropped one term

       geom%Bfield(:,:,k) = 1.0/( 1.0+trpeps*cos(zval(k)) )
       dBdx(:,:,k)=  -cos(zval(k))*geom%Bfield(:,:,k)**2/major_R
       geom%dBdz(:,:,k) = trpeps*sin(zval(k))*geom%Bfield(:,:,k)**2

       geom%jacobian(:,:,k) = q0/geom%Bfield(:,:,k)*major_R  !not the real jacobian for historic reasons

       geom%dxdR(:,k)=cos(zval(k))
       geom%dxdZ(:,k)=sin(zval(k))
       geom%R_hat(:,k)=major_R*(1.+trpeps*cos(zval(k)))
       geom%Z_hat(:,k)=major_R*(1.+trpeps*sin(zval(k)))
    end do


    gxz = 0.
    dBdy= 0.
    if (trpeps.ne.0.) then       
       gyz = 1/(trpeps*major_R)
       geom%gzz = (1/(trpeps*major_R))**2
    else !if trpeps.eq.0 the parallel coordinate is not really defined, but for test purposes, 
       !we want finite K_y
       if (write_pe) write(*,"(a)") &
         'trpeps is 0, THE COORDINATE SYSTEM MAKES NO SENSE! Continuing with wrong but finite curvature.'
       gyz = 0.
       geom%gzz = 0.
    endif

    C_y=trpeps*major_R/q0    ! used to compute n0_global

    if (yx_order) then
       geom%gii=gyy
       geom%giz=gyz
       geom%gjj=gxx
       geom%gjz=gxz
       geom%dBdj=dBdx
       geom%dBdi=dBdy
    else
       geom%gii=gxx
       geom%giz=gxz
       geom%gjj=gyy
       geom%gjz=gyz
       geom%dBdi=dBdx
       geom%dBdj=dBdy
    end if
    geom%gij=gxy

  end subroutine set_s_alpha

  !>Sets geometric quantities to a sheared slab geometry
  subroutine set_slab(geom,with_curv)
    TYPE(geomtype) :: geom
    logical, intent(in) :: with_curv
    real:: arg
    integer:: i,j,k
    
    ! initialise geometric coefficients etc.
    do k=lk1,lk2
       arg = zval(k) * shat
       if(yx_order) then
          geom%gjj(:,:,k) = 1.0
          geom%gii(:,:,k) = 1.0+arg**2
       else
          geom%gii(:,:,k) = 1.0
          geom%gjj(:,:,k) = 1.0+arg**2
       end if
       geom%gij(:,:,k) = arg
       geom%giz(:,:,k) = 0.
       geom%gjz(:,:,k) = 0.
       geom%gzz(:,:,k) = 0.

       geom%Bfield(:,:,k) = 1.0

       geom%dBdi(:,:,k) = 0.
       geom%dBdj(:,:,k) = 0.
       geom%dBdz(:,:,k) = 0.

       geom%jacobian(:,:,k) = parscale
       geom%dxdR=1.0
       geom%dxdZ=0.0
       geom%R_hat=1.0
       geom%Z_hat=0.0

       if (with_curv) then
          if(yx_order) then
             do j=pj1,pj2
                do i=pi1,pi2
                   geom%K_j(i,j,k)=-sin(zval(k))
                   geom%K_i(i,j,k)=-(cos(zval(k))+arg*sin(zval(k)))
                end do
             end do
          else
             do j=pj1,pj2
                do i=pi1,pi2
                   geom%K_i(i,j,k)=-sin(zval(k))
                   geom%K_j(i,j,k)=-(cos(zval(k))+arg*sin(zval(k)))
                end do
             end do
          endif
          geom%K_i = geom%K_i*q0/parscale  !replaced s-a major_R by parscale/q0
          geom%K_j = geom%K_j*q0/parscale
       endif

    end do

  end subroutine set_slab

  !>Writes geometric information to the output directory.
  Subroutine write_geometry
#ifdef WITHFUTILS
    use futils
#endif
    integer:: fileunit, k, ierr
    real, dimension(pmi1gl:pmi2gl,pj1:pj2,0:nz0-1)::&
         ggxx, ggxy, ggyy, gBfield
    real, dimension(pmi1gl:pmi2gl,pj1:pj2,0:nz0-1)::&
         ggxz, ggyz, ggzz, gdBdx, gdBdy, gdBdz, gjacobian
    real, dimension(pmi1gl:pmi2gl,0:nz0-1)::write_arr
    real, dimension(pmi1gl:pmi2gl,0:nz0-1):: gl_R, gl_z, gl_phi, &
         gl_dxdR, gl_dxdZ
    character(len=FILENAME_MAX):: gfilename, hfilename
#ifdef WITHFUTILS
    integer :: fidgeometry_h5
    character(len=100) :: dset_name_geomcoeff
#endif
    
    !Warning: Once we really have y-dep. metrics the gathering in y has to be done as well!
    if(yx_order) then
       call mpi_gather(geom%gjj(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxx(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gij(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gjz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gii(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggyy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%giz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggyz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%dBdj(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            gdBdx(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%dBdi(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            gdBdy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    else
       call mpi_gather(geom%gii(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxx(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gij(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%giz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggxz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gjj(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggyy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%gjz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            ggyz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%dBdi(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            gdBdx(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
       call mpi_gather(geom%dBdj(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
            gdBdy(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    end if
    
    call mpi_gather(geom%gzz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
         ggzz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%Bfield(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
         gBfield(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%dBdz(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
         gdBdz(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%jacobian(pmi1gl,pj1,lk1),pmx0*pj0*lk0, MPI_REAL_TYPE,&
         gjacobian(pmi1gl,pj1,0),pmx0*pj0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)

    call mpi_gather(geom%R(pmi1gl,lk1),pmx0*lk0, MPI_REAL_TYPE,&
         gl_R(pmi1gl,0),pmx0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%z(pmi1gl,lk1),pmx0*lk0, MPI_REAL_TYPE,&
         gl_z(pmi1gl,0),pmx0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%phi(pmi1gl,lk1),pmx0*lk0, MPI_REAL_TYPE,&
         gl_phi(pmi1gl,0),pmx0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%dxdR(pmi1gl,lk1),pmx0*lk0, MPI_REAL_TYPE,&
         gl_dxdR(pmi1gl,0),pmx0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)
    call mpi_gather(geom%dxdZ(pmi1gl,lk1),pmx0*lk0, MPI_REAL_TYPE,&
         gl_dxdZ(pmi1gl,0),pmx0*lk0,MPI_REAL_TYPE,0, mpi_comm_z, ierr)

    if (magn_geometry.eq.'tracer') then
       gfilename=trim(geomfile)
    else
       gfilename=trim(magn_geometry)
    endif

    if ((index(gfilename,'.dat').eq.0).or.(trim(file_extension).ne.'.dat')) &
         & gfilename=trim(gfilename)//trim(file_extension)
    if (index(gfilename,'.h5').eq.0) hfilename=trim(gfilename)//'.h5'

    if (mype.eq.0) then
       if(write_std) then
          call get_unit_nr(fileunit)
          open(fileunit,file=trim(diagdir)//'/'//trim(gfilename),Action='write')

          write(fileunit,"(A)")    "&parameters"
          write(fileunit,"(A,I5)") "gridpoints = ",nz0 
          write(fileunit,"(A,ES24.16)") "q0 = ", q0          
          write(fileunit,"(A,ES24.16)") "shat = ", shat
          write(fileunit,"(A,ES24.16)") "s0   = ", x0**2.0
          write(fileunit,"(A,ES24.16)") "minor_r= ", minor_r
          write(fileunit,"(A,ES24.16)") "trpeps = ", trpeps
          write(fileunit,"(A,ES24.16)") "beta = ", beta
          write(fileunit,"(A,ES24.16)") "Lref = ", Lref
          write(fileunit,"(A,ES24.16)") "Bref = ", Bref
          if (n_pol.ne.1) write(fileunit,"(A,I3)") "n_pol  = ", n_pol
          if ((dpdx_pm_arr(pi1).ne.0.).and.(pi0.eq.1)) &
               & write(fileunit,"(A,ES24.16)") "my_dpdx = ", dpdx_pm_arr(pi1)
          if ((x_local).and.(C_y(pi1).ne.1.0))  &
               & write(fileunit,"(A,ES24.16)") "Cy = ", C_y(pi1)
          if ((x_local).and.(C_xy(pi1).ne.1.0)) &
               & write(fileunit,"(A,ES24.16)") "Cxy = ", C_xy(pi1)
          if (magn_geo_geofile.ne."") then
             write(fileunit,"(3A)") "magn_geometry = '", TRIM(magn_geo_geofile),"'"
          else
             write(fileunit,"(3A)") "magn_geometry = '", TRIM(magn_geometry),"'"
          endif
          Write(fileunit,"(A)")    "/"

          if (xy_local) then
             do k=0,nz0-1
                write(fileunit,"(16ES24.16)") ggxx(pi1,pj1,k), ggxy(pi1,pj1,k), &
                     &ggxz(pi1,pj1,k), ggyy(pi1,pj1,k), ggyz(pi1,pj1,k), ggzz(pi1,pj1,k),&
                     &gBfield(pi1,pj1,k), gdBdx(pi1,pj1,k), gdBdy(pi1,pj1,k), gdBdz(pi1,pj1,k),&
                     &gjacobian(pi1,pj1,k), gl_R(pi1,k), gl_phi(pi1,k),gl_z(pi1,k),&
                     &gl_dxdR(pi1,k), gl_dxdZ(pi1,k)
             end do
          else
             call write_out_1d("q",q_prof(:),px0,fileunit)
             write_arr = ggxx(:,pj1,:)
             call write_out_2d("gxx",write_arr,pmx0,nz0,fileunit)
             write_arr = ggxy(:,pj1,:)
             call write_out_2d("gxy",write_arr,pmx0,nz0,fileunit)
             write_arr = ggxz(:,pj1,:)
             call write_out_2d("gxz",write_arr,pmx0,nz0,fileunit)
             write_arr = ggyy(:,pj1,:)
             call write_out_2d("gyy",write_arr,pmx0,nz0,fileunit)
             write_arr = ggyz(:,pj1,:)
             call write_out_2d("gyz",write_arr,pmx0,nz0,fileunit)
             write_arr = ggzz(:,pj1,:)
             call write_out_2d("gzz",write_arr,pmx0,nz0,fileunit)
             write_arr = gBfield(:,pj1,:)
             call write_out_2d("Bfield",write_arr,pmx0,nz0,fileunit)
             write_arr = gdBdx(:,pj1,:)
             call write_out_2d("dBdx",write_arr,pmx0,nz0,fileunit)
             write_arr = gdBdy(:,pj1,:)
             call write_out_2d("dBdy",write_arr,pmx0,nz0,fileunit)
             write_arr = gdBdz(:,pj1,:)
             call write_out_2d("dBdz",write_arr,pmx0,nz0,fileunit)
             write_arr = gjacobian(:,pj1,:)
             call write_out_2d("jacobian",write_arr,pmx0,nz0,fileunit)
             call write_out_1d("C_y",C_y(:),px0,fileunit)
             call write_out_1d("C_xy",C_xy(:),px0,fileunit)
             call write_out_2d("geo_R",gl_R(:,:),pmx0,nz0,fileunit)
             call write_out_2d("geo_Z",gl_z(:,:),pmx0,nz0,fileunit)
             call write_out_2d("geo_c1",gl_dxdR(:,:),pmx0,nz0,fileunit)
             call write_out_2d("geo_c2",gl_dxdZ(:,:),pmx0,nz0,fileunit)
          endif
          close(fileunit)
       endif
       
#ifdef WITHFUTILS
       if(write_h5) then
          if(xy_local) then
             call creatf(trim(diagdir)//'/'//trim(hfilename), &
                  fidgeometry_h5, "geometry related coefficients", 'd')
             call creatg(fidgeometry_h5, '/metric')
             call creatg(fidgeometry_h5, '/Bfield_terms')
             call creatg(fidgeometry_h5, '/shape')
             call creatg(fidgeometry_h5, '/parameters')
             
             write(dset_name_geomcoeff, "(A)") "/metric/g^xx"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxx(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^xy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxy(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^xz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxz(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^yy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggyy(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^yz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggyz(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^zz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggzz(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/C_y"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, C_y)
             write(dset_name_geomcoeff, "(A)") "/metric/C_xy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, C_xy)

             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/Bfield"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gBfield(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdx"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdx(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdy(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdz(pi1,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/Jacobian"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gjacobian(pi1,pj1,:))
                
             write(dset_name_geomcoeff, "(A)") "/shape/R"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_R(pi1,:))
             write(dset_name_geomcoeff, "(A)") "/shape/Z"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_z(pi1,:))
             write(dset_name_geomcoeff, "(A)") "/shape/phi"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_phi(pi1,:))
             write(dset_name_geomcoeff, "(A)") "/shape/dxdR"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_dxdR(pi1,:))
             write(dset_name_geomcoeff, "(A)") "/shape/dxdZ"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_dxdZ(pi1,:))

             write(dset_name_geomcoeff, "(A)") "/parameters/gridpoints"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ nz0 /))
             write(dset_name_geomcoeff, "(A)") "/parameters/shat"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ shat /))
             write(dset_name_geomcoeff, "(A)") "/parameters/trpeps"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ trpeps /))
             write(dset_name_geomcoeff, "(A)") "/parameters/beta"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ beta /))
             write(dset_name_geomcoeff, "(A)") "/parameters/q0"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ q0 /))

          else

             call creatf(trim(diagdir)//'/'//trim(hfilename), &
                  fidgeometry_h5, "geometric related coefficients", 'd')
             call creatg(fidgeometry_h5, '/profile')
             call creatg(fidgeometry_h5, '/metric')
             call creatg(fidgeometry_h5, '/Bfield_terms')
             call creatg(fidgeometry_h5, '/shape')
             call creatg(fidgeometry_h5, '/parameters')
  
             write(dset_name_geomcoeff, "(A)") "/profile/q_prof"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, q_prof)
             
             write(dset_name_geomcoeff, "(A)") "/metric/g^xx"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxx(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^xy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxy(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^xz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggxz(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^yy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggyy(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^yz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggyz(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/g^zz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, ggzz(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/metric/C_y"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, C_y)
             write(dset_name_geomcoeff, "(A)") "/metric/C_xy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, C_xy)

             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/Bfield"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gBfield(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdx"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdx(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdy"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdy(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/dBdz"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gdBdz(:,pj1,:))
             write(dset_name_geomcoeff, "(A)") "/Bfield_terms/Jacobian"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gjacobian(:,pj1,:))

             write(dset_name_geomcoeff, "(A)") "/shape/R"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_R(:,:))
             write(dset_name_geomcoeff, "(A)") "/shape/Z"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_z(:,:))
             write(dset_name_geomcoeff, "(A)") "/shape/dxdR"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_dxdR(:,:))
             write(dset_name_geomcoeff, "(A)") "/shape/dxdZ"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, gl_dxdZ(:,:))

             write(dset_name_geomcoeff, "(A)") "/parameters/gridpoints"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ nz0 /))
             write(dset_name_geomcoeff, "(A)") "/parameters/shat"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ shat /))
             write(dset_name_geomcoeff, "(A)") "/parameters/beta"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ beta /))
             write(dset_name_geomcoeff, "(A)") "/parameters/q0"
             call putarr(fidgeometry_h5, dset_name_geomcoeff, (/ q0 /))
          end if
          call closef(fidgeometry_h5)
       end if
#endif
    endif
    
  End Subroutine write_geometry
  

  !>Calculates the curvature from the metrics.
  subroutine set_curvature
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2):: ga1, ga2, ga3
    integer:: k,j,i

    if (trim(magn_geometry).eq.'slab_curv') return !already set
    
    if ((trim(magn_geometry).eq.'gist').or.(magn_geo_geofile.eq.'gist').or.trim(magn_geometry).eq.'miller_b') then
       !note that the GIST interface does not provide the individual dBdx/dBdy derivatives
       !or the gxz, gyz, gzz metric elements but directly the curvature terms as evaluated
       !in [P. Xanthopolous et al., Physics of Plasmas 16, 082303 (2009)]
       !which are here called dBdj and dBdi to avoid additional arrays
       do k=lk1,lk2
          do j=pj1,pj2
             do i=pi1,pi2
                geom%K_i(i,j,k)=1/C_xy(i)*geom%dBdj(i,j,k)
                geom%K_j(i,j,k)=1/C_xy(i)*geom%dBdi(i,j,k)
             end do
          end do
       end do
    else
       do k=lk1,lk2
          do j=pj1,pj2
             do i=pi1,pi2
                ga1(i,j,k)=geom%gii(i,j,k)*geom%gjj(i,j,k)-geom%gij(i,j,k)*geom%gij(i,j,k)
                ga2(i,j,k)=geom%gii(i,j,k)*geom%gjz(i,j,k)-geom%gij(i,j,k)*geom%giz(i,j,k)
                ga3(i,j,k)=geom%gij(i,j,k)*geom%gjz(i,j,k)-geom%gjj(i,j,k)*geom%giz(i,j,k)
                geom%K_i(i,j,k)=1/C_xy(i)*(geom%dBdj(i,j,k)+ga2(i,j,k)/ga1(i,j,k)*geom%dBdz(i,j,k))
                geom%K_j(i,j,k)=1/C_xy(i)*(geom%dBdi(i,j,k)-ga3(i,j,k)/ga1(i,j,k)*geom%dBdz(i,j,k))
             end do
          end do
       end do
    endif

    !fix the sign
    if (yx_order) then
       geom%K_j=-geom%K_j
    else
       geom%K_i=-geom%K_i
    end if

    if ((magn_geometry.eq.'s_alpha').or.(magn_geo_geofile.eq.'s_alpha')) then
       geom%K_i=geom%K_i/geom%Bfield(pi1:pi2,:,:)
       geom%K_j=geom%K_j/geom%Bfield(pi1:pi2,:,:)
    endif
    
  end subroutine set_curvature
  

  subroutine shift_metric
    real, dimension(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxy, gyy, gyz, dBdx
    
    !only works with yx_order=.false. at the moment
    !swap variables
    gxy=geom%gij
    gyy=geom%gjj
    gyz=geom%gjz
    dBdx=geom%dBdi
    geom%shift(:,:,lk1:lk2)=gxy/geom%gii

    !shift metric coefficients
    geom%gij=0.
    geom%gjj=gyy-gxy**2/geom%gii
    geom%gjz=gyz-geom%giz*gxy/geom%gii
    geom%dBdi=dBdx+geom%shift(:,:,lk1:lk2)*geom%dBdj

  end subroutine shift_metric


  !>Sets the safety factor profile according to the coefficients given by the user.
  subroutine set_q_profile(xval_a)
    implicit none
    real, dimension(pi1gl:pi2gl) :: xval_a
    integer :: o, i
    
    if (mag_prof) then
       if (.not.allocated(in_qprof)) then
          if (q_coeffs(0).eq.-1111.0) then !read from file
             call read_qprof_from_file
          else
             !analytical (polynomial) q profile
             do i=pi1gl,pi2gl
                q_prof(i) = q_coeffs(0)
                dqdx(i)   = 0.0
                do o=1,SIZE(q_coeffs)-1
                   q_prof(i)=q_prof(i)+q_coeffs(o)*xval_a(i)**o
                   dqdx(i)=dqdx(i)+o*q_coeffs(o)*xval_a(i)**(o-1)
                enddo
                if (x_local.or.lilo) then
                   if (magn_geometry.eq.'circular') then
                      q0 = q_prof(i)
                      trpeps = xval_a(i)*minor_r/major_R
                      shat = xval_a(i)/q0*dqdx(i)
                   else
                      stop "stop: mag_prof not allowed for non-circular local runs"
                   endif
                endif
             end do
          endif
       else
          q_prof = in_qprof
          call lag3deriv(q_prof,xval,nx0,dqdx,xval,nx0)
          dqdx = dqdx/rhostar !change normalization to minor_r
       endif

       eps_p=xval_a*minor_r/major_R
       dqdx_prof=dqdx/minor_r
    else
       q_prof = q0
       if (x_local.or.lilo) THEN
          if (trpeps.ne.0.0) then
             eps_p=trpeps
             dqdx=q0/(trpeps*major_R)*minor_r*shat !! dq/dx is normalized to a
          elseif (INDEX(magn_geometry,'slab').ne.1) then
             stop "trpeps must be given!"
          endif
       else
          eps_p=x0*minor_r/major_R
          dqdx=q0/(major_R*eps_p)*shat*minor_r
          ! dqdx(i)=minor_r*q0/(trpeps*major_R)*shat/sqrt(1-trpeps**2) !for qbar profiles?
          trpeps = x0*minor_r/major_R
       endif
       dqdx_prof=dqdx
             
    endif
  end subroutine set_q_profile
  
  
  !>Creates an output file containing the safety factor profile.
  subroutine write_q_profile
    INTEGER :: i, fileunit
    LOGICAL :: op
    
    ! find free fileunit
    fileunit=100
    DO 
       INQUIRE(unit=fileunit,opened=op)
       IF (.NOT.op) EXIT
       fileunit = fileunit+1
    END DO
    OPEN(fileunit,file=trim(diagdir)//'/q_profile.dat')
    ! write file header
    write(fileunit,"(A)") "#     x          xrhostar         q              dqdx"     
    Do i=pi1gl,pi2gl
       write(fileunit,"(4ES16.6)")  xval(i), xval(i)/rhostar, q_prof(i), dqdx(i)
    end Do
    close(fileunit)
  end subroutine write_q_profile


  !>Initialize additional geometry related variables
  subroutine compute_aux_geom_quantities
    integer :: j,k,ierr
    real, dimension(pg1gl:pg2gl) :: tmparr_loc

    !------------------------------------------------------------------
    ! Calculate average Jacobian
    if (y_local) then
       tmparr_loc = 0.0
       do k=lk1,lk2
          do j=pj1,pj2
             tmparr_loc = tmparr_loc + geom%jacobian(:,j,k) 
          end do
       end do
       tmparr_loc = tmparr_loc/real(pj0*nz0)
    else
       tmparr_loc = sum(geom%jacobian)/real(nky0*pj0*nz0)
    endif

    if (pj0.eq.1) then
       Call mpi_allreduce(tmparr_loc, geom%avg_jaco_yz, pmg0,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_z, ierr)
    else
       Call mpi_allreduce(tmparr_loc, geom%avg_jaco_yz, pmg0,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_yz, ierr)
    endif

    geom%avg_jaco = SUM(geom%avg_jaco_yz)/real(pmg0)


    !------------------------------------------------------------------
    !dVdx = Sum_(y,z) jacobian(y,z) ly/ny0 lz/nz0
    !     = ly*n0 * 2*pi*n_pol/nz0 * sum_z jacobian(z)
    !     = ly*n0 * 2*pi*n_pol * geom%avg_jaco_yz
    !in unit of Lref**2 using ly*n0= 2*pi*Cy
    dVdx = (2.0*pi)**2*C_y(pg1gl:pg2gl)*n_pol * geom%avg_jaco_yz(pg1gl:pg2gl)
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    !Area_FS = n0*Sum_(y,z) sqrt(geom%gxx(x,z)) jacobian(x,z) * ly/ny0 lz/nz0
    !in units of L_ref**2
    if (y_local) then
       tmparr_loc = 0.0
       do k=lk1,lk2
          do j=pj1,pj2
             if (yx_order) then
                tmparr_loc = tmparr_loc + sqrt(geom%gjj(:,j,k))*geom%jacobian(:,j,k) 
             else
                tmparr_loc = tmparr_loc + sqrt(geom%gii(:,j,k))*geom%jacobian(:,j,k) 
             endif
          end do
       end do
       tmparr_loc = tmparr_loc/real(pj0*nz0)
    else
       tmparr_loc = sum(sqrt(geom%gjj)*geom%jacobian)/real(nky0*pj0*nz0)
    endif

    if (pj0.eq.1) then
       Call mpi_allreduce(tmparr_loc, area_FS, pmg0,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_z, ierr)
    else
       Call mpi_allreduce(tmparr_loc, area_FS, pmg0,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_yz, ierr)
    endif

    area_FS = (2.0*pi)**2*C_y(pg1gl:pg2gl)*n_pol*area_FS
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    !sqrtgxx_FS = 1/dVdx * Sum_(y,z) sqrt(geom%gxx(x,z)) jacobian(x,z) * &
    ! & n0*ly/ny0 lz/nz0
    !           = n0*ly*lz/dVdx * Sum(z) sqrt(geom%gxx)*jacobian(x)/(nz0)
    sqrtgxx_FS = area_FS/dVdx
    !------------------------------------------------------------------

    if (magn_geometry=='tracer'.and.major_R==1..and.Lref.ne.0) then
       major_R=sum(geom%R(pi1,:)*geom%jacobian(pi1,pj1,:))
       call my_sum_to_all_real_0d(major_R,mpi_comm_z)
       major_R=major_R/(nz0*geom%avg_jaco)/Lref
    endif

    do k=lk1,lk2
       do j=pj1,pj2
          flux_geomfac(:,j,k) = 1.0/(C_xy(pi1:pi2))
          if (norm_flux_projection) then
             if (yx_order) then 
                flux_geomfac(:,j,k) = flux_geomfac(:,j,k)/&
                     sqrt(geom%gjj(pi1:pi2,j,k))
             else
                flux_geomfac(:,j,k) = flux_geomfac(:,j,k)/&
                     sqrt(geom%gii(pi1:pi2,j,k))
             endif
          endif
       enddo
    enddo
    
  end subroutine compute_aux_geom_quantities
  
  subroutine get_from_geom_in
    !geom_in is not distributed
    if (yx_order) then
       geom%gjj(pi1,pj1,lk1:lk2) = geom_in%gjj(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gij(pi1,pj1,lk1:lk2) = geom_in%gij(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gjz(pi1,pj1,lk1:lk2) = geom_in%gjz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gii(pi1,pj1,lk1:lk2) = geom_in%gii(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%giz(pi1,pj1,lk1:lk2) = geom_in%giz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gzz(pi1,pj1,lk1:lk2) = geom_in%gzz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%dBdj(pi1,pj1,lk1:lk2)= geom_in%dBdi(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%dBdi = 0.0    
    else
       geom%gii(pi1,pj1,lk1:lk2) = geom_in%gii(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gij(pi1,pj1,lk1:lk2) = geom_in%gij(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%giz(pi1,pj1,lk1:lk2) = geom_in%giz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gjj(pi1,pj1,lk1:lk2) = geom_in%gjj(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gjz(pi1,pj1,lk1:lk2) = geom_in%gjz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%gzz(pi1,pj1,lk1:lk2) = geom_in%gzz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%dBdi(pi1,pj1,lk1:lk2)= geom_in%dBdi(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
       geom%dBdj = 0.0    
    end if

    geom%jacobian(pi1,pj1,lk1:lk2) = geom_in%jacobian(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
    geom%Bfield(pi1,pj1,lk1:lk2) = geom_in%Bfield(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)
    geom%dBdz(pi1,pj1,lk1:lk2)= geom_in%dBdz(1,1,1+my_pez*lk0:1+(my_pez+1)*lk0)

    geom%R(pi1,lk1:lk2)    =    geom_in%R(1,1+my_pez*lk0:1+(my_pez+1)*lk0)
    geom%Z(pi1,lk1:lk2)    =    geom_in%Z(1,1+my_pez*lk0:1+(my_pez+1)*lk0)
    geom%dxdR(pi1,lk1:lk2) =    geom_in%dxdR(1,1+my_pez*lk0:1+(my_pez+1)*lk0)
    geom%dxdZ(pi1,lk1:lk2) =    geom_in%dxdZ(1,1+my_pez*lk0:1+(my_pez+1)*lk0)

  end subroutine get_from_geom_in


  !>Reads ASCII files containing the safety factor profile
  !!The format should be the following:
  !!1.) Two header lines beginning with the # symbol
  !!    where the first line should contain as many string entries
  !!    (without space characters) as data columns are present.
  !!    The second line is left for comments (e.g., time stamp)
  !!2.) The data should be arranged as follows:
  !!    1.column: x-axis using the same definition as GENE
  !!    2.column: arbitrary x-axis (e.g. for plotting) which is 
  !!              ignored
  !!    3.column: q profile
  SUBROUTINE read_qprof_from_file
    Character(Len=32) :: filename
    INTEGER :: ncols, nx0_file, fileunit, iostat, count
    REAL, dimension(:,:),allocatable :: data_arr
    REAL, dimension(:), allocatable :: q_file, xaxis_file
    Character(Len=128) :: header

    filename = 'profile_geom'

    call get_unit_nr(fileunit)
    count=10
    iostat=1
    do while (iostat.ne.0.and.count.gt.0)
       OPEN(fileunit,file=trim(filename),iostat=iostat)
       count=count-1
    enddo

    IF (iostat.ne.0) then
       WRITE(*,"(3A)") "ERROR reading q profile file ",trim(filename)
       STOP
    ENDIF

    nx0_file = get_nr_of_lines(fileunit)
    !skip header lines
    nx0_file = nx0_file - 2
    if (nx0_file.lt.1) then
       WRITE(*,"(3A)") "ERROR: found less than one radial grid point in ",trim(filename) 
       STOP
    endif

    rewind(fileunit)

    !read and ignore header lines which just contain general
    !information
    read(fileunit,"(A)",iostat=iostat) header
    !determine number of columns
    ncols = get_nr_of_substr(header)
    !ignore '#' symbol:
    ncols = ncols - 1

    if (ncols.lt.3) stop 'Less than 3 columns in profile file'

    allocate(data_arr(1:ncols,1:nx0_file))
    read(fileunit,"(A)",iostat=iostat) header
    read (fileunit,*,iostat=iostat) data_arr
    close(fileunit)

    if (iostat.ne.0) then
       WRITE(*,"(3A)") "Could not read profile from ",trim(filename)
       STOP
    ENDIF

    allocate(xaxis_file(1:nx0_file),q_file(1:nx0_file))

    xaxis_file = data_arr(1,:)
    q_file  = data_arr(3,:)
    deallocate(data_arr)

    call lag3interp(q_file,xaxis_file,nx0_file,&
         q_prof,xval_a,px0)
    call lag3deriv(q_file,xaxis_file,nx0_file,&
         dqdx,xval_a,px0)

    deallocate(xaxis_file, q_file)

  END SUBROUTINE read_qprof_from_file


 !>Deallocates geometry-related fields
 subroutine finalize_geometry
    implicit none

    deallocate(q_prof,dqdx_prof,eps_p,dqdx)
    deallocate(dpdx_pm_arr)
    deallocate(C_xy, C_y, C_j)
    deallocate(flux_geomfac)

    call finalize_geomtype(geom)

    if (par_in_dir.ne.'skip_parfile') &
         &deallocate(dVdx, area_FS, sqrtgxx_FS)

    !cleanup circular model, tracer, CHEASE??
    
  end subroutine finalize_geometry

end module geometry
