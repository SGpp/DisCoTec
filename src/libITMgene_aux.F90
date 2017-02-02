#include "redef.h"
!>Auxiliary module file for the ITM interface library
Module libITMgene_aux
  Use parameters_IO
  Use par_geom
  Use ITM_TYPES
  USE ITM_Constants
  USE Euitm_schemas

  implicit none

  public :: check_params, read_xml_parameters, &
       &nrho_transp, nion, turb_constructor, turb_destructor, &
       &fill_coretransp, set_metric_coeffs, &
       &mpi_rank, chigb, Qgb, Gamma_gb

  !  private

  INTEGER(ITM_I4), SAVE :: nrho_transp = 0, nion = 0

  REAL(r8) :: twopi=2.0_r8*itm_pi

  REAL(r8) :: mu_0 = itm_mu0
  REAL(r8) :: kb=itm_ev
  REAL(r8) :: ee=itm_qe
  REAL(r8) :: cc=1.0_r8, lcoul=14.0_r8

  REAL(r8) :: me=itm_me
  REAL(r8) :: md=itm_md
  
  REAL(R8) :: chigb,Qgb,Gamma_gb

  REAL(R8), DIMENSION(:), ALLOCATABLE :: rho_tor,rho_tor_norm, &
       nnex,ttex,nnix,ttix,qqx,zeffx,rlnex,rltex,rltix,shatx

  INTEGER :: mpi_rank

contains

  !-----------------------------------------------------------------------
  subroutine read_xml_parameters(code_parameters, return_status)

    !-----------------------------------------------------------------------
    ! calls the XML parser for the code parameters and assign the
    ! resulting values to the corresponding variables
    !-----------------------------------------------------------------------

    USE parameters_IO
    use euitm_schemas
    use euitm_xml_parser  

    implicit none

    type (type_param), intent(in) :: code_parameters
    integer(itm_i4), intent(out) :: return_status 

    type(tree) :: parameter_list
    type(element), pointer :: temp_pointer
    integer(itm_i4) :: i, nparm, n_vals
    character(len = 132) :: cname

    integer :: n

    !...  initialisations

    nparm = 0
    n_vals = 0
    return_status = 0      ! no error
    n = -1

    !some defaults for the ITM interface
    chptdir = '' !the xml lib has problems with an empty string as default value
    istep_nrg = 10
    istep_field = 100
    istep_mom = 500
    istep_energy = 200
    istep_g1 = 0
    istep_vsp = 1000
    istep_schpt = 1000
    
    x_local = .true.
    y_local = .true.
    ntimesteps = 1e9
    arakawa_zv = .true.
    collision_op = 'landau'
    magn_geometry = 'eq_cpo'
    hyp_z = 0.5
    hyp_v = 0.2

    !! transport definition !! CHECK
    norm_flux_projection = .false.


    !-- parse xml-string code_parameters%parameters using W3C XML schema in
    !   code_parameters%schema
    call euitm_xml_parse(code_parameters, nparm, parameter_list)


    !-- assign variables

    temp_pointer => parameter_list%first

    outer: do
       cname = char2str(temp_pointer%cname)   ! necessary for AIX

       select case (cname)
          !--   parameters overall
       case ("parameters")
          temp_pointer => temp_pointer%child
          cycle
          !------------------------------------------------------
       case ("parallelization")
          temp_pointer => temp_pointer%child
          cycle
          !--   individual parameters
       case ("n_procs_s")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_s)
       case ("n_procs_v")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_v)
       case ("n_procs_w")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_w)
       case ("n_procs_x")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_x)
       case ("n_procs_y")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_y)
       case ("n_procs_z")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_z)
       case ("n_procs_sim")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_procs_sim)
          !--------------------------------------------------------
       case ("box")
          temp_pointer => temp_pointer%child
          cycle
          !--   individual parameters
       case ("nx0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nx0)
       case ("nky0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nky0)
       case ("nz0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nz0)
       case ("nv0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nv0)
       case ("nw0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nw0)
       case ("n_spec")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, n_spec)
       case ("lx")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, lx)
       case ("kymin")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, kymin)
       case ("lv")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, lv)
       case ("lw")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, lw)
       case ("ky0_ind")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, ky0_ind)
       case ("kx_center")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, kx_center)
       case ("nexc")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nexc)
       case ("adapt_lx")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, adapt_lx)
       case ("x0")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, x0)
       case ("mu_grid_type")
          if (allocated(temp_pointer%cvalue)) &
               mu_grid_type = char2str(temp_pointer%cvalue)
       case ("n0_global")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n0_global)
          !--------------------------------------------------------
       case ("in_out")
          temp_pointer => temp_pointer%child
          cycle
       case ("diagdir")
          if (allocated(temp_pointer%cvalue)) &
               diagdir = char2str(temp_pointer%cvalue)
       case ("chptdir")
          if (allocated(temp_pointer%cvalue)) &
               chptdir = char2str(temp_pointer%cvalue)
       case ("read_checkpoint")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, read_checkpoint)
       case ("write_checkpoint")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, write_checkpoint)
       case ("istep_schpt")          
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_schpt)
       case ("istep_field")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_field)
       case ("istep_mom")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_mom)
       case ("istep_nrg")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_nrg)
       case ("istep_vsp")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_vsp)
       case ("istep_omega")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_omega)
       case ("istep_g1")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_g1)
       case ("istep_energy")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, istep_energy)
       !--------------------------------------------------------
       case ("general")
          temp_pointer => temp_pointer%child
          cycle      
       case ("nonlinear")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, nonlinear)
       case ("calc_dt")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, calc_dt)
       case ("dt_max")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, dt_max)
       case ("x_local")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, x_local)
!       case ("perf_vec")
          !scan_str2num is currently not available in the "old" libraries 
          !
          !        if (allocated(temp_pointer%cvalue)) &
          !             call scan_str2num(temp_pointer%cvalue, perf_vec, n_vals)
       case ("ntimesteps")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, ntimesteps)
       case ("timelim")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, timelim)
       case ("simtimelim")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, simtimelim)
       case ("collision_op")
          if (allocated(temp_pointer%cvalue)) &
               collision_op = char2str(temp_pointer%cvalue)
       case ("coll_cons_model")
          if (allocated(temp_pointer%cvalue)) &
               coll_cons_model = char2str(temp_pointer%cvalue)
       case ("avgflux_stime")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, avgflux_stime)
       case ("avgprof_stime")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, avgprof_stime)
       case ("omega_prec")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, omega_prec)
       case ("comp_type")
          if (allocated(temp_pointer%cvalue)) &
               comp_type=char2str(temp_pointer%cvalue)
       case ("which_ev")
          if (allocated(temp_pointer%cvalue)) &
               which_ev=char2str(temp_pointer%cvalue)
       case ("n_ev")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, n_ev)
       case ("arakawa_zv")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, arakawa_zv)
       case ("courant")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, courant)
       case ("hyp_x")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, hyp_x)
       case ("hyp_y")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, hyp_y)
       case ("hyp_z")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, hyp_z)
       case ("hyp_v")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, hyp_v)
          !--------------------------------------------------------
       case ("geometry")
          temp_pointer => temp_pointer%child
          cycle    
       case ("magn_geometry")
          if (allocated(temp_pointer%cvalue)) &
               magn_geometry = char2str(temp_pointer%cvalue)
       case ("shat")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, shat)
       case ("q0")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, q0)
       case ("major_R")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, major_R)
       case ("minor_r")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, minor_r)
       case ("trpeps")
          !read from CPO
          !            if (allocated(temp_pointer%cvalue)) &
          !                 call char2num(temp_pointer%cvalue, trpeps)
       case ("rhostar")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, rhostar)
       case ("flux_pos")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, flux_pos)
          !--------------------------------------------------------
       case ("species")
          !read from CPO
          !            if (n.eq.-1) then
          !               if (n_spec.lt.1) then
          !                  write (*,"(A)") "number of species has to be "//&
          !                       &"defined before first species namelist"
          !                  stop
          !               endif
          !               allocate(spec(0:n_spec-1))
          !            endif
          n = n + 1
          !            if (n.lt.n_spec) then
          !               temp_pointer => temp_pointer%child
          !               cycle    
          !            endif
       case ("name")
          if (allocated(temp_pointer%cvalue)) &
               spec(n)%name = char2str(temp_pointer%cvalue)
       case ("passive")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%passive)
       case ("omn")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%omn)
       case ("omt")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%omt)
       case ("mass")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%mass)
       case ("charge")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%charge)
       case ("dens")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%dens)
       case ("temp")
          if (allocated(temp_pointer%cvalue)) &
               call char2num(temp_pointer%cvalue, spec(n)%temp) 
       case default
          write(*, *) 'ERROR: invalid parameter', cname
          return_status = 1
          exit
       end select
       do
          if (associated(temp_pointer%sibling)) then
             temp_pointer => temp_pointer%sibling
             exit
          end if
          if (associated(temp_pointer%parent, parameter_list%first )) &
               exit outer
          if (associated(temp_pointer%parent)) then
             temp_pointer => temp_pointer%parent
          else
             write(*, *) 'ERROR: broken list.'
             return
          end if
       end do
    end do outer

    !-- destroy tree
    call destroy_xml_tree(parameter_list)

    auto_parall = ((n_procs_s.LT.1).OR.&
         &(n_procs_w.LT.1).OR.(n_procs_v.LT.1).OR.&
         &(n_procs_x.LT.1).OR.(n_procs_y.LT.1).OR.&
         &(n_procs_z.LT.1))

  end subroutine read_xml_parameters


  !--------------------------------------------------------------------------
  SUBROUTINE Turb_Constructor(coretransp, nrho0, nrho, nion)

    IMPLICIT NONE

    type (type_coretransp) :: coretransp
    integer(ITM_I4) :: nrho0, nrho, nion, nang0, nang

    !...  allocate transport structure

    ALLOCATE(coretransp%rho_tor(nrho0:nrho))
    ALLOCATE(coretransp%rho_tor_norm(nrho0:nrho))
    ALLOCATE(coretransp%ne_transp%flux(nrho0:nrho))
    ALLOCATE(coretransp%te_transp%flux(nrho0:nrho))
    ALLOCATE(coretransp%ni_transp%flux(nrho0:nrho,nion))
    ALLOCATE(coretransp%ti_transp%flux(nrho0:nrho,nion))
    ALLOCATE(coretransp%ne_transp%diff_eff(nrho0:nrho,3))
    ALLOCATE(coretransp%te_transp%diff_eff(nrho0:nrho))
    ALLOCATE(coretransp%ni_transp%diff_eff(nrho0:nrho,nion,3))
    ALLOCATE(coretransp%ti_transp%diff_eff(nrho0:nrho,nion))
    ALLOCATE(coretransp%ne_transp%vconv_eff(nrho0:nrho,3))
    ALLOCATE(coretransp%te_transp%vconv_eff(nrho0:nrho))
    ALLOCATE(coretransp%ni_transp%vconv_eff(nrho0:nrho,nion,3))
    ALLOCATE(coretransp%ti_transp%vconv_eff(nrho0:nrho,nion))

    coretransp%ne_transp%flux=0._R8
    coretransp%te_transp%flux=0._R8
    coretransp%ni_transp%flux=0._R8
    coretransp%ti_transp%flux=0._R8
    coretransp%ne_transp%diff_eff=0._R8
    coretransp%te_transp%diff_eff=0._R8
    coretransp%ni_transp%diff_eff=0._R8
    coretransp%ti_transp%diff_eff=0._R8
    coretransp%ne_transp%vconv_eff=0._R8
    coretransp%te_transp%vconv_eff=0._R8
    coretransp%ni_transp%vconv_eff=0._R8
    coretransp%ti_transp%vconv_eff=0._R8

  END SUBROUTINE Turb_Constructor

  !--------------------------------------------------------------------------
  SUBROUTINE Turb_Destructor(coretransp)
    !--------------------------------------------------------------------------

    IMPLICIT NONE

    TYPE (type_coretransp) :: coretransp

    !...  deallocate transport structure

    DEALLOCATE(coretransp%rho_tor)
    DEALLOCATE(coretransp%rho_tor_norm)
    DEALLOCATE(coretransp%ne_transp%flux)
    DEALLOCATE(coretransp%te_transp%flux)
    DEALLOCATE(coretransp%ni_transp%flux)
    DEALLOCATE(coretransp%ti_transp%flux)
    DEALLOCATE(coretransp%ne_transp%diff_eff)
    DEALLOCATE(coretransp%te_transp%diff_eff)
    DEALLOCATE(coretransp%ni_transp%diff_eff)
    DEALLOCATE(coretransp%ti_transp%diff_eff)
    DEALLOCATE(coretransp%ne_transp%vconv_eff)
    DEALLOCATE(coretransp%te_transp%vconv_eff)
    DEALLOCATE(coretransp%ni_transp%vconv_eff)
    DEALLOCATE(coretransp%ti_transp%vconv_eff)

  END SUBROUTINE Turb_Destructor

  !-----------------------------------------------------------------------------
  SUBROUTINE fill_coretransp(coretransp,i,gene_flux)
    IMPLICIT NONE
    
    type (type_coretransp) :: coretransp
    real, dimension(0:n_spec-1,2), intent(in) :: gene_flux

    integer :: i, n, pe, max, normspec, io, OUT

    !local variables
    REAL(R8) :: ff,gg,diff,chi, gxx_FS, lperp
    character(len=128):: results_file_extension

    write(results_file_extension,'(A,I3.3)') '_',i

    !fill coretransp CPO with results
    !... first the fluxes
    coretransp%ne_transp%flux(i) = gene_flux(n_spec-1,1)*Gamma_gb
    coretransp%te_transp%flux(i) = gene_flux(n_spec-1,2)*Qgb
    
    do n = 0, n_spec-2
       coretransp%ni_transp%flux(i,n+1) = gene_flux(n,1)*Gamma_gb
       coretransp%ti_transp%flux(i,n+1) = gene_flux(n,2)*Qgb
    enddo


    !... now the D's and V's
    lperp=1.0_R8/MAX( ABS(rltex(i)), ABS(rltix(i)), ABS(rlnex(i)))
    
    ff = coretransp%ne_transp%flux(i)/nnex(i)
    gg = coretransp%te_transp%flux(i)/(nnex(i)*kb*ttex(i))
    
    gxx_FS = sqrtgxx_FS(pi1gl)**2
    
    chi=ABS(gg)*lperp/gxx_FS
    diff=ABS(ff)*lperp/gxx_FS
    
    IF (diff < 0.2_R8*chi) diff=0.2_R8*chi
    chi = MAX(chi,0.2_R8*diff)
    
    coretransp%ne_transp%diff_eff(i,2) = diff
    coretransp%te_transp%diff_eff(i) = chi
    coretransp%ne_transp%vconv_eff(i,2) = (ff - diff*gxx_FS*rlnex(i))/gxx_FS
    coretransp%te_transp%vconv_eff(i) = (gg - chi*gxx_FS*rltex(i))/gxx_FS
    
    
    do n=1,n_spec-1
       ff = coretransp%ni_transp%flux(i,n)/nnix(i)
       gg = coretransp%ti_transp%flux(i,n)/(nnix(i)*kb*ttix(i))        
       chi = MAX(ABS(gg)*lperp,0.2_R8*diff)
       
       coretransp%ni_transp%diff_eff(i,n,2) = diff
       coretransp%ti_transp%diff_eff(i,n) = chi
       coretransp%ni_transp%vconv_eff(i,n,2) = (ff - diff*gxx_FS*rlnex(i))/gxx_FS
       coretransp%ti_transp%vconv_eff(i,n) = (gg - chi*gxx_FS*rltix(i))/gxx_FS
       
    enddo
    
    do pe = 0,0 !7
       if (mpi_rank.eq.pe) then
          do io=0,1
             if (io.eq.0) then
                OUT=6
             else
                call get_unit_nr(OUT)
                OPEN(OUT,file=TRIM(diagdir)//"/results"//trim(results_file_extension))
             endif

             write(OUT,'(4(A,(ES12.4,X)))') '#rho_tor_norm = ', coretransp%rho_tor_norm(i),&
                  &', Gamma_gb = ', Gamma_gb, ', Qgb = ', Qgb,&
                  &', <sqrt(gxx)>_FS = ',sqrtgxx_FS

             write(OUT,'(7(A12))') '#spec       ', 'nj_flux/gB', 'nj_flux', 'tj_flux/gB', 'tj_flux', 'nj_diff', 'tj_chi'
             write(OUT,'(A12,6(ES12.4,X))') 'electrons', coretransp%ne_transp%flux(i)/Gamma_gb,&
                  & coretransp%ne_transp%flux(i), coretransp%te_transp%flux(i)/Qgb,&
                  & coretransp%te_transp%flux(i), coretransp%ne_transp%diff_eff(i,2),&
                  & coretransp%te_transp%diff_eff(i)
             
             do n=1,n_spec-1
                write(OUT,'(A12,6(ES12.4,X))') spec(n-1)%name, coretransp%ni_transp%flux(i,n)/Gamma_gb,&
                     &coretransp%ni_transp%flux(i,1), coretransp%ti_transp%flux(i,n)/Qgb,&
                     &coretransp%ti_transp%flux(i,1), coretransp%ni_transp%diff_eff(i,n,2),&
                     &coretransp%ti_transp%diff_eff(i,n)
             enddo

             if (io.gt.0) CLOSE(OUT)
          enddo
          write(*,*)
       endif
    enddo

  END SUBROUTINE fill_coretransp


  !---------------------------------------------------------------------------------
  !------- GEOMETRY related stuff ----------
  !---------------------------------------------------------------------------------
  SUBROUTINE set_metric_coeffs(eq, coretransp,coreprof, geocont, Lref)
    TYPE (type_equilibrium) :: eq
    TYPE (type_coretransp) :: coretransp
    TYPE (type_coreprof) :: coreprof
    TYPE(geomtype) :: geocont
    real(r8), intent(in) :: Lref
    Real(R8),dimension(:,:),allocatable ::g11,g12,g13,g22,g23,g33
    Real(R8),dimension(:,:),allocatable ::B_2d,dBdi,dBdz,R_2d,Z_2d
    Real(R8),dimension(:,:),allocatable ::dRdx_2d,dZdx_2d,tmp_2d
    Real(R8),dimension(:),allocatable ::ageom,Rgeom,ageom_cpo,Rgeom_cpo
    Real(R8),dimension(:),allocatable ::dpsidrho_tor_norm,rho_tor_norm
    Real(R8),dimension(:),allocatable ::q,dqdrho_tor_norm,p,dpdrho_tor_norm
    Real(R8),dimension(:),allocatable ::rho_tor,shatx
    real(r8)::Bref,r0,p0,dqdpsi0,Cy,a00,Lperp

    REAL(R8), DIMENSION(1:nz0) :: theta
    !allocate conversion factors
    INTEGER(ITM_I4) :: nrho_eq, nrho_transp, nang_eq
    INTEGER(ITM_I4) :: i,j,k,kpi

    IF (trim(magn_geometry) .eq. 'eq_cpo') THEN
       IF (mpi_rank.eq.0) print*, "reading the equilibrium CPO"

       nrho_eq = SIZE(eq%profiles_1d%rho_vol)
       nang_eq = SIZE(eq%coord_sys%grid%dim2)
       nrho_transp = SIZE(coretransp%rho_tor_norm)

       !Allocation

       !initialize conversion factors

       !initialize temporary fields on angle=0..2pi grid
       ALLOCATE(g11(1:nrho_transp,1:nz0),g12(1:nrho_transp,1:nz0))
       ALLOCATE(g13(1:nrho_transp,1:nz0))
       ALLOCATE(g22(1:nrho_transp,1:nz0),g23(1:nrho_transp,1:nz0),g33(1:nrho_transp,1:nz0))
       ALLOCATE(R_2d(1:nrho_transp,1:nz0),Z_2d(1:nrho_transp,1:nz0))
       ALLOCATE(dRdx_2d(1:nrho_transp,1:nz0),dZdx_2d(1:nrho_transp,1:nz0))
       ALLOCATE(tmp_2d(1:nrho_transp,1:nang_eq))
       ALLOCATE(B_2d(1:nrho_transp,1:nz0))
       ALLOCATE(dBdi(1:nrho_transp,1:nz0),dBdz(1:nrho_transp,1:nz0))
       ALLOCATE(ageom(1:nrho_transp),Rgeom(1:nrho_transp))
       ALLOCATE(ageom_cpo(1:nrho_eq),Rgeom_cpo(1:nrho_eq))
       ALLOCATE(rho_tor_norm(1:nrho_transp))
       ALLOCATE(rho_tor(1:nrho_transp))
       ALLOCATE(shatx(1:nrho_transp))
       ALLOCATE(dpsidrho_tor_norm(1:nrho_transp))
       ALLOCATE(q(1:nrho_transp),dqdrho_tor_norm(1:nrho_transp))
       !IF (pressure_term) THEN
       ALLOCATE(p(1:nrho_transp),dpdrho_tor_norm(1:nrho_transp))
       !END IF

       !--- CPO bug check (to be erased in the future)
       i = nrho_eq/2
       if ((abs(eq%coord_sys%grid%dim2(nang_eq)-twopi-eq%coord_sys%grid%dim2(1)).le.&
            epsilon(eq%coord_sys%grid%dim2(1))).and.&
            &(abs(eq%coord_sys%g_11(i,1)-eq%coord_sys%g_11(i,nang_eq))&
            &.gt.epsilon(eq%coord_sys%g_11(i,1)))) then
          print*, 'WARNING: The equilibirum CPO probably uses a wrong theta axis'
          print*, '         Changing to coord_sys%grid%dim2 = 0, .., 2*pi - dz'
          do k=1,nang_eq
             eq%coord_sys%grid%dim2(k) = twopi/nang_eq*(k-1)
          enddo
       endif
       !---

       !read reference values from CPO: 
       Bref =coreprof%toroid_field%b0
       a00=eq%eqgeometry%a_minor

       rho_tor_norm = coretransp%rho_tor_norm
       rho_tor = coretransp%rho_tor

       !interpolate to target resolution
       !for the imp4 benchmark we assume that eq%rho_vol and rho_tor_norm is the same
       !TODO question: we use rho_vol and rho_tor_norm for interpolating metrics. maybe changed in future 

       !new grid (ranging here from [0,2*pi) for the interpolation 
       !          - later it will be shifted )
       do k=1,nz0
          theta(k) = twopi/nz0*(k-1)
       enddo

       call convert_metric_coeff(eq%coord_sys%g_11, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g11,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%g_12, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g12,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%g_13, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g13,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%g_22, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g22,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%g_23, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g23,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))         

       call convert_metric_coeff(eq%coord_sys%g_33, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            g33,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%position%R, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            R_2d,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       call convert_metric_coeff(eq%coord_sys%position%Z, &
            eq%profiles_1d%rho_vol,nrho_eq,&
            eq%coord_sys%grid%dim2,nang_eq,&
            Z_2d,rho_tor_norm,&
            nrho_transp,theta,SIZE(theta))

       !use rho_vol and rho_tor_norm as axes
       !dBdi = dBdrho_tor_norm
       call get_magnetic_field(eq%profiles_1d%F_dia,eq%coord_sys%position%R,&
            eq%coord_sys%g_11,eq%coord_sys%g_33,eq%profiles_1d%rho_vol,&
            nrho_eq,&
            eq%coord_sys%grid%dim2, size(eq%coord_sys%grid%dim2),&
            B_2d,dBdi,dBdz,&
            rho_tor_norm,nrho_transp,theta,size(theta) )

       ageom_cpo   = 0.5*(eq%coord_sys%position%R(:,1)- &
            &      eq%coord_sys%position%R(:,nang_eq/2) )
       Rgeom_cpo   = 0.5*(eq%coord_sys%position%R(:,1)+ &
            &      eq%coord_sys%position%R(:,nang_eq/2) )

       !bring ageom and Rgeom onto transp rho_tor grid
       call L3interp(ageom_cpo,eq%profiles_1d%rho_vol,nrho_eq,ageom,rho_tor_norm,nrho_transp)
       call L3interp(Rgeom_cpo,eq%profiles_1d%rho_vol,nrho_eq,Rgeom,rho_tor_norm,nrho_transp)
       CALL L3deriv ( eq%coord_sys%grid%dim1,&
            eq%profiles_1d%rho_vol, nrho_eq,&
            dpsidrho_tor_norm, rho_tor_norm, nrho_transp)
       CALL L3interp(eq%profiles_1D%q,eq%profiles_1D%rho_vol,nrho_eq,q,rho_tor_norm,nrho_transp)
       CALL L3deriv (eq%profiles_1D%q,eq%profiles_1D%rho_vol,nrho_eq,dqdrho_tor_norm,rho_tor_norm,nrho_transp)
       IF (pressure_term) THEN
          !get pressure profile and derivative on rho_tor_transp axis
          CALL L3interp(eq%profiles_1D%pressure,eq%profiles_1D%rho_vol,nrho_eq,p,rho_tor_norm,nrho_transp)
          CALL L3deriv (eq%profiles_1D%pressure,eq%profiles_1D%rho_vol,nrho_eq,dpdrho_tor_norm,rho_tor_norm,nrho_transp)
       END IF

       !for diagnostic purposes: dRdx, dZdx
       do k=1,nang_eq
          CALL L3deriv ( eq%coord_sys%position%R(:,k),&
               eq%profiles_1d%rho_vol, nrho_eq,&
               tmp_2d(:,k), rho_tor, nrho_transp)
       enddo
       do i=1,nrho_transp
          CALL L3interp(tmp_2d(i,:),eq%coord_sys%grid%dim2,nang_eq,&
               &dRdx_2d(i,:), theta, SIZE(theta))
       enddo
       !--
       do k=1,nang_eq
          CALL L3deriv ( eq%coord_sys%position%Z(:,k),&
               eq%profiles_1d%rho_vol, nrho_eq,&
               tmp_2d(:,k), rho_tor, nrho_transp)
       enddo
       do i=1,nrho_transp
          CALL L3interp(tmp_2d(i,:),eq%coord_sys%grid%dim2,nang_eq,&
               &dZdx_2d(i,:), theta, SIZE(theta))
       enddo       



       !metric is given as a function of poloidal flux coord and straight field line
       !angle -> transform to rho_tor and poloidal angle and interpolate
       !to simulation grid, the result is stored in 2D structures

       !coordinates in the eq. CPO: 
       !(poloidal flux label, field line label, straight field line angle)
       !this is the same as in chease
       !this part follows the cheas module

       !transform from rho_tor_norm to rho_tor whith a00 (not always correct!!)
       !\todo generalize!
       !dqdrho_tor_norm = dqdrho_tor_norm!/ageom_cpo(nrho_eq) 
       !x definition from chease: rho_tor only
       !\todo chease uses x= a * rho_tor_norm
       !!     here we want x= rho_tor_norm directly.
       !! should the ageom_cpo be taken here? I guess so.
       !         x_o_a = rho_tor_norm   ! equalize rho_tor_norm whith x/a
       !dxdpsi(:)=ageom_cpo(nrho_eq)*drho_tor_normdpsi(:)
       !         dxdpsi(:)=ageom_cpo(nrho_eq)/dpsidrho_tor_norm(:)


       !shift theta such that we later use the -pi to +pi-dz definition
       theta(nz0/2+1:nz0) = theta(nz0/2+1:nz0) - twopi

       ! loop over each flux surface
       Do i=1,nrho_transp

          !some reference quantities (from chease set xbox) 
          !NPSI0(on chease mesh) corresponds to index i of flux surface simulated i guess
          !r0=ageom(i)   !CHEASE like
!          r0=rho_tor_norm(i)*minor_r  !we do not use ageom, but minor_r to give x==rho_tor_norm a length
          r0 = rho_tor_norm(i)*a00
          q0=q(i)
          if (pressure_term) then 
             p0=p(i)
          end if
          dqdpsi0=dqdrho_tor_norm(i)/dpsidrho_tor_norm(i)
          !y def
          Cy=r0/q0
!          ! setting shear and trpeps (for local code, or for info in global code)
!          !shat     = Cy/dxdpsi*dqdpsi  
!          != Cy dqdrho_tor_norm drho_tor_normdx= r0/q0 dqdrho_tor_norm /ageom(NPSI)
!          !shatx = Cy*dqdrho_tor_norm/ageom_cpo(nrho_eq) !CHEASE like definition (we ues shat from cpo in libITMgene.F90)
!          shatx = Cy*dqdrho_tor_norm/minor_r !cpo definition like in libITMgene.F90)
          shatx = Cy*dqdrho_tor_norm/a00
          shat = shatx(i) !cpo definition like in libITMgene.F90)
 !         trpeps   = r0/major_R
          trpeps = r0/(major_R*Lref)

!          print*,"r0",r0
!          print*,"major_R", major_R
!          print*,"Lref", Lref
!          print*,"q",q(i)
!          print*,"shat",shatx(i)
!          print*,"trpeps", trpeps

          !*******************************
          ! Calculation of values needed in GENE, 
          ! note all output quantites are normalized
          ! we use x=rho_tor_norm z= straight field line angle, y=binormal coordinate
          !!todo normalization -> Lref and stuff. jacobian unclear (hkd)

          j=1 !pj1 not known during initialization, hence "1"
          DO k=1,nz0
             !indices are moved so that gab(i,1) is calculated at chi=-pi (check this...)
             kpi=MODULO(nz0/2+k-1,nz0)+1

             ! Calculation of the metric (normalized to Lnorm)
             geocont%gii(i,j,k) = g11(i,kpi)/dpsidrho_tor_norm(i)**2   
             geocont%gij(i,j,k) = Cy/dpsidrho_tor_norm(i)*(dqdpsi0 * theta(kpi)*g11(i,kpi)+ &
                  q(i)*g12(i,kpi)-g13(i,kpi))
             geocont%gjj(i,j,k) = Cy**2 * (dqdpsi0**2 * theta(kpi)**2 * g11(i,kpi)   + &
                  2*q(i)*dqdpsi0*theta(kpi)*g12(i,kpi)+q(i)**2*g22(i,kpi) &
                  -2*dqdpsi0*theta(kpi)*g13(i,kpi) - 2*q(i)*g23(i,kpi) + g33(i,kpi) )
             geocont%giz(i,j,k) = Lref / dpsidrho_tor_norm(i) * g12(i,kpi)
             geocont%gjz(i,j,k) = Lref * Cy * ( theta(kpi) * dqdpsi0* g12(i,kpi) + &
                  q(i) * g22(i,kpi) - g23(i,kpi))
             geocont%gzz(i,j,k) = Lref**2 * g22(i,kpi)
             geocont%Bfield(i,j,k)   = B_2d(i,kpi)/Bref
             geocont%dBdi(i,j,k) = Lref/Bref*dBdi(i,kpi)
             geocont%dBdz(i,j,k) = dBdz(i,kpi)/Bref
             ! jacobian is contravariant: J= (grad_psi x grad_chi . grad_phi)   (no ^-1) we need J_2d^-1
             !geocont%jacobian(i,j,k) = dpsidrho_tor_norm(i)/Cy*J_2d(i,kpi)/Lref
             !geocont%jacobian(i,j,k) = dpsidrho_tor_norm(i)/Cy*J_2d(i,kpi)**(-1)/Lref

             !self-consistent: J=1/sqrt(det(g^ij))
             geocont%jacobian(i,j,k)=sqrt(geocont%gii(i,j,k)*(geocont%gjj(i,j,k)*geocont%gzz(i,j,k)-&
                  &geocont%gjz(i,j,k)*geocont%gjz(i,j,k)) - &
                  &geocont%gij(i,j,k)*(geocont%gij(i,j,k)*geocont%gzz(i,j,k)-&
                  &geocont%giz(i,j,k)*geocont%gjz(i,j,k)) + &
                  &geocont%giz(i,j,k)*(geocont%gij(i,j,k)*geocont%gjz(i,j,k)-&
                  &geocont%giz(i,j,k)*geocont%gjj(i,j,k)))**(-1)

             !variables just for visualization (may not work right now)
             geocont%R(i,k) = R_2d(i,kpi)
             geocont%Z(i,k) = Z_2d(i,kpi)
             geocont%dxdR(i,k) = 1./dRdx_2d(i,kpi) !dxdR 
             geocont%dxdZ(i,k) = 1./dZdx_2d(i,kpi) !dxdz

          END DO ! k loop
       end Do  ! i loop

       DEALLOCATE(g11,g12,g13)
       DEALLOCATE(g22,g23,g33)
       DEALLOCATE(B_2d,R_2d,Z_2d)
       DEALLOCATE(dRdx_2d,dZdx_2d,tmp_2d)
       DEALLOCATE(dBdi,dBdz)
       DEALLOCATE(ageom,Rgeom)
       DEALLOCATE(ageom_cpo,Rgeom_cpo)
       DEALLOCATE(rho_tor_norm,rho_tor)
       DEALLOCATE(shatx)
       DEALLOCATE(dpsidrho_tor_norm)
       !IF (pressure_term) THEN
       DEALLOCATE(p,dpdrho_tor_norm)
       !END IF

    ELSE
       print*,"use circular geometry"
       magn_geometry = 'circular'
    ENDIF


  Contains

    Subroutine convert_metric_coeff(metric_in, axis1_in, naxis1_in, &
         axis2_in, naxis2_in, metric_out, axis1_out, naxis1_out, &
         axis2_out, naxis2_out)
      INTEGER(ITM_I4), INTENT(IN) :: naxis1_in, naxis2_in, naxis1_out,&
           naxis2_out
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_in) :: metric_in
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_out) :: metric_out
      REAL(R8), DIMENSION(1:naxis1_in) :: axis1_in
      REAL(R8), DIMENSION(1:naxis2_in) :: axis2_in
      REAL(R8), DIMENSION(1:naxis1_out) :: axis1_out
      REAL(R8), DIMENSION(1:naxis2_out) :: axis2_out          
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_in) :: temp_metric

      INTEGER :: i
      !step 1 radial interpolation
      Do i=1, naxis2_in
         CALL L3interp( metric_in(:,i), axis1_in, naxis1_in,&
              temp_metric(:,i), axis1_out, naxis1_out)
      End do

      !step 2 angular interpolation 
      Do i=1, naxis1_out
         CALL L3interp( temp_metric(i,:), axis2_in, naxis2_in,&
              metric_out(i,:), axis2_out, naxis2_out)
      End do

    End Subroutine Convert_metric_coeff

    Subroutine get_magnetic_field(F_dia_in,R_2d_in,g11_in,g33_in, axis1_in, naxis1_in, &
         axis2_in, naxis2_in, B_2d_out,dBdi_out,dBdz_out,axis1_out, naxis1_out, &
         axis2_out, naxis2_out)
      INTEGER(ITM_I4), INTENT(IN) :: naxis1_in, naxis2_in, naxis1_out,&
           naxis2_out
      REAL(R8), DIMENSION(1:naxis1_in) :: F_dia_in
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_in) :: R_2d_in
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_in) :: g11_in
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_in) :: g33_in
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_out) :: B_2d_out!on out,out grid
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_out) :: dBdi_out!on out,out grid
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_out) :: dBdz_out!on out,out grid
      REAL(R8), DIMENSION(1:naxis1_in) :: axis1_in
      REAL(R8), DIMENSION(1:naxis2_in) :: axis2_in
      REAL(R8), DIMENSION(1:naxis1_out) :: axis1_out
      REAL(R8), DIMENSION(1:naxis2_out) :: axis2_out          
      !temporary arrays
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_in) :: B_2d_cpo  !on in,in grid
      REAL(R8), DIMENSION(1:naxis1_out,1:naxis2_in) :: tmp_out_in !on out,in grid
      REAL(R8), DIMENSION(1:naxis1_in,1:naxis2_out) :: tmp_in_out !on in,out grid
      !REAL(R8) :: twopi=2.0_R8*itm_pi
      INTEGER(ITM_I4) :: i
      !magnetic field
      !F_dia = eq%profiles_1D%F_dia 
      !R_2d=eq%coord_sys%position%R
      !replace R_2d by g33 possibly
      do i=1,naxis2_in
         B_2d_cpo(:,i) = sqrt((F_dia_in*F_dia_in+g11_in(:,i)/twopi**2)/(R_2d_in(:,i)*R_2d_in(:,i)))
      end do

      do i=1,naxis2_in !CHI loop
         call L3interp(B_2d_cpo(:,i),axis1_in,naxis1_in,tmp_out_in(:,i),axis1_out,naxis1_out)
      end do

      do i=1,naxis1_out !x loop
         call L3interp(tmp_out_in(i,:),axis2_in,naxis2_in,B_2d_out(i,:),axis2_out,naxis2_out)
      end do
      !do the B derivative on in_grid and interpolate to out_grid in a second step
      do i=1,naxis2_in !CHI loop for x derivative
         call L3deriv(B_2d_cpo(:,i),axis1_in,naxis2_in,tmp_out_in(:,i),axis1_out,naxis1_out)
      end do
      do i=1,naxis1_out !Interpolate x derivative on chi out grid
         call L3interp(tmp_out_in(i,:),axis2_in,naxis2_in,dBdi_out(i,:),axis2_out,naxis2_out)
      end do
      do i=1,naxis1_in !PSI loop for CHI derivative
         call L3deriv(B_2d_cpo(i,:),axis2_in,naxis2_in,tmp_in_out(i,:),axis2_out,naxis2_out)
      end do
      do i=1,naxis2_out !Interpolate CHI derivative on x out grid
         call L3deriv(tmp_in_out(:,i),axis1_in,naxis1_in,dBdz_out(:,i),axis1_out,naxis1_out)
      end do
    End Subroutine get_magnetic_field

  END SUBROUTINE set_metric_coeffs

  !>Copies information from i'th radial entry 
  !!from geom_container to geom_in (see geometry.F90)
  SUBROUTINE assign_local_metric(geocont,i)
    TYPE(geomtype) :: geocont
    INTEGER :: i

    !assuming local code (px0=1) and nzb=2
    IF (.NOT.ALLOCATED(geom_in%gii)) &
         call initialize_geomtype(geom_in,1,1,1,1,1,1,&
         &1,1,1,1,1,nz0,2)
    
    geom_in%gii(1,1,:) = geocont%gii(i,1,:)
    geom_in%gij(1,1,:) = geocont%gij(i,1,:)
    geom_in%giz(1,1,:) = geocont%giz(i,1,:)
    geom_in%gjj(1,1,:) = geocont%gjj(i,1,:)
    geom_in%gjz(1,1,:) = geocont%gjz(i,1,:)
    geom_in%gzz(1,1,:) = geocont%gzz(i,1,:)
    geom_in%jacobian(1,1,:) = geocont%jacobian(i,1,:)
    geom_in%Bfield(1,1,:)   = geocont%Bfield(i,1,:)
    geom_in%dBdi(1,1,:)= geocont%dBdi(i,1,:)
    geom_in%dBdz(1,1,:)= geocont%dBdz(i,1,:)
    geom_in%R(1,:)= geocont%R(i,:)
    geom_in%Z(1,:)= geocont%Z(i,:)
    geom_in%dxdR(1,:)= geocont%dxdR(i,:)
    geom_in%dxdZ(1,:)= geocont%dxdZ(i,:)

  END SUBROUTINE assign_local_metric

end Module libITMgene_aux
