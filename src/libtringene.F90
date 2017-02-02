!#include "epik_user.inc"
#include "redef.h"

!>GENE-TRINITY interface
!!
!!\TODO in GENE:
!!* momentum flux?
!!* global ExB
!!
!! \param mpi_comm (sub-) mpi_comm_world provided by Trinity for this GENE instance
!! \param n_fluxtube_in Flux tube index
!! \param n_transp_it_in Transport code iteration number
!! \param electrostatic Switch to enforce electrostatic runs
!! \param px0_in number of radial grid points used for temperature/density profiles (=1 for local code)
!! \param px0_out number of bins for radially averaged output profiles
!! \param n_spec_in number of species (has to be <= number of species namelist in parameter file)
!! \param temp_io Array (rad. grid, n_spec_in) containing the temperature
!! \param dens_io Array (rad. grid, n_spec_in) containing the densities
!! \param mass_in Array containing species masses
!! \param vrot_in local code: ExBrate ExB shearing rate; Definition: 
!! \f$-\frac{c_{\rm ref}}{L_{\rm ref}}\cdot \frac{\rho_{{\rm tor},0}}{q_0}
!! \frac{\partial{\Omega_{\rm tor}}}{\partial \rho_{\rm tor}}\f$, 
!! global code: vrot in rad/s
!! \param charge_in Array containing species charges
!! \param omt_io Array containing the temp. gradients
!! \param omn_io Array containing the dens. gradients
!! \param q_in Safety factor
!! \param shat_in Magnetic shear
!! \param aspr_in aspect ratio to renormalize to major R for s_alpha or circular model
!! \param Lref_in reference length only required for analytical models else used for consistency check
!! \param Bref_in reference magnetic field only required for analytical models else used for consistency check
!! \param o_dVdx Radial derivative of volume in predefined bins in units of [Lref^2]
!! \param o_sqrtgxx_FS Flux surface avg. of sqrt(gxx) radially averaged in predefined bins
!! \param o_avgpflx Time avg'd.particle flux radially averaged in predefined bins
!! \param o_avgqflx Time avg'd heat flux radially averaged in predefined bins
!! \param o_temp_bin Time avg'd temperature profile radially averaged in predefined bins
!! \param o_dens_bin Time avg'd density profile radially averaged in predefined bins
subroutine tringene (mpi_comm_in,n_fluxtube_in,n_transp_it_in, electrostatic, & !administrative stuff
     &px0_in,px0_out,rho_in, rho_out,& !grids
     &n_spec_in, temp_io, dens_io, mass_in, charge_in, vrot_in, & !profiles
     &omt_io, omn_io, & !only required in local sims
     &q_in, shat_in, aspr_in, Lref_in, Bref_in,& !only required for analytical geometries
     &o_dVdx, o_sqrtgxx_FS, o_avgpflx, o_avgqflx, o_temp_bin, o_dens_bin) ! return values
  !OBSOLETE: !! \param rad_pos_in local code: radial position rho, global code: center positions of bins in rho
  Use mpi
  Use gene_subroutine
  Use parameters_IO
  Use communications
  use profiles_mod, only: set_spec_profiles, unset_spec_profiles, set_ref_vals
  Use external_contr, only: ExBrate
  use lagrange_interpolation
  
  Implicit None
  
  integer :: mpi_comm_in ! mpi communicator associated with the gene simulation
  integer, intent (in) :: n_fluxtube_in, n_transp_it_in ! for GENE output files index
  integer, intent (in) :: px0_in ! number of radial points used for input profiles (1 for local)
  integer, intent (in) :: px0_out ! number of radial points at which each GENE realization should return fluxes, etc. (1 for local)
  logical, intent (in) :: electrostatic
  real, dimension (px0_in), intent (in) :: rho_in   !grid for input profiles (r/a for analytical geom., rho_tor else)
  real, dimension (px0_out), intent (in) :: rho_out !grid for output profiles (r/a for analytical geom., rho_tor else)
  integer, intent (in) :: n_spec_in !number of species
  real, dimension (px0_in,0:n_spec_in-1), intent (inout) :: temp_io, dens_io !temperature/density profile in keV and 1E18m^-3
  real, dimension (0:n_spec_in-1), intent (in) :: mass_in, charge_in !mass in proton mass, charge in elementary charge
  real, dimension (px0_in), intent (in) :: vrot_in !tor. angular vel. vrot (global) or -ExBrate (local)
  real, dimension (px0_in,0:n_spec_in-1), intent (inout) :: omt_io, omn_io !input will only be used in local sims.
  real, dimension (px0_in), intent (in) :: q_in, shat_in !for analytical geometries *only*
  real, intent(in) :: aspr_in, Lref_in, Bref_in !for analytical geometries and consistency checks *only*
  ! following are binned (averaged) mean quantities and gradients from global GENE sim on rad_out_gene grid
  real, dimension (px0_out), intent (out) :: o_dVdx, o_sqrtgxx_FS
  real, dimension (px0_out,0:n_spec_in-1), intent (out) :: o_avgpflx, o_avgqflx
  real, dimension (px0_out,0:n_spec_in-1), intent (out) :: o_temp_bin, o_dens_bin

  !local variables
  character(len=128):: TRIN_par_in_dir='',TRIN_checkpoint_in=''
  character(len=20):: TRIN_file_extension
  integer :: file_unit = -1, is, rank, ierr
  call mpi_comm_rank (mpi_comm_in,rank, ierr)

  my_sim = n_fluxtube_in
  dens_io = 10.0*dens_io  !trinity uses density in 1E20m^-3

  if (n_spec_in.gt.1) then
     !enforce quasineutrality
     dens_io(:,n_spec_in-1) = 0.0
     do is = 0, n_spec_in-2
        dens_io(:,n_spec_in-1) = dens_io(:,n_spec_in-1)-&
             &(charge_in(is)/charge_in(n_spec_in-1))*dens_io(:,is)
     enddo
  endif

  !change file extension .dat to numbers
  WRITE(TRIN_file_extension,"(2(A,I3.3))") "_",n_fluxtube_in,"_",&
       &n_transp_it_in

  !set all n_procs to 1 to "fool" the check_parameters
  !subroutine; the actual values will be read later in the program
  n_procs_s = 1; n_procs_v = 1; n_procs_w = 1
  n_procs_x = 1; n_procs_y = 1; n_procs_z = 1

  TRIN_par_in_dir = ''
  print_ini_msg = .false.

  x0 = 0.5 ! --- need to look for a better solution (tbg)
  call read_parameters('')
  if (magn_geometry.eq.'tracer') then
     file_unit = -1
     call read_tracer_namelist(file_unit,.true.,.false.)
  endif

  !check compatibility with input parameters
  if (n_spec.ne.n_spec_in) STOP &
       'Number of species in parameters file and interface do not match!'

  call Modify_parameters
  TRIN_par_in_dir = 'skip_parfile'
  call check_params(TRIN_par_in_dir)

  call read_parall_nml('')

  call rungene(mpi_comm_in,TRIN_par_in_dir,TRIN_file_extension,&
       &TRIN_checkpoint_in)

  call Prepare_Output

  if (allocated(in_profiles)) deallocate(in_profiles)
  if (allocated(in_qprof)) deallocate(in_qprof)
  if (allocated(spec)) deallocate(spec)
  if (allocated(dVdx)) deallocate(dVdx)
  if (allocated(sqrtgxx_FS)) deallocate(sqrtgxx_FS)  
  if (allocated(area_FS)) deallocate(area_FS)

  
  !*************************************************************************!
  !************************* End of main program  **************************!
  !*************************************************************************!

Contains
  
  !>Replace parameters read from the parameters input file
  !!by values passed from Trinity
  Subroutine Modify_parameters
    IMPLICIT NONE
    Integer :: n
    real :: mp_o_e = 1.043968E-8 !kg/C

    x_local = (px0_in.eq.1)

    if (x_local) then
       x0 = rho_in(1)
    else
       x0 = 0.5
       allocate(in_profiles(1:px0_in,0:n_spec-1,0:4))
       mag_prof = .true.
       if (init_cond.eq.'ppj') init_cond='db'
    endif

    if ((trim(magn_geometry).ne.'tracer').and.&
         &(trim(magn_geometry).ne.'tracer_efit')) then
       if (x_local) then
          trpeps = rho_in(1)/aspr_in
          q0 = q_in(1)
          shat = shat_in(1)
       else
          allocate(in_qprof(1:px0_in))
          in_qprof = q_in
       endif
    endif

    !basic species settings
    do n=0,n_spec-1
       if (x_local) then
          spec(n)%dens = dens_io(1,n)
          spec(n)%temp = temp_io(1,n)
          spec(n)%omt  = omt_io(1,n)
          spec(n)%omn  = omn_io(1,n)
       else
          spec(n)%dens = 1.0
          spec(n)%temp = 1.0
          spec(n)%omt = 0.0
          spec(n)%omn = 0.0
          spec(n)%prof_type = 10
          in_profiles(:,n,0) = rho_in
          in_profiles(:,n,1) = temp_io(:,n)
          in_profiles(:,n,2) = omt_io(:,n) !gradients are calculated self-consist.
          in_profiles(:,n,3) = dens_io(:,n)
          in_profiles(:,n,4) = omn_io(:,n) !gradients are calculated self-consist.
       endif
       spec(n)%mass = mass_in(n)
       spec(n)%charge = charge_in(n)
    enddo


    !evaluate reference values!
    Tref = -1
    nref = -1
    mref = -1

    Lref = Lref_in !will be overwritten for non-analytic equil.
    Bref = abs(Bref_in) !will be overwritten for non-analytic equil.
    norm_index = 0

!if (rank.eq.0) then
!do n=1,px0_in
!   write(*,'(I3,5F12.6)') my_sim,rho_in(n), temp_io(n,0), dens_io(n,0),&
!        temp_io(n,1), dens_io(n,1)
!enddo
!write(*,*)
!endif

    if (x_local) then
       call set_ref_vals(rho_in,px0_in,spec,n_spec)
    else
       !temporarily set species profiles to trinity profile
       !to get reference values as we need them here to compute rhostat
       call set_spec_profiles(spec,rho_in,n_spec,0,px0_in-1,.false.)
       Call unset_spec_profiles(spec,n_spec)
    endif

    if (electrostatic) then
       beta = 1E-4
    else
       beta = -1
    endif
    coll = -1
    debye2 = -1
    rhostar = (SQRT(mref*1000*Tref*mp_o_e) / Bref) / (minor_r*Lref)
    
!if (rank.eq.0) write(*,'(I2,A,6F12.6)') my_sim, 'Tref,nref,mref,Lref,Bref,rho*: ', &
!         &Tref, nref, mref, Lref, Bref, rhostar

    if (.not.x_local) then
       lx = min(max(0.9,1.0-0.99*rho_out(1),1.01*rho_out(px0_out)),0.99)/rhostar
!       if (nx0.lt.floor(lx)) nx0 = 2.0*lx
       if (istep_prof.le.0) istep_prof=50
       if (ck_heat.le.0.0) ck_heat = 0.05
       if ((n_spec.gt.1).and.(ck_part.le.0.0)) ck_part = 0.05
       !ExBrate = -1111.0       
    else
       limit_shat = .true.
       if (ExBrate.ne.0) ExBrate = -vrot_in(1)
       !IMPORTANT: trinity assumes some q profile shape if q is
       !not found in the database file and hence the ExBrate might 
       !significantly differ from the values one would gain with the q from
       !EFIT; this should maybe fixed in the future (low priority as one
       !can alternatively take q from efit and add it to the iterdb file)
    endif

    !switch on calculation of time averaged fluxes if not set 
    !in parameters
    if (avgflux_stime.lt.0.0) avgflux_stime = 0.0
    if (.not.x_local) then
       if (avgprof_stime.lt.0.0) avgprof_stime = 0.0
    endif

    !check for checkpoints: first, try to take those with the same file extension
    !(i.e. continue a GENE run)
    !if not present, try to take those from the previous iteration
    !(i.e. resume with slightly modified gradients -> avoid big overshoots??)
    if (read_checkpoint) then
       reset_chpt_time = .false.
       TRIN_checkpoint_in = TRIM(chptdir)//'/checkpoint'//TRIN_file_extension
       if (.not.valid_chpt(mpi_comm_in,TRIN_checkpoint_in)) then
          !try secure checkpoint
          TRIN_checkpoint_in = TRIM(chptdir)//'/s_checkpoint'//TRIN_file_extension
          if (.not.valid_chpt(mpi_comm_in,TRIN_checkpoint_in)) then
             !try checkpoint from previous run
             WRITE(TRIN_checkpoint_in,"(2A,(I2.2,A,I2.2))") TRIM(chptdir),'/checkpoint_',&
                  & n_fluxtube_in,"_",(n_transp_it_in-1)
             if (.not.valid_chpt(mpi_comm_in,TRIN_checkpoint_in)) then
                TRIN_checkpoint_in = 'no' !i.e. switch off checkpoint read
             else
                reset_chpt_time = .true.
             endif
          endif
       endif
    else
       TRIN_checkpoint_in='no'
    endif

    print_ini_msg = .false.
    multiple_tracer_files = .true.

  End Subroutine Modify_parameters
  
  !> A routine setting the fluxes to
  !! pflux = dfac*a/Ln
  !! qflux = dfac*(p_s/p_i)*1.5*a/Lp_s
  !! Lp_s ~ L_n + L_T (?)
  !! for comparison with trinity
  !! flux_option = 'test1'
  !! grad_option = 'ntgrads'
  !! dfac can be set in trinity's fluxes
  !! namelist
  subroutine test_flux(n)
    integer :: n
    real :: dfac=0.3

    if (x_local) then
       do n=0, n_spec-1
          o_avgpflx(1,n) = dfac*spec(n)%omn
          o_avgqflx(1,n) = dfac*(spec(n)%dens*spec(n)%temp)/&
               &(spec(0)%dens*spec(0)%temp)*1.5*&
               &(spec(n)%omn+spec(n)%omt)
       enddo
    else
       !set avgprof(3,4) to omn/omt_prof in diag_df!
       avgprof(:,n,5) = dfac*avgprof(:,n,3)
       avgprof(:,n,6) = dfac*(avgprof(:,n,2)*avgprof(:,n,1))/&
            &(avgprof(:,0,2)*avgprof(:,0,1))*1.5*&
            &(avgprof(:,n,4)+avgprof(:,n,3))
       !this model is local hence change to global:
       avgprof(:,n,5) = avgprof(:,n,5)*(avgprof(:,0,1)*&
            &(avgprof(:,0,2))**1.5)
       avgprof(:,n,6) = avgprof(:,n,6)*(avgprof(:,0,1)*&
            &(avgprof(:,0,2))**2.5)
    endif
  end subroutine test_flux

  !>Prepare return values/arrays for Trinity
  Subroutine Prepare_Output
    IMPLICIT NONE

    integer :: n, ix
    Real, dimension (px0_out,0:n_spec-1) :: o_omn_bin, o_omt_bin
    
    o_avgpflx = 0.0
    o_avgqflx = 0.0

    if (x_local) then
       if (comp_type.eq.'NC') then
          do n=0, n_spec-1
             o_avgpflx(1,n) = neoflux(pi1,n,1)
             o_avgqflx(1,n) = neoflux(pi1,n,2)
          enddo
       elseif ((time-simtime_start).gt.0) then
          do n=0, n_spec-1
             o_avgpflx(1,n) = avgfluxes(n,1)
             o_avgqflx(1,n) = avgfluxes(n,2)
          enddo
       endif
       !call test_flux(n)

       do n=0,n_spec-1
          o_temp_bin(:,n) = spec(n)%temp
          o_dens_bin(:,n) = spec(n)%dens
          o_omn_bin(:,n) = spec(n)%omn
          o_omt_bin(:,n) = spec(n)%omt
       enddo
       o_dVdx = dVdx
       o_sqrtgxx_FS=sqrtgxx_FS
    else !nonlocal
       if ((time-avgprof_stime).ge.0.0) then
          do n=0,n_spec-1
             avgprof(:,n,2) = avgprof(:,n,2)*spec(n)%temp
             avgprof(:,n,1) = avgprof(:,n,1)*spec(n)%dens

!return current values instead of trinity values?
!maybe, in the future - for now, keep input values
!             call lag3interp(avgprof(:,n,2),avgprof(:,n,0),&
!                  &size(avgprof(:,n,0)),temp_io(:,n),rho_in,&
!                  &px0_in)
!             call lag3interp(avgprof(:,n,1),avgprof(:,n,0),&
!                  &size(avgprof(:,n,0)),dens_io(:,n),rho_in,&
!                  &px0_in)
!add units (see below)

             call lag3interp(avgprof(:,n,4),avgprof(:,n,0),&
                  &size(avgprof(:,n,0)),omt_io(:,n),rho_in,&
                  &px0_in)
             call lag3interp(avgprof(:,n,3),avgprof(:,n,0),&
                  &size(avgprof(:,n,0)),omn_io(:,n),rho_in,&
                  &px0_in)

             call calc_avg_in_bins(avgprof(:,n,2),avgprof(:,0,0),&
                  size(avgprof(:,0,0)),&
                  &o_temp_bin(:,n),rho_out,px0_out)
             call calc_avg_in_bins(avgprof(:,n,1),avgprof(:,0,0),&
                  size(avgprof(:,0,0)),&
                  &o_dens_bin(:,n),rho_out,px0_out)
 
             call calc_avg_in_bins(avgprof(:,n,4),avgprof(:,0,0),&
                  size(avgprof(:,0,0)),&
                  &o_omt_bin(:,n),rho_out,px0_out)
             call calc_avg_in_bins(avgprof(:,n,3),avgprof(:,0,0),&
                  size(avgprof(:,0,0)),&
                  &o_omn_bin(:,n),rho_out,px0_out)
             !call test_flux(n)

             if (comp_type.ne.'NC') then
                call calc_avg_in_bins(avgprof(:,n,5),avgprof(:,0,0),&
                     size(avgprof(:,0,0)),&
                     &o_avgpflx(:,n),rho_out,px0_out)
                call calc_avg_in_bins(avgprof(:,n,6),avgprof(:,0,0),&
                     size(avgprof(:,0,0)),&
                     &o_avgqflx(:,n),rho_out,px0_out)
             else
                call calc_avg_in_bins(neoflux(:,n,1),avgprof(:,0,0),&
                     size(avgprof(:,0,0)),&
                     &o_avgpflx(:,n),rho_out,px0_out)
                call calc_avg_in_bins(neoflux(:,n,2),avgprof(:,0,0),&
                     size(avgprof(:,0,0)),&
                     &o_avgqflx(:,n),rho_out,px0_out)
             endif

             !change normalization from center value to 'local' 
             !normalization in each bin (always w.r.t. to first species)
             o_avgpflx(:,n) = o_avgpflx(:,n)/&
                  &(o_dens_bin(:,0)*(o_temp_bin(:,0))**1.5)
             o_avgqflx(:,n) = o_avgqflx(:,n)/&
                  &(o_dens_bin(:,0)*(o_temp_bin(:,0))**2.5)
          enddo
       endif

       call calc_avg_in_bins(dVdx,avgprof(:,0,0),&
            size(avgprof(:,0,0)),&
            &o_dVdx,rho_out,px0_out)
       call calc_avg_in_bins(sqrtgxx_FS,avgprof(:,0,0),&
            size(avgprof(:,0,0)),&
            &o_sqrtgxx_FS,rho_out,px0_out)
    endif !x_local

    avgflux_stime = -1.0
    avgprof_stime = -1.0

    !add trinity units to profiles (except fluxes)
    o_temp_bin = o_temp_bin * Tref
    o_dens_bin = o_dens_bin * 0.1*nref
    dens_io = 0.1*dens_io

#if 0
if (rank==0) then
!   do ix = 1, px0_out
!      write (*,'(A,I3,7(F12.6))') "OUT-BIN",my_sim,rho_out(ix),&
!           &o_dens_bin(ix,0),o_temp_bin(ix,0), &
!           &o_dens_bin(ix,1),o_temp_bin(ix,1)
!   end do
!   write(*,*)
!   write(*,*)
   do ix = 1, px0_out
      write (*,'(A,I3,5(F8.3))') "GRD-BIN",my_sim,rho_out(ix),&
           &o_omn_bin(ix,0),o_omt_bin(ix,0), &
           &o_omn_bin(ix,1),o_omt_bin(ix,1)
   end do
   write(*,*)
   write(*,*)

   do ix = 1, px0_out
      write (*,'(A,I3,5(F8.3))') "FLX-BIN",my_sim,rho_out(ix),&
           &o_avgpflx(ix,0),o_avgqflx(ix,0), &
           &o_avgpflx(ix,1),o_avgqflx(ix,1)
   end do
   write(*,*)
   write(*,*)
endif
#endif

End Subroutine Prepare_Output
  
  !>Calculates radial averages in bins being
  !!centered between neighbouring output grid points
  Subroutine calc_avg_in_bins(inarr,rad_in,n_in,outarr,rad_out,n_out)
    
    Implicit None

    Integer, Intent(in) :: n_in, n_out

    Real, dimension(n_in) :: inarr, rad_in
    Real, dimension(n_out) :: outarr, rad_out

    Integer :: i,j,count
    Real :: xstart, xend
    
    outarr = 0.0

    Do i=1, n_out
       if (i.eq.1) then
          xstart = rad_in(1)
       else
          xstart = 0.5*(rad_out(i-1)+rad_out(i))
       endif
       if (i.eq.n_out) then
          xend = rad_in(n_in)
       else
          xend = 0.5*(rad_out(i)+rad_out(i+1))
       endif
       count = 0
       Do j=1,n_in
          if ((rad_in(j).ge.xstart).and.(rad_in(j).le.xend)) then
             count = count+1
             outarr(i) = outarr(i) + inarr(j)
          else
             cycle
          endif
       enddo
       if (count.eq.0) stop 'binning failed - use more input grid points'
       outarr(i) = outarr(i)/count
    Enddo

  End Subroutine calc_avg_in_bins


  logical function valid_chpt(mpi_comm_in,chpt_str)
    Implicit None
    integer,intent(in) :: mpi_comm_in
    character(len=*),intent(in) :: chpt_str
    integer :: handle=MPI_FILE_NULL, ierr
    integer(MPI_OFFSET_KIND) :: offset

    inquire(file=TRIM(chpt_str),exist=valid_chpt)

    if ((.not.chpt_read_h5).and.(valid_chpt)) then !check for file size
       call mpi_file_open(mpi_comm_in,TRIM(chpt_str),&
            &MPI_MODE_RDONLY,MPI_INFO_NULL,handle,ierr)
       !check file size
       call mpi_file_get_size(handle,offset,ierr)
       call mpi_file_close(handle,ierr)
       if (abs(offset).le.6) valid_chpt = .false.
       !using abs(offset) for now as OpenMPI has an issue with files sizes >2GB
       !and returns neg. numbers: https://svn.open-mpi.org/trac/ompi/ticket/2145
    endif
    
  end function valid_chpt



End Subroutine tringene
