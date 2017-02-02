#include "intrinsic_sizes.h"
#include "redef.h"
!>Read profiles from ASCII files and compute, e.g., the
!!normalized logarithmic gradients
MODULE profiles_mod
  use par_mod
  use geometry
  use coordinates, only: xval, xval_a
  use lagrange_interpolation
  use x_derivatives, only: calc_gradient
  use profile_io
  use file_io, only: get_unit_nr, get_nr_of_lines,&
       &get_nr_of_substr

  IMPLICIT NONE

  PUBLIC :: set_spec_profiles, unset_spec_profiles
  PUBLIC :: write_tempdens_profiles
  PUBLIC :: calc_derived_ref_vals, calc_dpdx_pm
  PUBLIC :: Tref, nref, mref, cref, drive_buffer, set_prof_defaults
  PUBLIC :: set_in_prof_defaults, unset_in_prof_defaults
  PUBLIC :: ldrive_buffer_size, udrive_buffer_size, ldrive_coef
  PUBLIC :: udrive_coef, ldrive_pow, udrive_pow
  PUBLIC :: iterdb_file, iterdb_time
  PUBLIC :: norm_index, set_ref_vals
  PUBLIC :: omegator, omegatorprime, omegatorref, zeff_prof, zeff, zeff0

  REAL :: cref=0.0
  REAL :: zeff  !input parameter (can be negative to determine zeff0 from input profile)
  REAL :: zeff0 !zeff at reference position, used in collisions, can not be negative
  REAL :: in_Tref=0, in_nref=0, in_beta=-10, in_coll=-10, in_debye2=-10, in_rhostar=-10
  REAL, dimension(:),allocatable :: in_temp, in_dens
  REAL, dimension(:),allocatable :: omegator, omegatorprime, zeff_prof
  logical:: drive_buffer=.false.
  REAL:: ldrive_buffer_size, udrive_buffer_size, ldrive_coef
  REAL:: udrive_coef, ldrive_pow, udrive_pow
  integer::init_status=0
  
  PRIVATE



CONTAINS

  subroutine set_prof_defaults
     Tref=0.
     nref=0.
     in_Tref=0.
     in_nref=0.
     in_beta=-10
     in_coll=-10
     in_debye2=-10
     in_rhostar=-10
     mref=0.0
     drive_buffer=.false.
     iterdb_file = ''
     iterdb_time = -1.0
     omegatorref = 0.0
     zeff = 1.0
  end subroutine set_prof_defaults
  
  SUBROUTINE set_in_prof_defaults
    
    if (init_status==0) then
       allocate(in_temp(0:n_spec-1), in_dens(0:n_spec-1))
       init_status = 1
    endif
    in_temp=-10
    in_dens=-10
  
  END SUBROUTINE set_in_prof_defaults

  SUBROUTINE unset_in_prof_defaults

    deallocate(in_temp,in_dens)
    init_status = 0

  END SUBROUTINE unset_in_prof_defaults

  subroutine set_spec_profiles(spec, rad_out, n_spec, pi1gl, pi2gl, write_pe)
    INTEGER, intent(in) :: n_spec, pi1gl, pi2gl
    TYPE(spectype) :: spec(0:n_spec-1)
    real, dimension(pi1gl:pi2gl), intent(in) :: rad_out
    logical :: prof_from_file, write_pe
    integer :: n, ppx0
    real, dimension(:), allocatable:: tmparr

    px0 = pi2gl-pi1gl+1

    prof_from_file = .false.

    zeff0=zeff
    do n=0, n_spec-1
       Allocate(spec(n)%temp_prof(pi1gl:pi2gl),spec(n)%dens_prof(pi1gl:pi2gl),&
            &spec(n)%omt_prof(pi1gl:pi2gl),spec(n)%omn_prof(pi1gl:pi2gl))
       if (allocated(in_profiles)) then
          allocate(tmparr(pi1gl:pi2gl))
          !interpolate from external to gene grid
          ppx0 = size(in_profiles(:,n,0))
          !--- temperature
          call lag3interp(in_profiles(:,n,1), in_profiles(:,n,0),ppx0,&
               spec(n)%temp_prof, rad_out, px0)
          tmparr=-LOG(spec(n)%temp_prof)
          CALL lag3deriv(tmparr, rad_out, px0, &
               &spec(n)%omt_prof, rad_out, px0)
          !          call lag3interp(in_profiles(:,n,2), in_profiles(:,n,0),ppx0,&
          !               spec(n)%omt_prof, rad_out, px0)          
          !--- density
          call lag3interp(in_profiles(:,n,3), in_profiles(:,n,0),ppx0,&
               spec(n)%dens_prof, rad_out, px0)
          tmparr=-LOG(spec(n)%dens_prof)
          CALL lag3deriv(tmparr, rad_out, px0, &
               &spec(n)%omn_prof, rad_out, px0)
          !          call lag3interp(in_profiles(:,n,4), in_profiles(:,n,0),ppx0,&
          !               spec(n)%omn_prof, rad_out, px0)          
          deallocate(tmparr)
          prof_from_file = .true.
       else
          Call set_spec_profiles_n(spec(n),minor_r,write_pe)
       endif
       prof_from_file = prof_from_file.or.(spec(n)%prof_type.lt.0)
    end do

    if (spec(0)%prof_type.lt.0) call enforce_quasineutrality(spec,write_pe)

    if ((prof_from_file).or.(x_local.and.spec(0)%prof_type>0)) then
       call set_ref_vals(rad_out,px0,spec,n_spec)
       if (minor_r.ne.1..and.write_pe) then
          write (*,'(A)') 'WARNING: Assuming input profiles to be given in&
               & x/a coordinates!'
       endif
    endif
    if ((mype==0).and.(.not.x_local)) then
       do n=0, n_spec-1
          Call write_tempdens_profiles(spec(n),.false.)
       enddo
    endif
    
  END SUBROUTINE set_spec_profiles

  SUBROUTINE set_spec_profiles_n(sp, minor_r,write_pe)
    type(spectype) :: sp
    REAL, intent(in) :: minor_r
    logical::write_pe
    real, dimension(:), allocatable:: tmparr
    INTEGER :: prof_type, px0

    prof_type = sp%prof_type
    px0 = SIZE(sp%temp_prof)

    IF (prof_type.ne.0) THEN
       IF (px0.NE.SIZE(xval_a)) &
            & stop 'inconsistent array lengths in set_spec_profiles'
    ENDIF

    !temperature and density profiles
    select case (prof_type)
    case(1)
       IF (sp%LT_center.ne.0.5) STOP 'currently only LT_center=0.5 allowed for prof_type=1'
       IF (sp%Ln_center.ne.0.5) STOP 'currently only Ln_center=0.5 allowed for prof_type=1'
       !Mathematica result for integrate(omt(x)) provided by Lausanne
       sp%temp_prof = exp(sp%kappa_T*minor_r/(sinh(sp%LT_center/sp%LT_width)**2)*(&
            xval_a-sp%LT_center-sp%LT_width*cosh(sp%LT_center/sp%LT_width)**2*&
            tanh((xval_a-sp%LT_center)/sp%LT_width)))
       sp%omt_prof = sp%kappa_T*(cosh((xval_a-sp%LT_center)/sp%LT_width)**(-2)-&
            &cosh(sp%LT_center/sp%LT_width)**(-2))/(1.-cosh(sp%LT_center/sp%LT_width)**(-2))
       
       sp%dens_prof = exp(sp%kappa_n*minor_r/(sinh(sp%Ln_center/sp%Ln_width)**2)*(&
            xval_a-sp%Ln_center-sp%Ln_width*cosh(sp%Ln_center/sp%Ln_width)**2*&
            tanh((xval_a-sp%Ln_center)/sp%Ln_width)))
       sp%omn_prof = sp%kappa_n*(cosh((xval_a-sp%Ln_center)/sp%Ln_width)**(-2)-&
               &cosh(sp%Ln_center/sp%Ln_width)**(-2))/(1-cosh(sp%Ln_center/sp%Ln_width)**(-2))
    case(2)
       sp%temp_prof = exp(-sp%kappa_T*minor_r*sp%LT_width*tanh((xval_a-sp%LT_center)/sp%LT_width))
       sp%dens_prof = exp(-sp%kappa_n*minor_r*sp%Ln_width*tanh((xval_a-sp%Ln_center)/sp%Ln_width))
       sp%omt_prof = sp%kappa_T*cosh((xval_a-sp%LT_center)/sp%LT_width)**(-2)
       sp%omn_prof = sp%kappa_n*cosh((xval_a-sp%Ln_center)/sp%Ln_width)**(-2)
    case(3)
       sp%temp_prof = (cosh((xval_a-sp%LT_center+sp%delta_x_T)/sp%LT_width)/&
            &cosh((xval_a-sp%LT_center-sp%delta_x_T)/sp%LT_width))**&
            &(-0.5*sp%kappa_T*minor_r*sp%LT_width)
       sp%omt_prof = 0.5*sp%kappa_T*(tanh((xval_a-sp%LT_center+sp%delta_x_T)/sp%LT_width)-&
            &tanh((xval_a-sp%LT_center-sp%delta_x_T)/sp%LT_width))
       sp%dens_prof = (cosh((xval_a-sp%Ln_center+sp%delta_x_n)/sp%Ln_width)/&
            &cosh((xval_a-sp%Ln_center-sp%delta_x_n)/sp%Ln_width))**(-0.5*sp%kappa_n*minor_r*sp%Ln_width)
       sp%omn_prof = 0.5*sp%kappa_n*(tanh((xval_a-sp%Ln_center+sp%delta_x_n)/sp%Ln_width)-&
            &tanh((xval_a-sp%Ln_center-sp%delta_x_n)/sp%Ln_width))
    case(4) !Falchetto et al., PPCF 50, 124015
       !The above cited paper mixes the definitions of r0 and ra
       !here, r0=center is the reference radius (usually 0.5a), ri=r0-0.5*delta and ra=r0+0.5*delta;
       !the values ri=0.1 and ra=0.9 therefore translate to center=0.5 and delta=0.8
       !the width is set to 0.04 in the paper
       !furthermore, cosh^(-2) = sech^2 has been used
       sp%temp_prof = exp(-sp%kappa_T*minor_r*(xval_a-sp%LT_center-sp%LT_width*&
            &(tanh((xval_a-sp%LT_center-0.5*sp%delta_x_T)/sp%LT_width)+&
           &tanh((xval_a-sp%LT_center+0.5*sp%delta_x_T)/sp%LT_width))))
       sp%omt_prof = sp%kappa_T*(1.-cosh((xval_a-sp%LT_center-0.5*sp%delta_x_T)/sp%LT_width)**(-2)-&
            &cosh((xval_a-sp%LT_center+0.5*sp%delta_x_T)/sp%LT_width)**(-2))
       sp%dens_prof = exp(-sp%kappa_n*minor_r*(xval_a-sp%Ln_center-sp%Ln_width*&
            &(tanh((xval_a-sp%Ln_center-0.5*sp%delta_x_n)/sp%Ln_width)+&
            &tanh((xval_a-sp%Ln_center+0.5*sp%delta_x_n)/sp%Ln_width))))
       sp%omn_prof = sp%kappa_n*(1.-cosh((xval_a-sp%Ln_center-0.5*sp%delta_x_n)/sp%Ln_width)**(-2)-&
            &cosh((xval_a-sp%Ln_center+0.5*sp%delta_x_n)/sp%Ln_width)**(-2))
    case(5) !exponential profiles
       sp%temp_prof = exp(-sp%kappa_T*minor_r*(xval_a-sp%LT_center))
       sp%omt_prof  = sp%kappa_T
       sp%dens_prof = exp(-sp%kappa_n*minor_r*(xval_a-sp%Ln_center))
       sp%omn_prof  = sp%kappa_n          
    case(-1) !read from file
       if (trim(sp%prof_file).eq.'') &
            &sp%prof_file='./profiles_'//trim(sp%name)
       call read_numerical_profiles(sp%prof_file,px0,sp,write_pe)
    case(-2) !read from iterdb file
       call read_numerical_profiles(iterdb_file,px0,sp,write_pe)
    case(-3) !read from d3d_iterdb file
       call read_numerical_profiles(iterdb_file,px0,sp,write_pe)
    case(-4) !read from CHEASE ogyropsi.h5
       call read_numerical_profiles('',px0,sp,write_pe)
    case default !local limit
       sp%temp_prof(:) = 1.0
       sp%dens_prof(:) = 1.0
       sp%omt_prof(:) = sp%omt
       sp%omn_prof(:) = sp%omn
    end select

    if (drive_buffer.and..not.x_local) then
       allocate(tmparr(px0))
       call profile_buffer(sp%temp_prof)
       call profile_buffer(sp%dens_prof)
       tmparr=-LOG(sp%temp_prof)
       CALL lag3deriv(tmparr, xval, nx0, &
            &sp%omt_prof, xval, nx0)
       tmparr=-LOG(sp%dens_prof)
       CALL lag3deriv(tmparr, xval, nx0, &
            &sp%omn_prof, xval, nx0)
       sp%omt_prof=sp%omt_prof/(rhostar*minor_r)
       sp%omn_prof=sp%omn_prof/(rhostar*minor_r)
       deallocate(tmparr)
    endif

  END SUBROUTINE set_spec_profiles_n


  SUBROUTINE unset_spec_profiles(spec,n_spec)
    integer, intent(in) :: n_spec
    type(spectype),intent(inout) :: spec(0:n_spec-1)
    integer :: n

    do n=0,n_spec-1
       deallocate(spec(n)%temp_prof,spec(n)%dens_prof,&
            spec(n)%omt_prof,spec(n)%omn_prof) 
    enddo 
    if (allocated(omegator)) deallocate(omegator,omegatorprime)
    if (allocated(zeff_prof)) deallocate(zeff_prof)

  END SUBROUTINE unset_spec_profiles
  

  SUBROUTINE profile_buffer(in_profile)
    real, dimension(pi1gl:pi2gl), intent(inout):: in_profile
    integer:: l_nxbuf, l_ibuf, u_nxbuf, u_ibuf, i
    real, dimension(pi1gl:pi2gl):: nu_K
    
    if (ldrive_coef.lt.10.) ldrive_coef=10.
    if (udrive_coef.lt.10.) udrive_coef=10.

    l_nxbuf=int(ldrive_buffer_size*nx0)
    l_ibuf =pi1gl+l_nxbuf
    u_nxbuf=int(udrive_buffer_size*nx0)
    u_ibuf =pi2gl-u_nxbuf
    
    if (l_nxbuf.gt.0) then
       do i=pi1gl, pi1gl+l_nxbuf
          nu_K(i)=ldrive_coef*(abs(i-l_ibuf)/real(l_nxbuf))**ldrive_pow
       end do
    endif
    
    if (u_nxbuf.gt.0) then
       do i=pi2gl-u_nxbuf, pi2gl
          nu_K(i)=udrive_coef*(abs(i-u_ibuf)/real(u_nxbuf))**udrive_pow
       end do
    endif

    do i=pi1gl+l_nxbuf,pi1gl,-1
       in_profile(i)=in_profile(i+1)-(in_profile(i+1)-in_profile(i))*exp(-nu_K(i))
    enddo
    do i=pi2gl-u_nxbuf,pi2gl
       in_profile(i)=in_profile(i-1)+(in_profile(i)-in_profile(i-1))*exp(-nu_K(i))
    enddo

  END SUBROUTINE profile_buffer


  !this routine should only be called on one (write) process
  SUBROUTINE write_tempdens_profiles(sp,append_1)
#ifdef WITHFUTILS
    use futils
#endif

    TYPE(spectype),intent(in) :: sp
    Character(Len=FILENAME_MAX) :: filename
    Logical, intent(in) :: append_1
    INTEGER :: i, fileunit
    REAL :: my_Tref, my_nref
#ifdef WITHFUTILS
    integer :: fidspecprof_h5
    character(len=FILENAME_MAX) :: dset_name_specprof
#endif

    if (.not.allocated(xval)) return

    my_Tref = Tref
    my_nref = nref

    if (my_Tref.eq.0) my_Tref=1.0
    if (my_nref.eq.0) my_nref=1.0

    if(write_std) then
       filename =TRIM(trim(diagdir)//'/profiles_'//trim(sp%name)//&
            &TRIM(file_extension))

       call get_unit_nr(fileunit)

       if (append_1) then
          OPEN(fileunit,file=trim(filename),position='APPEND',&
               &status='old',action='write')
          write(fileunit,"(A)") ""
          write(fileunit,"(A)") ""      
       else
          OPEN(fileunit,file=trim(filename))
          ! write file header
          write(fileunit,"(A)") "#   x/a             x/rho_ref           T"//&
               &"                   n                   omt                 omn"
       endif

       write(fileunit,"(A,F14.6)") '#', time
       Do i=pi1gl,pi2gl
          write(fileunit,"(6ES20.10)") xval_a(i), xval(i), &
               &my_Tref*sp%temp*sp%temp_prof(i), &
               &my_nref*sp%dens*sp%dens_prof(i), sp%omt_prof(i), sp%omn_prof(i)
       ENDDO
       close(fileunit)
    end if

#ifdef WITHFUTILS
    if(write_h5) then
       if(mype.eq.0) then
         call creatf(trim(diagdir)//'/profiles_'//trim(sp%name)//trim(file_extension)//'.h5', &
               fidspecprof_h5, "species profiles infos", 'd')
          
          call creatg(fidspecprof_h5, '/position')
          call creatg(fidspecprof_h5, '/temp')
          call creatg(fidspecprof_h5, '/density')
          call creatg(fidspecprof_h5, '/src')
          
          write(dset_name_specprof, "(A)") "/position/x_o_a"
          call putarr(fidspecprof_h5, dset_name_specprof, xval_a)
          write(dset_name_specprof, "(A)") "/position/x_o_rho_ref"
          call putarr(fidspecprof_h5, dset_name_specprof, xval_a/rhostar)          

          write(dset_name_specprof, "(A)") "/temp/T"
          call putarr(fidspecprof_h5, dset_name_specprof, my_Tref*sp%temp*sp%temp_prof)
          write(dset_name_specprof, "(A)") "/temp/omt"
          call putarr(fidspecprof_h5, dset_name_specprof, sp%omt_prof)
          
          write(dset_name_specprof, "(A)") "/density/n"
          call putarr(fidspecprof_h5, dset_name_specprof, my_nref*sp%dens*sp%dens_prof)
          write(dset_name_specprof, "(A)") "/density/omn"
          call putarr(fidspecprof_h5, dset_name_specprof, sp%omn_prof)
       
          call closef(fidspecprof_h5)
       end if
    end if
#endif

  END SUBROUTINE write_tempdens_profiles

  !>Generic read/interpolation routine for numerical profile input.
  subroutine read_numerical_profiles(filename,px0,sp,write_pe)
    Character(Len=*),intent(in)  :: filename
    integer, intent(inout):: px0
    type(spectype),intent(inout) :: sp    
    logical,intent(in):: write_pe
    integer:: nx0_file
    integer, dimension(1) :: i1
    real, dimension(:),allocatable :: temp_file, dens_file, xaxis_file, &
         &tmp_arr, omegator_file
    !for chease
    real, dimension(:),allocatable :: te_prof, ne_prof, ti_prof, ni_prof

    !For all profile types except IterDB (-2), we read the profiles directly from
    !file and interpolate them to the GENE grid in this routine. 
    select case(sp%prof_type)
    case(-1)
       call get_entrysize_gene(filename,nx0_file)
       allocate(xaxis_file(1:nx0_file),temp_file(1:nx0_file),&
            dens_file(1:nx0_file),tmp_arr(1:nx0_file),omegator_file(1:nx0_file))
       call read_profiles_gene(filename,nx0_file,temp_file,dens_file,&
            xaxis_file)
       !interpolate file data to GENE grid
       call lag3interp(temp_file,xaxis_file,nx0_file,&
            sp%temp_prof,xval_a,px0)
       tmp_arr = -log(temp_file)
       call lag3deriv(tmp_arr,xaxis_file,nx0_file,&
            sp%omt_prof,xval_a,px0)
       call lag3interp(dens_file,xaxis_file,nx0_file,&
            sp%dens_prof,xval_a,px0)
       tmp_arr = -log(dens_file)
       call lag3deriv(tmp_arr,xaxis_file,nx0_file,&
            sp%omn_prof,xval_a,px0)
       deallocate(xaxis_file,temp_file,dens_file,tmp_arr)
    case(-2)
       !NOTE: IterDB files can have different grids for each data set, and therefore have a special
       !treatment: If xval_a is not zero everywhere, the data from file will be interpolated to the
       !xval_a grid directly in the read_profiles_iterdb routine.
       if (sp%charge.lt.0.0) then !electrons
          call read_profiles_iterdb(iterdb_file, iterdb_time, 'TE', &
               &px0, sp%temp_prof, sp%omt_prof, xval_a)
          call read_profiles_iterdb(iterdb_file, iterdb_time, 'NE', &
               &px0, sp%dens_prof, sp%omn_prof, xval_a)
       else
          call read_profiles_iterdb(iterdb_file, iterdb_time, 'TI', &
               &px0, sp%temp_prof, sp%omt_prof, xval_a)
          !!! IMPROVE: 2nd ION SPECIES (helium) OR IMPURITIES, FAST IONS?       
          if (sp%charge.eq.1) then !should be improved in future
             call read_profiles_iterdb(iterdb_file, iterdb_time, 'NM1', &
                  &px0, sp%dens_prof, sp%omn_prof, xval_a)
          else
             call read_profiles_iterdb(iterdb_file, iterdb_time, 'NM2', &
                  &px0, sp%dens_prof, sp%omn_prof, xval_a)
          endif
       endif
       if (.not.allocated(omegator)) then
          allocate(omegator(pi1gl:pi2gl),omegatorprime(pi1gl:pi2gl))
          omegator = 0.0
          call read_profiles_iterdb(iterdb_file, iterdb_time,'VROT',&
               &px0,omegator,omegatorprime,xval_a,.true.)
          omegatorprime=-omegatorprime*omegator
       endif

       if (zeff.LT.0.0) then
          if (.not.allocated(zeff_prof)) then
             allocate(zeff_prof(pi1gl:pi2gl))
             allocate(temp_file(pi1gl:pi2gl))
             call read_profiles_iterdb(iterdb_file, iterdb_time,'ZEFFR',&
                  &px0,zeff_prof,temp_file,xval_a)
             if (px0.gt.1) then
                call lag3interp(zeff_prof,xval_a,px0,zeff0,x0)
                if (write_pe) Write(*,"(A,1ES13.4)") "read zeff profile: zeff(x0) = ",zeff0
             else
                zeff0 = zeff_prof(pi1)
                if (write_pe) Write(*,"(A,1ES13.4,A)") "read zeff = ", zeff0, " from profile"
             endif
          deallocate(temp_file)
          endif
       endif

       !convert to GENE reference units
       sp%temp_prof = 1E-3 * sp%temp_prof  !eV to keV
       sp%dens_prof = 1E-19 * sp%dens_prof !m^{-3} to 1E19 m^{-3}
       
    case(-3)
       call get_entrysize_iterdb_d3d(iterdb_file, nx0_file)
       allocate(xaxis_file(1:nx0_file),temp_file(1:nx0_file),&
            dens_file(1:nx0_file),tmp_arr(1:nx0_file),omegator_file(1:nx0_file))
       if (sp%charge.lt.0.0) then !electrons
          call read_profiles_iterdb_d3d(iterdb_file, &
               &'electron temperature, keV', nx0_file, temp_file, xaxis_file)
          call read_profiles_iterdb_d3d(iterdb_file, &
               &'electron density', nx0_file, dens_file, xaxis_file)
       else
          !note: there's indeed a typo in these files
          !may need special treatment once fixed at d3d
          call read_profiles_iterdb_d3d(iterdb_file, 'ion temperatue, keV', &
               &nx0_file, temp_file, xaxis_file)
          if (sp%charge.eq.1) then !should be improved in future
             call read_profiles_iterdb_d3d(iterdb_file, 'primary ion density', &
                  &nx0_file, dens_file, xaxis_file)
          else
             call read_profiles_iterdb_d3d(iterdb_file, &
                  &'impurity ion density', nx0_file, dens_file, xaxis_file)
          endif
       endif
       if (.not.allocated(omegator)) then
          allocate(omegator(pi1gl:pi2gl),omegatorprime(pi1gl:pi2gl))
          call read_profiles_iterdb_d3d(iterdb_file, &
               &'angular rotation speed', nx0_file, omegator_file, xaxis_file)
          call lag3interp(omegator_file,xaxis_file,nx0_file,&
               omegator,xval_a,px0)
          call lag3deriv(omegator_file,xaxis_file,nx0_file,&
               omegatorprime,xval_a,px0)
       endif
       if (zeff.LT.0.0) then
          if (.not.allocated(zeff_prof)) then
             allocate(zeff_prof(pi1gl:pi2gl))
             call read_profiles_iterdb_d3d(iterdb_file, &
                  &'zeff profile',nx0_file,tmp_arr,xaxis_file)
             call lag3interp(tmp_arr,xaxis_file,nx0_file,zeff_prof,xval_a,px0)
             if (px0.gt.1) then
                call lag3interp(tmp_arr,xaxis_file,nx0_file,zeff0,x0)
                if (write_pe) Write(*,"(A,1ES13.4)") "read zeff profile: zeff(x0) = ",zeff0
             else
                zeff0 = zeff_prof(pi1)
                if (write_pe) Write(*,"(A,1ES13.4,A)") "read zeff = ", zeff0, " from profile"
             endif
          endif
       endif

       !convert to GENE reference units
       dens_file = 1E-19 * dens_file !m^{-3} to 1E19 m^{-3}
       !interpolate file data to GENE grid
       call lag3interp(temp_file,xaxis_file,nx0_file,&
            sp%temp_prof,xval_a,px0)
       tmp_arr = -log(temp_file)
       call lag3deriv(tmp_arr,xaxis_file,nx0_file,&
            sp%omt_prof,xval_a,px0)
       call lag3interp(dens_file,xaxis_file,nx0_file,&
            sp%dens_prof,xval_a,px0)
       tmp_arr = -log(dens_file)
       call lag3deriv(tmp_arr,xaxis_file,nx0_file,&
            sp%omn_prof,xval_a,px0)
       deallocate(xaxis_file,temp_file,dens_file,tmp_arr,omegator_file)

    case(-4)
#ifdef WITHFUTILS
       call get_entrysize_chease(nx0_file)
       allocate(xaxis_file(1:nx0_file),te_prof(1:nx0_file),&
            &ne_prof(1:nx0_file),ti_prof(1:nx0_file),ni_prof(1:nx0_file),&
            &temp_file(1:nx0_file),dens_file(1:nx0_file),&
            &tmp_arr(1:nx0_file),omegator_file(1:nx0_file))       
       call read_profiles_chease(nx0_file,xaxis_file,te_prof,ti_prof,&
            &ne_prof,ni_prof)
       if (sp%charge.lt.0.0) then !electrons
          temp_file=te_prof*1.e-3 !convert to GENE convention (keV)
          dens_file=ne_prof*1e-19 
       else
          temp_file=ti_prof*1.e-3 !convert to GENE convention (keV)
          dens_file=ni_prof*1e-19
       endif
       !interpolate file data to GENE grid
       call lag3interp(temp_file,xaxis_file/maxval(xaxis_file),nx0_file,&
            sp%temp_prof,xval_a,px0)
       tmp_arr = -log(temp_file)
       call lag3deriv(tmp_arr,xaxis_file/maxval(xaxis_file),nx0_file,&
            sp%omt_prof,xval_a,px0)
       call lag3interp(dens_file,xaxis_file/maxval(xaxis_file),nx0_file,&
            sp%dens_prof,xval_a,px0)
       tmp_arr = -log(dens_file)
       call lag3deriv(tmp_arr,xaxis_file/maxval(xaxis_file),nx0_file,&
            sp%omn_prof,xval_a,px0)
       deallocate(xaxis_file,temp_file,dens_file,tmp_arr)
#else 
       print*, "You have to compile with FUITLS switch set to yes to use CHEASE profiles"  
       stop
#endif
    case default
       stop 'Error: wrong profile type'
    end select


    !assume input profiles to be given in coordinates of x/a
    sp%omt_prof = sp%omt_prof/minor_r
    sp%omn_prof = sp%omn_prof/minor_r

    if (px0.eq.1) then
       i1 = LBOUND(sp%temp_prof)
       sp%temp = sp%temp_prof(i1(1))
       sp%dens = sp%dens_prof(i1(1))
       sp%temp_prof = 1.0
       sp%dens_prof = 1.0
       sp%omt = sp%omt_prof(i1(1))
       sp%omn = sp%omn_prof(i1(1))
    elseif (lilo) then
       call lag3interp(sp%temp_prof,xval_a,nx0,sp%temp,x0)
       call lag3interp(sp%dens_prof,xval_a,nx0,sp%dens,x0)
       call lag3interp(sp%omt_prof,xval_a,nx0,sp%omt,x0)
       call lag3interp(sp%omn_prof,xval_a,nx0,sp%omn,x0)
    else
       sp%temp=1.0
       sp%dens=1.0
    endif

    if (allocated(omegator).and.(omegatorref.eq.-1).or.&
         (omegatorref.eq.-1111.0)) then
       omegatorref = omegator((pi1gl+pi2gl)/2)
       if (write_pe) Write(*,'(A,E13.6)') 'omegatorref set to    ',omegatorref
    endif

  end subroutine read_numerical_profiles


  SUBROUTINE set_ref_vals(rho,nrho,spec,nspec)
    INTEGER, intent(in) :: nrho, nspec!<number of profile grid points/species

    real, dimension(1:nrho) :: rho    !< radial grid for profiles
    real, dimension(1) :: ref_pos, dummy
    real, dimension(1:nrho) :: dummy_arr
    type(spectype) :: spec(0:nspec-1)
    INTEGER :: n, i, nt_norm_index
    logical :: write_pe
    write_pe = (mype.eq.0).and.print_ini_msg

    ref_pos = x0

    if (norm_index.eq.-1) then
       nt_norm_index = 0 !normalize to first species if 
                         !electrons are not found
       
       do n=0, nspec-1
          if (spec(n)%charge.lt.0) nt_norm_index = n
       enddo
    else
       nt_norm_index = norm_index
    endif

    if (spec(0)%prof_type > 0 .and. x_local) then
       !if parameter file is not re-read (after perf_opt), use the stored values
       if (in_Tref==0) in_Tref=Tref
       if (in_nref==0) in_nref=nref
       Tref=in_Tref
       nref=in_nref
       do n=0, n_spec-1
          i=lbound(spec(n)%temp_prof,1)
          if (in_temp(n)==-10) in_temp(n)=spec(n)%temp
          if (in_dens(n)==-10) in_dens(n)=spec(n)%dens
          spec(n)%temp = in_temp(n)
          spec(n)%dens = in_dens(n)
          !convert to dimensional qprofiles evaluated at x0
          !(Tref and nref are the values at the (LT,Ln)center positions )
          !(temp(n) can be used to adjust temperature ratio, dens(n) should fulfill quasineutrality)
          if (Tref>0) spec(n)%temp = spec(n)%temp_prof(i)*spec(n)%temp*Tref
          if (nref>0) spec(n)%dens = spec(n)%dens_prof(i)*spec(n)%dens*nref
          !gradients are dimensionless
          spec(n)%omn = spec(n)%omn_prof(i)
          spec(n)%omt = spec(n)%omt_prof(i)
          !move profile information to reference values in the local code
          spec(n)%temp_prof = 1
          spec(n)%dens_prof = 1
       enddo
       Tref=-1  !recalibrate T just as for prof_from_file
       nref=-1  !recalibrate n just as for prof_from_file
    endif

    if (Tref.lt.0.0) then
       if (x_local) then
          Tref = spec(nt_norm_index)%temp
       else
          !Tref = spec(nt_norm_index)%temp_prof((pi1gl+pi2gl)/2)
          call lag3interp(spec(nt_norm_index)%temp_prof,rho,nrho,&
               dummy,ref_pos,1)
          Tref = dummy(1)
       endif
       if (write_pe) Write(*,'(A,E13.6)') 'Tref/keV set to       ',Tref
    endif

    if (nref.lt.0.0) then
       if (x_local) then
          nref = spec(nt_norm_index)%dens
       else
          !nref = spec(nt_norm_index)%dens_prof((pi1gl+pi2gl)/2)
          call lag3interp(spec(nt_norm_index)%dens_prof,rho,nrho,&
               dummy,ref_pos,1)
          nref = dummy(1)
       endif
       if (write_pe)  Write(*,'(A,E13.6)') 'nref/1E19m^-3 set to ',nref
    endif

    if (mref.lt.0.0) then
       if (norm_index.ge.0) then
          mref = spec(norm_index)%mass
          do n=0, n_spec-1
             spec(n)%mass = spec(n)%mass/mref
          enddo
       else
          mref = m_deuteron/m_proton
       endif
       if (mype==0)  Write(*,'(A,E13.6)') 'mref/mproton set to  ',mref
    endif

    if (x_local) then
       do n=0, nspec-1
          if (nref.gt.0) spec(n)%dens = spec(n)%dens/nref
          if (Tref.gt.0) spec(n)%temp = spec(n)%temp/Tref
       enddo
    else
       do n=0, nspec-1
!          spec(n)%temp=spec(n)%temp_prof((pi1gl+pi2gl)/2)/&
          !   &spec(nt_norm_index)%temp_prof((pi1gl+pi2gl)/2)
!          spec(n)%dens=spec(n)%dens_prof((pi1gl+pi2gl)/2)/&
          !   &spec(nt_norm_index)%dens_prof((pi1gl+pi2gl)/2)
          dummy_arr = (spec(n)%temp_prof/spec(nt_norm_index)%temp_prof)
          call lag3interp(dummy_arr,rho,nrho,dummy,ref_pos,1)
          spec(n)%temp=dummy(1)
          dummy_arr = (spec(n)%dens_prof/spec(nt_norm_index)%dens_prof)
          call lag3interp(dummy_arr,rho,nrho,dummy,ref_pos,1)
          spec(n)%dens=dummy(1)
       enddo
       do n=0, nspec-1
          if (nref.gt.0) spec(n)%dens_prof = spec(n)%dens_prof/spec(n)%dens/nref
          if (Tref.gt.0) spec(n)%temp_prof = spec(n)%temp_prof/spec(n)%temp/Tref
          !for lilo, we enforce constant profiles and gradient profiles
          if (lilo) then
             spec(n)%dens_prof=1.0
             spec(n)%temp_prof=1.0
             dummy_arr = spec(n)%omn_prof
             call lag3interp(dummy_arr,rho,nrho,dummy,ref_pos,1)
             spec(n)%omn_prof=dummy(1)
             dummy_arr = spec(n)%omt_prof
             call lag3interp(dummy_arr,rho,nrho,dummy,ref_pos,1)
             spec(n)%omt_prof=dummy(1)
          endif
       enddo
    endif

  END SUBROUTINE set_ref_vals

  SUBROUTINE calc_derived_ref_vals(beta, coll, debye2, &
       & rhostar, Omega0_tor, write_pe)
    real, intent(inout) :: beta, coll, debye2, rhostar, Omega0_tor
    logical,intent(in):: write_pe

    if (spec(0)%prof_type > 0 .and. x_local) then
       !if parameter file is not re-read (after perf_opt), use the stored values
       if (in_coll==-10) in_coll = coll
       if (in_beta==-10) in_beta = beta
       if (in_debye2==-10) in_debye2 = debye2
       if (in_rhostar==-10) in_rhostar = rhostar
       coll = in_coll
       beta = in_beta
       debye2 = in_debye2
       rhostar = in_rhostar
    endif

    !Tref in keV, Lref in m, Bref in T, 
    !nref in 1E19 1/m^3, mref in m_proton
    if ((Tref.gt.0).and.(nref.gt.0).and.(mref.gt.0).and.&
         &(Lref.gt.0).and.(Bref.gt.0)) then
       cref = 9.7871518E3*sqrt((Tref*1E3)/mref) !see, e.g., NRL

       if (beta.lt.0.0.and.abs(beta).gt.epsilon(beta)) then
            beta = 403.0E-5*nref*Tref/(Bref*Bref)
            if (write_pe) &
                 Write(*,"(A,1ES13.4,A)") "computed beta    = ",beta,"  from profile"
       endif
       
       if (coll.lt.0.0.and.abs(coll).gt.epsilon(coll)) then
            coll = 2.3031E-5*Lref*(nref)/(Tref)**2*&
            (24.-log(sqrt(nref*1.0E13)/Tref*0.001))
            if (write_pe) &
                 Write(*,"(A,1ES13.4,A)") "computed coll    = ",coll,"  from profile"
       endif
       
       if (debye2.lt.0.0.and.abs(debye2).gt.epsilon(debye2)) then
            debye2 = 5.2936E-4 * (Bref)**2/(nref) * (1./mref)
            if (write_pe) &
                 Write(*,"(A,1ES13.4,A)") "computed debye2  = ",debye2,"  from profile"
       endif
       
       if (rhostar.lt.0.0.and.abs(rhostar).gt.epsilon(rhostar)) then
            rhostar = 3.2255E-3 * sqrt((mref)*Tref)/Bref/(minor_r*Lref)
            if (write_pe) &
                 Write(*,"(A,1ES13.4,A)") "computed rhostar = ",rhostar,"  from profile"
       endif
       
       if (Omega0_tor.eq.-1111.0) then
          Omega0_tor = omegatorref * Lref / cref
          if (write_pe) Write(*,"(A,1ES13.4,A)") "computed Omega0_tor = ", &
               Omega0_tor,"  from profile"
       endif

    else
       if ((beta.lt.0).or.(coll.lt.0).or.(debye2.lt.0).or.&
            &(Omega0_tor.eq.-1111.0)) &
            stop 'Missing reference values to compute derived input parameters'
    endif

  END SUBROUTINE calc_derived_ref_vals




  subroutine enforce_quasineutrality(spec,write_pe)
    TYPE(spectype) :: spec(0:n_spec-1)
    logical :: write_pe
    integer :: n, iind
    real, dimension(pi1gl:pi2gl) :: qn, qnomn

    if (n_spec.eq.1) return

    !find last ion species and compute chargedens
    qn = 0.0
    qnomn = 0.0
    iind = 0
    do n=0, n_spec-1
       if (spec(n)%passive) cycle
       if (spec(n)%charge.gt.0) iind=n
       qn = qn+spec(n)%charge*&
            &spec(n)%dens*spec(n)%dens_prof
       if (x_local) then
          qnomn = qnomn+spec(n)%charge*&
               &spec(n)%dens*spec(n)%omn
       else
          qnomn = qnomn+spec(n)%charge*&
            &spec(n)%dens*spec(n)%dens_prof*&
            &spec(n)%omn_prof
       endif
    enddo

    if (abs(sum(qn)).gt.epsilon(qn(pi1gl))) then
       if (write_pe) then
          write (*,'(A,ES8.1)') 'avg. quasineutrality violation in profiles: ' ,&
               sum(qn)/(pi2gl-pi1gl+1)
          write (*,'(A)') 'enforcing quasineutrality on last active ion species'
       endif

       !n_j = 1/q_j * sum(s<>j) q_s n_s
       qn = qn-spec(iind)%charge*spec(iind)%dens*spec(iind)%dens_prof
       if (x_local) then
          qnomn = qnomn-spec(iind)%charge*spec(iind)%dens*spec(iind)%omn
          spec(iind)%dens = -qn(pi1gl)/real(spec(iind)%charge)
          spec(iind)%omn = - qnomn(pi1gl)/(spec(iind)%dens*spec(iind)%charge)
       else
          qnomn = qnomn-spec(iind)%charge*spec(iind)%dens*spec(iind)%dens_prof*&
            &spec(iind)%omn_prof
          !assuming spec%dens=1 which should still hold here
          spec(iind)%dens_prof = -qn/real(spec(iind)%charge)
          spec(iind)%omn_prof = -qnomn/(spec(iind)%dens_prof*spec(iind)%charge)
       endif
    endif
  end subroutine enforce_quasineutrality


  !>calculate the negative pressure gradient from n and T gradients 
  !!and normalize to reference magnetic pressure
  Subroutine calc_dpdx_pm(dpdx_pm_arr)
    Real, Dimension(pi1gl:pi2gl),intent(out) :: dpdx_pm_arr
    Integer :: i, n
    
    do i=pi1,pi2
       dpdx_pm_arr(i) = 0.0
       do n=0, n_spec-1
         IF (.NOT. spec(n)%passive) &
           dpdx_pm_arr(i) = dpdx_pm_arr(i)+(spec(n)%omn_prof(i)+spec(n)%omt_prof(i))*&
               & spec(n)%temp*spec(n)%temp_prof(i)*&
               & spec(n)%dens*spec(n)%dens_prof(i)
       enddo
    enddo
    
    !switch normalization from pref to reference magnetic pressure
    dpdx_pm_arr = dpdx_pm_arr * beta

    if (n_spec.eq.1) then
       if ( (mype.eq.0).and.print_ini_msg ) WRITE(*,'(A)') &
            &'assuming same pressure gradient for adiabatic species in dpdx_pm computation'
       dpdx_pm_arr = dpdx_pm_arr * 2.0
    endif

    if (x_local.and.mype.eq.0) then !.and.print_ini_msg) then
       write(*,"(A,F12.6)") 'dpdx_pm taken from grad n/T: ', dpdx_pm_arr(pi1gl)
    endif

#if 0
    if (mype==0) then
       write(*,*)
       write(*,'(A)') 'from profiles:'
       write(*,'(A)') ' x0           pres         beta_p       dpdx'
       write(*,'(4(ES12.5,1X))') x0, nref*Tref*1.6022*1000*SUM(spec(:)%dens*spec(:)%temp), &
            &beta*SUM(spec(:)%dens*spec(:)%temp), dpdx_pm_arr(pi1gl)
       write(*,*)   
    endif
#endif
  End Subroutine calc_dpdx_pm

END MODULE profiles_mod
  
