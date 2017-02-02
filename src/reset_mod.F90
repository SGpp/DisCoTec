#include "redef.h"
!>Module containing subroutines to reset the distribution function
!!and all profile information throughout the whole GENE code
Module reset_mod
  Use par_mod
!  Use parameters_IO, only: write_parameters
  Use spectype_mod
  Use vel_space
  use prefactors
  use geometry, only: geom, rhostar,minor_r, write_geometry
  use fieldsolver_df
  use calc_rhs
  use diagnostics
  use diagnostics_df
  use coordinates, only: calc_center_value
  use profiles_mod, only: write_tempdens_profiles
  use f0_term, only: f0_contr
  use aux_fields, only: initialize_calc_aux_fields, &
       &finalize_calc_aux_fields
  use lagrange_interpolation, only: lag3deriv
  use profile_smoothing

  Implicit None

  !The basic reset interface
  PUBLIC :: reset_g1_prof, set_reset_mod_defaults

  PRIVATE
  
  integer :: reset_count = 0

Contains

  subroutine set_reset_mod_defaults
    call set_profile_smoothing_defaults
  end subroutine set_reset_mod_defaults

  Subroutine reset_g1_prof(g_1,new_profs)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(inout):: g_1
    real,dimension(0:nx0-1,0:n_spec-1,1:4),intent(inout) :: new_profs
    real,dimension(0:nx0-1,0:n_spec-1,1:4) :: smooth_new_profs
    integer :: j,k,l,m,n,dummy_int
    Real, Dimension(pi1:pi2) :: fm_new
    Logical :: dummy1, dummy2,dummy3
!    Real :: T_center, n_center
    real, dimension(pi1gl:pi2gl):: dens_electrons, omn_electrons

    !for now: reset equilibrium densities for all species to the 
    !electron density profile (times spec(n)%dens)
    do n=0,n_spec-1
       if (spec(n)%charge==-1) then
          dens_electrons=new_profs(:,n,1)
          omn_electrons=new_profs(:,n,3)
       endif
    enddo

    do n=0,n_spec-1
       if (n_spec.gt.1) then
          new_profs(:,n,1)=dens_electrons
          new_profs(:,n,3)=omn_electrons
       endif
       if (smooth_reset) then
          call smooth_profiles(new_profs(:,n,:), smooth_new_profs(:,n,:),smooth_reset_width)
       else
          smooth_new_profs(:,n,:) = new_profs(:,n,:)
       endif
       !If we only have one ion species and adiabatic electrons, leave density
       !alone
       if (n_spec .eq. 1) then
          smooth_new_profs(:,n,1) = spec(n)%dens_prof
          smooth_new_profs(:,n,3) = spec(n)%omn_prof
       endif
       !normalize profile to center value
       !call calc_center_value(smooth_new_profs(:,n,1),n_center)
       spec(n)%dens_prof=smooth_new_profs(:,n,1) !/n_center
       !call calc_center_value(smooth_new_profs(:,n,2),T_center)
       spec(n)%temp_prof=smooth_new_profs(:,n,2) !/T_center
       spec(n)%omn_prof=smooth_new_profs(:,n,3)
       spec(n)%omt_prof=smooth_new_profs(:,n,4)
    enddo

    !for now: g_1^new = g_1^old - Delta f_0 = (f_1^old-(f_0^old-f_0^new)) + q/c Aparbar^old vpar/T0 f_0^old
    !however, we might need to fix the second term as well in the future
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                fm_new = spec(n)%dens_prof(pi1:pi2) *&
                     (pi*spec(n)%temp_prof(pi1:pi2))**(-1.5) *&
                     & exp(-(mu(m)*geom%Bfield(pi1:pi2,pj1,k)+vp(l)**2)/&
                     spec(n)%temp_prof(pi1:pi2))
!                do j=lj1,lj2
                if (p_has_0_mode) then
                   j=lj1
                   g_1(:,j,k,l,m,n) = g_1(:,j,k,l,m,n)+(fm(pi1:pi2,pj1,k,l,m,n)&
                        &-fm_new(pi1:pi2))/(rhostar*minor_r)
                endif
                fm(pi1:pi2,pj1,k,l,m,n)=fm_new(pi1:pi2)
!                enddo
             enddo
          enddo
       enddo
    enddo

    call finalize_all_diags
    call finalize_CalFullRhs
    call finalize_calc_aux_fields

    if (mype==0) print*, 'old file_extension= ',file_extension
    read (file_extension(2:),'(I19)') dummy_int
    dummy_int = dummy_int+1
    write(file_extension,"(A,I4.4)") '_',dummy_int
    if (mype==0) print*, 'new file_extension= ', file_extension
    itime = 0
    
    call set_prefactors
    call initialize_calc_aux_fields
    call initialize_CalFullRhs
    call initialize_all_diags

    !adapt pressure gradient? (would require additional logs whether
    !pressure gradient originated from MHD file or from gradients)

    if (mype==0) then
       print*, 'Resetted g_1 and all n and T dep. profiles'
       do n=0,n_spec-1
          Call write_tempdens_profiles(spec(n),.false.)
       enddo
    endif

    !write parameters & geometry
!    If (mype == 0) Call write_parameters !Introduces circular reference  :-(
    if(write_geom) call write_geometry

    !write current time step
    call exec_all_diags(0,time,dummy1,dummy2,dummy3)

  End Subroutine reset_g1_prof

     
End Module reset_mod
