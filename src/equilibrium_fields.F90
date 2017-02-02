#include "redef.h"
#include "intrinsic_sizes.h"
!>Contains and computes equilibrium field
!!currently, only present/relevant in co-moving frames

Module equilibrium_fields
  use par_in, only: spec, Omega0_tor
  use discretization
  use geometry, only: geom, major_R, C_y
  use external_contr, only: ExBrate
  use communications, only: my_real_bcast


  implicit none

  public:: initialize_equilibrium_fields, finalize_equilibrium_fields, &
       & mem_est_equilibrium_fields, with_comoving_other, phi0, dens_co, &
       & set_equilibrium_fields_defaults, check_par_equilibrium_fields,&
       & dphi0dx_part1, dphi0_fac, R0_tor, dR0_tor_dx

  private

  !rotation/with_comoving_other frame 
  logical :: with_comoving_other
  real, dimension(:,:), allocatable :: phi0, dphi0dx_part1, dphi0_fac
  real:: R0_tor, dR0_tor_dx

Contains
  
  Real Function mem_est_equilibrium_fields(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    !phi0
    if (with_comoving_other) mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*lk0
    
    mem_est_equilibrium_fields=mem_req_in+mem_loc
  end Function mem_est_equilibrium_fields

  subroutine set_equilibrium_fields_defaults

    with_comoving_other=.true.
    R0_tor=-1
    dR0_tor_dx=-1

  end subroutine set_equilibrium_fields_defaults

  subroutine check_par_equilibrium_fields
    !with finite toroidal angular velocity, go to co-moving frame if
    !not explicitly set to false by user
    with_comoving_other = with_comoving_other.and.&
         ((abs(Omega0_tor).ge.epsilon(Omega0_tor)).and.x_local)
  end subroutine check_par_equilibrium_fields

  subroutine initialize_equilibrium_fields

    if (with_comoving_other) then

       !set reference position for background density:
       !R0_tor=<value>, dR0_tor_dx=<value>: arbitrary settings
       !R0_tor=-1, dR0_tor_dx set by GENE: use magnetic axis (default!)
       !R0_tor=-2, dR0_tor_dx set by GENE: use low field side R (useful for 
       !   benchmarks, e.g. to Casson PoP 17, 102305 (2010))

       if (R0_tor.eq.-1) then
          R0_tor=major_R
          dR0_tor_dx=0.
       else if (R0_tor.eq.-2.) then
          if (nz0/2.ge.lk1.and.nz0/2.le.lk2) then
             R0_tor=geom%R_hat(pi1,nz0/2)
             dR0_tor_dx=1./geom%dxdR(pi1,nz0/2)
          endif
          call my_real_bcast(R0_tor,nz0/2/lk0)
          call my_real_bcast(dR0_tor_dx,nz0/2/lk0)
       else if (R0_tor.gt.0.) then
          if (dR0_tor_dx.eq.-1.) then
             if (mype==0) print *, 'WARNING: dR0_tor_dx unspecified, setting to 0'
             dR0_tor_dx=0.
          endif
       else
          stop 'negative R0_tor makes no sense!'
       endif

       
       allocate(phi0(pi1:pi2,lk1:lk2))
       call phi0_solve_ff(phi0)
       allocate(dphi0dx_part1(pi1:pi2,lk1:lk2))
       allocate(dphi0_fac(pi1:pi2,lk1:lk2))
       call set_dphi0_parts
    end if

  end subroutine initialize_equilibrium_fields

  subroutine finalize_equilibrium_fields
    if (with_comoving_other) then
       deallocate(phi0)
       deallocate(dphi0dx_part1,dphi0_fac)
    endif
  end subroutine finalize_equilibrium_fields

  !>Solve quasineutrality equation for the equilibrium potential
  !!in co-moving frame with simple newton method
  Subroutine phi0_solve_ff(phi0)
    Real, Dimension(pi1:pi2,lk1:lk2) :: phi0
    Real :: qn, qn_prime, tmp, phi
    integer :: i, k, n, it, maxit=10000

    tmp = 1.0
    do n=0, n_spec-1
       if (spec(n)%passive) cycle
       tmp = tmp*spec(n)%charge
    enddo
    if ((n_spec.lt.2).or.(tmp.gt.0)) &
         & stop 'phi0_solve_ff: adiabatic species not implemented yet'

    phi = 0.1
    do k=lk1,lk2
       do i=pi1,pi2
          do it = 0, maxit
             qn = 0.0
             qn_prime = 0.0
             do n = 0, n_spec-1
                if (spec(n)%passive) cycle
                tmp = spec(n)%charge * spec(n)%dens * &
                     & EXP((-spec(n)%charge*phi + &
                     & 0.5*spec(n)%mass*Omega0_tor**2*&
                     & (geom%R_hat(i,k)**2-R0_tor**2))/spec(n)%temp)
                qn = qn + tmp
                qn_prime = qn_prime - spec(n)%charge/spec(n)%temp*tmp
             enddo

             if (abs(qn)<10.*epsilon(phi)) exit
             
             phi=phi-qn/qn_prime

          end do
          if (it.ge.maxit) stop 'reached max. iterations in phi0 computation'

          phi0(i,k) = phi

!benchmark with analytic 2-species result
!          phi = 0.5*spec(0)%temp*spec(n_spec-1)%temp/&
!               (spec(0)%temp+spec(n_spec-1)%temp)*&
!               (spec(0)%mass/spec(0)%temp-spec(n_spec-1)%mass/&
!               spec(n_spec-1)%temp)*Omega0_tor**2*&
!               (geom%R_hat(i,k)**2-R0_tor**2)
!          if (mype.eq.0) print*, k, phi0(i,k), phi

       enddo
    enddo

  End Subroutine phi0_solve_ff

  !>Returns *dens*ity correction in *co*moving frame
  real function dens_co(i,k,n)
    integer, intent(in) :: i,k,n

    if (with_comoving_other) then
       dens_co = exp((-spec(n)%charge*phi0(i,k)+&
            &0.5*spec(n)%mass*Omega0_tor**2*&
            &(geom%R_hat(i,k)**2-R0_tor**2))/spec(n)%temp)
    else
       dens_co = 1.0
    endif

  end function dens_co

  subroutine set_dphi0_parts
    real :: norm, part1, part2_fac
    integer :: i,k,n

    do k=lk1,lk2
       do i=pi1,pi2
          norm = 0.0
          part1 = 0.0
          part2_fac = 0.0
          do n=0, n_spec-1
             if (spec(n)%passive) cycle
             norm = norm + spec(n)%charge**2/&
                  &spec(n)%temp*spec(n)%dens*dens_co(i,k,n)
             part1 = part1 + spec(n)%charge*spec(n)%dens*dens_co(i,k,n)*&
                  &(-spec(n)%omn+spec(n)%omt*(-spec(n)%charge*phi0(i,k)+&
                  &0.5*spec(n)%mass*Omega0_tor**2*(geom%R_hat(i,k)**2-&
                  &R0_tor**2))/spec(n)%temp-spec(n)%mass/spec(n)%temp*&
                  &Omega0_tor*ExBrate/C_y(i)*(geom%R_hat(i,k)**2-R0_tor**2)-&
                  spec(n)%mass/spec(n)%temp*Omega0_tor**2*(R0_tor*dR0_tor_dx))
             part2_fac = part2_fac + spec(n)%charge*spec(n)%dens*&
                  &dens_co(i,k,n)*spec(n)%mass/spec(n)%temp
          enddo
          dphi0dx_part1(i,k) = part1/norm
          dphi0_fac(i,k) = part2_fac/norm
       enddo
    enddo

  end subroutine set_dphi0_parts

end Module equilibrium_fields
