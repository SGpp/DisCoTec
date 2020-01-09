#include "redef.h"
#include "intrinsic_sizes.h"
  
module quasilinear_model
  use discretization
  use coordinates
  use communications
  use par_in
  use par_other
  use geometry
  use aux_fields
  use time_averages, only: avgfluxes, avgflux_stime

  implicit none

  public :: set_quasilin,quasilin_model, ql_diffusivity, quasilin_flux

  private

  integer :: quasilin_model=0 !0=off,1=simple,>1 not implemented yet
  real :: ql_diffusivity = 0.0 
  real, dimension(:,:), allocatable :: quasilin_flux

Contains
  
  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_quasilin(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    mem_est_quasilin=mem_req_in+mem_loc
  end Function mem_est_quasilin


  !>Check parameters for quasilinear model
  subroutine check_par_quasilin
    if (quasilin_model.gt.0) then
!       if (istep_nrg.le.0) istep_nrg = 10
       if (avgflux_stime.lt.0.0) avgflux_stime = 0.0
    endif
  end subroutine check_par_quasilin


  subroutine set_quasilin(gamma)
    real, intent(in) :: gamma

    allocate(quasilin_flux(0:n_spec-1,2))

    if (quasilin_model.eq.1) &
         &call set_simple_quasilin(gamma,quasilin_flux)

  end subroutine set_quasilin

  
  !>A simple quasilinear model being based on the linear 
  !!growth rate devided by the normalized <k_perp^2 * |phi|^2>
  !!
  subroutine set_simple_quasilin(gamma,quasilin_flux)
    real, intent(in) :: gamma
    real, dimension(0:n_spec-1,2), intent(out) :: quasilin_flux
    !local variables
    real :: kperp2_av, normfield, max_flux, normflux
    integer :: i,j,n,o,normspec

    kperp2_av = 0.0
    normfield = 0.0

    if (.not.xy_local) stop 'nonlocal quasilinear model not implemented yet!'

    if (yx_order) then
       do j=lj1,lj2
          do i=li1,li2
             if ((evenx.eq.1).and.(i.eq.hkx+1)) cycle
             kperp2_av = kperp2_av + &
                  & sum((geom%gii(pi1,pj1,:)*ki(i)**2+2.0*geom%gij(pi1,pj1,:)*&
                  &ki(i)*kx_center+geom%gjj(pi1,pj1,:)*kx_center**2)*&
                  &geom%jacobian(pi1,pj1,:)*abs(emfields(i,j,lk1:lk2,1))**2)
             normfield = normfield + &
                  sum(geom%jacobian(pi1,pj1,:)*abs(emfields(i,j,lk1:lk2,1))**2)
          enddo
       enddo
    else
       do j=lj1,lj2
          do i=li1,li2
             if ((evenx.eq.1).and.(i.eq.hkx+1)) cycle
             kperp2_av = kperp2_av + &
                  & sum((geom%gii(pi1,pj1,:)*kx_center**2+2.0*geom%gij(pi1,pj1,:)*&
                  &kx_center*kj(j)+geom%gjj(pi1,pj1,:)*kj(j)**2)*&
                  &geom%jacobian(pi1,pj1,:)*abs(emfields(i,j,lk1:lk2,1))**2)
             normfield = normfield + &
                  sum(geom%jacobian(pi1,pj1,:)*abs(emfields(i,j,lk1:lk2,1))**2)
          enddo
       enddo
    endif

    call my_sum_to_all_real_0d(kperp2_av,mpi_comm_z)
    call my_sum_to_all_real_0d(normfield,mpi_comm_z)

    kperp2_av = kperp2_av/normfield

    ql_diffusivity = gamma/kperp2_av

    normspec = maxloc(avgfluxes(:,2),1)-1 !-1 to correct maxloc index
    max_flux = maxval(avgfluxes(:,2),1)
    
    if ((max_flux.gt.0).and.(ql_diffusivity.gt.0)) then
       normflux = spec(normspec)%dens*spec(normspec)%temp*&
            &ql_diffusivity*spec(normspec)%omt/max_flux
    else
       normflux = 0.0
    endif

    do o=1,2
       do n=0, n_spec-1
          quasilin_flux(n,o) = avgfluxes(n,o)*normflux
       enddo
    enddo

    if (mype.eq.0) then
       write(*,"(a)") 'Results of quasilinear transport model'
       write(*,"(a,f6.3,a,f6.3,a,f10.2)") 'ky= ',ky(lj1),'  <k_perp^2>= ',kperp2_av,&
                  '  gamma/<k_perp>^2= ',ql_diffusivity
       do n=0,n_spec-1
          write(*,"(2a,f6.3,a,f6.3,a,f10.2)") spec(n)%name,': G= ',&
               &quasilin_flux(n,1),', Q= ',quasilin_flux(n,2)
       enddo
       write(*,*)
    endif
    
  end subroutine set_simple_quasilin

  subroutine get_quasilin_flux(quasilin_flux_out)
    real, dimension(0:n_spec-1,2), intent(out) :: quasilin_flux_out

    quasilin_flux_out = quasilin_flux
    deallocate(quasilin_flux)

  end subroutine get_quasilin_flux

end module quasilinear_model
