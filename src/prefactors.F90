#include "redef.h"
#include "intrinsic_sizes.h"
!>Provides prefactors for, e.g., the Vlasov eq.
module prefactors
  Use par_mod
  use discretization
  use coordinates
  Use external_contr
  Use geometry
  use communications
  use vel_space
  use equilibrium_fields, only: with_comoving_other,phi0, &
       &dphi0dx_part1,dphi0_fac, R0_tor, dR0_tor_dx

  implicit none

  public:: initialize_prefactors, set_prefactors, finalize_prefactors, mem_est_prefactors
  public:: vTvpar, pdchibardi, pdchibardj, edr_for_energy,&
       pdg1di, pdg1dj, mu_tjqj
  public:: edr, curv, B_Bstar
  public:: with_coriolis, with_centrifugal, with_bxphi0
  public:: coriolis_i, centrifugal_i, bxphi0_i
  public:: set_prefactors_defaults, check_par_prefactors 

  private

!--------------------------------------------------------------------------
!   definition of prefactors 
!--------------------------------------------------------------------------

  !derived prefactors
  Real, Dimension(:,:), Allocatable :: vTvpar

  Real, Dimension(:,:,:,:,:,:), Allocatable:: pdchibardi, pdchibardj,edr_for_energy
  ! edr_for_energy2
  Real, Dimension(:,:,:,:,:,:), Allocatable:: pdg1di, pdg1dj

  ! variables for B_parallel
  REAL, DIMENSION(:,:), ALLOCATABLE :: mu_tjqj

  !switches for rotation effects
  logical :: with_coriolis, with_centrifugal, with_bxphi0

Contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_prefactors(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0
    
    !quantities from allocate_prefactors

    !pdg1di, pdg1dj, pdchibardi, pdchibardj
    mem_loc=4*SIZE_OF_REAL_MB*pi0*pj0*lklmn0

    !edr_for_energy
    if ((istep_energy.gt.0).or.diag_GyroLES.or.&
         (istep_energy3d.gt.0)) &
         &mem_loc = mem_loc + SIZE_OF_REAL_MB*pi0*pj0*lklmn0

    mem_est_prefactors=mem_req_in+mem_loc
  End Function mem_est_prefactors 

  Subroutine set_prefactors_defaults
    with_coriolis = .true.
    with_centrifugal = .true.
    with_bxphi0 = .true.
  end Subroutine set_prefactors_defaults

  Subroutine check_par_prefactors
    !with finite toroidal angular velocity, switch on rotation terms if
    !not explicitly deactivated by user
    with_coriolis = with_coriolis.and.&
         ((abs(Omega0_tor).ge.epsilon(Omega0_tor)).and.x_local)
    with_centrifugal = with_centrifugal.and.&
         ((abs(Omega0_tor).ge.epsilon(Omega0_tor)).and.x_local)
    with_bxphi0 = with_bxphi0.and.with_comoving_other
!         ((abs(Omega0_tor).ge.epsilon(Omega0_tor)).and.x_local)
  end Subroutine check_par_prefactors
  
  !>Allocates prefactors, e.g., for Vlasov eq.  
  Subroutine allocate_prefactors

    allocate(vTvpar(ll1:ll2, ln1:ln2))
    if (n_fields.gt.2) allocate(mu_tjqj(lm1:lm2,ln1:ln2))

    Allocate(pdg1di(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(pdg1dj(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(pdchibardi(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(pdchibardj(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    if (istep_energy.gt.0.or.istep_energy3d.gt.0.or.diag_GyroLES) &
         Allocate(edr_for_energy(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    
    equil_par_curr = .NOT. ((beta .EQ. 0.0) .OR. (ALL(currdens_par .LT. 1e-5)))

  end Subroutine allocate_prefactors


  !>Initializes prefactors, e.g., for Vlasov eq.
  !!\param vp v_parallel grid (explicitly passed to avoid cross reference gyro_average <-> prefactors <-> vel_space)
  !!\param mu mu grid (explicitly passed to avoid cross reference gyro_average <-> prefactors <-> vel_space)
  Subroutine initialize_prefactors

    Call allocate_prefactors
    Call set_prefactors

  end Subroutine initialize_prefactors


  !>Sets the prefactors, e.g., for Vlasov eq.
  !!\param vp v_parallel grid (explicitly passed to avoid cross reference gyro_average <-> prefactors <-> vel_space)
  !!\param mu mu grid (explicitly passed to avoid cross reference gyro_average <-> prefactors <-> vel_space)
  Subroutine set_prefactors
    Integer:: j,k, l, m, n, pni
    Real, Dimension(pi1:pi2) :: vabs2
    real, dimension(pi1:pi2) :: B_B_star, curv_, qjTjF0, press, edr_fe_loc, tmp,&
        & coriolis_i_, coriolis_j_, centrifugal_i_, centrifugal_j_,&
        & bxphi0_i_, bxphi0_j_

    press = 0.0D0
    maxclimit = 0.0

    ! ... 
    pni = pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          if (n_fields.gt.2) mu_Tjqj(m,n) = spec(n)%temp/spec(n)%charge * mu(m)
          do l=ll1,ll2
             !v_Tj(x0) * v_||
             vTvpar(l,n) = sqrt(2.*spec(n)%temp/spec(n)%mass)*vp(l)
             do k=lk1,lk2
                do j=pj1,pj2
                   !Squared velocity (temporary)
                   vabs2(:) =&
                        (mu(m)*geom%Bfield(pi1:pi2,j,k)+vp(l)**2)/&
                        spec(n)%temp_prof(pi1:pi2)

                   !calculate B0/B0||* prefactor if needed
                   B_B_star = B_Bstar(j,k,l,n)
                   
                   !calculate pressure term if needed
                   if (pressure_term) then
                      press = (vp(l)/geom%Bfield(pi1:pi2,j,k))**2
                      if (dpdx_term.eq.'gradB_eq_curv') & !full drift velocity minus grad p term in grad B drift
                           & press = press + 0.5*mu(m)/geom%Bfield(pi1:pi2,j,k)
                      press = press*B_B_star*spec(n)%temp/spec(n)%charge*&
                           &dpdx_pm_arr(pi1:pi2)/C_xy(pi1:pi2)
                   endif

                   
                   ! --- locally used prefactors ---
                   ! Prefactor to curvature term
                   curv_ = curv(j,k,l,m,n)

                   ! Coriolis drift terms
                   if (with_coriolis) then
                      coriolis_i_=coriolis_i(j,k,l,n)
                      coriolis_j_=coriolis_j(j,k,l,n)
                   endif

                   ! Centrifugal drift terms
                   if (with_centrifugal) then
                      centrifugal_i_=centrifugal_i(j,k,l,m,n)
                      centrifugal_j_=centrifugal_j(j,k,l,m,n)
                   endif

                   ! Equilibrium ExB drift due to comoving frame
                   if (with_bxphi0) then
                      bxphi0_i_ = bxphi0_i(j,k,l,m,n)
                      bxphi0_j_ = bxphi0_j(j,k,l,m,n)
                   endif

                   ! Another temporary prefactor
                   qjTjF0 = spec(n)%charge/(spec(n)%temp*spec(n)%temp_prof(pi1:pi2))*fm(:,j,k,l,m,pni)

                   ! --- prefactors of parallel dynamics terms ---
                   ! d/dz F_1 prefactor
                   maxclimit(3) = max(maxval(abs(-C_xy(pi1:pi2)*&
                        &vTvpar(l,n)/(geom%jacobian(pi1:pi2,j,k)*&
                        &geom%Bfield(pi1:pi2,j,k)))),maxclimit(3))

                   ! --- prefactors of (gyro averaged) chi derivatives ---
                   pdchibardi(:,j,k,l,m,n) = -qjTjF0*curv_*geom%K_i(:,j,k)
                   pdchibardj(:,j,k,l,m,n) = -qjTjF0*curv_*geom%K_j(:,j,k)

                   if (yx_order) then
                      pdchibardi(:,j,k,l,m,n)= pdchibardi(:,j,k,l,m,n)- edr(j,k,l,m,n)
                   else
                      pdchibardj(:,j,k,l,m,n)= pdchibardj(:,j,k,l,m,n)- edr(j,k,l,m,n)
                   end if

                   if (with_coriolis) then
                      pdchibardi(:,j,k,l,m,n) = pdchibardi(:,j,k,l,m,n) + qjTjF0*coriolis_i_
                      pdchibardj(:,j,k,l,m,n) = pdchibardj(:,j,k,l,m,n) + qjTjF0*coriolis_j_
                   endif

                   if (with_centrifugal) then
                      pdchibardi(:,j,k,l,m,n) = pdchibardi(:,j,k,l,m,n) + qjTjF0*centrifugal_i_
                      pdchibardj(:,j,k,l,m,n) = pdchibardj(:,j,k,l,m,n) + qjTjF0*centrifugal_j_
                   endif

                   if (with_bxphi0) then
                      pdchibardi(:,j,k,l,m,n) = pdchibardi(:,j,k,l,m,n) + qjTjF0*bxphi0_i_
                      pdchibardj(:,j,k,l,m,n) = pdchibardj(:,j,k,l,m,n) + qjTjF0*bxphi0_j_
                   endif                  

                   ! prefactors for energy diagnostics
                   if (istep_energy.gt.0.or.istep_energy3d.gt.0.or.diag_GyroLES) then
                      if (n_spec.eq.1) then
                         ! Prefactor to electrostatic drive term in the energy equation (for energy diagnostics)
                         edr_fe_loc = ( spec(n)%omt_prof(pi1:pi2) * &
                              &vabs2(:) ) * fm(:,j,k,l,m,pni)*B_B_star/C_xy(pi1:pi2)

                         edr_for_energy(:,j,k,l,m,n) = -edr_fe_loc
                      else
                         edr_for_energy(:,j,k,l,m,n) = -edr(j,k,l,m,n)
                      endif
                   endif


                   ! dchibary prefactor with pressure term
                   if(pressure_term) then
                      if(yx_order) then
                         pdchibardi(:,j,k,l,m,n) = pdchibardi(:,j,k,l,m,n) +&
                              press*qjTjF0
                      else
                         pdchibardj(:,j,k,l,m,n) = pdchibardj(:,j,k,l,m,n) +&
                              press*qjTjF0
                      end if
                   end if

                   ! --- prefactors of g1 derivatives ---
                   pdg1di(:,j,k,l,m,n)=-curv_(:)*geom%K_i(:,j,k)
                   pdg1dj(:,j,k,l,m,n)=-curv_(:)*geom%K_j(:,j,k)


                   !additional contributions from the pressure terms 
                   if(pressure_term) then
                      if(yx_order) then
                         pdg1di(:,j,k,l,m,n)  = pdg1di(:,j,k,l,m,n) + press
                      else
                         pdg1dj(:,j,k,l,m,n)  = pdg1dj(:,j,k,l,m,n) + press
                      endif
                   endif

                   IF ((with_phi_ext.or.with_apar_ext).AND.(.NOT.x_local)) THEN
                      ! prefactor containing external fields
                      tmp(pi1:pi2) = B_B_star/C_xy(pi1:pi2)
                      CALL add_phi_ext(pdg1dj(:,j,k,l,m,n),tmp,l,n)
                   END IF

                   if (with_coriolis) then
                      pdg1di(:,j,k,l,m,n)=pdg1di(:,j,k,l,m,n)+coriolis_i_
                      pdg1dj(:,j,k,l,m,n)=pdg1dj(:,j,k,l,m,n)+coriolis_j_
                   endif

                   if (with_centrifugal) then
                      pdg1di(:,j,k,l,m,n)=pdg1di(:,j,k,l,m,n)+centrifugal_i_
                      pdg1dj(:,j,k,l,m,n)=pdg1dj(:,j,k,l,m,n)+centrifugal_j_
                   endif
                   if (with_bxphi0) then
                      pdg1di(:,j,k,l,m,n)=pdg1di(:,j,k,l,m,n)+bxphi0_i_
                      pdg1dj(:,j,k,l,m,n)=pdg1dj(:,j,k,l,m,n)+bxphi0_j_
                   endif
                   ! --- external stuff (only interesting for tertiary instab. analysis) ---
                   ! set prefactor for zonal temperature/density  gradients (if needed)
                   ! currently without pressure contribution !!!
                   ! The reason why the following contributions cannot be combined i.e.
                   ! with edr is the additional convolution which is required due to the 
                   ! x dependence of omt_ext, omn_ext
                   if ((with_omt_ext).and.(xy_local)) &
                        & call build_ext_contr_prefac(1,(-(vabs2(pi1) - 1.5) * &
                        &  fm(pi1,pj1,k,l,m,pni)*B_B_star(pi1)/C_xy(pi1)),k,l,m,n)
                   if ((with_omn_ext).and.(xy_local)) &
                        & call build_ext_contr_prefac(2,(-fm(pi1,pj1,k,l,m,pni)*B_B_star(pi1)/&
                        &   C_xy(pi1)),k,l,m,n)
                   ! --- end external stuff ---
                enddo !j
             end do !k
          end do !l
          do k=lk1,lk2
             do j=pj1,pj2
                !trp term
                maxclimit(4) = max(maxval(abs(C_xy(pi1:pi2)*mu(m)/&
                     &sqrt(2.*spec(n)%mass/spec(n)%temp)*geom%dBdz(pi1:pi2,j,k)/&
                     &(geom%jacobian(pi1:pi2,j,k)*geom%Bfield(pi1:pi2,j,k)))),&
                     &maxclimit(4))
             end do
          end do
       end do
    end do
    !for arakawa_zv: use the same estimate, but with +20% safety factor.
    !because v/z discretization is actually not separated
    if (arakawa_zv) maxclimit(3)=maxclimit(3)*1.2
    if (arakawa_zv) maxclimit(4)=maxclimit(4)*1.2
    
    call maxclimits
    
  End Subroutine set_prefactors

  function B_Bstar(j,k,l,n)
    integer, intent(in) :: j,k,l,n
    real, dimension(pi1:pi2) :: B_Bstar
    
    if (equil_par_curr) then 
       B_Bstar = 1.0D0/(1.0D0+beta*&
            &sqrt(0.5*spec(n)%mass*spec(n)%temp)*vp(l)*&
            &currdens_par(pi1:pi2)/(spec(n)%charge*geom%Bfield(pi1:pi2,j,k)**2))
    else
       B_Bstar = 1.0
    endif    
    
  end function B_Bstar

  
  !> Return electrostatic drive term
  function edr(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2) :: edr

    integer :: pni
    if (pn0.gt.1) then
       pni=n
    else
       pni=pn1
    endif

    ! Prefactor to electrostatic drive term
    edr = (spec(n)%omn_prof(pi1:pi2) + spec(n)%omt_prof(pi1:pi2) * &
         &((mu(m)*geom%Bfield(pi1:pi2,j,k)+vp(l)**2)/spec(n)%temp_prof(pi1:pi2)- &
         & 1.5)) * fm(pi1:pi2,j,k,l,m,pni)*B_Bstar(j,k,l,n)/C_xy(pi1:pi2)
    
    
    if (xy_local.and.pfsrate.ne.0.) then
       if (pfsrate.eq.-1111.0) pfsrate=ExBrate
       if (magn_geometry=='s_alpha'.or.magn_geometry=='s_alpha_B') then
          !due to the approximations in s_alpha geometry, the expression for general 
          !geometry breaks down, and we use a simplified one
          edr=edr+pfsrate*q0/trpeps*B_Bstar(j,k,l,n)*fm(:,j,k,l,m,pni)/&
               spec(n)%temp*spec(n)%mass*vTvpar(l,n)/C_xy(pi1:pi2)
       else
          edr=edr+pfsrate*B_Bstar(j,k,l,n)*fm(:,j,k,l,m,pni)/&
               spec(n)%temp*spec(n)%mass*vTvpar(l,n)*q_prof(pi1:pi2)*&
               geom%R_hat(pi1:pi2,k)**2/C_y(pi1:pi2)/geom%jacobian(pi1:pi2,j,k)&
               /geom%Bfield(pi1:pi2,j,k)
       endif
    endif
    if (with_comoving_other) then
       edr=edr + &
            !addition to omt-term
            spec(n)%omt_prof(pi1:pi2) * (spec(n)%charge * phi0(pi1:pi2,k) - 0.5 * spec(n)%mass * Omega0_tor**2 * &
            &(geom%R_hat(pi1:pi2,k)**2 - R0_tor**2)) / spec(n)%temp &
            * fm(pi1:pi2,j,k,l,m,pni) * B_Bstar(j,k,l,n) / C_xy(pi1:pi2) + &
            !additional flow shear term
            ExBrate / C_y(pi1:pi2) * B_Bstar(j,k,l,m) * fm(:,j,k,l,m,pni) / spec(n)%temp * spec(n)%mass * Omega0_tor *&
            (geom%R_hat(pi1:pi2,k)**2 - R0_tor**2) / C_xy(pi1:pi2) +&
            !additional term due to radial variation of R0 (if applicable)
            B_Bstar(j,k,l,n) / C_xy(pi1:pi2) * fm(pi1:pi2,j,k,l,m,pni) / &
            spec(n)%temp * spec(n)%mass * Omega0_tor**2 *&
            R0_tor * dR0_tor_dx
    endif

  end function edr

  !> Return curvature term
  function curv(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2) :: curv

    curv = B_Bstar(j,k,l,n)*spec(n)%temp*(mu(m)*&
         &geom%Bfield(pi1:pi2,j,k)+2.0*vp(l)**2)/&
         (spec(n)%charge*geom%Bfield(pi1:pi2,j,k))
    
  end function curv

  !> Return coriolis term for 'i' coordinate
  function coriolis_i(j,k,l,n)
    integer, intent(in) :: j,k,l,n
    real, dimension(pi1:pi2):: angle_gx, angle_gz, gt_dot_gz, coriolis_i
    integer:: i

    !angle between grad(x) and grad(R)
    do i=pi1,pi2
       angle_gx(i)=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
       if (yx_order) then
          !angle between grad(z) and grad(R)
          angle_gz(i)=angle_gx(i)+acos(geom%gjz(i,j,k)/sqrt(geom%gjj(i,j,k)*geom%gzz(i,j,k)))
       else
          angle_gz(i)=angle_gx(i)+acos(geom%giz(i,j,k)/sqrt(geom%gii(i,j,k)*geom%gzz(i,j,k)))
       endif
       !grad(theta)*grad(Z)
       gt_dot_gz(i)=sqrt(geom%gzz(i,j,k))*sin(angle_gz(i))
    enddo

    coriolis_i = 2*spec(n)%mass*vTvpar(l,n)*Omega0_tor/spec(n)%charge/&
        &geom%Bfield(pi1:pi2,j,k)*B_Bstar(j,k,l,n)
    if (yx_order) then
       coriolis_i=coriolis_i * C_y(pi1:pi2)*(q_prof(pi1:pi2)*gt_dot_gz+zval(k)*&
        &dqdx_prof(pi1:pi2)*geom%dxdZ(pi1:pi2,k))
    else
       coriolis_i=coriolis_i * geom%dxdZ(pi1:pi2,k)
    endif

    
  end function coriolis_i

  !> Return coriolis term for 'j' coordinate
  function coriolis_j(j,k,l,n)
    integer, intent(in) :: j,k,l,n
    real, dimension(pi1:pi2):: angle_gx, angle_gz, gt_dot_gz, coriolis_j
    integer:: i

    !alignment of grad(x)
    do i=pi1,pi2
       angle_gx(i)=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
       if (yx_order) then
          angle_gz(i)=angle_gx(i)+acos(geom%gjz(i,j,k)/sqrt(geom%gjj(i,j,k)*geom%gzz(i,j,k)))
       else
          angle_gz(i)=angle_gx(i)+acos(geom%giz(i,j,k)/sqrt(geom%gii(i,j,k)*geom%gzz(i,j,k)))
       endif
       !grad(theta)*grad(Z)
       gt_dot_gz(i)=sqrt(geom%gzz(i,j,k))*sin(angle_gz(i))
    enddo

    coriolis_j = 2*spec(n)%mass*vTvpar(l,n)*Omega0_tor/spec(n)%charge/&
        &geom%Bfield(pi1:pi2,j,k)*B_Bstar(j,k,l,n)
    if (yx_order) then
       coriolis_j=coriolis_j * geom%dxdZ(pi1:pi2,k)
    else
       coriolis_j=coriolis_j * C_y(pi1:pi2)*(q_prof(pi1:pi2)*gt_dot_gz+zval(k)*&
        &dqdx_prof(pi1:pi2)*geom%dxdZ(pi1:pi2,k))
    endif
  end function coriolis_j

  !> Return centrifugal drift term for 'i' coordinate
  function centrifugal_i(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2):: dzdR, dydR, centrifugal_i
    integer:: i 
    real:: angle_gx, angle_gx_gz

    !angle_gx=angle between grad(x) and grad(R)
    !angle_gx_gz=angle between grad(x) and grad(z)
    !angle_gx+angle_gx_gz= angle between grad(z) and grad(R) and thus dz/dR 
    do i=pi1,pi2
       if (yx_order) then
          angle_gx=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
          angle_gx_gz=acos(geom%gjz(i,j,k)/sqrt(geom%gjj(i,j,k)*geom%gzz(i,j,k)))
       else
          angle_gx=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
          angle_gx_gz=acos(geom%giz(i,j,k)/sqrt(geom%gii(i,j,k)*geom%gzz(i,j,k)))
       endif
       dzdR(i)=cos(angle_gx+angle_gx_gz)*sqrt(geom%gzz(i,j,k))
       dydR(i)=C_y(i)*(dqdx_prof(i)*geom%dxdR(i,k)*zval(k)+q_prof(i)*dzdR(i))
    enddo

    centrifugal_i = - spec(n)%mass*Omega0_tor**2*geom%R_hat(pi1:pi2,k)*C_xy(pi1:pi2)/spec(n)%charge/&
        &geom%Bfield(pi1:pi2,j,k)**2* B_Bstar(j,k,l,n)
    if (yx_order) then
       centrifugal_i=centrifugal_i * (dydR(pi1:pi2)*geom%gij(pi1:pi2,j,k)-geom%dxdR(pi1:pi2,k)*geom%gii(pi1:pi2,j,k))
    else
       centrifugal_i=centrifugal_i * (dydR(pi1:pi2)*geom%gii(pi1:pi2,j,k)-geom%dxdR(pi1:pi2,k)*geom%gij(pi1:pi2,j,k))
    endif
    
  end function centrifugal_i


  !> Return centrifugal drift term for 'j' coordinate
  function centrifugal_j(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2):: dzdR, dydR, centrifugal_j
    integer:: i
    real:: angle_gx, angle_gx_gz

    do i=pi1,pi2
       if (yx_order) then
          angle_gx=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
          angle_gx_gz=acos(geom%gjz(i,j,k)/sqrt(geom%gjj(i,j,k)*geom%gzz(i,j,k)))
       else
          angle_gx=atan2(geom%dxdZ(i,k),geom%dxdR(i,k))
          angle_gx_gz=acos(geom%giz(i,j,k)/sqrt(geom%gii(i,j,k)*geom%gzz(i,j,k)))
       endif
       dzdR(i)=cos(angle_gx+angle_gx_gz)*sqrt(geom%gzz(i,j,k))
       dydR(i)=C_y(i)*(dqdx_prof(i)*geom%dxdR(i,k)*zval(k)+q_prof(i)*dzdR(i))
    enddo

    centrifugal_j = - spec(n)%mass*Omega0_tor**2*geom%R_hat(pi1:pi2,k)*C_xy(pi1:pi2)/spec(n)%charge/&
        &geom%Bfield(pi1:pi2,j,k)**2*B_Bstar(j,k,l,n)
    if (yx_order) then
       centrifugal_j=centrifugal_j * (dydR(pi1:pi2)*geom%gjj(pi1:pi2,j,k)-geom%dxdR(pi1:pi2,k)*geom%gij(pi1:pi2,j,k))
    else
       centrifugal_j=centrifugal_j * (dydR(pi1:pi2)*geom%gij(pi1:pi2,j,k)-geom%dxdR(pi1:pi2,k)*geom%gjj(pi1:pi2,j,k))
    endif
    
  end function centrifugal_j
  
  function bxphi0_i(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2):: bxphi0_i

    if (yx_order) then
       bxphi0_i = -(B_Bstar(j,k,l,n)/C_xy(pi1:pi2)*dphi0dx_part1(pi1:pi2,k)+&
            &dphi0_fac(pi1:pi2,k)*centrifugal_i(j,k,l,m,n)*spec(n)%charge/spec(n)%mass)
    else
       bxphi0_i = -dphi0_fac(pi1:pi2,k)*centrifugal_i(j,k,l,m,n)*spec(n)%charge/spec(n)%mass
    endif

  end function bxphi0_i

  function bxphi0_j(j,k,l,m,n)
    integer, intent(in) :: j,k,l,m,n
    real, dimension(pi1:pi2):: bxphi0_j

    if (yx_order) then    
       bxphi0_j = -dphi0_fac(pi1:pi2,k)*centrifugal_j(j,k,l,m,n)*spec(n)%charge/spec(n)%mass
    else
       bxphi0_j = -(B_Bstar(j,k,l,n)/C_xy(pi1:pi2)*dphi0dx_part1(pi1:pi2,k)+&
            &dphi0_fac(pi1:pi2,k)*centrifugal_j(j,k,l,m,n)*spec(n)%charge/spec(n)%mass)
    endif

  end function bxphi0_j


  !>Sets the global variables maxclimitx/maxclimity, which used to estimate the
  !!constraints on dt from perpendicular advection (also used for nonlinear time
  !!step adaption)
  subroutine maxclimits
    integer :: o

    ! Perpendicular time step limit
       
    maxclimit(1)=maxval(abs(pdg1di))
    maxclimit(2)=maxval(abs(pdg1dj))
    do o=1,4
       call my_real_max_to_all(maxclimit(o))
    enddo

  end subroutine maxclimits


  Subroutine finalize_prefactors

    deallocate(vTvpar)
    if (n_fields.gt.2) deallocate(mu_tjqj)
    deallocate(pdg1di,pdg1dj)
    deallocate(pdchibardi,pdchibardj)
    if (allocated(edr_for_energy)) deallocate(edr_for_energy) 
         !,edr_for_energy2)

  end subroutine finalize_prefactors


end module prefactors
