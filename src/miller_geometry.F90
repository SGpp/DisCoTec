#include "redef.h"
!>Implementation of the local equilibrium model of [R.L. Miller et al., PoP 5, 973 (1998)
!>Eventually will be expanded to include a 3-D generalization. 
MODULE miller_mod
  use discretization
  use coordinates
  use par_geom
  use fourier
  use par_in, only: fourier2D, beta
  use par_other, only: print_ini_msg
  use lagrange_interpolation

  implicit none
  public:: get_miller
  public:: rho, kappa, delta, s_kappa, s_delta, drR, drZ, zeta, s_zeta
  public:: major_Z
  public:: npts_pol,npts_tor, set_miller_defaults

  private

  INTEGER:: n_pol,thisunit
  REAL:: r0, a, Lnorm, Psi0                     ! r0/Lref where r0 radius at x=0
  INTEGER :: NPSI,NCHI, NPSI0, NCHI0, npts_tor, npts_pol, nmodespol, nmodestor
  INTEGER :: npol_max, pol
  REAL, DIMENSION(:), ALLOCATABLE :: theta, zzeta, theta2 
  REAL :: B0, R_0, major_Z, rho, kappa, delta, s_kappa, s_delta, dPsidrho, drR, drZ, PI, V_p
  REAL :: zeta, s_zeta
  REAL :: dp_dpsi, s_param, diota_dpsi, q_sigma, q_tau, q_p, mu_0, iota
  REAL, DIMENSION(:,:), ALLOCATABLE :: kappa_n,kappa_g,tau_n
  REAL, DIMENSION(:), ALLOCATABLE :: R_s,Z_s,dxdR_s,dxdZ_s,R_s_theta,Z_s_theta
  REAL, DIMENSION(:), ALLOCATABLE :: R_s_theta_theta,Z_s_theta_theta, sqrtg_int
  REAL, DIMENSION(:,:), ALLOCATABLE :: B2, gPsi2, sqrtg
  COMPLEX, PARAMETER :: imag=(0.0,1.0)
  LOGICAL :: fourier2D_sav

  ! SUBROUTINES DEFINITIONS
CONTAINS

  subroutine set_miller_defaults
    rho = -1.0
    kappa = 1.0
    s_kappa = 0.0
    delta = 0.0
    s_delta = 0.0
    drR = 0.0
    drZ = 0.0
    major_Z = 0.0
    zeta = 0.0
    s_zeta = 0.0
  end subroutine set_miller_defaults


  subroutine get_miller(geom,C_y,C_xy,dpdx_pm_geom,q_prof,dqdx_prof,Cyq0_x0,edge_opt)

    type(geomtype),intent(inout):: geom
    real, dimension(pi1gl:pi2gl), intent(out):: C_y, C_xy, dpdx_pm_geom, q_prof,&
         dqdx_prof
    real, intent(out):: Cyq0_x0
    real, intent(in):: edge_opt

    if (magn_geometry.eq.'miller') then
       call get_miller_a(geom,C_y,C_xy,dpdx_pm_geom,q_prof,dqdx_prof,Cyq0_x0,edge_opt)
    elseif (magn_geometry.eq.'miller_b') then
       call get_miller_b(geom,C_y,C_xy,dpdx_pm_geom)
       q_prof=q0
       dqdx_prof=shat*q0/x0
       Cyq0_x0=1.0
    endif
  end subroutine get_miller

  FUNCTION inv_fft_gen(coefs,eta,neta,alpha_0,pol_indices,nmodes,iota,flag) RESULT(output)

    INTEGER, INTENT(IN) :: nmodes, neta, flag
    REAL, DIMENSION(1:neta), INTENT(IN) :: eta
    REAL, INTENT(IN) :: alpha_0, iota
    REAL, DIMENSION(nz0+1) :: output
    INTEGER, DIMENSION(nmodes), INTENT(IN) :: pol_indices
    COMPLEX, DIMENSION(:,:), INTENT(IN) :: coefs
    INTEGER :: i,pol
    COMPLEX, DIMENSION(nz0+1) :: output_temp

    output_temp(:) = 0     
    output_temp = coefs(1,1)

    DO i=2,nmodes
       pol = pol_indices(i)
       output_temp(1:nz0+1) = output_temp(1:nz0+1) + 2.0*coefs(i,1)*EXP(imag*REAL(pol)*alpha_0) &
            &  * EXP(imag*REAL(pol)*iota*eta(:))
    END DO

    IF( flag .EQ. 1 ) THEN
    output = REAL(output_temp)
    ELSE
    output = AIMAG(output_temp)
    END IF
  END FUNCTION inv_fft_gen

  FUNCTION inv_fft_db(coefs,eta,neta,alpha_0,pol_indices,nmodes,iota) RESULT(output)
    INTEGER, INTENT(IN) :: nmodes, neta
    REAL, DIMENSION(neta), INTENT(IN) :: eta
    REAL, INTENT(IN) :: alpha_0, iota
    REAL, DIMENSION(neta) :: output
    INTEGER, DIMENSION(nmodes), INTENT(IN) :: pol_indices
    COMPLEX, DIMENSION(:,:), INTENT(IN) :: coefs
    INTEGER :: i,npts,pol
    COMPLEX, DIMENSION(neta) :: output_temp

    npts = SIZE(eta,1)
    output_temp(:) = 0

    DO i=1,nmodes
       pol = pol_indices(i)
       output_temp(:) = output_temp(:) + imag*(iota*REAL(pol))* &
            &  2.0*coefs(i,1)*EXP(imag*REAL(pol)*alpha_0) &
            &  * EXP(imag*iota*REAL(pol)*eta(:))
    END DO

    output = REAL(output_temp)

  END FUNCTION inv_fft_db

  FUNCTION inv_fft_arc_q(coefs,eta,neta,alpha_0,pol_indices,nmodes,iota) RESULT(output)
    INTEGER, INTENT(IN) :: nmodes, neta
    REAL, DIMENSION(neta), INTENT(IN) :: eta
    REAL, INTENT(IN) :: alpha_0, iota
    REAL, DIMENSION(neta) :: output
    INTEGER, DIMENSION(nmodes), INTENT(IN) :: pol_indices
    COMPLEX, DIMENSION(:,:), INTENT(IN) :: coefs
    INTEGER :: i, npts, pol
    COMPLEX, DIMENSION(neta) :: output_temp

    npts = SIZE(eta,1)
    output_temp(:) = 0
    output_temp(:) = coefs(1,1)*eta(:)

    DO i=2,nmodes
       pol = pol_indices(i)
       output_temp(:) = output_temp(:) -2.0*imag*coefs(i,1)*  &
            &   EXP(iota*imag*REAL(pol)*eta(:)) / (iota*REAL(pol)+10**(-10))
    END DO

    output = REAL(output_temp)

  END FUNCTION inv_fft_arc_q


  ! Top level subroutine which calculates the Miller equilibrium.
  subroutine get_miller_b(geom, C_y_out,C_xy_out,my_dpdx_pm)
    Implicit none

    type(geomtype),intent(inout):: geom
    real, dimension(pi1gl:pi2gl),intent(out):: C_y_out,C_xy_out,my_dpdx_pm

    REAL    :: int_tor, int_pol, int_eta, r_avg, r_ab
    INTEGER :: i, k, nz0n

    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: gB2_fft, q_sigma_fft,  &
         & q_lambda_fft, lambda_coef, &
         & q_tau_fft, q_p_fft, temp_fft
    REAL, DIMENSION(:), ALLOCATABLE :: modB_fl, modgradpsi_fl, kappa_g_fl, kappa_n_fl, s_fl, &
         & lambda_fl, tau_n_fl, modB_deriv_fl, sqrtg_fl, arc_length, q_LL_fl, &
         & eta, dB_deriv_fl, LL, R_s_fl, Z_s_fl, dxdR_s_fl, dxdZ_s_fl
    INTEGER, DIMENSION(:), ALLOCATABLE :: pol_indices
    REAL, DIMENSION(:,:), ALLOCATABLE :: gB, gB2, q_sigma_vec, q_tau_vec, q_p_vec, q_LL, &
         &lambda, sigma, s_int, modB, modgPsi, q_lambda
    real, dimension(lk1:lk2):: gxx,gxz,gyy,gyz,dBdx,dBdy
    character(len=6) :: Lref_str
    logical :: consistency_check = .false.
    
    !need to enforce 2d fourier transforms
    fourier2D_sav = fourier2D
    fourier2D=.true.

    nz0n = nz0+1
    npts_pol = 500
    IF(y_local) THEN
       npts_tor = 1
    ELSE
       npts_tor = 500
    ENDIF

    npol_max = npts_pol / 2
    nmodespol = npts_pol/2 + 1
    nmodestor = npts_tor

    R_0 = major_R
    iota = 1.0 / q0
    PI = ACOS(-1.0)
    mu_0 = 4.0 * PI * (10**(-7.0))  
    Lref_str = 'L_ref'
    if (minor_r.eq.1.0) Lref_str = 'a'
    if (major_R.eq.1.0) Lref_str = 'R_0'
    if (rho.lt.0.0) rho = trpeps*major_R
    if (rho.le.0.0) stop 'flux surface radius not defined'
    trpeps = rho/major_R

    ALLOCATE(R_s(npts_pol),Z_s(npts_pol))
    ALLOCATE(B2(npts_pol,npts_tor),gPsi2(npts_pol,npts_tor))
    ALLOCATE(dxdR_s(npts_pol),dxdZ_s(npts_pol))
    ALLOCATE(R_s_theta(npts_pol),Z_s_theta(npts_pol))
    ALLOCATE(R_s_theta_theta(npts_pol),Z_s_theta_theta(npts_pol)) 
    ALLOCATE(sqrtg_int(npts_pol) )

    ALLOCATE(gB2_fft(nmodespol,npts_tor), q_sigma_fft(nmodespol,npts_tor), &
         & q_lambda_fft(nmodespol,npts_tor), lambda_coef(nmodespol,npts_tor), &
         & q_p_fft(nmodespol,npts_tor), &
         & q_tau_fft(nmodespol,npts_tor), temp_fft(nmodespol,npts_tor))

    ALLOCATE(pol_indices(nmodespol), gB(npts_pol,npts_tor), &
         & gB2(npts_pol,npts_tor), q_sigma_vec(npts_pol,npts_tor), q_tau_vec(npts_pol,npts_tor), &
         & q_p_vec(npts_pol,npts_tor), q_LL(npts_pol,npts_tor), lambda(npts_pol,npts_tor), &
         & sigma(npts_pol,npts_tor), s_int(npts_pol,npts_tor), q_lambda(npts_pol,npts_tor), &
         & modB(npts_pol,npts_tor), modgPsi(npts_pol,npts_tor), sqrtg(npts_pol,npts_tor))

    ALLOCATE(modB_fl(nz0n), modgradpsi_fl(nz0n), kappa_g_fl(nz0n), kappa_n_fl(nz0n), &
         & s_fl(nz0n), lambda_fl(nz0n), tau_n_fl(nz0n), &
         & modB_deriv_fl(nz0n), sqrtg_fl(nz0n), arc_length(nz0n), q_LL_fl(nz0n), &
         & eta(nz0n), dB_deriv_fl(nz0n), LL(nz0n), &
         & R_s_fl(nz0n), Z_s_fl(nz0n),dxdR_s_fl(nz0n), dxdZ_s_fl(nz0n))

    ALLOCATE(theta(npts_pol),zzeta(npts_tor))

!    CALL finalize_fourier
    CALL initialize_fourier(npts_tor,npts_pol)

    DO i=1,npol_max+1
       pol_indices(i) = (i-1)
    ENDDO

    int_tor = (2*PI)/(REAL(npts_tor))
    int_pol = (2*PI)/(REAL(npts_pol))
    int_eta = (2*PI)/(iota*REAL(nz0n-1)) 

    DO i=1,npts_pol
       theta(i) = 0.0 + (i-1)*int_pol
    ENDDO

    DO i=1,nz0n
       eta(i) = -PI/iota + (i-1)*int_eta
    ENDDO

    DO i=1,npts_tor
       zzeta(i) = 0 + (i-1)*int_tor
    ENDDO

    B0 = 1.0 !1 / iota

    CALL mapping_miller

    gB  = sqrtg*(B2**(0.5))
    gB2 = sqrtg*(B2)

    modB    = B2**(0.5)
    modgPsi = gPsi2**(0.5)

    r_avg = sum(sqrt((R_s-major_R)**2.0+Z_s**2.0)*sqrtg_int)/sum(sqrtg_int)
    r_ab = abs(R_s(1)-R_s(npts_pol/2))*0.5

    Lnorm = r_avg !minor_r
    x0 = Lnorm
    Psi0 = B0 * Lnorm**2

    !definition of radial coordinate is here:
    !x = Psi/Psi0 * (q0 * Lnorm)

    if ((mype.eq.0).and.(print_ini_msg)) then
       WRITE(*,*)
       WRITE(*,'(A)') "*** Miller parameters ***"
       WRITE(*,"(1A,ES20.10)") "* rho:    ",rho
       WRITE(*,"(1A,ES20.10)") "* r_avg:  ",r_avg
       WRITE(*,"(1A,ES20.10)") "* r_ab:   ",r_ab
       WRITE(*,"(1A,ES20.10)") "* Psi0:   ",Psi0
       WRITE(*,"(1A,ES20.10)") "* V':     ",V_p
       WRITE(*,*)
       WRITE(*,'(A)') "Gradient conversion factor (x=Psi/Psi0*q0*r_avg)"
!       WRITE(*,"(A,1(ES20.10))") "* dr/dPsi:            ", 1./dPsidrho
       WRITE(*,"(A,1(ES20.10))") "* dr/dx:            ", 1./dPsidrho*Psi0/Lnorm*iota
!       WRITE(*,"(A,1(ES20.10))") "* Psi0*grad Psi_ab:   ", &
!            &((( (modgPsi(1,1)**(-1.0))+(modgPsi(npts_pol/2,1)**(-1.0)))/2.0))*Psi0/Lnorm
!       WRITE(*,"(A,1(ES20.10))") "* Psi0*grad Psi(0):   ", 1.0/modgPsi(1,1)*Psi0/Lnorm
!       WRITE(*,"(A,1(ES20.10))") "* Psi0*grad Psi(pi):  ", 1.0/modgPsi(npts_pol/2,1)*Psi0/Lnorm
       WRITE(*,'(A)') "Magnetic and ExB shear conversion factor"
       WRITE(*,"(A,1(ES20.10))") "* Psi0/rho*dr/dPsi: ", 1./dPsidrho*Psi0/rho*iota       
       WRITE(*,'(A)') "B_unit/Bref conversion factor"
       WRITE(*,"(A,1(ES20.10))") "* q/rho*dPsi/drho:  ", q0/rho*dPsidrho
       WRITE(*,*)
    endif

    !!NOTE: As for the gradients, the user has to make sure that shat
    !!is computed in the psi coordinate 
    !!(we cannot simply multiply the gradient factor on top of shat
    !! as this would be applied several times during the autoparallelization;
    !! adding in the diota_dpsi definition yields different shat in-/output 
    !! values)
    diota_dpsi = -iota*shat/Psi0 !*(Lnorm/rho) !internal r normalization is Lnorm

    !using for a 1st test: amhd = -q^2 R0 8 pi / Bref^2 dp/drho (circular geometries)
    !user is required to apply conversion factor to psi coordinate by hand
    !in input parameters
    dp_dpsi = -amhd * (iota)**2.0 / ( 2.0 * mu_0) * B0**2 / R_0 * & 
         & (Lnorm / Psi0) !for some reason, the agreement with other codes is worse
                          !if we add a factor of q0 here (see conversion factor)

    !neg. dpdx normalized to magnetic pressure for pressure term
    my_dpdx_pm = iota**2.0 * amhd / R_0

    ! remember, needs to be imag
    q_lambda = 2.0*mu_0*kappa_g*modgPsi * sqrtg / (modB)
    q_sigma_vec = sqrtg* ( (modB/modgPsi)**2.0)

    q_tau_vec = q_sigma_vec * tau_n
    q_p_vec   = q_sigma_vec
    q_LL      = q_sigma_vec
   
    gB2_fft = 0.0
    q_sigma_fft = 0.0
    q_lambda_fft = 0.0

    CALL fft_xy_to_ff(gB2,gB2_fft)
    CALL fft_xy_to_ff(q_sigma_vec,q_sigma_fft)
    CALL fft_xy_to_ff(q_lambda,q_lambda_fft)
    CALL fft_xy_to_ff(q_tau_vec,q_tau_fft)

    q_sigma = REAL(q_sigma_fft(1,1))
    q_tau = REAL(q_tau_fft(1,1))
    lambda_coef = 0.0

    DO i=2,npol_max
       pol = pol_indices(i)
       lambda_coef(i,1) = (q_lambda_fft(i,1)) / (-1.0*(iota*REAL(pol)))
    END DO

    lambda_coef(:,:) = lambda_coef(:,:) * imag    
    lambda_coef(1,1) = 0.0
    lambda(:,1) = 0.0

    DO i=2,npts_pol/2 + 1
       pol = (i-1)
       lambda(:,1) = lambda(:,1) + 2.0*REAL(lambda_coef(i,1))*COS(REAL(pol)*theta(:))
    END DO

    DO i=2,npts_tor
       lambda(:,i) = lambda(:,1)
    END DO

    ! We want to set lambda_coef(1,1) such that surf. avg of lambda*gB2 = 0
    ! Here, we set lambda_coef(1,1) = -<lambda*gB2> / <gB2> where <> is surf. avg
    lambda_coef(1,1) = -SUM(SUM(lambda*gB2,1)) / (REAL(npts_pol)*REAL(npts_tor)*gB2_fft(1,1))
    lambda = lambda + lambda_coef(1,1)

    q_p_vec = q_p_vec * lambda

    q_p_fft = 0.0

    CALL fft_xy_to_ff(q_p_vec,q_p_fft)

    q_p = REAL(q_p_fft(1,1))

    sigma = (diota_dpsi - q_p*dp_dpsi + 2*q_tau) / (q_sigma) 

    s_int = sigma + dp_dpsi*lambda - 2*tau_n
    q_LL = q_LL*s_int

    temp_fft = 0.0
    CALL fft_xy_to_ff(q_LL,temp_fft)
    q_LL_fl = inv_fft_arc_q(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota)

    temp_fft = 0.0
    CALL fft_xy_to_ff(modB,temp_fft)
    modB_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)
    modB_deriv_fl = inv_fft_db(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota)
    
    temp_fft = 0.0
    CALL fft_xy_to_ff(modgPsi,temp_fft)
    modgradpsi_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(kappa_g,temp_fft)
    kappa_g_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(kappa_n,temp_fft)
    kappa_n_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(s_int,temp_fft)
    s_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(lambda,temp_fft)
    lambda_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(tau_n,temp_fft)
    tau_n_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(sqrtg,temp_fft)
    sqrtg_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(R_s,temp_fft)
    R_s_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(Z_s,temp_fft)
    Z_s_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)   
    
    temp_fft = 0.0
    CALL fft_xy_to_ff(dxdR_s,temp_fft)
    dxdR_s_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)

    temp_fft = 0.0
    CALL fft_xy_to_ff(dxdZ_s,temp_fft)
    dxdZ_s_fl = inv_fft_gen(temp_fft,eta,nz0n,0.0,pol_indices,npol_max,iota,1)   

    arc_length = inv_fft_arc_q(gB2_fft,eta,nz0n,0.0,pol_indices,npol_max,iota)
    arc_length = arc_length / Lnorm

    dB_deriv_fl = R_0 * modB_deriv_fl / (sqrtg_fl * modB_fl)   
    LL = -(modgradpsi_fl**(2.0))*q_LL_fl / (modB_fl) ! / ( Lnorm**2.0 / V_p)

    !x = C_x(flux coordinate) = Psi/Psi0 * (q0 Lnorm)
    !y = C_y * field line label coord.
    !C_xy = dPsi/dPsi / (C_x' * C_y)
    !From metric: C_xy = Bfield / sqrt(gxx*gyy-gxy**2) = B0
    C_xy_out = B0

    !--> C_y = C_x'/C_xy = dPsi/dx/B0 = 1/(dx/dPsi)/B0 
    !        = Psi0/(q0*Lnorm*B0) = Lnorm/q0
    C_y_out = Lnorm/q0

    !Attention: index lk1 starts at 0 instead of 1
    do k = lk1,lk2
       gxx(k) = (modgradpsi_fl(k+1)*Lnorm/Psi0)**(2.0)
       geom%gij(pi1gl:pi2gl,pj1,k) = modB_fl(k+1)*LL(k+1)/B0
       gyy(k) = (modB_fl(k+1)/modgradpsi_fl(k+1))**2.0 * (1 + &
            &LL(k+1)**2.0) * (Psi0/(Lnorm*B0))**2.0
       gxz(k) = 0.0
       gyz(k) = 0.0
       geom%Bfield(pi1gl:pi2gl,pj1,k) = modB_fl(k+1)/B0
       geom%jacobian(pi1gl:pi2gl,pj1,k) = sqrtg_fl(k+1)/(V_p*iota) !/(Lnorm/rho)
       dBdx(k) = (modB_fl(k+1) / modgradpsi_fl(k+1)) * (Psi0/(Lnorm*B0)) * &
            &( kappa_n_fl(k+1) - kappa_g_fl(k+1)*LL(k+1))*R_0 - mu_0*dp_dpsi*(Lnorm*R_0)/(modB_fl(k+1))
       dBdy(k) = (R_0*Lnorm/Psi0)*modgradpsi_fl(k+1) * kappa_g_fl(k+1)
       geom%dBdz(pi1,pj1,k) = modB_deriv_fl(k+1) / (B0*iota)
    enddo

    if (consistency_check.and.print_ini_msg.and.(mype==0)) then
       write (*,'(A)') 'Miller consistency checks: '
       !shat consistency check: gxy/gxx(right)-gxy/gxx(left)=2pi*Cy*q0/x0*shat
       write (*,'(A,F12.6)') 'shat from metric: ',(Psi0/Lnorm)**2/B0* &
            &(modB_fl(nz0n)*LL(nz0n)/modgradpsi_fl(nz0n)**(2.0)-&
            & modB_fl(1)*LL(1)/modgradpsi_fl(1)**(2.0))/&
            &(2.0*pi*C_y_out*q0/x0)
       write (*,'(A,F12.6)') 'shat input:       ',shat
       write (*,*)       
       write (*,'(A,F12.6)') 'C_xy from metric: ',geom%Bfield(pi1gl,pj1,lk1)/&
            sqrt(gxx(lk1)*gyy(lk1)-geom%gij(pi1gl,pj1,lk1)**2.0)
       write (*,'(A,F12.6)') 'C_xy output     : ',C_xy_out
       write (*,*)
    endif


    !thisunit=333
    !OPEN(thisunit,file='./miller_out.dat')
    !DO k= lk1,lk2
    !      WRITE(thisunit,"(9(ES20.10))") &
    !           & gxx(k), gxy(k), gyy(k), modB_fl(k+1), sqrtg_fl(k+1), &
    !           & dBdx(k), dBdy(k), geom%dBdz(pi1,pj1,k), q_LL_fl(k+1)
    !END DO
    !CLOSE(thisunit)

    do k = lk1,lk2
       if (xy_local) then
          if (yx_order) then
             geom%gii(pi1gl:pi2gl,pj1,k)=gyy(k)
             geom%giz(pi1gl:pi2gl,pj1,k)=gyz(k)
             geom%gjj(pi1gl:pi2gl,pj1,k)=gxx(k)
             geom%gjz(pi1gl:pi2gl,pj1,k)=gxz(k)
             geom%dBdj(pi1gl:pi2gl,pj1,k)=dBdx(k)
             geom%dBdi(pi1gl:pi2gl,pj1,k)=dBdy(k)
          else
             geom%gii(pi1gl:pi2gl,pj1,k)=gxx(k)
             geom%giz(pi1gl:pi2gl,pj1,k)=gxz(k)
             geom%gjj(pi1gl:pi2gl,pj1,k)=gyy(k)
             geom%gjz(pi1gl:pi2gl,pj1,k)=gyz(k)
             geom%dBdi(pi1gl:pi2gl,pj1,k)=dBdx(k)
             geom%dBdj(pi1gl:pi2gl,pj1,k)=dBdy(k)
          end if
          geom%gzz(pi1gl:pi2gl,pj1,lk1:lk2) = 1.0 !???
       else
          if (yx_order) then
             do i=pi1gl,pi2gl
                geom%gii(i,pj1,k)=gyy(k)
                geom%giz(i,pj1,k)=gyz(k)
                geom%gjj(i,pj1,k)=gxx(k)
                geom%gjz(i,pj1,k)=gxz(k)
                geom%gzz(i,pj1,lk1:lk2) = 1.0
                geom%dBdj(i,pj1,k)=dBdx(k)
                geom%dBdi(i,pj1,k)=dBdy(k)
             enddo
          else
             do i=pi1gl,pi2gl
                geom%gii(i,pj1,k)=gxx(k)
                geom%giz(i,pj1,k)=gxz(k)
                geom%gjj(i,pj1,k)=gyy(k)
                geom%gjz(i,pj1,k)=gyz(k)
                geom%gzz(i,pj1,lk1:lk2) = 1.0
                geom%dBdi(i,pj1,k)=dBdx(k)
                geom%dBdj(i,pj1,k)=dBdy(k)
             enddo
          end if
       endif
    enddo
    geom%gzz(pi1gl:pi2gl,pj1,lk1:lk2) = 1.0 !???

    if (Bref.le.0) Bref = B0
    if (Lref.le.0) Lref = 1.0

    do k=lk1,lk2
       geom%R(pi1gl:pi2gl,k)     = R_s_fl(k+1)*Lref
       geom%PHI(pi1gl:pi2gl,k)   = 0.0 !FILL?
       geom%Z(pi1gl:pi2gl,k)     = Z_s_fl(k+1)*Lref
       geom%R_hat(pi1gl:pi2gl,k) = R_s_fl(k+1)
       geom%Z_hat(pi1gl:pi2gl,k) = Z_s_fl(k+1)
       geom%dxdR(pi1gl:pi2gl,k)  = dxdR_s_fl(k+1)*Lnorm/Psi0
       geom%dxdZ(pi1gl:pi2gl,k)  = dxdZ_s_fl(k+1)*Lnorm/Psi0
    enddo

    DEALLOCATE(pol_indices,gB,gB2,q_sigma_vec,q_tau_vec,q_p_vec,q_LL,lambda,sigma)
    DEALLOCATE(s_int, q_lambda, modB, modgPsi)
    DEALLOCATE(modB_deriv_fl, sqrtg_fl, arc_length, q_LL_fl, eta, dB_deriv_fl, LL)
    DEALLOCATE(theta,zzeta,modB_fl,kappa_n,kappa_g,tau_n,R_s_fl, Z_s_fl,dxdR_s_fl,dxdZ_s_fl)
    DEALLOCATE(R_s,Z_s,B2,gPsi2,sqrtg_int,sqrtg,dxdR_s,dxdZ_s,R_s_theta,Z_s_theta)

    DEALLOCATE(gB2_fft, q_sigma_fft, q_lambda_fft, lambda_coef, &
         & q_tau_fft, q_p_fft)

    CALL finalize_fourier  
    fourier2D = fourier2D_sav
!    CALL initialize_fourier(nx0,nky0)

  end subroutine get_miller_b



  ! Constructs the Miller equilibrium and converts to straight field line coordinates
  SUBROUTINE mapping_miller

    REAL, DIMENSION(:,:), ALLOCATABLE :: R_g, Z_g, R_theta_g, Z_theta_g, &
         &  dX_dtheta_mag_g, Bp_top_g, Bp_bottom_g, &
         &  Bp_g_unnorm, Bp_g, ode_rhs
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ode_rhs_fft, sqrtg_fft 
    REAL, DIMENSION(:), ALLOCATABLE :: q_integrand
    REAL :: inter, xint, int_pol, inter_top
    INTEGER :: i
    REAL, DIMENSION(1:npts_pol) :: Theta_g, & !dtt_dtheta_g, dt_dTheta_g,
         & dt_dTheta_s, d2t_dTheta2_s
    REAL, DIMENSION(1:npts_pol) :: newterm, theta_s, & !costerm,
         & gradPsi1,gradPsi2,gradPsi3,Br,Bphi,Bz,B2_int,gPsi2_int

    ALLOCATE(R_g(npts_pol,npts_tor), Z_g(npts_pol,npts_tor), &
         & R_theta_g(npts_pol,npts_tor), Z_theta_g(npts_pol,npts_tor), &
         & dX_dtheta_mag_g(npts_pol,npts_tor), Bp_top_g(npts_pol,npts_tor), &
         & Bp_bottom_g(npts_pol,npts_tor), Bp_g_unnorm(npts_pol,npts_tor),  &
         & Bp_g(npts_pol,npts_tor), ode_rhs(npts_pol,npts_tor), &
         & ode_rhs_fft(npts_pol/2+1,npts_tor), q_integrand(npts_pol), &
         & sqrtg_fft(npts_pol/2+1,npts_tor) )

    int_pol = (2*PI)/(npts_pol-1)
    xint = ASIN(delta)

    !flux surface shape 
    R_g(:,1) = R_0 + rho*COS( theta + xint*SIN(theta) )
    Z_g(:,1) = major_Z + kappa*rho*SIN( theta + zeta*SIN(2.0*theta) ) 

    !theta derivative
    R_theta_g(:,1) = -rho*SIN(theta + xint*SIN(theta) ) * ( 1 + xint*COS(theta) )
    Z_theta_g(:,1) = kappa*rho*(1.0+2.0*zeta*COS(2.0*theta))*COS(theta + zeta*SIN(2.0*theta))


    dX_dtheta_mag_g = (R_theta_g**2 + Z_theta_g**2)**(0.5)

    !B_poloidal nominator/denominator
    Bp_top_g(:,1) = ( (SIN(theta + xint*sin(theta)) * ( 1.0 + xint*COS(theta)))**2 + &
         & (kappa*(1.0+2.0*zeta*COS(2.0*theta))*COS(theta+zeta*SIN(2.0*theta)))**2 )**(0.5)

    Bp_bottom_g(:,1) = kappa*R_g(:,1)* &
         & ( (drR + COS(theta+xint*SIN(theta)) - s_delta*SIN(theta)*SIN(theta+xint*SIN(theta)))*&
         &   (1.0+2.0*zeta*COS(2.0*theta))*COS(theta+zeta*SIN(2.0*theta))&
         &  +(s_zeta*COS(theta+zeta*SIN(2.0*theta))*SIN(2.0*theta)+(1.0+s_kappa)*&
         &    SIN(theta+zeta*SIN(2.0*theta))+drZ/kappa)*(1.0+xint*COS(theta))*SIN(theta+xint*SIN(theta)) )

    !B_poloidal/(dPsi/dr)  =  |grad r|/R
    Bp_g_unnorm = Bp_top_g / Bp_bottom_g

    q_integrand = dX_dtheta_mag_g(:,1) / ( Bp_g_unnorm(:,1) * (R_g(:,1)**2))

    inter_top = sum(q_integrand)

    inter = inter_top / REAL(npts_pol)

    !Bref set to magnetic axis value if major_R = 1.0?
    dPsidrho = (R_0 * B0 / q0 ) * inter
    Bp_g = Bp_g_unnorm * dPsidrho

    !rhs of ODE in Eq.(7) in Tom's Nucl. Fusion 53, 013004
    !assuming f = 1 on given flux surface
    ode_rhs = ((R_0*B0*dX_dtheta_mag_g)) / (q0 * (R_g**2) * Bp_g)
    CALL fft_xy_to_ff(ode_rhs,ode_rhs_fft)

    !now solve the ODE and construct the straight field line theta
    !using the (odd) Fourier representation
!    costerm = 0.0
    newterm = 0.0
!    costerm = costerm + REAL(ode_rhs_fft(1,1))
    DO i=2,npol_max
       pol = (i-1)
!       costerm(:) = costerm(:) + REAL(2.0*ode_rhs_fft(i,1))*COS(REAL(pol)*theta(:))
       newterm(:) = newterm(:) + (REAL(2.0*ode_rhs_fft(i,1))/REAL(pol))*SIN(REAL(pol)*theta(:))
    ENDDO
    
    !straight field line angle:
    Theta_g = theta + newterm

    !derivatives
!    dt_dTheta_g = 1.0 / (1.0 + costerm)
!    dtt_dtheta_g = 1.0 + costerm

    !create geometrical angle array which is equidistant in straight field line angle
    CALL lag3interp(theta(:),Theta_g(:),npts_pol,theta_s,theta,npts_pol)

    !derivatives (replace with lag3deriv?)
    call lag3deriv(theta_s,theta,npts_pol,dt_dTheta_s,theta,npts_pol)
    call lag3deriv(dt_dTheta_s,theta,npts_pol,d2t_dTheta2_s,theta,npts_pol)

    !flux surface coordinates as functions of theta_s
    R_s(:) = R_0 + rho * COS( theta_s(:) + xint*SIN(theta_s(:)))
    Z_s(:) = kappa * rho * SIN ( theta_s(:) + zeta*SIN(2.0*theta_s(:)) )

!    R_s_rho(:) = drR + COS( theta_s(:) + xint*SIN(theta_s(:))) - &
!         & s_delta * SIN(theta_s(:)) * SIN(theta_s(:) + xint*SIN(theta_s(:))) 
!    Z_s_rho(:) = kappa * ( (1.0 + s_kappa) * SIN ( theta_s(:) + zeta*SIN(2.0*theta_s(:))) + &
!         & s_zeta * COS(theta_s(:) + zeta*SIN(2.0*theta_s(:)))*SIN(2.0*theta_s(:)) )

    R_s_theta(:) = -rho*SIN(theta_s(:) + xint*SIN(theta_s(:)))*(1.0 + xint*COS(theta_s(:))) * (dt_dTheta_s(:))
    Z_s_theta(:) = kappa*rho*(1.0+2.0*zeta*COS(2.0*theta_s(:)))*COS(theta_s(:) + zeta*SIN(2.0*theta_s(:))) * (dt_dTheta_s(:))

    R_s_theta_theta(:) = -rho * cos( theta_s(:) + xint*sin(theta_s(:))) * ( (1+xint*cos(theta_s(:)))**(2.0)) * &
         &                                                                 (dt_dTheta_s(:)**(2.0)) &
         &                   +rho * sin( theta_s(:) + xint*sin(theta_s(:)))*xint*sin(theta_s(:))*( dt_dTheta_s(:)**(2.0)) &
         &                   -rho * sin( theta_s(:) + xint*sin(theta_s(:)))*(1+xint*cos(theta_s(:)))*(d2t_dTheta2_s)

    Z_s_theta_theta(:) = -kappa * rho * (4.0*zeta*sin(2.0*theta_s(:))*cos(theta_s(:)+zeta*sin(2.0*theta_s(:))) + &
         & (1.0+2.0*zeta*cos(2.0*theta_s(:)))**2.0*sin(theta_s(:)+zeta*sin(2.0*theta_s(:)))) * (dt_dTheta_s(:)**(2.0)) &
         & + kappa*rho*(1.0+2.0*zeta*COS(2.0*theta_s(:)))*COS(theta_s(:) + zeta*SIN(2.0*theta_s(:))) * d2t_dTheta2_s

    sqrtg_int(:) = (R_s(:)**(2.0))/(B0*R_0)

    CALL curvatures

    Br = (iota*R_s_theta)/sqrtg_int
    Bphi = (-R_s) / sqrtg_int
    Bz = (iota*Z_s_theta)/sqrtg_int

    B2_int = Br**2.0 + Bphi**2.0 + Bz**2.0

    gradPsi1 = (Z_s_theta*R_s) / sqrtg_int
    gradPsi2 = 0
    gradPsi3 = (-R_s_theta*R_s) / sqrtg_int

    dxdR_s = gradPsi1/iota
    dxdZ_s = gradPsi3/iota

    gPsi2_int = gradPsi1**2.0 + gradPsi2**2.0 + gradPsi3**2.0

    DO i=1,npts_tor
       B2(:,i) = B2_int(:)
       gPsi2(:,i) = gPsi2_int(:)
       sqrtg(:,i) = sqrtg_int(:)
    ENDDO

    CALL fft_xy_to_ff(sqrtg,sqrtg_fft)
    V_p = R_0/B0 !sqrtg_fft(1,1)
!    V_p = sqrtg_fft(1,1)

    DEALLOCATE(R_g,Z_g,R_theta_g,Z_theta_g,Bp_top_g,Bp_bottom_g,dX_dtheta_mag_g, Bp_g_unnorm, &
         &   q_integrand, Bp_g, ode_rhs, ode_rhs_fft, sqrtg_fft)  


  END SUBROUTINE mapping_miller

  ! Constructs the geometric curvature quantities

  SUBROUTINE curvatures

    REAL, DIMENSION(:,:), ALLOCATABLE :: R_int, &
         &              R_zeta,Z_zeta,Phi_zeta,Lcurv,Q, &
         &              Phi_zeta_zeta,kappa_g_top_1, &
         &              kappa_g_top_2,kappa_g_top_3,tau_n_top_1, &
         &              tau_n_top_2,tau_n_top_3,tau_n_top_4, &
         &              R_theta_zeta,Z_theta_zeta, &
         &              Phi_theta, Phi_theta_theta, Phi_theta_zeta, R_zeta_zeta, Z_zeta_zeta, &
         &              R_theta_int,Z_theta_int,R_theta_theta_int,Z_theta_theta_int,kg_t,kn_t,tn_t
   
    INTEGER :: i

    ALLOCATE(kappa_n(npts_pol,npts_tor),kappa_g(npts_pol,npts_tor), &
         &         tau_n(npts_pol,npts_tor), &
         &         Phi_theta(npts_pol,npts_tor), &
         &         R_zeta(npts_pol,npts_tor), & 
         &         Phi_zeta(npts_pol,npts_tor), &
         &         Z_zeta(npts_pol,npts_tor), & 
         &         Lcurv(npts_pol,npts_tor),Q(npts_pol,npts_tor), & 
         &         Phi_theta_theta(npts_pol,npts_tor), & 
         &         Phi_zeta_zeta(npts_pol,npts_tor), & 
         &         Phi_theta_zeta(npts_pol,npts_tor), & 
         &         R_theta_zeta(npts_pol,npts_tor), &
         &         R_zeta_zeta(npts_pol,npts_tor), &
         &         Z_theta_zeta(npts_pol,npts_tor), &
         &         Z_zeta_zeta(npts_pol,npts_tor), &
         &         kappa_g_top_1(npts_pol,npts_tor), & 
         &         kappa_g_top_2(npts_pol,npts_tor), & 
         &         kappa_g_top_3(npts_pol,npts_tor), & 
         &         tau_n_top_1(npts_pol,npts_tor), &
         &         tau_n_top_2(npts_pol,npts_tor), &
         &         tau_n_top_3(npts_pol,npts_tor), &
         &         tau_n_top_4(npts_pol,npts_tor), &
         &         R_int(npts_pol,npts_tor), R_theta_int(npts_pol,npts_tor), &
         &         Z_theta_int(npts_pol,npts_tor), R_theta_theta_int(npts_pol,npts_tor), &
         &         Z_theta_theta_int(npts_pol,npts_tor),kg_t(npts_pol,npts_tor), &
         &         kn_t(npts_pol,npts_tor), tn_t(npts_pol,npts_tor) ) 

    DO i=1,npts_tor
       R_int(:,i) = R_s(:)
       Z_theta_int(:,i) = Z_s_theta(:)
       R_theta_int(:,i) = R_s_theta(:)
       Z_theta_theta_int(:,i) = Z_s_theta_theta(:)
       R_theta_theta_int(:,i) = R_s_theta_theta(:)
    ENDDO
    Phi_theta(:,:)= 0.0
    R_zeta(:,:) = 0.0
    Phi_zeta(:,:)=-1.0
    Z_zeta(:,:)= 0.0 
    Phi_theta_theta(:,:) = 0
    Phi_theta_zeta(:,:) = 0
    Phi_zeta_zeta(:,:) = 0.0
    R_zeta_zeta(:,:) = 0
    Z_zeta_zeta(:,:) = 0
    R_theta_zeta(:,:) = 0
    Z_theta_zeta(:,:) = 0


!    iota = 1.0 / 3.0

    Lcurv = sqrt( (R_zeta + iota*R_theta_int)**2 &
         &         + (Z_zeta + iota*Z_theta_int)**2 &
         &         + (R_int**2) * (Phi_zeta + iota*Phi_theta)**2) 

    Q = ( (R_int**2)*(Phi_theta*Z_zeta - &
         &           Phi_zeta*Z_theta_int)**2 &
         &         + (R_int**2)*(R_theta_int*Phi_zeta - &
         &           R_zeta*Phi_theta)**2 &
         &         + (Z_theta_int*R_zeta - &
         &           Z_zeta*R_theta_int)**2)**(1.0/4.0)


    kn_t=((Phi_theta*Z_zeta - &
         &                  Phi_zeta*Z_theta_int)* &
         &               (R_int*R_zeta_zeta + &
         &                 iota*iota*R_int*R_theta_theta_int &
         &                 +2.0*iota*R_int*R_theta_zeta &
         &                 -R_int*R_int*( &
         &                  (Phi_zeta+iota*Phi_theta)**2)) &
         &             +(Z_theta_int*R_zeta-Z_zeta*R_theta_int)* &
         &               (R_int*Phi_zeta_zeta & 
         &                 +iota*iota*Phi_theta_theta*R_int &
         &                 +2.0*iota*R_int*Phi_theta_zeta & 
         &                 +2.0*(R_zeta + iota*R_theta_int)* &
         &                   (Phi_zeta + iota*Phi_theta)) &
         &             +(R_theta_int*Phi_zeta - &
         &                  R_zeta*Phi_theta)* &
         &               (R_int*Z_zeta_zeta &
         &                 +iota*iota*R_int*Z_theta_theta_int &
         &                 +2.0*iota*R_int*Z_theta_zeta) ) &
         &             / ( (Lcurv**2.0) * (Q**2.0) ) 

    kappa_n = kn_t

    kappa_g_top_1= &
         &      (R_theta_int*Phi_zeta-R_zeta*Phi_theta) &
         &       *( R_int*R_int *  & 
         &                (Phi_zeta+iota*Phi_theta)* &
         &                (R_zeta_zeta + iota*iota*R_theta_theta_int &
         &                + 2.0*iota*R_theta_zeta) &
         &         -R_int*R_int*R_int*  & 
         &                (Phi_zeta+iota*Phi_theta)**(3.0) &
         &         -R_int*R_int *  & 
         &                (R_zeta + iota*R_theta_int)* &
         &                (Phi_zeta_zeta +  &
         &                 iota*iota*Phi_theta_theta + & 
         &                 2.0*iota*Phi_theta_zeta) &
         &         -2.0*R_int* & 
         &                ( (R_zeta+iota*R_theta_int)**(2.0))* &
         &                (Phi_zeta+iota*Phi_theta) )



    kappa_g_top_2= &
         &      (Z_theta_int*R_zeta - Z_zeta*R_theta_int) &
         &       *( -(Z_zeta+iota*Z_theta_int)* & 
         &           (R_zeta_zeta+iota*iota*R_theta_theta_int &
         &            +2.0*iota*R_theta_zeta) & 
         &          +R_int*(Z_zeta+iota*Z_theta_int) &
         &            *((Phi_zeta+iota*Phi_theta)**2.0) &
         &          +(R_zeta+iota*R_theta_int) & 
         &            *(Z_zeta_zeta+iota*iota*Z_theta_theta_int &
         &            +2.0*iota*Z_theta_zeta))


    kappa_g_top_3= & 
         &      ((Phi_theta*Z_zeta - Phi_zeta*Z_theta_int)) &
         &       *( R_int*R_int* & 
         &           (Z_zeta+iota*Z_theta_int)* & 
         &           (Phi_zeta_zeta+iota*iota*Phi_theta_theta & 
         &            +2.0*iota*Phi_theta_zeta) & 
         &        + 2.0*R_int* & 
         &           (R_zeta+iota*R_theta_int)* & 
         &           (Z_zeta+iota*Z_theta_int)* & 
         &           (Phi_zeta+iota*Phi_theta) &
         &        - R_int*R_int* & 
         &            (Z_zeta_zeta + iota*iota*Z_theta_theta_int &
         &              + 2.0*iota*Z_theta_zeta)* &
         &            (Phi_zeta + iota*Phi_theta))

    kg_t= (kappa_g_top_1+kappa_g_top_2 &
         &               +kappa_g_top_3) / ( (Q**(2.0))* &
         &                                        (Lcurv**(3.0)))

!    tau_n_top_1= & 
!         &      (R_theta_zeta*Phi_zeta+ & 
!         &            R_theta_int*Phi_zeta_zeta+ & 
!         &            iota*R_theta_theta_int*Phi_zeta+ & 
!         &            iota*R_theta_int*Phi_theta_zeta- & 
!         &            R_zeta_zeta*Phi_theta- & 
!         &            R_zeta*Phi_theta_zeta- & 
!         &            iota*R_theta_zeta*Phi_theta- & 
!         &            iota*R_zeta*Phi_theta_theta)* & 
!         &      (-(R_int**3.0)*(Phi_zeta+iota*Phi_theta)* &
!         &                        (Phi_theta*Z_zeta- & 
!         &                         Phi_zeta*Z_theta_int) & 
!         &        +R_int*(R_zeta+iota*R_theta_int)* & 
!         &                (Z_theta_int*R_zeta- & 
!         &                 Z_zeta*R_theta_int)) 

    tau_n_top_1 = -(R_int**3.0)*iota*R_theta_theta_int*Z_theta_int

    tau_n_top_2= &
         &      (Z_theta_zeta*R_zeta+ & 
         &            Z_theta_int*R_zeta_zeta+ & 
         &            iota*Z_theta_theta_int*R_zeta+ & 
         &            iota*Z_theta_int*R_theta_zeta- & 
         &            Z_zeta_zeta*R_theta_int- & 
         &            Z_zeta*R_theta_zeta- & 
         &            iota*Z_theta_zeta*R_theta_int- & 
         &            iota*Z_zeta*R_theta_theta_int)* &
         &      (R_int*(Z_zeta+iota*Z_theta_int)* & 
         &              (Phi_theta*Z_zeta - & 
         &               Phi_zeta*Z_theta_int) & 
         &      -R_int*(R_zeta+iota*R_theta_int)* &
         &              (R_theta_int*Phi_zeta- & 
         &               R_zeta*Phi_theta))

    tau_n_top_3= & 
         &      (Phi_theta_zeta*Z_zeta+ &
         &              Phi_theta*Z_zeta_zeta+ & 
         &              iota*Phi_theta_theta*Z_zeta+ &
         &              iota*Phi_theta*Z_theta_zeta- &
         &              Phi_zeta_zeta*Z_theta_int- & 
         &              Phi_zeta*Z_theta_zeta- & 
         &              iota*Phi_theta_zeta*Z_theta_int- & 
         &              iota*Phi_zeta*Z_theta_theta_int)* & 
         &      (-R_int*(Z_zeta+iota*Z_theta_int)* & 
         &               (Z_theta_int*R_zeta- & 
         &                Z_zeta*R_theta_int) & 
         &      +(R_int**3.0)*(Phi_zeta+iota*Phi_theta)* & 
         &               (R_theta_int*Phi_zeta- & 
         &                R_zeta*Phi_theta)) 


    tau_n_top_4= & 
         &      ((Phi_theta*Z_zeta- &
         &        Phi_zeta*Z_theta_int)**2.0)* &
         &      (R_int**2.0)*(Phi_zeta+iota*Phi_theta)* &
         &                      (Z_zeta+iota*Z_theta_int)+ &
         &      ((Z_theta_int*R_zeta- & 
         &        Z_zeta*R_theta_int)**2.0)* &
         &      (Phi_zeta+iota*Phi_theta)* & 
         &      (Z_zeta+iota*Z_theta_int)- & 
         &      (Phi_theta*Z_zeta- & 
         &       Phi_zeta*Z_theta_int)* & 
         &      (R_theta_int*Phi_zeta- & 
         &       R_zeta*Phi_theta)* & 
         &      (R_int**2.0)* & 
         &      (R_zeta+iota*R_theta_int)* & 
         &      (Phi_zeta+iota*Phi_theta)- &
         &      (R_theta_int*Phi_zeta- & 
         &       R_zeta*Phi_theta)* & 
         &      (Z_theta_int*R_zeta- & 
         &       Z_zeta*R_theta_int)* & 
         &      ( (R_int**2.0)*( (Phi_zeta+iota*Phi_theta) & 
         &         **2.0 ) - ( (R_zeta+iota*R_theta_int)**2.0))- &
         &      (Phi_theta*Z_zeta- & 
         &       Phi_zeta*Z_theta_int)* & 
         &      (Z_theta_int*R_zeta- & 
         &       Z_zeta*R_theta_int)* & 
         &      (R_zeta+iota*R_theta_int)* & 
         &      (Z_zeta+iota*Z_theta_int) 

    tn_t = (tau_n_top_1+tau_n_top_2+tau_n_top_3  &
         &         + tau_n_top_4)/( (Lcurv**2.0)*(Q**4.0))

    tau_n = tn_t
    kappa_g = kg_t

    DEALLOCATE(R_zeta,Z_zeta,Phi_zeta,Lcurv,Q,Phi_zeta_zeta,kappa_g_top_1, &
         &   kappa_g_top_2,kappa_g_top_3,tau_n_top_1,tau_n_top_2,tau_n_top_3, &
         &   tau_n_top_4, R_s_theta_theta, Z_s_theta_theta, Phi_theta, Phi_theta_theta, &
         &   Phi_theta_zeta, R_zeta_zeta, kn_t,kg_t,tn_t, R_theta_zeta, Z_zeta_zeta,&
         &   Z_theta_zeta ) 

    DEALLOCATE(R_theta_int,Z_theta_int,R_theta_theta_int,Z_theta_theta_int,R_int)

  END SUBROUTINE curvatures

  !> Miller with d/dr input and explicit metric tensor, following Candy PPCF 51 (2009) 105009
  subroutine get_miller_a(geom, C_y_out,C_xy_out,my_dpdx_pm,q_prof,dqdx_prof,Cyq0_x0,edge_opt)
    type(geomtype),intent(inout):: geom
    real, dimension(pi1gl:pi2gl), intent(out):: C_y_out, C_xy_out, my_dpdx_pm, q_prof,&
         dqdx_prof
    real, intent(out):: Cyq0_x0
    real, intent(in):: edge_opt
    integer,parameter:: np=500
    real, dimension(np):: R,Z,R_theta,Z_theta,R_theta_theta,Z_theta_theta,dlp,Rc,cosu,sinu,Bphi
    real, dimension(np):: theta, Bp_numerator, Bp_denom, Bp_n, Bp, D0, D1, D2, D3, nu, chi, theta_s
    real, dimension(np):: R_s, Z_s, R_theta_s, Z_theta_s, Rc_s, cosu_s, sinu_s, Bphi_s, B, B_s, Bp_s, dlp_s
    real, dimension(np):: nu1, nu1_s, dchidx, dchidx_s, dB_drho, dB_drho_s, dB_dl_s, dnu_drho_s, dnu_dl_s, grad_nu
    real, dimension(np):: gxx, gxy, gxz, gyy, gyz, gzz, dtheta_dchi_s, dBp_dchi_s, jacobian, dBdx, dBdz
    real, dimension(np):: g_xx, g_xy, g_xz, g_yy, g_yz, g_zz, tmp_arr, dxdR_s, dxdZ_s, K_x, K_y !tmp_arr2
    real, dimension(0:nz0-1):: gxx_out,gxy_out,gxz_out,gyy_out,gyz_out,gzz_out,Bfield_out,jacobian_out, dBdx_out, dBdz_out, chi_out
    real, dimension(0:nz0-1):: R_out, Z_out, dxdR_out, dxdZ_out
    real:: d_inv, drPsi, dxPsi, dq_dx, dq_dpsi, R0, Z0, B0, F, D0_full, D1_full, D2_full, D3_full
    !real :: Lnorm, Psi0 ! currently module-wide defined anyway
    real:: pprime, ffprime, D0_mid, D1_mid, D2_mid, D3_mid, dx_drho, pi, mu_0, dzprimedz
    integer:: i, k

    if (rho.lt.0.0) rho = trpeps*major_R
    if (rho.le.0.0) stop 'flux surface radius not defined'
    trpeps = rho/major_R
    x0=rho

    R0=major_R
    B0=1.0
    F=R0*B0
    Z0=major_Z
    pi = acos(-1.0)
    mu_0=4*pi

    theta=linspace(-pi,pi,np)
    d_inv=asin(delta)

    !flux surface parametrization
    R=R0+rho*cos(theta+d_inv*sin(theta))
    Z=Z0+kappa*rho*sin(theta+zeta*sin(2*theta))
    R_theta=-rho*(d_inv*cos(theta)+1)*sin(d_inv*sin(theta)+theta)
    Z_theta=rho*kappa*(2*zeta*cos(2*theta)+1)*cos(zeta*sin(2*theta)+theta)
    R_theta_theta=rho*d_inv*sin(theta)*sin(d_inv*sin(theta)+theta)-rho*(d_inv*cos(theta)+1)**2*cos(d_inv*sin(theta)+theta)
    Z_theta_theta=-rho*kappa*(2*zeta*cos(2*theta)+1)**2*sin(zeta*sin(2*theta)+theta)-&
         4*rho*kappa*zeta*sin(2*theta)*cos(zeta*sin(2*theta)+theta)

    !dl/dtheta
    dlp=(R_theta**2+Z_theta**2)**0.5
    !curvature radius
    Rc=dlp**3*(R_theta*Z_theta_theta-Z_theta*R_theta_theta)**(-1)
    ! some useful quantities (see papers for definition of u)
    cosu=Z_theta/dlp
    sinu=-R_theta/dlp
    Bphi=F/R

    !poloidal field
    Bp_numerator=sqrt(rho**2*kappa**2*(2*zeta*cos(2*theta)+1)**2*cos(zeta*sin(2*theta)+theta)**2+&
         rho**2*(d_inv*cos(theta)+1)**2*sin(d_inv*sin(theta)+theta)**2)

    Bp_denom=R*(rho*(d_inv*cos(theta)+1)*sin(d_inv*sin(theta)+theta)*(kappa*s_kappa*sin(zeta*sin(2*theta)+theta)&
         +kappa*sin(zeta*sin(2*theta)+theta)+kappa*s_zeta*sin(2*theta)*cos(zeta*sin(2*theta)+theta)+drZ)+&
         rho*kappa*(2*zeta*cos(2*theta)+1)*(-s_delta*sin(theta)*sin(d_inv*sin(theta)+theta)+&
         cos(d_inv*sin(theta)+theta)+drR)*cos(zeta*sin(2*theta)+theta))

    Bp_n=Bp_numerator/Bp_denom

    tmp_arr=1./Bp_n/R**2
    drPsi=F/2./pi/q0*dlp_int(tmp_arr,dlp)
    Bp=Bp_n*drPsi
    B=sqrt(Bphi**2+Bp**2)

!for comparison with miller_b
!    tmp_arr=sqrt((R-R0)**2+Z**2)
!    tmp_arr2=q0*R**2/F
!    Lnorm=average(tmp_arr,tmp_arr2)
!    Psi0=B0*Lnorm**2

    !definition of radial coordinate! dx_drho=1 --> x = r
    dx_drho=1. !drPsi/Psi0*Lnorm*q0
    if (mype==0.and.print_ini_msg) write(*,"(A,ES12.4)") 'Using radial coordinate with dx/dr = ',dx_drho
    dxPsi=drPsi/dx_drho
    C_y_out=dxPsi
    Cyq0_x0=C_y_out(pi1gl)*q0/rho
    C_xy_out=B0*dxPsi/C_y_out
    if (mype==0.and.print_ini_msg) then 
       write(*,"(A,ES12.4,A,ES12.4,A,ES12.4)") &
            "Setting C_xy = ",C_xy_out(pi1gl),' C_y = ', C_y_out(pi1gl)," C_x' = ", 1./dxPsi
       write(*,'(A,ES12.4)') "B_unit/Bref conversion factor = ", q0/rho*drPsi
    endif



    !--------shear is expected to be defined as rho/q*dq/drho--------!
    dq_dx=shat*q0/rho/dx_drho
    dq_dpsi=dq_dx/dxPsi
    pprime=-amhd/q0**2/R0/(2*mu_0)*B0**2/drPsi

    !neg. dpdx normalized to magnetic pressure for pressure term
    my_dpdx_pm=amhd/q0**2/R0/dx_drho

    !integrals for ffprime evaluation
    do i=1,np
       tmp_arr=(2./Rc/R/Bp-2*cosu/R**2/Bp)/R**2/Bp
       D0(i)=-F*dlp_int_ind(tmp_arr,dlp,i)
       tmp_arr=B**2/Bp**2/R**2/Bp
       D1(i)=-dlp_int_ind(tmp_arr,dlp,i)/F
       tmp_arr=mu_0/Bp**2/R**2/Bp
       D2(i)=-dlp_int_ind(tmp_arr,dlp,i)*F
       tmp_arr=1./R**2/Bp*F
       D3(i)=-dlp_int_ind(tmp_arr,dlp,i)
    enddo
    tmp_arr=(2./Rc/R/Bp-2*cosu/R**2/Bp)/R**2/Bp
    D0_full=-F*dlp_int(tmp_arr,dlp)
    tmp_arr=B**2/Bp**2/R**2/Bp
    D1_full=-dlp_int(tmp_arr,dlp)/F
    tmp_arr=mu_0/Bp**2/R**2/Bp
    D2_full=-dlp_int(tmp_arr,dlp)*F
    tmp_arr=1./R**2/Bp*F
    D3_full=-dlp_int(tmp_arr,dlp)
    D0_mid=0.5*(D0(np/2)+D0(np/2+1))
    D1_mid=0.5*(D1(np/2)+D1(np/2+1))
    D2_mid=0.5*(D2(np/2)+D2(np/2+1))
    D3_mid=0.5*(D3(np/2)+D3(np/2+1))
    ffprime=-(dq_dpsi*2.*pi+D0_full+D2_full*pprime)/D1_full

    D0=D0-D0_mid
    D1=D1-D1_mid
    D2=D2-D2_mid
    nu=D3-D3_mid
    !straight field line angle defined on equidistant theta grid
    chi=-nu/q0
    !correct small scaling error (<0.5%, due to finite integration resolution)
    chi=chi*2*pi/(maxval(chi)-minval(chi))
    
    call lag3interp(theta,chi,np,theta_s,theta,np)
    dtheta_dchi_s=deriv_fd(theta_s,theta)

    !arrays equidistant in straight field line angle
    R_s=R0+rho*cos(theta_s+d_inv*sin(theta_s))
    Z_s=Z0+kappa*rho*sin(theta_s+zeta*sin(2*theta_s))
    R_theta_s=-rho*(d_inv*cos(theta_s)+1)*sin(d_inv*sin(theta_s)+theta_s)*dtheta_dchi_s
    Z_theta_s=rho*kappa*(2*zeta*cos(2*theta_s)+1)*cos(zeta*sin(2*theta_s)+theta_s)*dtheta_dchi_s

    call lag3interp(Bp,theta,np,Bp_s,theta_s,np)
    call lag3interp(dlp,theta,np,dlp_s,theta_s,np)
    call lag3interp(Rc,theta,np,Rc_s,theta_s,np)
    dBp_dchi_s=deriv_fd(Bp_s,theta)
    Bphi_s=F/R_s
    B_s=sqrt(Bphi_s**2+Bp_s**2)
    cosu_s=Z_theta_s/dlp_s/dtheta_dchi_s
    sinu_s=-R_theta_s/dlp_s/dtheta_dchi_s

    nu1=R*Bp*(D0+D1*ffprime+D2*pprime)
    call lag3interp(nu1,theta,np,nu1_s,theta_s,np)
    
    !radial derivative of straight field line angle
    dchidx=-1./q0*nu1/R/Bp*dxPsi+nu/q0**2*dq_dx
    dchidx_s=-1./q0*nu1_s/R_s/Bp_s*dxPsi-theta/q0*dq_dx

    !Bfield derivatives in Mercier-Luc coordinates
    dB_drho=-1/B*(F**2/R**3*cosu+Bp**2/Rc+mu_0*R*Bp*pprime)
    dB_drho_s=-1/B_s*(F**2/R_s**3*cosu_s+Bp_s**2/Rc_s+mu_0*R_s*Bp_s*pprime)
    dB_dl_s=1./B_s*(Bp_s*dBp_dchi_s/dtheta_dchi_s/dlp_s+F**2/R_s**3*sinu_s)

    dnu_drho_s=nu1_s
    dnu_dl_s=-F/R_s**2/Bp_s
    grad_nu=sqrt(dnu_drho_s**2+dnu_dl_s**2)

    !contravariant metric coefficients
    gxx=R_s**2*Bp_s**2/dxPsi**2
    gxy=-R_s*Bp_s/dxPsi*C_y_out(pi1gl)*nu1_s
    gxz=-R_s*Bp_s/dxPsi/q0*nu1_s-R_s**2*Bp_s**2/dxPsi/q0*dq_dpsi*theta
    gyy=C_y_out(pi1gl)**2*(grad_nu**2+1/R_s**2)
    gyz=C_y_out(pi1gl)/q0*grad_nu**2+C_y_out(pi1gl)/q0*dq_dpsi*nu1_s*R_s*Bp_s*theta
    gzz=1/q0**2*grad_nu**2+2/q0**2*dq_dpsi*nu1_s*R_s*Bp_s*theta+1/q0**2*dq_dpsi**2*R_s**2*Bp_s**2*theta**2

    jacobian=1./sqrt(gxx*gyy*gzz + 2.*gxy*gyz*gxz - gxz**2*gyy - gyz**2*gxx - gzz*gxy**2)
    !covariant metric coefficients
    g_xx=jacobian**2*(gyy*gzz-gyz**2)
    g_xy=jacobian**2*(gxz*gyz-gxy*gzz)
    g_xz=jacobian**2*(gxy*gyz-gxz*gyy)
    g_yy=jacobian**2*(gxx*gzz-gxz**2)
    g_yz=jacobian**2*(gxz*gxy-gxx*gyz)
    g_zz=jacobian**2*(gxx*gyy-gxy**2)
    
    !Bfield derivatives
    dBdx=jacobian*C_y_out(pi1gl)/q0*(F/R_s**3/Bp_s*dB_drho_s+(nu1_s/R_s+dq_dpsi*theta*Bp_s)*dB_dl_s)
    dBdz=1./B_s*(Bp_s*dBp_dchi_s-F**2/R_s**3*R_theta_s)
    !curvature terms (these are just local and will be recalculated in geometry.F90)
    K_x = (0.-g_yz/g_zz*dBdz)
    K_y = (dBdx-g_xz/g_zz*dBdz)

    !(R,Z) derivatives for visualization
    dxdR_s = dx_drho/drPsi*R_s*Bp_s*cosu_s
    dxdZ_s = dx_drho/drPsi*R_s*Bp_s*sinu_s

    if (edge_opt==0.0) then
       !gene z-grid
       chi_out=linspace(-pi,pi-2*pi/nz0,nz0)
    else
       !new parallel coordinate chi_out==zprime
       !see also tracer_aux.F90
       do k=0,nz0-1
          chi_out(k)=sinh((-pi+k*2.*pi/nz0)*log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/pi)/edge_opt
       enddo
       !transform metrics according to chain rule
       do k=1,np
          dzprimedz=edge_opt*pi/log(edge_opt*pi+sqrt(edge_opt**2*pi**2+1))/&
               sqrt(edge_opt**2*theta(k)**2+1)
          gxz(k)=gxz(k)*dzprimedz
          gyz(k)=gyz(k)*dzprimedz
          gzz(k)=gzz(k)*dzprimedz**2
          jacobian(k)=jacobian(k)/dzprimedz
          dBdz(k)=dBdz(k)/dzprimedz
       enddo
    endif !edge_opt

    !interpolate down to GENE z-grid
    call lag3interp(gxx,theta,np,gxx_out,chi_out,nz0)
    call lag3interp(gxy,theta,np,gxy_out,chi_out,nz0)
    call lag3interp(gxz,theta,np,gxz_out,chi_out,nz0)
    call lag3interp(gyy,theta,np,gyy_out,chi_out,nz0)
    call lag3interp(gyz,theta,np,gyz_out,chi_out,nz0)
    call lag3interp(gzz,theta,np,gzz_out,chi_out,nz0)
    call lag3interp(B_s,theta,np,Bfield_out,chi_out,nz0)
    call lag3interp(jacobian,theta,np,jacobian_out,chi_out,nz0)
    call lag3interp(dBdx,theta,np,dBdx_out,chi_out,nz0)
    call lag3interp(dBdz,theta,np,dBdz_out,chi_out,nz0)
    call lag3interp(R_s,theta,np,R_out,chi_out,nz0)
    call lag3interp(Z_s,theta,np,Z_out,chi_out,nz0)
    call lag3interp(dxdR_s,theta,np,dxdR_out,chi_out,nz0)
    call lag3interp(dxdZ_s,theta,np,dxdZ_out,chi_out,nz0)

    !select local k range
    do i=pi1gl,pi2gl
       if (yx_order) then
          geom%gii(i,pj1,lk1:lk2)=gyy_out(lk1:lk2)
          geom%gjj(i,pj1,lk1:lk2)=gxx_out(lk1:lk2)
          geom%giz(i,pj1,lk1:lk2)=gyz_out(lk1:lk2)
          geom%gjz(i,pj1,lk1:lk2)=gxz_out(lk1:lk2)
          geom%dBdi(i,pj1,lk1:lk2)=0.
          geom%dBdj(i,pj1,lk1:lk2)=dBdx_out(lk1:lk2)
       else
          geom%gii(i,pj1,lk1:lk2)=gxx_out(lk1:lk2)
          geom%gjj(i,pj1,lk1:lk2)=gyy_out(lk1:lk2)
          geom%giz(i,pj1,lk1:lk2)=gxz_out(lk1:lk2)
          geom%gjz(i,pj1,lk1:lk2)=gyz_out(lk1:lk2)
          geom%dBdi(i,pj1,lk1:lk2)=dBdx_out(lk1:lk2)
          geom%dBdj(i,pj1,lk1:lk2)=0.
       endif
       geom%gij(i,pj1,lk1:lk2)=gxy_out(lk1:lk2)
       geom%gzz(i,pj1,lk1:lk2)=gzz_out(lk1:lk2)
       geom%Bfield(i,pj1,lk1:lk2)=Bfield_out(lk1:lk2)
       geom%jacobian(i,pj1,lk1:lk2)=jacobian_out(lk1:lk2)
       geom%dBdz(i,pj1,lk1:lk2)=dBdz_out(lk1:lk2)
       geom%R(i,lk1:lk2)=R_out(lk1:lk2)*Lref
       geom%Z(i,lk1:lk2)=Z_out(lk1:lk2)*Lref
       geom%R_hat(i,lk1:lk2)=R_out(lk1:lk2)
       geom%Z_hat(i,lk1:lk2)=Z_out(lk1:lk2)
       geom%dxdR(i,lk1:lk2)= dxdR_out(lk1:lk2)
       geom%dxdZ(i,lk1:lk2)= dxdZ_out(lk1:lk2)
    enddo
    q_prof=q0
    dqdx_prof=shat*q0/x0
    
#if 0
    if (mype==0) then
       write (*,'(A,F12.6)') 'C_xy from metric: ',geom%Bfield(pi1gl,pj1,lk1)/&
            sqrt(geom%gii(pi1,pj1,lk1)*geom%gjj(pi1,pj1,lk1)-geom%gij(pi1gl,pj1,lk1)**2.0)
       write (*,'(A,F12.6)') 'C_xy output     : ',C_xy_out(pi1gl)
    endif
#endif
  contains

    !> Generate an equidistant array from min to max with n points
    function linspace(min,max,n) result(out)
      real:: min, max
      integer:: n
      real, dimension(n):: out
      
      do i=1,n
         out(i)=min+(i-1)*(max-min)/(n-1)
      enddo
    end function linspace

    !> Weighted average
    real function average(var,weight)
      real, dimension(np):: var, weight
      average=sum(var*weight)/sum(weight)
    end function average

    !> full theta integral with weight function dlp
    real function dlp_int(var,dlp)
      real, dimension(np):: var, dlp
      dlp_int=sum(var*dlp)*2*pi/(np-1)
    end function dlp_int

    !> theta integral with weight function dlp, up to index 'ind'
    real function dlp_int_ind(var,dlp,ind)
      real, dimension(np):: var, dlp
      integer:: ind

      dlp_int_ind=0.
      if (ind.gt.1) then
         dlp_int_ind=dlp_int_ind+var(1)*dlp(1)*pi/(np-1)
         dlp_int_ind=dlp_int_ind+(sum(var(2:ind-1)*dlp(2:ind-1)))*2*pi/(np-1)
         dlp_int_ind=dlp_int_ind+var(ind)*dlp(ind)*pi/(np-1)
      endif
    end function dlp_int_ind

    !> 1st derivative with 2nd order finite differences
    function deriv_fd(y,x) result(out)
      real, dimension(np):: x,y,out,dx

      out=0.
      do i=2,np-1
         out(i)=out(i)-y(i-1)/2
         out(i)=out(i)+y(i+1)/2
      enddo
      out(1)=y(2)-y(1)
      out(np)=y(np)-y(np-1)
      dx=x(2)-x(1)
      out=out/dx
    end function deriv_fd


  end subroutine get_miller_a
  

END MODULE miller_mod
