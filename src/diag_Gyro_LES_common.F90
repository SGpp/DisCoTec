#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
!Routines used by the various Gyro_LES diags
Module diag_Gyro_LES_common
  
  Use par_mod
  use antenna, only: antenna_type
  Use vel_space, only: fm, mat_00
  use geometry
  use communications
!  use coordinates
  use aux_fields
  Use gyro_average_ff_mod, only: gyro_average_ff, gyro_average_ff_bpar
!  Use gyro_average_df_mod, only: gyro_average_df
!  use blockindex
  use aux_func, only: gamma0

  Implicit None

  public:: filter, get_cfgamma, get_cfgamma_mod, get_electrostatic
  public:: integral_vp, integral_w, integral_vpw
  public:: calc_n_shells_ky, calc_n_shells_log_kperp_1, calc_n_shells_log_kperp, calc_n_shells_kx

  !***** Common parameters to all GyroLES diags *****

  Complex, dimension(:,:,:,:,:,:),allocatable :: g_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable :: gc_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable :: d_v_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable :: d_z_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable  :: d_kperp_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable  :: coll_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable  :: c_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable  :: p_rhs
  Complex, dimension(:,:,:,:,:,:),allocatable  :: cfgamma, cfgamma1, cfgamma2, cfgamma_ant, cfgamma_ant_ls1
  Complex, Dimension(:,:,:,:), Allocatable :: cfgamma3
  Complex, Dimension(:,:,:,:,:,:,:), Allocatable :: v7d_spec
  real, dimension(:,:,:), allocatable :: kperp

  ! Allocatable variables for the 3d transfer
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: g_ls1
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: g_ls2 
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: g_ls3
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_mod_ls1
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_mod_ls2
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_mod_ls3
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_ls2_ls3_g
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_ls2_ls3_f
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_ls2_ls3_h
  Complex, Dimension(:,:,:,:,:,:,:), Allocatable :: fe_tf_term
  Complex, Dimension(:,:,:,:,:,:,:), Allocatable :: fe_tf_term_p
  Complex, Dimension(:,:,:,:,:,:,:), Allocatable :: fe_tf_term_n
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: f_ls1
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: f_ls2
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: f_ls3
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_ls1
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_ls2
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_ls3
  Complex, Dimension(:,:,:,:), Allocatable :: fields_ls1
  Complex, Dimension(:,:,:,:), Allocatable :: fields_ls2
  Complex, Dimension(:,:,:,:), Allocatable :: fields_ls3

  integer :: istep_fe_twoD = 0, avg_window = 1
  logical :: tky = .false.
  logical :: tkx = .false.
  logical :: tkperp = .false.
  logical :: t_log_kperp = .false.
  logical :: t_log_kperp_h = .false.
  logical :: t_log_kperp_fb = .false.
  logical :: t_log_kperp_1 = .false.
  logical :: t_log_ky = .false.
  logical :: tk3_log_kperp = .false.
  logical :: tk3_log_kperp_1 = .false.
  logical :: tk3_log_kperp_h = .false.
  logical :: tk3_ky = .false.
  logical :: tk3_kx = .false.
  integer :: n_shells_ky = 25 ! For security in the allocatable variable
  integer :: n_shells_kx = 25 ! For security in the allocatable variable
  integer :: n_shells_kperp = 25 ! For security in the allocatable variable
  integer :: n_shells_log_kperp = 25 ! For security in the allocatable variable
  integer :: n_shells_log_ky = 25 ! For security in the allocatable variable
  integer :: n_shells_log_kperp_1 = 25 ! For security in the allocatable variable
  logical :: with_nvp = .false.
  real ::  delta_ky, delta_kx, delta_kperp, lambda1, lambda2
  real :: kperp_max, k0, k0_y
  integer :: n_transfer_models = 7

contains

  subroutine init_kperp
    real :: kperp_maxloc
    integer :: i, j, k, ierr

    if (.not.allocated(kperp)) allocate(kperp(li1:li2,lj1:lj2,lk1:lk2))

    Do k=lk1,lk2
       Do j=lj1,lj2
          Do i=li1,li2
             kperp(i,j,k) = sqrt(geom%gii(pi1,pj1,k)*kx(i)**2+2.0*geom%gij(pi1,pj1,k)*&
                  &kx(i)*ky(j)+geom%gjj(pi1,pj1,k)*ky(j)**2)
          Enddo
       Enddo
    Enddo

    kperp_maxloc = maxval(kperp)
    call mpi_allreduce(kperp_maxloc,kperp_max,1,MPI_REAL_TYPE,MPI_MAX,mpi_comm_xyz,ierr)

    ! Initialization parameters for the shells and transfer quantities
    n_shells_ky = nky0
    delta_ky = kymin
    n_shells_kx = real(nx0/2)
    delta_kx = real(2.*pi/lx)
    delta_kperp = real(kperp_max)/real(n_shells_kperp)
    k0 = kperp_max/(2**(real(n_shells_log_kperp-1)/5.))
    k0_y = kymax/(2**(real(n_shells_log_ky-1)/6.))

  end subroutine init_kperp

!!!**************************************************************************!!!
!!!*********************** Calculate the filter  ************************
  subroutine filter(f_in,f_out)

    implicit none

    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: f_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: f_out
    Integer :: fkx, fky

    ! Filter in the kx direction, it gives you the indice of the maximum value of the kx mode you keep
    fkx = int((fracx*hkx))
    ! Filter in the ky direction, it gives you the indice of the maximum value of the ky mode you keep 
    !fky = int((fracy*lj2))
    fky = int(fracy*(nky0-1))

    !Copy the input function to the output function FGS
    f_out = f_in
    !Filter FGS in the kx direction for modes higher than fkx and lower than -fkx 

    f_out( li1 + fkx + 1 : li2 - fkx,:,:,:,:,:) = CMPLX(0.,0.)

    IF (fky.ge.lj1 .and. fky.ge.lj2) then
    ELSE IF (fky.ge.lj1.and.fky.lt.lj2) then
       f_out(:,fky + 1 : lj2,:,:,:,:) = CMPLX(0.,0.)
    ELSE IF (fky.ge.lj1.and.fky.eq.lj2) then
    ELSE
       f_out(:,:,:,:,:,:)= CMPLX(0.,0.)
    ENDIF

  End Subroutine filter

!!!*********************************************************************************

  SUBROUTINE get_cfgamma(f_in,fields_in,cfgamma1,cfgamma2)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: fields_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: f_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma2,cfgamma1
    complex, dimension(li1:li2, lj1:lj2):: phi_bar, bpar_bar
    integer :: k,l,m,n,pni
    pni=pn1
    DO n=ln1,ln2
       if (pn0.gt.1) pni=n
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                CALL gyro_average_ff(fields_in(:,:,k,1),phi_bar,k,m,n)
                cfgamma2(:,:,k,l,m,n)=spec(n)%temp*spec(n)%dens*&
                     (spec(n)%charge*CONJG(phi_bar(:,:))/spec(n)%temp)
                if (n_fields.gt.2) then
                   !bpar contribution
                   CALL gyro_average_ff_bpar(fields_in(:,:,k,3),bpar_bar,k,m,n)
                   cfgamma2(:,:,k,l,m,n) = cfgamma2(:,:,k,l,m,n)+spec(n)%temp*spec(n)%dens*mu(m)*conjg(bpar_bar)
                endif
                cfgamma1(:,:,k,l,m,n)=spec(n)%temp*spec(n)%dens/fm(pi1,pj1,k,l,m,pni)*&
                     (CONJG(f_in(:,:,k,l,m,n)))
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE get_cfgamma

  !WARNING: f_in and fields_in are not used in this subroutine
  !does this make sense?????
  SUBROUTINE get_cfgamma_mod(f_in,fields_in,cfgamma1m)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: fields_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: f_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out):: cfgamma1m
!    complex, dimension(li1:li2, lj1:lj2):: phi_bar
    integer :: k,l,m,n,pni
    pni=pn1
    DO n=ln1,ln2
       if (pn0.gt.1) pni=n
       DO m=lm1,lm2
          DO l=ll1,ll2
             DO k=lk1,lk2
                !not used...
                !CALL gyro_average_ff(fields_in(:,:,k,1),phi_bar,k,m,n)
                cfgamma1m(:,:,k,l,m,n)=fm(pi1,pj2,k,l,m,pni)/(spec(n)%temp*spec(n)%dens) 
             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE get_cfgamma_mod


  SUBROUTINE get_electrostatic(fields_in,cfgamma3)
    Implicit None
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: fields_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2), intent(out):: cfgamma3
    Real :: b2spec, k2, b2
    integer :: i,j,k,n

    Do k=lk1,lk2
       Do j = lj1,lj2
          Do i = li1,li2
!\todo could be precomputed during initialization
             k2 = geom%gii(pi1,pj1,k)*ki(i)**2+2.0*geom%gij(pi1,pj1,k)*&
                  & ki(i)*kj(j)+geom%gjj(pi1,pj1,k)*kj(j)**2
             b2 = k2/(geom%Bfield(pi1,pj1,k)**2)

!\todo not performant: innermost loop over last index
             do n=ln1,ln2
                b2spec = b2*spec(n)%mass*spec(n)%temp/(spec(n)%charge**2)
                cfgamma3(i,j,k,n)=(spec(n)%dens*spec(n)%charge**2&
                     /(2.0*spec(n)%temp))*(1.0 - gamma0(b2spec))*&
                     CONJG(fields_in(i,j,k,1))*fields_in(i,j,k,1)
             enddo
          enddo
       enddo
    enddo

  END SUBROUTINE get_electrostatic

!!!**************************************************************************!!!
!!!********************* Subroutine that the module needs ************************!!! 
!!!**************************************************************************!!!

  subroutine integral_vpw(v5d,v3d)
    !This subroutine performs the mat_00 integral and doesnt have explicit species
    !dependence

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(out) :: v3d
    integer :: k,m,l

    v3d=cmplx(0.0,0.0)
    do k=lk1,lk2
       do m=lm1,lm2
          do l=ll1,ll2
             call saxpy(2*li0*lj0,mat_00(pi1,pj1,k,l,m),v5d(li1,lj1,k,l,m),1,v3d(li1,lj1,k),1)
          enddo
       enddo
    enddo

    call my_complex_sum_vw(v3d,size(v3d)) 


  end subroutine integral_vpw

  subroutine integral_vp(v5d,v4d)
    !This subroutine performs the v parallel integral and doesnt have explicit species
    !dependence

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2), intent(in) :: v5d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2),intent(out) :: v4d
    integer :: m,l

 
    v4d=cmplx(0.0,0.0)
    
    do m=lm1,lm2
        do l=ll1,ll2
              call saxpy(2*li0*lj0*lk0,vp_weight(l),v5d(li1,lj1,lk1,l,m),1,v4d(li1,lj1,lk1,m),1)
        enddo
    enddo

    call my_complex_sum_v(v4d,size(v4d)) 

  end subroutine integral_vp

  subroutine integral_w(v4d,v3d)
    !This subroutine perform the mu integral, which it is
    !2*pi v_perp d_vperp = pi \hat{B} dmu,  and doesnt have explicit species
    !dependence

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2), intent(in) :: v4d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2),intent(out) :: v3d
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2,lm1:lm2) :: mat_00_mu
    integer :: k,m,j

    v3d=cmplx(0.0,0.0)
    
    Do m=lm1,lm2
        Do k=lk1,lk2
            Do j=pj1,pj2
                mat_00_mu(pi1:pi2,j,k,m) = pi*mu_weight(m)*&
                    &geom%Bfield(pi1:pi2,j,k)
            Enddo
        Enddo
    Enddo

    do m=lm1,lm2
       do k=lk1,lk2
             call saxpy(2*li0*lj0,mat_00_mu(pi1,pj1,k,m),v4d(li1,lj1,k,m),1,v3d(li1,lj1,k),1)
       enddo
    enddo

    call my_complex_sum_w(v3d,size(v3d)) 


  end subroutine integral_w

!!!******************************************************************************************

  Subroutine calc_n_shells_ky(ls1,g_1,g_ls1)
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    integer :: j
    g_ls1=cmplx(0.0,0.0)

    do j=lj1,lj2
       if (ls1.eq.1) then
          if ((ky(j).ge.(0)).and.(ky(j).lt.(delta_ky*ls1))) then
             g_ls1(:,j,:,:,:,:) = g_1(:,j,:,:,:,:)
          endif
       else 
          if ((ky(j).ge.(delta_ky*(ls1-1))).and.(ky(j).lt.(delta_ky*ls1))) then
             g_ls1(:,j,:,:,:,:) = g_1(:,j,:,:,:,:)
          endif
       endif
    enddo

  End Subroutine calc_n_shells_ky

!!!******************************************************************************************

  Subroutine calc_n_shells_kx(ls1,g_1,g_ls1)
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    integer :: i

    !Calculate the shell splitting
    !delta_kx of the shells
    g_ls1=cmplx(0.0,0.0)
    do i=li1,li2
       if (ls1.eq.1) then
          if ((abs(kx(i)).eq.(0))  ) then
             g_ls1(i,:,:,:,:,:) = g_1(i,:,:,:,:,:)
          endif
       else if ((ls1.gt.1).and.(ls1.le.(n_shells_kx)) ) then
          if ((abs(kx(i)).ge.(delta_kx*(ls1-1))).and.(abs(kx(i)).lt.(delta_kx*ls1))  ) then
             g_ls1(i,:,:,:,:,:) = g_1(i,:,:,:,:,:)
          endif
       endif
    enddo
  End Subroutine calc_n_shells_kx


!!!******************************************************************************************

  Subroutine calc_n_shells_kperp(ls1,g_1,g_ls1)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    integer :: i,j,k

    !Calculate the shell splitting
    !delta_k of the shells       
    g_ls1=cmplx(0.0,0.0)
    do k=lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
             if (ls1.eq.1) then
                if ((kperp(i,j,k).ge.(0)).and.(kperp(i,j,k).le.(delta_kperp*real(ls1)))) then
                   g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                endif
             else if ((ls1.gt.1).and.(ls1.le.(n_shells_kperp)) ) then
                if ((kperp(i,j,k).gt.(real(ls1-1)*delta_kperp)).and.(kperp(i,j,k).le.(delta_kperp*real(ls1)))) then
                   g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                endif
             endif
          enddo
       enddo
    enddo

  End Subroutine calc_n_shells_kperp

!!!******************************************************************************************

  Subroutine calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    real, intent(in) :: k0
    integer :: i,j,k

    !Calculate the shell splitting
    !delta_k of the shells       
    g_ls1=cmplx(0.0,0.0)
    do k=lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
            if (ls1.eq.1) then
                if ((kperp(i,j,k).ge.(0)).and.(kperp(i,j,k).le.(k0*2**real((ls1-1)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                endif
            else if ((ls1.gt.1).and.(ls1.lt.n_shells_log_kperp)) then
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-2)/5.))).and.(kperp(i,j,k).le.(k0*2**real((ls1-1)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                endif
            else
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-2)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                endif
            endif
         enddo
       enddo
    enddo

  End Subroutine calc_n_shells_log_kperp

!**********************************************************************************************

  Subroutine calc_n_shells_log_kperp_pn(k0,ls1,g_1,g_ls1_p,g_ls1_n)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2),intent(out) :: g_ls1_p, g_ls1_n
    integer, intent(in) :: ls1
    real, intent(in) :: k0
    integer :: i,j,k,l,m

    !Calculate the shell splitting
    !delta_k of the shells       
    g_ls1_p=cmplx(0.0,0.0)
    g_ls1_n=cmplx(0.0,0.0)
    do k=lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
            if (ls1.eq.1) then
                if ((kperp(i,j,k).ge.(0)).and.(kperp(i,j,k).le.(k0*2**real((ls1-1)/5.)))) then
                    do m=lm1,lm2
                        do l=ll1,ll2
                            if (real(g_1(i,j,k,l,m)).ge.0) then
                                g_ls1_p(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            else 
                                g_ls1_n(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            endif
                        enddo
                    enddo
                endif
            else if ((ls1.gt.1).and.(ls1.lt.n_shells_log_kperp)) then
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-2)/5.))).and.(kperp(i,j,k).le.(k0*2**real((ls1-1)/5.)))) then
                    do m=lm1,lm2
                        do l=ll1,ll2
                            if (real(g_1(i,j,k,l,m)).ge.0) then    
                                g_ls1_p(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            else 
                                g_ls1_n(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            endif
                        enddo
                    enddo  
                endif    
            else
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-2)/5.)))) then
                    do m=lm1,lm2
                        do l=ll1,ll2
                            if (real(g_1(i,j,k,l,m)).ge.0) then
                                g_ls1_p(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            else 
                                g_ls1_n(i,j,k,l,m) = cmplx(real(g_1(i,j,k,l,m)),0.0)
                            endif
                        enddo
                    enddo
                endif
            endif
         enddo
       enddo
    enddo

  End Subroutine calc_n_shells_log_kperp_pn


!!!******************************************************************************************

Subroutine calc_n_shells_log_kperp_1(ls1,g_1,g_ls1)
    
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    integer :: i,j,k
    
    !Calculate the shell splitting
    !Lambda 1 and Lambda 2, define the size of the first 2 shells.
    !delta_k of the shells (excluding the first two shells)        
    g_ls1=cmplx(0.0,0.0)
        do k=lk1,lk2
            do j=lj1,lj2
                do i=li1,li2
                    if (ls1.eq.1) then
                        if ((kperp(i,j,k).ge.(0)).and.(kperp(i,j,k).le.(lambda1))) then
                            g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                        endif
                    else if (ls1.eq.2) then
                        if ((kperp(i,j,k).gt.(lambda1)).and.(kperp(i,j,k).le.(lambda1+lambda2))) then
                             g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                        endif
                    else if (ls1.eq.3) then
                        if ((kperp(i,j,k).gt.(lambda1+lambda2)).and.(kperp(i,j,k).le.((lambda1+lambda2)*2**real((ls1-2)/5.)))) then
                            g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                        endif
                    else if ((ls1.gt.3).and.(ls1.lt.n_shells_log_kperp_1)) then
                        if ((kperp(i,j,k).gt.((lambda1+lambda2)*2**real((ls1-3)/5.))).and.&
                            (kperp(i,j,k).le.((lambda1+lambda2)*2**real((ls1-2)/5.)))) then
                            g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                        endif
                    else
                        if ((kperp(i,j,k).gt.((lambda1+lambda2)*2**real((ls1-3)/5.)))) then
                            g_ls1(i,j,k,:,:,:) = g_1(i,j,k,:,:,:)
                        endif 
                    endif
                enddo
            enddo
        enddo   
            
End Subroutine calc_n_shells_log_kperp_1


!!!******************************************************************************************
  Subroutine calc_n_shells_log_ky(k0_y,ls1,g_1,g_ls1)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    integer, intent(in) :: ls1
    real, intent(in) :: k0_y
    integer :: j

    !Calculate the shell splitting
    !delta_k of the shells       
    g_ls1=cmplx(0.0,0.0)
    do j=lj1,lj2
        if (ls1.eq.1) then
            if ((ky(j).ge.(0)).and.(ky(j).le.(k0_y*2**real((ls1-1)/6.)))) then
                g_ls1(:,j,:,:,:,:) = g_1(:,j,:,:,:,:)
            endif
        else if ((ls1.gt.1).and.(ls1.lt.n_shells_log_ky)) then
            if ((ky(j).gt.(k0_y*2**real((ls1-2)/6.))).and.(ky(j).le.(k0_y*2**real((ls1-1)/6.)))) then
                g_ls1(:,j,:,:,:,:) = g_1(:,j,:,:,:,:)
            endif
        else
            if ((ky(j).gt.(k0_y*2**real((ls1-2)/6.)))) then
                g_ls1(:,j,:,:,:,:) = g_1(:,j,:,:,:,:)
            endif
        endif
    enddo

  End Subroutine calc_n_shells_log_ky


!!!******************************************************************************************

  Subroutine calc_n_shells_log_kperp_zonal(k0,ls1,g_1,g_ls1)
    !check why the 2 procedures dont give the same result
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(out) :: g_ls1
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)  :: g_1_phi
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)  :: g_1_del
    complex, dimension(li1:li2,ll1:ll2, lm1:lm2, ln1:ln2)  :: g_phi
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2) :: f_del_1
    complex, dimension(li1:li2)  :: fields_phi
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz)  :: fields_del
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields)  :: fields_del_1
    integer, intent(in) :: ls1
    real, intent(in) :: k0
    integer :: i,j,k,l,m,n

    !Calculate the shell splitting
    !delta_k of the shells       
    g_ls1=cmplx(0.0,0.0)
    g_1_phi = cmplx(0.0,0.0)
    g_phi = cmplx(0.0,0.0)
    fields_phi = cmplx(0.0,0.0)
    fields_del(:,:,:) = emfields(:,:,:,1)
    f_del_1 = cmplx(0.0,0.0)
    
    g_1_del = g_1
    if (p_has_0_mode) then
        do  n=ln1,ln2
            do m = lm1,lm2
                do l =ll1,ll2
                   do i=li1,li2
                      g_phi(i,l,m,n) = Sum(g_1(i,lj1,lk1:lk2,l,m,n)*geom%jacobian(pi1,pj1,lk1:lk2)) / (nz0*geom%avg_jaco)
                   Enddo
                Enddo
            Enddo
        enddo
       do n=ln1,ln2
            do m =lm1,lm2
                do l =ll1,ll2
                    Call sum_to_all_complex(g_phi(:,l,m,n), Size(g_phi(:,l,m,n)), mpi_comm_z)
                enddo
            enddo
        enddo
        Do n=ln1,ln2
            do m=lm1,lm2
                do l=ll1,ll2         
                    Do k=lk1,lk2
                        g_1_del(:,lj1,k,l,m,n) = g_1_del(:,lj1,k,l,m,n) - g_phi(:,l,m,n)
                        g_1_phi(:,lj1,k,l,m,n) = g_phi(:,l,m,n)
                    End Do
                enddo
            enddo
        enddo
  
        do i = li1,li2
            fields_phi(i) = Sum(emfields(i,lj1,lk1:lk2,1)*geom%jacobian(pi1,pj1,lk1:lk2)) / (nz0*geom%avg_jaco)
        enddo
        Call sum_to_all_complex(fields_phi, Size(fields_phi), mpi_comm_z)
       do k=lk1,lk2
            fields_del(:,lj1,k) = fields_del(:,lj1,k) - fields_phi
        End Do
    endif

    call calc_aux_fields(g_1_del, fields_del_1, f_del_1,.false.)
    ! Compare fields_del, with fields_del_1
!    counter = 0
!    if (counter.eq.0) then
!    Do k=lk1,lk2
!        do j=lj1,lj2
!            do i=li1,li2
!        if (mype.eq.0.and.j.eq.0) write(*,*) "checking zonal flow components",&
!             & i,j,k,fields_del_1(i,j,k,1) + fields_phi(i) - emfields(i,j,k,1),&
!             & fields_del(i,j,k) - emfields(i,j,k,1) + fields_phi(i)
!if (mype.eq.0.and.j.ne.0) write(*,*) "checking non zonal flow components",&
!             & i,j,k,fields_del_1(i,j,k,1) - emfields(i,j,k,1), &
!             & fields_del(i,j,k) -emfields(i,j,k,1)
!            enddo
!        enddo
!    enddo
!    counter = 1
!    endif


    do k=lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
            if (ls1.eq.1) then
                g_ls1(i,j,k,:,:,:) = g_1_phi(i,j,k,:,:,:)
            else if (ls1.eq.2) then
                if ((kperp(i,j,k).ge.(0)).and.(kperp(i,j,k).le.(k0*2**real((ls1-2)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1_del(i,j,k,:,:,:)
                endif
            else if ((ls1.gt.2).and.(ls1.lt.n_shells_log_kperp)) then
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-3)/5.))).and.(kperp(i,j,k).le.(k0*2**real((ls1-2)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1_del(i,j,k,:,:,:)
                endif
            else
                if ((kperp(i,j,k).gt.(k0*2**real((ls1-3)/5.)))) then
                    g_ls1(i,j,k,:,:,:) = g_1_del(i,j,k,:,:,:)
                endif
            endif
         enddo
       enddo
    enddo

  End Subroutine calc_n_shells_log_kperp_zonal


end Module diag_Gyro_LES_common
