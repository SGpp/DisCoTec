#include "redef.h"
!>Provides a circular concentric geometry
MODULE circular_mod
  use discretization
  use coordinates
  use par_geom, only: geomtype, Bprof_coeffs

  implicit none
  public:: init_circ, get_circ
  private

  INTEGER:: n_pol
  REAL:: r0                     ! r0/Lref where r0 radius at x=0
  REAL:: maj_R                  ! R0/Lref
  REAL:: minor_r                ! r0/Lref
  INTEGER:: NPSI,NCHI, NPSI0, NCHI0
  REAL, DIMENSION(:),   ALLOCATABLE   :: q                         ! safety factor
  REAL, DIMENSION(:),   ALLOCATABLE   :: dqdxa                     ! dqdxa : dq/d(x/a)
  REAL, DIMENSION(:),   ALLOCATABLE   :: r                         ! r/Lref
  REAL, DIMENSION(:),   ALLOCATABLE   :: eps                       ! r/R0
  
! SUBROUTINES DEFINITIONS
CONTAINS
  
  SUBROUTINE init_circ(pn_pol,pq,peps,pdqdxa,pmaj_R,pminor_r)
    INTEGER :: pn_pol 
    REAL    :: pmaj_R, pminor_r
    REAL, dimension(pi1gl:) :: pq,peps,pdqdxa       ! use assumed-shape arrays
    INTEGER :: i

    n_pol=pn_pol
    maj_R=pmaj_R
    minor_r=pminor_r
    allocate(q(pi1gl:pi2gl),r(pi1gl:pi2gl),eps(pi1gl:pi2gl),dqdxa(pi1gl:pi2gl))
    do i=pi1gl,pi2gl
       q(i)=pq(i);  eps(i)=peps(i); r(i)=eps(i)*maj_R; dqdxa(i)=pdqdxa(i);
    end do
  END SUBROUTINE init_circ

  
  SUBROUTINE get_circ(geom,Cy,C_xy,q0)
    ! Using assumed-shape arrays instead of explicit-shape arrays:
    ! When in DEBUG mode, we get an error "subscript out of bounds"
    ! if the actual array is too small.
    TYPE(geomtype), INTENT(INOUT) :: geom
    REAL, DIMENSION(pi1gl:), INTENT(OUT):: Cy, C_xy
    REAL, INTENT(OUT) :: q0
    REAL    :: g11,g12,g22,g33,&! normalized metric in (PSI,CHI,PHI)
         B,dBdr_c,dBdchi        ! J_PCP : J_PSI,CHI,P
    ! intermediate variables, t stands for theta the poloidal angle. 
    REAL    :: qbar,dpsidr,dchidr,dchidt,dBdr,dBdt,cost,sint   

    REAL    :: fac1, fac2, dCydr_Cy
    REAL    :: Bprof_fit, ddr_Bprof_fit
    INTEGER :: i, j, k
    REAL    :: dxdr, dqdr
    REAL    :: dqbardr
    LOGICAL :: use_Bprof_corr = .false.
    REAL, DIMENSION(pi1gl:pi2gl,pj1:pj2,lk1:lk2):: gxx, gxz, gyy, gyz, dBdx

    if (mod((pi1gl+pi2gl),2).eq.0) THEN !odd number of radial points
       r0=r((pi2gl+pi1gl)/2)
       q0=q((pi2gl+pi1gl)/2)
    else !even number of radial points
       !estimate center by linear interpolation of neighboring grid points
       r0=0.5*(r((pi1gl+pi2gl)/2)+r((pi1gl+pi2gl+1)/2))
       q0=0.5*(q((pi1gl+pi2gl)/2)+q((pi1gl+pi2gl+1)/2))
    end if

    do i=pi1gl,pi2gl
       qbar=q(i)*sqrt(1-eps(i)**2)
       dxdr=1.0                ! x=r
       dpsidr=r(i)/qbar                ! 1/(Bref*Lref)*dpsi/dr
       dqdr=dqdxa(i)*dxdr/minor_r
       dqbardr=sqrt(1-eps(i)**2)*(dqdr-q(i)*eps(i)/maj_R/(1-eps(i)**2))
       Cy(i) = r0/q0 !r(i)/q(i)
       dCydr_Cy = 0.0 !1./r(i)+dqdr
       C_xy(i)=dpsidr/(dxdr*Cy(i))

       ! x dep. intermediate variables
       fac1 = (1+eps(i))/(1-eps(i))
       fac2 = (dCydr_Cy+dqdr/q(i))

       ! compute correction function and derivative
       ! if Bprof_coeffs are specified
       if (abs(sum(Bprof_coeffs)) .gt. epsilon(1.0)) use_Bprof_corr = .true.
       if (use_Bprof_corr) then
         Bprof_fit = 1.0
         ddr_Bprof_fit = 0.0
         do k = 1, 8
           Bprof_fit = Bprof_fit + &
             eps(i)**k * (Bprof_coeffs(k) / Bprof_coeffs(0))
           ddr_Bprof_fit = ddr_Bprof_fit + &
             k * eps(i)**(k-1) * (Bprof_coeffs(k) / (Bprof_coeffs(0) * maj_R))
         end do
       end if

       do k=lk1,lk2
          ! (x,z) dep. intermediate variables
          cost = (cos(zval(k)) - eps(i))/(1-eps(i)*cos(zval(k)))
          sint = 2.0*sqrt(fac1)*tan(zval(k)/2.0)/(1+fac1*tan(zval(k)/2.0)**2)

          ! derivatives with respect to r,theta
          dchidr =  -1.0/maj_R*sin(zval(k))/(1-eps(i)**2)
          dchidt = sqrt(1-eps(i)**2)/(1+eps(i)*cost)

          B = 1/(1+eps(i)*cost)*sqrt(1+(eps(i)/qbar)**2)
          if (.not. use_Bprof_corr) then
            dBdr = B*(- cost / (maj_R*(1+eps(i)*cost)) + &
              eps(i)/(maj_R*qbar**2)/(1+(eps(i)/qbar)**2))
          else
            dBdr = B*(- cost / (maj_R*(1+eps(i)*cost)) + &
              (eps(i)/(maj_R*qbar**2)/(1+(eps(i)/qbar)**2)) * &
              (1 - maj_R*eps(i)*dqbardr/qbar) + &
              (dqbardr/qbar) + (ddr_Bprof_fit/Bprof_fit))
          end if
          dBdt = eps(i)*sint / (1+eps(i)*cost)**2 * &
              sqrt(1+(eps(i)/qbar)**2)

          ! metric in (r,CHI,PHI)
          g11=1
          g22=dchidr**2+1/(maj_R*eps(i))**2*dchidt**2
          g12=dchidr
          g33=1/maj_R**2*1/(1+eps(i)*cost)**2

          ! magnetic field derivatives in
          dBdchi=dBdt/dchidt
          dBdr_c=1/g11*(dBdr-g12*dBdt/dchidt)

          Do j=pj1,pj2
             ! metric in (x,y,z)
             gxx(i,j,k)=(dxdr)**2*g11
             gyy(i,j,k)=(Cy(i)*q(i))**2*((fac2*zval(k))**2*g11+&
                  &2*fac2*zval(k)*g12+g22)+Cy(i)**2*g33
             geom%gij(i,j,k)=dxdr*Cy(i)*q(i)*(fac2*zval(k)*g11+g12)
             gxz(i,j,k)=dxdr*g12
             gyz(i,j,k)=Cy(i)*q(i)*(fac2*zval(k)*g12+g22)
             geom%gzz(i,j,k)=g22
             
             ! jacobian
             geom%jacobian(i,j,k)=C_xy(i)*q(i)*maj_R*(1+eps(i)*cost)**2
             
             ! Bfield
             geom%Bfield(i,j,k)=B
             dBdx(i,j,k)=dBdr_c/dxdr
             geom%dBdz(i,j,k)=dBdchi
          end do
          !derivatives with respect to cylindrical coordinates
          geom%dxdR(i,k)=cost
          geom%dxdZ(i,k)=sint
          geom%R_hat(:,k)=maj_R*(1.+eps(i)*cost)
          geom%Z_hat(:,k)=maj_R*(1.+eps(i)*sint)
       end do
    end do

    if (yx_order) then
       geom%gii=gyy
       geom%giz=gyz
       geom%gjj=gxx
       geom%gjz=gxz
       geom%dBdi=0.
       geom%dBdj=dBdx
    else
       geom%gii=gxx
       geom%giz=gxz
       geom%gjj=gyy
       geom%gjz=gyz
       geom%dBdi=dBdx
       geom%dBdj=0.
    end if

    deallocate(q,dqdxa,r,eps)

  END SUBROUTINE get_circ

  !-----------------------------------------------------------
END MODULE circular_mod
