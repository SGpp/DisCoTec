#include "redef.h"
!>Provides Gauss-Legendre weights and knots for mu integration
module GaussQuadrature
  use discretization

  implicit none
  public:: mu_grid_type, GetMuWeightsAndKnots, set_Gauss_defaults

  private
  
  character(len=20) :: mu_grid_type

contains
  
  subroutine set_Gauss_defaults
     mu_grid_type='gau_lag'
  end subroutine set_Gauss_defaults

  subroutine GetMuWeightsAndKnots(muweights,muknots,lw,nw_in)
    integer,intent(in)::nw_in
    real,intent(in):: lw
    real,dimension(1:nw_in),intent(inout):: muknots,muweights

    real:: deltamu, sub_lw_norm
    integer:: i, ndomain, nw_per_block=4
#ifdef COMBI
    integer:: m,filippi_count,cc_count,pEnd,p
    real:: filippi_sum,cc_sum,lastElement,thetaS,stretchFactor,quadratureShift
    REAL, parameter   :: pi= 3.141592653589793239d0
#endif
    real(8),dimension(1:nw_in):: dblweights, dblknots
    real(8) :: dbl_lb=0.0, dbl_ub

    select case(mu_grid_type)
    case('equidist')
       deltamu=lw/nw_in
       do i=1,nw_in
#ifdef COMBI
          muknots(i)=(i-1)*deltamu+WSHIFT
#else
          muknots(i)=(i-0.5)*deltamu
#endif
          muweights(i)=deltamu
       enddo       
#ifdef COMBI
    case('filippi')
        write (*,*) 'Use filippi integration'
        do i=1,nw_in
            muknots(i) = (-cos(pi*i/(nw_in+1))+1)*lw/2.0
            m = floor((nw_in-1.0)/2.0)
            filippi_sum = 0.0
            do filippi_count=0,m
                filippi_sum = filippi_sum + sin((2*filippi_count+1)*i*pi/(nw_in+1))/(2.0*filippi_count+1.0)
            enddo
            muweights(i) = 4.0/(nw_in+1) * sin(i*pi/(nw_in+1)) * filippi_sum

            if (mype.eq.0) then
                write (*,*) i,muknots(i),muweights(i)
            endif

        enddo
    case('clenshaw_curtis')
        if (mype .eq. 0) then
            write (*,*) 'using Clenshaw-Curtis for integration'
        endif

        stretchFactor = 2.7
        quadratureShift = 0.00

        do i = 1, nw_in
            if (mod(nw_in,2).eq.0) then
                lastElement = 0.5
            else
                lastElement = 1.0
            endif
            thetaS = pi*(i-1)/nw_in
            muknots(i) = log(1.0/(cos(thetaS/2.0)**2))

            muweights(i)= 0.5
            pEnd = floor(nw_in/2.0)
            do p = 1, pEnd-1
                muweights(i) = muweights(i) + cos(2.0*p*thetaS)/(1.0-4*p**2)
!                if (mype .eq. 0) then
!                    write (*,*) i,p,cos(2.0*p*thetaS)/(1.0-4*p**2)
!                end if
            end do
            muweights(i) = muweights(i) + lastElement*cos(2.0*pEnd*thetaS)/(1.0-4*p**2)
            muweights(i) = muweights(i) * 2.0 / nw_in * exp(muknots(i)) * stretchFactor

            muknots(i) = muknots(i)*stretchFactor + quadratureShift

            if (mype .eq. 0) then
                write (*,*) i,muknots(i),muweights(i)
            end if

        end do
        muweights(1)=muweights(1)*0.5
#endif
    case('eq_vperp')
       deltamu=lw/nw_in**2
       do i=1,nw_in
          muknots(i)=(i-0.5)**2 *deltamu
          muweights(i)= (2*i-1)*deltamu
       enddo 
    case('gau_leg')
       dbl_ub=lw
       call gauleg(dbl_lb, dbl_ub, dblknots, dblweights, nw_in)
       muweights=real(dblweights)
       muknots=real(dblknots)

       ! rest of the integration is added to the last weight
       muweights(nw_in) = muweights(nw_in)+Exp(muknots(nw_in)-lw)

    case('gau_leg_block')
       !lw subblocks with equal widths containing nw_per_block
       !grid points
       ndomain = nw_in/nw_per_block
       if (modulo(nw_in,nw_per_block).ne.0) &
            &stop 'nw_in not divisible by nw_per_block'
       dbl_ub = 0.0
       do i=0,ndomain-1
          dbl_lb=dbl_ub
          dbl_ub=lw/ndomain*(i+1)
          call gauleg(dbl_lb, dbl_ub, dblknots, dblweights, &
               &nw_per_block)
          muweights(i*nw_per_block+1:(i+1)*nw_per_block)=real(dblweights)
          muknots(i*nw_per_block+1:(i+1)*nw_per_block)=real(dblknots)
       enddo
    case('gau_leg_block_lin')
       !lw subblocks with lin. increasing widths containing 
       !nw_per_block grid points
       ndomain = nw_in/nw_per_block
       if (modulo(nw_in,nw_per_block).ne.0) &
            &stop 'nw_in not divisible by nw_per_block'
       sub_lw_norm = 0.0
       do i=0,ndomain-1
          sub_lw_norm = sub_lw_norm+real(i+1)
          !for exp. increase of box size (does not improve conv.):
          !sub_lw_norm = sub_lw_norm+exp(real(i))
       enddo
       dbl_ub = 0.0
       do i=0,ndomain-1
          dbl_lb=dbl_ub
          dbl_ub=dbl_ub+lw*real(i+1)/sub_lw_norm
          !exp. width: dbl_ub=dbl_ub+lw*exp(real(i))/sub_lw_norm
          call gauleg(dbl_lb, dbl_ub, dblknots, dblweights, &
               &nw_per_block)
          muweights(i*nw_per_block+1:(i+1)*nw_per_block)=real(dblweights)
          muknots(i*nw_per_block+1:(i+1)*nw_per_block)=real(dblknots)
       enddo     

        ! rest of the integration is added to the last weight
       muweights(nw_in) = muweights(nw_in)+Exp(muknots(nw_in)-lw)    
    case('gau_lag')
       if (nw_in.gt.124) stop 'Gauss-Laguerre integration is currently not possible for nw0>124.'
       dbl_ub=lw
       call gaulag(dbl_lb, dbl_ub, dblknots, dblweights, nw_in)
       muweights=real(dblweights)
       muknots=real(dblknots)

    case default
       print*, 'no valid mu_grid_type selected'
       stop
    end select

  end subroutine GetMuWeightsAndKnots

  !>Compute Gauss-Legendre abscissas and weights in interval [x1, x2]
  !!Originally written by T-M Tran (CRPP, Lausanne) but slightly 
  !!modified
  !!\param n size of weights/knots array
  !!\param x1 lower bound
  !!\param x2 upper bound
  !!\param x  Gauss-Legendre knots
  !!\param w  Gauss-Legendre weights
  SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER, INTENT(in) :: n
    REAL(kind=8),INTENT(in) :: x1,x2
    REAL(kind=8),INTENT(out) :: x(n),w(n)
    REAL(kind=8) :: EPS
    INTEGER i,j,m
    REAL(kind=8) :: p1,p2,p3,pp,xl,xm,z,z1

    eps=EPSILON(eps)
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z=COS(3.14159265358979323846d0*(i-.25d0)/(n+.5d0))
       DO
          p1=1.d0; p2=0.d0
          DO j=1,n
             p3=p2; p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
          END DO
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          IF( ABS(z-z1) .LT. EPS ) EXIT
       END DO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    END DO
  END SUBROUTINE gauleg


  !> This routine is adapted from a F77 specfun routine.
  !!\param n size of weights/knots array
  !!\param x1 lower bound
  !!\param x2 upper bound
  !!\param x  Gauss-Legendre knots
  !!\param w  Gauss-Legendre weights
  subroutine gaulag(x1,x2,X,W,N)
!
!       =========================================================
!       Purpose : Compute the zeros of Laguerre polynomial Ln(x)
!                 in the interval [0,oo], and the corresponding
!                 weighting coefficients for Gauss-Laguerre
!                 integration
!       Input :   n    --- Order of the Laguerre polynomial
!                 X(n) --- Zeros of the Laguerre polynomial
!                 W(n) --- Corresponding weighting coefficients
!       =========================================================
!
    implicit none
    integer :: N
    REAL(kind=8),INTENT(in) :: x1,x2
    REAL(kind=8),INTENT(out) :: x(N),w(N)
    integer :: nr, i, k, j, it
    REAL(kind=8) :: hn, pf, pd, z, z0, p, f0, &
         & f1, fd, q, wp, gd, fac
    
    HN=1.0D0/N
    PF=0.0D0
    PD=0.0D0
    do  NR=1,N
       Z=HN
       if (NR.gt.1) Z=X(NR-1)+HN*NR**1.27D0
       IT=0
10     IT=IT+1
       Z0=Z
       P=1.0D0
       do I=1,NR-1
          P=P*(Z-X(I))
       enddo
       F0=1.0D0
       F1=1.0D0-Z
       do K=2,N
          PF=((2.0D0*K-1.0D0-Z)*F1-(K-1.0D0)*F0)/K
          PD=K/Z*(PF-F1)
          F0=F1
          F1=PF
       enddo
       FD=PF/P
       Q=0.0D0
       do I=1,NR-1
          WP=1.0D0
          do J=1,NR-1
             if (J.ne.I) WP=WP*(Z-X(J))
          enddo
          Q=Q+WP
       enddo
       GD=(PD-Q*FD)/P
       Z=Z-FD/GD
       if (IT.le.40.and.ABS((Z-Z0)/Z).gt.1.0D-15) GO TO 10
       X(NR)=Z
       W(NR)=1.0D0/(Z*PD*PD)
    enddo
    w=w*exp(x)
    fac=x2/sum(w)
    x=x*fac
    w=w*fac


  end subroutine gaulag

END MODULE GaussQuadrature

