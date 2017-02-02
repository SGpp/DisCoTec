!>Finds the complex roots of a polynomial (only used if the corresponding 
!!NAG routines are not available)
MODULE polyroots
  implicit none
  public:: zroots
  private

  INTERFACE zroots
     MODULE PROCEDURE zroots_complex, zroots_real
  END INTERFACE

CONTAINS
  !> zroots finds the complex roots of a polynomial

  !>zroots for real input
  SUBROUTINE zroots_real(a,m,roots,polish)
    INTEGER :: m  !<degree of the polynomial
    !>the m+1 complex coefficients of the polynomial in ascending
    !! order. a1+a2*z+a3*z^2+...+a(m+1)*z^(m)
    REAL    :: a(m+1)
    COMPLEX :: roots(m) !<result
    LOGICAL :: polish !<enables polishing (whatsoever that is...)
    COMPLEX :: a_complex(m+1)
    
    a_complex = a
    CALL zroots_complex(a_complex,m,roots,polish)

  END SUBROUTINE zroots_real

  !>zroots for complex input
  SUBROUTINE zroots_complex(a,m,roots,polish)
    INTEGER :: m  !<degree of the polynomial
    !>the m+1 complex coefficients of the polynomial in ascending
    !! order. a1+a2*z+a3*z^2+...+a(m+1)*z^(m)
    COMPLEX :: a(m+1)
    COMPLEX :: roots(m) !<result
    LOGICAL :: polish !<enables polishing (whatsoever that is...)
    REAL,PARAMETER :: EPS=1.0e-6
    INTEGER,parameter :: MAXM=101

    INTEGER :: i,j,jj,its
    COMPLEX :: ad(MAXM),x,b,c

    DO j=1,m+1
       ad(j)=a(j)
    END DO

    DO j=m,1,-1
       x=CMPLX(0.,0.)
       CALL laguer(ad,j,x,its)
       IF(ABS(AIMAG(x)).LE.2.*EPS**2*ABS(REAL(x))) x=CMPLX(REAL(x),0.)
       roots(j)=x
       b=ad(j+1)
       DO jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
       END DO
    END DO
    
    IF (polish) THEN
       DO j=1,m
          CALL laguer(a,m,roots(j),its)
       END DO
    ENDIF
    
    DO j=2,m
       x=roots(j)
       DO i=j-1,1,-1
          IF(REAL(roots(i)).LE.REAL(x)) GOTO 10
          roots(i+1)=roots(i)
       END DO
       i=0
10     roots(i+1)=x
     END DO
   END SUBROUTINE zroots_complex

   SUBROUTINE laguer(a,m,x,its)
     INTEGER :: m,its
     COMPLEX :: a(m+1),x

     REAL,PARAMETER :: EPSS=2.0e-7
     INTEGER,PARAMETER :: MR=8,MT=10,MAXIT=MT*MR

     INTEGER :: iter,j
     REAL :: abx,abp,abm,err,frac(MR)=(/.5,.25,.75,.13,.38,.62,.88,1./)
     COMPLEX :: dx,x1,b,d,f,g,h,sq,gp,gm,g2
     SAVE frac

     DO iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=ABS(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        DO j=m,1,-1
           f=x*f+d
           d=x*d+b
           b=x*b+a(j)
           err=ABS(b)+abx*err
        END DO
        err=EPSS*err
        IF(ABS(b).LE.err) THEN
           RETURN
        ELSE
           g=d/b
           g2=g*g
           h=g2-2.*f/b
           sq=SQRT((m-1)*(m*h-g2))
           gp=g+sq
           gm=g-sq
           abp=ABS(gp)
           abm=ABS(gm)
           IF(abp.LT.abm) gp=gm
           IF (MAX(abp,abm).GT.0.) THEN
              dx=m/gp
           ELSE
              dx=EXP(CMPLX(LOG(1.+abx),real(iter)))
           ENDIF
        ENDIF
        x1=x-dx
        IF(x.EQ.x1) RETURN
        IF (MOD(iter,MT).NE.0) THEN
           x=x1
        ELSE
           x=x-dx*frac(iter/MT)
        ENDIF
     END DO
     STOP 'too many iterations in polyroots/laguer'
   END SUBROUTINE laguer

 END MODULE polyroots
