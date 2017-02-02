!> Fifth order Runge-Kutta scheme used by the tracer module. 
module tracer_rk5
  use tracer_rk5_util
  use tracer_IO
  implicit none
  
  public:: odeint, rhs
  private

contains

  !> Integrate ordinary differential equations.
  subroutine odeint(ystart,x1,x2,eps,h1,hmin)
       
    real, dimension(:), intent(INOUT) :: ystart
    real, intent(IN) :: x1,x2,eps,h1,hmin
    real, parameter :: TINY=1.0e-30
    integer, parameter :: MAXSTP=1000000
    integer :: nstp
    real :: h,hdid,hnext,x
    real, dimension(size(ystart)) :: dydx,y,yscal
    
    x=x1
    h=sign(h1,x2-x1)
    y(:)=ystart(:)
    
    do nstp=1,MAXSTP
       call rhs(y,dydx)
       yscal(:)=abs(y(:))+abs(h*dydx(:))+TINY
       
       if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x
       call rkqs(y,dydx,x,h,eps,yscal,hdid,hnext)
       
       if ((x-x2)*(x2-x1) >= 0.0) then
          ystart(:)=y(:)
          return
       end if
       if (abs(hnext) < hmin)&
            call nrerror('stepsize smaller than minimum in odeint')
       h=hnext
    end do
    call nrerror('too many steps in odeint')
    
  end subroutine odeint

  !================================================================================!

  subroutine rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext)

    real, dimension(:), intent(INOUT) :: y
    real, dimension(:), intent(IN) :: dydx,yscal
    real, intent(INOUT) :: x
    real, intent(IN) :: htry,eps
    real, intent(OUT) :: hdid,hnext
        integer :: ndum
    real :: errmax,h,htemp,xnew
    real, dimension(size(y)) :: yerr,ytemp
    real, parameter :: SAFETY=0.9,PGROW=-0.2,PSHRNK=-0.25,&
         ERRCON=1.89e-4
    
    ndum=assert_eq(size(y),size(dydx),size(yscal),'rkqs')
    h=htry
    do
       call rkck(y,dydx,h,ytemp,yerr)
       errmax=0.
       errmax=maxval(abs(yerr(:)/yscal(:)))/eps
       if (errmax <= 1.0) exit
       htemp=SAFETY*h*(errmax**PSHRNK)
       h=sign(max(abs(htemp),0.1*abs(h)),h)
       xnew=x+h
       if (xnew == x) call nrerror('stepsize underflow in rkqs')
    end do
    if (errmax > ERRCON) then
       hnext=SAFETY*h*(errmax**PGROW)
    else
       hnext=5.0*h
    end if
    hdid=h
    x=x+h
    y(:)=ytemp(:)
    
  end subroutine rkqs
  
  !================================================================================!
  
  subroutine rkck(y,dydx,h,yout,yerr)

    real, dimension(:), intent(IN) :: y,dydx
    real, intent(IN) :: h
    real, dimension(:), intent(OUT) :: yout,yerr
    integer :: ndum
    real, dimension(size(y)) :: ak2,ak3,ak4,ak5,ak6,ytemp
    real, parameter :: A2=0.2,A3=0.3,A4=0.6,A5=1.0,&
         A6=0.875,B21=0.2,B31=3.0/40.0,B32=9.0/40.0,&
         B41=0.3,B42=-0.9,B43=1.2,B51=-11.0/54.0,&
         B52=2.5,B53=-70.0/27.0,B54=35.0/27.0,&
         B61=1631.0/55296.0,B62=175.0/512.0,&
         B63=575.0/13824.0,B64=44275.0/110592.0,&
         B65=253.0/4096.0,C1=37.0/378.0,&
         C3=250.0/621.0,C4=125.0/594.0,&
         C6=512.0/1771.0,DC1=C1-2825.0/27648.0,&
         DC3=C3-18575.0/48384.0,DC4=C4-13525.0/55296.0,&
         DC5=-277.0/14336.0,DC6=C6-0.25
    ndum=assert_eq(size(y),size(dydx),size(yout),size(yerr),'rkck')
    ytemp=y+B21*h*dydx
    call rhs(ytemp,ak2)
    ytemp=y+h*(B31*dydx+B32*ak2)
    call rhs(ytemp,ak3)
    ytemp=y+h*(B41*dydx+B42*ak2+B43*ak3)
    call rhs(ytemp,ak4)
    ytemp=y+h*(B51*dydx+B52*ak2+B53*ak3+B54*ak4)
    call rhs(ytemp,ak5)
    ytemp=y+h*(B61*dydx+B62*ak2+B63*ak3+B64*ak4+B65*ak5)
    call rhs(ytemp,ak6)
    yout=y+h*(C1*dydx+C3*ak3+C4*ak4+C6*ak6)
    yerr=h*(DC1*dydx+DC3*ak3+DC4*ak4+DC5*ak5+DC6*ak6)
  end subroutine rkck
  


  subroutine rhs(y,yp)
    implicit none
    
    !Input arguments:
    !    
    !phi     
    !y(1)=r
    !y(2)=z
    !y(3)=c(1,1)
    !y(4)=c(1,2)
    !y(5)=c(1,3)
    !y(6)=c(2,1)
    !y(7)=c(2,2)
    !y(8)=c(3,3)
    
    real, dimension(:), intent(OUT) :: yp
    real, dimension(:), intent(IN)  :: y
    
    real :: Bpm1,df1_dy1,df2_dy1,df1_dy2,df2_dy2
    real :: df3_dy1,df3_dy2,df1_dy3,df2_dy3,df3_dy3
    real :: rrr,zzz
    real:: Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,&
         &dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
    !Coords 
    
    rrr = y(1)
    zzz = y(2)
    
    call get_from_efit(rrr,zzz,Br,Bp,Bz,dBrdR,dBrdp,dBrdZ,&
         &dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
    
    Bpm1 = 1.d0/Bp
    
    yp(1) = Br*rrr*Bpm1
    yp(2) = Bz*rrr*Bpm1
    
    !Auxilliaries
    
    df1_dy1 = dBrdR*rrr * Bpm1
    df2_dy1 = dBzdR*rrr * Bpm1
    df1_dy2 = rrr*dBrdZ * Bpm1
    df2_dy2 = rrr*dBzdZ * Bpm1
    df3_dy1 = Bpm1*dBpdR - 1.d0/rrr
    df3_dy2 = Bpm1*dBpdZ
    df1_dy3 = rrr*Bpm1*dBrdp
    df2_dy3 = rrr*Bpm1*dBzdp
    df3_dy3 = dBpdp *Bpm1
    
    !RHS
    
    yp(3) = -(  y(3)*df1_dy1 + y(4)*df2_dy1 + y(5)*df3_dy1 )
    yp(4) = -(  y(3)*df1_dy2 + y(4)*df2_dy2 + y(5)*df3_dy2 )
    yp(5) = -(  y(3)*df1_dy3 + y(4)*df2_dy3 + y(5)*df3_dy3 )
    yp(6) = -(  y(6)*df1_dy1 + y(7)*df2_dy1 + y(8)*df3_dy1 )
    yp(7) = -(  y(6)*df1_dy2 + y(7)*df2_dy2 + y(8)*df3_dy2 )
    yp(8) = -(  y(6)*df1_dy3 + y(7)*df2_dy3 + y(8)*df3_dy3 )
    
    
  end subroutine rhs
  

end module tracer_rk5


