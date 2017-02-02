#include "redef.h"
#include "intrinsic_sizes.h"
!>Contains the coefficients and initialization routines related to the Runge-Kutta time stepping schemes
!!and a routine to compute numerical stability criteria related to these schemes.
module RK_coefficients
  use par_mod, only: timescheme, direct_rhs, print_ini_msg, coll_split_scheme, coll_split
  use discretization, only: mype
  use polyroots

  implicit none
  public:: initialize_RK_coefficients, finalize_RK_coefficients, compute_dt_from_ev
  public:: explicit_RK, a_rk, b_rk, c_rk, rkstages, low_mem, rk_corr, implicit_scheme
  public:: a_rk_coll, b_rk_coll, c_rk_coll, rkstages_coll,&
       switch_coll_scheme, max_rkstages_coll
  private

  integer :: init_status_std = 0
  integer :: init_status_coll = 0

  !Runge-Kutta-coefficients
  logical:: explicit_RK=.true., implicit_scheme=.false.
  real,dimension(:),allocatable:: a_rk, b_rk, c_rk, RK_Taylor
  real,dimension(:),allocatable:: a_rk_coll, b_rk_coll, c_rk_coll, RK_Taylor_coll
  integer:: rkstages, rkstages_coll=1, max_rkstages_coll=4
  real,dimension(2):: rk_corr
  !real::rk_corr_coll
  logical:: low_mem
contains

  !>Initializes the Runge-Kutta coefficients for the various schemes. The coefficients of the corresponding 
  !!expansion in dt (RK_Taylor) can in principle be computed from the Butcher coefficients, but this is not 
  !!implemented, so that they have to be given explicitly.
  subroutine initialize_RK_coefficients
    implicit none
    integer::rkstages_coll_loc

    !if rkstages_coll should be ill defined, use 1 for initialization
    rkstages_coll_loc = rkstages_coll
    if (rkstages_coll.le.0) rkstages_coll_loc = 1

    call initialize_RK_coefficients_std
    if (init_status_coll.ne.rkstages_coll_loc) call finalize_RK_coefficients_coll
    call initialize_RK_coefficients_coll(rkstages_coll_loc)

  end subroutine initialize_RK_coefficients

  !>initializes RK coefficients for the timescheme
  subroutine initialize_RK_coefficients_std
    
    if (init_status_std.eq.0) then
       select case (timescheme)
       case('RHS')
          allocate(b_rk(rkstages))
          b_rk=(/1.0/)
          implicit_scheme=.true.
          rkstages=1
          low_mem=.false.
          explicit_RK=.false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "only RHS computation"
       case('IE1s','IE1p','IE1f')
          implicit_scheme=.true.
          if(direct_rhs) then
             !if direct_rhs is set, this is only the performance optimization -> switch to RK4
             rkstages=4
             allocate(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
             a_rk=(/0.5,0.5,1.0/)
             b_rk=(/1.0,2.0,2.0,1.0/)/6.
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "4th order Runge-Kutta for performance optimization"
          else
             rkstages=1
             explicit_RK=.false.
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "1st order implicit Euler"
          endif
          low_mem=.false.
       case('RK3','RK3f')
          rkstages=3
          allocate(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          a_rk=(/1.0,2.0/)/3.
          b_rk=(/0.25,0.0,0.75/)
          RK_Taylor=(/1., 1., 1./2, 1./6/)
          low_mem=.false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "3rd order Runge-Kutta"
       case('RK4','RK4f')
          rkstages=4
          allocate(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          a_rk=(/0.5, 0.5, 1.0/)
          b_rk=(/1.0, 2.0, 2.0, 1.0/)/6.
          RK_Taylor=(/1., 1., 1./2, 1./6, 1./24/)
          low_mem=.false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "4th order Runge-Kutta"
       case('RK4M','RK4Mf')
          rkstages=6
          allocate(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          a_rk=(/0.16791846623918, 0.48298439719700, 0.70546072965982,0.09295870406537, 0.76210081248836/)
          b_rk=(/-0.15108370762927, 0.75384683913851, -0.36016595357907, 0.52696773139913, 0., 0.23043509067071/)
          RK_Taylor=(/1., 1., 1./2, 1./6, 1./24, 0.005562334, 0.0009340186/)
          low_mem=.false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "modified 4th order Runge-Kutta"
       case('RK4T')
          rkstages=4
          allocate(a_rk(2:rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          a_rk=(/0.311526414979283, 0.422649730810374, 1.2658238542103426/)
          !       a_rk=(/1.2658238542103426,0.422649730810374,0.311526414979283/)
          b_rk=(/0.25, 0.25, 0.25, 0.25/)
          RK_Taylor=(/1., 1., 1./2, 1./6, 1./24/)
          low_mem=.false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "equal weights 4th order Runge-Kutta"
       case('RK3lm')
          rkstages=4
          !compared to standard Butcher notation, the last entry of b has been moved to a 
          !to fit into the structure of the implementation
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          !NASA report Kennedy,Carpenter,Lewis: Low-storage, Explicit Runge-Kutta Schemes 
          !for the Compressible Navier-Stokes Equation (they claim it's 3rd order, 
          !but at least linearly, it seems to be 4th order)
          a_rk=(/11847461282814./36547543011857.,&
               3943225443063./7078155732230.,&
               -346793006927./4029903576067.,& 
               -101169746363290./37734290219643./)
          b_rk=(/1017324711453./9774461848756.,&
               8237718856693./13685301971492.,&
               57731312506979./19404895981398.,&
               0./)
          RK_Taylor=(/1., 1., 1./2, 1./6, 4.166666666666667E-002/)
          low_mem=.true.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "3rd order low memory Runge-Kutta (NASA)"
       case('RK3ssp')
          rkstages=4
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          !Spiteri, Ruuth: A new class of optimal high-order strong-stability-preserving time
          !discretization methods. The last 4 digits seem incorrect, but I still copied them..
          a_rk=(/1.03216665875130, 0.18793881263711, 0.15215751854315, 0.65675174856653/)
          b_rk=(/0.10250480379728, 0.18793881266399, 0.05280467502407, 0./)
          RK_Taylor=(/1., 1., 1./2, 1./6, 1.938478371506826E-002/)
          low_mem=.true.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "3rd order low memory Runge-Kutta (SSP)"
       case('turbo')
          rkstages=4
          allocate(a_rk(rkstages),b_rk(rkstages),c_rk(rkstages),RK_Taylor(rkstages+1))
          a_rk=(/1./3, -3./8, 1./4,  1./4/)
          b_rk=(/1./3, -3./4, 1./2, -1./6/)
          c_rk=(/  1., -3./2,   0.,    0./)
          RK_Taylor=(/1., 1., 1./2, 1./6, 1./32/)
          low_mem=.true.
       case('EE1','RKC1')
          rkstages=1
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          call initialize_adaptive_RKC(rkstages,a_rk,b_rk,RK_Taylor)
          low_mem = .false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "explicit Euler (1st order, 1 stage)"
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "WARNING: unconditionally unstable for advection problems"
       case('RKC2')
          rkstages=2
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          call initialize_adaptive_RKC(rkstages,a_rk,b_rk,RK_Taylor)
          low_mem = .false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages, " stage Runge-Kutta-Chebychev"
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "WARNING: unconditionally unstable for advection problems"
       case('RKC3')
          rkstages=3
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          call initialize_adaptive_RKC(rkstages,a_rk,b_rk,RK_Taylor)
          low_mem = .false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages, " stage Runge-Kutta-Chebychev"
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "WARNING: unconditionally unstable for advection problems"
       case('RKC4')
          rkstages=4
          allocate(a_rk(rkstages),b_rk(rkstages),RK_Taylor(rkstages+1))
          call initialize_adaptive_RKC(rkstages,a_rk,b_rk,RK_Taylor)
          low_mem = .false.
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages, " stage Runge-Kutta-Chebychev"
          if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "WARNING: unconditionally unstable for advection problems"
       case default
          if (mype.eq.0) write(*,"(3A)") "the timescheme '",timescheme,"' is unknown, please choose a valid scheme"
          stop
       end select

       !this computes the stability interval on the imaginary / negative real axis (i.e. Re=0/Im=0) for dt=1.
       if (explicit_RK) then
          call compute_stability_criterion((0.,1.),rk_corr(1),.false.)
          call compute_stability_criterion((-1.,0.),rk_corr(2),.false.)
       endif
    endif 
    init_status_std = 1

  end subroutine initialize_RK_coefficients_std
   
  !>initializes RK coefficients for collisions, 
  !!if operator splitting is applied
  subroutine initialize_RK_coefficients_coll(s_in)
    implicit none
    integer,intent(in)::s_in !input for adaptive scheme, default: 1

    if (coll_split) then
       if (init_status_coll.ne.rkstages_coll) then
          select case (coll_split_scheme)
          case('RK2')
             rkstages_coll=2
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             !a_rk_coll=(/0.0,0.5/)  !generic: /(0,x/)
             !b_rk_coll=(/0.0,1.0/)            /(1-1/(2x),1/(2x)/)
             a_rk_coll=(/0.0,0.8/)
             b_rk_coll=(/0.375,0.625/)
             RK_Taylor_coll=(/1., 1., 1./2/)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions: RK2"
          case('RK3')
             rkstages_coll=3
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             a_rk_coll=(/0.0,1.0,2.0/)/3.
             b_rk_coll=(/0.25,0.0,0.75/)
             RK_Taylor_coll=(/1., 1., 1./2, 1./6/)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions: RK3"
          case('RK4')
             rkstages_coll=4
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             a_rk_coll=(/0.0,0.5, 0.5, 1.0/)
             b_rk_coll=(/1.0, 2.0, 2.0, 1.0/)/6.
             RK_Taylor_coll=(/1., 1., 1./2, 1./6, 1./24/)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions: RK4"
          case('EE1','RKC1')
             rkstages_coll=1
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             call initialize_adaptive_RKC(rkstages_coll,a_rk_coll,b_rk_coll,RK_Taylor_coll)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions:"
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "explicit Euler (= 1st order, 1 stage Runge-Kutta-Chebychev)"
          case('RKC2')
             rkstages_coll=2
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             call initialize_adaptive_RKC(rkstages_coll,a_rk_coll,b_rk_coll,RK_Taylor_coll)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions:"
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages_coll, " stage Runge-Kutta-Chebychev"
          case('RKC3')
             rkstages_coll=3
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             call initialize_adaptive_RKC(rkstages_coll,a_rk_coll,b_rk_coll,RK_Taylor_coll)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions:"
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages_coll, " stage Runge-Kutta-Chebychev"
          case('RKC4')
             rkstages_coll=4
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             call initialize_adaptive_RKC(rkstages_coll,a_rk_coll,b_rk_coll,RK_Taylor_coll)
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A)") "Using operator splitting for collisions:"
             if ((mype.eq.0).and.print_ini_msg) write(*,"(A,I1,A)") "1st order", rkstages_coll, " stage Runge-Kutta-Chebychev"
          case('RKCa')
             rkstages_coll=s_in
             allocate(a_rk_coll(rkstages_coll),b_rk_coll(rkstages_coll),RK_Taylor_coll(rkstages_coll+1))
             call initialize_adaptive_RKC(rkstages_coll,a_rk_coll,b_rk_coll,RK_Taylor_coll)
             max_rkstages_coll = 4  !we have 4, 3, 2 and 1 stages implemented
          case ('none')
          case default
             if (mype.eq.0) write(*,"(3A)") "coll_split_scheme = '",coll_split_scheme,"' is unknown, please choose a valid scheme"
             stop
          end select

          !call compute_stability_criterion((-1.,0.),rk_corr_coll,.true.)
          !Write(*,"(A,G12.4)") 'the maximum dt for a collisional damping rate of 1 is ', rk_corr_coll
          !Write(*,"(A)") '[also the maximum (negative real) stable eigenvalue at dt=1]'
       endif 
       init_status_coll=rkstages_coll
    endif 

  end subroutine initialize_RK_coefficients_coll

  !s stage version of first order Runge Kutta Chebychev Method (Verver Applied Numerical Mathematics 22 1996 359-379)
  !max. dt ~ s**2, but s>2 require extra entries in the a_rk matrix (extra temporary arrays in the RK scheme) 
  subroutine initialize_adaptive_RKC(s,a,b,t)
    !s: number of stages (in)
    !a: a_rk
    !b: b_rk
    !t: coefficients of the taylor expansion of the RK scheme
    implicit none
    integer,intent(in)::s
    real,dimension(s),intent(inout)::a,b
    real,dimension(s+1),intent(inout)::t
    select case (s)
    case (1)
       !explicit Euler method
       !a_rk(1) is not used..
       a=(/0.0/)
       b=(/1.0/)
       t=(/1., 1./)
    case (2)
       !a_i=sum(a_ij): diagonal RK scheme in Butcher notation.
       !the damped chebyshev stability polynomial requires
       !t_1=1, t_2=sum(b_i)=1, t_3=b_3*a_3+b_2*a_2, t_4=b_3*a_3*a_2
       !we make the following assumptions considering internal stability of the temporary g's
       !that are computed at lower order, i.e.
       !a_1 < max_dt(s=3)/max_dt(s=1), a_2 < max_dt(s=3)/max_dt(s=2)

       !first version
       !a=(/0.0,0.256134735558603887639606000448/)
       !b=(/0.5,0.5/)

       !second version: sub-polynomials are stable (with abs and |Pj|==1-0.05)
       a=(/0.0,0.1706953512053041/)
       b=(/0.2497313671696385,0.7502686328303615/)
       t=(/1., 1.,0.128067367779301943819803000224/)
    case (3)
       !first version
       !a=(/0.0,0.1142497465633991,0.4460533098497248/)
       !b=(/0.0002730531845217814,0.886078706826887,0.11364823998859093/)

       !second version: sub-polynomials are stable (with abs and |Pj|==1-0.05)
       b=(/0.1115050079072208,-0.2891490316411831,1.177644023733962/)
       a=(/0.0,0.03569626261410427,0.1377742151226115/)
       t=(/1., 1.,0.1519274412957031,0.005791682236923549/)
    case (4)
       !first version
       !a=(/0.0,0.06436995847949907,0.2513128816328438,0.5628510386702287/)
       !b=(/-0.2663166985235108,0.869253213935325,0.3824639680427712,0.01459951654541513/)

       !second version: sub-polynomials are stable (with abs and |Pj|==1-0.05)
       a=(/0.0,0.01569528165214393,0.04916684563864518,0.1484202501429531/)
       b=(/0.001610859821807129,0.1193627217893438,-0.2816076282040107,1.16063404659286/)
       t=(/1., 1.,0.1602892682704187,0.00825224619254291,0.0001329321183124034/)
    endselect
  end subroutine initialize_adaptive_RKC


  subroutine switch_coll_scheme(s_in)
     implicit none
     integer,intent(in)::s_in

     if (init_status_coll.ne.s_in) call finalize_rk_coefficients_coll
     call initialize_rk_coefficients_coll(s_in)

  end subroutine switch_coll_scheme


  subroutine finalize_RK_coefficients

       call finalize_rk_coefficients_std
       call finalize_rk_coefficients_coll

  end subroutine finalize_RK_coefficients
  
  subroutine finalize_RK_coefficients_std

    if (init_status_std.gt.0) then
       if(explicit_RK) deallocate(a_rk,b_rk,RK_Taylor)
       if(timescheme.eq.'RHS') deallocate(b_rk)
       if(timescheme.eq.'turbo') deallocate(c_rk)
       rkstages=0
    endif
    init_status_std=0

  end subroutine finalize_RK_coefficients_std

  subroutine finalize_RK_coefficients_coll

    if (init_status_coll.gt.0) then
       deallocate(a_rk_coll,b_rk_coll,RK_Taylor_coll)
    endif

    init_status_coll=0

  end subroutine finalize_RK_coefficients_coll

  !>This routine is only a wrapper to avoid circular calls of the initialize_RK_coefficients routine 
  !!and compute_stability_criterion
  subroutine compute_dt_from_ev(ev,dtmax,for_coll_scheme)
    complex,intent(in):: ev
    real,intent(out):: dtmax
    logical,intent(in)::for_coll_scheme

    call initialize_RK_coefficients
    call compute_stability_criterion(ev,dtmax,for_coll_scheme)

  end subroutine compute_dt_from_ev

  !>Computes the criterion for neutral stability
  !!
  !!The linearized GK equation for an eigenvector \f$|\lambda\rangle\f$ can be written as 
  !!\f$\partial_t |\lambda\rangle=L|\lambda\rangle=\lambda |\lambda\rangle\f$ (L is the linear operator, 
  !!\f$\lambda\f$) the eigenvalue. The time evolution is approximated by a Runge-Kutta scheme, 
  !!\f$|\lambda(t_0+\Delta t)\rangle=\sum_n c_n \Delta t^n (\partial_t^n |\lambda\rangle)|_{t_0}=
  !!\sum_n c_n \Delta t^n \lambda^n |\lambda(t_0)\rangle\f$, where the \f$c_n\f$ are the coefficients of the 
  !!dt expansion corresponding to the RK scheme.
  !!The system is neutrally stable if the amplitude of \f$|\lambda\rangle\f$ remains constant, i.e. if
  !!\f$|\sum_n c_n \Delta t^n \lambda^n|^2=1\f$. This is a polynomial equation in \f$(\Delta t \lambda)\f$. 
  !!The 1 on the right cancels with the n=0 contribution on the left, the remaining equation can be devided by
  !!\f$(\Delta t \lambda)\f$ because the trivial solution 0 is not of interest. The remainder is a polynomial of 
  !!order (2*rkstages-1) with three (real) degrees of freedom: \f$\Delta t, Re(\lambda), Im(\lambda)\f$. If 
  !!\f$\lambda\f$ is specified, the (maximal) \f$\Delta t\f$ can be computed, if \f$\Delta t=1,Re(\lambda)=0\f$ are
  !!specified, the root of the polynomial corresponds to the stability interval on the imaginary axis
  !!(i.e. for advection problems) which is used for (approximate) dt estimates.
  subroutine compute_stability_criterion(cnum,res,for_coll_scheme)
    complex,intent(in):: cnum
    real,intent(out):: res
    logical,intent(in)::for_coll_scheme
    Character(len=6):: scheme
    integer::stages
    real:: re,im
    integer:: nroots, rnum
    real, dimension(:), allocatable:: poly_ev
    complex(8), dimension(:), allocatable:: droots
    complex, dimension(:), allocatable:: roots
    real,dimension(0:11):: arr
    real,dimension(0:6)::b    

    b=0.
    if (for_coll_scheme) then
       scheme=coll_split_scheme
       stages=rkstages_coll
       b(0:stages)= RK_Taylor_coll
    else
       scheme=timescheme
       stages=rkstages
       b(0:stages)= RK_Taylor
    endif

    select case (scheme)    
    case('IE1s','IE1p','IE1f','RHS')
       res=1.
       return
    end select

    nroots=2*stages-1
    allocate(poly_ev(2*stages), droots(2*stages-1), roots(2*stages-1))
    res=-1.
    re=real(cnum)
    im=Aimag(cnum)
    !this is stability polynomial for up to six stages
    arr=(/2*b(1)*re, &
         2*b(2)*(-im**2+re**2)+b(1)**2*(im**2+re**2),& 
         2*re*(b(3)*(-3*im**2+re**2)+b(1)*b(2)*(im**2+re**2)),& 
         b(2)**2*(im**2+re**2)**2+2*(b(1)*b(3)*(-im**4+re**4)+b(4)*(im**4-6*im**2*re**2+re**4)),& 
         2*re*(b(2)*b(3)*(im**2+re**2)**2+b(5)*(5*im**4-10*im**2*re**2+re**4)+& 
         b(1)*b(4)*(-3*im**4-2*im**2*re**2+re**4)), &
         b(3)**2*(im**2+re**2)**3+2*(-(b(2)*b(4)*(im**2-re**2)*(im**2+re**2)**2)+&
         b(6)*(-im**6+15*im**4*re**2-15*im**2*re**4+re**6)+b(1)*b(5)*(im**6-5*im**4*re**2-5*im**2*re**4+re**6)),&
         2*re*(im**2+re**2)*(b(3)*b(4)*(im**2+re**2)**2+b(1)*b(6)*(5*im**4-10*im**2*re**2+re**4)+&
         b(2)*b(5)*(-3*im**4-2*im**2*re**2+re**4)), &
         (im**2+re**2)**2*(b(4)**2*(im**2+re**2)**2+2*(b(3)*b(5)*(-im**4+re**4)+&
         b(2)*b(6)*(im**4-6*im**2*re**2+re**4))),&
         2*re*(im**2+re**2)**3*(b(3)*b(6)*(-3*im**2+re**2)+b(4)*b(5)*(im**2+re**2)), &
         (im**2+re**2)**4*(2*b(4)*b(6)*(-im**2+re**2)+b(5)**2*(im**2+re**2)),& 
         2*b(5)*b(6)*re*(im**2+re**2)**5, &
         b(6)**2*(im**2+re**2)**6/)
    poly_ev=arr(0:nroots)

    !find the roots of the polynomial
    call zroots(poly_ev,nroots,roots,.true.)

    !return only the maximal real+positive solutions. For cnum=\lambda, this corresponds to the maximal dt 
    !allowed for the input eigenvector; for cnum=(0.,1.)/(-1.,0.) to the extent of the stable region along 
    !the imaginary /negative realaxis for dt=1.
    do rnum=1,nroots
       if (aimag(roots(rnum)).le.1e-8.and.real(roots(rnum)).gt.0.0) res=max(res,real(roots(rnum)))
    enddo
    deallocate(poly_ev, roots, droots)

  end subroutine compute_stability_criterion

  
end module RK_coefficients
