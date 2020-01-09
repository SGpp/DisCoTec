#include "redef.h"
#include "intrinsic_sizes.h"
!>Module for antenna drive of A_par in astrophysical applications. 
Module antenna
  Use par_mod
  Use communications
  USE coordinates
  use mtrandom
  use blockindex
  use vel_space
  use gyro_average_ff_mod, only: jfac

  implicit none

  public:: Apar0_antenna, omega0_antenna, lv_antenna_amp, lv_antenna_modes, lv_antenna_freq
  public:: Apar_pre_antenna, mem_est_antenna, check_antenna, add_apar_antenna, initialize_antenna
  public:: finalize_antenna, set_antenna_defaults, antenna_type, evolve_antenna_amplitudes, n_antennae
  public:: lv_antenna_initamp, amp, antenna_contrib, add_dApar_dt_antenna, calc_antenna_evolution
  public:: get_dApar_dt_antenna, get_Apar_antenna
  
  private
  
  integer:: init_status=0
  integer:: antenna_type !0=Off; 1=single antenna applied to all modes; 2=Langevin antenna; 3=mod. Langevin antenna 
  !Ref. for Langevin antenna: J.M. TenBarge et al., CPC 185, 578 (2014)

  !array for storing Apar without the antenna contribution (this will be written to field.dat!)
  COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: Apar_pre_antenna

  !parameters for type 1 antenna
  real:: Apar0_antenna, omega0_antenna
  logical:: antenna_contrib

  !parameters for Langevin antenna
  integer, parameter:: nmax=10 !maximum number of antennae
  integer:: n_antennae
  integer, dimension(nmax*3):: lv_antenna_modes
  complex, dimension(nmax):: lv_antenna_amp, lv_antenna_initamp, lv_antenna_freq
  real, dimension(nmax):: k_par
  integer, dimension(nmax):: i_ant, j_ant
  real, dimension(:,:,:,:), allocatable:: prefac

  !arrays for random number generator (mtrandom.F90)
  integer, dimension(:,:), allocatable, target:: mt
  integer, dimension(:), allocatable,target:: mti

  complex, dimension(:), allocatable:: amp, d_amp_dt
contains


  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_antenna(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

 
    !Apar_pre_antenna
    IF (antenna_type.ne.0) mem_loc = mem_loc + SIZE_OF_COMPLEX_MB * li0*lj0*lz0

    mem_est_antenna=mem_req_in+mem_loc
  End Function mem_est_antenna


  subroutine set_antenna_defaults
    antenna_type=0
    Apar0_antenna=0.
    omega0_antenna=0.
    antenna_contrib=.false.

    lv_antenna_modes=0
    lv_antenna_amp=0.
    lv_antenna_initamp=0.
    lv_antenna_freq=0.
  end subroutine set_antenna_defaults

  subroutine check_antenna
    if (Apar0_antenna.ne.0.0) antenna_type=1
    if (antenna_type==1.and.Apar0_antenna.eq.0.0) antenna_type=0
    if (antenna_type==2.and.all(lv_antenna_amp.eq.0.0).and.all(lv_antenna_initamp.eq.0.0)) antenna_type=0
    if (antenna_type==3.and.all(lv_antenna_amp.eq.0.0).and.all(lv_antenna_initamp.eq.0.0)) antenna_type=0
!    IF ((Apar0_antenna .NE. 0.0) .AND. (ABS(omega0_antenna) .GT. 25.13274 / dt_max)) &
!      STOP "set lower antenna frequency or lower dt_max"
    if (antenna_type==3) nonlin_h=.true.

    IF ((antenna_type.ne.0) .AND. (beta.EQ.0)) &
         STOP "cannot use apar antenna in electrostatic runs"
    if (.not.xy_local) stop "antenna implemented only for local runs"
    if (antenna_type.ne.0) antenna_contrib=.true.
  end subroutine check_antenna

  subroutine initialize_antenna

    select case(antenna_type)
    case(1)
       if (init_status.eq.0) &
            ALLOCATE(Apar_pre_antenna(li1:li2,lj1:lj2,lbz:ubz))
    case(2)
       call initialize_langevin_antenna
    case(3) 
       call initialize_prefac
       call initialize_langevin_antenna
    end select

    init_status=1

  end subroutine initialize_antenna

  subroutine initialize_langevin_antenna
    integer:: p

    if (init_status.eq.0) then
       ALLOCATE(Apar_pre_antenna(li1:li2,lj1:lj2,lbz:ubz))
       !Apar_pre_antenna is set to zero to avoid uninitialized boundary points
       Apar_pre_antenna=0.
       n_antennae=0
       do p=1,nmax
          !if no specific initial amplitude is specified, we use lv_antenna_amp
          if (abs(lv_antenna_initamp(p)).eq.0.) lv_antenna_initamp(p)=lv_antenna_amp(p)
          
          !count antennae
          if (abs(lv_antenna_amp(p)).ne.0..and..not.all(lv_antenna_modes(p*3-2:p*3).eq.0)) &
               n_antennae=n_antennae+1
          
          !initialize antenna modes
          i_ant(p)=lv_antenna_modes(p*3-2)
          do while (i_ant(p)<0) 
             i_ant(p)=i_ant(p)+nx0
          end do
          j_ant(p)=lv_antenna_modes(p*3-1)
          k_par(p)=real(lv_antenna_modes(p*3))
       enddo
       allocate(mt(1:n_antennae,nrnd),mti(1:n_antennae),amp(1:n_antennae),d_amp_dt(1:n_antennae))
    endif
    do p=1,n_antennae
       !initialize separate random number generator for each antenna
       call seed_rnd(p,lv_antenna_amp(p))
       !initialize all antennae with complex amplitudes (possibly from an earlier run)
       amp(p)=lv_antenna_amp(p)
    enddo
    call calc_antenna_evolution

  end subroutine initialize_langevin_antenna

  subroutine finalize_antenna
    select case(antenna_type)
    case(1)
       deallocate(Apar_pre_antenna)
    case(2)
       !we write out the amplitude status at the end of the simulation to allow smooth continuation
       if (.not.in_perf_opt) lv_antenna_amp(1:n_antennae)=amp
       deallocate(mt,mti,Apar_pre_antenna,amp,d_amp_dt)
    case(3)
       !we write out the amplitude status at the end of the simulation to allow smooth continuation
       if (.not.in_perf_opt) lv_antenna_amp(1:n_antennae)=amp
       call finalize_prefac
       deallocate(mt,mti,Apar_pre_antenna,amp,d_amp_dt)
    end select
    init_status=0
  end subroutine finalize_antenna

  subroutine add_Apar_antenna(emfields,a11det_inv_antenna)
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,n_fields),intent(inout):: emfields
    real,dimension(li1:li2,lj1:lj2,lk1:lk2,1:2),intent(in),optional:: a11det_inv_antenna
    complex:: A_ant
    real:: k_par1
    integer:: i,j,k,p

    Apar_pre_antenna = emfields(li1:li2,lj1:lj2,lbz:ubz,2)
    select case(antenna_type)
    case(1)
       k_par1 = 1.0
       DO k = lk1, lk2
          DO j = lj1, lj2
             DO i = li1, li2
                emfields(i,j,k,2) = emfields(i,j,k,2) - &
                     a11det_inv_antenna(i,j,k,1) * &
                     Apar0_antenna * EXP(imag*(k_par1*zval(k)-omega0_antenna*time))
             END DO
          END DO
       END DO
    case(2)
       do p=1,n_antennae
          if (j_ant(p).ge.lj1.and.j_ant(p).le.lj2.and..not.(i_ant(p).eq.0.and.j_ant(p)+ky0_ind.eq.0)) then
             DO k = lk1, lk2
                A_ant=amp(p)*exp(imag*k_par(p)*zval(k))
                Apar_pre_antenna(i_ant(p),j_ant(p),k)=emfields(i_ant(p),j_ant(p),k,2)+&
                     A_ant*a11det_inv_antenna(i_ant(p),j_ant(p),k,2)
                emfields(i_ant(p),j_ant(p),k,2) = emfields(i_ant(p),j_ant(p),k,2) + &
                     a11det_inv_antenna(i_ant(p),j_ant(p),k,1)*A_ant
                !reality condition
                if (j_ant(p)+ky0_ind.eq.0.and.i_ant(p).ne.0) then
                   Apar_pre_antenna(nx0-i_ant(p),j_ant(p),k)=emfields(nx0-i_ant(p),j_ant(p),k,2) + &
                        conjg(A_ant)*a11det_inv_antenna(nx0-i_ant(p),j_ant(p),k,2)
                   emfields(nx0-i_ant(p),j_ant(p),k,2) = emfields(nx0-i_ant(p),j_ant(p),k,2) + &
                        a11det_inv_antenna(nx0-i_ant(p),j_ant(p),k,1)*conjg(A_ant)
                endif
             END DO
          endif

       enddo
    case(3)
       do p=1,n_antennae
          if (j_ant(p).ge.lj1.and.j_ant(p).le.lj2.and..not.(i_ant(p).eq.0.and.j_ant(p)+ky0_ind.eq.0)) then
             DO k = lk1, lk2
                A_ant=amp(p)*exp(imag*k_par(p)*zval(k))
                emfields(i_ant(p),j_ant(p),k,2) = emfields(i_ant(p),j_ant(p),k,2) + A_ant
                !reality condition
                if (j_ant(p)+ky0_ind.eq.0.and.i_ant(p).ne.0) then
                   emfields(nx0-i_ant(p),j_ant(p),k,2) = emfields(nx0-i_ant(p),j_ant(p),k,2) + conjg(A_ant)
                endif
             END DO
          endif
       enddo
       
    end select
    
  end subroutine add_Apar_antenna
  
  subroutine calc_antenna_evolution
    complex:: omega_cplx
    real:: sigma
    integer:: p

    do p=1,n_antennae
       omega_cplx=lv_antenna_freq(p)
       
       !we need to use the initial amplitude to set the size of the random contribution,
       !and thus need to know the initial amplitude also for continuation runs
       if (in_perf_opt) then 
          !in perf_opt we drop the random part, since this requires a properly initialized dt,
          !which is not there when running this for the first time
          sigma=0.
       else
          sigma=abs(lv_antenna_initamp(p))*sqrt(12*abs(aimag(omega_cplx))/dt)
       end if
       
       !amplitude evolution consists of damped oscillation plus random part
       if (.not.in_perf_opt) then
          if (dt.ge.epsilon(dt)) then
             d_amp_dt(p)=amp(p)*(exp(-imag*omega_cplx*dt)-1.)/dt+rnd(p)*sigma
          else
             d_amp_dt(p)=0.
          endif
       else
          d_amp_dt(p)=0.
       endif
    enddo

  end subroutine calc_antenna_evolution

  subroutine evolve_antenna_amplitudes
    integer:: p

    call calc_antenna_evolution
    do p=1,n_antennae
       amp(p)=amp(p)+d_amp_dt(p)*dt
    enddo

  end subroutine evolve_antenna_amplitudes

  subroutine get_dApar_dt_antenna(p,dApar_dt,i,j)
    complex,dimension(lk1:lk2):: dApar_dt
    integer,intent(in):: p 
    integer,optional:: i,j
    integer:: k

    do k=lk1,lk2
       dApar_dt(k)=d_amp_dt(p)*exp(imag*k_par(p)*zval(k))
    enddo
    if (present(i).and.present(j)) then
       i=i_ant(p)
       j=j_ant(p)
    endif

  end subroutine get_dApar_dt_antenna

  subroutine get_Apar_antenna(p,A_ant,i,j)
    integer,intent(in):: p
    complex, dimension(lk1:lk2):: A_ant
    integer,optional:: i,j
    integer:: k

    DO k = lk1, lk2
       A_ant(k)=amp(p)*exp(imag*k_par(p)*zval(k))
    END DO
    if (present(i).and.present(j)) then
       i=i_ant(p)
       j=j_ant(p)
    endif
    
  end subroutine get_Apar_antenna


  subroutine initialize_prefac
    integer:: k,l,m,n,pni

    if (init_status.eq.0) then
       allocate(prefac(lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   prefac(k,l,m,n)=spec(n)%charge*vp(l)*sqrt(2./spec(n)%temp/spec(n)%mass)*fm(pi1,pj1,k,l,m,pni)
                end do
             end do
          end do
       end do
    endif

  end subroutine initialize_prefac

  subroutine finalize_prefac
    deallocate(prefac)
  end subroutine finalize_prefac

  subroutine add_dApar_dt_antenna(p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    integer:: klmn, p
    complex, dimension(lk1:lk2):: dApar_dt

    do p=1,n_antennae
       call get_dApar_dt_antenna(p,dApar_dt)
       if (j_ant(p).ge.lj1.and.j_ant(p).le.lj2.and..not.(i_ant(p).eq.0.and.j_ant(p)+ky0_ind.eq.0)) then
          do klmn=lb1,lb2
             p_rhs(i_ant(p),j_ant(p),klmn)=p_rhs(i_ant(p),j_ant(p),klmn)-&
                  prefac(sk(klmn),sl(klmn),sm(klmn),sn(klmn))*dApar_dt(sk(klmn))*&
                  jfac(i_ant(p),j_ant(p),sk(klmn),sm(klmn),sn(klmn))
          end do
       endif
       !reality condition
       if (j_ant(p)+ky0_ind.eq.0.and.i_ant(p).ne.0.and.p_has_0_mode) then
          do klmn=lb1,lb2
             p_rhs(nx0-i_ant(p),j_ant(p),klmn)=p_rhs(nx0-i_ant(p),j_ant(p),klmn)-&
                  prefac(sk(klmn),sl(klmn),sm(klmn),sn(klmn))*conjg(dApar_dt(sk(klmn)))*&
                  jfac(nx0-i_ant(p),j_ant(p),sk(klmn),sm(klmn),sn(klmn))
          end do
       endif
    enddo
    
  end subroutine add_dApar_dt_antenna


  function rnd(p)
    complex:: rnd
    integer, intent(in):: p
    integer, dimension(0:nrnd-1):: tmparr
    integer:: tmpint

    tmparr=mt(p,:)
    tmpint=mti(p)
    rnd=grnd(tmparr,tmpint)-0.5+imag*(grnd(tmparr,tmpint)-0.5)
    mt(p,:)=tmparr
    mti(p)=tmpint
    return
  end function rnd
    
  subroutine seed_rnd(p,initamp)
    integer, intent(in):: p
    complex, intent(in):: initamp
    integer:: seed, tmpint
    integer, dimension(0:nrnd-1):: tmparr

    tmparr=mt(p,:)
    tmpint=mti(p)
    seed=int(abs(initamp)*1000*p)
    call sgrnd(seed,tmparr,tmpint)
    mt(p,:)=tmparr
    mti(p)=tmpint
  end subroutine seed_rnd

end module antenna
