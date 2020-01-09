!>Sets initial value of g_1
#include "redef.h"
module initcond 
  use par_mod
  use communications
  use aux_func, only: gamma0, j0
  use geometry
  USE gyro_average_ff_mod, only: get_jfac
  USE aux_fields
  use fourier, only: initialize_fourier_x_1d, to_real_x_1d, to_fourier_x_1d, &
    finalize_fourier_x_1d, initialize_fourier, finalize_fourier
  USE parallel_nonlin,only: nl_to_fourier_xy
  USE discretization
  USE vel_space
  use neo_equilibrium
  use antenna, only: get_Apar_antenna, n_antennae, antenna_type
  use mtrandom
  ! adaptivity module necessary for nullify_g1_outside_adptv subroutine
  use discretization_adptv_module

  Implicit None

  public:: init_g_1, check_init_cond
  public:: init_cond, init_aux_x, init_aux_y, init_aux_z, init_aux_amp
  public:: init_vpar_fluct

  ! block structured grid subroutines
  public :: nullify_g1_outside_adptv
  
  private
  
  character(len=20) :: init_cond = ''
  real   :: init_aux_x = -100., init_aux_y = -100., init_aux_z = -100.
  real   :: init_aux_amp = -100.0, init_vpar_fluct=-100.0

  !local variables
  LOGICAL:: cond_local_in_i=.true., write_pe

contains 

  subroutine check_init_cond
    ! set initial condition if not found in parameter file
    IF (init_cond.EQ.'') THEN
       IF (xy_local) THEN
          IF (nonlinear) THEN
             init_cond = 'ppj'
          ELSE
             init_cond = 'alm'
          END IF
       ELSE
          IF ((.NOT.y_local).AND.(.not.nonlinear)) THEN
             init_cond = 'mmj'
          ELSE
             init_cond = 'db'
          ENDIF
       ENDIF
    END IF


  end subroutine check_init_cond

  Subroutine init_g_1

    !Local variables
    REAL :: local_sum, global_sum
    integer :: ierr
    LOGICAL :: OUTPUT=.FALSE.

    write_pe = (mype==0).and.(print_ini_msg)

    If (write_pe) Then 
       Write(*,*)
       Write(*,"(A,A)") "initialization: ", init_cond
    End If
    Select case (init_cond)
      case('sw') ! single mode in kx, ky and z
         Call StartWave(init_aux_x, init_aux_y, init_aux_z, init_aux_amp)
      case('fb') ! gaussians in kx, ky and z
         Call fourier_blob(init_aux_x, init_aux_y, init_aux_z)
      case('db') ! gaussians in x, ky and z
         Call density_blob(init_aux_x, init_aux_y, init_aux_z, init_aux_amp)
      case('db1')
         Call density_blob(init_aux_x, init_aux_y, init_aux_z, init_aux_amp)
      case('alm') ! all Fourier modes in x- and y-direction, Gaussian in parallel direction
         Call all_modes(init_aux_z,init_aux_amp,.false.,.false.)
      case('almmt') !microtearing mode
         Call all_modes(init_aux_z,init_aux_amp,.false.,.true.)
      case('almnf') ! all modes, but phi zero
         Call all_modes(init_aux_z,init_aux_amp,.true.,.false.)
      case('ppg') ! Power law in kx and ky; gaussian distribution in z
         Call powgauss(init_aux_x, init_aux_y, init_aux_z,init_aux_amp)
      case('ppj') ! Power law in kx and ky; powers of jacobian in z
         Call powjac(init_aux_x, init_aux_y, init_aux_z, init_aux_amp,.false.,.false.)
      case('ppjrn') ! Power law in kx and ky; powers of jacobian in z; 
         !random phases in kx and ky
         Call powjac(init_aux_x, init_aux_y, init_aux_z, init_aux_amp,.true.,.false.)
      case('ppjmt') ! Power law in kx and ky; powers of jacobian in z, odd in z
         Call powjac(init_aux_x, init_aux_y, init_aux_z, init_aux_amp,.false.,.true.)
      case('pmtrn') ! Power law in kx and ky; powers of jacobian in z, odd in z
         !random phases in kx and ky
         Call powjac(init_aux_x, init_aux_y, init_aux_z, init_aux_amp,.true.,.true.)
      case('mmj') ! Single modes in kx and ky; powers of jacobian in z
         Call modejac(init_aux_x, init_aux_y, init_aux_z,init_aux_amp)
      case('cosn') ! Cosine in x and y + white noise
         call cosnoise(init_aux_x, init_aux_y, init_aux_z,init_aux_amp)
      case('gmj') ! Single modes in kx and ky; powers of jacobian in z
         Call globalmodejac(init_aux_x, init_aux_y, init_aux_z)
      case('gaminit_y')
         CALL gaminit_y(init_aux_x,init_aux_z)  
      case('gam') ! Single modes in kx and ky; powers of jacobian in z
         Call gaminit(init_aux_x,init_aux_z)
      case('cs') ! current sheets
         Call current_sheets(init_aux_x,init_aux_y,init_aux_z,init_aux_amp,.false.)
      case('ffcs') ! force free current sheet setup
         Call current_sheets(init_aux_x,init_aux_y,init_aux_z,init_aux_amp,.true.)
      case('plasmoid') ! one or two reconnection plasmoids
         CALL plasmoids(init_aux_x,init_aux_y,init_aux_z)
      case('zero') ! everything zero (use with antenna currect)
         Call antenna_zero
      case default
       If(mype==0) Write(*,"(3(A))") "!! initial condition ", &
            &TRIM(init_cond)," does not exist !!"
       Stop
    End select
    If (write_pe) Write(*,*)

    time = 0.0

    IF (OUTPUT) THEN
       ! Output of the square sum of the initialized g_1
       local_sum = REAL(SUM(g_1*CONJG(g_1)))
       WRITE(*,"(I3,A,ES20.12)") mype,": local_sum(g_1) = ", local_sum
       CALL mpi_reduce(local_sum, global_sum,1, MPI_REAL_TYPE,MPI_SUM,0,MY_MPI_COMM_WORLD,ierr)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "square sum of g_1 after initialization is ",global_sum
    END IF
    
    !avoid underflow exit due to g=0
    if ((init_cond.eq.'zero').or.&
          &((maxval(ki).le.epsilon(maxval(ki))).and.(kjmax.le.epsilon(kjmax)))) then
       if (underflow_limit.ge.0) then
          underflow_limit = -1
          !if (mype==0) Write(*,"(A)") "disabled underflow exit by setting underflow_limit=-1"
       endif
    endif

  end Subroutine init_g_1
  
  !*************************************************************************!

  !> Fourier blob initialization:
  !! Gaussians in kx, ky and z
  Subroutine fourier_blob(aux_x, aux_y, aux_z)
    REAL, intent(IN) :: aux_x, aux_y, aux_z
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real::   amplitude,sigma_i,sigma_j,sigma_z
    Integer::    i, j, k, n

    sigma_i = maxval(ki)/4.0
    sigma_j = kjmax/4.0
    sigma_z = pi/4.0

    !if only 0 modes: avoid 0./0.
    if (maxval(ki).le.epsilon(maxval(ki))) sigma_i = 1.0
    if (kjmax.le.epsilon(kjmax)) sigma_j = 1.0

    IF (.not.yx_order) THEN
       If ((aux_x.ne.-100.).And.(aux_x.ne.0.0)) sigma_i=sigma_i/aux_x
       If ((aux_y.ne.-100.).And.(aux_y.ne.0.0)) sigma_j=sigma_j/aux_y
    ELSE
       If ((aux_x.ne.-100.).And.(aux_x.ne.0.0)) sigma_i=sigma_i/aux_y
       If ((aux_y.ne.-100.).And.(aux_y.ne.0.0)) sigma_j=sigma_j/aux_x
    ENDIF
    If ((aux_z.ne.-100.).And.(aux_z.ne.0.0)) sigma_z=sigma_z/aux_z

    If (write_pe) Then
       Write(*,"(A)") "  fourier blob (gaussians in kx, ky and z)"
       Write(*,"(3(A,ES10.3))") "  sigma_i=",sigma_i,", sigma_j=",sigma_j,&
            & ", sigma_z=", sigma_z
    End If
    
    ! initial amplitude
    amplitude=0.1

    Do k = lk1,lk2
       Do j = lj1,lj2
          Do i = li1,li2
             dens(i,j,k) = amplitude*exp(-(ki(i)/sigma_i)**2/2.0&
                  &-(kj(j)/sigma_j)**2/2.0-(zval(k)/sigma_z)**2/2.0)
          End Do
       End Do

       if (p_has_00_mode) dens(li1,lj1,k) = 0.0
    End Do

    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  End Subroutine fourier_blob

  !*************************************************************************!

  !> Start Wave initialization:
  !! Single mode in kx, ky and kpar (i.e., shat = 0 only)
  Subroutine StartWave(aux_x,aux_y,aux_z,aux_amp)
    Real,Intent(IN) :: aux_x,aux_y,aux_z,aux_amp
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Integer :: osi,osj,k,n, minus_osi
    Real:: p,kpar

    osi = 1
    osj = 1
    kpar = 0.0
    p = 0.1
    
    IF (.not.yx_order) THEN
       IF (aux_x.ne.-100.) osi = INT(aux_x)
       IF (aux_y.ne.-100.) osj = ABS(INT(aux_y))
    ELSE
       IF (aux_y.ne.-100.) osi = INT(aux_y)
       IF (aux_x.ne.-100.) osj = ABS(INT(aux_x))
    ENDIF
    IF (osi.lt.0) osi = li2+osi+1

    If (aux_z.ne.-100.) kpar = aux_z
    IF (aux_amp.ne.-100.0) p = aux_amp

    If (write_pe) Then       
       Write(*,"(A)") "  single mode with"
       IF (.not.yx_order) then
          Write(*,"(3(A,ES10.3))") "  kx=",kx(osi),", ky=",ky(osj),&
               &", kpar= ",kpar
       else
          Write(*,"(3(A,ES10.3))") "  kx=",kx(osj),", ky=",ky(osi),&
                         &", kpar= ",kpar
       endif
       If (shat.ne.0.0) Then 
          Write(*,"(A)") "  !! makes no sense with shat <> 0 !!"
          STOP
       End If
    End If

    dens = 0.
    If ((my_pey == osj/lj0)) then
       Do k = lk1,lk2
          dens(osi, lj1+Modulo(osj,lj0), k) = p*&
               &exp(cmplx(0.0,1.0)*kpar*zval(k))
       End Do
    Endif

    IF ((osj.eq.0).and.(abs(osi).gt.0)) THEN !fulfill reality condition
       minus_osi = li2+1-osi
       dens(minus_osi,lj1+Modulo(osj,lj0),:) = &
            CONJG(dens(osi,lj1+Modulo(osj,lj0),:))
    ENDIF
    
    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  End Subroutine Startwave

  !*************************************************************************!

  !> All modes initialization
  !! All kx, ky modes are equally excited, the jacobian to the power of 
  !! aux_z is taken in the z direction
  !! odd_z_parity sets additional z-dependence (kpar=1)
  Subroutine all_modes(aux_z,aux_amp,no_field,odd_z_parity)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_z, aux_amp
    Logical, INTENT(in) :: no_field !< divide g_1 by J0 if true
    logical, intent(in) :: odd_z_parity
    Integer::   i, k, n
    Real :: p, es

    if (xy_local) then
       es = 2.0
    else
       es = 0.0 !avoid jacobian which is not in Fourier space in 1st dim.
    endif
    p = 0.1

    If (aux_z.ne.-100.) es=aux_z
    If (aux_amp.ne.-100.) p = aux_amp

    If (write_pe) Then
       Write(*,"(A)") "  all kx, ky modes equally excited,"
       if (odd_z_parity) then
          Write(*,"(A,ES10.3,A)") "sin(z) * jacobian to the power of ", es, " in z"
       else
          Write(*,"(A,ES10.3,A)") "  jacobian to the power of ", es, " in z"
       endif
    End If


    Do k = lk1,lk2
       if (pmi0.eq.1) then
          dens(:,:,k) = p*cmplx(1.0,0.)*(geom%jacobian(pi1,pj1,k)/geom%avg_jaco)**es
       else
          Do i=li1,li2
             dens(i,:,k) = p*cmplx(1.0,0.)*(geom%jacobian(i,pj1,k)/geom%avg_jaco)**es
          End do
       endif
       if (odd_z_parity) dens(:,:,k) = dens(:,:,k)*cmplx(sin(zval(k)-epsilon(zval(lk1))),0.D0)
       !shift by eps to avoid exact zero (problem for diag_omega)

       !set 00 fourier mode to zero(if this is not the only mode in the system)
       if (p_has_00_mode.and.(nx0*nky0.gt.1)) dens(li1,lj1,k)=0.0
    End do

    !weaken zonal flow:
    !    if (p_has_0_mode) dens(:,lj1,:) = 0.01 * dens(:,lj1,:)

    
    do n=ln1,ln2
       if (no_field) then
          Call add_fm(dens,n,-1.0)
       else
          Call add_fm(dens,n)
       endif
    end do

  End Subroutine all_modes




  !*************************************************************************!

  !> Powgauss initialization
  !! Power laws in kx and ky, Gaussian in z direction
  SUBROUTINE powgauss(aux_x,aux_y,aux_z,aux_amp)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_x,aux_y,aux_z,aux_amp
    Real ::    amplitude, ex, ey, sigma_z
    Integer :: i, j, k, n, ikx, iky

    ex = -1.0 
    ey = -1.0
    sigma_z = pi/4.0
    amplitude = 0.1

    If (aux_x.ne.-100.) ex=-Abs(aux_x) 
    If (aux_y.ne.-100.) ey=-Abs(aux_y)
    If ((aux_z.ne.-100.).AND.(aux_z==0.0)) sigma_z = pi/(4.0*aux_z)
    if (aux_amp.ne.-100.) amplitude = abs(aux_amp)

    If (write_pe) Then
       Write(*,"(A)") "  power law in kx and ky, gaussian in z"
       Write(*,"(3(A,ES10.3))") "  kx to the power of ", ex, &
               & ", ky to the power of ", ey, &
               & ", sigma_z = ", sigma_z
    End If
    
    do k = lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
             if (yx_order) then
                iky = i
                ikx = j
             else
                ikx = i
                iky = j
             endif
             if (ky(iky).lt.epsilon(ky(iky))) then
                if (abs(kx(ikx)).lt.epsilon(kx(ikx))) then
                   dens(i,j,k) = 0.0
                else
                   dens(i,j,k) = 0.5*amplitude*((abs(kx(ikx))+kxmin)/kxmin)**ex&
                        & *(0.5+0.5*exp(-(zval(k)/sigma_z)**2/2.0) &
                        & *sqrt(2.0*sqrt(pi)/sigma_z))
                endif
             else
                dens(i,j,k) = amplitude*((abs(kx(ikx))+kxmin)/kxmin)**ex &
                     & *(ky(iky)/kymin)**ey&
                     & *(0.5+0.5*exp(-(zval(k)/sigma_z)**2/2.0) &
                     & *sqrt(2.0*sqrt(pi)/sigma_z))
                if (abs(kx(ikx)).lt.epsilon(kx(ikx))) &
                     & dens(i,j,k) = 0.5*amplitude*(ky(iky)/kymin)**ey&
                     & *(0.5+0.5*exp(-(zval(k)/sigma_z)**2/2.0) &
                     & *sqrt(2.0*sqrt(pi)/sigma_z)) 
             end if
          end do
       end do
    end do
    
    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  END SUBROUTINE powgauss

  !*************************************************************************!

  !> Powjac initialization
  !! Power laws in kx and ky, Jacobian to the power of aux_z in z direction
  SUBROUTINE powjac(aux_x,aux_y,aux_z, aux_amp,rnd_phase,odd_z_parity)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_x,aux_y,aux_z,aux_amp
    Logical,Intent(IN) :: rnd_phase,odd_z_parity
    complex, dimension(0:ni0-1,lj1:lj2):: phasefac
    Real ::    amplitude, ex, ey, es
    Integer :: i, j, k, n, pii, ikx, iky

    ex = -1.0 
    ey = -1.0
    es =  2.0
    amplitude = 1.0

    If (aux_x.ne.-100.) ex=-Abs(aux_x) 
    If (aux_y.ne.-100.) ey=-Abs(aux_y)
    If (aux_z.ne.-100.) es=Abs(aux_z)
    if (aux_amp.ne.-100.) amplitude = abs(aux_amp)

    If (write_pe) Then       
       Write(*,"(A)") "  power law in kx and ky, powers of jacobian in z"
       Write(*,"(3(A,ES10.3))") "  kx to the power of ", ex, &
            & ", ky to the power of ", ey, &
            & ", jacobian to the power of ", es
       if (rnd_phase) write(*,"(A)") "  random phases in kx, ky"
    End If
 
    pii = pi1gl

    if (rnd_phase) then
       call sgrnd(1)
       if (yx_order) then
          do i=0,ni0-1
             do j=lj1,lj2
                if ((abs(ki(i)).lt.epsilon(amplitude)).and.&
                     &(kj(j).lt.0.0)) then !fulfill reality condition (ki half space)
                   phasefac(i,j) = conjg(phasefac(i,nj0-j))
                elseif ((abs(kj(j)).lt.epsilon(amplitude)).and.&
                     &(ki(i).lt.0.0)) then !fulfill reality condition (kj half space <- y global)
                   phasefac(i,j) = conjg(phasefac(ni0-i,j))
                else
                   phasefac(i,j)=exp(imag*grnd()*2*pi)
                endif
             enddo
          enddo
       else
          do j=lj1,lj2
             do i=0,ni0-1
                if ((abs(kj(j)).lt.epsilon(amplitude)).and.&
                     &(ki(i).lt.0.0)) then !fulfill reality condition
                   phasefac(i,j) = conjg(phasefac(ni0-i,j))
                else
                   phasefac(i,j) = exp(imag*grnd()*2*pi)
                endif
             enddo
          enddo
       endif
    else
       phasefac = 1.0
    endif

    do k = lk1,lk2
       do j=lj1,lj2
          do i=li1,li2
             if (yx_order) then
                iky = i
                ikx = j
             else
                ikx = i
                iky = j
             endif
             if (pmi0.gt.1) pii=i
             if (abs(ky(iky)).lt.epsilon(ky(iky))) then !ky=0 mode
                if (abs(kx(ikx)).lt.epsilon(kx(ikx))) then !kx=ky=0 mode
                   dens(i,j,k) = 0.0
                else
                   dens(i,j,k) = 0.5*amplitude*((abs(kx(ikx))+kxmin)/kxmin)**ex*phasefac(i,j)
                endif
             else
                dens(i,j,k) = amplitude*((abs(kx(ikx))+kxmin)/kxmin)**ex &
                     & *(ky(iky)/kymin)**ey*phasefac(i,j)
                if (abs(kx(ikx)).lt.epsilon(kx(ikx))) &
                     & dens(i,j,k) = 0.5*amplitude*(ky(iky)/kymin)**ey*phasefac(i,j)
             end if
             !add k dependence
             dens(i,j,k) = dens(i,j,k)*(geom%jacobian(pii,pj1,k)/geom%avg_jaco)**es
             if (odd_z_parity) &
                  dens(i,j,k) = dens(i,j,k)*cmplx(sin(zval(k)-epsilon(zval(0))),0.D0)
             !shift by eps to avoid exact zero (problem for diag_omega)
          end do
       end do
    end do

    do n=ln1,ln2
       call add_fm(dens,n)
    end do

  END SUBROUTINE powjac




  !*************************************************************************!

  !> Tatsuno initialization
  !! White noise plus cosines in x and y
  subroutine cosnoise(aux_x,aux_y,aux_z, aux_amp)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_x,aux_y,aux_z,aux_amp
    Real ::    amplitude
    integer :: i, j, k, n, ex, ey, es
    !noise array is not parallelized in y to guarantee identical results
    !also with y-parallelization
    complex, dimension(li1:li2,0:nky0-1):: noise

    ex = 2
    ey = 2
    !es is currently not used
    es =  1
    amplitude = 1.0e-4

    If (aux_x.ne.-100.) ex=Abs(aux_x) 
    If (aux_y.ne.-100.) ey=Abs(aux_y)
    If (aux_z.ne.-100.) es=Abs(aux_z)
    if (aux_amp.ne.-100.) amplitude = abs(aux_amp)

    If (write_pe) Then       
       Write(*,"(A)") " cosine in x and y plus white noise "
    End If
 
    
    call sgrnd(1)
    noise=0.
    !reality condition for noise
    do i=1,hkx-li1
       noise(li1+i,0)=amplitude*grnd()+amplitude*imag*grnd()
       noise(li2-i+1,0)=conjg(noise(i,0))
    enddo
    do j=1,nky0-1
       do i=li1,li2
          noise(i,j)=amplitude*grnd()+amplitude*imag*grnd()
       enddo
    enddo

    dens=0.
    do k = lk1,lk2
       if (p_has_0_mode) then
          dens(li1+ex,lj1,k)=1.0
          dens(li2-ex+1,lj1,k)=1.0
       endif
       do j=lj1,lj2
          if (j==ey) dens(0,j,k)=1.0
          do i=li1,li2
             dens(i,j,k)=dens(i,j,k)+noise(i,j)
          enddo
       enddo
    enddo
    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  end subroutine cosnoise




  !*************************************************************************!

  !> modejac initialization
  !! Single mode (default: kx=0, ky=kymin) excitation;
  !! power of Jacobian in z direction
  Subroutine modejac(aux_x,aux_y,aux_z,aux_amp)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_x,aux_y,aux_z,aux_amp
    Integer :: ikx,iky, iki, ikj
    Real:: amplitude, es
    Integer::   k, n, pii

    ikx = 0
    iky = 1
    If ((ky0_ind.gt.0).or.(.not.y_local)) iky=0
    es  = 2.0
    amplitude = 0.1

    If (aux_x.ne.-100.) ikx=INT(aux_x)
    If (aux_y.ne.-100.) iky=INT(aux_y)
    If (aux_z.ne.-100.) es=INT(aux_z)
    If (aux_amp.ne.-100.) amplitude=aux_amp

    if (yx_order) then
       iki = iky
       ikj = ikx
    else
       iki = ikx
       ikj = iky
    endif

    if (iki.lt.0) iki=li2+iki+1
    
    dens = 0.0 !1E-5 * amplitude

    If (write_pe) Then
       Write(*,"(A)") "  single mode in kx and ky, power of jacobian in z"
       WRITE(*,"((A,I3),(A,ES10.3),(A,I3),2(A,ES10.3))") "  kx(",ikx,")=",&
               & kx(ikx),", ky(",iky,")=",ky(iky),&
               &", jacobian to the power of ",es
    End If

    If ((my_pey == ikj/lj0)) then
       if (pmi0.eq.1) then
          pii=pi1
       else
          pii=li1+iki
       endif
       Do k = lk1,lk2
          dens(li1+iki,lj1+MODULO(ikj,lj0), k) = amplitude*geom%jacobian(pii,pj1,k)**es
          !Fulfill reality condition
          IF ((ikj.eq.0).and.(iki.ne.0)) dens(li2-iki+1,lj1, k) = CONJG(dens(li1+iki,lj1, k))
       End Do
    Endif

    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  End Subroutine modejac

  !*************************************************************************!

  !> globalmodejac initialization
  !! Single mode (default: kx=0, ky=kymin) excitation;
  !! power of Jacobian in z direction
  !! similar to modejac but implemented in direct space in the first dimension
  Subroutine globalmodejac(aux_x,aux_y,aux_z)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    Real,Intent(IN) :: aux_x,aux_y,aux_z
    Integer :: ikx, iky, kii,kji
    Real:: amplitude, es, xvar !, kxmin
    Integer::   i, k, n

    cond_local_in_i = .false.

    ikx = 0
    iky = 1
    If ((ky0_ind.gt.0).or.(.not.y_local)) iky=0
    es  = 2.0
    amplitude = 1.0

    If (aux_x.ne.-100.) ikx=INT(aux_x)
    If (aux_y.ne.-100.) iky=INT(aux_y)
    If (aux_z.ne.-100.) es = aux_z

    if (yx_order) then
       kii = iky
       kji = ikx
    else
       kii = ikx
       kji = iky
    endif

    if (kii.lt.0) kii=li2+kii+1

    dens = 0.0 !1E-5 * amplitude
!    kxmin=2.0*pi/lx

    If (write_pe) Then
       Write(*,"(A)") "  single mode in i and kj, power of jacobian in z"
       WRITE(*,"((A,I3),(A,ES10.3),(A,I3),2(A,ES10.3))") "  kx(",ikx,")=",&
               & kx(ikx),", ky(",iky,")=",ky(iky),&
               &", jacobian to the power of ",es
    End If

    If ((my_pey == kji/lj0)) then
       Do k=lk1,lk2
          Do i=li1,li2
             xvar = -0.5*lp1+(i-li1)*deli
             dens(i, lj1+MODULO(kji,lj0), k) = amplitude*geom%Bfield(i,pj1,k)**es*&
                  !geom%jacobian(i,pj1,k)**es*&
                  &COS(ki(kii)*xvar)
          Enddo
       Enddo
    Endif

    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  End Subroutine globalmodejac

  !*************************************************************************!

  !> density blob initialization
  !! Gaussian in direct space in first and third dimension, Fourier
  !! space in second dimension
  Subroutine density_blob(aux_x,aux_y,aux_z,aux_amp)
    REAL, intent(IN) :: aux_x, aux_y, aux_z, aux_amp
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    REAL::   amplitude,sigma_i,sigma_j,sigma_z,xvar
    INTEGER:: i, j, k, n

    cond_local_in_i = .false.

    sigma_i = lp1/8.0
    if(abs(kjmax).gt.epsilon(kjmax)) then
       sigma_j = kjmax/4.0
    else
       sigma_j = 1.
    endif
    sigma_z = pi/4.0
    amplitude = 1.0D0

    If ((aux_x.ne.-100.).And.(aux_x.ne.0.0)) sigma_i=sigma_i/aux_x
    If ((aux_y.ne.-100.).And.(aux_y.ne.0.0)) sigma_j=sigma_j/aux_y
    If ((aux_z.ne.-100.).And.(aux_z.ne.0.0)) sigma_z=sigma_z/aux_z
    If (aux_amp.ne.-100.) amplitude = aux_amp
    
    If (write_pe) Then
       Write(*,"(A)") "  density blob (gaussians in x, ky<>0 and z)"
       Write(*,"(3(A,ES10.3))") "  sigma_i=",sigma_i,", sigma_j=",sigma_j,&
            & ", sigma_z=", sigma_z
    End If
    
    dens = CMPLX(0.0,0.0)

    Do k = lk1,lk2
       Do j = lj1,lj2
          Do i = li1,li2
             xvar = -0.5*lp1 + i*deli
             dens(i,j,k) = amplitude*exp(-(xvar/sigma_i)**2/2.0&
                  &-(kj(j)/sigma_j)**2/2.0-(zval(k)/sigma_z)**2/2.0)
          End Do
       End Do
    End Do
    if (y_local) then
       If (p_has_0_mode.and.(abs(kjmax).gt.epsilon(kjmax))) then
          do k=lk1,lk2
             dens(:,lj1,k) = 0.01*dens(:,lj1,k)/lp1**2
          enddo
       endif
    endif
    do n=ln1,ln2
       Call add_fm(dens,n)
    end do

  End Subroutine density_blob


  !*************************************************************************!

  !> gaminit_y initialization
  !! Initializes only a certain kx and ky=0 mode and a power of the jacobian
  !! the density is given by (1-Gamma0)*amplitude and hence the electrostatic
  !! potential is initialized to amplitude which is desirable for GAM
  !! investigations.  This is intended for benchmarking of the y-global code. 

  Subroutine gaminit_y(aux_x,aux_z)
    REAL, intent(IN) :: aux_x,aux_z
    real :: b2spec
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2,lm1:lm2,ln1:ln2) :: jfac_cmplx
    INTEGER:: i, j, k, l,m, n, kxind, pii, pji, pni
    REAL :: amplitude, k_perp2, z_exp, rho
    LOGICAL :: adiabatic_ions = .true.

    amplitude = 1.0

    If (aux_x.eq.-100.) then
       if (kx_center.ne.0.0) then
          kxind = 0
       else
          kxind = 1
       endif
    else
       kxind = int(aux_x)
    endif

    if (aux_z.eq.-100) then
       z_exp=0.0
    else
       z_exp = aux_z
    endif

    If (write_pe) Then
       Write(*,"(A)") "  GAM init. cond. ((1-Gamma0(kperp2,z,species)*amplitude))"
       Write(*,"(2(A,ES10.3))") "  kx=",kx(kxind), &
            &' jacobian to the power of ',z_exp
    End If
    
    do n=0,n_spec-1
       if (spec(n)%charge.gt.0) adiabatic_ions = .false.
    enddo

    pii = pi1
    pji = pj1
    pni = pn1

    do n=ln1,ln2
       if (pn0.gt.1) pni=n

       dens = cmplx(0.0,0.0)

       do i=li1,li2
          j = kxind
          do k=lk1,lk2
             k_perp2 = geom%gjj(pi1,pj1,k)*kj(j)**2
             b2spec = k_perp2*spec(n)%mass*spec(n)%temp/&
                  &(spec(n)%charge**2*geom%Bfield(pi1,pj1,k)**2)
             dens(i,j,k) = spec(n)%charge/spec(n)%temp*(1.0-gamma0(b2spec))*amplitude*&
                  (geom%jacobian(pi1,pj1,k)/geom%avg_jaco)**z_exp
             if (adiabatic_ions) dens(i,j,k) = dens(i,j,k) - tau*&
                  &amplitude*(geom%jacobian(pi1,pj1,k)/geom%avg_jaco)**z_exp
          enddo
       enddo
       
       Do m=lm1,lm2
          Do l= ll1,ll2
             Do k = lk1, lk2
                Do j = lj1,lj2
                   if (pj0.gt.1) pji=j
                   Do i = li1, li2
                      if (pi0.gt.1) pii=i
                      k_perp2 = geom%gjj(pi1,pj1,k)*kj(j)**2
                      rho = sqrt(2.0*spec(n)%mass*spec(n)%temp*mu(m)/&
                           &spec(n)%charge**2/geom%Bfield(i,pj1,k))
                      jfac_cmplx(i,j,k,m,n)=cmplx(j0(sqrt(k_perp2)*rho),0.0)
                      g_1(i,j,k,l,m,n) = dens(i,j,k)*fm(pii,pji,k,l,m,pni)/jfac_cmplx(i,j,k,m,n)
                   enddo
                enddo
             enddo
          endDo
       endDo
    enddo
    
  End Subroutine gaminit_y



  !*************************************************************************!

  !> gaminit initialization (for GAM investigations)
  !! Initializes only ky=0 and a certain kx mode (either just a single kx 
  !! grid point if nx0=1 or the plus/minus components) and a power of the jacobian;
  !!
  !! the quasi-neutrality equation is employed to initialize the density such
  !! that phi=amplitude:
  !! phi = sum_j n0j qj int J0 g1j dv dmu / sum_s n0s qs^2/Ts (1-Gamma0)
  !!     = amplitude
  !! While previous version of this routine considered
  !! g1j = qj/Tj (1-Gamma0) F0j / J0
  !! which caused problems at higher k where J0~0, we now use 
  !! g1j = qj/Tj (1-Gamma0)/Gamma0 F0j * J0
  Subroutine gaminit(aux_x,aux_z)
    REAL, intent(IN) :: aux_x,aux_z
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    INTEGER:: i, j, k, n, kxind, pii, pji, ii, uii=0, i1, j1
    REAL :: amplitude, k_perp2, b2spec, z_exp
    LOGICAL :: adiabatic_ions = .true.

    z_exp = 0.0
    amplitude = 1.0

    If (aux_x.eq.-100.) then
       if (kx_center.ne.0.0) then
          kxind = 0
       else
          kxind = 1
       endif
    else
       kxind = int(aux_x)
       if (nx0.lt.(2*abs(kxind)+1)) &
            &stop 'nx0>=2*init_aux_x+1 required in gaminit'
    endif
    if (aux_z.ne.-100) z_exp = aux_z


    if (nx0.gt.2) uii=1

    If (write_pe) Then
       Write(*,"(A)") "  GAM IC (phi=amplitude)"
       Write(*,"(2(A,ES10.3))") "  kx=",kx(kxind), &
            &', jacobian to the power of ',z_exp
    End If    

    pii = pi1
    pji = pj1

    do n=0,n_spec-1
       if (spec(n)%charge.gt.0) adiabatic_ions = .false.
    enddo

    do n=ln1,ln2
       dens = cmplx(0.0,0.0)

       if (p_has_0_mode) then
          i = kxind
          j = lj1
          do ii=0,uii
             if (ii.gt.0) then
                i=nx0-kxind
             endif
             if (yx_order) then
                i1 = j; j1 = i
             else
                i1 = i; j1 = j
             endif

             do k=lk1,lk2
                k_perp2 = geom%gii(pi1,pj1,k)*ki(i1)**2+2.0*geom%gij(pi1,pj1,k)*&
                     &ki(i1)*kj(j1)+geom%gjj(pi1,pj1,k)*kj(j1)**2
                
                b2spec = k_perp2*spec(n)%mass*spec(n)%temp/&
                     &(spec(n)%charge*geom%Bfield(pi1,pj1,k))**2

                dens(i1,j1,k) = spec(n)%charge/spec(n)%temp*(1.0-gamma0(b2spec))
                
                if (adiabatic_ions) dens(i1,j1,k) = dens(i1,j1,k) - tau
                
                dens(i1,j1,k) = dens(i1,j1,k)/gamma0(b2spec)*&
                     &amplitude*(geom%jacobian(pi1,pj1,k)/geom%avg_jaco)**z_exp
             enddo
          enddo
       endif

       Call add_fm(dens,n,1.0)

    enddo !n loop
    
  End Subroutine gaminit

  !*************************************************************************!

  !> current sheets initialization
  !! (counter-flowing species)
  !! note: only implemented for local code!
  SUBROUTINE current_sheets(aux_x,aux_y,aux_z,aux_amp,force_free)

    IMPLICIT NONE

    REAL, INTENT(IN) :: aux_x, aux_y, aux_z, aux_amp
    logical, intent(in):: force_free
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: dens
    COMPLEX, DIMENSION(li1:li2) :: temp_x
    REAL :: amplitude, relative_pert, vel_shift, vel_shift_l, lmod, Gauss_offset_corr
    REAL :: fm_lmod, dv, nx0_inv, pi_n15, ifrac, bitsum, ffcs_width
    INTEGER :: i, j, k, l, m, n, x_mode
    LOGICAL :: vmult_Maxw, x_v_Gaussian, nonrandom_phases
    INTEGER :: rseed_size
    INTEGER, DIMENSION(:), ALLOCATABLE :: rseed
    REAL, DIMENSION(li1:li2) :: phases

    ! check: code needs to be run in local mode
    IF (.NOT. xy_local) THEN
      IF (write_pe) WRITE(*,"(A)") &
        "ERROR: current sheet initial condition requires local code"
      STOP
    END IF

    !two different interpretations of aux_x, depending on choice of current sheet setup!
    if (.not.force_free) then
       ! aux_x contains three bit-wise switches, plus 2-digit x_mode parameter
       !   to set the latter to, e.g., 5, set aux_x = [].05
       vmult_Maxw = .false. ! bit 1; 0 = shifted MW, 1 = v_par-multiplied MW
       x_v_Gaussian = .false. ! bit 2, for vmult_Maxw = 0; 0 = consine, 1 = Gaussian
       nonrandom_phases = .false. ! bit 3, for vmult_Maxw = 0; 0 = random, 1 = nonrandom
       bitsum = 0
       x_mode = 1 ! >1 for >1 cosine; for vmult_Maxw = 0, x_v_Gaussian = 0; see above
       IF ((aux_x .NE. -100.0) .AND. (aux_x .NE. 0.0)) bitsum = aux_x
       IF (bitsum .GE. 4.0) THEN
          bitsum = bitsum - 4.0
          nonrandom_phases = .true.
       END IF
       IF (bitsum .GE. 2.0) THEN
          bitsum = bitsum - 2.0
          x_v_Gaussian = .true.
       END IF
       IF (bitsum .GE. 1.0) THEN
          bitsum = bitsum - 1.0
          vmult_Maxw = .true.
       END IF
       IF (bitsum .GT. 0.015) x_mode = NINT(100.0*bitsum)
    else
       vmult_Maxw = .false.
       x_v_Gaussian = .false.
       nonrandom_phases = .false.
       !with force-free current sheet, aux_x simply defines the width of the current sheet
       ffcs_width=1.0
       if (aux_x.ne.-100.0) ffcs_width=aux_x
    endif
    
    ! initial amplitude perturbation
    amplitude = 1.0
    if (aux_amp.ne.-100.0) amplitude=aux_amp
    relative_pert = 1e-10
    IF ((aux_y .NE. -100.0) .AND. (aux_y .NE. 0.0)) relative_pert = aux_y

    CALL initialize_fourier_x_1d

#ifdef COMBI
    dv = 2.0 * lv / (nv0 - 1.0)
#else
    dv = 2.0 * lv / (nv0)
#endif

    nx0_inv = 1.0 / nx0
    pi_n15 = pi**(-1.5)

    IF (.NOT. vmult_Maxw) THEN ! cosine or Gaussian
       ! parallel electron bulk velocity in units of v_Tj
       vel_shift = 0.5 ! default value
       IF ((aux_z .NE. -100.0) .AND. (aux_z .NE. 0.0)) vel_shift = aux_z
       IF (vel_shift .GT. lv / 2) THEN
          vel_shift = 0.5 * lv
          IF (write_pe) WRITE(*,"(A)") "  current bulk speed limited to Lv/2"
       END IF
       ! from v_Tj shift, calculate shift in terms of the grid
       vel_shift_l = vel_shift / lv * 0.5 * nv0

       IF (write_pe) THEN
          if (.not.force_free) then
             WRITE(*,"(A)") "  current sheets (counter-flowing species)"
             IF (x_v_Gaussian) THEN
                WRITE(*,"(A,L3,A,L3,A,L3,A,ES10.3,A,ES10.3)") &
                     "  vmult_Maxw=", vmult_Maxw, ", x_v_Gaussian=", x_v_Gaussian, &
                     ", nonrandom_phases=", nonrandom_phases, &
                     ", relative_pert=", relative_pert, ", vel_shift=", vel_shift
             ELSE
                WRITE(*,"(A,L3,A,L3,A,L3,A,ES10.3,A,ES10.3,A,I3)") &
                     "  vmult_Maxw=", vmult_Maxw, ", x_v_Gaussian=", x_v_Gaussian, &
                     ", nonrandom_phases=", nonrandom_phases, &
                     ", relative_pert=", relative_pert, ", vel_shift=", vel_shift, &
                     ", x_mode=", x_mode
             END IF
          else
             WRITE(*,"(A)") "  force-free current sheets (only electron flow)"
             WRITE(*,"(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)") &
                  "  ffcs_width=", ffcs_width, ", amplitude=", amplitude,&
                  ", relative_pert=", relative_pert, ", vel_shift=", vel_shift
          endif
       END IF

       
       ! prepare random phases for x modes
       ! note: yields different results on different machines!
       IF (.NOT. nonrandom_phases) THEN
          CALL random_seed(SIZE=rseed_size)
          ALLOCATE(rseed(rseed_size))
          rseed(:) = 1
          CALL random_seed(PUT=rseed)
          DEALLOCATE(rseed)
          
          !default for standard current sheet: identical phases for each ky
          if(.not.force_free) then
             CALL random_number(phases)
             phases(:) = real(EXP(-2.0*pi*CMPLX(0.0,1.0,kind(phases))*phases(:)))
          endif
       END IF
       
       ! fill spatial modes
       DO j = lj1, lj2
          !default for force_free current sheet: different phases for each ky
          if (force_free) then
             CALL random_number(phases)
             phases(:) = real(EXP(-2.0*pi*CMPLX(0.0,1.0,kind(phases))*phases(:)))
          endif
          DO i = li1, li2
             IF ((i .EQ. 0) .AND. (j .EQ. 0)) THEN
                dens(i,j,:) = amplitude
             ELSE
                IF (.NOT. nonrandom_phases) THEN
                   if (force_free) then
                      !force-free cs is initialized with noise modeling the Apar spectrum of the PIC code ACRONYM
                      dens(i,j,:)=phases(i)*amplitude*relative_pert*(min(abs(ki(i)),abs(ki(12)))*&
                           min(abs(kj(j)),kymin*12))**0.8
                      !an additional perturbation, similar to that of the GEM challenge is added
                      if (i.eq.0.and.j.eq.1) dens(i,j,:)=amplitude*1e-3
                   else
                      !standard current sheet: uniform noise
                      dens(i,j,:) = phases(i) * relative_pert * amplitude
                   endif
                ELSE
                   ifrac = REAL(MODULO(i,2)) - 0.5
                   dens(i,j,:) = CMPLX(2.0*ifrac,1.0-2.0*ifrac) * &
                        relative_pert * amplitude
                END IF
             END IF
          END DO
          
          ! ensure complex conjugate condition for ky=0
          IF (j .EQ. 0) THEN
             k = hkx
             DO i = lkx, li2
                dens(i,j,:) = CONJG(dens(k,j,:))
                k = k - 1
             END DO
          END IF
          
          DO k = lk1, lk2
             CALL to_real_x_1d(dens(:,j,k),temp_x)
             dens(:,j,k) = temp_x
          END DO
       END DO
       
       ! for x_v_Gaussian, ensure that both Gaussians are shifted by a
       ! constant value to reach zero at the box center/edge
       IF (x_v_Gaussian) Gauss_offset_corr = EXP(-(10.0*(-nx0*0.25)*nx0_inv)**2)
       
       DO n = ln1, ln2
          DO m = lm1, lm2
             DO l = ll1, ll2
                DO k = lk1, lk2
                   if (.not.force_free) then
                      DO i = li1, li2
                         ! shifted Maxwellians (different bulk velocity sign for ions, electrons)
                         IF (spec(n)%charge .EQ. -1) THEN
                            IF (x_v_Gaussian) THEN
                               IF (i .LT.nx0 / 2) THEN
                                  lmod = l - vel_shift_l * &
                                       (EXP(-(10.0*(i-nx0*0.25)*nx0_inv)**2) - Gauss_offset_corr)
                               ELSE
                                  lmod = l + vel_shift_l * &
                                       (EXP(-(10.0*(i-nx0*0.75)*nx0_inv)**2) - Gauss_offset_corr)
                               END IF
                            ELSE
                               lmod = l - vel_shift_l * COS(-2.0*pi*i*nx0_inv*x_mode)
                            END IF
                         ELSE
                            IF (x_v_Gaussian) THEN
                               IF (i .LT. nx0 / 2) THEN
                                  lmod = l + vel_shift_l * &
                                       (EXP(-(10.0*(i-nx0*0.25)*nx0_inv)**2) - Gauss_offset_corr)
                               ELSE
                                  lmod = l - vel_shift_l * &
                                       (EXP(-(10.0*(i-nx0*0.75)*nx0_inv)**2) - Gauss_offset_corr)
                               END IF
                            ELSE
                               lmod = l + vel_shift_l * COS(-2.0*pi*i*nx0_inv*x_mode)
                            END IF
                         END IF
                         ! construct fractional index (lmod) shifted Maxwellian
#ifdef COMBI
                         fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv+SHIFT)**2))
#else
                         fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv)**2))
#endif
                         g_1(i,:,k,l,m,n) = dens(i,:,k) * fm_lmod
                      END DO
                      
                      DO j = lj1, lj2
                         CALL to_fourier_x_1d(g_1(:,j,k,l,m,n),temp_x)
                         g_1(:,j,k,l,m,n) = temp_x
                      END DO
                   else
                      !now force-free current sheet
                      !we need to discern different ky modes
                      do j=lj1,lj2
                         if (j.eq.0) then
                            DO i = li1, li2
                               ! shifted Maxwellians (bulk velocity only for electrons)
                               IF (spec(n)%charge .EQ. -1) THEN
                                  lmod = l - vel_shift_l * xfunc(i,ffcs_width) 
                               ELSE
                                  !stationary ions
                                  lmod=l
                               END IF
                               
                               ! construct fractional index (lmod) shifted Maxwellian
#ifdef COMBI
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv+SHIFT)**2))
#else
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv)**2))
#endif
                               g_1(i,j,k,l,m,n) = dens(i,j,k) * fm_lmod
                            END DO
                            !GEM perturbation
                         else if (j.eq.1) then
                            DO i = li1, li2
                               ! shifted Maxwellians (different bulk velocity sign for ions, electrons)
                               IF (spec(n)%charge .EQ. -1) THEN
                                  lmod = l - vel_shift_l * cos(2.0*pi*i/nx0-pi/2.) 
                               ELSE
                                  !stationary ions
                                  lmod=l
                               END IF
                               
                               ! construct fractional index (lmod) shifted Maxwellian
#ifdef COMBI
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv+SHIFT)**2))
#else
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv)**2))
#endif
                               g_1(i,j,k,l,m,n) = dens(i,j,k) * fm_lmod
                            END DO
                         else
                            do i=li1,li2
                               if (spec(n)%charge.eq.-1) then
                                  lmod=l-vel_shift_l 
                               else
                                  lmod=l
                               endif
#ifdef COMBI
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv+SHIFT)**2))
#else
                               fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv)**2))
#endif
                               g_1(i,j,k,l,m,n) = dens(i,j,k) * fm_lmod
                            enddo
                         endif
                         CALL to_fourier_x_1d(g_1(:,j,k,l,m,n),temp_x)
                         g_1(:,j,k,l,m,n) = temp_x
                         !we zero out the (0,0) mode for the force-free current to match AstroGK more closely
                         if (p_has_0_mode) g_1(0,0,:,:,:,:)=0.
                      enddo
                   endif
                END DO
             END DO
          END DO
       END DO
       
    ELSE ! alternate implementation: v_par multiplication
       
       IF (write_pe) THEN
          WRITE(*,"(A)") "  current sheets (counter-flowing species)"
          WRITE(*,"(A,L3,A,ES10.3)") &
               "  vmult_Maxw=", vmult_Maxw, ", relative_pert=", relative_pert
       END IF
       
       DO n = ln1, ln2
          DO m = lm1, lm2
             DO l = ll1, ll2
                DO k = lk1, lk2
                   fm_lmod = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+l*dv)**2))
                   DO i = li1, li2
                      IF (spec(n)%charge .EQ. -1) THEN
#ifdef COMBI
                         g_1(i,:,k,l,m,n) = 30.0 * amplitude * relative_pert * &
                              fm_lmod * ((-lv+l*dv+SHIFT) * SIN(2.0*pi*i*nx0_inv))! + 1.0)
                         g_1(i,lj1,k,l,m,n) = 30.0 * amplitude * &
                              fm_lmod * ((-lv+l*dv+SHIFT) * SIN(2.0*pi*i*nx0_inv))! + 1.0)
#else
                         g_1(i,:,k,l,m,n) = 30.0 * amplitude * relative_pert * &
                              fm_lmod * ((-lv+l*dv) * SIN(2.0*pi*i*nx0_inv))! + 1.0)
                         g_1(i,lj1,k,l,m,n) = 30.0 * amplitude * &
                              fm_lmod * ((-lv+l*dv) * SIN(2.0*pi*i*nx0_inv))! + 1.0)
#endif
                      ELSE
                         g_1(i,:,k,l,m,n) = 30.0 * amplitude * relative_pert * fm_lmod
                         g_1(i,lj1,k,l,m,n) = 30.0 * amplitude * fm_lmod
                      END IF
                   END DO
                   
                   DO j = lj1, lj2
                      CALL to_fourier_x_1d(g_1(:,j,k,l,m,n),temp_x)
                      g_1(:,j,k,l,m,n) = temp_x
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF
    
    CALL finalize_fourier_x_1d
    
    ! the variable g_1 above is actually the f distribution; thus, need to
    ! compute the fields from f and evaluate the actual g_1 distribution
    CALL calc_aux_fields_from_f(g_1,emfields,f_)
    g_1 = f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    
  contains
    real function xfunc(i,width)
      integer, intent(in):: i
      real:: x, width
      
      x=(deli*i-lx/4.)/width
      xfunc=1./cosh(x)**2 
      x=(deli*i+lx/4.)/width
      xfunc=xfunc-1./cosh(x)**2 
      x=(deli*i-3.*lx/4.)/width
      xfunc=xfunc-1./cosh(x)**2 
      x=(deli*i-5.*lx/4.)/width
      xfunc=xfunc+1./cosh(x)**2 


    end function xfunc
    
  END SUBROUTINE current_sheets

  !*************************************************************************!

  !> plasmoids initialization
  !! note: only implemented for local code!
  SUBROUTINE plasmoids(aux_x,aux_y,aux_z)

    IMPLICIT NONE

    REAL, INTENT(IN) :: aux_x, aux_y, aux_z
    REAL, DIMENSION(:,:), ALLOCATABLE :: temp_sxsy, expDr2, expDr2_2, &
      expDr2_ring, expDr2_ring_2
    COMPLEX, DIMENSION(li1:li2,lj1:lj2) :: temp_kxky
    REAL :: pl_sep, pl_width, pl_width_2, vel_shift, vel_shift_2, &
      vel_shift_l, vel_shift_2_l, neg_ring_amp, neg_ring_amp_2, &
      pi_n15, lmod, dv
    INTEGER :: n_plasmoids, i, j, k, l, m, n, nlj_lo, nlj_hi
    LOGICAL :: use_extra_pars

    ! for additional options/parameters, see below
    use_extra_pars = .true. ! .false.

    ! check: code needs to be run in local mode
    IF (.NOT. xy_local) THEN
      IF (write_pe) WRITE(*,"(A)") &
        "ERROR: plasmoid initial condition requires local code"
      STOP
    END IF

    ! number of plasmoids (1 or 2) and
    ! initial separation in units of L_x (for 2 plasmoids)
    n_plasmoids = 1
    pl_sep = 0.0
    IF ((aux_x .NE. -100.0) .AND. (aux_x .NE. 0.0)) THEN
      IF (aux_x .GE. 2.0) THEN
        n_plasmoids = 2
        pl_sep = aux_x - REAL(n_plasmoids)
        IF (pl_sep .LE. 0.0) pl_sep = 0.5
      ELSE
        n_plasmoids = 1
      END IF
    END IF

    ! plasmoid width in units of L_x
    pl_width = 0.1
    IF ((aux_y .NE. -100.0) .AND. (aux_y .NE. 0.0)) pl_width = aux_y

    ! plasmoid amplitude (electron v_par shift in units of v_Te)
    vel_shift = 0.5
    IF ((aux_z .NE. -100.0) .AND. (aux_z .NE. 0.0)) vel_shift = aux_z
    ! from v_Te shift, calculate shift in terms of the grid
    vel_shift_l = vel_shift / lv * 0.5 * nv0

    IF (write_pe) THEN
      WRITE(*,"(A)") "  plasmoids initialization"
      WRITE(*,"(A,I1,A,ES10.3,A,ES10.3,A,ES10.3)") &
        "  n_plasmoids=", n_plasmoids, ", pl_sep=", pl_sep, &
        ", pl_width=", pl_width, ", vel_shift=", vel_shift
    END IF

    IF (use_extra_pars) THEN
      pl_width_2 = pl_width ! 0.2 ! width of plasmoid 2

      ! parallel velocity shift of plasmoid 2
      vel_shift_2 = vel_shift ! 0.25
      vel_shift_2_l = vel_shift_2 / lv * 0.5 * nv0

      ! relative amplitude of negative current ring for plasmoids 1 and 2
      neg_ring_amp = 0.0 ! 0.5
      neg_ring_amp_2 = 0.0 ! 0.5

      IF (write_pe) THEN
        WRITE(*,"(A)") "  with additional parameters"
        WRITE(*,"(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)") &
          "  pl_width_2", pl_width_2, ", vel_shift_2=", vel_shift_2, &
          ", vel_shift_2=", vel_shift_2, &
          ", neg_ring_amp=", neg_ring_amp, ", neg_ring_amp_2=", neg_ring_amp_2
      END IF
    ELSE
      pl_width_2 = pl_width
      vel_shift_2 = vel_shift
      vel_shift_2_l = vel_shift_2 / lv * 0.5 * nv0

      neg_ring_amp = 0.143 ! roughly zero total current
      neg_ring_amp_2 = neg_ring_amp
    END IF

    CALL initialize_fourier(li0da,ly0da)

    nlj_lo = my_pey * li0da / n_procs_y
    nlj_hi = (my_pey + 1) * li0da / n_procs_y - 1

    ALLOCATE(temp_sxsy(0:ly0da-1,nlj_lo:nlj_hi),expDr2(0:ly0da-1,nlj_lo:nlj_hi),&
      expDr2_2(0:ly0da-1,nlj_lo:nlj_hi),expDr2_ring(0:ly0da-1,nlj_lo:nlj_hi),&
      expDr2_ring_2(0:ly0da-1,nlj_lo:nlj_hi))

    ! exp(-Delta r^2) for 1 or 2 plasmoids (nl_to_fourier_xy expects reverse index order!)
    DO j = nlj_lo, nlj_hi
      DO i = 0, ly0da - 1
        expDr2(i,j) = EXP(-((((j / REAL(li0da) - 0.5) + 0.5 * pl_sep) * lx)**2 + &
          ((i / REAL(ly0da) - 0.5) * ly)**2)/pl_width)
        expDr2_ring(i,j) = EXP(-((SQRT((((j / REAL(li0da) - 0.5) + 0.5 * pl_sep) * lx)**2 + &
          ((i / REAL(ly0da) - 0.5) * ly)**2) - pl_width * lx)**2)/pl_width)

        IF (n_plasmoids .EQ. 2) THEN
          expDr2_2(i,j) = EXP(-((((j / REAL(li0da) - 0.5) - 0.5 * pl_sep) * lx)**2 + &
            ((i / REAL(ly0da) - 0.5) * ly)**2)/pl_width_2)
          expDr2_ring_2(i,j) = EXP(-((SQRT((((j / REAL(li0da) - 0.5) - 0.5 * pl_sep) * lx)**2 + &
            ((i / REAL(ly0da) - 0.5) * ly)**2) - pl_width_2 * lx)**2)/pl_width_2)
        END IF
      END DO
    END DO

#ifdef COMBI
    dv = 2.0 * lv / (nv0 )
#else
    dv = 2.0 * lv / (nv0 - 1.0)
#endif
    pi_n15 = pi**(-1.5)

    DO n = ln1, ln2
      DO m = lm1, lm2
        DO l = ll1, ll2
          DO k = lk1, lk2
            DO j = nlj_lo, nlj_hi
              DO i = 0, ly0da - 1
                ! shifted Maxwellians (only for electrons)
                IF (spec(n)%charge .EQ. -1) THEN
                  lmod = l - vel_shift_l * (expDr2(i,j) - neg_ring_amp * expDr2_ring(i,j))
                  IF (n_plasmoids .EQ. 2) lmod = lmod - &
                    vel_shift_2_l * (expDr2_2(i,j) - neg_ring_amp_2 * expDr2_ring_2(i,j))
                ELSE
                  lmod = l
                END IF

#ifdef COMBI
                temp_sxsy(i,j) = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv+SHIFT)**2))
#else
                temp_sxsy(i,j) = pi_n15 * EXP(-(mu(m)*geom%Bfield(pi1,pj1,k)+(-lv+lmod*dv)**2))
#endif
              END DO
            END DO

            CALL nl_to_fourier_xy(temp_sxsy,temp_kxky)

            g_1(li1:li2,lj1:lj2,k,l,m,n) = temp_kxky(li1:li2,lj1:lj2)
          END DO
        END DO
      END DO
    END DO

    CALL finalize_fourier

    DEALLOCATE(temp_sxsy,expDr2,expDr2_2,expDr2_ring,expDr2_ring_2)

    ! the variable g_1 above is actually the f distribution; thus, need to
    ! compute the fields from f and evaluate the actual g_1 distribution
    CALL calc_aux_fields_from_f(g_1,emfields,f_)
    g_1 = f_(:,:,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)

  END SUBROUTINE plasmoids

  !*************************************************************************!

  !> zero initialization
  !! initializes a zeroed out distribution
  !! intended use: antenna current
  SUBROUTINE antenna_zero
    
    IF (write_pe) WRITE(*,"(A)") "  zero init. cond."

    g_1(:,:,:,:,:,:) = 0.0
                      
  End Subroutine antenna_zero

  !*************************************************************************!
  !*************************************************************************!

  !> Add the velocity space distribution
  !! Typically, we use a Maxwellian for the initialization
  !! However, minor modification, e.g. division by J0 can be switched on 
  !! additionally
  subroutine add_fm(dens_in,n,J0_pow_in)
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2),INTENT(IN) :: dens_in
    INTEGER, INTENT(IN) :: n
    real, optional, intent(IN) :: J0_pow_in
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lk1:lk2) :: my_dens
    COMPLEX, DIMENSION(li1:li2):: temp_x
    complex,dimension(lk1:lk2):: A_ant
    INTEGER :: i,j,k,l,m, pii,pji,pni, hki, p
    real :: J0_pow
    real, dimension(li1:li2,lj1:lj2):: vshift

    IF (present(J0_pow_in)) THEN
       J0_pow = J0_pow_in
    ELSE
       J0_pow = 0
    ENDIF
    
    pii=pi1; pji=pj1; pni=pn1 

    my_dens = dens_in
    hki = (ni0-1)/2
    
    !set Nyquist mode to zero
    IF (cond_local_in_i) THEN
       IF ((yx_order).AND.(xy_local)) THEN
          IF ((evenx.eq.1).AND.((hkx+1.GE.lg1).AND.(hkx+1.LE.lg2))) &
               &my_dens(:,hkx+1,:)=(0.,0.)
       ELSE
          IF ((even_ni0.eq.1).AND.(hki+1.GE.lg1).AND.(hki+1.LE.lg2)) &
               my_dens(hki+1,:,:)=(0.,0.)
       ENDIF
    END IF
    
    IF (cond_local_in_i.AND.(.NOT.xy_local)) THEN
       CALL initialize_fourier_x_1d
       DO k = lk1, lk2
          DO j=lj1,lj2
             CALL to_real_x_1d(my_dens(:,j,k), temp_x)
             my_dens(li1:li2,j,k) = temp_x(li1:li2)
          END DO
       END DO
       CALL finalize_fourier_x_1d
    ELSEIF ((.NOT.cond_local_in_i).AND.(xy_local)) THEN
       CALL initialize_fourier_x_1d
       DO k = lk1, lk2
          DO j=lj1,lj2
             CALL to_fourier_x_1d(my_dens(:,j,k),temp_x)
             my_dens(:,j,k) = temp_x
          ENDDO
       END DO
       CALL finalize_fourier_x_1d
    ENDIF

    if (pn0.gt.1) pni=n
    if (init_vpar_fluct.gt.0.) then
       call sgrnd(1)
       do j=lj1,lj2
          do i=li1,li2
             vshift(i,j)=init_vpar_fluct*grnd()
          enddo
       enddo
    endif

    Do m=lm1,lm2
       Do l= ll1,ll2
          Do k = lk1, lk2          
             Do j = lj1,lj2
                if (pj0.gt.1) pji=j
                Do i = li1, li2
                   if (pmi0.gt.1) pii=i
                   if (init_vpar_fluct.gt.0.) then
                      g_1(i,j,k,l,m,n) = my_dens(i,j,k)*spec(n)%dens_prof(pii)*&
                           (pi*spec(n)%temp_prof(pii))**(-1.5)*&
                           exp(-(mu(m)*geom%Bfield(pii,pji,k)+(vp(l)+vshift(i,j))**2)/&
                           spec(n)%temp_prof(pii))
                   else
                      g_1(i,j,k,l,m,n) = my_dens(i,j,k)*fm(pii,pji,k,l,m,pni)           
                   endif
                   if (J0_pow.ne.0) then
                      if ((J0_pow.lt.0).and.(abs(get_jfac(i,j,k,m,n)).lt.epsilon(abs(get_jfac(i,j,k,m,n))))) then
                         g_1(i,j,k,l,m,n) = 0.0
                      else
                         g_1(i,j,k,l,m,n) = g_1(i,j,k,l,m,n)*get_jfac(i,j,k,m,n)**J0_pow
                      endif
                   endif
                End Do
             End Do
             if (precomp_nc) then
                if(x_local.and.p_has_00_mode) g_1(0,0,k,l,m,n) = g1_nc(pi1,pj1,k,l,m,n)
                if((.not.x_local).and.p_has_0_mode) g_1(:,0,k,l,m,n) = g1_nc(:,pj1,k,l,m,n)
             end if
          End Do
          if (antenna_type.eq.2) then
             do p=1,n_antennae
                call get_Apar_antenna(p,A_ant,i,j)
                do k=lk1,lk2
                   g_1(i,j,k,l,m,n)=g_1(i,j,k,l,m,n)+spec(n)%charge*sqrt(2./spec(n)%temp/spec(n)%mass)*&
                        vp(l)*fm(pi1,pj1,k,l,m,pni)*A_ant(k)*get_jfac(i,j,k,m,n)
                   if (j+ky0_ind.eq.0.and.p_has_0_mode) &
                        g_1(nx0-i,j,k,l,m,n)=g_1(nx0-i,j,k,l,m,n)+spec(n)%charge*sqrt(2./spec(n)%temp/spec(n)%mass)*&
                        vp(l)*fm(pi1,pj1,k,l,m,pni)*conjg(A_ant(k))*get_jfac(nx0-i,j,k,m,n)
                enddo
             enddo
          endif
       End Do
    End Do

  End Subroutine add_fm
  
  !=========================================================================================
  ! subroutines necessary for block structured grids

  ! > set g1 to zero outside block structured grid
  subroutine nullify_g1_outside_adptv
    integer :: i, j, k, l, m, n

    do n = ln1, ln2
       do m = lm1, lm2
          do l = ll1, ll2
             do k = lk1, lk2
                do j = lj1, lj2
                   do i = li1, li1_vwadp(l,m)-1
                      g_1(i,j,k,l,m,n) = 0.
                   enddo
                   do i = li2_vwadp(l,m)+1, li2
                      g_1(i,j,k,l,m,n) = 0.
                   enddo
                enddo
             enddo
          endDo
       endDo
    enddo

  end subroutine nullify_g1_outside_adptv
  
end module initcond
