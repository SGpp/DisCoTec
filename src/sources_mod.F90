!> Module containing sources and sinks subroutines
#include "redef.h"
#include "intrinsic_sizes.h"
Module sources_mod
  use communications
  use discretization
  use discretization_adptv_module
  use coordinates
  use prefactors
  use geometry
  use par_mod, only: pi, in_timestep_comp, time, g_1
  use par_in, only: spec, n0_global, n_pol, diagdir, &
       &file_extension, rad_bc_type
  use vel_space, only : mat_00, mat_10, fm
  use par_other, only : print_ini_msg, p_has_0_mode,dt
  use file_io, only: get_unit_nr
  use blockindex
  use profiles_mod, only: drive_buffer, ldrive_buffer_size, udrive_buffer_size, &
       &ldrive_coef, udrive_coef, ldrive_pow, udrive_pow

  Implicit None

  PUBLIC :: initialize_f1_sources,  finalize_f1_sources, &
       &compute_f_vabs_allspec, add_f1_sources,&
       &add_f1_sources_block, mem_est_f1_sources,&
       &compute_f_vabs_n, ck_heat, ck_part, intsource_heat_coeff, &
       &intsource_part_coeff, intsource_time,&
       &coeff_g1_en_source, g1_init, delivered_g1ens_msg, &
       &heatsrc_on, partsrc_on, simtime_start,&
       !relative low pass filter sizes for heat/part. krook op.
       &ck_heat_smooth, ck_part_smooth, ck_filter_type,&
       !buffer zone definitions for krook operator
       &l_buffer_size,lcoef_krook,buffer_on_ky0,&
       &u_buffer_size, ucoef_krook, lpow_krook, upow_krook,&
       &initialize_krookbuffer_operator, add_krookbuffer, &
       &finalize_krookbuffer_operator,add_kBuffer_explicit,&
       !buffer zone definitions for heat/particle krook operator
       &lck_buffer_size, lckcoef_krook,&
       &uck_buffer_size, uckcoef_krook, lckpow_krook, uckpow_krook,&
       !heat pulses
       &puls_start, puls_amp, puls_on_duration, puls_off_duration,&
       &puls_func, add_ck_heat, add_ck_part_1, add_ck_part_2, &
       &initialize_f0_heat_source,f0_heat_src_active,&
       &finalize_f0_heat_source, check_sources, &
       &initialize_g1_source, finalize_g1_source, explicit_buffer, with_sources,&
       &psource_type
  
  ! --- routines for the adaptive grids
  public:: add_kBuffer_explicit_adptv, add_f1_sources_block_adptv, add_krookBuffer_adptv

  PRIVATE
  logical:: with_sources=.true.

  Integer :: init_status = 0, init_status_krookbuffer_operator = 0, &
       &init_status_heat_source = 0
  logical :: f0_heat_src_active = .false., explicit_buffer=.false.

  ! Krook buffer zone operator
  Real, Dimension(:), Allocatable :: nu_K
  real :: simtime_start=-1.0

  ! Variables for the Krook buffer zone operator l:lower; u:upper
  ! l/u_buffer_size between 0 and 1
  Real :: l_buffer_size=0.0, lcoef_krook=0.0, & 
       u_buffer_size=0.0, ucoef_krook=0.0
  Integer :: lpow_krook=4,upow_krook=4

  ! Input variable for the Krook-Type heat and particle sources
  ! to set base level
  Real :: ck_heat = 0.0, ck_part = 0.0
  ! Relative (wrt lx) low pass filter sizes:
  Real :: ck_heat_smooth = 0.0, ck_part_smooth = 0.0
  integer:: ck_filter_type=0
  ! Additional profiles for these coefficients to, e.g., introduce
  ! buffer layers
  Real, dimension(:), Allocatable :: ck_heat_prof, ck_part_prof
  Real :: lck_buffer_size=0.0, lckcoef_krook=-1.0, & 
       uck_buffer_size=0.0, uckcoef_krook=-1.0
  Integer :: lckpow_krook=4,uckpow_krook=4

  ! Variables for g1 energy source
  REAL :: coeff_g1_en_source = 0.0
  logical :: allocated_g1_en_source = .false., delivered_g1ens_msg = .false., &
       heatsrc_on=.false., partsrc_on=.false., buffer_on_ky0=.false.

  integer,dimension(:),allocatable:: map_to_f
  integer:: psource_type=2 !type of particle source: 1=Told, 2=McMillan
  real:: intsource_heat_coeff=0.0, intsource_part_coeff=0.0     !coefficients to integrating sources
  real:: intsource_time=1000.0     !characteristic decay time of integrating source, in normalized units
  real:: norm_time=1000.0 !normalization time for integrator; is set to the minimum
  !of integration time and time since (re)start of the simulation

  ! Flux surface average fm, and its integral over vel space 
  ! for the krook type heat src
  Real, Dimension(:,:,:,:), Allocatable :: k_heat_prefac, k_part_prefac


  Complex, Dimension(:,:,:,:), Allocatable :: f_vabs         !< <f(x,|vpar|,mu)>_FS 
                                                             !!   = <0.5*(f(x,vpar,mu) + f(x,- vpar,mu))>_FS
  Complex, Dimension(:,:), Allocatable :: int_vel_f_vabs     !< int <f(x,|vpar|,mu)>_FS dv
  Complex, Dimension(:,:), Allocatable :: int_vel_vsq_f_vabs !< int v**2 <f(x,|vpar|,mu)>_FS dv
  Real, Dimension(:,:,:,:), Allocatable :: fm_fsa, heat_prof, part_prof
  real, dimension(:,:), allocatable :: int_vel_fm_fsa
  Real, Dimension(:,:), Allocatable :: int_vel_vsq_fm_fsa
  complex, dimension(:), allocatable:: sum_vint_f_vabs !< species sum of int_vel_f_vabs
  complex, dimension(:), allocatable:: sum_q_vint_f_vabs !< species sum of int_vel_f_vabs*charge
  Real, Dimension(:,:), Allocatable :: coefficient !dynamically adapted coefficient
  !to eliminate heat contribution from particle source
!  Complex, Dimension(:,:,:,:), Allocatable :: krook_heat_cont !< krook heat source contribution
  COMPLEX, DIMENSION(:,:,:,:,:), ALLOCATABLE :: g1_init !< initial g1 for g1_en_source

  Logical :: write_pe

  !Variables for heat pulses
  Real :: puls_start=100000.    !starting time for heat pulses
  Real :: puls_amp=1.0          !relative heat pulse amplitude
  Real :: puls_on_duration=50.  !time duration of heat pulse
  Real :: puls_off_duration=50. !time between two subsequent pulses


  INTERFACE add_f1_sources
     MODULE PROCEDURE add_f1_sources_klmn, add_f1_sources_block
  END INTERFACE  
  
Contains

  subroutine check_sources
    integer:: n

    f0_heat_src_active=.false.
    do n=0, n_spec-1
       f0_heat_src_active=f0_heat_src_active.or.(ANY(spec(n)%src_amp.GT.0))
    enddo
    if (f0_heat_src_active.and.n0_global.eq.0.and.kymin.eq.0)  then 
       if (mype==0) then
          write (*,'(A)') "STOP f0_heat_source does not make sense with n0_global=0 and kymin=0"
          write (*,'(A)') "     (cannot define ly for flux surface averaging)!!"
          if(nky0==1) write (*,'(A)') "     use n0_global=1 if only ky=0 mode is computed"
       endif
       stop
    endif
    !Heaviside buffer zones on ky=0 currently require explicit buffers
    if (buffer_on_ky0) explicit_buffer=.true.

    if (x_local) then
       if (ck_heat.ne.0.0) then
          if (mype.eq.0) print*, 'Setting ck_heat=0.0 in x local simulations'
          ck_heat = 0.0
       endif
       if (ck_part.ne.0.0) then
          if (mype.eq.0) print*, 'Setting ck_part=0.0 in x local simulations'
          ck_part = 0.0
       endif
    endif

    if (ck_part.gt.0.0.or.ck_heat.gt.0.0.or.intsource_heat_coeff.gt.0.0.or.intsource_part_coeff.gt.0.0) &
         heatsrc_on=.true.
    if (ck_part.gt.0.0.or.intsource_part_coeff.gt.0.0) partsrc_on=.true.

    if (drive_buffer) then
       ldrive_buffer_size=l_buffer_size
       udrive_buffer_size=u_buffer_size
       ldrive_coef=lcoef_krook
       udrive_coef=ucoef_krook
       ldrive_pow=lpow_krook
       udrive_pow=upow_krook
    endif

    if ((coeff_g1_en_source.lt.epsilon(coeff_g1_en_source)).and.&
         &.not.heatsrc_on.and.&
         &.not.partsrc_on) then 
       with_sources=.false.
    else
       with_sources=.true.
    endif

  end subroutine check_sources

  Real function mem_est_f1_sources(mem_req_in)
    real:: mem_req_in, mem_loc

    !krook buffer zone operator
    mem_loc = pi0*SIZE_OF_REAL_MB
    
    if (ck_heat.GT.0.or.ck_part.GT.0) then
       !krook type heat source: fm_fsa, int_vel_fm_sfa, f_vabs, int_vel_f_vabs, k_heat_prefac
       mem_loc = mem_loc + (3.*li0*ll0*lm0*ln0+li0*ln0)*SIZE_OF_COMPLEX_MB
       !ck_heat_prof, ck_part_prof
       mem_loc = mem_loc + 2.*nx0*SIZE_OF_REAL_MB
    endif

    !buffer_f
    mem_loc = mem_loc + li0*ll0*SIZE_OF_COMPLEX_MB

    ! g1 energy source
    IF (coeff_g1_en_source .GT. 0.0) mem_loc = mem_loc + &
      (li0*lk0*ll0*lm0*ln0) * SIZE_OF_COMPLEX_MB

    mem_est_f1_sources=mem_req_in + mem_loc
  end function mem_est_f1_sources

  !>Initialize krook type heat and particle sources/sinks
  Subroutine initialize_f1_sources

    if (init_status==1) then
       heat_prof=0.
       part_prof=0.
       return
    endif 

    write_pe = ((mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec))&
         &.AND.(print_ini_msg))
    
    ! Damping Krook operator
!    call init_krookbuffer_operator !is now called from time_scheme

    ! Heat/Particle source krook using f1
    if (heatsrc_on.or.partsrc_on) CALL init_k_heat_src

    CALL initialize_g1_source

    init_status = 1
  end Subroutine initialize_f1_sources

  !>Initialize Krook-type g1 source
  SUBROUTINE initialize_g1_source
    IF (coeff_g1_en_source .GT. 0.0) THEN
      IF (.NOT. allocated_g1_en_source) THEN
        ALLOCATE(g1_init(li1:li2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
        allocated_g1_en_source = .true.
      END IF
    END IF
  END SUBROUTINE initialize_g1_source

  !>Initialize source mod
  Subroutine finalize_f1_sources

    if (heatsrc_on.or.partsrc_on) then
       deallocate(ck_heat_prof, ck_part_prof)
       deallocate(f_vabs, coefficient)
       deallocate(fm_fsa, int_vel_f_vabs,k_heat_prefac,k_part_prefac,int_vel_vsq_f_vabs, int_vel_fm_fsa, &
            int_vel_vsq_fm_fsa,sum_vint_f_vabs,sum_q_vint_f_vabs,heat_prof,part_prof)
    endif

    CALL finalize_g1_source

    init_status = 0
  end Subroutine finalize_f1_sources

  !>Finalize Krook-type g1 source
  SUBROUTINE finalize_g1_source
    IF (coeff_g1_en_source .GT. 0.0) THEN
      IF (allocated_g1_en_source) THEN
        DEALLOCATE(g1_init)
        allocated_g1_en_source = .false.
      END IF
    END IF
  END SUBROUTINE finalize_g1_source

  !--------------------------------------------------------------------
  ! Routines used for the Krook-type Heat src

  !>Initializes Krook type heat src, e.g., for Vlasov eq.
  Subroutine init_k_heat_src
    INTEGER:: i,k,l,m,n, nxbuf, ibuf
    real,dimension(li1:li2):: sum_vint_fm_fsa
    
    PERFON('in_krook')
    ALLOCATE(ck_heat_prof(pi1:pi2))
    ALLOCATE(ck_part_prof(pi1:pi2))
    ALLOCATE(fm_fsa(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2))
    ALLOCATE(f_vabs(li1:li2, ll1:ll2, lm1:lm2, ln1:ln2))
    fm_fsa=0.0
    ALLOCATE(int_vel_fm_fsa(li1:li2,ln1:ln2))
    ALLOCATE(coefficient(li1:li2,ln1:ln2))
    ALLOCATE(int_vel_vsq_fm_fsa(li1:li2,ln1:ln2))
    ALLOCATE(k_heat_prefac(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2))
    ALLOCATE(k_part_prefac(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2))
    ALLOCATE(int_vel_f_vabs(li1:li2,ln1:ln2))
    ALLOCATE(int_vel_vsq_f_vabs(li1:li2,ln1:ln2))
    allocate(sum_q_vint_f_vabs(li1:li2))
    allocate(sum_vint_f_vabs(li1:li2))
    allocate(heat_prof(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2),part_prof(li1:li2,ll1:ll2,lm1:lm2,ln1:ln2))
    int_vel_fm_fsa = 0.0
    int_vel_vsq_fm_fsa = 0.0
    heat_prof=0.0
    part_prof=0.0
    
    ck_heat_prof = ck_heat
    ck_part_prof = ck_part

    if ((lckcoef_krook.ge.0.0).and.(lck_buffer_size.gt.0)) then
       nxbuf=int(lck_buffer_size*nx0)
       ibuf =pi1gl+nxbuf
       do i=pi1gl, pi1gl+nxbuf
          If (i<pi1 .Or. i>pi2) Cycle
          ck_heat_prof(i)=ck_heat_prof(i)+(lckcoef_krook-ck_heat_prof(i))*&
               &(abs(i-ibuf)/real(nxbuf))**lckpow_krook
          ck_part_prof(i)=ck_part_prof(i)+(lckcoef_krook-ck_part_prof(i))*&
               &(abs(i-ibuf)/real(nxbuf))**lckpow_krook
       end do
    endif
    
    if ((uckcoef_krook.ge.0.0).and.(uck_buffer_size.gt.0)) then
       nxbuf=int(uck_buffer_size*nx0)
       ibuf =pi2gl-nxbuf
       do i=pi2gl-nxbuf, pi2gl
          If (i<pi1 .Or. i>pi2) Cycle
          ck_heat_prof(i)=ck_heat_prof(i)+(uckcoef_krook-ck_heat_prof(i))*&
               &(abs(i-ibuf)/real(nxbuf))**uckpow_krook
          ck_part_prof(i)=ck_part_prof(i)+(uckcoef_krook-ck_part_prof(i))*&
               &(abs(i-ibuf)/real(nxbuf))**uckpow_krook
       end do
    endif

    IF (p_has_0_mode) then
       ! Compute < fm >_yz

       do n=ln1,ln2
          Do m=lm1,lm2
             Do l=ll1,ll2
                DO k=lk1,lk2
                   fm_fsa(li1:li2,l,m,n)=fm_fsa(li1:li2,l,m,n)+ & 
                        fm(li1:li2,lj1,k,l,m,n)*geom%jacobian(li1:li2,lj1,k)
                end Do
                call  my_sum_to_all_real(fm_fsa(li1:li2,l,m,n),li0,mpi_comm_z)
                fm_fsa(li1:li2,l,m,n)=fm_fsa(li1:li2,l,m,n)/(Real(nz0)*geom%avg_jaco_yz(li1:li2))
             end Do
          end Do
       end do
       ! Compute < int < fm >_yz dv >_yz


          do n=ln1,ln2
             Do m=lm1,lm2
                Do l=ll1,ll2
                   do k=lk1,lk2
                      int_vel_fm_fsa(li1:li2,n) = int_vel_fm_fsa(li1:li2,n) + &
                           fm_fsa(li1:li2,l,m,n)*mat_00(li1:li2,lj1,k,l,m)*geom%jacobian(li1:li2,pj1,k)
                   end do
                end Do
             end Do
             call my_sum_to_all_real(int_vel_fm_fsa(li1:li2,n),li0,mpi_comm_vw)
             call my_sum_to_all_real(int_vel_fm_fsa(li1:li2,n),li0,mpi_comm_z)
             int_vel_fm_fsa(li1:li2,n)=int_vel_fm_fsa(li1:li2,n)/(Real(nz0)*geom%avg_jaco_yz(li1:li2))
          end Do

          do n=ln1,ln2
             Do m=lm1,lm2
                Do l=ll1,ll2
                   do k=lk1,lk2
                      int_vel_vsq_fm_fsa(li1:li2,n) = int_vel_vsq_fm_fsa(li1:li2,n) + &
                           & (vp(l)**2+mu(m)*geom%Bfield(li1:li2,pj1,k))*&
                           & fm_fsa(li1:li2,l,m,n)*mat_00(li1:li2,lj1,k,l,m)*geom%jacobian(li1:li2,pj1,k)
                   end do
                end Do
             end Do
             call my_sum_to_all_real(int_vel_vsq_fm_fsa(li1:li2,n),li0,mpi_comm_vw)
             call my_sum_to_all_real(int_vel_vsq_fm_fsa(li1:li2,n),li0,mpi_comm_z)
             int_vel_vsq_fm_fsa(li1:li2,n)=int_vel_vsq_fm_fsa(li1:li2,n)/(Real(nz0)*geom%avg_jaco_yz(li1:li2))
          end Do

          !combine fm_fsa and int_vel_fm_fsa in a new prefactor
          Do n=ln1,ln2
             Do m=lm1,lm2
                Do l=ll1,ll2
                   k_heat_prefac(li1:li2,l,m,n) = fm_fsa(li1:li2,l,m,n)/int_vel_fm_fsa(li1:li2,n)
                End Do
             End Do
          End Do
          
          sum_vint_fm_fsa=sum(int_vel_fm_fsa(li1:li2,:),2)
          call my_real_sum_spec(sum_vint_fm_fsa,li0)
         
          !combine fm_fsa and int_vel_fm_fsa in a new prefactor
          Do n=ln1,ln2
             Do m=lm1,lm2
                Do l=ll1,ll2
                   k_part_prefac(li1:li2,l,m,n) = fm_fsa(li1:li2,l,m,n)/sum_vint_fm_fsa
                End Do
             End Do
          End Do

       end IF

!    Deallocate(fm_fsa)

    PERFOFF

  end Subroutine init_k_heat_src

 !> Compute f(x,|vpar|,mu) = f(x,vpar,mu) + f(x,- vpar,mu) for all species
  subroutine compute_f_vabs_allspec(p_f_) 
    COMPLEX, DIMENSION(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),INTENT(IN):: p_f_ !< distribution function
    real, dimension(li1:li2) :: int_vel_vsq_krookheat, int_vel_vsq_krookpart

    integer :: i,n

    PERFON('f_vabs')

    if(.not. is_grid_adptv) then
       do n=ln1,ln2
          call compute_f_vabs_n(p_f_,n)
       enddo
    else
       do n=ln1,ln2
          call compute_f_vabs_n_adptv(p_f_,n)
       enddo
    endif

    !compute total perturbed gyrocenter charge density
    do i=li1,li2
       sum_q_vint_f_vabs(i)=sum(spec(ln1:ln2)%charge*int_vel_f_vabs(i,ln1:ln2))
       sum_vint_f_vabs(i)=sum(int_vel_f_vabs(i,ln1:ln2))
    enddo
    call my_complex_sum_spec(sum_q_vint_f_vabs,li0)
    call my_complex_sum_spec(sum_vint_f_vabs,li0)

    !if psource_type=1: 
    !compute coefficient to Krook heat source term, correcting for any heat added by the
    !Krook particle source
    !if psource_type=2: no corrections to the heat source prefactor
    if (psource_type==1) then
       do n=ln1,ln2
          int_vel_vsq_krookheat = (int_vel_vsq_f_vabs(li1:li2,n) - int_vel_f_vabs(li1:li2,n) * &
               int_vel_vsq_fm_fsa(li1:li2,n) / int_vel_fm_fsa(li1:li2,n))
          int_vel_vsq_krookpart = (int_vel_vsq_f_vabs(li1:li2,n) - sum_q_vint_f_vabs(i) /&
                     n_spec / spec(n)%charge * int_vel_vsq_fm_fsa(li1:li2,n) / int_vel_fm_fsa(li1:li2,n))
          if (ck_heat_smooth.gt.0.0) call smooth_x(int_vel_vsq_krookheat, ck_heat_smooth)
          if (ck_part_smooth.gt.0.0) call smooth_x(int_vel_vsq_krookpart, ck_part_smooth)

          do i=li1,li2
             !denominator has a lower limit of 0.1 to avoid numerical problems
             if (abs(int_vel_vsq_krookheat(i)).gt.0.1) then
                coefficient(i,n)= ck_heat_prof(i) - &
                     !v^2 integrated particle source
                     ck_part_prof(i) * int_vel_vsq_krookpart(i) /&
                     !v^2 integrated heat source 
                     int_vel_vsq_krookheat(i)
             else 
                !if denominator<0.1 (typically the case for very weak fluctuations),
                !we don't correct for the heat introduced by the particle source
                coefficient(i,n)=ck_heat_prof(i)
             endif
          enddo
       enddo
    else
       do n=ln1,ln2
          do i=li1,li2       
             coefficient(i,n)=ck_heat_prof(i)
          enddo
       enddo
    endif
    PERFOFF

  end subroutine compute_f_vabs_allspec

  !> Compute f(x,|vpar|,mu) = f(x,vpar,mu) + f(x,- vpar,mu) for species n
  subroutine compute_f_vabs_n(p_f_,n) 
    COMPLEX, DIMENSION(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),INTENT(IN):: p_f_ !< distribution function
    INTEGER, INTENT(IN) :: n                                                                   !< species index
    INTEGER :: m,l,k,ierr,ind_v,ind_loc,dest_pev,stat(MPI_STATUS_SIZE),tag
    COMPLEX, DIMENSION(li1:li2,nv0-1-ll2:nv0-1-ll1) :: buffer_f

    f_vabs(:,:,:,n) = 0.0


    if (n_procs_v.NE.1) then
       IF (MODULO(n_procs_v,2).eq.1) THEN
          print*, '! The krook type heat source is only implemented for even nproc_v !'
          STOP
       END IF
       do m=lm1,lm2
          do k=lk1,lk2
             do l=ll1,ll2
                ! compute flux surface average
                ! index of -vpar point
                dest_pev=(nv0-1-ll1)/ll0
                tag=359
                ind_loc=nv0-1-ll1+l-ll2
                Call mpi_sendrecv(&
                     p_f_(li1:li2,lj1,k,l,m,n), li0, MPI_COMPLEX_TYPE, dest_pev, tag,&
                     buffer_f(li1:li2,ind_loc), li0, MPI_COMPLEX_TYPE, MPI_ANY_SOURCE, tag,&
                     mpi_comm_v, stat, ierr)
             end do  ! l loop         
             do l=ll1,ll2
                ind_v=nv0-1-l
                f_vabs(li1:li2,l,m,n) =  f_vabs(li1:li2,l,m,n) + 0.5*(p_f_(li1:li2,lj1,k,l,m,n) + &
                     buffer_f(li1:li2,ind_v))*geom%jacobian(li1:li2,pj1,k)/ (REAL(nz0)*geom%avg_jaco_yz(li1:li2))
             end do ! l loop
          end do     ! k loop   
          Call sum_to_all_complex(f_vabs(li1:li2,ll1:ll2,m,n),li0*ll0,mpi_comm_z)
       end do        ! m loop
    else
       do m=lm1,lm2
          do l=ll1,ll2
             ! compute flux surface average
             do k=lk1,lk2
                ! index of -vpar point
                ind_v=nv0-1-l
                f_vabs(li1:li2,l,m,n) =  f_vabs(li1:li2,l,m,n) + 0.5*(p_f_(li1:li2,lj1,k,l,m,n) + &
                     & p_f_(li1:li2,lj1,k,ind_v,m,n))*geom%jacobian(li1:li2,pj1,k)
             end do  ! k loop
             Call sum_to_all_complex(f_vabs(li1:li2,l,m,n),li0,mpi_comm_z)
             f_vabs(li1:li2,l,m,n) = f_vabs(li1:li2,l,m,n) / (REAL(nz0)*geom%avg_jaco_yz(li1:li2))
          end do     ! l loop   
       end do        ! m loop
    end if

    int_vel_f_vabs(:,n)=0.0
    do m=lm1,lm2
       do l=ll1,ll2
          do k=lk1,lk2
             int_vel_f_vabs(li1:li2,n) = int_vel_f_vabs(li1:li2,n) + &
                  & f_vabs(li1:li2,l,m,n)*mat_00(li1:li2,pj1,k,l,m)*geom%jacobian(li1:li2,pj1,k) 
          end do
       end do
    end do
    call my_complex_sum_vw(int_vel_f_vabs(li1:li2,n),li0)
    call sum_to_all_complex(int_vel_f_vabs(li1:li2,n),li0,mpi_comm_z)
    int_vel_f_vabs(li1:li2,n)=int_vel_f_vabs(li1:li2,n)/&
         &(Real(nz0)*geom%avg_jaco_yz(li1:li2))


    int_vel_vsq_f_vabs(:,n)=0.0
    do m=lm1,lm2
       do l=ll1,ll2
          do k=lk1,lk2
             int_vel_vsq_f_vabs(li1:li2,n) = int_vel_vsq_f_vabs(li1:li2,n) + &
                  & (vp(l)**2+mu(m)*geom%Bfield(li1:li2,pj1,k))*f_vabs(li1:li2,l,m,n)*&
                  & mat_00(li1:li2,pj1,k,l,m)*geom%jacobian(li1:li2,pj1,k) 
          end do
       end do
    end do
    call my_complex_sum_vw(int_vel_vsq_f_vabs(li1:li2,n),li0)
    call sum_to_all_complex(int_vel_vsq_f_vabs(li1:li2,n),li0,mpi_comm_z)
    int_vel_vsq_f_vabs(li1:li2,n)=int_vel_vsq_f_vabs(li1:li2,n)/&
         &(Real(nz0)*geom%avg_jaco_yz(li1:li2))

  end subroutine compute_f_vabs_n

  !> Update the profiles for the integrator parts of the particle and heat
  !! sources
  subroutine update_integrator(l,m,n)
    integer, intent(in):: l,m,n
    real, dimension(li1:li2) :: krookheatpart_src

    norm_time=min(intsource_time,time-simtime_start)
    if (norm_time.le.0.) norm_time=dt

    if (heatsrc_on) then
       krookheatpart_src = (f_vabs(li1:li2,l,m,n) - &
            k_heat_prefac(li1:li2,l,m,n) * int_vel_f_vabs(li1:li2,n) )
       if (ck_heat_smooth.gt.0.0) call smooth_x(krookheatpart_src, ck_heat_smooth)
       heat_prof(:,l,m,n)=(1-dt/intsource_time)*heat_prof(:,l,m,n)+krookheatpart_src
    endif

    if (partsrc_on) then
       if (psource_type==2) then
          krookheatpart_src = sum_vint_f_vabs(li1:li2) * k_part_prefac(li1:li2,l,m,n)
       else
          krookheatpart_src = (f_vabs(li1:li2,l,m,n) - sum_q_vint_f_vabs(li1:li2) / n_spec * &
               k_heat_prefac(li1:li2,l,m,n) / spec(n)%charge)
       endif
       if (ck_part_smooth.gt.0.0) call smooth_x(krookheatpart_src, ck_part_smooth)
       part_prof(:,l,m,n)=(1-dt/intsource_time)*part_prof(:,l,m,n)+krookheatpart_src       
    endif
  end subroutine update_integrator

  !-----------------------------------------------
  !> Add krook-type heat source
  subroutine add_ck_heat(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs
    real, dimension(li1:li2) :: krookheat_src
    
    krookheat_src = (f_vabs(li1:li2,l,m,n) - &
         k_heat_prefac(li1:li2,l,m,n) * int_vel_f_vabs(li1:li2,n) )
    if (ck_heat_smooth.gt.0.0) call smooth_x(krookheat_src, ck_heat_smooth)

    p_rhs(li1:li2) = p_rhs(li1:li2) - coefficient(li1:li2,n) * krookheat_src(li1:li2) 
    
    if (intsource_heat_coeff.gt.0.0) p_rhs(li1:li2) = p_rhs(li1:li2) - &
         &intsource_heat_coeff / norm_time  * heat_prof(:,l,m,n)

  end subroutine add_ck_heat
  
  !> Add krook-type particle source
  subroutine add_ck_part_1(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    real, dimension(li1:li2) :: krookpart_src
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs
    
    krookpart_src = (f_vabs(li1:li2,l,m,n) - sum_q_vint_f_vabs(li1:li2) / n_spec * &
         k_heat_prefac(li1:li2,l,m,n) / spec(n)%charge)
    if (ck_part_smooth.gt.0.0) call smooth_x(krookpart_src, ck_part_smooth)

    p_rhs(li1:li2) = p_rhs(:) - ck_part * krookpart_src
   
    if (intsource_part_coeff.gt.0.0) p_rhs(li1:li2) = p_rhs(:) - &
         &intsource_part_coeff / norm_time * part_prof(:,l,m,n)

  end subroutine add_ck_part_1

  !> Add krook-type particle source
  subroutine add_ck_part_2(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    real, dimension(li1:li2) :: krookpart_src
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs

    krookpart_src = sum_vint_f_vabs(li1:li2) * k_part_prefac(li1:li2,l,m,n)
    if (ck_part_smooth.gt.0.0) call smooth_x(krookpart_src, ck_part_smooth)
    
    p_rhs(li1:li2) = p_rhs(:) - ck_part * krookpart_src

    if (intsource_part_coeff.gt.0.0) p_rhs(li1:li2) = p_rhs(:) - &
         &intsource_part_coeff / norm_time * part_prof(li1:li2,l,m,n)

  end subroutine add_ck_part_2

  
  !> Add krook-type g1 energy source (for reconnection/current sheets)
  SUBROUTINE add_g1_en_source(k,l,m,n,p_rhs)
    INTEGER, INTENT(IN) :: k, l, m, n !< z, vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2), INTENT(INOUT) :: p_rhs !< local rhs

    ! if first run in series, use initial distribution function
    ! follow-up runs: g1_init set in initial_value_comp before checkpoint is read
    IF (time .EQ. 0.0) THEN
      IF (.NOT. delivered_g1ens_msg) THEN
        IF (mype .EQ. 0) PRINT*, 'setting g1_init for g1_en_source'
        delivered_g1ens_msg = .true.
      END IF
      IF (p_has_0_mode) g1_init(li1:li2,k,l,m,n) = g_1(li1:li2,lj1,k,l,m,n)
    ELSE
      IF (p_has_0_mode) p_rhs(li1:li2) = p_rhs(li1:li2) - coeff_g1_en_source * &
        (g_1(li1:li2,lj1,k,l,m,n) - g1_init(li1:li2,k,l,m,n))
    END IF

  END SUBROUTINE add_g1_en_source
  

  !> Add heat sources
  subroutine add_f1_sources_klmn(k,l,m,n,p_rhs,stage)
    !real, intent(in) :: time
    integer, intent(in) :: k,l,m,n,stage   !< z, vpar, mu, species indices; Runge-Kutta stage
    COMPLEX, DIMENSION(li1:li2,lj1:lj2),intent(INOUT) :: p_rhs !< local rhs     

    PERFON('add_src')
    if (p_has_0_mode) then
       if (stage==1.and.k==lk1) call update_integrator(l,m,n)
       if (heatsrc_on) &
            call add_ck_heat(l,m,n,p_rhs(li1:li2,lj1))
       if (partsrc_on.and.psource_type==1) &
            call add_ck_part_1(l,m,n,p_rhs(li1:li2,lj1))
       if (partsrc_on.and.psource_type==2) &
            call add_ck_part_2(l,m,n,p_rhs(li1:li2,lj1))
!       if (f0_heat_src_active) call add_heat_src_term(time,k,l,m,n,p_rhs(li1:li2,lj1))
       IF (coeff_g1_en_source .GT. 0.0) CALL add_g1_en_source(k,l,m,n,p_rhs(li1:li2,lj1))
    endif
    PERFOFF

  end subroutine add_f1_sources_klmn


  !-----------------------------------------------
  !> Add heat sources (block version)
  !! \todo this routine could be better incorporated in the block structures
  subroutine add_f1_sources_block(p_rhs, lb1, lb2,stage)
    !real, intent(in) :: time
    integer, intent(in) :: lb1,lb2,stage
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lb1:lb2),intent(INOUT) :: p_rhs

    integer :: klmn

    PERFON('add_src')
    if (p_has_0_mode) then
       do klmn=lb1,lb2
          if (stage==1.and.sk(klmn)==lk1) call update_integrator(sl(klmn),sm(klmn),sn(klmn))    
          if (heatsrc_on) call add_ck_heat(sl(klmn),sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          if (partsrc_on.and.psource_type==1) call add_ck_part_1(sl(klmn),&
               & sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          if (partsrc_on.and.psource_type==2) call add_ck_part_2(sl(klmn),&
               & sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          IF (coeff_g1_en_source .GT. 0.0) CALL add_g1_en_source(&
            sk(klmn),sl(klmn),sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
       enddo
    endif
    
    PERFOFF

  end subroutine add_f1_sources_block


!*****************************************************************************
! Damping Krook subroutines
  
  !>Initialize the Krook operator
  subroutine initialize_KrookBuffer_operator
    integer :: i,k,l,m,n,l_nxbuf,l_ibuf,u_nxbuf,u_ibuf,klmn

    if (init_status_krookbuffer_operator.eq.1) return
    if (explicit_buffer) then
       allocate(map_to_f(lklmn0))
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2 
                do k=lk1,lk2
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                enddo
             enddo
          enddo
       enddo
    endif
    ALLOCATE(nu_K(pi1:pi2))
    nu_k=0.0
    
    ! check buffer size
    if ((abs(l_buffer_size).ge.1).or.((abs(u_buffer_size).ge.1))) then
       print*, l_buffer_size, u_buffer_size
       stop 'buffer_size should be between 0 and 1'
    end if
    
    ! find index of the buffer region's end
    l_nxbuf=int(l_buffer_size*nx0)
    l_ibuf =pi1gl+l_nxbuf
    u_nxbuf=int(u_buffer_size*nx0)
    u_ibuf =pi2gl-u_nxbuf
    
    if (l_nxbuf.gt.0) then
       do i=pi1gl, pi1gl+l_nxbuf
          If (i<pi1 .Or. i>pi2) Cycle
          nu_K(i)=lcoef_krook*(abs(i-l_ibuf)/real(l_nxbuf))**lpow_krook
       end do
    endif
    
    if (u_nxbuf.gt.0) then
       do i=pi2gl-u_nxbuf, pi2gl
          If (i<pi1 .Or. i>pi2) Cycle
          nu_K(i)=ucoef_krook*(abs(i-u_ibuf)/real(u_nxbuf))**upow_krook
       end do
    endif

    !    do i=pi1,pi2
    !       print*, i, nu_k(i)
    !    end do
    
    init_status_krookbuffer_operator = 1

  end subroutine initialize_KrookBuffer_operator

  subroutine add_kBuffer_explicit(p_f,p_rhs,lb1,lb2)
    integer,intent(in):: lb1,lb2
    complex,dimension(li1:li2,lj1:lj2,lb1:lb2):: p_rhs
    complex,dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    integer:: j,klmn


    do klmn=lb1,lb2
       do j=lj1,lj2
          if (.not.buffer_on_ky0.or.(j==lj1.and.p_has_0_mode)) &
          p_rhs(pi1:pi2,j,klmn)=p_rhs(pi1:pi2,j,klmn)-nu_K(pi1:pi2)*p_f(pi1:pi2,j,map_to_f(klmn))
       enddo
    enddo

  end subroutine add_kBuffer_explicit

  !>add Krook operator (used in time scheme)
  subroutine add_krookBuffer(p_g1)
    complex,dimension(li1:li2,0:ljklmn0-1),intent(inout):: p_g1
    integer(KIND=8) :: ind

    do ind=0,ljklmn0-1
       p_g1(pi1:pi2,ind) = exp(-nu_K(pi1:pi2)*dt)*p_g1(pi1:pi2,ind)
    enddo

  end subroutine add_krookBuffer
  
  subroutine finalize_krookBuffer_operator
    if (explicit_buffer) deallocate(map_to_f)
    deallocate(nu_K)
    init_status_krookbuffer_operator = 0
  end subroutine finalize_krookBuffer_operator

!*****************************************************************************

  subroutine initialize_f0_heat_source(mp_f0_contr)
    real, dimension(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), intent(inout) :: mp_f0_contr
    integer :: s, i,j,k,l,m,n, prof_type, ierr
    real, dimension(:),allocatable :: src_prof, src_prof_unnorm
    real :: x, src_x_norm_loc, src_x_norm, value, area_center_loc
    real :: area_center, chi
    real, dimension(0:2) :: src_mom, src_mom_tmp

    if (init_status_heat_source.eq.1) return

    if ((puls_amp.ne.1.0).and.(mype.eq.0)) then
       write (*,'(A)') "Warning: the pulsed heat source is switched on "
       write (*,'(A)') "         which would also affect the F0 term in the"
       write (*,'(A)') "         Vlasov eq. if activated (check this!)"
    endif

    allocate(src_prof(pi1:pi2), src_prof_unnorm(pi1:pi2))

    Do n=ln1,ln2
       !skip species which have no source/sink
       IF (SUM(spec(n)%src_prof_type).eq.0.0) CYCLE

       src_prof = 0.0
       !construct sources & sinks profile
       !loop over up to five functions
       Do s=0,SIZE(spec(n)%src_prof_type)-1
          prof_type = spec(n)%src_prof_type(s)          
          if (prof_type.eq.0) CYCLE

          src_x_norm_loc = 0.0
          src_prof_unnorm = 0.0
          
          Do i=pi1,pi2
             x=x0+rhostar*(real(i-pi1gl)*deli-0.5D0*lp1)
             select case (prof_type)
             case(1) !Gaussian heat source profile
                value = exp(-4.0*LOG(2.)*((x-spec(n)%src_x0(s))/spec(n)%src_width(s))**2)
             case(2) !GYSELA type heat source profile, see Y. Sarazin, accepted for pub. in Nucl. Fusion (2010)
                value = -0.5*(tanh((x-spec(n)%src_x0(s)-3*spec(n)%src_width(s))/&
                     &spec(n)%src_width(s))+tanh((spec(n)%src_x0(s)-3*spec(n)%src_width(s)-x)/&
                     &spec(n)%src_width(s)))
             case default !no sources & sinks
                stop 'this source profile type is not implemented'
             end select
             
             src_prof_unnorm(i) = value
             src_x_norm_loc = src_x_norm_loc + value*geom%avg_jaco_yz(i)
          End Do

          !now we need to sum over all x processes
          Call mpi_allreduce(src_x_norm_loc,src_x_norm, 1,&
               MPI_REAL_TYPE, MPI_SUM, mpi_comm_x, ierr)
          src_x_norm = lx*abs(n0_global)*ly*2.0*pi*n_pol*src_x_norm/nx0
          !src_norm now contains the full "int dx dy dz S(x) J(x,z)"

          ! ... and add the normalized heat source profile multiplied by 
          ! the desired power amplitude to the total source profile
          src_prof(pi1:pi2) = src_prof(pi1:pi2) + &
               &spec(n)%src_amp(s)*src_prof_unnorm(pi1:pi2)/src_x_norm

       End Do !loop over different profiles
       
       !set values smaller than machine precision to zero 
       !(to avoid problems with ASCII file I/O)
       Do i=pi1,pi2
          if (ABS(src_prof(i)).LT.EPSILON(src_prof(i))) src_prof(i)=0.0
       Enddo

       if (print_ini_msg) call write_f0_heat_source_profile(src_prof,n)

       !now add the normalization of the energy term (function of x)

       src_prof(pi1:pi2) = src_prof(pi1:pi2)*&
            &2./(3.*spec(n)%dens*spec(n)%dens_prof(pi1:pi2)*&
            &spec(n)%temp*spec(n)%temp_prof(pi1:pi2))

       src_mom = 0.0

       !in principle, one could construct a new prefactor to increase readability
       !however, it would be high-dimensional and thus memory consuming
       !thus, we simply add the source to the F0 contribution

       Do m=lm1,lm2
          Do l=ll1,ll2
             Do k=lk1,lk2
                Do j=pj1,pj2
                   Do i=pi1,pi2
                      value = src_prof(i) * ((mu(m)*geom%Bfield(i,j,k)+vp(l)**2)/&
                           & spec(n)%temp_prof(i)-1.5)*fm(i,j,k,l,m,n)
                      mp_f0_contr(i,j,k,l,m,n) = mp_f0_contr(i,j,k,l,m,n) + &
                           & value
                      
                      !calculate some moments to check whether the numerical results
                      !match the analytical predictions
                      value = value*geom%jacobian(i,j,k)
                      src_mom(0) = src_mom(0) + value*mat_00(i,j,k,l,m)
                      src_mom(1) = src_mom(1) + value*mat_10(i,j,k,l,m,n)
                      src_mom(2) = src_mom(2) + value*(mu(m)*geom%Bfield(i,j,k)+vp(l)**2)*&
                           &mat_00(i,j,k,l,m)
                   End Do
                End Do
             End Do
          End Do
       End Do

       src_mom = src_mom * lx / nx0 * ly * abs(n0_global) * 2.0 * pi * n_pol / nz0

       do s=0,2
          src_mom_tmp(s) = sum_vw(src_mom(s))
          Call mpi_allreduce(src_mom_tmp(s),src_mom(s),1,&
               MPI_REAL_TYPE, MPI_SUM, mpi_comm_z, ierr)
          call mpi_allreduce(src_mom(s),src_mom_tmp(s),1,&
               MPI_REAL_TYPE, MPI_SUM, mpi_comm_x, ierr)
       end do

       area_center_loc = ly*abs(n0_global)*2.0*pi*n_pol*&
            & SUM(geom%jacobian((pi1gl+pi2gl)/2,pj1,:))/real(nz0)
       Call mpi_allreduce(area_center_loc,area_center,1,&
            MPI_REAL_TYPE, MPI_SUM, mpi_comm_z, ierr)

       If (write_pe) then
          write (*,*) ' '
          write (*,'(2A)') 'heat source for ', spec(n)%name
          write(*,'(A,ES12.4)') &
               & '0th phase space moment (full flux surface):    ',src_mom_tmp(0)
          write(*,'(A,ES12.4)') &
               & '1st phase space moment (full flux surface):    ',src_mom_tmp(1)
          write(*,'(A,ES12.4)') &
               & 'energy phase space moment (full flux surface): ',src_mom_tmp(2)
          write (*,*) ' '

          write(*,'(A,ES12.4)') &
               & 'projected heat flux at center of sim. domain:  ',&
               & src_mom_tmp(2)/area_center
          chi = src_mom_tmp(2)/(area_center*spec(n)%dens*&
               & spec(n)%dens_prof((pi1gl+pi2gl)/2)*spec(n)%temp*&
               & spec(n)%temp_prof((pi1gl+pi2gl)/2)*spec(n)%omt_prof((pi1gl+pi2gl)/2))
          write(*,'(A,ES12.4)') &
               & 'estimated heat diff. at center of sim. domain: ',chi
          IF ((chi*(rhostar/minor_r)**2).GT.1e-12) THEN
             write(*,'(A,ES12.4)') &
                  & 'estimated energy confinement time:             ',&
                  & (minor_r/rhostar)**2/chi
          end IF

       endif

    End Do !n loop

    deallocate(src_prof, src_prof_unnorm)

    init_status_heat_source = 1

  end subroutine initialize_f0_heat_source


  subroutine write_f0_heat_source_profile(src_prof,n)
    real, dimension(pi1:pi2), intent(in) :: src_prof
    integer, intent(in) :: n
    !local variables
    real, dimension(pi1gl:pi2gl) :: src_prof_glob
    real :: x
    integer :: SRC_PROF_FILE, i, ierr

    !construct full array in x on my_pex=0
    call mpi_gather(src_prof,pi0, MPI_REAL_TYPE,&
         src_prof_glob,pi0,MPI_REAL_TYPE,0, mpi_comm_x, ierr)

    if (write_pe) then
       call get_unit_nr(SRC_PROF_FILE)
       OPEN(SRC_PROF_FILE,file=trim(diagdir)//'/src_profiles_'//&
            &trim(spec(n)%name)//trim(file_extension))
       write(SRC_PROF_FILE,"(A)") "#     x           xrhostar        src"
       Do i=pi1gl,pi2gl
          x=x0+rhostar*(real(i-pi1gl)*deli-0.5*lp1)
          write(SRC_PROF_FILE,"(3ES16.6)") x,x/rhostar,src_prof_glob(i)
       End Do
       close(SRC_PROF_FILE)
    endif
    
  end subroutine write_f0_heat_source_profile

  subroutine finalize_f0_heat_source
    init_status_heat_source = 0
  end subroutine finalize_f0_heat_source


  !> Function modelling puls func
  !! Currently, a step function is implemented
  real function puls_func(time)
    real, intent(in) :: time
    
    if (in_timestep_comp) then !just return the maximum of the pulse
                               !for the linear time step computation
       puls_func=max(1.0,puls_amp)
    else
       if (mod((time-puls_start),(puls_on_duration+puls_off_duration))&
            .lt.puls_on_duration) then
          puls_func = puls_amp
       else
          puls_func = 1.0
       endif
    endif

  end function puls_func
     
  subroutine smooth_x (data, relwidth)
    real, dimension(li1:li2), intent(inout) :: data !< in/out array
    real, intent(in) :: relwidth
    real, dimension(0:nx0-1) :: full_data
    integer :: ierr

    if (n_procs_x.eq.1) then
       call smooth(data, nx0, relwidth, rad_bc_type, ck_filter_type)
    else
       CALL mpi_allgather(data,li0,MPI_REAL_TYPE,full_data,li0,&
            MPI_REAL_TYPE,mpi_comm_x,ierr)
       call smooth(full_data, nx0, relwidth, rad_bc_type, ck_filter_type)
       data = full_data(li1:li2)      
    endif

  end subroutine smooth_x

  !> Smoothing routine
  !> \param data array containing data to be smoothed
  !> \param ndata number of elements in (local part of the) array
  !> \param relwidth rel. fraction of the global array to be smoothed
  !> \param bc_type boundary condition
  !> \param filter_type: 0=boxcar average, 1=Gaussian
  subroutine smooth(data, ndata, relwidth, bc_type, filter_type)
    integer, intent(in) :: ndata
    real, dimension(1:ndata), intent(inout) :: data !< in/out array
    real, intent(in) :: relwidth
    integer, intent(in) :: bc_type
    integer, intent(in):: filter_type 
    real, dimension(1:ndata) :: temp, weight
    integer :: i, o, w, locw

    w = int(relwidth*ndata)/2

    select case (bc_type)
    case (1) !Dirichlet b.c.: zero value outside box
       if (filter_type==0) then
          do i = 1, ndata
             locw = w - max(-(i-w-1),i+w-ndata,0)
             temp(i) = sum(data(i-locw:i+locw))/real(2*w+1)
          enddo
       else
          do i = 1, ndata
             locw = w - max(-(i-w-1),i+w-ndata,0)
             do o=1,2*locw+1
                weight(o)=exp(-(o-locw-1)**2/0.3/real(relwidth*ndata/2)**2)
             enddo
             weight=weight/sum(weight(1:2*locw+1))
             temp(i) = sum(data(i-locw:i+locw)*weight(1:2*locw+1))
          enddo
       endif
    case (2) !v Neumann b.c.
       stop 'not implemented yet'
    case default
       !reduce smoothing window width when hitting the boundary
       do i = 1, ndata
          locw = w - MAX(-(i-w-1),i+w-ndata,0)
          temp(i) = sum(data(i-locw:i+locw))/real(2*locw+1)
       enddo
    end select

    data = temp

  end subroutine smooth

  ! ---- modified routines for adaptive grids

  subroutine compute_f_vabs_n_adptv(p_f_,n) 
    COMPLEX, DIMENSION(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2),INTENT(IN):: p_f_ !< distribution function
    INTEGER, INTENT(IN) :: n                                                                   !< species index
    INTEGER :: m,l,k,ierr,ind_v,ind_loc,dest_pev,stat(MPI_STATUS_SIZE),tag
    COMPLEX, DIMENSION(li1:li2,nv0-1-ll2:nv0-1-ll1) :: buffer_f
    integer :: i1, i2

    f_vabs(:,:,:,n) = 0.0

    if (n_procs_v.NE.1) then
       IF (MODULO(n_procs_v,2).eq.1) THEN
          print*, '! The krook type heat source is only implemented for even nproc_v !'
          STOP
       END IF
       do m=lm1,lm2
          do k=lk1,lk2
             do l=ll1,ll2
                i1 = li1_vwadp(l,m)
                i2 = li2_vwadp(l,m)
                ! compute flux surface average
                ! index of -vpar point
                dest_pev=(nv0-1-ll1)/ll0
                tag=359
                ind_loc=nv0-1-ll1+l-ll2
                Call mpi_sendrecv(&
                     p_f_(i1:i2,lj1,k,l,m,n), li0_vwadp(l,m), MPI_COMPLEX_TYPE, &
                     & dest_pev, tag, buffer_f(i1:i2,ind_loc), li0_vwadp(l,m), &
                     & MPI_COMPLEX_TYPE, MPI_ANY_SOURCE, tag, &
                     & mpi_comm_v, stat, ierr)
             end do  ! l loop         
             do l=ll1,ll2
                i1 = li1_vwadp(l,m)
                i2 = li2_vwadp(l,m)
                ind_v=nv0-1-l
                f_vabs(i1:i2,l,m,n) =  &
                     &f_vabs(i1:i2,l,m,n) + 0.5*(p_f_(i1:i2,lj1,k,l,m,n) + &
                     &buffer_f(i1:i2,ind_v))*geom%jacobian(i1:i2,pj1,k)/ &
                     &(REAL(nz0)*geom%avg_jaco_yz(i1:i2))
             end do ! l loop
          end do     ! k loop   
          Call sum_to_all_complex(f_vabs(li1:li2,ll1:ll2,m,n),li0*ll0,mpi_comm_z)
       end do        ! m loop
    else
       do m=lm1,lm2
          do l=ll1,ll2
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             ! compute flux surface average
             do k=lk1,lk2
                ! index of -vpar point
                ind_v=nv0-1-l
                f_vabs(i1:i2,l,m,n) =  f_vabs(i1:i2,l,m,n) + &
                     &0.5*(p_f_(i1:i2,lj1,k,l,m,n) + &
                     &p_f_(i1:i2,lj1,k,ind_v,m,n))*geom%jacobian(i1:i2,pj1,k)
             end do  ! k loop
             Call sum_to_all_complex(f_vabs(li1:li2,l,m,n),li0,mpi_comm_z)
             f_vabs(i1:i2,l,m,n) = f_vabs(i1:i2,l,m,n) / &
                  &(REAL(nz0)*geom%avg_jaco_yz(i1:i2))
          end do     ! l loop   
       end do        ! m loop
    end if

    int_vel_f_vabs(:,n)=0.0
    do m=lm1,lm2
       do l=ll1,ll2
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do k=lk1,lk2
             int_vel_f_vabs(i1:i2,n) = int_vel_f_vabs(i1:i2,n) + &
                  &f_vabs(i1:i2,l,m,n)*mat_00(i1:i2,pj1,k,l,m)*&
                  &geom%jacobian(i1:i2,pj1,k) 
          end do
       end do
    end do
    call my_complex_sum_vw(int_vel_f_vabs(li1:li2,n),li0)
    call sum_to_all_complex(int_vel_f_vabs(li1:li2,n),li0,mpi_comm_z)
    int_vel_f_vabs(li1:li2,n)=int_vel_f_vabs(li1:li2,n)/&
         &(Real(nz0)*geom%avg_jaco_yz(li1:li2))


    int_vel_vsq_f_vabs(:,n)=0.0
    do m=lm1,lm2
       do l=ll1,ll2
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do k=lk1,lk2
             int_vel_vsq_f_vabs(i1:i2,n) = int_vel_vsq_f_vabs(i1:i2,n) + &
                  & (vp(l)**2+mu(m)*geom%Bfield(i1:i2,pj1,k))*f_vabs(i1:i2,l,m,n)*&
                  & mat_00(i1:i2,pj1,k,l,m)*geom%jacobian(i1:i2,pj1,k) 
          end do
       end do
    end do
    call my_complex_sum_vw(int_vel_vsq_f_vabs(li1:li2,n),li0)
    call sum_to_all_complex(int_vel_vsq_f_vabs(li1:li2,n),li0,mpi_comm_z)
    int_vel_vsq_f_vabs(li1:li2,n)=int_vel_vsq_f_vabs(li1:li2,n)/&
         &(Real(nz0)*geom%avg_jaco_yz(li1:li2))

  end subroutine compute_f_vabs_n_adptv

  subroutine add_kBuffer_explicit_adptv(p_f,p_rhs,lb1,lb2)
    integer,intent(in):: lb1,lb2
    complex,dimension(li1:li2,lj1:lj2,lb1:lb2):: p_rhs
    complex,dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    integer:: j, klmn, l, m, pi1_a, pi2_a


    do klmn=lb1,lb2
       l = sl(klmn)
       m = sm(klmn)
       pi1_a = pi1_vwadp(l,m)
       pi2_a = pi2_vwadp(l,m)
       do j=lj1,lj2
          p_rhs(pi1_a:pi2_a,j,klmn)=p_rhs(pi1_a:pi2_a,j,klmn)-nu_K(pi1_a:pi2_a)*p_f(pi1_a:pi2_a,j,map_to_f(klmn))
       enddo
    enddo

  end subroutine add_kBuffer_explicit_adptv

  !-----------------------------------------------
  !> Add heat sources (block version)
  !! \todo this routine could be better incorporated in the block structures
  subroutine add_f1_sources_block_adptv(p_rhs,lb1,lb2,stage)
    !real, intent(in) :: time
    integer, intent(in) :: lb1,lb2,stage
    COMPLEX, DIMENSION(li1:li2,lj1:lj2,lb1:lb2),intent(INOUT) :: p_rhs

    integer :: klmn

    PERFON('add_src')
    if (p_has_0_mode) then
       do klmn=lb1,lb2
          if (stage==1.and.sk(klmn)==lk1) call update_integrator_adptv(sl(klmn),sm(klmn),sn(klmn))    
          if (heatsrc_on) call add_ck_heat_adptv(sl(klmn),sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          if (partsrc_on.and.psource_type==1) call add_ck_part_1_adptv(sl(klmn),&
               & sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          if (partsrc_on.and.psource_type==2) call add_ck_part_2_adptv(sl(klmn),&
               & sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
          IF (coeff_g1_en_source .GT. 0.0) CALL add_g1_en_source_adptv(&
            sk(klmn),sl(klmn),sm(klmn),sn(klmn),p_rhs(li1:li2,lj1,klmn))
       enddo
    endif
    
    PERFOFF

  end subroutine add_f1_sources_block_adptv

  !> Update the profiles for the integrator parts of the particle and heat
  !! sources
  subroutine update_integrator_adptv(l,m,n)
    integer, intent(in):: l,m,n
    integer :: i1, i2
    real, dimension(:), allocatable :: krookheatpart_src

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    allocate(krookheatpart_src(i1:i2))

    if (heatsrc_on) then
       krookheatpart_src = (f_vabs(i1:i2,l,m,n) - &
            k_heat_prefac(li1:li2,l,m,n) * int_vel_f_vabs(i1:i2,n) )
       if (ck_heat_smooth.gt.0.0) call smooth(krookheatpart_src,i2-i1+1, ck_heat_smooth,rad_bc_type,ck_filter_type)
       heat_prof(i1:i2,l,m,n)=(1-dt/intsource_time)*heat_prof(i1:i2,l,m,n)+krookheatpart_src
    endif

    if (partsrc_on) then
       if (psource_type==2) then
          krookheatpart_src = sum_vint_f_vabs(i1:i2) * k_part_prefac(i1:i2,l,m,n)
       else
          krookheatpart_src = (f_vabs(i1:i2,l,m,n) - sum_q_vint_f_vabs(i1:i2) / n_spec * &
               k_heat_prefac(i1:i2,l,m,n) / spec(n)%charge)
       endif
       if (ck_part_smooth.gt.0.0) call smooth(krookheatpart_src,i2-i1+1, ck_part_smooth,rad_bc_type,ck_filter_type)
       part_prof(i1:i2,l,m,n)=(1-dt/intsource_time)*part_prof(i1:i2,l,m,n)+krookheatpart_src       
    endif
    
    deallocate(krookheatpart_src)
  end subroutine update_integrator_adptv

  !-----------------------------------------------

  ! modified routines for adaptive grids
  !>add Krook operator (used in time scheme)
  subroutine add_krookBuffer_adptv(p_g1)
    complex,dimension(li1:li2,lj1:lj2,1:lklmn0),intent(inout):: p_g1
    integer(KIND=8) :: ind, j
    integer :: l, m, p1, p2

    do ind=1,lklmn0
       l = sl(ind)
       m = sm(ind)
       p1 = pi1_vwadp(l,m)
       p2 = pi2_vwadp(l,m)
       do j = lj1, lj2
          p_g1(p1:p2,j,ind) = exp(-nu_K(p1:p2)*dt)*p_g1(p1:p2,j,ind)
       end do
    end do

  end subroutine add_krookBuffer_adptv
  
  !> Add krook-type heat source
  subroutine add_ck_heat_adptv(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs
    real, dimension(:), allocatable :: krookheat_src
    integer :: i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    allocate(krookheat_src(i1:i2))

    krookheat_src = (f_vabs(i1:i2,l,m,n) - &
         k_heat_prefac(i1:i2,l,m,n) * int_vel_f_vabs(i1:i2,n) )
    if (ck_heat_smooth.gt.0.0) call smooth(krookheat_src,i2-i1+1, ck_heat_smooth,rad_bc_type,ck_filter_type)

    p_rhs(i1:i2) = p_rhs(i1:i2) - coefficient(i1:i2,n) * krookheat_src(i1:i2) 
    
    if (intsource_heat_coeff.gt.0.0) p_rhs(i1:i2) = p_rhs(i1:i2) - &
         &intsource_heat_coeff / norm_time  * heat_prof(:,l,m,n)

    deallocate(krookheat_src)
  end subroutine add_ck_heat_adptv

    !> Add krook-type particle source
  subroutine add_ck_part_1_adptv(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs
    real, dimension(:), allocatable :: krookpart_src
    integer :: i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    allocate(krookpart_src(i1:i2))    

    krookpart_src = (f_vabs(i1:i2,l,m,n) - sum_q_vint_f_vabs(i1:i2) / n_spec * &
         k_heat_prefac(i1:i2,l,m,n) / spec(n)%charge)
    if (ck_part_smooth.gt.0.0) call smooth(krookpart_src,i2-i1+1, ck_part_smooth,rad_bc_type,ck_filter_type)
   
    p_rhs(i1:i2) = p_rhs(i1:i2) - ck_part * krookpart_src

    if (intsource_part_coeff.gt.0.0) p_rhs(i1:i2) = p_rhs(i1:i2) - &
         &intsource_part_coeff / norm_time * part_prof(i1:i2,l,m,n)
    
    deallocate(krookpart_src)

  end subroutine add_ck_part_1_adptv

  subroutine add_ck_part_2_adptv(l,m,n,p_rhs)
    integer, intent(in) :: l,m,n   !< vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2),intent(INOUT) :: p_rhs !< local rhs
    real, dimension(:), allocatable :: krookpart_src
    integer :: i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    allocate(krookpart_src(i1:i2))

    krookpart_src = sum_vint_f_vabs(i1:i2) * k_part_prefac(i1:i2,l,m,n)
    if (ck_part_smooth.gt.0.0) call smooth(krookpart_src,i2-i1+1, ck_part_smooth,rad_bc_type,ck_filter_type)
   
    p_rhs(i1:i2) = p_rhs(i1:i2) - ck_part * krookpart_src

    if (intsource_part_coeff.gt.0.0) p_rhs(i1:i2) = p_rhs(i1:i2) - &
         &intsource_part_coeff / norm_time * part_prof(i1:i2,l,m,n)
    
    deallocate(krookpart_src)

  end subroutine add_ck_part_2_adptv

  !> Add krook-type g1 energy source (for reconnection/current sheets)
  SUBROUTINE add_g1_en_source_adptv(k,l,m,n,p_rhs)
    INTEGER, INTENT(IN) :: k, l, m, n !< z, vpar, mu, species indices
    COMPLEX, DIMENSION(li1:li2), INTENT(INOUT) :: p_rhs !< local rhs

    integer :: i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)

    ! if first run in series, use initial distribution function
    ! follow-up runs: g1_init set in initial_value_comp before checkpoint is read
    IF (time .EQ. 0.0) THEN
      IF (.NOT. delivered_g1ens_msg) THEN
        IF (mype .EQ. 0) PRINT*, 'setting g1_init for g1_en_source'
        delivered_g1ens_msg = .true.
      END IF
      IF (p_has_0_mode) g1_init(i1:i2,k,l,m,n) = g_1(i1:i2,lj1,k,l,m,n)
    ELSE
      IF (p_has_0_mode) p_rhs(i1:i2) = p_rhs(i1:i2) - coeff_g1_en_source * &
        (g_1(i1:i2,lj1,k,l,m,n) - g1_init(i1:i2,k,l,m,n))
    END IF

  END SUBROUTINE add_g1_en_source_adptv


end Module sources_mod
