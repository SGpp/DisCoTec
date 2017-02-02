#include "redef.h"

!>Computes or estimates the maximal timestep that can be used in explicit initial 
!!value computations 
!!
!!The estimate is also used to set dt for the performance optimization
!!\todo write a maxclimits routine which doesn't need full initialize_current_parall
module compute_dt
  use par_mod
  use discretization
  use coordinates
  use communications
  use collisions, only: ev_coll_est
  use phys_ini
  use aux_fields
  use eigen_parameters
  use RK_coefficients
  use calc_rhs, only: rhs_nl, rhs_f0, rhs_only_coll
  use external_contr, only: ExBrate, ExB_stime
  use numerical_damping
  use x_derivatives, only: max_ddi
#ifdef WITHSLEPC
  use petsc_aux
  use slepc_aux
#endif

  implicit none
  public:: estimate_dt_max, set_dt_max, check_compute_dt
  public:: initialize_adapt_dt, adapt_dt, finalize_adapt_dt
#ifdef WITHSLEPC
  public:: dt_from_ev
#endif
  private

  !quantities needed for nonlinear time step adaption
  integer:: win_dts=100, win_est=25
  real,dimension(:),allocatable:: last_dts, last_estimates
  real:: ltav
  real,dimension(4):: dt_coll=0.
  logical:: write_pe, coll_scheme_adaptive
  !integer:: TIMEFILE=78

contains

  subroutine check_compute_dt
    if (comp_type=='IV') then
       if (calc_dt) then
          if(n_ev.eq.0) then
#ifdef oldparbc
             n_ev=1
#else
             n_ev=20
#endif
          end if
       else !calc_dt=.f.
          if ((collision_op.ne.'none').and.coll_split) then 
             if (dt_vlasov.le.0.) dt_vlasov=dt_max
             if ((ev_coll==-1.).or.(dt_vlasov.le.0.)) then 
                if (mype==0) Write(*,'(A)') "ERROR: specify ev_coll and dt_vlasov, or enable automatic calculation of time step"
                stop
             endif
          else
             !without collisions, or with collisions but no splitting
             if (dt_max.le.0.) stop "ERROR: specify dt_max or enable automatic calculation of time step"
          endif
       end if
    end if

    if (nonlinear) then
       if (courant.le.0.0) then
          if(mype==0) Write(*,'(A)') "WARNING: courant<=0 does not make sense, setting courant=1.0"
          courant = 1.0
       endif
       if (courant.lt.1.0) then
          if(mype==0) Write(*,'(A)') "WARNING: courant >= 1.0 is recommended! (if no instability is detected)"
       endif
    endif

  end subroutine check_compute_dt
  
  !>If calc_dt=.true., this routine either computes or estimates the maximum dt for 
  !>numerical stability of the time stepping scheme
  subroutine set_dt_max
    
    if (calc_dt) then   
#ifdef WITHSLEPC
       call compute_dt_max
#else
       call estimate_dt_max(.true.,dt_max,dt_vlasov,ev_coll)
#endif
    else 
       write_pe=((mype.eq.0).and.(print_ini_msg))
       if ((collision_op.ne.'none').and.(coll_split_scheme.ne.'none')) then
          !with coll_split

          call check_dt_vlasov(dt_max,dt_vlasov,write_pe)
          if (ev_coll.le.0.and.write_pe) then
             Write(*,"(a)") "ill defined ev_coll --exit"
             stop
          endif

          call get_dt_coll(ev_coll,dt_vlasov)
          dt_max = min(dt_vlasov,dt_coll(rkstages_coll))

          if (write_pe) then
             write(*,"(a,es10.3)") 'using time step limit from parameters:'
             write(*,"(a,es10.3)") 'collisions: dt_coll = C_RK/ev_coll = ',dt_coll(rkstages_coll)
             write(*,"(a,es10.3)") 'collisionless part: dt_vlasov = ',dt_vlasov
             write(*,"(a,es10.3)") 'total: dt< ',dt_max
             write(*,*)
          endif
       else
          !without collisions, or with collisions but no splitting
          call check_dt_vlasov(dt_max,dt_vlasov,write_pe)
          if (write_pe) &
               write(*,"(a,es10.3)") 'using time step limit from parameters: ',dt_max
       endif
    endif

  end subroutine set_dt_max

  subroutine check_dt_vlasov(dtm,dtv,write_pe)
     real::dtm,dtv
     logical::write_pe
     !set dt_vlasov, in case it is not set
     if (dtv.le.0) then
        if (dtm .gt.0) then
           dtv = dtm
        else
           if(write_pe) Write(*,"(a)") "ill defined dt_max(dt_vlasov) --exit"
           stop
        endif
     endif
  end subroutine check_dt_vlasov


#ifdef WITHSLEPC
  !> Computes dt_max (dt_coll and dt_vlasov) from the eigenvalues of the linear operator
  !!
  !!The eigenvalues with the largest absolute value are computed with SLEPc
  !!(the number can be changed with n_ev). Then, for each eigenvalue the maximal dt 
  !!necessary to keep this eigenvalue in the stability region of the time 
  !!stepping scheme is computed, the minimum of these values is then written to dt_max.
  !!The most popular time stepping schemes supported at the moment are RK3, RK4, RK4M.
  !!Note that this method assumes that the eigenvalues are mostly aligned with the imaginary
  !!axis, so that the ones with the largest magnitude are actually the most critical for stability
  !!(see flm's PhD thesis).
  !!Note: if the spectrum is distorted, e.g. by large hyperdiffusion terms
  !!or dissipative boundary conditions, the evs with the largest magnitude
  !!might not be the most critical for stability, in these cases 
  !!n_ev may have to be set considerably >1 (maybe 50 or so).
  !!
  !!When arakawa_zv (with hyp_on_h) is used for parallel dynamics,
  !!extremely damped eigenmodes may appear that limit the timestep.
  !!In this case, a set of smallest_real eigenvalue runs
  !!clarifies, if the hypz_compensation terms need to be added to avoid a small timestep.
  !!
  !!When operator splitting is used, also the n_ev largest magnitude collisional eigenvalues
  !!are computed with SLEPc, from which the collisional timestep is computed, so that
  !!it falls into the stability region of the chosen collision time scheme.
  !!Reasonable schemes are RKC1, RKC2, RKC3 and RKC4.
  !!The default (adaptive) RKCa selects the most appropriate of the above RKC schemes.
  !!
  !!\todo ky-parallelization (if any) is not exploited
  subroutine compute_dt_max

    double precision:: time_1, time_2
    integer:: ky0_ind_in, n_procs_y_in, n_procs_sim_in
    integer:: j, blocknum, num_ev, ierr
    integer:: mpi_comm_world_save
    integer,dimension(:),allocatable:: kyind_for_ev
    logical:: sav_print_ini_msg, sav_antenna_contrib
    real:: dt_max_in

    PERFON('comp_dt')

    in_timestep_comp = .true.

    write_pe=((mype.eq.0).and.(print_ini_msg))

    sav_print_ini_msg=print_ini_msg
    print_ini_msg = .false.

    if(write_pe) then
       write(*,"(a)") '******************************************************'
       write(*,"(a)") 'Starting time step computation with PETSC/SLEPc'
       write(*,"(a)") '******************************************************'
    endif

    if ((.not.slepc_restartable()).and.(n_procs_y.gt.1).and.&
         &slepc_initialized_from_outside) then
       if (mype.eq.0) then
          write(*,"(a)") 'ERROR: the installed slepc version does not support'
          write(*,"(A)") '       repeated initialization as is required for '
          write(*,"(A)") '       n_procs_y > 1.'
          write(*,"(A)") '       Either install slepc > 3.1-p4 (recommended),'
          write(*,"(a)") '       use n_procs_y=1 or precompute dt_max.'
       endif
       stop
    endif

    call get_systime(time_1)
    !save present values
    ky0_ind_in=ky0_ind
    dt_max_in = dt_max
    n_procs_y_in=n_procs_y
    n_procs_sim_in=n_procs_sim

    !always use linear version of CalFullRhs for eigenvalue computation and neglect constant
    !f0 contribution
    rhs_nl=.false.
    rhs_f0=.false.
    sav_antenna_contrib=antenna_contrib
    antenna_contrib=.false.
    rhs_only_coll=.false.
    num_ev = 0

    !for y_local, the linear operator is block diagonal in the ky coordinate;
    !the time step limit for the whole operator is therefore the minimum of the
    !time step limits of the individual blocks (which can be computed much faster 
    !than the full system). 
    !Optimization and parallelization is not altered EXCEPT the n_procs_y,
    !which means that only 1/n_procs_y processors are used for the computation,
    !the rest idles
    if(y_local.and.(nky0.gt.1)) then
       n_procs_sim=n_procs_sim_in/n_procs_y_in
       n_procs_y=1
       nky0=1
       
       mpi_comm_world_save = MY_MPI_COMM_WORLD
       CALL set_my_mpi_comm_world(mpi_comm_world_save)

       if(mype.lt.n_procs_sim) then
          ! to further speed up the search for eigenvalues, we only check
          ! a maximum of 3 blocks. ky=ky0_ind, ky=ky0_ind+1 and ky=kymax.
          ! If nky0_in<=3, we check all ky modes given in the system
          
          if (nky0_in.gt.3) then
             blocknum=3
             allocate(kyind_for_ev(blocknum))
             do j=1,2
                kyind_for_ev(j)=ky0_ind_in+j-1
             enddo
             kyind_for_ev(3)=ky0_ind_in+nky0_in-1
          else
             blocknum=nky0_in
             allocate(kyind_for_ev(blocknum)) 
             do j=1,nky0_in
                kyind_for_ev(j)=ky0_ind_in+j-1
             enddo
          endif

          if (write_pe) write(*,"(a,i3,a)") "selecting ", blocknum,"  ky modes"
          
          call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)

          call compute_dt_vlasov(dt_vlasov,blocknum,kyind_for_ev,num_ev)

          if (arakawa_zv.and.hyp_on_h.and.hypz_opt) &
               call compute_dt_hypz(dt_vlasov,blocknum,kyind_for_ev,num_ev)

          if (coll_split) &
               call compute_dt_coll(dt_vlasov,dt_coll,blocknum,kyind_for_ev,num_ev)
          
          if (slepc_restartable()) &
               call my_SlepcFinalize(.false.)
       endif

       !restore the values from the parameters file
       n_procs_sim=n_procs_sim_in

       call set_my_mpi_comm_world(mpi_comm_world_save)

       nky0=nky0_in
       ky0_ind=ky0_ind_in
       n_procs_y=n_procs_y_in

    else !y_local=.false. .or. nky0.eq.1

       if(mype.lt.n_procs_sim) then
          
          call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)
          
          blocknum=1
          allocate(kyind_for_ev(blocknum))
          kyind_for_ev=ky0_ind

          call compute_dt_vlasov(dt_vlasov,blocknum,kyind_for_ev,num_ev)

          if (arakawa_zv.and.hyp_on_h.and.hypz_opt) &
               call compute_dt_hypz(dt_vlasov,blocknum,kyind_for_ev,num_ev)

          if (coll_split) &
               call compute_dt_coll(dt_vlasov,dt_coll,blocknum,kyind_for_ev,num_ev)

          if (slepc_restartable()) &
               call my_SlepcFinalize(.false.)
       end if
    end if

    in_timestep_comp = .false.

    !broadcasting the results: 
    call mpi_bcast(dt_max,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)
    call mpi_bcast(dt_vlasov,1,MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)

    if (arakawa_zv.and.hyp_on_h.and.hypz_opt) &
       &call mpi_bcast(hypz_compensation,1,MPI_LOGICAL,0,MY_MPI_COMM_WORLD,ierr)

    if (coll_split) then
       !the chosen RKCa coll_split_scheme (found on my_pey=0) is broadcasted
       call mpi_bcast(rkstages_coll,1,MPI_INTEGER,0,MY_MPI_COMM_WORLD,ierr)
       call mpi_bcast(dt_coll,size(dt_coll),MPI_REAL_TYPE,0,MY_MPI_COMM_WORLD,ierr)
       if (coll_split_scheme=='RKCa') then 
          call switch_coll_scheme(rkstages_coll)
       endif
    endif

    call get_systime(time_2)
    time_ev=time_2-time_1
    
    print_ini_msg = sav_print_ini_msg
    if(write_pe) then
       write(*,"(a,es12.3,a)") 'dt_max = ',dt_max ,' found with SLEPc'
       write(*,"(a,i4,a)") 'calculated ',num_ev,' eigenvectors'
       write(*,"(a,f10.3,a)") 'time for computation of dt_max:',time_ev,' sec'
       write(*,"(a)") '******************************************************'
    endif

    call my_barrier()
    antenna_contrib=sav_antenna_contrib

    PERFOFF
  end subroutine compute_dt_max

  subroutine compute_dt_vlasov(dt_vlasov,blocknum,kyind_for_ev,num_ev_out)
    implicit none
    complex,dimension(:),allocatable:: eigenvalues
    real,intent(inout)::dt_vlasov
    integer,intent(in)::blocknum
    integer,dimension(blocknum),intent(in)::kyind_for_ev
    integer,intent(inout):: num_ev_out
    integer:: num_ev
    integer::j
    
   ! if (allocated(eigenvalues)) deallocate(eigenvalues)
    num_ev=0
    dt_vlasov = -1
       
    do j=1,blocknum
       ky0_ind=kyind_for_ev(j)  
       call compute_relevant_evs(eigenvalues,num_ev)
    end do

    call dt_from_ev(eigenvalues,num_ev,dt_vlasov,.false.)

    if (dt_vlasov.gt.0) then
       dt_max = dt_vlasov
       if (write_pe) write(*,"(a,es10.3)") 'time step limit set to ',dt_max
    else
       if (mype==0)  write(*,"(a)") ' problem in timestep computation: dt_max not adapted!!'
    endif

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    num_ev_out=num_ev_out+num_ev

  end subroutine compute_dt_vlasov

  subroutine compute_dt_hypz(dt_vlasov,blocknum,kyind_for_ev,num_ev_out)
    !when z hyperdiffusion is aplied on h, it is also applied on phi
    !at small ky, this can sometimes lead to a large negative real eigenvalue
    !hypz_compensation substracts the hyp_on_phi term, to possibly allow a larger timestep
    !this routine searches this extreme smallest_real eigenvalue, determines dt restrictions
    !and switches off hypz_compensation when it is not needed
    implicit none
    complex,dimension(:),allocatable:: eigenvalues
    real,intent(inout)::dt_vlasov
    integer,intent(in)::blocknum
    integer,intent(inout):: num_ev_out
    integer,dimension(blocknum),intent(in)::kyind_for_ev
    character(len=20):: which_ev_in
    real::dt_hyp
    integer::j
    integer:: num_ev

    dt_vlasov = dt_max  !stores previous value (with hypz_compensation)
    dt_hyp = -1

    num_ev=0
    which_ev_in = which_ev

    !switching off hypz_compensation
    hypz_compensation=.false.
    !switch ev solver to smallest real (this is where the critical eigenvalue is)
    which_ev = 'smallest_real'
   
    if (write_pe) write(*,"(a)") 'testing time step limit without hypz_compensation:'
    do j=1,blocknum
       ky0_ind=kyind_for_ev(j)          
       call compute_relevant_evs(eigenvalues,num_ev)
    end do

    which_ev = which_ev_in
    call dt_from_ev(eigenvalues,num_ev,dt_hyp)
    if (dt_hyp.gt.0) then
       if (write_pe) write(*,"(a,es10.3)") '  small real eigenvalue found. time step limit: ',dt_hyp
       !if a negative real eigenvalue is present: is its timestep smaller?
       dt_max = min(dt_vlasov,dt_hyp)
    else
       if (write_pe)  write(*,"(a,i1)") '  no additional restriction due to small real eigenvalue found.'
    endif
    
    if (dt_vlasov.gt.1.05*dt_max) then 
       !restoring previous setting with hypz_compensation
       hypz_compensation=.true.
       dt_max=dt_vlasov
       if (write_pe) write(*,"(A)") '  using hyp_z compensation for larger timestep.'
    else
       !large negative ev is not found: it is safe to swith off hyp_z_compensation
       if (write_pe) write(*,"(A)") '  not using hyp_z compensation.'
       !we adapt the possibly new value for dt_vlasov
       dt_vlasov=dt_max
    endif
    if (allocated(eigenvalues)) deallocate(eigenvalues)
    num_ev_out=num_ev_out+num_ev

  end subroutine compute_dt_hypz



  subroutine compute_dt_coll(dt_vlasov,dt_coll,blocknum,kyind_for_ev,num_ev_out)
    implicit none
    complex,dimension(:),allocatable:: eigenvalues
    real,intent(inout)::dt_vlasov
    real,dimension(max_rkstages_coll),intent(inout)::dt_coll
    integer,intent(in)::blocknum
    integer,intent(inout)::num_ev_out
    integer,dimension(blocknum),intent(in)::kyind_for_ev
    integer::s, num_ev, ndigits

    dt_coll   = -1

    num_ev=0
    if (write_pe) write(*,"(A)") 'testing time step limit for collision scheme:'
    rhs_only_coll = .true.
    !take maximum ky0_ind (assuming that this is relevant, if spatial diffusion is on)
    ky0_ind=kyind_for_ev(blocknum)
    call compute_relevant_evs(eigenvalues,num_ev)
    rhs_only_coll = .false.

    ev_coll = maxval(abs(real(eigenvalues)))

    !round to available precision + 1 digit
    ndigits= int(log10(1.0/ev_prec))+1
    ev_coll=round_dt(ev_coll,ndigits)

    call my_real_max_to_all(ev_coll)

    if (coll_split_scheme=='RKCa') then
       !theoretical limits for eps=0.05 damped shifted Chebychev Polynomials
       !dt_coll(2) = 3.90419518*dt_coll(1)
       !dt_coll(3) = 8.74400189*dt_coll(1)
       !dt_coll(4) = 15.6167807*dt_coll(1)

       !compute the dt_coll array
       do rkstages_coll=1,max_rkstages_coll
          call dt_from_ev(eigenvalues,num_ev,dt_coll(rkstages_coll),.true.)
          if (dt_coll(rkstages_coll).gt.0) then
             if (write_pe) write(*,"(a,i1,a,es10.3)") '  RKC',rkstages_coll,' : time step limit is ',dt_coll(rkstages_coll)
          else
             if (mype==0)  write(*,"(a,i1)") '  could not determine time step for scheme RKC',rkstages_coll
          endif
       enddo

       !find appropriate rkstages_coll (default:1)
       !goal is to minimize stages and obtain a larger timestep than dt_vlasov/1.2
       !this 1.2 tries to account for the extra numerical effort
       rkstages_coll=1
       if (dt_coll(1).gt.0.and.dt_vlasov.gt.0) then
          rkstages_coll=max_rkstages_coll
          do s=max_rkstages_coll,1,-1
             if ((dt_coll(s).gt.dt_vlasov/1.2).and.(dt_coll(s).gt.0))&
                  &rkstages_coll=s
          enddo
       endif
       if (write_pe) write(*,"(a,i1)") '  adaptive coll_split_scheme: we chose RKC',rkstages_coll
    else
       call dt_from_ev(eigenvalues,num_ev,dt_coll(1),.true.)
       dt_coll = dt_coll(1)
       if (dt_coll(1).gt.0) then
          if (write_pe) write(*,"(a,i1,a,es10.3)") &
                '  RKC',rkstages_coll,' : time step limit is ',dt_coll(1)
       else
          if (mype==0)  write(*,"(a)") '  could not determine timestep for collision scheme!'
       endif
    endif

    !decide on total timestep, if possible
    if (dt_vlasov.gt.0.and.dt_coll(rkstages_coll).gt.0) dt_max = min(dt_vlasov,dt_coll(rkstages_coll))
    if (dt_vlasov.lt.0.and.dt_coll(rkstages_coll).gt.0) dt_max = dt_coll(rkstages_coll)
    if (dt_vlasov.gt.0.and.dt_coll(rkstages_coll).lt.0) dt_max = dt_vlasov
    !else (both .lt.0): do not change dt_max 

    if (allocated(eigenvalues)) deallocate(eigenvalues)
    num_ev_out=num_ev_out+num_ev

  end subroutine compute_dt_coll

  subroutine compute_relevant_evs(eigenvalues,num_ev)
    implicit none
    complex,dimension(:),allocatable:: eigenvalues
    integer:: num_ev
    integer:: iev, new_evs
    complex,dimension(:),allocatable:: old_eigenvalues
    call initialize_current_parall
    call initialize_calc_aux_fields
    
    call initialize_petsc_mat(.false.)
    call initialize_slepc_eps(.false.)
    call EPSSetFromOptions(eps_obj,globerr)
    PERFON('evsolver')
    call EPSSolve(eps_obj,globerr)
    PERFOFF
    
    call EPSGetConverged(eps_obj,new_evs,globerr)

    if (new_evs.eq.0) then
       if (num_ev.eq.0) then 
          !for some compilers we need to allocate a zero-size array here.
          if (.not.allocated(eigenvalues)) allocate(eigenvalues(0:-1))
       endif
       if (y_local) then 
          if (mype.eq.0) write(*,"(a,G12.4,a)") '   ***no eigenvalues computed for ky=',ky,'!!***'
       else
          if (mype.eq.0) write(*,"(a)") '  ***no eigenvalues computed !!***'
       endif
       call EPSGetIterationNumber(eps_obj,it_ev,globerr)
       if (mype.eq.0) write(*,"(a,i6)") '  number of ev iterations:',it_ev
    else
       if(num_ev.eq.0) then
          if (allocated(eigenvalues)) deallocate(eigenvalues)
          allocate(eigenvalues(0:new_evs-1))
       else
          allocate(old_eigenvalues(0:num_ev-1))
          old_eigenvalues=eigenvalues
          deallocate(eigenvalues)
          allocate(eigenvalues(0:num_ev+new_evs-1))
          eigenvalues(0:num_ev-1)=old_eigenvalues
          deallocate(old_eigenvalues)
       endif
       do iev=0,new_evs-1
          call EPSGetEigenpair(eps_obj,iev,eigenvalues(iev+num_ev),&
               &PETSC_NULL_OBJECT,glob_g,PETSC_NULL_OBJECT,globerr)
       enddo
       num_ev=num_ev+new_evs
       call EPSGetIterationNumber(eps_obj,it_ev,globerr)
       !if (write_pe) write(*,"(a,i6)") ' number of ev iterations:',it_ev
    endif
    
    call finalize_slepc_eps
    call finalize_petsc_mat
             
    call finalize_current_parall
    call finalize_calc_aux_fields

  end subroutine compute_relevant_evs

#endif  

  !>Calculates the dt_max for all eigenvalues in ev for the chosen time stepping
  !!scheme writes the minimum to the global variable dt_max
  subroutine dt_from_ev(ev,num_ev,dt_out,for_coll_scheme_in)
    implicit none
    integer,intent(in):: num_ev !<number of eigenvectors that is passed
    real,intent(inout):: dt_out
    complex,dimension(0:num_ev-1),intent(in):: ev !<array of the complex eigenvalues
    logical,intent(in),optional::for_coll_scheme_in
    logical::for_coll_scheme
    real:: dt_final, dt_ev
    integer:: nr

    if (present(for_coll_scheme_in))then
       for_coll_scheme=for_coll_scheme_in
    else
       for_coll_scheme=.false.
    endif

    !set initial value
    dt_final=100.
    
    !compute dt_max for each eigenvalue and find the smallest
    do nr=0,num_ev-1
       call compute_dt_from_ev(ev(nr),dt_ev,for_coll_scheme)
       if (dt_ev.ne.-1.) dt_final=min(dt_final,dt_ev)
    enddo

    !round to 3 significant digits
    dt_final=round_dt(0.99*dt_final,3)
    if (dt_final.ne.99.) then 
       dt_out=dt_final
       if (write_pe .and. which_ev.ne.'none') write(*,"(a)") '**********************************'
    !else: do not adapt dt!
    endif

  end subroutine dt_from_ev

  !>Rounds dt to ndigits significant digits (always decreases, never increases 
  !!dt to maintain stability)
  real function round_dt(dt_ev,ndigits)
    real:: dt_ev !<value to be rounded
    integer::ndigits !<significant digits
    real:: dt_aux
    integer::digits

    digits = 0

    dt_aux=dt_ev
    do while (dt_aux.lt.10**(ndigits-1))
       dt_aux=dt_aux*10
       digits=digits+1
    end do
    round_dt=real(floor(dt_aux))/10**digits

  end function round_dt

  !!Estimates without actual computation

  !>Prints the estimates for the time step limits in the various directions used
  !!by estimate_dt_max
  subroutine print_estimates
    real:: corrfac
    integer:: s

    if (write_pe) then
       write(*,*)
       write(*,"(A)") &
            &"=========== estimates for time step limit ==========="
       select case (parscheme)
       case('c4th')
          corrfac=1.4
       case('c6th')
          corrfac=1.6
       case default
          corrfac=1.0
       end select
       write(*,"(2(a,es10.3))") "CFL for parallel dynamics: dt<",&
            & dz/corrfac/maxclimit(3)*rk_corr(1)
       if (maxclimit(4).gt.0.0) then
          write(*,"(A,es10.3)")  "CFL for parallel velocity: dt<",&
               dv/1.4/maxclimit(4)*rk_corr(1)
       endif
       write(*,*)
       
       if ((maxclimit(1)>0).and.(maxval(kx).gt.0)) &
            & write(*,"(a,es10.3)") "linear x advection limit:  dt<",&
            & rk_corr(1)/(maxclimit(1)*maxval(kx))

       if ((maxclimit(2)>0).and.(kjmax.ne.0)) &
            & write(*,"(a,es10.3)") "linear y advection limit:  dt<",&
            & rk_corr(1)/(maxclimit(2)*kjmax)

       if(hyp_x.gt.0.) then
          if (yx_order) then
             write(*,"(a,es10.3)") "radial diffusion limit:    dt<",&
                  &rk_corr(2)/(hyp_x*djdiff_max)
          else
             write(*,"(a,es10.3)") "radial diffusion limit:    dt<",&
                  &rk_corr(2)/(hyp_x*didiff_max)
          endif
       endif
       
       if(hyp_y.gt.0.) then
          if (yx_order) then
             write(*,"(a,es10.3)") "binormal diffusion limit:  dt<",&
                  &rk_corr(2)/(hyp_y*didiff_max)
          else
             write(*,"(a,es10.3)") "binormal diffusion limit:  dt<",&
                  &rk_corr(2)/(hyp_y*djdiff_max)
          endif
       endif
       
       if(hyp_perp.gt.0.) then
          write(*,"(a,es10.3)") "perpend. diffusion limit:  dt<", &
               &rk_corr(2)/(hyp_perp*dperpdiff_max)
       endif
       
       if ((hyp_z_order.eq.4).and.(hyp_z.gt.0.)) &
            & write(*,"(a,es10.3)") "parallel diffusion limit:  dt<",&
            &rk_corr(2)/abs(hyp_z)

       if (hyp_v.gt.0.) &
            & write(*,"(a,es10.3)") "v_par diffusion limit:     dt<",&
            &rk_corr(2)/abs(hyp_v)

       if (collision_op.ne.'none') then
          if (coll_split_scheme=='RKCa')then 
             write(*,*)
             write(*,"(a,es10.3)") 'max. collisional eigenvalue = ',ev_coll_est
             write(*,"(a)") '!!WARNING: rough estimate without SLEPc!!'
             do s=1,max_rkstages_coll
                write(*,"(a,i1,a,es10.3)") 'coll. v-diffusion, RKC',s,':      dt<',dt_coll(s)
             enddo
          else
             write(*,"(a,es16.6,a)") 'coll. v-diffusion     dt<',dt_coll(1)
          endif
       endif

    endif
 
  end subroutine print_estimates

  !>Rough estimate of the maximum time step, the different advection terms are analyzed independently
  !!and the minimum is returned
  !!
  !!maxclimit must have been called befors so that maxclimitx/maxclimity are set
  !!\todo Routine doesn't take into account the stability regions of the different timeschemes
  !!\todo Adapt to nonlocal code
  subroutine estimate_dt_max(withoutput,edt,dtv,evc)
    logical,intent(in):: withoutput !<switch whether estimates are printed to standard output
    real, intent(out):: edt !<estimate dt_max
    real, intent(out):: dtv !<estimate dt_vlasov
    real, intent(out):: evc !<estimate ev_coll
    real::    cour_min=100.
    logical:: finalize_on_exit
    real:: vxmax, vymax, vzmax, vvmax

    write_pe=((mype==0).and.(print_ini_msg))

    !initialize RK scheme to get rk_corr
    call initialize_RK_coefficients

    !initialize maxclimitx/maxclimity if it hasn't been done already in perf_opt
    finalize_on_exit=.false.
    if (.not.allocated(kx)) then
       call initialize_current_parall
       finalize_on_exit=.true.
    endif


    !compute the advection limits in the different directions
    vxmax=0.
    vymax=0.
    vzmax=0.
    vvmax=0.
    
    !x
    IF (maxclimit(1).gt.0) then
       if (xy_local) then
          if (maxval(ki).gt.0) &
               vxmax=(maxclimit(1)*maxval(ki))
       else
          !vxmax = 1.4/deli*maxclimit(1)
          vxmax =  max_ddi*maxclimit(1)
       endif
    endif

    !y
    if ((maxclimit(2).gt.0).and.(kjmax.ne.0)) then
       vymax = maxclimit(2)*kjmax
    endif

    !z 
    !additional factor 1.4/1.6 for 4th/6th order finite differencing schemes
    !has been found experimentally but should be somewhere in the literature/
    !can probably derived after some thinking..
    select case (parscheme)
    case('c4th')
       vzmax=1.4/dz*maxclimit(3)
    case('c6th')
       vzmax=1.6/dz*maxclimit(3)
    case default
       vzmax=1./dz*maxclimit(3)
    end select
    
    !v_parallel
    if (maxclimit(4).gt.0.0) vvmax=1.4/dv*maxclimit(4)
    
    !in case vmax is zero, estimate dt ~ 1e-4
    if ((vxmax+vymax+vzmax+vvmax) .le. epsilon(vxmax)) vxmax=10000.0

    cour_min = rk_corr(1)/(vxmax+vymax+vzmax+vvmax)

    !diffusive limits
    if(hyp_x.gt.0.) then
       if (yx_order) then
          cour_min = min(cour_min, rk_corr(2)/(hyp_x*djdiff_max))
       else
          cour_min = min(cour_min, rk_corr(2)/(hyp_x*didiff_max))
       endif
    endif

    if(hyp_y.gt.0.) then
       if (yx_order) then
          cour_min = min(cour_min, rk_corr(2)/(hyp_y*didiff_max))
       else
          cour_min = min(cour_min, rk_corr(2)/(hyp_y*djdiff_max))
       endif
    endif
    
    if(hyp_perp.gt.0.) then
       cour_min = min(cour_min, rk_corr(2)/(hyp_perp*dperpdiff_max))
    endif     

    if ((hyp_z_order.eq.4).and.(hyp_z.gt.0.)) then
       cour_min = min(cour_min, rk_corr(2)/abs(hyp_z))
    endif

    if (hyp_v.gt.0.) then
       cour_min = min(cour_min, rk_corr(2)/abs(hyp_v))
    endif
    

    Call my_real_min_to_all(cour_min)
    
    if (collision_op .ne. 'none') then
       !set the timestep output parameters
       dtv = cour_min
       evc = ev_coll_est
       Call get_dt_coll(evc,dtv)
       edt = min(dtv, dt_coll(rkstages_coll))
    else
       edt = cour_min
       evc = -1
       dtv = edt
    endif

    if (withoutput) then
       call print_estimates      
       if (mype.eq.0) then
          if (coll_split) then
             write(*,"(a,es10.3)") 'set collisionless time step limit to ',dtv
             write(*,"(a,es10.3)") 'set collision time step limit to ',dt_coll(rkstages_coll)
             write(*,"(a,es10.3)") 'set total time step limit to ',edt
          else
             write(*,"(a,es10.3)") 'set time step limit to ',edt
          endif
       endif
    endif

    if(finalize_on_exit) call finalize_current_parall

  end subroutine estimate_dt_max

  !this routine sets dt_coll from a given maximum collisional eigenvalue.
  !for the RKCa coll_split_scheme, also the number of RKC stages is determined
  subroutine get_dt_coll(ev_coll,dtv)
    real,intent(in)::ev_coll,dtv
    complex,dimension(1):: evc
    real:: dt_vl
    integer:: s

    evc=cmplx(-ev_coll)

    !a very generous amplification accounts for nonlinear eigenvalue shifts.
    !this affects the initial choice of the number of RKC stages (in performance optimization e.g.)
    !during a nonlinear run, the number of RKC stages is adapted using more precise information
    if (nonlinear) then
       dt_vl=dtv*5
    else
       dt_vl=dtv
    endif

    if (coll_split_scheme=='RKCa') then
       !compute the dt_coll array
       do rkstages_coll=1,max_rkstages_coll
          call dt_from_ev(evc,1,dt_coll(rkstages_coll),.true.)
       enddo
       rkstages_coll=1
       !determine rkstages_coll
       if (dt_coll(1).gt.0.and.(dt_vl.gt.0)) then
          rkstages_coll=max_rkstages_coll
          do s=max_rkstages_coll,1,-1
             if ((dt_coll(s).gt.dt_vl/1.2).and.(dt_coll(s).gt.0))&
                  &rkstages_coll=s
          enddo
       endif
       if (write_pe) write(*,"(a,i1)") '  adaptive coll_split_scheme: we chose RKC',rkstages_coll
    else
       !use standard timescheme for collisional EV
       call initialize_RK_coefficients
       dt_coll = rk_corr(2)/ev_coll
    endif
  end subroutine get_dt_coll

  subroutine initialize_adapt_dt
    allocate(last_dts(0:win_dts-1),last_estimates(0:win_est-1))
    last_dts=dt
    last_estimates=dt
    ltav=0.

    !in perf_opt, one should not adapt the splitting scheme:
    !since dt_coll and dt_vlasov are only estimates, dt is set very small.
    coll_scheme_adaptive = ((coll_split_scheme=='RKCa').and.(.not.in_perf_opt))

    !if(mype.eq.0) open(TIMEFILE, file=trim(diagdir)//'/timeest.dat', form='formatted',&
    !        status='replace', position='rewind')
    !if(mype.eq.0) Write(TIMEFILE, "(A16,E16.8,A13,E16.8)") "#  dt_max=", dt_max, ", dt_vlasov=",dt_vlasov
    !if(mype.eq.0) Write(TIMEFILE, "(5A16)") "#     time", "nl_dt","dt_all","normvar","dt"

  end subroutine initialize_adapt_dt

  !!\todo the definition of nl_dt should be generalized for the global case
  subroutine adapt_dt
    real:: nl_dt, dt_inst, dt_est, dt_all, dt_GyroLES, nlev_max
    real, dimension(win_est):: sortdts
    real, dimension(win_est-2):: deri
    real:: av, var, normvar
    !threshold values for nonlinear adaption:
    real,parameter:: thresh_down=0.95,thresh_up=1.05

    !time step restriction due to maximized  chi x B term
    if (xy_local) then
       nlev_max = ve_max(1)*ki(hki)+ve_max(2)*kjmax
    else
       nlev_max = ve_max(1)*max_ddi+ve_max(2)*kjmax
    endif
    if (parallel_nl) &
         &nlev_max = nlev_max+1.4/dv*ve_pnl_max !~factor for 4th order centered differencing
    nl_dt=rk_corr(1)/nlev_max

    last_estimates(mod(itime-1,win_est)) = nl_dt

    !some statistical analysis to detect instabilities
    sortdts=cshift(last_estimates,mod(itime-1,win_est)-win_est+1)
    !2nd derivative
    deri=2*sortdts(2:win_est-1)-sortdts(3:win_est)-sortdts(1:win_est-2)
   
    !variance of the sample
    av=sum(deri)/(win_est-2)
    var=sum((deri-av)**2)/(win_est-3)

    !normalize variance to the old average
    if (((itime.le.win_est+1).or.(ltav.lt.(1000*epsilon(ltav)))).or. &
         &((abs(ExBrate).gt.epsilon(ExBrate)).and.(time.ge.ExB_stime))) then
       normvar=0. 
    else
       normvar=var/ltav
    end if

    !limit the reduction factor of the variance term to 0.3
    if (normvar.gt.50) normvar=50

    !update long term average
    if (itime.gt.win_est) then
       if(itime.gt.win_est+500) then
          ltav=0.998*ltav+0.002*var
       else
          ltav=((itime-win_est-1)*ltav+var)/(itime-win_est)
       end if
    end if
    
    !variance analysis: reduce dt when nl_dt decreases quickly (likely an instability),
    !courant>1 enables to compensate overestimation of chi x B effect
    dt_est = nl_dt*courant/cosh(normvar/4)**0.1
    !add linear eigenvalues and chi x B shift (worst case scenario)
    !1/dt_vlasov is a good estimate for eigenvalue due to dominating immag. part
    !without coll_split: dt_vlasov=dt_max is used
    !never exceed dt_max, which might be set by collisions for coll_split=T
    dt_all = min(1.0/(1.0/dt_vlasov+1.0/dt_est),dt_max)

    last_dts(mod(itime-1,win_dts))=dt_all

    dt_inst=minval(last_dts)

    if (dt_inst/dt.lt.thresh_down) then
       !ramp down only gradually 
       dt=dt_inst
       last_dts(mod(itime-1,win_dts))=dt
    elseif(dt_inst/dt.gt.thresh_up) then
       !ramp up dt only gradually
       dt=thresh_up*dt
       !modify last_dts to store this value for win_dts timesteps
       last_dts(mod(itime-1,win_dts))=dt
    endif

    !Modification to GyroLES 
    if (GyroLES) then
        dt_GyroLES = min(rk_corr(2)/(maxval(hyp_x_spec)*didiff_max),rk_corr(2)/(maxval(hyp_y_spec)*djdiff_max))
        dt=min(dt,dt_max,dt_GyroLES)
        !if(mype.eq.0) Write(TIMEFILE, "(6ES16.4)") time,nl_dt,dt_all,normvar,dt, dt_GyroLES
    else
        !if(mype.eq.0) Write(TIMEFILE, "(5ES16.4)") time,nl_dt,dt_all,normvar,dt
        dt=min(dt,dt_max) ! required in some cases to bridge the gap
                          ! between thresh_down and 1.0 if dt_inst=dt_max<dt
    endif
    
    if (coll_scheme_adaptive) &
         &call adapt_coll_scheme(dt_all,dt_est,thresh_up)


  end subroutine adapt_dt


  subroutine adapt_coll_scheme(dt_curr,dt_est,thresh_up)
    implicit none
    real,intent(in)::dt_curr,dt_est,thresh_up
    real::dt_s  !timestep for different number of stages
    real::dt_max_new
   
    !if current nonlinear dt is significantly (by 1.05*thresh_up) smaller than the maximum
    !allowed timestep of a lower stages coll_split_scheme,
    !we reduce the number of stages
    !to be safe, the threshold must be larger than the one for dt_inst to dt 
    if (rkstages_coll.gt.1) then
       if (dt_coll(rkstages_coll-1).gt.(1.05*thresh_up*dt_curr)) then
          rkstages_coll=rkstages_coll-1
          call switch_coll_scheme(rkstages_coll)
          dt_max=min(dt_vlasov,dt_coll(rkstages_coll))
          if (mype.eq.0) then
             write(*,"(a,i1,a,es10.3)") 'adapting coll_split_scheme to ', rkstages_coll,' stages at time=',time
             write(*,"(a,es10.3)") '   reset dt_max to ', dt_max
          endif
       endif
    endif

    !if nonlinear dt would be significantly (1.1*thresh_up times) larger with a higher stages coll_split_scheme,
    !we increase the number of stages (this should be rare..)
    if (rkstages_coll.lt.max_rkstages_coll) then
       !linear dt for higher scheme
       dt_max_new=min(dt_vlasov,dt_coll(rkstages_coll+1))

       !combine nonlinear, higher scheme linear and variance constraints
       dt_s=min(1.0/(1.0/dt_vlasov + 1.0/dt_est),dt_max_new)

       if (dt_s.gt.1.1*thresh_up*dt_curr) then
          rkstages_coll=rkstages_coll+1
          call switch_coll_scheme(rkstages_coll)
          dt_max=dt_max_new
          !dt is not adapted (this is left for next timestep)
          if (mype.eq.0) then
             write(*,"(a,i1,a,es10.3)") 'adapting coll_split_scheme to ', rkstages_coll,' stages at time =',time
             write(*,"(a,es10.3)") '   reset dt_max to ', dt_max
             write(*,"(a,es10.3,a,es10.3)") '   nonlinear estimate increases from dt=', dt_curr, ' to dt=',dt_s 
          endif
       endif
    endif

  end subroutine adapt_coll_scheme

  subroutine finalize_adapt_dt
    !if(mype.eq.0) close(TIMEFILE)
    deallocate(last_dts,last_estimates)
  end subroutine finalize_adapt_dt


end module compute_dt
