#include "redef.h"
#include "switches.h"
#include "intrinsic_sizes.h"
!>Contains routines to diagnose the eigenvalue spectrum in a nonlinear run
!!(linearizes the right hand side by taking a snapshot of the electromagnetic fields)
!/todo use this module in diag.F90 (reqires to remove circular dependency with slepc_aux.F90)
!/todo implementation for global runs (requires interface to global nonlinear term in calc_rhs.F90)
Module diag_nl_eigenvalues
  use calc_rhs, only: rhs_nl, rhs_f0, rhs_only_coll, rhs_nlev
  Use communications
  Use coordinates, only: kx, kjmax, dv
  use discretization, only: hkx
  use eigen_parameters
  Use file_io, only: get_unit_nr
  Use par_mod
  !Use phys_ini
  use RK_coefficients, only: rk_corr
  !use time_scheme
#ifdef WITHSLEPC
  use compute_dt, only: dt_from_ev
  use petsc_aux
  use slepc_aux
#endif

  Implicit None

  !The basic diag interface
  PUBLIC :: initialize_all_diags_nlev, exec_all_diags_nlev, &
       &finalize_all_diags_nlev,mem_est_all_diags_nlev, &
       &check_diag_nlev,&
       &istep_nlev, nlev_stime, n_nlev

  PRIVATE 

  integer:: istep_nlev=-1,n_nlev = 0
  real:: nlev_stime = -1
  integer:: NLEVFILE
  Character(Len=8):: filestat='replace', filepos='rewind'
  complex,dimension(:),allocatable::eigenvalues_lin
  complex::max_ev_lin

Contains

!!************** public subroutines *********************!!

  !> check input parameters for this module
  Subroutine check_diag_nlev

#ifndef WITHSLEPC
    if (istep_nlev.gt.0) then
       if (mype==0) write(*,'(a)') 'no SLEPC/PETSc linked: we set istep_nlev = 0'
       istep_nlev=0
    endif
#endif
    if ((istep_nlev.gt.0).and.(.not.nonlinear)) then
       if (mype==0) Write(*,'(A)') 'switching off diag_nlev (only possible for nonlinear runs)'
       istep_nlev=0
    endif
    
    if ((istep_nlev.gt.0).and.(.not.xy_local)) then
       if (mype==0) Write(*,'(A)') 'switching off diag_nlev (currently not implemented for global runs)'
       istep_nlev=0
    endif

    if (istep_nlev.gt.0) then
#ifndef with_extended_diags
       if (mype==0) Write(*,'(a)') 'ERROR: with_extended_diags must be set in switches.h for using diag_nlev'
       stop
#endif
       !setting the number of nonlinear eigenvalues
       if(n_nlev.le.0) then
          if (n_ev.gt.0) then 
             n_nlev = n_ev
          else
             n_nlev = 20
          endif
          if (mype==0) Write(*,'(a,i6,a)') 'using ',n_nlev,' eigenvalues for diag_nlev'
       end if
    end if

  End Subroutine check_diag_nlev


  !>Give an estimate of the memory requirements of this module
  !!(no big memory demand at the moment)
  Real Function mem_est_all_diags_nlev(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !if (istep_nlev.gt.0) mem_loc = mem_loc+0.0

    mem_est_all_diags_nlev=mem_req_in+mem_loc

  End Function mem_est_all_diags_nlev


  Subroutine initialize_all_diags_nlev

#ifdef with_extended_diags
    if (istep_nlev.gt.0) call initialize_diag_nlev
#endif

  End Subroutine initialize_all_diags_nlev

  Subroutine exec_all_diags_nlev(itime,time,first)
    integer, intent(in) :: itime
    real, intent(in) :: time
    logical,intent(in)::first
#ifdef WITHSLEPC
#ifdef with_extended_diags
    logical::rhs_nl_in
    
    if (istep_nlev.gt.0) then
       if(first) then
          !one set of eigenvalues for the linear case
          !the resulting timestep should coincide with dt_vlasov/dt_max from a ky subset
          rhs_nl_in =  rhs_nl
          rhs_nl =  .false.
          call diag_nlev(time,-1)
          rhs_nl=rhs_nl_in
       else
          if ((modulo(itime,istep_nlev) == 0).and.(time.ge.nlev_stime)) then
                call diag_nlev(time,0) !nlev
                !call diag_nlev(time,1) !model 1 : vemax
                call diag_nlev(time,5) !model 5 : lin+ve
          endif
       endif
    endif
#endif
#endif

  End Subroutine exec_all_diags_nlev

  Subroutine finalize_all_diags_nlev

#ifdef with_extended_diags
    if (istep_nlev.gt.0) call finalize_diag_nlev
    if (allocated(eigenvalues_lin)) deallocate(eigenvalues_lin)
#endif

  End Subroutine finalize_all_diags_nlev

!**************************internal routines***************************************

#ifdef with_extended_diags

  Subroutine initialize_diag_nlev
    Implicit None
    
    if (write_std) then
       if(mype==0) call get_unit_nr(NLEVFILE)
       if(mype==0) print*,NLEVFILE
       if(mype==0) open(NLEVFILE, file=trim(diagdir)//'/eigenvalues_nl'//&
            &''//trim(file_extension),&
            form='formatted', status=filestat, position=filepos)  
    endif

  End Subroutine initialize_diag_nlev

#ifdef WITHSLEPC
  Subroutine write_evfile(eigenvalues,num_ev,dt_nl,dt_ExB,time,what)
    implicit none
    integer,intent(in)::num_ev
    complex,dimension(0:num_ev-1),intent(in)::eigenvalues
    real,intent(in)::dt_nl, dt_ExB, time
    integer,intent(in):: what
    integer::i
    character(len=16)::whatstr

    !output to file
    if(write_std) then
       if (mype.eq.0) then
          select case(what)
          case(0)
             whatstr='nlev'
          case(1)
             whatstr='model vemax'
          case(5)
             whatstr='model lin+ve'
          endselect

          select case(what)
          case(-1)
             write(NLEVFILE, "(a,ES12.3,a,ES12.3,a)") &
                  &'#linear eigenvalues (re/im part),  (dt_max = ',dt_max,' , dt_vlasov = ',dt_vlasov,')'
          case(0,1,5)
             !             #...  time ()=( dt   ,  dt    , dt     ) what
             write(NLEVFILE, "(a,f10.3,a,ES10.3,a,ES10.3,a,ES10.3,a,a)") &
                &'#nonlinear eigenvalues (re/im part) at time ', time,&
                &' (dt_nl,dt_ExB,dt) = (',dt_nl,',',dt_ExB,',',dt,')',whatstr
          endselect
          do i=0,num_ev-1
             write(NLEVFILE, "(8ES20.8)") eigenvalues(i)
          enddo
          write(NLEVFILE, "(a)") NEW_LINE('(a)')
       end if
       call flush(NLEVFILE)
    end if

  End Subroutine write_evfile
#endif
!*************************************************

  !>computes and writes nonlinear eigenvalues
  !!this is basically the same as compute_dt, but using all ky modes
  !!Runge-Kutta timestep is computed(for comparison), but not adapted
  Subroutine diag_nlev(time,what)
    implicit none
    real,intent(in):: time
    integer,intent(in):: what
#ifdef WITHSLEPC
    double precision:: time_1, time_2
    complex,dimension(:),allocatable:: eigenvalues
    integer:: num_ev, j
    logical:: sav_print_ini_msg, sav_rhs_f0
    real:: dt_nl, dt_ExB
    real:: max_ev_nl

    PERFON('diag_nl_eigenvalues')

    if (what.ne.-1) then
!   from compute_dt: adapt_dt: nonlinear estimate from maximum ExB vel.
       if (parallel_nl) then
          dt_ExB=rk_corr(1)/(ve_max(1)*kx(hkx)+ve_max(2)*kjmax+ve_pnl_max/dv)
       else
          dt_ExB=rk_corr(1)/(ve_max(1)*kx(hkx)+ve_max(2)*kjmax)
       endif
    else
       dt_ExB = 0.0
    endif

    in_timestep_comp = .true.
    
    sav_print_ini_msg=print_ini_msg
    print_ini_msg = .false.
    !alway neglect constant f0 contribution to CalFullRhs for eigenvalue computation
    sav_rhs_f0 = rhs_f0
    rhs_f0=.false.
    

    call get_systime(time_1)

    select case(what)
    case(-1,0) !really compute ev's (total, only nl, or with model nl)
       num_ev=0
       if((mype.eq.0).and.print_ini_msg) then
          write(*,"(a)") '******************************************************'
          write(*,"(a)") 'Starting eigenvalue computation with PETSC/SLEPc'
          write(*,"(a)") '******************************************************'
       endif
       call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)
       call compute_nonlinear_evs(eigenvalues,num_ev)
       call my_SlepcFinalize(.false.)

       if (what==-1) then 
          allocate(eigenvalues_lin(size(eigenvalues)))
          eigenvalues_lin=eigenvalues
          !maximum imaginary part of linear eigenvalue (absolute value)
          !to combine with max. ExB drift for nonlinear dt estimate)
          max_ev_lin = maxval(abs(aimag(eigenvalues)))
       endif
    case(5) 
       !nl_model by adding some ve_max to linear ev's
       !currently, only what=5 is implemented
       if (allocated(eigenvalues)) deallocate(eigenvalues)
       num_ev = size(eigenvalues_lin)
       allocate(eigenvalues(1:num_ev))
       !max. (positive) ExB drift
       max_ev_nl = (ve_max(1)*kx(hkx)+ve_max(2)*kjmax)
       do j=1,num_ev
          !combine with linear ev (sign of linear ev is used!!)
          eigenvalues(j)=eigenvalues_lin(j)+imag*sign(max_ev_nl,aimag(eigenvalues_lin(j)))
       enddo
    case(1) 
       !just the nonlinear shift
       if (allocated(eigenvalues)) deallocate(eigenvalues)
       allocate(eigenvalues(1:1))
       eigenvalues(1) = imag*(ve_max(1)*kx(hkx)+ve_max(2)*kjmax)
    endselect

    call dt_from_ev(eigenvalues,num_ev,dt_nl)
    
    call write_evfile(eigenvalues,num_ev,dt_nl,dt_ExB,time,what)
    if (allocated(eigenvalues)) deallocate(eigenvalues)

    in_timestep_comp = .false.
    
    call get_systime(time_2)
    time_ev=time_2-time_1
    
    print_ini_msg = sav_print_ini_msg
    rhs_f0 = sav_rhs_f0
    if (mype.eq.0) write(*,"(a,f10.3,a)") 'Time for computation of nonlinear eigenvalues:',time_ev,' sec'

    PERFOFF
#endif
  End Subroutine diag_nlev

!*************************************************
#ifdef WITHSLEPC
  subroutine compute_nonlinear_evs(eigenvalues,num_ev)
    implicit none
    complex,dimension(:),allocatable:: eigenvalues
    integer:: num_ev
    integer:: iev
    
    !use the linearized version of the nonlinearity in CalFulRhs
    rhs_nlev = .true.

    call initialize_petsc_mat(.false.,for_nlev=.true.)
    call initialize_slepc_eps(.false.,n_nlev)
    call EPSSetFromOptions(eps_obj,globerr)
    !PERFON('nl_evsol')
    call EPSSolve(eps_obj,globerr)
    !PERFOFF
    
    call EPSGetConverged(eps_obj,num_ev,globerr)

    if (num_ev.eq.0) then
       if (mype.eq.0) write(*,"(a)") '***no eigenvalues computed !!***'
       call EPSGetIterationNumber(eps_obj,it_ev,globerr)
       if (mype.eq.0) write(*,"(a,i6)") 'number of iterations:',it_ev
    else
       allocate(eigenvalues(0:num_ev-1))
       do iev=0,num_ev-1
          call EPSGetEigenpair(eps_obj,iev,eigenvalues(iev),&
               &PETSC_NULL_OBJECT,glob_g,PETSC_NULL_OBJECT,globerr)
       enddo
    endif
    
    call finalize_slepc_eps
    call finalize_petsc_mat(.false.)

    rhs_nlev = .false.

  end subroutine compute_nonlinear_evs
#endif

  Subroutine finalize_diag_nlev

    if (write_std.and.mype.eq.0) close(NLEVFILE)

  End Subroutine finalize_diag_nlev
  
#endif
!******************************************************************!!!

End Module diag_nl_eigenvalues
