#include "redef.h"

!>Solver for neoclassical equilibrium.
!!
!!At present, only an iterative solver based on PETSc is implemented; this could be extended by a direct
!!solver based on scalapack, similar to eigenvalue_comp
module neo_equilibrium
  use par_mod
  use discretization
  use f0_term  
  use phys_ini
  use aux_fields
  use diagnostics
  use calc_rhs, only: rhs_nl, rhs_f0
  use communications, only: MY_MPI_COMM_WORLD, get_systime, my_barrier
  use checkpoint

#ifdef WITHSLEPC 
  !slepc is not used, only the my_SlepcInitialize, which is used instead of PetscInitialize to
  !avoid problems in scans (should eventually be removed!!)
  use slepc_aux
  use petsc_aux
  use petsc_precond
#endif

  implicit none
  public:: initialize_g1_nc, finalize_g1_nc, comp_nc
  public:: g1_nc
  public:: nc_prec, nc_max_it, set_nc_defaults
  private

#ifdef WITHSLEPC 
#include "finclude/slepc.h"
#include "petscversion.h"
#if (PETSC_VERSION_MINOR<1)&&(PETSC_VERSION_MAJOR<4)
#define KSPLGMRES "lgmres"
#endif
#endif

  complex,dimension(:,:,:,:,:,:),allocatable:: g1_nc
  real:: nc_prec = 5.e-7
  integer:: nc_max_it=100000

contains

  !>Computes the neoclassical equilibrium, runs the diagnostics and cleans up afterwards.
  subroutine comp_nc
    call initialize_g1_nc
    call finalize_g1_nc    
  end subroutine comp_nc


  !>Sets the input parameter dafaults for neoclassical equilibrium computation
  subroutine set_nc_defaults
    nc_max_it = 100000
    nc_prec = 5.e-7
  end subroutine set_nc_defaults


  !>Initializes g1_nc, i.e. the converged g1(kx=0,ky=0)(t) due to the f0 contribution.
  !!
  !!The routine temporarily changes the resolution and/or parallelization of gene, 
  !!it therefore (re-) initializes and finalizes GENE and should therefore NOT be called when
  !!GENE is already initialized. 
  subroutine initialize_g1_nc
#ifdef WITHSLEPC 
    logical:: dummy, only_zonal_in, nc_not_converged
    integer:: nx0_in, nky0_in, ky0_ind_in, istep_neoclass_in, istep_schpt_in
    double precision :: time_1, time_2
    complex,dimension(:),allocatable:: rhs

    if (n_procs_y.gt.1) then
       print*,'The computation of the neoclassical equilibrium works only for n_procs_y=1'
       stop
    end if
    if (.not.y_local) then
       print*,'The computation of the neoclassical equilibrium only works for y_local=T'
       stop
    endif


    include_f0_contr=.true.

    if (xy_local) then
       nx0_in=nx0
       nx0=1
       evenx=0
    endif

    istep_neoclass_in=istep_neoclass
    istep_neoclass=1
    istep_schpt_in=istep_schpt
    istep_schpt=0
    
    nky0_in=nky0
    nky0=1
    ky0_ind_in=ky0_ind
    ky0_ind=0
    only_zonal_in=only_zonal
    only_zonal=.true.
    rhs_nl=.false.
    nc_not_converged = .false.

    final_init=.true.
    call initialize_current_parall
    call initialize_calc_aux_fields

    if(mype.eq.0) then
       write(*,"(a)") '***********************************************************'
       write(*,"(a)") 'Starting computation of neoclassical equilibrium with PETSc'
       write(*,"(a)") '***********************************************************'
    end if

    !time measurement
    call get_systime(time_1)

    allocate(g1_nc(pi1:pi2,pj1:pj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    
    !initialize PETSc, use wrapper to avoid errors at the end of scans
    PETSC_COMM_WORLD=MY_MPI_COMM_WORLD
    !call PetscInitialize(PETSC_NULL_CHARACTER,globerr)
    call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)

    rhs_f0=.true.
    !here, the rhs comp. (CalFullRhs) is given to PETSc
    !(this includes the solve for the fields, as well)
    call initialize_petsc_mat(.false.) 
    call initialize_petsc_vec
    call initialize_petsc_flags

    !Now, the f0_term is switched off in CalFullRhs
    !and set as the right-hand-side vector instead,
    !see KSPsolve below
    allocate(rhs(lijklmn0))
    rhs=0.
    call add_f0_term(rhs,1,lk0*ll0*lm0*ln0,0.)
    rhs=-rhs

    in_timestep_comp=.true.
    rhs_f0=.false.

    order_rcm=.true.
    call KSPCreate(MY_MPI_COMM_WORLD,ksp_obj,globerr) 
    call KSPSetType(ksp_obj,KSPLGMRES,globerr)
    call KSPSetTolerances(ksp_obj,1.e-30,nc_prec,300.,nc_max_it,globerr)
    call KSPGMRESSetRestart(ksp_obj,60,globerr)

    if(pc_type.eq.'none') then
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>4)       
       call KSPSetOperators(ksp_obj,shellmat,shellmat,globerr)
#else
       call KSPSetOperators(ksp_obj,shellmat,shellmat,DIFFERENT_NONZERO_PATTERN,globerr)
#endif
    else
       if(mype.eq.0) write(*,"(a)") 'computing preconditioner'
       ksp_prec=nc_prec
       call initialize_L_g 
#if (PETSC_VERSION_MAJOR>3) || (PETSC_VERSION_MINOR>4)
       call KSPSetOperators(ksp_obj,shellmat,L_g_mat,globerr)
#else
       call KSPSetOperators(ksp_obj,shellmat,L_g_mat,DIFFERENT_NONZERO_PATTERN,globerr)
#endif
       call KSPGetPC(ksp_obj,pc_obj,globerr)
       call set_pc_obj
       call my_barrier()
       call get_systime(time_2)
       time_nc=time_2-time_1
       if(mype.eq.0) write(*,"(a,f10.3,a)") 'time to compute the preconditioner:',time_nc,' sec'
    end if
    if (read_checkpoint.and.comp_type.eq.'NC') then
       !!tell PETSC to not zero out the initial guess glob_g passed to KSPSolve 
       call KSPSetInitialGuessNonzero(ksp_obj,petsc_t,globerr)
       call initialize_checkpoint_read
       call checkpoint_read(g1_nc)
       call finalize_checkpoint_read
       call fo2pe(g1_nc,glob_g)
    endif

    !read additional run time options
    call KSPSetFromOptions(ksp_obj,globerr)

    call fo2pe(rhs,glob_rhs)

!precond_approx=.true.
!call MatMult(shellmat,glob_rhs,glob_g,globerr)
!call pe2fo(glob_g,g1_nc)

!call MatMult(L_g_mat,glob_rhs,glob_g,globerr)
!call pe2fo(glob_g,g_1)
!print*,mype,sum(abs(g_1-g1_nc)),sum(abs(g_1))
!stop

    call KSPSolve(ksp_obj,glob_rhs,glob_g,globerr)

    call KSPGetConvergedReason(ksp_obj,kspreason,globerr)
    if (kspreason.lt.0) then
        if(mype.eq.0) print *, 'PETSc did not converge! No diagnostic output will be generated!'
        nc_not_converged = .true.
    endif
!    rhs_f0=.true.
!    call MatMult(shellmat,glob_g,glob_rhs,globerr)
!    call pe2fo(glob_rhs,g1_nc)
!    print*, mype,sum(abs(rhs)),sum(abs(g1_nc)),sum(abs(rhs))-sum(abs(g1_nc))

    call pe2fo(glob_g,g1_nc)

    call KSPDestroy(ksp_obj,globerr)
    if(pc_type.ne.'none') call finalize_L_g

    call get_systime(time_2)
    time_nc=time_2-time_1

    if(mype.eq.0) then
       write(*,"(a,f10.3,a)") 'time for computation of neoclassical equilibrium:',time_nc,' sec'
       write(*,"(a)") '******************************************************'
    endif
    
    deallocate(rhs)
    in_timestep_comp=.false.

    !output results
    call initialize_all_diags
    g_1=g1_nc
    if (.not. nc_not_converged) call exec_all_diags(0,time,dummy,dummy,dummy)
    call finalize_all_diags
    rhs_f0=.true.
    call finalize_petsc_vec
    call finalize_petsc_mat

    !call PetscFinalize(globerr)
    call my_SlepcFinalize(.false.)

    call finalize_current_parall
    call finalize_calc_aux_fields

    if (xy_local) then
       nx0=nx0_in
       if (mod(nx0,2).eq.0) then
          evenx=1
       else
          evenx=0
       endif
    end if

    istep_neoclass=istep_neoclass_in
    istep_schpt=istep_schpt_in

    nky0=nky0_in
    ky0_ind=ky0_ind_in
    only_zonal=only_zonal_in
#else
    print *, 'The neoclassical equilibrium can only be computed with PETSc at the moment.'
    print *, 'Please install PETSc/SLEPc and recompile GENE with the appropriate linker command.'
    stop
#endif

  end subroutine initialize_g1_nc
  
  !>Deletes g1_nc
  subroutine finalize_g1_nc
    deallocate(g1_nc)
  end subroutine finalize_g1_nc

end module neo_equilibrium
