#include "redef.h"

!>Iterative eigenvalue solver based on SLEPc
!!
!!All input variables like the number of eigenvalues n_ev, specification of the method 
!!and spectral transform which_ev etc. can be set in the parameters file and/or can be 
!!overwritten at runtime with the appropriate SLEPc runtime options (see SLEPc manual) 
module eigen_iterative
  use par_mod, only: g_1, time, time_ev
  use eigen_parameters
  use par_in, only: diagdir, ev_out,file_extension, write_std, write_h5
  use discretization
  use file_io, only: get_unit_nr
  use communications, only: MY_MPI_COMM_WORLD, get_systime
  use rhs_computations
  use diagnostics
  use diagnostics_energy, only: get_energy_total
  use petsc_aux
  use slepc_aux
  use checkpoint
  use initcond
  use quasilinear_model
#ifdef COMBI
  use petsc_precond
#endif
#ifdef WITHFUTILS
  use futils
#endif

  implicit none
  public:: ev_iterative, n_ev, evfile
#ifdef COMBI
  public:: ev_operator_only, solve_for_residual
#endif
  private

  integer:: evfile
#ifdef COMBI
  complex:: local_shift
#endif

#include "petscversion.h"
contains
#ifdef COMBI
  subroutine ev_operator_only
    double precision :: time_1, time_2
    integer :: num_ev, i
    logical:: eof
    integer :: ev_handle
    real ::norm_energy
    logical:: dummy
!    character(len=10)::pc_string,ksp_string,eps_string

    PERFON('ev_it')
    if(mype.eq.0) then
       write(*,"(a)") '******************************************************'
       write(*,"(a)") 'Start to compute L*g only using PETSC/SLEPc'
       write(*,"(a)") '******************************************************'
    endif


    ! always use linear version of CalFullRhs for eigenvalue computation
    call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)

    !initialize_petsc_mat can't be used because derived rhs-computations
    !are allowed in this routine (e.g. scalapack interface)
    call MatCreateShell(MY_MPI_COMM_WORLD,vlen,vlen,vlen*n_procs_sim,&
         vlen*n_procs_sim,PETSC_NULL_INTEGER,shellmat,globerr)

    call initialize_calc_k
    call MatShellSetOperation(shellmat,MATOP_MULT,matmul_k,globerr)


    call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,glob_g,globerr)

    !initial space
    if (read_checkpoint) Then
        if(mype.eq.0) then
           write(*,"(a)") 'read checkpoint*********************************'
        endif
       call initialize_checkpoint_read
       eof=.false.
       i=0
       do while ((i.lt.5).and.(.not.eof))
          call checkpoint_read(g_1,eof)
          if(.not.eof) then
             i=i+1
             call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(i),globerr)
             call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0_out(i),globerr)
             call fo2pe(g_1,v0(i))
          end if
       end do
       n_ch=i
       call finalize_checkpoint_read
    else
       n_ch=1
       call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(1),globerr)
       call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0_out(1),globerr)
       call init_g_1
       call fo2pe(g_1,v0(1))
    end if


    call initialize_slepc_eps(.true.)


    !time measurement only for solver (avoiding file I/O)
    call get_systime(time_1)

!    call initialize_all_diags
    PERFON('ev_compute matvec')

!!! here one can insert the matrix vector multiplication
    do i = 1,n_ch
      call MatMult(shellmat,v0(i),v0_out(i),globerr)
    end do

    PERFOFF
!    call finalize_all_diags

    call get_systime(time_2)
    time_ev=time_2-time_1

   !!! here begin analyze

    call initialize_all_diags

    if(ev_out) call ev_out_open(ev_handle,1)

    do i=1,n_ch
       call pe2fo(v0_out(i),g_1)

       time=i

       !normalize eigenvectors so energy is 1.0
!       if (xy_local) then
!          call get_energy_total(g_1,norm_energy)
!          g_1=g_1/sqrt(norm_energy)
!       endif

       if(ev_out) call ev_out_write(g_1,i,n_ch,ev_handle)
       call exec_all_diags(0,time,dummy,dummy,dummy)
    enddo

    time=0.0

    !output to file and standard out
    if(write_std) then
       if(ev_out) call ev_out_close(ev_handle)
    end if

    call finalize_all_diags

   !!! here end analyze

    !Destroy initial vectors
    do i=1,n_ch
       call VecDestroy(v0(i),globerr)
       call VecDestroy(v0_out(i),globerr)
    end do

    call finalize_calc_k

    call finalize_slepc_eps
    call MatDestroy(shellmat,globerr)
    call VecDestroy(glob_g,globerr)
    call my_SlepcFinalize(.false.)

    if(mype.eq.0) then
       write (*,*) n_ch
       write(*,"(a,f10.3,a)") 'time for matrix-vector computation:',time_ev,' sec'
       write(*,"(a)") '******************************************************'
    endif
    PERFOFF

  end subroutine ev_operator_only

subroutine solve_for_residual
    double precision :: time_1, time_2
    integer :: num_ev, i
    logical:: eof
    integer :: ev_handle
    real ::norm_energy
    logical:: dummy


!    character(len=10)::pc_string,ksp_string,eps_string

    PERFON('ev_it')
    local_shift = ev_shift_ref
    if(mype.eq.0) then
       write(*,"(a)") '******************************************************'
       write(*,*) 'solve system for iterative refinement with lambda=', local_shift
       write(*,"(a)") '******************************************************'
    endif

    ! always use linear version of CalFullRhs for eigenvalue computation
    call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)

    !initialize_petsc_mat can't be used because derived rhs-computations
    !are allowed in this routine (e.g. scalapack interface)
    call MatCreateShell(MY_MPI_COMM_WORLD,vlen,vlen,vlen*n_procs_sim,&
         vlen*n_procs_sim,PETSC_NULL_INTEGER,shellmat,globerr)

    call initialize_calc_k
!    call MatShellSetOperation(shellmat,MATOP_MULT,matmul_k_shift,globerr)
    call MatShellSetOperation(shellmat,MATOP_MULT,matmul_k,globerr)


    call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,glob_g,globerr)


    !initial space
    if (read_checkpoint) Then
        if(mype.eq.0) then
           write(*,"(a)") 'read checkpoint*********************************'
        endif
       call initialize_checkpoint_read
       eof=.false.
       i=0
       do while ((i.lt.5).and.(.not.eof))
          call checkpoint_read(g_1,eof)
          if(.not.eof) then
             i=i+1
             call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(i),globerr)
             call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0_out(i),globerr)
             call fo2pe(g_1,v0(i))
             call fo2pe(g_1,v0_out(i))
             write (*,*) g_1(1,1,1,1,1,1)
          end if
       end do
       n_ch=i
       call finalize_checkpoint_read
    else
       n_ch=1
       call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(1),globerr)
       call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0_out(1),globerr)
       call init_g_1
       call fo2pe(g_1,v0(1))
    end if


    call initialize_slepc_eps(.true.)


    !time measurement only for solver (avoiding file I/O)
    call get_systime(time_1)

!    call initialize_all_diags
    PERFON('solve system')

!!! the solution of the system


    call KSPCreate(MY_MPI_COMM_WORLD,ksp_test_obj,globerr)
    call KSPSetFromOptions(ksp_test_obj,globerr)
    call initialize_L_g


    petsc_shift = -local_shift
    if(mype.eq.0) then
           write(*,"(a)") 'L_g is set up'
    endif
    call MatShift(L_g_mat, petsc_shift,globerr)
    call MatShift(shellmat,petsc_shift,globerr)

    call KSPSetOperators(ksp_test_obj,shellmat,L_g_mat,DIFFERENT_NONZERO_PATTERN,globerr)
    do i = 1,n_ch
        call KSPSolve(ksp_test_obj,v0(i),v0_out(i),globerr)
    end do

    call KSPDestroy(ksp_test_obj,globerr)

    PERFOFF
!    call finalize_all_diags

    call get_systime(time_2)
    time_ev=time_2-time_1

   !!! here begin analyze

    call initialize_all_diags

    if(ev_out) call ev_out_open(ev_handle,1)

    do i=1,n_ch
       call pe2fo(v0_out(i),g_1)

       time=i

       !normalize eigenvectors so energy is 1.0
!       if (xy_local) then
!          call get_energy_total(g_1,norm_energy)
!          g_1=g_1/sqrt(norm_energy)
!       endif

       if(ev_out) call ev_out_write(g_1,i,n_ch,ev_handle)
       call exec_all_diags(0,time,dummy,dummy,dummy)
    enddo

    time=0.0

    !output to file and standard out
    if(write_std) then
       if(ev_out) call ev_out_close(ev_handle)
    end if

    call finalize_L_g
    call finalize_all_diags

   !!! here end analyze

    !Destroy initial vectors
    do i=1,n_ch
       call VecDestroy(v0(i),globerr)
       call VecDestroy(v0_out(i),globerr)
    end do

    call finalize_calc_k

    call finalize_slepc_eps
    call MatDestroy(shellmat,globerr)
    call VecDestroy(glob_g,globerr)
    call my_SlepcFinalize(.false.)

    if(mype.eq.0) then
       write(*,"(a,f10.3,a)") 'time for solving_operation:',time_ev,' sec'
       write(*,"(a)") '******************************************************'
    endif
    PERFOFF

  end subroutine solve_for_residual

#endif

  !>Computes the desired eigenvalues
  subroutine ev_iterative
    
    double precision :: time_1, time_2
    integer :: num_ev, i
    logical:: eof
!    character(len=10)::pc_string,ksp_string,eps_string

    PERFON('ev_it')
    if(mype.eq.0) then
       write(*,"(a)") '******************************************************'
       write(*,"(a)") 'Starting eigenvalue computation with PETSC/SLEPc'
       write(*,"(a)") '******************************************************'
    endif

    
    ! always use linear version of CalFullRhs for eigenvalue computation
    call my_SlepcInitialize(MY_MPI_COMM_WORLD,.false.)

    !initialize_petsc_mat can't be used because derived rhs-computations
    !are allowed in this routine (e.g. scalapack interface)
    call MatCreateShell(MY_MPI_COMM_WORLD,vlen,vlen,vlen*n_procs_sim,&
         vlen*n_procs_sim,PETSC_NULL_INTEGER,shellmat,globerr)

    call initialize_calc_k
    call MatShellSetOperation(shellmat,MATOP_MULT,matmul_k,globerr)
    call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,glob_g,globerr)

    !initial space
    if (read_checkpoint) Then
       call initialize_checkpoint_read
       eof=.false.
       i=0
       do while ((i.lt.5).and.(.not.eof))
          call checkpoint_read(g_1,eof)
          if(.not.eof) then
             i=i+1
             call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(i),globerr)
             call fo2pe(g_1,v0(i))
          end if
       end do
       n_ch=i
       call finalize_checkpoint_read
    else
       n_ch=1
       call VecCreateMPI(MY_MPI_COMM_WORLD,vlen,vlen*n_procs_sim,v0(1),globerr)
       call init_g_1
       call fo2pe(g_1,v0(1))
    end if

    call initialize_slepc_eps(.true.)

    do i=1,n_ch
       call VecDestroy(v0(i),globerr)
    end do

    !time measurement only for solver (avoiding file I/O)
    call get_systime(time_1)

    if(which_ev=='jd_as'.or.which_ev=='ks_as') call initialize_all_diags 
    PERFON('evsolver')
    call EPSSolve(eps_obj,globerr)
    PERFOFF
    if(which_ev=='jd_as'.or.which_ev=='ks_as') call finalize_all_diags 

    call get_systime(time_2)
    time_ev=time_2-time_1

    call eps_info(num_ev)
    
    if (num_ev.gt.0) call analyze_eigenpairs(num_ev)
    call finalize_calc_k

    call finalize_slepc_eps
    call MatDestroy(shellmat,globerr)
    call VecDestroy(glob_g,globerr)
    call my_SlepcFinalize(.false.)

    if(mype.eq.0) then
       write(*,"(a,f10.3,a)") 'time for eigenvalue computation:',time_ev,' sec'
       write(*,"(a)") '******************************************************'
    endif
    PERFOFF
  end subroutine ev_iterative

  !>Analyzes the eigenpairs stored in eps_obj, i.e. writes the eigenvalues.dat file and
  !!applies the chosen diagnostics to the eigenvectors
  subroutine analyze_eigenpairs(num_ev)
    integer, intent(in):: num_ev !<number of eigenpairs in eps_obj
    complex,dimension(:),allocatable:: eigenvalue
    integer :: ev_handle
    real ::norm_energy
    integer:: i
    logical:: dummy

#ifdef WITHFUTILS
    integer :: fidev_h5
#endif

    allocate(eigenvalue(0:num_ev-1))
    call initialize_all_diags

    if(ev_out) call ev_out_open(ev_handle,1)
    do i=0,num_ev-1
       call EPSGetEigenpair(eps_obj,i,eigenvalue(i),PETSC_NULL_SCALAR,glob_g,&
            &PETSC_NULL_OBJECT,globerr)
       if(which_ev.eq.'shift_invert_s') then
          !backtransform by hand for shift_invert_s
          eigenvalue(i)=conjg(1./eigenvalue(i)+ev_shift)
       end if
       call pe2fo(glob_g,g_1) 
       
       time=i+1
       
       !normalize eigenvectors so energy is 1.0
       if (xy_local) then
          call get_energy_total(g_1,norm_energy)
          g_1=g_1/sqrt(norm_energy)
       endif

       if(ev_out) call ev_out_write(g_1,i+1,num_ev,ev_handle)
       call exec_all_diags(0,time,dummy,dummy,dummy)
       if ((i.eq.0).and.(quasilin_model.gt.0)) &
            & call set_quasilin(real(eigenvalue(i))) 
    enddo
    
    time=0.0
    
    !output to file and standard out
    if(write_std) then
       if(ev_out) call ev_out_close(ev_handle)

       if (mype.eq.0) then
          if(.not.cat_output) then
             call get_unit_nr(evfile)
             Open(evfile, file=Trim(diagdir)//'/eigenvalues'//trim(file_extension),&
                  &form='formatted', status='replace')
          end if
          write(*,"(a)") 'computed eigenvalues: '
          write(evfile, "(a,F7.3)") '#eigenvalues (re/im part) at ky=', ky(lg1)
          do i=0,num_ev-1
             write(*,"(2f15.8)") eigenvalue(i)
             write(evfile, "(8ES20.8)") eigenvalue(i)
          enddo
          if(.not.cat_output) close(evfile)
       endif
    end if


#ifdef WITHFUTILS
    if(write_h5) then
       if((mype).eq.0) then
          call creatf(trim(diagdir)//"/eigenvalues"//trim(file_extension)//".h5", &
               fidev_h5, "eigenvalues", 'd')
          call putarr(fidev_h5, "/eigenvalues", eigenvalue)
          call closef(fidev_h5)
       end if
    end if
#endif


    deallocate(eigenvalue)
    call finalize_all_diags
    
  end subroutine analyze_eigenpairs
  
  !>Defines the matrix-free matrix-vector multiplication by using calc_k
  !!
  !!This allows for the use of derived rhs-computations e.g. based on SCALAPACK
  !!(in contrast to imp_matmul provided in petsc_aux),
  !!which is important for the computation of left eigenvectors
  subroutine matmul_k(locmat,petsc_in,petsc_out,ierr)
#if !((PETSC_VERSION_MAJOR>2) && (PETSC_VERSION_MINOR>0))
#include "finclude/petsc.h"
#endif
#include "finclude/petscvec.h"   
#include "finclude/petscmat.h"
    Mat locmat
    Vec petsc_in,petsc_out
    PetscErrorCode ierr
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2):: fort_g,fort_rhs

    PERFON('rhsmult')

    call pe2fo(petsc_in,fort_g)
    Call calc_k(fort_g,fort_rhs,0)
    call fo2pe(fort_rhs,petsc_out)

    PERFOFF

  end subroutine matmul_k

#ifdef COMBI
    !>Defines the matrix-free matrix-vector multiplication by using calc_k
  !!
  !!This allows for the use of derived rhs-computations e.g. based on SCALAPACK
  !!(in contrast to imp_matmul provided in petsc_aux),
  !!which is important for the computation of left eigenvectors
  subroutine matmul_k_shift(locmat,petsc_in,petsc_out,ierr)
#if !((PETSC_VERSION_MAJOR>2) && (PETSC_VERSION_MINOR>0))
#include "finclude/petsc.h"
#endif
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
    Mat locmat
    Vec petsc_in,petsc_out
    PetscErrorCode ierr
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2):: fort_g,fort_rhs

    PERFON('rhsmult')

    call pe2fo(petsc_in,fort_g)
    Call calc_k(fort_g,fort_rhs,0)
    fort_rhs = fort_rhs - local_shift*fort_g

    call fo2pe(fort_rhs,petsc_out)

    PERFOFF

  end subroutine matmul_k_shift

#endif
  
end module eigen_iterative
