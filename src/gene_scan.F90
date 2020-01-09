#include "redef.h"
#include "intrinsic_sizes.h"
module gene_scan
  use communications
  use parameters_IO
  use file_io, only: get_unit_nr
  use mpi
  use discretization
  use par_in, only: diagdir
  use par_other, only: par_in_dir
  use gene_subroutine
  use diagnostics, only: cat_output, nrgfile
  use checkpoint, only: checkpoint_in, read_checkpoint, chptdir
#ifdef WITHSLEPC
  !PETSc/SLEPc 3.0 or older cannot be initialized/finalized more than once,
  !therefore the initialization routines are called here and not in the initialization of
  !EV and dt computations.
  use slepc_aux
  use eigen_iterative, only: EVFILE
#endif

  implicit none
  
  public:: check_for_scan, run_scan
  public:: SCANFILE, LOGFILE, idle_global
  private

  integer, dimension(10):: scan_dims=0
  character(len=FILENAME_MAX):: chpt_in
  integer:: n_scan_dims, n_sims
  character(len=7):: form
  integer,dimension(:),allocatable:: p_stat,sr_tags,r_requests,s_buf
  logical:: continuation
  integer :: comm_group, comm_parall
  integer:: nextp,sbi
  integer:: SCANFILE, LOGFILE
  real(8):: idle_global
  integer:: ierr

  namelist /scan/ scan_dims, par_in_dir, chpt_in, n_procs_sim

contains

  subroutine check_for_scan(mult_par)
    logical, intent(out):: mult_par
    integer :: n_procs, PARAMFILE

    !read relevant parameters
    !scan defaults: only one simulation
    scan_dims(1)=1
    n_parallel_sims=1

    call read_diagdir

    call read_parall_nml(par_in_dir)

    call get_unit_nr(PARAMFILE)
    open(PARAMFILE, file='parameters', form='formatted')
    read(PARAMFILE, nml=scan,iostat=ierr)
    if (ierr.lt.0) then
       mult_par=.false.
    elseif (ierr.gt.0) then
       stop 'on i/o error: error in scan namelist'
    endif
    close(PARAMFILE)

    call mpi_comm_size (mpi_comm_world, n_procs, ierr)
    if (n_procs_sim.eq.-1) then
       if(modulo(n_procs,n_parallel_sims).eq.0) then
          n_procs_sim=n_procs/n_parallel_sims
       else
          stop 'n_procs is not a multiple of n_procs_sim' 
       end if
    end if

    !compute number of simulations and dimensionality of the scan
    n_scan_dims=0.
    n_sims=1
    do while (scan_dims(n_scan_dims+1).gt.0)
       n_scan_dims=n_scan_dims+1
       n_sims=n_sims*scan_dims(n_scan_dims)
    end do

#ifndef COMBI_MGR
    !check input
    if (n_procs.ne.n_parallel_sims*n_procs_sim) then
       write (*,"(A,I4,A,I4,A,I6,A)") 'ERROR: nr. of parallel sims. (',n_parallel_sims,&
            ') x nr. of procs per sim. (',n_procs_sim,&
            ') <> nr. (',n_procs,')!!'
       stop
    end if
#endif
    
  end subroutine check_for_scan


  subroutine run_scan(gene_comm, comm_parall_in)
    integer:: gene_comm, comm_parall_in
    real(8):: start_idle, idle_time
    character(len=20):: f_ext
    character(len=FILENAME_MAX):: ch_in
    character(len=7):: strnumber
    integer:: pnum, ierr

    comm_group  = gene_comm
    comm_parall = comm_parall_in

#if defined(WITHSLEPC)
    !initialize SLEPc (should be removed once everyone is using SLEPc 3.1)
    if (.not.slepc_restartable()) &
         call my_SlepcInitialize(comm_group,.true.)
#endif
     !initialize while loop parameters
    
     !determine file extension format
     pnum=1
     !check if parameters_0001 exists (scanscript interface)
     form='(i4.4)'
     if(.not.par_exists(pnum,form)) then 
        !check if parameters_1 exists (SG++ interface)
        form='(i7)'
     end if

     if (cat_output.and.(mype.eq.0)) then
        write(strnumber,'(i7)') my_sim
        f_ext='_'//trim(adjustl(strnumber))
        SCANFILE=65
        open(SCANFILE, file=trim(diagdir)//'/scan_out'//f_ext, &
             form='formatted', status='replace', position='rewind')
        !for the moment, only eigenvalues.dat, nrg.dat and parameters.dat are collected in one file
#ifdef WITHSLEPC
        EVFILE=SCANFILE
#endif
        NRGFILE=SCANFILE
        PAR_OUT_FILE=SCANFILE
     end if

     call init_next_p(nextp)

     write(strnumber,form) nextp
     f_ext='_'//trim(adjustl(strnumber))
     do while (par_exists(nextp,form))
        if (cat_output.and.(mype.eq.0)) write(SCANFILE,'(a)') 'parameters'//f_ext
        call get_closest_checkpoint(nextp,ch_in)
        call update_logfile(nextp,'s')
        call rungene(comm_group, par_in_dir, f_ext,ch_in)

        if (cat_output.and.(mype.eq.0)) call flush(SCANFILE)        

        !compute next file extension
        call get_next_p(nextp)
        write(strnumber,form) nextp
        f_ext='_'//trim(adjustl(strnumber))
     end do

     start_idle= MPI_Wtime()
     call mpi_barrier(MPI_COMM_WORLD, ierr)

     idle_time=MPI_Wtime()-start_idle
     call mpi_reduce(idle_time, idle_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0,&
          MPI_COMM_WORLD, ierr)

     if (cat_output.and.(mype.eq.0)) then
        close(SCANFILE)
     end if

     call finalize_logfile

#if defined(WITHSLEPC)
     if (.not.slepc_restartable()) &
          call my_SlepcFinalize(.true.)
#endif

   end subroutine run_scan
  

  subroutine get_closest_checkpoint(pnum,checkpoint)
    integer, intent(in):: pnum      
    character(len=FILENAME_MAX), intent(out):: checkpoint
    character(len=20):: checkpoint_extension  
    real,dimension(pnum-1):: dist
    logical,dimension(pnum-1):: d_mask
    
    character(len=7):: strnumber
    logical:: valid
    integer:: chpt_size=0
    integer:: oldruns, nearest, PARAMFILE, ierr
    
    !default: no checkpoint
    chpt_in='no'
    write(strnumber,form) pnum
    call get_unit_nr(PARAMFILE)
    open(PARAMFILE, file=trim(par_in_dir)//'/parameters_'//&
         &trim(adjustl(strnumber)), form='formatted')
    read(PARAMFILE, nml=scan,iostat=ierr)
    close(PARAMFILE)

    if ((ierr.ne.0).or.(chpt_in.eq.'no')) then
       !assume that the parameter increment in each dimension has a similar impact on the 
       !eigenvector (a more sophisticted metric of the parameter space would require a distinction
       !of the possible scan directions, like, e.g., beta and omt)
       do oldruns=1,pnum-1
          dist(oldruns)=sqrt(sum(abs(real(map_to_cartesian(pnum)-map_to_cartesian(oldruns)))**2))
       enddo

       if(pnum.gt.1) then
          !look for the minimum of the distance
          valid=.false.
          d_mask=.true.
          !stop if a valid checkpoint has been found or none are left to check
          do while ((.not.valid).and.any(d_mask))
             nearest=minloc(dist,1,d_mask)
             write(strnumber,form) nearest
             checkpoint_extension  ='_'//trim(adjustl(strnumber))
             !check if the corresponding run has finished already
             inquire(file=trim(diagdir)//'/checkpoint'//trim(checkpoint_extension), exist=valid, size=chpt_size)
             valid=(valid).and.(chpt_size.gt.6)
             !make sure that all procs of one group found the same result
             call mpi_bcast(valid,1,MPI_LOGICAL,0,comm_group,ierr)
             !otherwise continue the search without this entry
             if(.not.valid) d_mask(nearest)=.false.
          end do
       else
          valid=.false.
       endif

       if (valid) then
          checkpoint=trim(diagdir)//'/checkpoint'//checkpoint_extension
       else
          checkpoint='no'
       end if
    else
       checkpoint=chpt_in
    endif

  end subroutine get_closest_checkpoint
  
  function map_to_cartesian(runnr)
    integer,dimension(n_scan_dims):: map_to_cartesian
    integer:: runnr
    integer:: dim, tmp
    
    tmp=runnr-1
    
    do dim=1,n_scan_dims
       map_to_cartesian(dim)=tmp/product(scan_dims(dim+1:n_scan_dims))
       tmp=tmp-map_to_cartesian(dim)*product(scan_dims(dim+1:n_scan_dims))
    enddo
    
  end function map_to_cartesian
  
  function par_exists(pnum,form)
    logical:: par_exists
    integer:: pnum
    character(len=7):: form

    character(len=20):: f_ext
    character(len=7):: strnumber
    logical:: p_exists

    write(strnumber,form) pnum

    f_ext='_'//trim(adjustl(strnumber))
    if(mype.eq.0) inquire(file=trim(par_in_dir)//'/parameters'//&
         &trim(f_ext), exist=p_exists)
    call mpi_bcast(p_exists,1,MPI_LOGICAL,0,comm_group,ierr)
    par_exists=p_exists
  end function par_exists
  
  subroutine init_next_p(nextp)
    integer,intent(out):: nextp
    integer:: src, prob,num_s,ind
    character, dimension(:), allocatable:: prob_stat
    integer(MPI_OFFSET_KIND):: n_old
    integer:: n_old_i
    
    !initialize log file
    inquire(file=trim(par_in_dir)//'/gene_status', exist=continuation)
   

    if(mype.eq.0) then
       call mpi_file_open(comm_parall,trim(par_in_dir)//'/gene_status',&
            &MPI_MODE_RDWR+MPI_MODE_CREATE,MPI_INFO_NULL,LOGFILE,ierr)
       allocate(p_stat(0:n_parallel_sims-1))
       allocate(sr_tags(0:n_parallel_sims-1))
       sr_tags=1
       allocate(r_requests(0:n_parallel_sims-1),s_buf(0:49))
       sbi=0
       do src=0,n_parallel_sims-1
          if (src.ne.my_sim) call mpi_irecv(p_stat(src), 1, MPI_INTEGER, src, sr_tags(src), &
               comm_parall, r_requests(src), ierr)
       end do
       
       !initialize p_stat
       if (continuation) then
          call mpi_file_get_size(LOGFILE,n_old,ierr)
          allocate(prob_stat(n_old))
          n_old_i=n_old
          call mpi_file_read(LOGFILE,prob_stat, n_old_i, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
          num_s=0
          do ind=1,n_old
             if (prob_stat(ind)=='s') then
                if (num_s.lt.n_parallel_sims) p_stat(num_s)=ind
                num_s=num_s+1
             end if
          end do
          deallocate(prob_stat)
          do ind=num_s,n_parallel_sims-1
             p_stat(ind)=n_old+ind-num_s+1
          end do

          if (my_sim.eq.0) write(*,"(A,I5)") 'restarting old scan at parameters file ',minval(p_stat)           
       else
          do prob=1,n_parallel_sims
             p_stat(prob-1)=prob
          end do
       end if
       nextp=p_stat(my_sim)
    end if
    call mpi_bcast(nextp,1,MPI_INTEGER,0,comm_group,ierr)
    
  end subroutine init_next_p

  subroutine get_next_p(nextp)
    integer,intent(inout):: nextp
    integer:: src, dest, request
    logical:: flag, valid
    integer:: status(MPI_STATUS_SIZE)
    integer:: short_wait
    integer:: ierr

    !mark previous run as finished
    call update_logfile(nextp,'f')
    if(mype.eq.0) then
       valid=.false.
       nextp=0
       do while(.not.valid)
          !update p_stat
          do src=0,n_parallel_sims-1
             if (src.ne.my_sim) then
                call mpi_test(r_requests(src), flag, status, ierr)
                do while(flag)
                   sr_tags(src)=sr_tags(src)+1
                   call mpi_irecv(p_stat(src), 1, MPI_INTEGER, src, sr_tags(src), comm_parall, r_requests(src), ierr)
                   call mpi_test(r_requests(src), flag, status, ierr)
                end do
             end if
          end do
          
          !check for validity and exit if valid
          if(nextp.gt.0) valid=.true.
          do src=0,my_sim-1
             if (p_stat(src).eq.nextp) valid=.false.
          end do

          if(.not.valid) then
             !choose next problem
             nextp=maxval(p_stat)+1
             
             !broadcast choice
             s_buf(sbi)=nextp       
             p_stat(my_sim)=nextp
             do dest=0,n_parallel_sims-1
                if (dest.ne.my_sim) call mpi_isend(s_buf(sbi), 1, MPI_INTEGER, &
                     &dest, sr_tags(my_sim), comm_parall, request, ierr) 
             end do
             sr_tags(my_sim)=sr_tags(my_sim)+1       
             sbi=mod(sbi+1,size(s_buf))

             !wait one millisecond (interface to c)
             ierr=short_wait()
          end if
       end do
       
    end if

    call mpi_bcast(nextp,1,MPI_INTEGER,0,comm_group,ierr)

  end subroutine get_next_p

  subroutine update_logfile(nextp,stat_char)
    integer, intent(in):: nextp
    character, intent(in):: stat_char
    integer(MPI_OFFSET_KIND):: offset
    integer:: ierr

    if(mype.eq.0) then
       offset=nextp-1
       call mpi_file_write_at(LOGFILE, offset, stat_char, 1,&
            &MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    end if
    
  end subroutine update_logfile

  subroutine finalize_logfile

    call mpi_file_close(LOGFILE,ierr)
  end subroutine finalize_logfile
end module gene_scan
