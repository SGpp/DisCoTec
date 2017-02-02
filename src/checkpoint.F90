#include "intrinsic_sizes.h"
#include "redef.h"
!>Routines for reading and writing the distribution function from/into checkpoints
MODULE checkpoint
  Use par_mod
  use file_io, only: get_unit_nr
  Use communications
  use lagrange_interpolation
  use GaussQuadrature 
  use boundaries
  Use mpi
  use par_other, only: p_has_0_mode
#ifdef WITHFUTILS
  USE futils
#endif
#ifdef WITHHAC
  USE hlst_adios_checkpoint, ONLY: &
       & HAC_KIND, HAC_VERBOSE, hac_info, hac_read, hac_write, hac_init, hac_exit
#endif
  Implicit None
  Public:: initialize_checkpoint_module, finalize_checkpoint_module, &
       & initialize_checkpoint_read, checkpoint_read, finalize_checkpoint_read, &
       & initialize_checkpoint_write, checkpoint_write, finalize_checkpoint_write, &
       & chptdir, read_checkpoint, write_checkpoint, many_chpts, checkpoint_in, &
       & chpt_h5, chpt_read_h5, chpt_write_h5, &
       & chpt_hac, chpt_read_hac, chpt_write_hac, &
#ifdef COMBI_MGR
       & check_par_checkpoint,reset_chpt_time, istep_schpt, istep_g1, istep_gav, time_cp
#else
       & check_par_checkpoint,reset_chpt_time, istep_schpt, istep_g1, istep_gav
#endif
  Private

  !set chpt format (0=multiple binaries, 1=mpiio, 2=h5, 3=hac)
  Integer :: chpt_fmt_out = 0, chpt_fmt_in = 0
  Integer :: init_status = 0
#ifndef COMBI
  Character(Len=128):: chptdir  
#else
  Character(Len=FILENAME_MAX):: chptdir  
#endif
  Character(len=FILENAME_MAX) :: checkpoint_in = ''
  Logical:: read_checkpoint=.false., write_checkpoint=.true.
  Logical:: reset_chpt_time=.false., many_chpts=.false.
  Integer, dimension(6) :: arr_of_sizes, arr_of_subsizes, arr_of_starts
  integer:: garr,garr_r
  integer,target:: ch_handle=MPI_FILE_NULL, s_ch_handle=MPI_FILE_NULL, &
       &g1_handle=MPI_FILE_NULL, gav_handle=MPI_FILE_NULL
  integer(MPI_OFFSET_KIND),target:: ch_offset, s_ch_offset, g1_offset, &
       &gav_offset
  !MPI_Info object for tuning IO 
  integer :: info 
  integer:: istep_schpt=0, istep_g1=0, istep_gav=0
  character(len=6):: chpt_prec=''
  Logical:: chpt_h5=.false., chpt_read_h5=.false., chpt_write_h5=.false.
  Logical:: chpt_hac=.false., chpt_read_hac=.false., chpt_write_hac=.false.
#ifdef COMBI_MGR
  double precision:: time_cp
#endif

CONTAINS

  subroutine check_par_checkpoint

    !default checkpoint format is MPI-IO binary
    chpt_fmt_out = 1
    chpt_fmt_in = 1

    if (many_chpts.and.chpt_h5) then
       if (mype==0) then
          print*, 'WARNING: Only one checkpoint format can be selected'
          print*, '         Using chpt_h5 = T ...'
       endif
       many_chpts =.false.
    endif

    if (chpt_hac.and.chpt_h5) then
       if (mype==0) then
          print*, 'WARNING: Only one checkpoint format can be selected'
          print*, '         Using chpt_h5 = T ...'
       endif
       chpt_hac =.false.
    endif    

    if (chpt_h5) then
       chpt_read_h5 = .true.
       chpt_write_h5 = .true.
    endif

    if (chpt_hac) then
       chpt_read_hac = .true.
       chpt_write_hac = .true.
    endif

    if (many_chpts) then
       chpt_fmt_in = 0
       chpt_fmt_out = 0
    endif

    if (chpt_read_hac) chpt_fmt_in = 3
    if (chpt_write_hac) chpt_fmt_out = 3

    if (chpt_read_h5) chpt_fmt_in = 2
    if (chpt_write_h5) chpt_fmt_out = 2

  end subroutine check_par_checkpoint

  !>Initializes checkpoint module
  !!required by some libraries for general start-up
  subroutine initialize_checkpoint_module
    IMPLICIT NONE
#ifdef WITHHAC
    INTEGER, PARAMETER :: HACGVLVL = HAC_VERBOSE
    INTEGER :: ierr
#endif
   
    if (init_status.gt.0) return
    
#ifdef WITHHAC
    if (chpt_read_hac.or.chpt_write_hac) &
         CALL hac_init (MY_MPI_COMM_WORLD, HACGVLVL,ierr)
#endif
    
  end subroutine initialize_checkpoint_module


  !*************************************************************************!
  !******************* Checkpoint Input  ***********************************!
  !*************************************************************************!
  Subroutine initialize_checkpoint_read
    implicit none
    character(Len=FILENAME_MAX):: filename
    character(len=4)::  checkpnum
    integer,pointer:: handle
    integer(MPI_OFFSET_KIND),pointer:: offset
    integer:: ierr

    arr_of_sizes=(/ni0,nj0,nz0,nv0,nw0,n_spec/)
    arr_of_subsizes=(/li0,lj0,lk0,ll0,lm0,ln0/)
    arr_of_starts=(/my_pex*li0,my_pey*lj0,my_pez*lk0,my_pev*ll0,&
         &my_pew*lm0,my_pespec*ln0/)

    handle => ch_handle
    offset => ch_offset

    filename = trim(checkpoint_in)
    
    select case( chpt_fmt_in )
    case(0) !binary chpt per process
       call get_unit_nr(handle)
       write(checkpnum,"(i4.4)") mype
       filename=trim(filename)//'.'//checkpnum

       Open(handle, file=trim(filename), Form='unformatted', Status='old', Action='read')
       
       !determine precision of checkpoint
       read (handle) chpt_prec
       if (mype.eq.0) write(*,"(3a)") 'reading ',chpt_prec,' precision checkpoints'    
    case(1) !read MPI-IO checkpoint 
       call mpi_info_create(info,ierr)
       call mpi_info_set(info,"romio_ds_write", "disable",ierr)
       call mpi_info_set(info,"romio_ds_read", "disable",ierr)
       
       !call mpi_file_open(MY_MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,handle,ierr)
       call mpi_file_open(MY_MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,info,handle,ierr)
#ifdef COMBI
       write(*,"(2a)") filename
#endif
       if (ierr/=0) then
          write(*,"(2a)") 'ERROR while reading checkpoint ',trim(filename)
          STOP
       endif
       
       !check file size
       call mpi_file_get_size(handle,offset,ierr)
       if (abs(offset).le.6) then !using abs() for now as OpenMPI has an issue with files sizes >2GB
          if (mype.eq.0) Write(*,'(2(A))') 'ERROR: Invalid checkpoint ',trim(filename)
          stop
       endif
       
       !determine precision of checkpoint
       call mpi_file_read(handle,chpt_prec,6,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
       if (mype.eq.0) write(*,"(4a)") 'reading ',chpt_prec,' precision checkpoint  ',&
            &trim(filename)     
       offset=6
    case(2) !hdf5 checkpoint
       !nothing to be done
    case(3) !adios checkpoint
       !nothing to be done
    case default
       stop 'invalid chpt_fmt_in'
    end select
    
  End Subroutine initialize_checkpoint_read

  Subroutine checkpoint_read(g_,eof)
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), intent(inout):: g_
    logical, optional, intent(out):: eof
    Integer:: count1, count2, count_rate,garr_r
    Integer,dimension(6):: resolution, lind, uind,npoints
    Real   :: oldtime
    real(4):: timearr_s(2)
    real(8):: timearr_d(2)
    integer,pointer:: handle
    integer(MPI_OFFSET_KIND),pointer:: offset
    integer,dimension(MPI_STATUS_SIZE):: status
    Complex(4),dimension(:,:,:,:,:,:),allocatable:: check_in_s
    Complex(8),dimension(:,:,:,:,:,:),allocatable:: check_in_d
    Complex,dimension(:,:,:,:,:,:),allocatable:: check_in
    Integer:: ierr, count, MPI_COMPLEX_CHPT_TYPE
#ifdef WITHFUTILS
    INTEGER :: fidchkptin_h5, ni0_old, nj0_old, nz0_old, nv0_old, nw0_old, n_spec_old
    integer,dimension(6):: resolution_old
    complex, dimension(:,:,:,:,:,:), allocatable:: g_old
#endif
#ifdef WITHHAC
    character :: ntc
!    HAC_FIELD(KIND=NT_K),DIMENSION(:,:,:,:,:,:) :: g_old
#endif

    Call system_clock (count1)
    Call my_barrier
    PERFON('chpt_in')

    handle => ch_handle
    offset => ch_offset

    select case( chpt_fmt_in )
    case(0) !binary chpt per process
       if (chpt_prec.eq.'SINGLE') then
          read (handle,iostat=ierr) timearr_s
          if (ierr.ne.0) then 
             if (present(eof)) eof=.true.
             PERFOFF
             return
          end if
          time=timearr_s(1)
          dt=timearr_s(2)
          read (handle) resolution
          if (sum(abs(resolution-arr_of_sizes)).ne.0) stop &
               "checkpoint doesn't match the box specifications in parameters file"
          allocate(check_in_s(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
          read (handle) check_in_s
          g_=check_in_s
          deallocate(check_in_s)    
       else
          read (handle,iostat=ierr) timearr_d
          if (ierr.ne.0) then
             if (present(eof)) eof=.true.
             PERFOFF
             return
          end if
          time=timearr_d(1)
          dt=timearr_d(2)
          read (handle) resolution
          if (sum(abs(resolution-arr_of_sizes)).ne.0) stop &
               "checkpoint doesn't match the box specifications in parameters file"
          allocate(check_in_d(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
          read (handle) check_in_d
          g_=check_in_d
          deallocate(check_in_d) 
       end if
       
       oldtime = time
       If (mype /= 0) time = 0
       Call my_real_max_to_all (time)
       If (time.NE.oldtime) STOP&
            "checkpoint times do not match!"
    case(1) !read MPI-IO checkpoint 
       !read time, resolution and g_
       !step 1) read time and resolution
       !call mpi_file_set_view(handle,offset,MPI_BYTE,MPI_BYTE,"native",MPI_INFO_NULL,ierr)
       call mpi_file_set_view(handle,offset,MPI_BYTE,MPI_BYTE,"native",info,ierr)
       if (chpt_prec.eq.'SINGLE') then
          MPI_COMPLEX_CHPT_TYPE = MPI_COMPLEX
          call mpi_file_read(handle,timearr_s(1),2,MPI_REAL,status,ierr)
          call mpi_get_count(status,MPI_REAL,count,ierr)
       else
          MPI_COMPLEX_CHPT_TYPE = MPI_DOUBLE_COMPLEX
          call mpi_file_read(handle,timearr_d(1),2,MPI_DOUBLE_PRECISION,status,ierr)
          call mpi_get_count(status,MPI_DOUBLE_PRECISION,count,ierr)
       endif
       
       if (count.eq.0) then
          if (present(eof)) eof=.true.
          PERFOFF
          return
       end if
       
       call mpi_file_read(handle,resolution(1),6,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
       
       if (chpt_prec.eq.'SINGLE') then
          time=timearr_s(1)
          dt=timearr_s(2)
          offset=offset+2*4+6*4  !2*real, 6*integer
       else
          time=timearr_d(1)
          dt=timearr_d(2)
          offset=offset+2*8+6*4  !2*real, 6*integer
       endif
       
       call set_chpt_indices(resolution,lind,uind,npoints,arr_of_starts)
       
       !step 2) set mpi view and subarray type garr_r for reading checkpoint
       call mpi_type_create_subarray(6,resolution,npoints,&
            arr_of_starts,MPI_ORDER_FORTRAN,MPI_COMPLEX_CHPT_TYPE,garr_r,ierr)
       call mpi_type_commit(garr_r,ierr)
       
       !call mpi_file_set_view(handle,offset,MPI_BYTE,garr_r,"native",MPI_INFO_NULL,ierr)
       call mpi_file_set_view(handle,offset,MPI_BYTE,garr_r,"native",info,ierr)
       
       !step 3)  adapt precision and resolution: 2 cases 
       !step 3a) precision and resolution match: read to g_ directly
       
       if ((sum(abs(resolution-arr_of_sizes)).eq.0).and.(prec.eq.chpt_prec)) then
          call mpi_file_read_all(handle,g_(li1,lj1,lk1,ll1,lm1,ln1),li0*lj0*lk0*ll0*lm0*ln0,&
               MPI_COMPLEX_CHPT_TYPE,MPI_STATUS_IGNORE,ierr)
       else 
          !3b) precision and/or resolution change: temporary arrays needed
          !second temporary array could in principle be avoided when only precision or resolution change.
          allocate(check_in(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),&
               &lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))
          if (chpt_prec.eq.'SINGLE') then
             allocate(check_in_s(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),&
                  &lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))
             call mpi_file_read_all(handle,check_in_s(lind(1),lind(2),lind(3),lind(4),lind(5),lind(6)),&
                  npoints(1)*npoints(2)*npoints(3)*npoints(4)*npoints(5)*npoints(6),&
                  MPI_COMPLEX,MPI_STATUS_IGNORE,ierr)
             check_in = check_in_s
             deallocate(check_in_s)
          else 
             allocate(check_in_d(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),&
                  &lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))
             call mpi_file_read_all(handle,check_in_d(lind(1),lind(2),lind(3),lind(4),lind(5),lind(6)),&
                  npoints(1)*npoints(2)*npoints(3)*npoints(4)*npoints(5)*npoints(6),&
                  MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
             check_in = check_in_d
             deallocate(check_in_d)
          endif
          call copy_chpt_to_g1(resolution,check_in,lind,uind,g_)
          deallocate(check_in)
       endif
       
       !step 4) set file offset
       if (chpt_prec.eq.'SINGLE') then
          offset=offset+product(resolution)*8
       else
          offset=offset+product(resolution)*16
       endif
    case(2) !hdf5 checkpoint
#ifdef WITHFUTILS
       !fullpath_chkptr=trim(checkpoint_in)//'.h5' !!trim(chptdir)//'/checkpoint.h5'
       call openf(trim(checkpoint_in)//'.h5', fidchkptin_h5, 'r', mpi_comm_xyzvwspec)
       if (yx_order) then
          call getatt(fidchkptin_h5, '/dist/g_', 'nx0' , nj0_old)
          call getatt(fidchkptin_h5, '/dist/g_', 'nky0', ni0_old)
       else
          call getatt(fidchkptin_h5, '/dist/g_', 'nx0' , ni0_old)
          call getatt(fidchkptin_h5, '/dist/g_', 'nky0', nj0_old)
       endif
       call getatt(fidchkptin_h5, '/dist/g_', 'nz0' , nz0_old)
       call getatt(fidchkptin_h5, '/dist/g_', 'nv0' , nv0_old)
       call getatt(fidchkptin_h5, '/dist/g_', 'nw0' , nw0_old)
       call getatt(fidchkptin_h5, '/dist/g_', 'n_spec' , n_spec_old)
       call getatt(fidchkptin_h5, '/dist/g_', 'time', time)
       call getatt(fidchkptin_h5, '/dist/g_', 'dt', dt)
       
       resolution_old = (/ni0_old, nj0_old, nz0_old, nv0_old, nw0_old, n_spec_old/)
       
       if ((sum(abs(resolution_old-arr_of_sizes)).ne.0).or.(prec.NE.'DOUBLE')) then
          dt = dt_max
          call set_chpt_indices(resolution_old,lind,uind,npoints,arr_of_starts)
          allocate(g_old(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),&
               &lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))
          call getarrnd(fidchkptin_h5, '/dist/g_', g_old, (/6, 5, 4, 3, 2, 1/))
          call copy_chpt_to_g1(resolution_old,g_old,lind,uind,g_)
          deallocate(g_old)
       else
          call getarrnd(fidchkptin_h5, '/dist/g_', g_, (/6, 5, 4, 3, 2, 1/))
       end if
       call closef(fidchkptin_h5)
#endif
    case(3) !adios checkpoint
#ifdef WITHHAC
       CALL hac_info (TRIM(checkpoint_in)//'.hac',GADA_=resolution,ntc_=ntc,&
            ts_=time, st_=dt)

       IF (((prec.eq.'DOUBLE').and.(ntc.ne.'Z')).or.((prec.eq.'SINGLE').and.(ntc.ne.'C'))) &
            stop 'checkpoint precision conversion not yey implemented for HAC interface'

       if (sum(abs(resolution-arr_of_sizes)).eq.0) then
          CALL hac_read (g_, LBOUND(g_), TRIM(checkpoint_in)//'.hac',&
               SHAPE(g_))
       else 
          !3) resolution change: temporary arrays needed
          call set_chpt_indices(resolution,lind,uind,npoints,arr_of_starts)

          allocate(check_in(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),&
               &lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))
          CALL hac_read (check_in, LBOUND(check_in), &
               TRIM(checkpoint_in)//'.hac',SHAPE(check_in))
          call copy_chpt_to_g1(resolution,check_in,lind,uind,g_)
          deallocate(check_in)
       endif
#endif
       !nothing to be done
    case default
       stop 'invalid chpt_fmt_in'
    end select

    IF (nonlinear) THEN
       ! use checkpoint dt if smaller than dt_max due to nonlinear time step adaption
       dt = min(dt,dt_max)
    ELSE
       ! for linear runs, ignore checkpoint dt
       dt = dt_max
    END IF
    
    if (reset_chpt_time) time = 0.0

    Call my_barrier
    
    Call system_clock (count2, count_rate)
    If (mype == 0)&
         print "('Time for checkpoint_in:',F10.3,' sec')",&
         Real(count2-count1)/count_rate

    PERFOFF
    
  End Subroutine checkpoint_read
  
  Subroutine finalize_checkpoint_read
    integer,pointer:: handle
    integer:: ierr

    handle => ch_handle

    select case( chpt_fmt_in )
    case(0) !binary chpt per process
       close(handle)
    case(1) !read MPI-IO checkpoint 
       call mpi_file_close(handle,ierr)
       !call tracebackqq()
       call mpi_info_free(info,ierr)
    case(2) !hdf5 checkpoint
       !nothing to be done
    case(3) !adios checkpoint
       !nothing to be done
    case default
       stop 'invalid chpt_fmt_in'
    end select

    handle = MPI_FILE_NULL

  End Subroutine finalize_checkpoint_read


  !>Sets the indices used to read in checkpoints (also with different resolution)
  Subroutine set_chpt_indices(res,lind,uind,npoints,arr_of_starts)
    Integer, Dimension(6), INTENT(IN):: res
    Integer, Dimension(6), INTENT(INOUT):: lind, uind, npoints,arr_of_starts
    !resolution    <-old checkpoint
    !npoints       <-dimension of local array to be read by mpi i/o
    !arr_of_starts <-where to start reading in the file
    !lind,uind     <-for allocation of the local array to be read into
     
    !standard setup for unchanged resolution
    lind = (/li1,lj1,lk1,ll1,lm1,ln1/)
    uind = (/li2,lj2,lk2,ll2,lm2,ln2/)
    npoints = (/li0,lj0,lk0,ll0,lm0,ln0/)

    arr_of_starts = (/my_pex*npoints(1),my_pey*npoints(2),my_pez*npoints(3),my_pev*npoints(4),&
         my_pew*npoints(5),my_pespec*npoints(6)/)
    
    !ni0 change: read global array on each proc, then interpolate(global) or pick modes (local)
    If (ni0.ne.res(1)) then
       if (mype.eq.0) then
          if (yx_order) then
             WRITE(*,"(a)") 'read checkpoint: changing y resolution...'
          else
             WRITE(*,"(a)") 'read checkpoint: changing x resolution...'
          endif
       endif
       npoints(1)=res(1)
       arr_of_starts(1)=0
       lind(1)=0
       uind(1)=lind(1)+res(1)-1
    End if

    !kj modes are read on new grid until old domain ends, then
    if ((nj0.ne.res(2)).and.(mype.eq.0)) then
       if (yx_order) then
          WRITE(*,"(a)") 'read checkpoint: changing x resolution...'
       else
          WRITE(*,"(a)") 'read checkpoint: changing y resolution...'
       endif
    endif
    if(lj2 .gt. (res(2)-1)) then
       if (lj1 .gt. (res(2)-1))then
          !outside old domain: read one (irrelevant) point, since size of subarray cannot be 0
          uind(2)=lind(2)
          npoints(2)=1
          arr_of_starts(2)=0
       else
          !partially inside old domain
          uind(2)=res(2)-1
          npoints(2)=uind(2)-lind(2)+1
          !arr_of_starts(2)= my_pey*lj0 !start is inside old domain.
       endif
    endif

    !nz0 change: read global array on each proc, then interpolate and pick subarray
    if (nz0.ne.res(3)) then
       !stop 'Changing z resolution currently not allowed for continuation runs.'
       if (mype.eq.0) WRITE(*,"(a)") 'read checkpoint: changing z resolution...'
       npoints(3)=res(3)
       arr_of_starts(3)=0
       lind(3)=0
       uind(3)=lind(3)+res(3)-1
    endif

    !nv0 change: read global array on each proc, then interpolate and pick subarray
    If (nv0.ne.res(4)) then
       if (mype.eq.0) WRITE(*,"(a)") 'read checkpoint: changing vp resolution...'
       npoints(4)=res(4)
       arr_of_starts(4)=0
       lind(4)=0
       uind(4)=lind(4)+res(4)-1
    endif
    
    !nw0 change: read global array on each proc, then interpolate and pick subarray
    If (nw0.ne.res(5)) then
       if (mype.eq.0) WRITE(*,"(a)") 'read checkpoint: changing mu resolution...'
       if (mype.eq.0) WRITE(*,"(a)") 'WARNING: extrapolation at mu boundaries might be used'
       npoints(5)=res(5)
       arr_of_starts(5)=0
       lind(5)=0
       uind(5)=lind(5)+res(5)-1
    endif

    If (n_spec.ne.res(6)) stop 'Changing number of species is not allowed for continuation runs.'
  End Subroutine set_chpt_indices


  !>Copies checkpoint data to the initial distribution function (with transformations in case of changes in resolution)
  !!\todo Find a smart way to use this routine also many_chpts=.t. (without just copying it)
  !res is a vector of the resolution of the input checkpoint
  !chpt is a complex array as read in by mpi i/o, call with check_in_s or check_in_d not possible (data type)
  !lind and uind is the upper and lower bounds of chpt
  !g_is a g-type array on the new grid that does not have to be the input grid
  Subroutine copy_chpt_to_g1(res,chpt,lind,uind,g_)
    Complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), intent(inout):: g_
    Integer, Dimension(6), INTENT(IN):: res, lind, uind
    Complex, Dimension(lind(1):uind(1),lind(2):uind(2),lind(3):uind(3),lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)), &
         INTENT(IN):: chpt
    !two arrays to swap around the checkpoint data as we cycle through the dimensions and adapt the resolution
    Complex, Dimension(:,:,:,:,:,:), Allocatable:: tmp1,tmp2
    Complex, Dimension(:,:,:), Allocatable:: tmp3
    Integer:: i, j, k, l, m, n, nind
    real,dimension(:),allocatable::axis_in
    complex,dimension(:),allocatable::data_out, data_in
    real,dimension(1:res(5))::old_mu_weight
    real::dz_old,dv_old,dx_old,x0_rhostar
    integer:: lk1_in, lk2_in, lbz_in, ubz_in, n_procs_z_in

    PERFON('chpt_trafo')
    if (sum(abs(res-(/ni0,nj0,nz0,nv0,nw0,n_spec/))).eq.0) then
       g_=chpt
    else
       !different ni0 --- add/remove modes in local mode, interpolate in global mode
       Allocate(tmp1(li1:li2,lind(2):uind(2),lind(3):uind(3),lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))       
       if (xy_local.and.ni0.ne.res(1)) then
          If (nx0.gt.res(1)) then
             nind=res(1)
          Else if (ni0.lt.res(1)) then
             nind=ni0
          End if
          tmp1 = 0.0 !initialize with zeros
          !paste kx.ge.0
          tmp1(li1:li1+(nind-1)/2+evenx,:,:,:,:,:) = chpt(lind(1):lind(1)+(nind-1)/2+evenx,:,:,:,:,:)
          !paste kx.lt.0
          tmp1(li2-(nind-1)/2+1:li2,:,:,:,:,:) = chpt(uind(1)-(nind-1)/2+1:uind(1),:,:,:,:,:)
       else if (.not.x_local.and.ni0.ne.res(1)) then
          !arrays for interpolation (without boundaries)
          allocate(axis_in(lind(1):uind(1)))
          allocate(data_in(lind(1):uind(1)))
          allocate(data_out(0:ni0-1))
          
          !initialise old radial coordinate
          if (rad_bc_type.eq.0) then
             if (mype.eq.0) print*,"WARNING: extrapolation in x for increasing resolution"
             !global periodic case: x=-lx/2,..,lx/2-dx
             dx_old = lx / res(1)
          else
             !global and nonperiodic: x=-lx/2,..,lx/2
             dx_old = lx / (res(1)-1) 
          endif
          do i=0,res(1)-1
             axis_in(i) = -0.5*lx + i*dx_old
          enddo
          x0_rhostar = (xval(0)+0.5*lx)
          axis_in = axis_in + x0_rhostar
          
          do n=lind(6),uind(6)
             do m=lind(5),uind(5)
                do l=lind(4),uind(4)
                   do k=lind(3),uind(3)
                      do j=lind(2),uind(2)
                         data_in = chpt(lind(1):uind(1),j,k,l,m,n)
                         call lag3interp_complex(data_in,axis_in,res(1),&
                              data_out,xval,nx0)
                         tmp1(li1:li2,j,k,l,m,n)=data_out(li1:li2)
                      enddo
                   enddo
                enddo
             enddo
          enddo
          deallocate(axis_in,data_in,data_out)
       else if (.not.y_local.and.ni0.ne.res(1)) then
          stop 'Changing resolution for y_local=.false. is not implemented'
       else
          tmp1=chpt
       endif

       !different nky0 --- set modes out of old domain to zero
       Allocate(tmp2(li1:li2,lj1:lj2,lind(3):uind(3),lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))           
       If (nj0.ne.res(2).and.nonlinear) then
          if (.not.yx_order) then
             !local indices for copy
             nind=lj0
             if (lj2.gt.res(2)-1) nind=res(2)-lj1
             if (lj1.gt.res(2)-1) nind=0
             !initialize with zeros
             tmp2=0.0
             !initialize with fraction of 0-mode when starting from a
             !neoclassical equilibrium
             if (res(2) .eq. 1 .and. p_has_0_mode) then
                do i=li1,li2
                    do j=lj1,lj2
                        tmp2(i,j,:,:,:,:) = 0.01 * tmp1(i,lind(2),:,:,:,:)
                    end do
                end do 
             end if
             Do i=li1,li2
                tmp2(i,lj1:lj1+nind-1,:,:,:,:)=tmp1(i,lind(2):lind(2)+nind-1,:,:,:,:)
             End do
          else
             stop 'Changing x resolution for y_local=.false. is not implemented yet'
          endif
       Else if (.not.nonlinear.and.nky0.ne.res(2)) then
          stop 'Changing number of ky modes is not allowed for linear continuation runs.'
       Else
          tmp2=tmp1
       End if
       Deallocate(tmp1)

       !different nz0 --- global array is read on each proc, interpolate and distribute
       Allocate(tmp1(li1:li2,lj1:lj2,lk1:lk2,lind(4):uind(4),lind(5):uind(5),lind(6):uind(6)))       
       If (nz0.ne.res(3)) then
          allocate(axis_in(lind(3)-nzb:uind(3)+nzb))
          allocate(data_in(lind(3)-nzb:uind(3)+nzb))
          allocate(data_out(0:nz0-1))
          allocate(tmp3(li1:li2,lj1:lj2,(lind(3)-nzb):(uind(3)+nzb)))

          !initialise old parallel coordinate
          dz_old = (2.0 * pi * n_pol) / res(3)
          do k=lind(3)-nzb,uind(3)+nzb
             axis_in(k) = -pi*n_pol + k*dz_old !+dz_old/2
          enddo
          if(mod(res(3),2).ne.0) then
             axis_in = axis_in+dz_old/2
          end if

          !temporarily redefine z discretization
          lk1_in=lk1; lk2_in=lk2; lbz_in=lbz; ubz_in=ubz
          lk1=lind(3)
          lk2=uind(3)
          lbz=lind(3)-nzb
          ubz=uind(3)+nzb
          !arrays are not distributed:
          !pretend to switch off z parallelization
          n_procs_z_in=n_procs_z;
          n_procs_z=1

          do n=lind(6),uind(6)
             do m=lind(5),uind(5)
                do l=lind(4),uind(4)
                   !xyz array with old resolution plus ghost cells
                   tmp3(:,:,lind(3):uind(3))=tmp2(:,:,lind(3):uind(3),l,m,n)
                   call exchange_z_noinner(tmp3) !exchange w/o z parallelization
                   do j=lj1,lj2
                      do i=li1,li2
                         data_in = tmp3(i,j,:)
                         call lag3interp_complex(data_in,axis_in,res(3)+2*nzb,&
                            data_out,zval,nz0)
                         tmp1(i,j,lk1_in:lk2_in,l,m,n)=data_out(lk1_in:lk2_in)
                      enddo
                   enddo
                enddo
             enddo
          enddo 
          deallocate(axis_in,data_in,data_out,tmp3)

          !restore z discretization and parallelization
          lk1=lk1_in; lk2=lk2_in; lbz=lbz_in; ubz=ubz_in
          n_procs_z=n_procs_z_in;

       else
          tmp1=tmp2
       endif
       Deallocate(tmp2)

       !different nv0 --- global array is read on each proc, interpolate and distribute
       Allocate(tmp2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lind(5):uind(5),lind(6):uind(6)))       
       If (nv0.ne.res(4)) then
          !arrays for interpolation (without boundaries)
          allocate(axis_in(lind(4):uind(4)))
          allocate(data_in(lind(4):uind(4)))
          allocate(data_out(0:nv0-1))
          
          !initialise old coordinate
#ifndef COMBI
          dv_old = 2.0*lv/(res(4)-1.0)
         !do l=-2,res(4)+1
          do l=0,res(4)-1
             axis_in(l) = -lv + l*dv_old
          end do
#else
          dv_old = 2.0*lv/(res(4))
         !do l=-2,res(4)+1
          do l=0,res(4)-1
             axis_in(l) = -lv + l*dv_old + SHIFT
          end do
#endif


          do n=lind(6),uind(6)
             do m=lind(5),uind(5)
                do k=lk1,lk2
                   do j=lj1,lj2
                      do i=li1,li2
                         data_in = tmp1(i,j,k,lind(4):uind(4),m,n)
                         call lag3interp_complex(data_in,axis_in,res(4),&
                              data_out,vp,nv0)
                         tmp2(i,j,k,ll1:ll2,m,n)=data_out(ll1:ll2)
                      enddo
                   enddo
                enddo
             enddo
          enddo

          deallocate(axis_in,data_in,data_out)
       else
          tmp2=tmp1
       endif
       Deallocate(tmp1)
       
       
       !different nw0 --- global array is read on each proc, interpolate and distribute
       Allocate(tmp1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,lind(6):uind(6)))       
       If (nw0.ne.res(5)) then
          !arrays for interpolation (without boundaries)
          allocate(axis_in(lind(5):uind(5)))
          allocate(data_in(lind(5):uind(5)))
          allocate(data_out(0:nw0-1))
          
          !initialise old coordinate

          call GetMuWeightsAndKnots(old_mu_weight,axis_in,lw,res(5))

          do n=lind(6),uind(6)
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      do i=li1,li2
                         data_in = tmp2(i,j,k,l,lind(5):uind(5),n)
                         call lag3interp_complex(data_in,axis_in,res(5),&
                              data_out,mu,nw0)
                         tmp1(i,j,k,l,lm1:lm2,n)=data_out(lm1:lm2)
                      enddo
                   enddo
                enddo
             enddo
          enddo

          deallocate(axis_in,data_in,data_out)
       else
          tmp1=tmp2
       endif
       Deallocate(tmp2)

       !no change of species allowed! 

       g_=tmp1
       Deallocate(tmp1)
    endif !resolution changed

    PERFOFF
  End Subroutine copy_chpt_to_g1
  
  
  !*************************************************************************!
  !******************* Checkpoint Output ***********************************!
  !*************************************************************************!

  subroutine initialize_checkpoint_write
    character(Len=16):: filename
    character(len=4)::  checkpnum
    character(Len=FILENAME_MAX):: fullpath
    integer,pointer:: handle
    integer(MPI_OFFSET_KIND),pointer:: offset
    integer:: f, ierr

    !print*,mype,": initialize_checkpoint_write start"
    arr_of_sizes=(/ni0,nj0,nz0,nv0,nw0,n_spec/)
    arr_of_subsizes=(/li0,lj0,lk0,ll0,lm0,ln0/)
    arr_of_starts=(/my_pex*li0,my_pey*lj0,my_pez*lk0,my_pev*ll0,&
         &my_pew*lm0,my_pespec*ln0/)

    !open all write_checkpoint based files if more efficient
    do f=1,4
       select case (f)
       case(1)   
          if (write_checkpoint) then
             filename='checkpoint' 
             handle => ch_handle
             offset => ch_offset
          else
             cycle
          end if
       case(2)
          if (istep_schpt.gt.0) then
             filename='s_checkpoint' 
             handle => s_ch_handle
             offset => s_ch_offset
          else
             cycle
          end if
       case(3)
          if(istep_g1.gt.0) then
             filename='g1'//trim(file_extension)
             handle => g1_handle
             offset => g1_offset
          else
             cycle
          endif
       case(4)
          if(istep_gav.gt.0) then
             filename='gav'//trim(file_extension)
             handle => gav_handle
             offset => gav_offset
          else
             cycle
          endif
       end select

       If ((f.lt.3).and.trim(file_extension).ne.'.dat') then
          fullpath=trim(chptdir)//'/'//trim(filename)//trim(file_extension)
       else
          fullpath=trim(chptdir)//'/'//trim(filename)
       endif

       select case( chpt_fmt_out )
       case(0) !binary chpt per process
          write(checkpnum,"(i4.4)") mype
          call get_unit_nr(handle)   
          fullpath=trim(fullpath)//'.'//checkpnum
          Open(handle, file=trim(fullpath), status='replace',&
               form='unformatted',action='write')
          if(f.ne.3) write (handle) prec
       case(1) !MPI-IO checkpoint 
          !create mpi type
          call mpi_type_create_subarray(6,arr_of_sizes,arr_of_subsizes,&
               &arr_of_starts,MPI_ORDER_FORTRAN,MPI_COMPLEX_TYPE,garr,ierr)
          call mpi_type_commit(garr,ierr)
          call mpi_info_create(info,ierr)
          call mpi_info_set(info,"romio_ds_write", "disable",ierr)
          call mpi_info_set(info,"romio_ds_read", "disable",ierr)
          ! call mpi_file_delete(fullpath,MPI_INFO_NULL,ierr)
          ! depending on system this (or workaround with dummy file ) caused problems
          ! because slower procs deleted the already created NEW checkpoint
          if (mype == 0) then
             open(99+my_sim,file=fullpath,status='replace')
             close(99+my_sim,status='delete')
          endif
          !open file 
          !call mpi_file_open(MY_MPI_COMM_WORLD,fullpath,MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,handle,ierr)
          call mpi_file_open(MY_MPI_COMM_WORLD,fullpath,MPI_MODE_WRONLY+MPI_MODE_CREATE,info,handle,ierr)
          if(f.ne.3) then
             !write header
             if (mype.eq.0) call mpi_file_write(handle,prec,6,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
             offset=6  !offset in byte: 6*character
          else
             offset=0
          end if
       case(2) !hdf5 checkpoint
#ifdef WITHFUTILS
          fullpath=trim(fullpath)//'.h5'
          call creatf(trim(fullpath),handle, "Checkpoint", 'd', mpi_comm_xyzvwspec)
#endif
       case(3) !adios checkpoint
          !nothing to be done
       case default
          stop 'invalid chpt_fmt_out'
       end select
    enddo
    !print*,mype,": initialize_checkpoint_write end"

  end subroutine initialize_checkpoint_write

  subroutine checkpoint_write (g_,filename)
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2), intent(in):: g_
    character(Len=*), intent(in):: filename
    integer,pointer:: handle
    integer(MPI_OFFSET_KIND),pointer:: offset
    integer:: count1, count2, count_rate, ierr, rkind
    real, dimension(2) :: time_arr
    logical:: write_res, frewind, append
    CHARACTER(len=FILENAME_MAX) :: fullname
    CHARACTER(len=16) :: gdesc

    PERFON('chpt_out')
    !print*,mype,": Starting to write checkpoint"

    Call system_clock (count1)

    time_arr = (/time, dt/)

    rkind = 4
    if (prec.eq.'DOUBLE') rkind=rkind*2
    gdesc = '/dist'

    frewind=.false.
    append=.false.
    write_res=.true.
    if (filename.eq.'checkpoint') then 
       handle => ch_handle
       offset => ch_offset
    elseif (filename.eq.'s_checkpoint') then
       handle => s_ch_handle
       offset => s_ch_offset
       frewind=.true.
    elseif (filename.eq.'g1') then
       handle => g1_handle
       offset => g1_offset
       write_res=.false.
       append=.true.
    elseif (filename.eq.'gav') then
       handle => gav_handle
       offset => gav_offset
    endif

    fullname=trim(filename)
    if (trim(file_extension).ne.'.dat') &
         &fullname=trim(fullname)//trim(file_extension)

    select case( chpt_fmt_out )
    case(0) !binary chpt per process
       !file is already open
       if (frewind) then
          rewind(handle)
          if (write_res) write (handle) prec
       endif
       if (write_res) then
          write (handle) time_arr
          write (handle) arr_of_sizes
       else
          write (handle) time
       endif
       write (handle) g_
    case(1) !MPI-IO checkpoint 
       !write header
       !call mpi_file_set_view(handle,offset,MPI_BYTE,MPI_BYTE,"native",MPI_INFO_NULL,ierr)
       call mpi_file_set_view(handle,offset,MPI_BYTE,MPI_BYTE,"native",info,ierr)
       if (mype.eq.0) then
          if(write_res) then 
             call mpi_file_write(handle,time_arr,2,MPI_REAL_TYPE,MPI_STATUS_IGNORE,ierr)
             call mpi_file_write(handle,arr_of_sizes,6,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
          else
             call mpi_file_write(handle,time,1,MPI_REAL_TYPE,MPI_STATUS_IGNORE,ierr)
          endif
       endif
       if(write_res) then
          !offset in byte: 2*real, 6*integer
          offset=offset+2*rkind+6*4
       else
          !offset: 1*real
          offset=offset+rkind
       endif
       !call mpi_file_set_view(handle,offset,MPI_BYTE,garr,"native",MPI_INFO_NULL,ierr)
       call mpi_file_set_view(handle,offset,MPI_BYTE,garr,"native",info,ierr)
       call mpi_file_write_all(handle,g_(li1,lj1,lk1,ll1,lm1,ln1),li0*lj0*lk0*ll0*lm0*ln0,&
            MPI_COMPLEX_TYPE,MPI_STATUS_IGNORE,ierr)
       if (.not.frewind) then
          !update offset to last position of the file
          offset=offset+li0*lj0*lk0*ll0*lm0*ln0*n_procs_sim*2*rkind
       else
          !rewind offset to position after header in order to overwrite with the next call
          offset=offset-2*rkind-6*4
       end if
    case(2) !hdf5 checkpoint
#ifdef WITHFUTILS
       if (append) then
          write(gdesc,"(a,i8.8)") '/dist',itime
       elseif (frewind) then !reset file
          fullname=trim(fullname)//'.h5'
          call closef(handle)
          call creatf(trim(chptdir)//trim(fullname),handle, &
               &"Checkpoint", 'd', mpi_comm_xyzvwspec)
       endif
       call creatg(handle, trim(gdesc))
       call putarrnd(handle, trim(gdesc)//'/g_', &
            & g_, (/6, 5, 4, 3, 2, 1/), desc='distribution for all species')
       call attach(handle, trim(gdesc)//'/g_', 'time', time)
       call attach(handle, trim(gdesc)//'/g_', 'dt', dt)
       call attach(handle, trim(gdesc)//'/g_', 'nx0' , nx0)
       call attach(handle, trim(gdesc)//'/g_', 'nky0', nky0)
       call attach(handle, trim(gdesc)//'/g_', 'nz0' , nz0)
       call attach(handle, trim(gdesc)//'/g_', 'nv0' , nv0)
       call attach(handle, trim(gdesc)//'/g_', 'nw0' , nw0)
       call attach(handle, trim(gdesc)//'/g_', 'n_spec', n_spec)
       call attach(handle, trim(gdesc)//'/g_', 'precision', prec)
       call flushh5(handle)
#endif
    case(3) !adios checkpoint
#ifdef WITHHAC
       fullname=trim(fullname)//'.hac'
       CALL hac_write (g_, LBOUND(g_), TRIM(chptdir)//TRIM(fullname),&
            time,dt)
#endif
    case default
       stop 'invalid chpt_fmt_out'
    end select

    call system_clock (count2, count_rate)
    if ((mype == 0).and.print_ini_msg) THEN
       print "('Time for ',a,' out:',F10.3,' sec')",&
            & trim(fullname), real(count2-count1)/count_rate          
    end if

#ifdef COMBI_MGR
    time_cp = Real(count2-count1)/count_rate
#endif

    PERFOFF

  end subroutine checkpoint_write

  subroutine finalize_checkpoint_write
    integer,pointer:: handle
    integer:: f, ierr
    integer(MPI_OFFSET_KIND) :: offset
    character(Len=FILENAME_MAX):: fullpath

    if (chpt_fmt_out.eq.1) then
       ! free the info object
       call mpi_info_free(info,ierr)
       ! free garr type
       call mpi_type_free(garr,ierr)
    endif

    !close files
    do f=1,4
       handle=>null()
       select case (f)
       case(1)   
          handle => ch_handle
       case(2)
          handle => s_ch_handle
       case(3)
          handle => g1_handle
       case(4)
          handle => gav_handle
       end select
       if (.not.associated(handle).or.(handle.eq.MPI_FILE_NULL)) cycle
       
       
       select case( chpt_fmt_out )
       case(0) !binary chpt per process
          close(handle)
       case(1) !MPI-IO checkpoint 

          !get file size
          call mpi_file_get_size(handle,offset,ierr)
          call mpi_file_close(handle,ierr)
          !erase checkpoint files if no distribution function is stored
          !!using abs(offset) for now as OpenMPI has an issue with files sizes >2GB
          !!and returns negative numbers: https://svn.open-mpi.org/trac/ompi/ticket/2145
          if ((f==1).and.(mype==0).and.(abs(offset).le.6)) then
             fullpath=trim(chptdir)//'/checkpoint'
             if (trim(file_extension).ne.'.dat') &
                  &fullpath=trim(fullpath)//trim(file_extension)
             open(99+my_sim,file=fullpath,status='replace')
             close(99+my_sim,status='delete')
          endif
       case(2) !hdf5 checkpoint
#ifdef WITHFUTILS
          call closef(handle)
#endif
       case(3) !adios checkpoint
          !nothing to be done
       case default
          stop 'invalid chpt_fmt_in'
       end select
    end do
    
  end subroutine finalize_checkpoint_write

  !>Finalize checkpoint module
  subroutine finalize_checkpoint_module
    implicit none
#ifdef WITHHAC
!    if (chpt_read_hac.or.chpt_write_hac) &
!         CALL hac_exit
    init_status = 0
#endif
  end subroutine finalize_checkpoint_module


END MODULE checkpoint
