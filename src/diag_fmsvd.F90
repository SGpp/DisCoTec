#include "redef.h"
!> This modules contains routines for SVD diagnostics 
Module diagnostics_fmsvd
  Use par_mod 
  Use file_io, only: get_unit_nr
  use geometry
  use communications
  use mpi
  use arrays
  use phys_ini

  Implicit None

  PUBLIC :: svd_projection_fmom  !########

  PRIVATE 

!For parallel svd:

#ifdef WITHSCAL
#ifdef with_extended_diags
  integer :: nrgcols
  !!!Matrix dimensions
  integer :: num_col,num_row
!  integer :: num_gridp

  !!!MPI and BLACS
  integer :: mpi_comm_mat
  integer :: mypeb,my_context,procs_blacs,myrow,mycol
  integer :: mype_mat(2)
  integer, dimension(:,:), allocatable :: blacs_map
  !integer :: MPI_REAL_TYPE
  integer :: n_row_procs, n_col_procs

  !!!!Arrays and Matrices
  complex, dimension(:,:), allocatable  :: Amat_loc
  complex, dimension(:,:), allocatable  :: Umat_loc
  complex, dimension(:,:), allocatable  :: VTmat_loc
  real, dimension(:), allocatable  :: svals
  complex, allocatable, dimension(:,:) :: mat_glob

  !n=global
  !l=local
  !r=row
  !c=column
  integer :: nr0A  !Total number of rows (global)
  integer :: nc0A  !Total number of columns (global)
  integer :: pr0A  !Processor block size
  integer :: pc0A  !Processor block size
  !Number of row/column elements for each processor
  integer :: lr0A
  integer :: lc0A

  integer :: nr0U  !Total number of rows (global)
  integer :: nc0U  !Total number of columns (global)
  integer :: pr0U  !Processor block size
  integer :: pc0U  !Processor block size
  !Number of row/column elements for each processor
  integer :: lr0U
  integer :: lc0U

  integer :: nr0VT  !Total number of rows (global)
  integer :: nc0VT  !Total number of columns (global)
  integer :: pr0VT  !Processor block size
  integer :: pc0VT  !Processor block size
  !Number of row/column elements for each processor
  integer :: lr0VT
  integer :: lc0VT

  !Min(num_row,num_col)
  integer :: size0

  !Descriptors
  integer :: desc_U(9),desc_VT(9),desc_A(9)
  integer :: block_size


  real, allocatable, dimension(:) :: time_arr
  !for integral weight
  real,allocatable, dimension(:) :: jac_glob
#endif 
#endif 

Contains

!!!******************************************************************!!!

  subroutine svd_projection_fmom 

    if(svd_field_mom) then
#ifdef with_extended_diags
#ifdef WITHSCAL
      call initialize_svd_fmom
      call svd_fmom
      call finalize_svd_fmom
#else
      if(mype==0) write(*,*) "Need scalapack to use parallel svd routine!"
      stop
#endif
#endif
    end if

  end subroutine svd_projection_fmom 

#ifdef with_extended_diags
#ifdef WITHSCAL

  subroutine initialize_svd_fmom

    real :: temp1d(0:nz0-1)
    integer :: ierr

    final_init=.true.
    call split_comm
    call initialize_discretization(.true.)
    call allocate_arrays
    call init_physical

    allocate(jac_glob(0:nz0-1))

    temp1d=0.0
    temp1d(lk1:lk2)=geom%jacobian(pi1,pj1,lk1:lk2)
    call MPI_ALLREDUCE(temp1d,jac_glob,nz0,MPI_REAL_TYPE,MPI_SUM,mpi_comm_z,ierr)

    nrgcols=10
    !if(.not.momentum_flux) nrgcols=8
    !if(mype==0)  write(*,*) "jac_glob",jac_glob

  end subroutine initialize_svd_fmom 

  subroutine finalize_svd_fmom
    call finalize_physical
    call svd_deallocate_arrays
    call free_comm
    deallocate(jac_glob)
  end subroutine finalize_svd_fmom

!This subroutine calculates the SVD of phi, Tpar and Tperp
!At this time this routine only analyzes ion moments
!Instructions: Use the parameters.dat file from a previous simulation and modify:
!SVD_field_mom=T
!Input parameters (in extended_diags namelist):
!SVD_df_n_time: specify number of time steps to analyze for each run: note that
!this is for the field file--number of mom steps is calculated from this
!SVD_start_time (default:0) : start time for analysis
!SVD_sparse_factor (default:1): analysis takes every SVD_sparse_factor'th time step
!SVD_df_file_path: default is diagdir
!SVD_n_restarts: number of restarts for df files (default = 1, i.e. one run w/o restarts)
!SVD_df_file_suffix (default:'.dat'): array of file endings for a series of runs
  subroutine svd_fmom
    implicit none

    !integer :: i,i0,p,q,lr,lc,pr,pc
    integer :: numxyz
    integer :: s_index_field(SVD_n_restarts)  !start index
    integer :: s_index_mom(SVD_n_restarts)  !start index
    integer :: n_time,n_time_mom

    if(mype==0) write(*,*) "Performing SVD of field and mom data."
    block_size=64
    numxyz=nx0*nky0*nz0
    num_row=3*numxyz

    if(mype==0) write(*,*) "Getting time data for field file."
    SVD_df_file_name='field'
    call read_fmom_time(SVD_df_file_name,SVD_sparse_factor*SVD_istep_mf_factor,1,s_index_field,n_time)
    if(mype==0) write(*,*) "Getting time data for mom file."
    SVD_df_file_name='mom_'//trim(spec(0)%name)
    call read_fmom_time(SVD_df_file_name,SVD_sparse_factor,2,s_index_mom,n_time_mom)
    if(n_time_mom.ne.n_time) then
      write(*,*) "n_time_mom.ne.n_time!!!",mype
    end if
    num_col=n_time
    allocate(time_arr(num_col))
    time_arr=0.0

    !call init_comm
    if(mype==0) write(*,*) "Initializing communicators."
    call comm_svd
    if(mype==0) write(*,*) "Allocating distributed matrices."
    call arrays_A

    if(mype==0) write(*,*) "Reading potential"
    SVD_df_file_name='field'
    call read_fmom(SVD_df_file_name,SVD_sparse_factor*SVD_istep_mf_factor,1,s_index_field,n_time)
    if(mype==0) write(*,*) "Reading Tpar"
    SVD_df_file_name='mom_'//trim(spec(0)%name)
    call read_fmom(SVD_df_file_name,SVD_sparse_factor,2,s_index_mom,n_time_mom)
    if(mype==0) write(*,*) "Reading Tperp"
    call read_fmom(SVD_df_file_name,SVD_sparse_factor,3,s_index_mom,n_time_mom)


   if(mype==0) write(*,*) "Done reading data from file."
   call arrays_UV
   if(mype==0) write(*,*) "Calculating SVD."
   call get_svd
   if(mype==0) write(*,*) "Outputting data."
   call output_data_parallel_fm_beta
   if(mype==0) write(*,*) "Deallocating arrays."
   call svd_deallocate_arrays
   if(mype==0) write(*,*) "Finalizing blacs."
   !call finalize_blacs
   if(mype==0) write(*,*) "Done with parallel svd."
   deallocate(time_arr)

  end subroutine svd_fmom

  subroutine read_fmom(file_name,this_sparse_factor,read_number,s_index,n_time)
    implicit none
    character(len=128), intent(in) :: file_name
    character(len=128)  :: file_path_name(20)
    integer, intent(in) :: this_sparse_factor
    integer, intent(in) :: read_number
    integer :: field_handle
    complex, allocatable, dimension(:) :: fm_in
    complex, allocatable, dimension(:) :: fm_in2
    integer :: numxyz,p2
    integer :: num_fmom,ifmom
    integer :: i,i0,p,q,lr,lc,pr,pc
    integer :: rnum
    integer, intent(in) :: s_index(SVD_n_restarts),n_time  !start index
    integer :: INDXG2L,  INDXG2P,ierr
    real :: dummy_r
    integer :: df_n_time(SVD_n_restarts)
!    real :: temp

    block_size=64
    numxyz=nx0*nky0*nz0
    if(read_number==1) then
      num_fmom=n_fields
      df_n_time=SVD_df_n_time
    end if
    if(read_number==2) then
      num_fmom=n_moms
      df_n_time=SVD_df_n_time/SVD_istep_mf_factor
    end if
    if(read_number==3) then
      num_fmom=n_moms
      df_n_time=SVD_df_n_time/SVD_istep_mf_factor
    end if
    do i=1,20
      file_path_name(i)= trim(SVD_df_file_path)//'/'//trim(file_name)//&
             trim(SVD_df_file_suffix(i))
    end do

    allocate(fm_in(numxyz))
    allocate(fm_in2(numxyz))


   i0=0
   do rnum=1,SVD_n_restarts
    call get_unit_nr(field_handle)
    if(mype==0) &
#ifdef F2003_NO_OPEN_CONVERT
         &open(field_handle,file=file_path_name(rnum),&
         &form='unformatted',status='unknown')
#else
         &open(field_handle,file=file_path_name(rnum),&
         &form='unformatted',status='unknown',convert=SVD_endian)
#endif
    i=0
    do while(i.lt.df_n_time(rnum))
      i=i+1 
      if(mod(i,this_sparse_factor)==0.and.i.ge.s_index(rnum).and.i0.lt.n_time) then
        i0=i0+1
        if(mype==0) read(field_handle) time_arr(i0)
        !if(mype==0) write(*,*) "tk",time_arr(i0)

        do ifmom=1, num_fmom
          if(mype==0) read(field_handle) fm_in
          if(ifmom==read_number) fm_in2=fm_in   !1 for potential in field file
                                                !2 for Tpar in mom file
                                                !3 for Tperp in mom file
        end do
        fm_in=fm_in2
        !!!!!!!!!!!!Testing
        !if(mype==0) write(*,*) abs(sum(fm_in))
        !if(mype==0) call fm_sum_3d(fm_in,temp)
        !if(mype==0) write(*,*) "test",read_number,temp
        !!!!!!!!!!!!Testing
        if(mype==0) call integral_weight_fm(fm_in,fm_in2)
        call MPI_BCAST(fm_in2,numxyz,MPI_COMPLEX_TYPE,0,my_mpi_comm_world,ierr) 
        q=i0
        do p=1,numxyz  !Get local processor and indices from global indices
          p2=p+numxyz*(read_number-1)
          !if(mype==0) write(*,*) q,p,p2
          !Find the local matrix indices corresponding to the global indices
          lr=INDXG2L(p2,pr0A,0,0,n_row_procs)
          lc=INDXG2L(q,pc0A,0,0,n_col_procs)
          !Find the processors corresponding to these global indices
          pr=INDXG2P(p2,pr0A,0,0,n_row_procs)
          pc=INDXG2P(q,pc0A,0,0,n_col_procs)
          if(myrow==pr.and.mycol==pc) then
            Amat_loc(lr,lc)=fm_in2(p)
          end if
        end do
        if(mype==0) write(*,*) "Done reading and distributing step:", i0

      else
        if(mype==0) read(field_handle) dummy_r
        !if(mype==0) write(*,*) "tdk", dummy_r
        do ifmom=1, num_fmom
          if(mype==0) read(field_handle) fm_in
        end do
      end if
    end do  !while i.lt.. . .
    close(field_handle)
   end do
   deallocate(fm_in)
   deallocate(fm_in2)

  end subroutine read_fmom

  subroutine read_fmom_time(file_name,this_sparse_factor,read_number,s_index,n_time)
    implicit none
    character(len=128), intent(in) :: file_name
    character(len=128)  :: file_path_name(20)
    integer, intent(in) :: this_sparse_factor
    integer, intent(in) :: read_number
    integer :: field_handle,mom_handle
    complex, allocatable, dimension(:) :: fm_in
    integer :: numxyz
    integer :: num_fmom,ifmom
    integer :: i
    integer :: rnum
    integer, intent(out) :: s_index(SVD_n_restarts)  !start index
    integer, intent(out) :: n_time
    integer :: n_time_all(SVD_n_restarts)
    integer ::   ierr
    real :: dummy_r
    integer :: df_n_time(SVD_n_restarts)

    numxyz=nx0*nky0*nz0
    if(read_number==1) then
      num_fmom=n_fields
      df_n_time=SVD_df_n_time
    end if
    if(read_number==2) then
      num_fmom=n_moms
      df_n_time=SVD_df_n_time/SVD_istep_mf_factor
    end if
    if(read_number==3) then
      num_fmom=n_moms
      df_n_time=SVD_df_n_time/SVD_istep_mf_factor
    end if
    do i=1,20
      file_path_name(i)= trim(SVD_df_file_path)//'/'//trim(file_name)//&
             trim(SVD_df_file_suffix(i))
    end do
    do i=1,20
      file_path_name(i)= trim(SVD_df_file_path)//'/'//trim(file_name)//&
             trim(SVD_df_file_suffix(i))
    end do

    allocate(fm_in(numxyz))

    if(mype==0) then  !Calculate number of 
      do rnum=1,SVD_n_restarts
        !write(*,*) file_path_name(rnum)
        call get_unit_nr(field_handle)
#ifdef F2003_NO_OPEN_CONVERT
        open(field_handle,file=file_path_name(rnum),form='unformatted',status='unknown')
#else
        open(field_handle,file=file_path_name(rnum),form='unformatted',&
             &status='unknown',convert=SVD_endian)
#endif
        call get_unit_nr(mom_handle)
        s_index(rnum)=-1
        i=0
        do while(s_index(rnum)==-1.and.i.lt.df_n_time(rnum))
          i=i+1
          read(field_handle) dummy_r
          if(dummy_r.ge.SVD_start_time) s_index(rnum)=i
          do ifmom=1,num_fmom
            read(field_handle) fm_in
          end do
          
        end do
        if(s_index(rnum)==-1) s_index(rnum)=df_n_time(rnum)
        close(field_handle)
        n_time_all(rnum)=(df_n_time(rnum)-s_index(rnum)+1)/this_sparse_factor
        close(field_handle)
      end do

      n_time=sum(n_time_all)

      write(*,*) "df_n_time",df_n_time(1:SVD_n_restarts)
      write(*,*) "n_time_all",n_time_all
      write(*,*) "n_time=", n_time
    end if

    call MPI_BCAST(n_time,1,MPI_INTEGER,0,my_mpi_comm_world,ierr) 
    call MPI_BCAST(s_index,SVD_n_restarts,MPI_INTEGER,0,my_mpi_comm_world,ierr) 

   deallocate(fm_in)

  end subroutine read_fmom_time



subroutine comm_svd

  implicit none

  integer, dimension(:,:), allocatable :: blacs_map_temp
  integer :: i,j,ierr,dummy_i
  integer :: mat_arr_dims(2)

  call BLACS_PINFO(mypeb,procs_blacs)
  !write(*,*) "mypeb",mypeb

  n_row_procs=n_procs_sim
  n_col_procs=1

  if(n_procs_sim.ne.procs_blacs) then
    write(*,*) "Error BLACS procs not equal to MPI procs."
    write(*,*) "BLACS procs:", procs_blacs
    write(*,*) "MPI procs:", n_procs_sim
    stop
  end if

  !if(mype==0) write(*,*) "Creating mpi proc grid."
  mat_arr_dims=(/n_row_procs,n_col_procs/)
  call MPI_CART_CREATE(my_mpi_comm_world,2,mat_arr_dims,(/.false.,.false./),.false.,mpi_comm_mat,ierr)
  call MPI_CART_COORDS(mpi_comm_mat,mype,2,mype_mat,ierr)
  !write(*,*) "mype_mat",mype_mat
  !if(mype==0) write(*,*) "Creating blacs proc grid."
  !if(mype==0) write(*,*) "n_row_procs,n_col_procs",n_row_procs,n_col_procs
  allocate(blacs_map(0:n_row_procs-1,0:n_col_procs-1))
  allocate(blacs_map_temp(0:n_row_procs-1,0:n_col_procs-1))
  blacs_map=0
  myrow=mype_mat(1)
  mycol=mype_mat(2)
  do i=0,n_row_procs-1
    do j=0,n_col_procs-1
      if(myrow==i.and.mycol==j) then
         blacs_map(i,j)=mype
      end if
    end do
  end do
  call MPI_ALLREDUCE(blacs_map,blacs_map_temp,n_row_procs*n_col_procs&
      ,MPI_INTEGER,MPI_SUM,my_mpi_comm_world,ierr)
  blacs_map=blacs_map_temp
  !write(*,*) "blacs_map",blacs_map

  call BLACS_GET(dummy_i,0,my_context)
  !write(*,*) "after blacs_get, before blacs_gridmap,mype",mype
  call BLACS_GRIDMAP(my_context,blacs_map,n_row_procs,n_row_procs,n_col_procs)

  deallocate(blacs_map_temp)

end subroutine comm_svd

subroutine finalize_blacs

  call BLACS_EXIT(0)
  !if(mype==0) write(*,*) "Done."

end subroutine finalize_blacs

subroutine arrays_A

  use par_mod
  use mpi
  implicit none

  integer :: NUMROC,ierr

!Note:
!*  Alignment requirements
!*  ======================
!*
!*  The routine PCGESVD inherits the same alignement requirement as
!*  the routine PCGEBRD, namely:
!*
!*  The distributed submatrix sub( A ) must verify some alignment proper-
!*  ties, namely the following expressions should be true:
!*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
!*          where NB = MB_A = NB_A,
!*          IROFFA = MOD( IA-1, NB ), ICOFFA = MOD( JA-1, NB ),
!Note: The above is satisfied by making row and column blocking factors equal 
!==>pr0A=pc0A

  size0=min(num_row,num_col)  !number of singular values, etc.
  allocate(svals(size0))

  !!!!!!!!!!!!!!Set up A matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up A matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up A matrix!!!!!!!!!!!!!!!!!!!!

  nr0A=num_row    !Total number of rows (global)
  nc0A=num_col    !Total number of columns (global)

  !Processor block size: reasonable values seem to be 32 or 64, google for more details 
  pr0A=block_size
  pc0A=pr0a  !See alignment requirements above
  
  !Number of row / column elements for each processor
  lr0A=max(1,NUMROC(nr0A,pr0A,myrow,0,n_row_procs))
  lc0A=max(1,NUMROC(nc0A,pc0A,mycol,0,n_col_procs))
  allocate(Amat_loc(lr0A,lc0A))
  call DESCINIT(desc_A,nr0A,nc0A,pr0A,pc0A,0,0,my_context,lr0A,ierr)


end subroutine arrays_A


subroutine arrays_UV

  use par_mod
  use mpi
  implicit none

  integer :: NUMROC,ierr

!Note:
!*  Alignment requirements
!*  ======================
!*
!*  The routine PCGESVD inherits the same alignement requirement as
!*  the routine PCGEBRD, namely:
!*
!*  The distributed submatrix sub( A ) must verify some alignment proper-
!*  ties, namely the following expressions should be true:
!*  ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
!*          where NB = MB_A = NB_A,
!*          IROFFA = MOD( IA-1, NB ), ICOFFA = MOD( JA-1, NB ),
!Note: The above is satisfied by making row and column blocking factors equal 
!==>pr0A=pc0A

  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!

  nr0U=num_row    !Total number of rows (global)
  nc0U=size0    !Total number of columns (global)

  pr0U=block_size
  pc0U=pr0U
  
  !Number of row / column elements for each processor
  lr0U=max(1,NUMROC(nr0U,pr0U,myrow,0,n_row_procs))
  lc0U=num_col  !max(1,NUMROC(nc0U,pc0U,mycol,0,n_col_procs))
  !Lower and upper row and column bounds
  allocate(Umat_loc(lr0U,lc0U))
  call DESCINIT(desc_U,nr0U,nc0U,pr0U,pc0U,0,0,my_context,lr0U,ierr)

  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up U matrix!!!!!!!!!!!!!!!!!!!!
 

  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!

  nr0VT=size0    !Total number of rows (global)
  nc0VT=num_col    !Total number of columns (global)

  pr0VT=block_size
  pc0VT=pr0VT
  
  !Number of row / column elements for each processor
  lr0VT=size0   !max(1,NUMROC(nr0VT,pr0VT,myrow,0,n_row_procs))
  lc0VT=max(1,NUMROC(nc0VT,pc0VT,mycol,0,n_col_procs))
  !Lower and upper row and column bounds
  allocate(VTmat_loc(lr0VT,lc0VT))
  call DESCINIT(desc_VT,nr0VT,nc0VT,pr0VT,pc0VT,0,0,my_context,lr0VT,ierr)

  !write(*,*) "VTmat_loc info:"
  !write(*,*) "mype,nr,nc,lr,lc",mype,nr0VT,nc0VT,lr0VT,lc0VT

  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!Set up VT matrix!!!!!!!!!!!!!!!!!!!!

end subroutine arrays_UV

subroutine svd_deallocate_arrays
  use par_mod

  if(allocated(Amat_loc)) deallocate(Amat_loc)
  if(allocated(Umat_loc)) deallocate(Umat_loc)
  if(allocated(VTmat_loc)) deallocate(VTmat_loc)
  if(allocated(mat_glob)) deallocate(mat_glob)
  if(allocated(svals)) deallocate(svals)
  if(allocated(blacs_map)) deallocate(blacs_map)

end subroutine svd_deallocate_arrays

subroutine get_svd

  !!!!!!!!!!P?GESVD parameters!!!!!!!!!!!!!!!!
  integer :: info
  integer :: df_handle
  character(len=1) :: jobu, jobvt
  !integer :: m,n = num_row,num_col
  !a=Aloc_mat
  integer :: ia,ja
  !desca=desc_A
  !s=svals
  !u=Umat_loc
  !vt=VTmat_loc
  integer :: iu,ju
  integer :: ivt,jvt
  complex, allocatable, dimension(:) :: work
  integer :: lwork
  real, allocatable, dimension(:) :: rwork 
  !!!!!!!!!!P?GESVD parameters!!!!!!!!!!!!!!!!
  complex :: dummy_work,Amat_dummy
  real :: dummy_rwork
  integer :: i
 !mp = number of local rows in A and U     ===> lp0A,lp0U
 !nq = number of local columns in A and VT ===>lq0A,lq0VT
 !size = min(m, n)                         ===>size0
 !sizeq = number of local columns in U     ===>lq0U
 !sizep = number of local rows in VT       ===>lp0VT

  jobu='V'
  jobvt='V'
  ia=1
  ja=1
  iu=1
  ju=1
  ivt=1
  jvt=1

 !Calculate lwork:
  !if(mype==0) write(*,*) "Calculating sufficient lwork"
  lwork=-1  !This singals the routine to calculate the appropriate lwork.
  if(mype==0) write(*,*) "Calculating lwork."
  if(prec=='double'.or.prec=='DOUBLE') then
    call PZGESVD(jobu,jobvt,nr0A,nc0A,Amat_dummy,ia,ja,desc_A,svals,Umat_loc,& 
         iu,ju,desc_U,VTmat_loc,ivt,jvt,desc_VT,dummy_work,lwork,dummy_rwork,info)
  else
    call PCGESVD(jobu,jobvt,nr0A,nc0A,Amat_dummy,ia,ja,desc_A,svals,Umat_loc,&
         iu,ju,desc_U,VTmat_loc,ivt,jvt,desc_VT,dummy_work,lwork,dummy_rwork,info)
  end if 
  if(info==0) then
     lwork=nint(real(dummy_work))
     if(mype==0) write(*,*) "lwork",lwork
     if(mype==0) write(*,*) "size of rwork",nint(dummy_rwork)
     allocate(work(lwork))
     allocate(rwork(nint(dummy_rwork)))
     if(mype==0) write(*,*) "Done allocating work arrays."
  else
     if(mype==0) write(*,*) "Error calculating necessary size of rwork and lwork."
     stop 
  endif
 
  if(prec=='double'.or.prec=='DOUBLE') then
    call PZGESVD(jobu,jobvt,nr0A,nc0A,Amat_loc,ia,ja,desc_A,svals,Umat_loc,&
         iu,ju,desc_U,VTmat_loc,ivt,jvt,desc_VT,work,lwork,rwork,info)
  else
    call PCGESVD(jobu,jobvt,nr0A,nc0A,Amat_loc,ia,ja,desc_A,svals,Umat_loc,&
         iu,ju,desc_U,VTmat_loc,ivt,jvt,desc_VT,work,lwork,rwork,info)
  end if 

  call get_unit_nr(df_handle)
  if(info==0) then
    if(mype==0) write(*,*) "SVD successfully calculated!"
    if(mype==0) write(*,*) "max,min,max / min",svals(1),svals(size0),svals(1)/svals(size0)
    if(mype==0) then
      open(df_handle,file=trim(diagdir)//'/svals_fmom',status='unknown')
      !open(unit=100,file=trim(diagdir)//'/svals.dat',status='unknown')
      do i=1,size0
        write(df_handle,*) i,svals(i)
      end do
      close(df_handle)
    end if
    call my_barrier()

  else
    if(mype==0) write(*,*) "Error in SVD routine, info=",info  
    if(mype==0.and.info==-19) write(*,*) "Illegal lwork.  Most likely a memory problem."
    if(mype==0) write(*,*) "max,min,max / min",svals(1),svals(size0),svals(1)/svals(size0)
    if(mype==0) then
      !open(unit=100,file=trim(diagdir)//'/svals.dat',status='unknown')
      !do i=1,size0
      !  write(100,*) i,svals(i)
      !end do
      !close(100)
    end if
    stop
  end if

  deallocate(work)
  deallocate(rwork)

end subroutine get_svd



  subroutine integral_weight_fm(fm_in,fm_out)
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1), intent(in) :: fm_in
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1), intent(out) :: fm_out
    integer :: k
    complex, dimension(0:nx0-1,0:nky0-1) :: one_array

    if(xy_local) then
      one_array=cmplx(1.0,0.0) 
      !For now g_1_out will be used to store the integral weighting
      do k=0,nz0-1
              fm_out(:,:,k)=one_array(:,:)*sqrt(jac_glob(k))
      enddo
      fm_out=fm_out/sqrt((real(nz0)*geom%avg_jaco))/2.0
      fm_out=fm_in*fm_out
      fm_out(:,0,:) = 0.5*fm_out(:,0,:)
    else
      if(mype==0)  write(*,*) "!!!!!!!Must implement integral_weight_fm for nonlocal runs."
      fm_out=fm_in
    end if

  end subroutine integral_weight_fm

  subroutine integral_unweight_fm(fm_in,fm_out)
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1), intent(in) :: fm_in
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1), intent(out) :: fm_out
    integer :: k
    complex, dimension(0:nx0-1,0:nky0-1) :: one_array

    if(xy_local) then
      one_array=cmplx(1.0,0.0) 
      !For now g_1_out will be used to store the integral weighting
      do k=0,nz0-1
        fm_out(:,:,k)=one_array(:,:)/sqrt(jac_glob(k))
      enddo
      fm_out=fm_out*sqrt((real(nz0)*geom%avg_jaco))*2.0
      fm_out=fm_in*fm_out
      fm_out(:,0,:) = 2.0*fm_out(:,0,:)
    else
      if(mype==0)  write(*,*) "!!!!!!!Must implement integral_weight_fm for nonlocal runs."
      fm_out=fm_in
    end if

  end subroutine integral_unweight_fm

  !subroutine get_indices(vi,i,j,k,l,m,n,n1,n2,n3,n4,n5)
!
!    implicit none
!
!    integer, intent(in) :: vi,n1,n2,n3,n4,n5
!    integer, intent(out) :: i,j,k,l,m,n 
!    integer :: lt6,lt5,lt4,lt3,lt2
!    integer :: v2
!
!    v2=vi-1
!    lt6=n5*n4*n3*n2*n1
!    lt5=n4*n3*n2*n1
!    lt4=n3*n2*n1
!    lt3=n2*n1
!    lt2=n1
!
!    n=v2/lt6+1
!    m=(v2-(n-1)*lt6)/lt5+1
!    l=(v2-(n-1)*lt6-(m-1)*lt5)/lt4+1
!    k=(v2-(n-1)*lt6-(m-1)*lt5-(l-1)*lt4)/lt3+1
!    j=(v2-(n-1)*lt6-(m-1)*lt5-(l-1)*lt4-(k-1)*lt3)/lt2+1
!    i=(v2-(n-1)*lt6-(m-1)*lt5-(l-1)*lt4-(k-1)*lt3-(j-1)*lt2)+1
!
!    n=n-1
!    m=m-1
!    l=l-1
!    k=k-1
!    j=j-1
!    i=i-1
!
!  end subroutine get_indices

  !subroutine get_df_proc(i1,i2,i3,i4,i5,i6,df_proc)
!
!    integer, intent(in) :: i1,i2,i3,i4,i5,i6
!    integer, intent(out) :: df_proc
!    integer :: ierr
!    integer :: proc_temp
!
!    proc_temp=0
!    if(i1.ge.li1.and.i1.le.li2&
!       .and.i2.ge.lj1.and.i2.le.lj2&
!       .and.i3.ge.lk1.and.i3.le.lk2&
!       .and.i4.ge.ll1.and.i4.le.ll2&
!       .and.i5.ge.lm1.and.i5.le.lm2&
!       .and.i6.ge.ln1.and.i6.le.ln2) then
!       proc_temp=mype
!    end if
!    
!    call MPI_ALLREDUCE(proc_temp,df_proc,1&
!          ,MPI_INTEGER,MPI_SUM,my_mpi_comm_world,ierr)
!
!    !write(*,*) "mype,df_proc",mype,df_proc
!
!  end subroutine get_df_proc


  subroutine output_data_parallel_fm

    complex :: field_out(nx0*nky0*nz0)
    complex :: tpar_out(nx0*nky0*nz0)
    complex :: tperp_out(nx0*nky0*nz0)
    complex :: fm_temp(nx0*nky0*nz0)
    complex :: zeros(nx0*nky0*nz0)
    integer :: q,p,lr,lc,pr,pc
    integer :: INDXG2L,INDXG2P,ierr
    integer :: field_handle,mom_handle,numxyz,i,df_handle
    complex :: vt_out(size0) 
    complex :: vt_out2(size0) 

    numxyz=nx0*nky0*nz0

    call get_unit_nr(field_handle)
    if(mype==0) open(field_handle,file=trim(diagdir)//'/field_svd',&
         &form='unformatted',status='unknown')
    call get_unit_nr(mom_handle)
    if(mype==0) open(mom_handle,file=trim(diagdir)//'/mom_'//trim(spec(0)%name)//'_svd',&
         &form='unformatted',status='unknown')

    zeros=0.0
    do q=1,size0
      field_out=cmplx(0.0,0.0)
      tpar_out=cmplx(0.0,0.0)
      tperp_out=cmplx(0.0,0.0)
      do p=1,num_row  !Get local processor and indices from global indices
        !Find the local matrix indices corresponding to the global indices
        lr=INDXG2L(p,pr0A,0,0,n_row_procs)
        lc=INDXG2L(q,pc0A,0,0,n_col_procs)
        !Find the processors corresponding to these global indices
        pr=INDXG2P(p,pr0A,0,0,n_row_procs)
        pc=INDXG2P(q,pc0A,0,0,n_col_procs)
        !call get_indices(p,i1,i2,i3,i4,i5,i6,nx0,nky0,nz0,nv0,nw0)
        if(myrow==pr.and.mycol==pc.and.p.le.numxyz) then
          field_out(p)=Umat_loc(lr,lc)
        elseif(myrow==pr.and.mycol==pc.and.(p.le.2*numxyz.and.p.gt.numxyz)) then
          tpar_out(p-numxyz)=Umat_loc(lr,lc)
        elseif(myrow==pr.and.mycol==pc.and.p.gt.2*numxyz) then
          tperp_out(p-2*numxyz)=Umat_loc(lr,lc)
        end if
      end do
      time=real(q)
      call MPI_ALLREDUCE(field_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,field_out)
      call MPI_ALLREDUCE(tpar_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,tpar_out)
      call MPI_ALLREDUCE(tperp_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,tperp_out)
      if(mype==0) write(field_handle) real(q)
      if(mype==0) write(field_handle) field_out
      do i=2,n_fields
        if(mype==0) write(field_handle) zeros
      end do

      if(mype==0) write(mom_handle) real(q)
      do i=1,n_moms
        if(mype==0.and.i==2) then
          write(mom_handle) tpar_out
        elseif(mype==0.and.i==3) then
          write(mom_handle) tperp_out
        else
          if(mype==0) write(mom_handle) zeros
        end if
      end do

    end do

    close(field_handle)
    close(mom_handle)

    call get_unit_nr(df_handle)
    if(mype==0) open(df_handle,file=trim(diagdir)//'/c_tot_svd.dat',&
         form='unformatted',status='unknown')
    do q=1,num_col
      vt_out=cmplx(0.0,0.0)
      do p=1,size0
        !Find the local matrix indices corresponding to the global indices
        lr=INDXG2L(p,pr0A,0,0,n_row_procs)
        lc=INDXG2L(q,pc0A,0,0,n_col_procs)
        !Find the processors corresponding to these global indices
        pr=INDXG2P(p,pr0A,0,0,n_row_procs)
        pc=INDXG2P(q,pc0A,0,0,n_col_procs)
        !call get_indices(p,i1,i2,i3,i4,i5,i6,nx0,nky0,nz0,nv0,nw0)
        if(myrow==pr.and.mycol==pc) then
          vt_out(p)=svals(p)*VTmat_loc(lr,lc)
        end if
      end do
      call MPI_ALLREDUCE(vt_out,vt_out2,size0&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      
      if(mype==0) write(df_handle) time_arr(q)
      if(mype==0) write(df_handle) vt_out2(:)

    end do
    close(df_handle)

  !nr0VT=size0    !Total number of rows (global)
  !nc0VT=num_col    !Total number of columns (global)

    !if(mype==0) open(df_handle,file=trim(diagdir)//'/c_tot_svd.dat',&
    !       &form='unformatted',status='unknown')
    !do i=1,n_time
    !  write(df_handle) time_arr(i)
    !  write(df_handle) vt(:,i)
    !end do
    !close(df_handle)
    !write(*,*) "Done."

  end subroutine output_data_parallel_fm

  subroutine output_data_parallel_fm_beta

    complex, allocatable :: field_out(:)
    complex, allocatable :: tpar_out(:)
    complex, allocatable :: tperp_out(:)
    !!!!!Testing purposes!!!!!!!!
    !!!!!Testing purposes!!!!!!!!
    !!!!!Testing purposes!!!!!!!!
    !complex :: field_sum(nx0*nky0*nz0,size0)
    !complex :: tpar_sum(nx0*nky0*nz0,size0)
    !complex :: tperp_sum(nx0*nky0*nz0,size0)
    !!!!!Testing purposes!!!!!!!!
    !!!!!Testing purposes!!!!!!!!
    !!!!!Testing purposes!!!!!!!!
    complex :: fm_temp(nx0*nky0*nz0)
    !complex :: zeros(nx0*nky0*nz0)
    integer :: q,p,lr,lc,pr,pc
    integer :: INDXG2L,INDXG2P,ierr
    integer :: field_handle,mom_handle,numxyz,i,df_handle
    integer :: nrg_handle
    !real :: nrg_test(num_col,nrgcols)
    real :: nrg_out(nrgcols)
!    real :: nrg_temp

    write(*,*) "In output_data_parallel_fm_beta"

    numxyz=nx0*nky0*nz0
    allocate(field_out(numxyz))
    allocate(tpar_out(numxyz))
    allocate(tperp_out(numxyz))
    nrg_out=0.0

    call get_unit_nr(field_handle)
    if(mype==0) open(field_handle,file=trim(diagdir)//'/field_svd',&
         &form='unformatted',status='unknown')
    call get_unit_nr(mom_handle)
    if(mype==0) open(mom_handle,file=trim(diagdir)//'/mom_'//trim(spec(0)%name)//'_svd',&
         &form='unformatted',status='unknown')
    call get_unit_nr(df_handle)
    if(mype==0) open(df_handle,file=trim(diagdir)//'/c_tot_svd.dat',&
         form='unformatted',status='unknown')
    call get_unit_nr(nrg_handle)
    if(mype==0) open(nrg_handle,file=trim(diagdir)//'/nrg_svd',status='unknown')
    if(mype==0.and.xy_local) write(nrg_handle,*) "#Warning: Not completely debugged/verified!"

    !nrg_test=0.0

    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!
    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!
    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!
    !field_sum=cmplx(0.0,0.0)
    !tpar_sum=cmplx(0.0,0.0)
    !tperp_sum=cmplx(0.0,0.0)
    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!
    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!
    !!!!!!!!!!!!Testing purposes!!!!!!!!!!!!

    do q=1,size0
      field_out=cmplx(0.0,0.0)
      tpar_out=cmplx(0.0,0.0)
      tperp_out=cmplx(0.0,0.0)
      do p=1,num_row  !Get local processor and indices from global indices
        !Find the local matrix indices corresponding to the global indices
        lr=INDXG2L(p,pr0A,0,0,n_row_procs)
        lc=INDXG2L(q,pc0A,0,0,n_col_procs)
        !Find the processors corresponding to these global indices
        pr=INDXG2P(p,pr0A,0,0,n_row_procs)
        pc=INDXG2P(q,pc0A,0,0,n_col_procs)
        !call get_indices(p,i1,i2,i3,i4,i5,i6,nx0,nky0,nz0,nv0,nw0)
        if(myrow==pr.and.mycol==pc.and.p.le.numxyz) then
          field_out(p)=Umat_loc(lr,lc)
        elseif(myrow==pr.and.mycol==pc.and.(p.le.2*numxyz.and.p.gt.numxyz)) then
          tpar_out(p-numxyz)=Umat_loc(lr,lc)
        elseif(myrow==pr.and.mycol==pc.and.p.gt.2*numxyz) then
          tperp_out(p-2*numxyz)=Umat_loc(lr,lc)
        end if
      end do

      time=real(q)
      call MPI_ALLREDUCE(field_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,field_out)
      call MPI_ALLREDUCE(tpar_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,tpar_out)
      call MPI_ALLREDUCE(tperp_out,fm_temp,numxyz&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
      call integral_unweight_fm(svals(q)*fm_temp,tperp_out)

      if(xy_local) then
        call get_QES_from_svd(field_out,tpar_out,tperp_out,nrg_out(7))
        call fm_sum_3d(tpar_out,nrg_out(3))
        call fm_sum_3d(tperp_out,nrg_out(4))
        if(mype==0) write(nrg_handle,'(1f12.6)') real(q)
        if(mype==0) write(nrg_handle,'(10es12.4)') nrg_out
      else
        if(mype==0) write(*,*) "Implement this for global."
      end if

      if(mype==0) write(field_handle) real(q)
      if(mype==0) write(field_handle) field_out
      field_out=0.0
      do i=2,n_fields
        if(mype==0) write(field_handle) field_out !zeros
      end do

      if(mype==0) write(mom_handle) real(q)
      do i=1,n_moms
        if(mype==0.and.i==2) then
          write(mom_handle) tpar_out
        elseif(mype==0.and.i==3) then
          write(mom_handle) tperp_out
        else
          if(mype==0) write(mom_handle) field_out !zeros
        end if
      end do

      if(mype==0) write(df_handle) real(q)
      if(mype==0) write(df_handle) VTmat_loc(q,:) 

      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !do p=1,num_col
      !  field_sum(:,p)=field_sum(:,p)+field_out(:)*VTmat_loc(q,p)
      !  tpar_sum(:,p)=tpar_sum(:,p)+tpar_out(:)*VTmat_loc(q,p)
      !  tperp_sum(:,p)=tperp_sum(:,p)+tperp_out(:)*VTmat_loc(q,p)
      !end do
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!

    end do

      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !if(mype==0) then
      !  write(*,*) "Testing SVD recomposition."
      !  open(unit=888,file=trim(SVD_df_file_path)//'/'//'mom_ions.dat',form='unformatted',status='unknown')
      !  do p=1,num_col
      !     read(888) nrg_temp 
      !     write(*,*) nrg_temp, time_arr(p)
      !     read(888) field_out
      !     read(888) tpar_out
      !     read(888) tperp_out
      !     read(888) field_out
      !     read(888) field_out
      !     read(888) field_out
      !     write(*,*) "tpar diff",abs(sum(tpar_sum(:,p)-tpar_out,1))/abs(sum(tpar_out))
      !     write(*,*) "tperp diff",abs(sum(tperp_sum(:,p)-tperp_out,1))/abs(sum(tperp_out))
      !  end do
      !  close(888)
      !end if
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!

      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!Now calculate heat fluxes as a test
      !nrg_test=0.0
      !do p=1,num_col
      !  if(xy_local) call get_QES_from_svd(field_sum(:,p),tpar_sum(:,p),tperp_sum(:,p),QES)
      !  nrg_test(p,7)=QES
      !  if(xy_local) call fm_sum_3d(tpar_sum(:,p),nrg_temp)
      !  nrg_test(p,3)=nrg_temp
      !  if(xy_local) call fm_sum_3d(tperp_sum(:,p),nrg_temp)
      !  nrg_test(p,4)=nrg_temp
      !end do
      !call get_unit_nr(nrg_handle)
      !if(mype==0) open(nrg_handle,file=trim(diagdir)//'/nrg_test.dat',status='unknown')
      !do i=1,num_col
      !  write(nrg_handle,'(1f12.6)') time_arr(i)
      !  write(nrg_handle,'(10es12.4)') nrg_test(i,:)
      !end do
      !close(nrg_handle)
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!
      !!!!!!!!!!Testing purposes!!!!!!!!!

    !if(mype==0) write(*,*) field_handle,mom_handle,df_handle,nrg_handle
    close(field_handle)
    close(mom_handle)
    close(df_handle)
    close(nrg_handle)

    deallocate(field_out)
    deallocate(tpar_out)
    deallocate(tperp_out)

  end subroutine output_data_parallel_fm_beta

  subroutine get_QES_from_svd(phi,tpar,tperp,QES_out)
    implicit none
    complex, intent(in) :: phi(0:nx0-1,0:nky0-1,0:nz0-1)
    complex, intent(in) :: tpar(0:nx0-1,0:nky0-1,0:nz0-1)
    complex, intent(in) :: tperp(0:nx0-1,0:nky0-1,0:nz0-1)
    complex  :: dphidy(0:nx0-1,0:nky0-1,0:nz0-1)
    complex  :: QES(0:nx0-1,0:nky0-1,0:nz0-1)
    integer :: j
    real, intent(out) :: QES_out

    do j=0,nky0-1
      dphidy(:,j,:)=phi(:,j,:)*imag*kymin*j
    end do
    QES=-conjg(tpar+tperp)*dphidy

    if (yx_order) then
       if ((xy_local).and.(evenx.eq.1)) QES(:,hkx+1,:) = 0.0
    else
       if ((xy_local).and.(evenx.eq.1)) QES(hkx+1,:,:) = 0.0
    endif

    QES(:,0,:) = 0.5*QES(:,0,:)

    QES_out = 2.0*real(sum(sum(sum(QES,1),1)*jac_glob(:))/&
            (real(nz0)*geom%avg_jaco))

  end subroutine get_QES_from_svd

  subroutine fm_sum_3d(fm_in,sum_out)
    implicit none
    complex, intent(in) :: fm_in(0:nx0-1,0:nky0-1,0:nz0-1)
    complex, dimension(:,:,:), allocatable :: fm_temp
    real, intent(out) :: sum_out

    allocate(fm_temp(0:nx0-1,0:nky0-1,0:nz0-1))
    fm_temp=fm_in
    if (yx_order) then
       stop "this is not suitable for yx_order."
    else
       !if ((xy_local).and.(evenx.eq.1)) fm_temp(hkx+1,:,:) = 0.0
       fm_temp(:,0,:) = 0.5*fm_temp(:,0,:)
    endif

    sum_out = real(2.0*sum(sum(sum(conjg(fm_temp)*fm_temp,1),1)*jac_glob(:))/&
            (real(nz0)*geom%avg_jaco))

    !sum_out = real(2.0*sum(sum(sum(conjg(fm_temp)*fm_temp,1),1)*geom%jacobian(pi1,pj1,:))/&
    !        (real(nz0)*geom%avg_jaco))

    deallocate(fm_temp)

  end subroutine fm_sum_3d

#endif
#endif

End Module diagnostics_fmsvd

