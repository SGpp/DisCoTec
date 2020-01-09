#include "redef.h"
#include "intrinsic_sizes.h"
!> This modules contains routines for SVD diagnostics 
Module diagnostics_svd
  Use par_mod 
  Use file_io, only: get_unit_nr
  Use vel_space, only: fm, mat_00
  use geometry
  use diagnostics, only: exec_all_diags 
  use communications
  use mpi
  use arrays
  
  Implicit None
  
  PUBLIC :: svd_projection
  
  PRIVATE 
  
  !For parallel svd:
  
  character(len=16) :: faccess='sequential'
#ifdef WITHSCAL
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
  real,allocatable, dimension(:,:,:) :: fm_glob
  real,allocatable, dimension(:,:,:) :: mat_00_glob
  real,allocatable, dimension(:) :: jac_glob
#endif 
  
Contains
  
!!!******************************************************************!!!
  
  subroutine svd_projection
    
    if (dfout_mpiio.and..not.SVD_field_mom) then
       faccess='stream'
    else
       faccess='sequential'
    endif
    
    if(svd_parallel) then
#ifdef WITHSCAL
       call initialize_svd_parallel
       call svd_projection_parallel
       call finalize_svd_parallel
#else
       if(mype==0) write(*,*) "Need scalapack to use parallel svd routine!"
       stop
#endif
    else
       call svd_projection_serial
    end if
    
  end subroutine svd_projection
  
  !>This subroutine calculates the SVD of distribution function data output by GENE.
  !This routine can only be run serially.
  !The input distribution functions come from the diag_df_ev routine (triggered by istep_dfout).
  !Instructions: Use the parameters.dat file from a previous simulation and modify:
  !diagdir, all n_procs=1, plus the following:
  !Input parameters (in extended_diags namelist):
  !SVD_proj=T
  !SVD_f0_flag (default:T): normalizes to 1/F0 to be more similar to free energy
  !SVD_kx_ind:  central kx index to analyze
  !SVD_ky_ind:  ky index to analyze
  !SVD_nkx0: number of kx modes including connections (overwrites nx0)
  !SVD_df_n_time: specify number of time steps to analyze for each run
  !SVD_start_time (default:0) : start time for analysis
  !SVD_sparse_factor (default:1): analysis takes every SVD_sparse_factor'th time step
  !SVD_df_file_name: default constructs file name from SVD_kx_ind and SVD_ky_ind
  !SVD_df_file_path: default is diagdir
  !SVD_n_restarts: number of restarts for df files (default = 1, i.e. one run w/o restarts)
  !SVD_df_file_suffix (default:'.dat'): array of file endings for a series of runs
  !SVD_parallel: uses parallel SVD for large datasets (see svd_projection_parallel)
  !Note: new kymin is constructed from old kymin and SVD_ky_ind
  !Note: kx_center is constructed from old lx and SVD_kx_ind
  !Note: parallel execution is now available with SCALAPACK
  subroutine svd_projection_serial
    complex, allocatable, dimension(:,:) :: df_matrix
    complex, allocatable, dimension(:,:) :: df_matrix2
    !complex :: df(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    complex :: df2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
    character(len=1) :: jobu,jobvt
    integer :: info,lda,ldvt,lwork,m,n ,i,ldu,df_handle,j
    real :: dummy_r
    real, allocatable, dimension(:) :: time_array
    real, allocatable, dimension(:) :: rwork
    complex, allocatable, dimension(:,:) :: U_mat, VT
    complex, allocatable, dimension(:) :: work
    complex :: dummy_c,dc2
    complex, allocatable, dimension(:,:) :: c_tot
    complex, allocatable, dimension(:,:) :: c_svd
    complex :: df_sum2(nx0*nz0*nv0*nw0*n_spec)
    !complex :: df3d(li1:li2,lj1:lj2,lk1:lk2)
    integer :: i0
    logical :: dummyl
    real :: maxsv!,stddev!,error_avg
    real, allocatable, dimension(:) :: singular_values
    logical :: error_check
    integer :: rnum
    integer :: s_index(SVD_n_restarts)
    integer :: n_time_all(SVD_n_restarts),n_time
    
    error_check=.false.
    
    do rnum=1,SVD_n_restarts
       call get_unit_nr(df_handle)
       open(df_handle,file=trim(SVD_df_file_path)//'/'//trim(SVD_df_file_name)//&
            trim(SVD_df_file_suffix(rnum)),form='unformatted',status='unknown',&
            access=faccess)      
       s_index(rnum)=-1
       i=0
       do while(s_index(rnum)==-1.and.i.lt.SVD_df_n_time(rnum))
          i=i+1
          read(df_handle) dummy_r
          if(dummy_r.ge.SVD_start_time) s_index(rnum)=i
          read(df_handle) df2
       end do
       if(s_index(rnum)==-1) s_index(rnum)=SVD_df_n_time(rnum)
       close(df_handle)
       n_time_all(rnum)=(SVD_df_n_time(rnum)-s_index(rnum)+1)/SVD_sparse_factor
       close(df_handle)
    end do
    
    n_time=sum(n_time_all)
    write(*,*) "SVD_df_n_time",SVD_df_n_time(1:SVD_n_restarts)
    write(*,*) "n_time_all",n_time_all
    write(*,*) "n_time=", n_time
    allocate(df_matrix(nx0*nz0*nw0*nv0*n_spec,n_time))
    allocate(df_matrix2(nx0*nz0*nw0*nv0*n_spec,n_time))
    allocate(time_array(n_time))
    
    allocate(c_tot(n_time,n_time))
    allocate(c_svd(n_time,n_time))
    allocate(singular_values(n_time))
    
    write(*,*) "Beginning SVD analysis."
    !The size of the left singular vectors
    m=nx0*nz0*nv0*nw0*n_spec
    !The size of the right singular vectors
    n=n_time
    
    i0=0
    do rnum=1,SVD_n_restarts
       call get_unit_nr(df_handle)
       open(df_handle,file=trim(SVD_df_file_path)//'/'//trim(SVD_df_file_name)// &
            trim(SVD_df_file_suffix(rnum)),form='unformatted',status='unknown',&
            access=faccess)
       i=0
       !Read in distribution function from earlier nonlinear simulation
       do while(i.lt.SVD_df_n_time(rnum))
          i=i+1 
          if(mod(i,SVD_sparse_factor)==0.and.i.ge.s_index(rnum).and.i0.lt.n_time) then
             i0=i0+1
             read(df_handle) time_array(i0)
             read(df_handle) df_matrix(:,i0)
          else
             read(df_handle) dummy_r
             read(df_handle) df2
          end if
       end do
       close(df_handle)
    end do
    
    if(i0.ne.n_time) write(*,*) "i0 does not equal n_time!!!!"
    write(*,*) "Weighting the vectors so that scalar product = integral."
    !This gives more physical meaning to the POD procedure
    do i=1,n_time
       call integral_weight_serial(df_matrix(:,i),df_sum2)
       df_matrix(:,i)=df_sum2
    end do
    write(*,*) "Done."
    
!!!!!!!!!!!!!!!!ZGESVD!!!!!!!!!!!!!!!!!
    !Decomposes A=U*SIGMA*conjugate-transpose(V)
    !A=df_matrix2
    !U=u_mat
    !S=singular_values
    !VT=vt
!!!!!!!!!!!!!!!!ZGESVD!!!!!!!!!!!!!!!!!
    
    allocate(u_mat(m,n)) 
    jobu='s'   !the first min(m,n) columns of U (the left singular vectors) are returned in U
    jobvt='s'  !the first min(m,n) rows of v**H (the right singular vectors) are returned in the array VT
    lda=m      !the leading dimension of the array A
    ldu=m      !the leading dimension of the array U
    ldvt=n     !the leading dimension of the array VT
    allocate(vt(n,n)) 
    lwork=100*( 2*n+m )  !>=2*min*(M,N)+MAX(M,N)
    allocate(work(lwork))   
    allocate(rwork(5*m))
    write(*,*) "Calculating svd."  
    write(*,*) "n=",n
    !put df_matrix into df_matrix2 since it will be destroyed in zgesvd
    df_matrix2=df_matrix
    !Need cgesvd for single precision and zgesvd for double precision
    
    if(prec.eq.'double'.or.prec.eq.'DOUBLE') then
       call zgesvd(jobu,jobvt,m,n,df_matrix2,lda,singular_values,u_mat,ldu,&
            vt,ldvt,work,lwork,rwork,info)
    else
       call cgesvd(jobu,jobvt,m,n,df_matrix2,lda,singular_values,u_mat,ldu,&
            vt,ldvt,work,lwork,rwork,info)
    end if
    
    deallocate(df_matrix2)
    if(info==0) then
       write(*,*) "Singular values successfully calculated!" 
       if(minval(singular_values).ne.0.0) write(*,*) "sigma_min/sigma_max",&
            minval(singular_values)/maxval(singular_values)
       !Adjust lwork in future if desired
       write(*,*) "lwork=",lwork,"Optimal lwork=",work(1)
    else
       write(*,*) "Error in cgesvd!"
       stop
    end if
    
    !calculate the standard gene output data for each of the POD modes
    close(df_handle)
    open(df_handle,file=trim(diagdir)//'/SVD_'//trim(SVD_df_file_name),&
         &form='unformatted',status='unknown')
    
    time=0.0
    do i=1,n
       !label the modes with integers just as in eigenvalue data output
       time=time+1.0
       !This puts u_mat into g_1 applying the appropriate parallel boundary condition, etc.
       call integral_unweight_serial(u_mat(:,i),df_sum2)
       call umat_to_g1(singular_values(i)*df_sum2)  !fills g_1
       call exec_all_diags(0,time,dummyl,dummyl,dummyl)
       write(df_handle) nint(time)
       write(df_handle) g_1
    end do
    close(df_handle)
    
    write(*,*) "Writing singular values to svals.dat"
    maxsv=maxval(singular_values)
    open(df_handle,file=trim(diagdir)//'/svals_'//trim(SVD_df_file_name),status='unknown')
    do i=1,n
       write(df_handle,*) i,singular_values(i)
    end do
    Write(*,*) "Done."
    
    write(*,*) "Weighting time vectors."
    do i=1,n
       do j=1,n
          vt(j,i)=singular_values(j)*vt(j,i)
       end do
    end do
    write(*,*) "Done."
    
    if(error_check) then
       write(*,*) "Checking errors."
       do i=1,n_time
          if(mod(i,50)==0) write(*,*) i, " of ", n_time
          df_sum2=cmplx(0.0,0.0)
          do j=1,n
             call scalar_product_serial(u_mat(:,j),df_matrix(:,i),c_tot(j,i))
             df_sum2=df_sum2+c_tot(j,i)*u_mat(:,j)
          end do
          
          call scalar_product_serial(df_matrix(:,i)-df_sum2,df_matrix(:,i)-df_sum2,dummy_c)
          call scalar_product_serial(df_matrix(:,i),df_matrix(:,i),dc2)
          if(abs(real(conjg(dc2)*dc2)).ge.1.0e-14) dummy_c=dummy_c/dc2
          write(*,*) "Error at time step ", i, " . . . . . ", sqrt(real(conjg(dummy_c)*dummy_c))
       end do
    end if
      
    close(df_handle)
    
    !open(df_handle,file=trim(diagdir)//'/svd_vectors.dat',form='unformatted',status='unknown')
    !write(*,*) "Unweighting vectors and writing vectors to svd_vectors.dat"
    !do i=1,n
    !  df_sum2=u_mat(:,i)
    !  call integral_unweight(df_sum2,u_mat(:,i))
    !  write(df_handle) i
    !  write(df_handle) u_mat(:,i)
    !end do
    !close(df_handle)
    !write(*,*) "Done."
    
    write(*,*) "Writing c_tot_svd.dat."
    open(df_handle,file=trim(diagdir)//'/c_tot_svd.dat',form='unformatted',status='unknown')
    do i=1,n_time
       write(df_handle) time_array(i)
       write(df_handle) vt(:,i)
    end do
    close(df_handle)
    write(*,*) "Done."
    
  end subroutine svd_projection_serial
  
  subroutine integral_weight_serial(g_1_in,g_1_out)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: g_1_out
    integer :: k,l,m,n,pni
    complex, dimension(li1:li2,lj1:lj2) :: one_array
    
    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    if(evenx.ne.1) then
       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do k=lk1,lk2
             do m=lm1,lm2
                do l=ll1,ll2
                   g_1_out(:,:,k,l,m,n)=one_array(:,:)*sqrt(mat_00(pi1,pj1,k,l,m)*geom%jacobian(pi1gl,pj1,k))
                   if(SVD_f0_flag) then
                      g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)/sqrt(fm(pi1,pj1,k,l,m,pni))
                   end if
                enddo
             enddo
          enddo
       enddo
    else
       write(*,*) "Error.  evenx==1."
       stop
    end if
    
    g_1_out=g_1_out/sqrt((real(nz0)*geom%avg_jaco))
    g_1_out=g_1_in*g_1_out
    
    
  end subroutine integral_weight_serial
  
  subroutine integral_unweight_serial(g_1_in,g_1_out)
    
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in) :: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(out) :: g_1_out
    integer :: k,l,m,n,pni
    complex, dimension(li1:li2,lj1:lj2) :: one_array
    
    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    if(evenx.ne.1) then
       pni=pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do k=lk1,lk2
             do m=lm1,lm2
                do l=ll1,ll2
                   g_1_out(:,:,k,l,m,n)=one_array(:,:)/(sqrt(mat_00(pi1,pj1,k,l,m)*geom%jacobian(pi1gl,pj1,k)))
                   !var(o,n) = sum( sum( Sum(mom(:,:,:,o,n),1),1 )*geom%jacobian(pi1,pj1,:) )  
                   if(SVD_f0_flag) then
                      g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)*sqrt(fm(pi1,pj1,k,l,m,pni))
                   end if
                enddo
             enddo
          enddo
       enddo
    else
       write(*,*) "Error.  evenx==1."
       stop
    end if
    
    g_1_out=g_1_out*sqrt((real(nz0)*geom%avg_jaco))
    g_1_out=g_1_in*g_1_out
    
    
  end subroutine integral_unweight_serial
  
  subroutine umat_to_g1(u_mat_n)
    
    complex, intent(in) :: u_mat_n(nx0,nky0,nz0,nv0,nw0,n_spec)
    integer :: i
    !logical :: apply_pbc
    
    ! apply_pbc=.true.
    g_1(:,:,:,:,:,:)=u_mat_n(:,:,:,:,:,:)
    if(SVD_pbc_phase) then
       do i = 1,nx0
          if(mod(i,2)==0.and.i.le.(nx0/2+1)) g_1(i-1,:,:,:,:,:)=-1.0*g_1(i-1,:,:,:,:,:)
          if(mod(i,2)==1.and.i.gt.(nx0/2+1)) g_1(i-1,:,:,:,:,:)=-1.0*g_1(i-1,:,:,:,:,:)
       end do
    end if
    
  end subroutine umat_to_g1
  
  subroutine scalar_product_serial(v1,v2,sp)
    
    complex, intent(in), dimension(li0*lj0*lk0*ll0*lm0*ln0) :: v1,v2
    complex, intent(out) :: sp
    
    sp=sum(conjg(v1)*v2)
    
  end subroutine scalar_product_serial

#ifdef WITHSCAL
!!!!!!!!!!!!!!!!!!!Parallel!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Parallel!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!Parallel!!!!!!!!!!!!!!!!!!!!!
  subroutine initialize_svd_parallel
    
    real :: temp3d(0:nz0-1,0:nv0-1,0:nw0-1)
    real :: temp1d(0:nz0-1)
    integer :: ierr
    
    allocate(fm_glob(0:nz0-1,0:nv0-1,0:nw0-1))
    allocate(mat_00_glob(0:nz0-1,0:nv0-1,0:nw0-1))
    allocate(jac_glob(0:nz0-1))
    
    temp3d=0.0
    temp3d(lk1:lk2,ll1:ll2,lm1:lm2)=fm(pi1,pj1,lk1:lk2,ll1:ll2,lm1:lm2,pn1)
    call MPI_ALLREDUCE(temp3d,fm_glob,nz0*nv0*nw0,MPI_REAL_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
    temp3d=0.0
    temp3d(lk1:lk2,ll1:ll2,lm1:lm2)=mat_00(pi1,pj1,lk1:lk2,ll1:ll2,lm1:lm2)
    call MPI_ALLREDUCE(temp3d,mat_00_glob,nz0*nv0*nw0,MPI_REAL_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
    temp1d=0.0
    temp1d(lk1:lk2)=geom%jacobian(pi1,pj1,lk1:lk2)
    call MPI_ALLREDUCE(temp1d,jac_glob,nz0,MPI_REAL_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
    
  end subroutine initialize_svd_parallel
  
  subroutine finalize_svd_parallel
    deallocate(fm_glob)
    deallocate(mat_00_glob)
    deallocate(jac_glob)
  end subroutine finalize_svd_parallel
  
  !Everything carries over from svd_projection_serial
  !Activate with SVD_projection_parallel=T
  subroutine svd_projection_parallel
    
    integer :: df_handle
    complex, allocatable, dimension(:) :: df_in
    complex, allocatable, dimension(:) :: df_in2
    integer :: rnum
    integer :: s_index(SVD_n_restarts)
    integer :: n_time_all(SVD_n_restarts),n_time
    integer :: i,i0,p,q,lr,lc,pr,pc
    real :: dummy_r
    integer :: INDXG2L,  INDXG2P,ierr
    
    block_size=64
    num_row=nx0*nky0*nz0*nv0*nw0*n_spec
    allocate(df_in(num_row))
    allocate(df_in2(num_row))
    
    if(mype==0) then  !Calculate number of 
       do rnum=1,SVD_n_restarts
          call get_unit_nr(df_handle)
          open(df_handle,file=trim(SVD_df_file_path)//'/'//trim(SVD_df_file_name)//&
               trim(SVD_df_file_suffix(rnum)),form='unformatted',status='unknown',&
               access=faccess)
          s_index(rnum)=-1
          i=0
          do while(s_index(rnum)==-1.and.i.lt.SVD_df_n_time(rnum))
             i=i+1
             read(df_handle) dummy_r
             if(dummy_r.ge.SVD_start_time) s_index(rnum)=i
             read(df_handle) df_in
          end do
          if(s_index(rnum)==-1) s_index(rnum)=SVD_df_n_time(rnum)
          close(df_handle)
          n_time_all(rnum)=(SVD_df_n_time(rnum)-s_index(rnum)+1)/SVD_sparse_factor
          close(df_handle)
       end do
       
       n_time=sum(n_time_all)
       
       write(*,*) "SVD_df_n_time",SVD_df_n_time(1:SVD_n_restarts)
       write(*,*) "n_time_all",n_time_all
       write(*,*) "n_time=", n_time
    end if
    call MPI_BCAST(n_time,1,MPI_INTEGER,0,my_mpi_comm_world,ierr) 
    call MPI_BCAST(s_index,SVD_n_restarts,MPI_INTEGER,0,my_mpi_comm_world,ierr) 
    num_col=n_time
    allocate(time_arr(num_col))
    time_arr=0.0
    
    !call init_comm
    call comm_svd
    call arrays_A
    
    i0=0
    do rnum=1,SVD_n_restarts
       call get_unit_nr(df_handle)
       if(mype==0) open(df_handle,file=trim(SVD_df_file_path)//'/'//trim(SVD_df_file_name)// &
            trim(SVD_df_file_suffix(rnum)),form='unformatted',status='unknown',&
            access=faccess)
       i=0
       !Read in distribution function from earlier nonlinear simulation
       do while(i.lt.SVD_df_n_time(rnum))
          i=i+1 
          if(mod(i,SVD_sparse_factor)==0.and.i.ge.s_index(rnum).and.i0.lt.n_time) then
             i0=i0+1
             if(mype==0) read(df_handle) time_arr(i0)
             if(mype==0) read(df_handle) df_in
             !if(mype==0) write(*,*) "Calling integral_weight_serial2",mype
             if(mype==0) call integral_weight_serial2(df_in,df_in2)
             !if(mype==0) write(*,*) "After integral_weight_serial2",mype
             !write(*,*) "num_row,mype",num_row,mype
             call MPI_BCAST(df_in2,num_row,MPI_COMPLEX_TYPE,0,my_mpi_comm_world,ierr) 
             !write(*,*) "Filling Amat_loc.",mype
             q=i0
             do p=1,num_row  !Get local processor and indices from global indices
                !Find the local matrix indices corresponding to the global indices
                lr=INDXG2L(p,pr0A,0,0,n_row_procs)
                lc=INDXG2L(q,pc0A,0,0,n_col_procs)
                !Find the processors corresponding to these global indices
                pr=INDXG2P(p,pr0A,0,0,n_row_procs)
                pc=INDXG2P(q,pc0A,0,0,n_col_procs)
                if(myrow==pr.and.mycol==pc) then
                   Amat_loc(lr,lc)=df_in2(p)
                end if
             end do
             if(mype==0) write(*,*) "Done reading and distributing step:", i0
             
          else
             if(mype==0) read(df_handle) dummy_r
             if(mype==0) read(df_handle) df_in
          end if
       end do  !while i.lt.. . .
       close(df_handle)
    end do
    deallocate(df_in)
    
    !Broadcast time array
    !call MPI_BCAST(time_arr,num_col,MPI_REAL_TYPE,0,my_mpi_comm_world,ierr) 
    
    call arrays_UV
    if(mype==0) write(*,*) "Calculating SVD."
    call get_svd
    if(mype==0) write(*,*) "Done."
    call output_data_parallel 
    if(mype==0) write(*,*) "Deallocating arrays."
    call svd_deallocate_arrays
    if(mype==0) write(*,*) "Finalizing blacs."
      !call finalize_blacs
    if(mype==0) write(*,*) "Done with parallel svd."
    deallocate(time_arr)
    
  end subroutine svd_projection_parallel
  
  
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
          open(df_handle,file=trim(diagdir)//'/svals_'//trim(SVD_df_file_name),status='unknown')
          !open(unit=100,file=trim(diagdir)//'/svals.dat',status='unknown')
          do i=1,size0
             write(df_handle,*) i,svals(i)
          end do
          close(df_handle)
       end if
       
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
    
  subroutine integral_weight_serial2(g_1_in,g_1_out)
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1), intent(in) :: g_1_in
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1), intent(out) :: g_1_out
    integer :: k,l,m,n
    complex, dimension(0:nx0-1,0:nky0-1) :: one_array
    
    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    if(evenx.ne.1) then
       do n=0,n_spec-1
          do k=0,nz0-1
             do m=0,nw0-1
                do l=0,nv0-1
                   g_1_out(:,:,k,l,m,n)=one_array(:,:)*sqrt(mat_00_glob(k,l,m)*jac_glob(k))
                   if(SVD_f0_flag) then
                      g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)/sqrt(fm_glob(k,l,m))
                   end if
                enddo
             enddo
          enddo
       enddo
    else
       write(*,*) "Error.  evenx==1."
       stop
    end if
    
    g_1_out=g_1_out/sqrt((real(nz0)*geom%avg_jaco))/2.0
    g_1_out=g_1_in*g_1_out
    
    
  end subroutine integral_weight_serial2
  
  subroutine integral_unweight_serial2(g_1_in,g_1_out)
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1), intent(in) :: g_1_in
    complex, dimension(0:nx0-1,0:nky0-1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1), intent(out) :: g_1_out
    integer :: k,l,m,n
    complex, dimension(0:nx0-1,0:nky0-1) :: one_array
    integer :: i
    
    one_array=cmplx(1.0,0.0) 
    !For now g_1_out will be used to store the integral weighting
    if(evenx.ne.1) then
       do n=0,n_spec-1
          do k=0,nz0-1
             do m=0,nw0-1
                do l=0,nv0-1
                   g_1_out(:,:,k,l,m,n)=one_array(:,:)/(sqrt(mat_00_glob(k,l,m)*jac_glob(k)))
                   if(SVD_f0_flag) then
                      g_1_out(:,:,k,l,m,n)=g_1_out(:,:,k,l,m,n)*sqrt(fm_glob(k,l,m))
                   end if
                enddo
             enddo
          enddo
       enddo
    else
       write(*,*) "Error.  evenx==1."
       stop
    end if
    
    g_1_out=g_1_out*sqrt((real(nz0)*geom%avg_jaco))*2.0
    g_1_out=g_1_in*g_1_out
    
    if(SVD_pbc_phase) then
       do i = 1,nx0
          if(mod(i,2)==0.and.i.le.(nx0/2+1)) g_1_out(i-1,:,:,:,:,:)=-1.0*g_1_out(i-1,:,:,:,:,:)
          if(mod(i,2)==1.and.i.gt.(nx0/2+1)) g_1_out(i-1,:,:,:,:,:)=-1.0*g_1_out(i-1,:,:,:,:,:)
       end do
    end if
    
  end subroutine integral_unweight_serial2
  
  subroutine output_data_parallel
    
    complex :: df_out(num_row)
    complex :: df_out2(num_row)
    complex :: g_out(0:nx0-1,0:nky0-1,0:nz0-1,0:nv0-1,0:nw0-1,0:n_spec-1)
    integer :: q,p,lr,lc,pr,pc
    logical :: dummyl
    integer :: INDXG2L,INDXG2P,ierr
    integer :: df_handle
    complex :: vt_out(size0) 
    complex :: vt_out2(size0) 
    
    
    call get_unit_nr(df_handle)
    if(mype==0) open(df_handle,file=trim(diagdir)//'/SVD_'//trim(SVD_df_file_name),&
         &form='unformatted',status='unknown')
    
    do q=1,size0
       df_out=cmplx(0.0,0.0)
       do p=1,num_row  !Get local processor and indices from global indices
          !Find the local matrix indices corresponding to the global indices
          lr=INDXG2L(p,pr0A,0,0,n_row_procs)
          lc=INDXG2L(q,pc0A,0,0,n_col_procs)
          !Find the processors corresponding to these global indices
          pr=INDXG2P(p,pr0A,0,0,n_row_procs)
          pc=INDXG2P(q,pc0A,0,0,n_col_procs)
          !call get_indices(p,i1,i2,i3,i4,i5,i6,nx0,nky0,nz0,nv0,nw0)
          if(myrow==pr.and.mycol==pc) then
             df_out(p)=Umat_loc(lr,lc)
          end if
       end do
       time=real(q)
       call MPI_ALLREDUCE(df_out,df_out2,num_row&
            ,MPI_COMPLEX_TYPE,MPI_SUM,my_mpi_comm_world,ierr)
       call integral_unweight_serial2(svals(q)*df_out2,g_out)
       !call df_to_g(df_out2,g_out) 
       g_1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)=&
            g_out(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)
       call exec_all_diags(0,time,dummyl,dummyl,dummyl)
       if(mype==0) write(df_handle) q
       if(mype==0) write(df_handle) g_out
    end do
    
    close(df_handle)
    
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

  end subroutine output_data_parallel
  
#endif

End Module diagnostics_svd

