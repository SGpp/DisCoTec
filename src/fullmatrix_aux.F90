#include "redef.h"
!>Contains routines to build and manipulate an explicit representation of the linear gyrokinetic 
!!operator used in GENE
!!
!!All routines (and arrays) in this module are fully parallelized, nevertheless having an explicit
!!representation of L is only feasible for relatively small linear problems due to memory constraints
!!(and it's never efficient)
module fullmatrix_aux
use par_mod
use eigen_parameters
use communications
use calc_rhs
use file_io
use phys_ini
use aux_fields
use par_other ,only: precond_approx

implicit none
public:: get_expl_sl_mat, manipulate_matrix, output_operator
public:: mm_format
private

character(len=10):: mm_format='coordinate'

contains
  !>Computes an explicit, parallel representation of the linear gyrokinetic operator
  subroutine get_expl_sl_mat(sl_mat)
    !>full matrix containing the linear operator, parallel in the first dimension
    complex,dimension(0:vlen-1,0:vlen_gl-1),intent(inout):: sl_mat 
    complex,dimension(:,:,:,:,:,:),allocatable:: g_loc, rhs_loc
    double precision :: time_1, time_2
    integer:: col, coord, writeproc, memory
    logical:: del_old

    PERFON('build_ex')
    
    memory=int(8./1024**2*vlen_gl*vlen)
    if (prec.eq.'DOUBLE') memory=2*memory
    
    if (mype.eq.0) then
       write(*,"(A)") '--------------------------------------------------------------'
       write(*,"(A)") ''
       write(*,"(A)") 'Creating explicit matrix representation of the linear operator ' 
       write(*,"(a,i8)") 'Rank of the linear operator: ',vlen_gl
       write(*,"(A,i5,A,F8.2,A,F8.2,A)") 'Memory requirements: ',memory,' MB per core, i.e. ',&
         1./1024*memory*n_procs_sim,' GB in total'
    endif
    
    call get_systime(time_1)
    
    allocate(g_loc(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),&
         rhs_loc(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    g_loc=0.
    del_old=.false.

    call initialize_CalFullRhs(for_petsc=.true.)

    do col=0,vlen_gl-1
       call set_g(col,g_loc,coord,del_old,writeproc)
       call CalFullRhs(g_loc,rhs_loc,0)
       call ccopy(vlen,rhs_loc(li1,lj1,lk1,ll1,lm1,ln1),1,sl_mat(0,col),1)
    enddo

    call finalize_CalFullRhs

    deallocate(g_loc,rhs_loc)
    
    call get_systime(time_2)
    if (mype.eq.0) then
       write(*,"(A,F8.2)") 'Time to build explicit L: ', time_2-time_1  
       write(*,"(A)") ''
       write(*,"(A)") '--------------------------------------------------------------'
       write(*,"(A)") ''
    endif
    PERFOFF
    
  end subroutine get_expl_sl_mat

  !>Wrapper to initialize an explicit representation of the linear operator and write it to a file
  subroutine output_operator
    complex, dimension(:,:), allocatable:: sl_mat
    
    final_init=.true.
    call initialize_current_parall
    call initialize_calc_aux_fields

    allocate(sl_mat(vlen,vlen_gl))
    call get_expl_sl_mat(sl_mat)
    call write_matrix(sl_mat)
    deallocate(sl_mat)

    call finalize_calc_aux_fields
    call finalize_current_parall   

  end subroutine output_operator

  !>Writes the linear operator to a file
  !!
  !!The two matrix market formats 'array' and 'coordinate' are possible. The routine is a bit slow,
  !!but this is mostly due to the formatted output and the I/O bandwidth. Separating the computation of the
  !!individual lines for the 'coordinates' output (which then runs in parallel on the individual procs) 
  !!and the actual write statement needs much more memory and doesn't help very much.
  subroutine write_matrix(sl_mat)
    complex,dimension(vlen,vlen_gl),intent(inout):: sl_mat 
    complex,dimension(:,:,:),allocatable:: mat_out
   
    character(len=120):: header
    real:: absval
    double precision :: time_1, time_2
    integer:: MATFILE, loc_nnz, nnz, ierr, memory
    integer:: ind1, ind2, pr_ind, bl_ind

    call get_systime(time_1)

    memory=int(8./1024**2*vlen_gl*vlen)
    if (prec.eq.'DOUBLE') memory=2*memory

    if(mype.eq.0) then
       write(*,"(3A)") 'Writing matrix market ',trim(mm_format),' file of the linear operator to file (this may take a while)'
       write(*,*) 'Additional memory needed for output: ',memory,' MB'
    endif

    allocate(mat_out(vlen,vlen,n_procs_sim))
    !reorder to get more convenient version for output
    call mpi_alltoall(sl_mat, vlen*vlen, MPI_COMPLEX_TYPE, &
         mat_out, vlen*vlen, MPI_COMPLEX_TYPE, MY_MPI_COMM_WORLD,ierr) 

    do pr_ind=0,n_procs_sim-1
       !only one process opens and writes at a time
       call get_unit_nr(MATFILE)
       !process 0 replaces the file and writes the header, the others append only
       if(pr_ind.eq.mype) then
          if (mype.eq.0) then
             if(mm_format.eq.'array') then
                write(header,fmt=*) '%%MatrixMarket matrix array complex general',new_line('A'),&
                     vlen_gl, vlen_gl
             elseif(mm_format.eq.'coordinate') then
                !the third entry of the second line should contain the number of nonzero entries; 
                !since this is not known at this point the entry has to be replaced later in this routine.
                write(header,fmt=*) '%%MatrixMarket matrix coordinate complex general',new_line('A'),&
                     vlen_gl, vlen_gl, vlen_gl*vlen_gl
             endif
             open(MATFILE, file=trim(diagdir)//'/mm_'//trim(mm_format)//trim(file_extension), &
                  form='formatted', status='replace', position='rewind')
             write(MATFILE,"(A)") trim(header)
          else
             open(MATFILE, file=trim(diagdir)//'/mm_'//trim(mm_format)//trim(file_extension), &
                  form='formatted', status='old', position='append')
          endif
          if(mm_format.eq.'array') then
             !output full matrix
             do ind2=1,vlen
                do bl_ind=1,n_procs_sim
                   do ind1=1,vlen
                      write(MATFILE,fmt=*)real(mat_out(ind1,ind2,bl_ind)),aimag(mat_out(ind1,ind2,bl_ind))
                   enddo
                enddo
             enddo
          elseif(mm_format.eq.'coordinate') then
             loc_nnz=0
             do ind2=1,vlen
                do bl_ind=1,n_procs_sim
                   do ind1=1,vlen
                      absval=abs(mat_out(ind1,ind2,bl_ind))
                      if(absval.gt.epsilon(absval)) then
                         write(MATFILE,fmt=*) ind1+vlen*(bl_ind-1),ind2+vlen*mype,&
                              real(mat_out(ind1,ind2,bl_ind)),aimag(mat_out(ind1,ind2,bl_ind))
                         loc_nnz=loc_nnz+1
                      end if
                   enddo
                enddo
             enddo
          end if
          close(MATFILE)
       endif
       call my_barrier
    enddo
    deallocate(mat_out)
 
    if(mm_format.eq.'coordinate') then
       !add information about the number of nonzero entries to the file
       call mpi_reduce(loc_nnz, nnz, 1,&
            MPI_INTEGER, MPI_SUM, 0, MY_MPI_COMM_WORLD, ierr)
       if(mype.eq.0) then
          call get_unit_nr(MATFILE)
          open(MATFILE, file=trim(diagdir)//'/mm_'//trim(mm_format)//trim(file_extension), &
               status='old', form='formatted', access='direct', recl=len_trim(header),iostat=ierr)
          !rewrite header
          write(header,fmt=*) '%%MatrixMarket matrix coordinate complex general',new_line('A'),&
               vlen_gl, vlen_gl, nnz
          write(MATFILE,"(A)",rec=1) trim(header)
          close(MATFILE,iostat=ierr)
       endif
    endif

    call get_systime(time_2)
    if (mype.eq.0) then
       write(*,"(A,F8.2)") 'Time to write matrix: ', time_2-time_1  
       write(*,"(A)") ' '
    endif


  end subroutine write_matrix


  !>Auxiliary routine that sets one entry of g_1 to 1, rest zero, to compute one column of
  !!the linear operator
  subroutine set_g(pos,loc_g,coord,del_old,writeproc)
    complex,dimension(0:vlen-1),intent(inout):: loc_g !<the distribution function g
    integer,intent(in):: pos
    integer,intent(inout):: coord, writeproc
    logical,intent(inout):: del_old
    if (del_old) then
       !remove old entry
       if (mype.eq.writeproc) loc_g(coord)=0.
    endif
    
    !compute new coordinates
    writeproc=pos/vlen
    coord=modulo(pos,vlen)
    if (mype.eq.writeproc) loc_g(coord)=1.
    
    del_old=.true.
  end subroutine set_g

  !>Routine to manipulate the explicit, full, parallel representation of the linear operator L.
  !!
  !!Replaces L with (L-ev_shift) for shift_invert_s or (1-dt*L) for implicit Euler schemes and 
  !!creates the Hermitian adjoint of the matrix if left eigenvectors are to be computed. The 
  !!Hermitian adjoint is taken instead of the simple transpose in order to get the same results
  !!as the scalapack solver for the left eigenvectors. The eigenvalues for the Hermitian adjoint
  !!matrix are complex conjugate to the initial ones, which is corrected in the back transform.
  subroutine manipulate_matrix(sl_mat)
    complex, dimension(0:vlen-1,0:vlen_gl-1),intent(inout):: sl_mat !<in:L, out: manipulated matrix 
    complex,dimension(:,:,:),allocatable:: aux_mat
    integer:: ind, ierr, proc

    !add diagonal entries
    if (which_ev.eq.'shift_invert_s') then
       !for shift-invert spectral transform
       do ind=0,vlen-1
          sl_mat(ind,ind+mype*vlen)=sl_mat(ind,ind+mype*vlen)-ev_shift
       enddo
    elseif((timescheme.eq.'IE1s').or.(timescheme.eq.'IE1f')) then
       !for implicit time stepping
       sl_mat=-dt_max*sl_mat
       do ind=0,vlen-1
          sl_mat(ind,ind+mype*vlen)=1.+sl_mat(ind,ind+mype*vlen)
       enddo
    end if

    if (ev_left) then
       !this creates the Hermitian adjoint linear operator
       allocate(aux_mat(0:vlen-1,0:vlen-1,0:n_procs_sim-1))
       do proc=0,n_procs_sim-1
          aux_mat(:,:,proc)=conjg(transpose(sl_mat(:,proc*vlen:(proc+1)*vlen-1)))
       enddo
       call mpi_alltoall(aux_mat, vlen*vlen, MPI_COMPLEX_TYPE, &
            sl_mat, vlen*vlen, MPI_COMPLEX_TYPE, MY_MPI_COMM_WORLD,ierr)
       deallocate(aux_mat)
       if(mype.eq.0) write(*,"(A)") 'Transposed linear operator matrix in order to obtain left eigenvectors'
    endif

  end subroutine manipulate_matrix

  !>Routine to gather the full matrix on one processor (not used)
  subroutine sl_2_gl(sl_mat_,gl_mat_)
    complex,dimension(0:vlen-1,0:vlen_gl-1), intent(in):: sl_mat_
    complex,dimension(0:vlen_gl-1,0:vlen_gl-1), intent(out):: gl_mat_
    integer:: ierr, col
    
    do col=0,vlen_gl-1
       call mpi_gather(sl_mat_(0,col),vlen, MPI_COMPLEX_TYPE,&
            gl_mat_(0,col),vlen,MPI_COMPLEX_TYPE,0,MY_MPI_COMM_WORLD,ierr)
    enddo
    
  end subroutine sl_2_gl

#ifdef F2003_NO_NEW_LINE
  CHARACTER(1) FUNCTION new_line(ch) 
    CHARACTER(len=*) :: ch
    new_line = achar(10)
  END FUNCTION new_line
#endif

end module fullmatrix_aux
