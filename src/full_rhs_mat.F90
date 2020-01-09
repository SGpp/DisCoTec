#include "redef.h"

!>Computes the RHS for linear problems with an explicit representation of the gyrokinetic
!!operator
!!
!!This module is useful for investigations where an explicit representation of L is analyzed 
!!and manipulated, however full matrix computations are very slow compared to the other methods
!!implemented in GENE
module full_rhs_mat
  use par_mod
  use communications
  use fullmatrix_aux
  use calc_rhs
  use RK_coefficients, only: implicit_scheme
#ifdef WITHSCAL
  use impl_scalapack
#endif
  implicit none
  public:: initialize_CalFullRhs_fullmat, CalFullRhs_fullmat, finalize_CalFullRhs_fullmat 
  private

  complex,dimension(:,:,:),allocatable:: lin_op_mat

contains
  !>Initializes the full matrix used for the RHS computation, including the explicit matrix 
  !!inversion (very slow!) if the implicit IE1f timescheme is chosen
  subroutine initialize_CalFullRhs_fullmat
    complex,dimension(:,:),allocatable:: lin_op_slices 
    integer:: ierr
    
    allocate(lin_op_slices(0:vlen-1,0:vlen_gl-1))
    call get_expl_sl_mat(lin_op_slices)
    call initialize_CalFullRhs
    
    call manipulate_matrix(lin_op_slices)

    allocate(lin_op_mat(0:vlen-1,0:vlen-1,0:n_procs_sim-1))
#ifdef WITHSCAL
    if(implicit_scheme) then
       call invert_matrix(lin_op_slices)
    endif
#endif

    !reorder to get more convenient version of the linear operator
    call mpi_alltoall(lin_op_slices, vlen*vlen, MPI_COMPLEX_TYPE, &
         lin_op_mat, vlen*vlen, MPI_COMPLEX_TYPE, MY_MPI_COMM_WORLD,ierr) 
    
    deallocate(lin_op_slices)
  
  end subroutine initialize_CalFullRhs_fullmat

  !>Computes the RHS from the full matrix
  subroutine CalFullRhs_fullmat(loc_g,rhs)
    complex,dimension(0:vlen-1),intent(in):: loc_g !<g_1 type array
    complex,dimension(0:vlen-1),intent(out):: rhs !<result of L*g_1
    complex,dimension(0:vlen-1):: loc_rhs
    integer:: proc, ierr
    
    PERFON('rhs_ex')
    do proc=0,n_procs_sim-1
       loc_rhs=matmul(lin_op_mat(:,:,proc),loc_g)
       call mpi_reduce(loc_rhs, rhs, vlen,&
            MPI_COMPLEX_TYPE, MPI_SUM, proc, MY_MPI_COMM_WORLD, ierr)
    enddo

    PERFOFF
    
  end subroutine CalFullRhs_fullmat

  !>Deletes the explicit representation of the linear operator
  subroutine finalize_CalFullRhs_fullmat
    deallocate(lin_op_mat)
    call finalize_CalFullRhs
  end subroutine finalize_CalFullRhs_fullmat

end module full_rhs_mat
