#include "redef.h"
module impl_scalapack
  use par_mod
  use communications
  use scalapack_aux
  use fullmatrix_aux
  use calc_rhs

  implicit none
  public:: initialize_CalFullRhs_LU, CalFullRhs_LU, finalize_CalFullRhs_LU, invert_matrix
  private

  integer,dimension(:),allocatable:: pivot
  complex,dimension(:,:), allocatable:: g_bl

contains

  subroutine initialize_CalFullRhs_LU
    complex,dimension(:,:),allocatable:: lin_op_slices 
    integer:: info

    allocate(lin_op_slices(0:vlen-1,0:vlen_gl-1))

    call get_expl_sl_mat(lin_op_slices)
    
    call manipulate_matrix(lin_op_slices)
    
    call initialize_sl_desc
    call initialize_bl

    call comp_LU(lin_op_slices,bl_mat)

    deallocate(lin_op_slices)

    !initialize the descriptors for the state vector
    call descinit(slvec_desc,vlen_gl,1,vlen,1,0,0,sl_con,vlen,info)  
    call descinit(blvec_desc,vlen_gl,1,blocksize,1,0,0,bl_con,loc_rows,info)  
    allocate(g_bl(loc_rows,1))

    call initialize_CalFullRhs

  end subroutine initialize_CalFullRhs_LU

  subroutine CalFullRhs_LU(loc_g,rhs)
    complex,dimension(0:vlen-1),intent(in):: loc_g
    complex,dimension(0:vlen-1),intent(out):: rhs
    integer::info

    call sl_2_bl_vec(loc_g,g_bl) 
    call pcgetrs('N',vlen_gl,1,bl_mat,1,1,bl_desc,pivot,g_bl,1,1,blvec_desc,info)
    call bl_2_sl_vec(g_bl,rhs)

  end subroutine CalFullRhs_LU

  subroutine finalize_CalFullRhs_LU
    
    call finalize_sl_desc
    deallocate(g_bl)
    call finalize_bl
    call finalize_CalFullRhs

  end subroutine finalize_CalFullRhs_LU


  subroutine invert_matrix(sl_mat)
    complex, dimension(0:vlen-1,0:vlen_gl-1),intent(inout):: sl_mat

    call initialize_sl_desc
    call initialize_bl

    call comp_LU(sl_mat,bl_mat)
    call comp_inverse(bl_mat,sl_mat)

    call finalize_sl_desc
    call finalize_bl

  end subroutine invert_matrix 

  subroutine comp_LU(sl_mat,bl_mat)
    complex,dimension(0:vlen-1,0:vlen_gl-1),intent(in):: sl_mat
    complex,dimension(loc_rows,loc_cols),intent(out):: bl_mat

    integer:: info
    double precision :: time1,time2
    
    call get_systime(time1)
    
    call sl_2_bl(sl_mat,bl_mat)
    
    !-------------------
    !LU decomposition  
    
    PERFON('calc_LU')
    !too much but I can't make sense of the manual
    allocate(pivot(vlen_gl))
    
    call pcgetrf(vlen_gl,vlen_gl,bl_mat,1,1,bl_desc,pivot,info)
    call get_systime(time2)
    if(mype.eq.0) print*, 'Time for LU factorization: ', time2-time1
  
  end subroutine comp_LU
  
  
  !-------------------
  !explicit computation of the inverse matrix 
  subroutine comp_inverse(bl_mat,sl_mat)
    complex, dimension(0:vlen-1,0:vlen_gl-1),intent(out):: sl_mat
    complex,dimension(loc_rows,loc_cols),intent(inout):: bl_mat
    complex,dimension(:),allocatable:: work 
    integer,dimension(:),allocatable:: iwork
    integer:: lwork, liwork, info
    double precision :: time1,time2

    call get_systime(time1) 
    !compute and adapt size of work buffer
    allocate(work(1),iwork(1))
    lwork = -1
    liwork= -1
    call pcgetri(vlen_gl,bl_mat,1,1,bl_desc,pivot,work,lwork,iwork,liwork,info)
    lwork=int(real(work(1)))
    liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),iwork(liwork))
    !actual computation

    call pcgetri(vlen_gl,bl_mat,1,1,bl_desc,pivot,work,lwork,iwork,liwork,info)
    deallocate(work,iwork,pivot)

    call get_systime(time2)  
    if(mype.eq.0) print*, 'Time for computation of the inverse matrix: ', time2-time1
    
    call bl_2_sl(bl_mat,sl_mat)
    
    PERFOFF

  end subroutine comp_inverse
    
end module impl_scalapack
