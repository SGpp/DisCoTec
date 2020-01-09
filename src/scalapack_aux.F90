#include "redef.h"

module scalapack_aux
  use par_mod
  implicit none
  integer:: loc_rows, loc_cols
  integer:: blocksize=64
  integer,dimension(9):: sl_desc, bl_desc
  integer,dimension(9):: blvec_desc, slvec_desc
  integer:: sl_con, bl_con
  integer:: nprow, npcol
  complex,dimension(:,:),allocatable:: bl_mat

contains
  subroutine comp_bl_dims(bl_rows,bl_cols,nprow,npcol)
    integer,intent(out):: bl_rows,bl_cols
    integer:: nblocks, iceil, nprow, npcol
    
    !NACHDENKEN!!
    npcol=int(sqrt(dble(n_procs_sim)))
    do while(mod(n_procs_sim,npcol).ne.0)
       npcol=npcol-1
    enddo
    nprow=n_procs_sim/npcol
    
    nblocks=ceiling(real(vlen_gl)/blocksize)
    
    ! array dimensions for 2D block cyclic distributions
    bl_cols=iceil(iceil(vlen_gl,blocksize),npcol)*blocksize
    bl_rows=iceil(iceil(vlen_gl,blocksize),nprow)*blocksize

  end subroutine comp_bl_dims

  subroutine initialize_bl
    integer:: ierr

    call comp_bl_dims(loc_rows,loc_cols,nprow,npcol)
    !block cyclic distribution of the linear operator matrix
    allocate(bl_mat(loc_rows,loc_cols))
    call sl_init(bl_con,nprow,npcol)
    call descinit(bl_desc,vlen_gl,vlen_gl,blocksize,blocksize,0,0,bl_con,loc_rows,ierr)

  end subroutine initialize_bl
  
  subroutine finalize_bl
    
    deallocate(bl_mat)
    call blacs_gridexit(bl_con)
    
  end subroutine finalize_bl
  
  
  subroutine initialize_sl_desc
    integer:: ierr
    
    !generate SCALAPACK descriptor for sl_mat
    call sl_init(sl_con,n_procs_sim,1)
    call descinit(sl_desc,vlen_gl,vlen_gl,vlen,vlen_gl,0,0,sl_con,vlen,ierr)  
    
  end subroutine initialize_sl_desc
  
  subroutine finalize_sl_desc
    
    call blacs_gridexit(sl_con)
    
  end subroutine finalize_sl_desc
  
  subroutine sl_2_bl(sl_mat_,bl_mat_)
    complex,dimension(0:vlen-1,0:vlen_gl-1), intent(in):: sl_mat_
    complex,dimension(loc_rows,loc_cols), intent(out):: bl_mat_

    call pcgemr2d(vlen_gl,vlen_gl,sl_mat_,1,1, &
         sl_desc,bl_mat_,1,1,bl_desc,bl_con)

  end subroutine sl_2_bl

  subroutine sl_2_bl_vec(sl_vec_,bl_vec_)
    complex,dimension(0:vlen-1), intent(in):: sl_vec_
    complex,dimension(loc_rows,1), intent(out):: bl_vec_

    call pcgemr2d(vlen_gl,1,sl_vec_,1,1, &
         sl_desc,bl_vec_,1,1,bl_desc,bl_con)

  end subroutine sl_2_bl_vec

  subroutine bl_2_sl_vec(bl_vec_,sl_vec_)
    complex,dimension(loc_rows), intent(in):: bl_vec_
    complex,dimension(0:vlen-1,1), intent(out):: sl_vec_
   
    call pcgemr2d(vlen_gl,1,bl_vec_,1,1, &
         blvec_desc,sl_vec_,1,1,slvec_desc,sl_con)

  end subroutine bl_2_sl_vec

  subroutine bl_2_sl(bl_mat_,sl_mat_)
    complex,dimension(loc_rows,loc_cols), intent(in):: bl_mat_
    complex,dimension(0:vlen-1,0:vlen_gl-1), intent(out):: sl_mat_
   
    call pcgemr2d(vlen_gl,vlen_gl,bl_mat_,1,1, &
         bl_desc,sl_mat_,1,1,sl_desc,sl_con)

  end subroutine bl_2_sl

end module scalapack_aux
