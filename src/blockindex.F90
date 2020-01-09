#include "redef.h"
#include "intrinsic_sizes.h"
!module containing arrays/functions to transform,
!e.g., (k,l,m,n) indices to more global indices (klmn)
!as it is required for strip mining
module blockindex
  use discretization
  implicit none

  public :: mem_est_klmn_conv_arrs, initialize_klmn_conv_arrs,&
       &finalize_klmn_conv_arrs
  public :: sk,sl,sm,sn,lbg0
  public:: blk,bll,blm,bln,dblk
  public:: blk1,blk2,bll1,bll2,blm1,blm2,bln1,bln2
  public :: compute_gy_av, initialize_blockindex

  private

  integer :: init_status = 0

  !arrays for conversion from global index klmn to k,l,m,n
  integer, dimension(:),allocatable :: sk, sl, sm, sn

  integer,dimension(:,:), allocatable:: blk, bll, blm, bln
  integer:: dblk, blk1,blk2,bll1,bll2,blm1,blm2,bln1,bln2
  !array to indicate when gyroaverages have to be performed
  !in klmn loops (see dchidxy, dchidz)
  logical, dimension(:),allocatable :: compute_gy_av

  !integer(C_INT),bind(c) :: lbg0
  integer :: lbg0

contains 

  function mem_est_klmn_conv_arrs(mem_req_in)
    real:: mem_req_in, mem_est_klmn_conv_arrs, mem_loc
    
    !sk,sl,sm,sn arrays
    mem_loc = 4.*lklmn0*SIZE_OF_INTEGER/(1024.)**2
    
    !compute_gy_av
    mem_loc = mem_loc + lklmn0*SIZE_OF_LOGICAL/(1024.)**2

    mem_est_klmn_conv_arrs=mem_req_in+mem_loc

  end function mem_est_klmn_conv_arrs

  subroutine initialize_blockindex
    lbg0 = lklmn0/nblocks
  end subroutine initialize_blockindex
    
  
  !>Initialize 4 integer arrays which allow for conversion from
  !!the global klmn index to k,l,m,n indices
  subroutine initialize_klmn_conv_arrs
    integer:: k,l,m,n,klmn,bl

    if (init_status==1) return

    allocate(sk(lklmn0),sl(lklmn0),sm(lklmn0),sn(lklmn0))
    allocate(compute_gy_av(lklmn0))

    compute_gy_av = .false.

    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                sk(klmn)=k
                sl(klmn)=l
                sm(klmn)=m
                sn(klmn)=n
                if(modulo(klmn-1,lkl0).eq.0) compute_gy_av(klmn)=.true.
             enddo
          enddo
       enddo
    enddo

    allocate(blk(2,nblocks),bll(2,nblocks),blm(2,nblocks),bln(2,nblocks))
    do bl=1,nblocks
       !first index of block
       klmn=lbg0*(bl-1)+1
       blk(1,bl)=sk(klmn)
       bll(1,bl)=sl(klmn)
       blm(1,bl)=sm(klmn)
       bln(1,bl)=sn(klmn)
       !last index of block
       klmn=lbg0*bl
       blk(2,bl)=sk(klmn)
       bll(2,bl)=sl(klmn)
       blm(2,bl)=sm(klmn)
       bln(2,bl)=sn(klmn)
    end do
    dblk=blk(2,1)-blk(1,1)+1

    init_status = 1

  end subroutine initialize_klmn_conv_arrs

  subroutine finalize_klmn_conv_arrs
    deallocate(sk,sl,sm,sn)
    deallocate(blk,bll,blm,bln)
    deallocate(compute_gy_av)
    init_status = 0
  end subroutine finalize_klmn_conv_arrs

end module blockindex

