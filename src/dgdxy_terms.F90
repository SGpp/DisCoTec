#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the terms containing x and y derivatives of chi and g_1 (i.e. the drive and curvature terms)
module dgdxy_terms
  use par_mod
  use discretization
  use discretization_adptv_module
  use blockindex
  use coordinates
  use prefactors
  use external_contr
  use x_derivatives, only: x_deriv_exc, x_deriv, x_deriv_1d_adptv, x_deriv_exc_2d_adptv

  implicit none
  public:: mem_est_dgdxy, initialize_dgdxy, add_dgdxy, finalize_dgdxy
  ! routines for adaptive grids
  public:: add_dgdxy_adptv
  private

  complex,dimension(:,:),allocatable:: iki, ikj
  complex, dimension(:,:,:), allocatable:: curv_pref

  integer:: lzero1, lzero2, init_status = 0
contains

  function mem_est_dgdxy(mem_req_in)
    real:: mem_req_in, mem_est_dgdxy
    real:: mem_loc
    
    if (xy_local) then
       if (perf_vec(1).eq.1) then
          mem_loc=mem_est_dgdxy_1()
       else 
          mem_loc=mem_est_dgdxy_2()
       end if
    else
       mem_loc=mem_est_dgdxy_df_1()
    endif

    mem_est_dgdxy=mem_req_in+mem_loc
    
  end function mem_est_dgdxy

  subroutine initialize_dgdxy

    if ((init_status.ne.perf_vec(1)).and.(init_status.gt.0)) &
         call finalize_dgdxy

    if (init_status.eq.0) then
       if (xy_local) then
          if (perf_vec(1).eq.1) then
             call initialize_dgdxy_1
          else
             call initialize_dgdxy_2
          end if
       else
          call initialize_dgdxy_df_1
       endif

       init_status = perf_vec(1)
    endif

  end subroutine initialize_dgdxy

  subroutine add_dgdxy(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    integer, intent(IN) :: lb1,lb2
    complex, dimension(lbi:ubi, lj1:lj2, 1:*), intent(in):: g_block
    !complex, dimension(lbi:ubi, lj1:lj2, 1:lklmn0), intent(in):: g_block
    real, dimension(pi1:pi2, pj1:pj2, 1:lklmn0), intent(in):: p_pdg1di,p_pdg1dj
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer:: p_dgdxy

    if (xy_local) then
       if (perf_vec(1).eq.1) then
          call add_dgdxy_1(g_block,p_rhs,p_pdg1di,p_pdg1dj,lb1,lb2)
       else
          call add_dgdxy_2(g_block(:,:,lb1:lb2),curv_pref(:,:,lb1:lb2),p_rhs,lb1,lb2)
       end if
    else
       call add_dgdxy_df_1(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    endif
    
  end subroutine add_dgdxy

  subroutine finalize_dgdxy

    if (xy_local) then
       if (init_status.eq.1) then
          call finalize_dgdxy_1
       else
          call finalize_dgdxy_2
       end if
    else
       call finalize_dgdxy_df_1
    endif

    init_status = 0

  end subroutine finalize_dgdxy

!---------------------------------------------------------------------------------

  function mem_est_dgdxy_1()
    real :: mem_est_dgdxy_1
    !iki, ikj, localrhs
    mem_est_dgdxy_1=3.*lij0*SIZE_OF_COMPLEX_MB
    
  end function mem_est_dgdxy_1

  subroutine initialize_dgdxy_1
    integer:: i,j

    allocate(iki(li1:li2,lj1:lj2),ikj(li1:li2,lj1:lj2))
    do j=lj1,lj2
       do i=li1,li2
          iki(i,j)=imag*ki(i)
          ikj(i,j)=imag*kj(j)
       enddo
    enddo
    
  end subroutine initialize_dgdxy_1

  subroutine add_dgdxy_1(p_g,p_rhs,p_pdg1di,p_pdg1dj,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, 1:lklmn0), intent(in):: p_g
    real, dimension(pi1:pi2, pj1:pj2, 1:lklmn0), intent(in):: p_pdg1di,p_pdg1dj
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs    

    complex, dimension(li1:li2, lj1:lj2) :: localrhs
    integer:: klmn

    PERFON_I('dgdxy1')

    do klmn=lb1,lb2
!       n = ln1+(klmn-1)/lklm0
!       m = lm1+Modulo(klmn-1,lklm0)/lkl0
!       l = ll1+Modulo(klmn-1,lkl0)/lk0
!       k = lk1+Modulo(klmn-1,lk0)

       ! g terms
       localrhs=p_pdg1di(pi1,pj1,klmn)*iki + p_pdg1dj(pi1,pj1,klmn)*ikj
       p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + localrhs* p_g(:,:,klmn)
       if (with_phi_ext) call add_external_contr(p_rhs(:,:,klmn),&
            &0,p_g(li1:li2,lj1:lj2,klmn), kxind_phi_ext,sk(klmn),sl(klmn),&
            &sm(klmn),sn(klmn))
       if (with_apar_ext) call add_external_contr(p_rhs(:,:,klmn),&
            &3,p_g(li1:li2,lj1:lj2,klmn), kxind_apar_ext,sk(klmn),sl(klmn),&
            &sm(klmn),sn(klmn))
    enddo

    PERFOFF_I

  end subroutine add_dgdxy_1

  subroutine finalize_dgdxy_1
    deallocate(iki,ikj)
  end subroutine finalize_dgdxy_1

!---------------------------------------------------------------------------------

  function mem_est_dgdxy_2()
    real:: mem_est_dgdxy_2
    
    mem_est_dgdxy_2=lijklmn0*SIZE_OF_COMPLEX_MB
    
  end function mem_est_dgdxy_2

  subroutine initialize_dgdxy_2
    integer:: klmn,i,j !,k,l,m,n

    allocate(curv_pref(li1:li2,lj1:lj2,1:lklmn0)) 
    do klmn=1,lklmn0
       do j=lj1,lj2
          do i=li1,li2
             curv_pref(i,j,klmn)=imag*(pdg1dj(pi1,pj1,sk(klmn),sl(klmn),&
                  &sm(klmn),sn(klmn))*kj(j)+&
                  &pdg1di(pi1,pj1,sk(klmn),sl(klmn),&
                  &sm(klmn),sn(klmn))*ki(i))
          enddo
       enddo
    enddo

  end subroutine initialize_dgdxy_2

  subroutine add_dgdxy_2(p_g,curv_pref,p_rhs,lb1,lb2)
    integer, intent(in):: lb1,lb2
    complex, dimension(lij0*lbg0), intent(in):: p_g, curv_pref
    complex, dimension(lij0*lbg0), intent(inout):: p_rhs    
    
    PERFON_I('dgdxy2')
    p_rhs=p_rhs+curv_pref*p_g

    if (with_phi_ext.or.with_apar_ext) then
       ! external contribution (radial derivative)
       call add_external_block(p_g,p_rhs,lb1,lb2)
    endif

    PERFOFF_I

  end subroutine add_dgdxy_2

  subroutine add_external_block(p_g,p_rhs,lb1,lb2)
    integer,intent(in):: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(in):: p_g
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs        
    integer:: klmn
    if (with_phi_ext) then
       do klmn=lb1,lb2
          call add_external_contr(p_rhs(li1:li2,lj1:lj2,klmn),0,&
               &p_g(li1:li2,lj1:lj2,klmn), kxind_phi_ext,sk(klmn),&
               &sl(klmn),sm(klmn),sn(klmn))
       enddo
    endif
    if (with_apar_ext) then
       do klmn=lb1,lb2
          call add_external_contr(p_rhs(li1:li2,lj1:lj2,klmn),3,&
               &p_g(li1:li2,lj1:lj2,klmn), kxind_apar_ext,sk(klmn),&
               &sl(klmn),sm(klmn),sn(klmn))
       enddo
    endif

  end subroutine add_external_block

  subroutine finalize_dgdxy_2
    deallocate(curv_pref)
  end subroutine finalize_dgdxy_2

!---------------------------------------------------------------------------------

  function mem_est_dgdxy_df_1()
    real:: mem_est_dgdxy_df_1

    !ikj, dgdx, g_func
    mem_est_dgdxy_df_1=(2.*lij0+lx0*lj0)*SIZE_OF_COMPLEX_MB
    
  end function mem_est_dgdxy_df_1

  subroutine initialize_dgdxy_df_1
    integer:: i,j

    allocate(ikj(li1:li2,lj1:lj2))
    do j=lj1,lj2
       do i=li1,li2
          ikj(i,j)=imag*kj(j)
       enddo
    enddo

  end subroutine initialize_dgdxy_df_1

  subroutine add_dgdxy_df_1(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    integer,intent(IN) :: lb1, lb2
    complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in):: g_block
    real, dimension(pi1:pi2, pj1:pj2, 1:lklmn0), intent(in):: p_pdg1di,p_pdg1dj
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs    
    complex, dimension(:,:,:,:), pointer:: p_dgdxy

    complex, dimension(li1:li2, lj1:lj2) :: dgdx
    !complex, dimension(lbi:ubi, lj1:lj2) :: g_func
    integer:: j,klmn
    logical:: save_dgdxy
    
    PERFON_I('dgdxy1')
    
    save_dgdxy=associated(p_dgdxy)

    do klmn=lb1,lb2
#ifdef prec_x
       !g_func(li1:li2,:) = p_g(:,:,klmn)
       !do j=lj1,lj2
       !   g_func(li1:li2,j) = p_g(:,j,klmn)
       !enddo
       call x_deriv(g_block(:,:,klmn),dgdx)
#else
       if(precond_approx) then
          dgdx=(0.,0.)
       else
          do j=lj1,lj2
             g_func(li1:li2,j) = p_g(:,j,klmn)
          enddo
          call x_deriv_exc(g_func,dgdx)
       end if
#endif
       if (save_dgdxy) then
          p_dgdxy(:,:,1,klmn-lb1+1) = dgdx
          p_dgdxy(:,:,2,klmn-lb1+1) = ikj*g_block(li1:li2,:,klmn)
       endif

       !take out the curvature drift term in y_global neoclassical computations
       if (.not.(x_local.and.only_neo)) then
          do j=lj1,lj2
             ! g terms
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + p_pdg1di(:,pj1,klmn)*dgdx(:,j) + &
                  &p_pdg1dj(:,pj1,klmn)*ikj(:,j)*g_block(li1:li2,j,klmn)
          enddo
       endif
    enddo

    PERFOFF_I

  end subroutine add_dgdxy_df_1

  subroutine finalize_dgdxy_df_1
    deallocate(ikj)
  end subroutine finalize_dgdxy_df_1

  subroutine add_dgdxy_adptv(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    integer, intent(IN) :: lb1,lb2
    complex, dimension(lbi:ubi, lj1:lj2, 1:*), intent(in):: g_block
    !complex, dimension(lbi:ubi, lj1:lj2, 1:lklmn0), intent(in):: g_block
    real, dimension(pi1:pi2, pj1:pj2, 1:lklmn0), intent(in):: p_pdg1di,p_pdg1dj
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer:: p_dgdxy

    if (xy_local) then
       stop "no adaptive grids for xy local version!"
    else
       call add_dgdxy_df_1_adptv(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    endif
    
  end subroutine add_dgdxy_adptv

  subroutine add_dgdxy_df_1_adptv(g_block,p_rhs,p_dgdxy,p_pdg1di,p_pdg1dj,lb1,lb2)
    integer,intent(IN) :: lb1, lb2
    complex, dimension(lbi:ubi, lj1:lj2, lb1:lb2), intent(in):: g_block
    real, dimension(pi1:pi2, pj1:pj2, 1:lklmn0), intent(in):: p_pdg1di,p_pdg1dj
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs    
    complex, dimension(:,:,:,:), pointer:: p_dgdxy

    complex, dimension(li1:li2, lj1:lj2) :: dgdx
    integer:: j, klmn, l, m, i1, i2
    logical:: save_dgdxy
    
    PERFON_I('dgdxy1')
    
    save_dgdxy=associated(p_dgdxy)

    do klmn=lb1,lb2
       l = sl(klmn)
       m = sm(klmn)
       i1 = li1_vwadp(l,m)
       i2 = li2_vwadp(l,m)
#ifdef prec_x
       call x_deriv_1d_adptv(g_block(:,:,klmn),dgdx,i1,i2)
#else
       if(precond_approx) then
          dgdx=(0.,0.)
       else
          do j=lj1,lj2
             g_func(i1:i2,j) = p_g(i1:i2,j,klmn)
          enddo
          call x_deriv_exc_2d_adptv(g_func,dgdx)
       end if
#endif
       if (save_dgdxy) then
          p_dgdxy(i1:i2,:,1,klmn-lb1+1) = dgdx(i1:i2,:)
          p_dgdxy(i1:i2,:,2,klmn-lb1+1) = ikj*g_block(i1:i2,:,klmn)
       endif
       
       !take out the curvature drift term in y_global neoclassical computations
       if (.not.(x_local.and.only_neo)) then
          do j=lj1,lj2
             ! g terms
             p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + p_pdg1di(i1:i2,pj1,klmn)*dgdx(i1:i2,j) + &
                  &p_pdg1dj(i1:i2,pj1,klmn)*ikj(i1:i2,j)*g_block(i1:i2,j,klmn)
          enddo
       endif
    enddo

    PERFOFF_I

  end subroutine add_dgdxy_df_1_adptv

end module dgdxy_terms
