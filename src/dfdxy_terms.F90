#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the terms containing x and y derivatives of chi and g_1 (i.e. the drive and curvature terms)
!!\todo move all perpendicular hyperdiffusion quantities here. 
!note: For the local code, this could also be included in dfdzv_1!
module dfdxy_terms
  use par_mod
  use discretization
  use discretization_adptv_module
!  use prefactors, only: f0_contr
  use numerical_damping
  use blockindex
  use boundaries, only: exchange_x, lx0_boundary
 
  implicit none
  public:: initialize_dfdxy, add_dfdxy, finalize_dfdxy, mem_est_dfdxy
  public:: add_dfdxy_adptv
  private

  integer :: init_status = 0
  integer,dimension(:),allocatable :: map_to_f
  real,dimension(:),allocatable:: deriv_sten

contains

  function mem_est_dfdxy(mem_req_in)
    real:: mem_req_in, mem_est_dfdxy
    
    mem_est_dfdxy=mem_req_in+lklmn0*SIZE_OF_INTEGER/(1024.)**2
    
  end function mem_est_dfdxy

  subroutine initialize_dfdxy
    integer :: klmn,k,l,m,n

    if (init_status==1) return

    allocate(map_to_f(1:lklmn0))
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                map_to_f(klmn)=(n-ln1)*lz0*lv0*lw0+(nwb+m-lm1)*lz0*lv0+&
                     &(nvb+l-ll1)*lz0+nzb+k-lk1+1
            enddo
          enddo
       enddo
    enddo

    if (Erad.ne.0..and..not.y_local) then
       allocate(deriv_sten(-2:2))
       deriv_sten=1/(12.*deli)*(/1,-8,0,8,-1/)
    endif

    init_status = 1

  end subroutine initialize_dfdxy


  !>Wrapper for dfdx & dfdy terms
  subroutine add_dfdxy(p_f,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs    
    
    if (xy_local) then
       call add_dfdxy_ff(p_f,p_rhs,lb1,lb2)
    else
       call add_dfdxy_df(p_f,p_rhs,lb1,lb2)
    endif

  end subroutine add_dfdxy


  !>Routine for the local code
  subroutine add_dfdxy_ff(p_f,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs    
    
    integer:: ind1,ind2,j,klmn

    PERFON('dfdxy')
    do klmn=lb1,lb2       
!       If (p_has_00_mode) Then
!           p_rhs(li1,0,klmn) = p_rhs(li1,0,klmn) + f0_contr(pi1,pj1,sk,sl,sm,sn)
!        Endif
       if (damp_i.and..not.GyroLES) then
          ind1=map_to_f(klmn)
          do j=lj1,lj2
             p_rhs(:,j,klmn)= p_rhs(:,j,klmn)+ hyp_i*didiff*p_f(:,j,ind1) 
          end do
       endif
       if (damp_j.and..not.GyroLES) then
          do j=lj1,lj2
             ! numerical damping in binormal direction
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) &
                  & + hyp_j * p_f(:,j,map_to_f(klmn)) * djdiff(j)
          enddo
       endif
       if (damp_perp.and..not.GyroLES) then
          ind1=map_to_f(klmn)
          ind2=Modulo(klmn-1,lk0)+lk1
          do j=lj1,lj2
             ! numerical damping in perpendicular plane
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) &
                  & - hyp_perp * p_f(:,j,ind1) * dperpdiff(:,j,ind2)
          enddo
       endif
       if (GyroLES) then
          do j=lj1,lj2
             p_rhs(:,j,klmn)= p_rhs(:,j,klmn)+ hyp_x_spec(sn(klmn))*didiff*p_f(:,j,map_to_f(klmn))& 
             + hyp_y_spec(sn(klmn)) * p_f(:,j,map_to_f(klmn)) * djdiff(j)
          end do
       endif
    enddo
    PERFOFF

  end subroutine add_dfdxy_ff


  !>Routine for radially nonlocal code
  subroutine add_dfdxy_df(p_f,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer:: i,j,sten,klmn
    complex, dimension(lbi:ubi, lj1:lj2) :: f_wb

    PERFON('dfdxy')
    do klmn=lb1,lb2
       if (damp_i.or.Erad.ne.0.) then
          do j=lj1,lj2
             f_wb(li1:li2,j) = p_f(li1:li2,j,map_to_f(klmn))
          enddo
          CALL exchange_x(lx0_boundary,f_wb)
          If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
             Do i=1,nib
                f_wb(li1-i,lj1) = f_wb(li1+i,lj1)
             Enddo
          Endif
       endif
#ifdef prec_x
       if (damp_i) then
#else
       if (damp_i.and.(.not.precond_approx)) then
#endif
          Do j=lj1,lj2
             Do i=li1,li2
                ! numerical damping of high k modes in first coordinate
                ! by adding hyp_i_order'th derivative of f1 (not g1)
                Do sten=-hyp_i_order/2,hyp_i_order/2
                   p_rhs(i,j,klmn) = p_rhs(i,j,klmn) &
                        & + hyp_i * f_wb(i+sten,j) * didiff(sten)
                End Do
             End Do
          End Do
       endif
       if (damp_j) then
          do j=lj1,lj2
             ! numerical damping in binormal direction
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) &
                  & + hyp_j * p_f(:,j,map_to_f(klmn)) * djdiff(j)
          enddo
       endif
       if (Erad.ne.0..and..not.y_local) then
             ! constant radial electric field
          do j=lj1,lj2
             do sten=-2,2
                p_rhs(:,j,klmn) = p_rhs(:,j,klmn) &
                     & + Erad * f_wb(li1+sten:li2+sten,j) * deriv_sten(sten)
             enddo
          enddo
       endif
    enddo

    PERFOFF

  end subroutine add_dfdxy_df


  subroutine finalize_dfdxy
    deallocate(map_to_f)
    init_status = 0
    if (Erad.ne.0..and..not.y_local) deallocate(deriv_sten)


  end subroutine finalize_dfdxy

  ! ---- routines modified (mirrored) for the adaptive grids -----

  subroutine add_dfdxy_adptv(p_f,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs    
    
    if (xy_local) then
       stop "no adaptive grids for xy local version!"
    else
       call add_dfdxy_df_adptv(p_f,p_rhs,lb1,lb2)
    endif

  end subroutine add_dfdxy_adptv

  subroutine add_dfdxy_df_adptv(p_f,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,1:lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer:: i, j, sten, klmn, l, m, i1, i2
    complex, dimension(lbi:ubi, lj1:lj2) :: f_wb

    PERFON('dfdxy')
    do klmn=lb1,lb2
       l = sl(klmn)
       m = sm(klmn)
       i1 = li1_vwadp(l,m)
       i2 = li2_vwadp(l,m)
       if (damp_i.or.Erad.ne.0.) then
          do j=lj1,lj2
             f_wb(i1:i2,j) = p_f(i1:i2,j,map_to_f(klmn))
          enddo
          CALL exchange_x(lx0_boundary,f_wb)
          If ((rad_bc_type.eq.2).and.(p_has_0_mode).and.(my_pex.eq.0)) then
             Do i=1,nib
                f_wb(i1-i,lj1) = f_wb(i1+i,lj1)
             Enddo
          Endif
       endif
#ifdef prec_x
       if (damp_i) then
#else
       if (damp_i.and.(.not.precond_approx)) then
#endif
          Do j=lj1,lj2
             Do i=i1,i2
                ! numerical damping of high k modes in first coordinate
                ! by adding hyp_i_order'th derivative of f1 (not g1)
                Do sten=-hyp_i_order/2,hyp_i_order/2
                   p_rhs(i,j,klmn) = p_rhs(i,j,klmn) &
                        & + hyp_i * f_wb(i+sten,j) * didiff(sten)
                End Do
             End Do
          End Do
       endif
       if (damp_j) then
          do j=lj1,lj2
             ! numerical damping in binormal direction
             p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) &
                  & + hyp_j * p_f(i1:i2,j,map_to_f(klmn)) * djdiff(j)
          enddo
       endif
       if (Erad.ne.0..and..not.y_local) then
             ! constant radial electric field
          do j=lj1,lj2
             do sten=-2,2
                p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) &
                     & + Erad * f_wb(i1+sten:i2+sten,j) * deriv_sten(sten)
             enddo
          enddo
       endif
    enddo

    PERFOFF

  end subroutine add_dfdxy_df_adptv

end module dfdxy_terms
