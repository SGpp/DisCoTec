#include "redef.h"
#include "intrinsic_sizes.h"

!>Computes the f0_contr contribution
module f0_term
  use par_mod
  use discretization
  use discretization_adptv_module
  use prefactors, only: edr, curv, B_Bstar, with_coriolis, with_centrifugal, &
       with_bxphi0, coriolis_i, centrifugal_i, bxphi0_i
  use blockindex
  use geometry
  use sources_mod, only: initialize_f0_heat_source, &
       &finalize_f0_heat_source, puls_func,f0_heat_src_active,&
       &puls_amp, puls_start
  !use boundaries, only: exchange_x, lx0_boundary
  use vel_space 

  implicit none
  public:: initialize_f0_term, add_f0_term, finalize_f0_term
  public:: mem_est_f0_term, add_f0_term_df_1d
  ! routines modified for adaptive grids
  public:: add_f0_term_adptv

  public:: f0_contr

  private

  Real, Dimension(:,:,:), Allocatable :: f0_contr

  integer :: init_status = 0

contains

  function mem_est_f0_term(mem_req_in)
    real:: mem_req_in, mem_est_f0_term, mem_loc
    
    !f0 term
    mem_loc=SIZE_OF_REAL_MB*pi0*pj0*lklmn0

    !localized heat_source
    !mem_loc = mem_est_localized_heat_src

    mem_est_f0_term=mem_req_in+SIZE_OF_REAL_MB*pi0*pj0*lklmn0
    
  end function mem_est_f0_term

  subroutine initialize_f0_term
    integer :: j,klmn,pni 
    real, dimension(pi1:pi2):: Erad_contr 
    if (init_status==1) return

    Allocate(f0_contr(pi1:pi2, pj1:pj2, 1:lklmn0))
    f0_contr = 0.0

    if (.not.x_local) then
       call initialize_f0_heat_source(f0_contr)
    endif

    pni =pn1
    if (include_f0_contr) then
       !inhomogeneity stemming from pure F0 terms
       if (yx_order) then
          do klmn=1,lklmn0
             if (pn0.gt.1) pni=sn(klmn)
             do j=pj1,pj2
                if(Erad.ne.0.) then 
                   Erad_contr(:) = B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) /  C_xy(pi1:pi2) * spec(sn(klmn))%charge * Erad /&
                     &spec(sn(klmn))%temp_prof(pi1:pi2) * fm(pi1:pi2,j,sk(klmn),sl(klmn),sm(klmn),pni)
                   f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                        &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                        &( edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn))  + Erad_contr)* &
                        &curv(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                        &geom%K_j(:,j,sk(klmn))
                else
                   f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                        &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                        &( edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)))* &
                        &curv(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                        &geom%K_j(:,j,sk(klmn))
                endif
             enddo
          enddo
       else
          do klmn=1,lklmn0
             if (pn0.gt.1) pni=sn(klmn)
             do j=pj1,pj2
                f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                     &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                     &edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                     &curv(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                     &geom%K_i(:,j,sk(klmn))
             enddo
          enddo
          if (with_coriolis) then
             do klmn=1,lklmn0
                if (pn0.gt.1) pni=sn(klmn)
                do j=pj1,pj2
                   f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                        &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                        &edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                        &coriolis_i(j,sk(klmn),sl(klmn),sn(klmn))
                enddo
             enddo
          endif
          if (with_centrifugal) then
             do klmn=1,lklmn0
                if (pn0.gt.1) pni=sn(klmn)
                do j=pj1,pj2
                   f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                        &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                        &edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                        &centrifugal_i(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
                enddo
             enddo
          endif
          if (with_bxphi0) then
             do klmn=1,lklmn0
                if (pn0.gt.1) pni=sn(klmn)
                do j=pj1,pj2
                   f0_contr(:,j,klmn)= f0_contr(:,j,klmn) + C_xy(pi1:pi2)/&
                        &B_Bstar(j,sk(klmn),sl(klmn),sn(klmn)) * &
                        &edr(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                        &bxphi0_i(j,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
                enddo
             enddo
          endif
       end if
    end if

    init_status = 1

  end subroutine initialize_f0_term


  !>Wrapper for f0_term
  subroutine add_f0_term(p_rhs,lb1,lb2,time)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs  
    real, intent(in) :: time

    if (xy_local) then
       call add_f0_term_ff(p_rhs,lb1,lb2)
    else
       call add_f0_term_df(p_rhs,lb1,lb2,time)
    endif
 
  end subroutine add_f0_term


  !>Routine adding f0_contribution for the local code
  subroutine add_f0_term_ff(p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs    
    
    integer:: klmn

    !PERFON('f0_contr')
    If (p_has_00_mode) Then
       do klmn=lb1,lb2       
          p_rhs(li1,lj1,klmn) = p_rhs(li1,lj1,klmn) + &
               &f0_contr(pi1,pj1,klmn)
       enddo
    endif
    !PERFOFF

  end subroutine add_f0_term_ff


  !>Routine for radially nonlocal code
  !!This includes the option of pulsed f0 heat source contributions
  !!(which should naturally not be used in the present form
  !! if the neoclassical term is switched on)
  subroutine add_f0_term_df(p_rhs,lb1,lb2,time)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    real, intent(in) :: time
    
    integer :: klmn

    !PERFON('f0')
    if (p_has_0_mode) then
       if ((puls_amp.eq.1.0).or.(time.lt.puls_start)) then
          do klmn=lb1,lb2
             p_rhs(:,lj1,klmn) = p_rhs(:,lj1,klmn) + &
                  &f0_contr(pi1:pi2,pj1,klmn)          
          enddo
       else
          do klmn=lb1,lb2
             p_rhs(:,lj1,klmn) = p_rhs(:,lj1,klmn) + &
                  &puls_func(time)*f0_contr(pi1:pi2,pj1,klmn)          
          enddo          
       endif
    endif
    !PERFOFF

  end subroutine add_f0_term_df

  !>Same routine as add_f0_term_df but different interface
  !!(p_rhs shape + indices) for diagnostic reasons
  subroutine add_f0_term_df_1d(time,k,l,m,n,p_rhs)
    integer,intent(in) :: k,l,m,n
    complex, dimension(li1:li2), intent(inout):: p_rhs
    real, intent(in) :: time

    integer :: klmn

    !PERFON('f0')
    if (p_has_0_mode) then
       klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
       if ((puls_amp.eq.1.0).or.(time.lt.puls_start)) then
          p_rhs = p_rhs + f0_contr(pi1:pi2,pj1,klmn) 
       else
          p_rhs(:) = p_rhs(:) + puls_func(time)*&
               &f0_contr(pi1:pi2,pj1,klmn)          
       endif
    endif
    !PERFOFF

  end subroutine add_f0_term_df_1d

  subroutine finalize_f0_term

    deallocate(f0_contr)

    if (.not.x_local) call finalize_f0_heat_source

    init_status = 0
  end subroutine finalize_f0_term

  ! ---- modified routines for adaptive grids -----

  subroutine add_f0_term_adptv(p_rhs,lb1,lb2,time)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs  
    real, intent(in) :: time

    if (xy_local) then
       stop "no adaptive grids for local version"
    else
       call add_f0_term_df_adptv(p_rhs,lb1,lb2,time)
    endif
 
  end subroutine add_f0_term_adptv

  subroutine add_f0_term_df_adptv(p_rhs,lb1,lb2,time)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    real, intent(in) :: time
    
    integer :: klmn, l, m, i1, i2, pi1_a, pi2_a

    !PERFON('f0')
    if (p_has_0_mode) then
       if ((puls_amp.eq.1.0).or.(time.lt.puls_start)) then
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             pi1_a = pi1_vwadp(l,m)
             pi2_a = pi2_vwadp(l,m)
             p_rhs(i1:i2,lj1,klmn) = p_rhs(i1:i2,lj1,klmn) + &
                  &f0_contr(pi1_a:pi2_a,pj1,klmn)          
          enddo
       else
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             pi1_a = pi1_vwadp(l,m)
             pi2_a = pi2_vwadp(l,m)
             p_rhs(i1:i2,lj1,klmn) = p_rhs(i1:i2,lj1,klmn) + &
                  &puls_func(time)*f0_contr(pi1_a:pi2_a,pj1,klmn)          
          enddo          
       endif
    endif
    !PERFOFF

  end subroutine add_f0_term_df_adptv

end module f0_term
