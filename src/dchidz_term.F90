#include "redef.h"
#include "intrinsic_sizes.h"

!>Computes the term containing the parallel derivative of chi
module dchidz_term
  use par_mod
  use discretization
  use discretization_adptv_module
  use blockindex
  use prefactors
!  use boundaries, only: exchange_z
  use gyro_average_ff_mod, only: jfac, I1_factor
  use gyro_average_df_mod, only: gyro_average_df_wb
  use axpy
  use geometry, only: geom, C_xy
  use vel_space, only: fm

  implicit none
  public:: initialize_dchidz, add_dchidz, finalize_dchidz,mem_est_dchidz, z_deriv
  ! --- routines for adaptive grids ---------
  public::  add_dchidz_adptv
  private

  integer :: init_status = 0
  complex, dimension(:,:,:), allocatable:: p_dxi_dz
  real, dimension(:,:,:), allocatable:: pdphidz

#ifdef WITHOMP_BLOCKLOOP
  !$OMP THREADPRIVATE(p_dxi_dz)
#endif
contains

  function mem_est_dchidz(mem_req_in)
    real:: mem_req_in, mem_est_dchidz
    
    mem_est_dchidz = mem_req_in + (lijk0+lij0*lz0)*SIZE_OF_COMPLEX/(1024.)**2

    !pdphidz
    mem_est_dchidz = mem_est_dchidz + SIZE_OF_REAL_MB*pi0*pj0*lklmn0

  end function mem_est_dchidz
    

  subroutine initialize_dchidz
    integer :: j,klmn,pni,my_thread
#ifdef WITHOMP
    integer :: omp_get_thread_num
#endif

    if (init_status==1) return

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) private(my_thread) shared(li1,li2,lj1,lj2,lk1,lk2)
    my_thread = omp_get_thread_num()
    !$OMP CRITICAL
#else
    my_thread = 0
#endif
    !print*,mype,": Allocating p_dxi_dz on thread ",my_thread
    allocate(p_dxi_dz(li1:li2, lj1:lj2, lk1:lk2))
    !print*,mype,": Allocated  p_dxi_dz on thread ",my_thread
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif

    allocate(pdphidz(pi1:pi2, pj1:pj2, 1:lklmn0))

    ! d/dz phi prefactor
    pni=pn1
    do klmn=1,lklmn0
       if (pn0.gt.1) pni=sn(klmn)
       do j=pj1,pj2
          pdphidz(:,j,klmn) = - C_xy(pi1:pi2)*vTvpar(sl(klmn),sn(klmn))/&
               &(geom%jacobian(pi1:pi2,j,sk(klmn))*geom%Bfield&
               &(pi1:pi2,j,sk(klmn)))*spec(sn(klmn))%charge/&
               &(spec(sn(klmn))%temp*spec(sn(klmn))%temp_prof(pi1:pi2))*&
               &fm(pi1:pi2,j,sk(klmn),sl(klmn),sm(klmn),pni)
       enddo
    enddo

    init_status = 1

  end subroutine initialize_dchidz


  subroutine add_dchidz(p_emfields,p_bar_emfields,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs

    if (par_sten_bound.eq.0) return

    if (xy_local) then
       call add_dchidz_ff(p_emfields,p_rhs,lb1,lb2)
    else
       call add_dchidz_df(p_bar_emfields,p_rhs,lb1,lb2)
    endif
  end subroutine add_dchidz
  

  subroutine add_dchidz_ff(p_emfields,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs

    complex, dimension(li1:li2, lj1:lj2, lbz:ubz):: gy_xi
    integer:: klmn

    PERFON_I('dchidz')

    do klmn=lb1,lb2
       !recompute p_dxi_dz only if the m or n index has changed
#ifdef WITHOMP_BLOCKLOOP
       if ((compute_gy_av(klmn)).or.(klmn.eq.lb1)) then
#else
       if (compute_gy_av(klmn)) then
#endif
          !gyroaverage
          gy_xi = p_emfields(:,:,:,1) * jfac(:,:,:,sm(klmn),sn(klmn))
          IF (n_fields .GT. 2) gy_xi = gy_xi + mu_tjqj(sm(klmn),sn(klmn)) * &
            p_emfields(:,:,:,3) * I1_factor(:,:,:,sm(klmn),sn(klmn))

          !z derivative
          call z_deriv(gy_xi,p_dxi_dz,sl(klmn))
       endif
       
       !compute the term          
       call axpy_ij(lij0,pdphidz(pi1,pj1,klmn),&
            &p_dxi_dz(:,:,sk(klmn)),p_rhs(:,:,klmn))

!       p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+p_dxi_dz(:,:,sk(klmn))*&
!            &pdphidz(pi1,pj1,klmn)

    enddo 

    PERFOFF_I

  end subroutine add_dchidz_ff

  subroutine add_dchidz_df(p_bar_emfields,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs

!    complex, dimension(li1:li2, lj1:lj2, lbz:ubz):: gy_xi
    integer:: j,klmn

    PERFON_I('dchidz')

    do klmn=lb1,lb2
       !recompute p_dxi_dz only if the m or n index has changed
       if (compute_gy_av(klmn) &
#ifdef WITHOMP_BLOCKLOOP
            &.or.(klmn.eq.lb1) &
#endif
            &) then
          !gyroaverage
 !         call gyro_average_df_wb(p_emfields(:,:,:,1),gy_xi,sm(klmn),sn(klmn))

          !B_parallel is not supported in the global version

          !z derivative
!          call exchange_z(p_bar_emfields(:,:,:,sm(klmn),sn(klmn),1))

          call z_deriv(p_bar_emfields(:,:,:,sm(klmn),sn(klmn),1),p_dxi_dz,sl(klmn))
       endif

       !compute the term  
       do j=lj1,lj2
          p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+p_dxi_dz(:,j,sk(klmn))*&
               &pdphidz(:,pj1,klmn)
       enddo
    enddo 

    PERFOFF_I

  end subroutine add_dchidz_df
  
  subroutine z_deriv(y,dydz,l)
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz), intent(in)  :: y
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2), intent(out) :: dydz
    real, dimension(-par_sten_bound:par_sten_bound) :: par_sten2
    integer :: k, l, sten, stensign

    !z derivative
    sten = -par_sten_bound
    if ((.not.x_local).and.shifted_metric) then
       dydz = par_sten(sten)*geom%phasefac(pi1:pi2,:,:,sten)*&
            &y(:,:,lk1+sten:lk2+sten)
       
       if (parscheme(1:1).eq.'c') then
          do sten=-par_sten_bound+1,-1
             dydz = dydz + par_sten(sten)*&
                  &geom%phasefac(pi1:pi2,:,:,sten)* &
                  &y(:,:,lk1+sten:lk2+sten)
          enddo
          !skip the center value (sten = 0) for centered schemes
          do sten=1,par_sten_bound
             dydz = dydz + par_sten(sten)*&
                  &geom%phasefac(pi1:pi2,:,:,sten)* &
                  &y(:,:,lk1+sten:lk2+sten)
          enddo
       else
          do sten=-par_sten_bound+1,par_sten_bound
             dydz = dydz + par_sten(sten)*&
                  &geom%phasefac(pi1:pi2,:,:,sten)* &
                  &y(:,:,lk1+sten:lk2+sten)
          enddo
       endif
    else
       if (parscheme(1:1).eq.'c') then
          dydz = par_sten(sten)*y(:,:,lk1+sten:lk2+sten)
          do k=lk1,lk2
             do sten=-par_sten_bound+1,-1
                call axpy_ij(lij0,par_sten(sten),y(:,:,k+sten),dydz(:,:,k))
             enddo
             !skip the center value (sten = 0) for centered schemes
             do sten=1,par_sten_bound
                call axpy_ij(lij0,par_sten(sten),y(:,:,k+sten),dydz(:,:,k))
             enddo
          enddo
       else 
#if 1
          !keep dphi/dz c4th if df/dz is treated upwind 3rd
          !needs to be adapted if higher order upwind schemes shall be tested
          par_sten2 = (/1,-8,0,8,-1/)/(12.0*dz)
          dydz = par_sten2(sten)*y(:,:,lk1+sten:lk2+sten)
          do k=lk1,lk2
             do sten=-par_sten_bound+1,par_sten_bound
                call axpy_ij(lij0,par_sten2(sten),y(:,:,k+sten),dydz(:,:,k))
             enddo
          enddo
#else 
          ! Upwind in parallel direction considering sign of parallel velocity
          stensign = sign(1.0,vp(l))
          dydz = stensign*par_sten(sten)*y(:,:,lk1+stensign*sten:lk2+stensign*sten)
          do k=lk1,lk2
             Do sten=-par_sten_bound+1,par_sten_bound
                if (par_sten(sten).eq.0.0) cycle
                call axpy_ij(lij0,stensign*par_sten(sten),y(:,:,k+stensign*sten),dydz(:,:,k))
             enddo
          End Do
#endif
       endif
    endif
    
  end subroutine z_deriv
  

  subroutine finalize_dchidz
#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL
    !$OMP CRITICAL
#endif
    deallocate(p_dxi_dz)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif
    deallocate(pdphidz)
    init_status = 0
  end subroutine finalize_dchidz

  subroutine add_dchidz_adptv(p_emfields,p_bar_emfields,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs

    if (par_sten_bound.eq.0) return

    if (xy_local) then
       stop "no adaptive grids for xy local version!"
    else
       call add_dchidz_df_adptv(p_bar_emfields,p_rhs,lb1,lb2)
    endif
  end subroutine add_dchidz_adptv

  subroutine add_dchidz_df_adptv(p_bar_emfields,p_rhs,lb1,lb2)
    integer,intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs

    integer:: j, klmn, l, m, i1, i2

    PERFON_I('dchidz')

    do klmn=lb1,lb2
       l = sl(klmn)
       m = sm(klmn)
       i1 = li1_vwadp(l,m)
       i2 = li2_vwadp(l,m)
       if (compute_gy_av(klmn) &
#ifdef WITHOMP_BLOCKLOOP
            &.or.(klmn.eq.lb1) &
#endif
            &) then
          ! \todo check if the following line is correct
          call z_deriv(p_bar_emfields(:,:,:,sm(klmn),sn(klmn),1),p_dxi_dz,l)
       endif

       !compute the term  
       do j=lj1,lj2
          p_rhs(i1:i2,j,klmn)=p_rhs(i1:i2,j,klmn)+p_dxi_dz(i1:i2,j,sk(klmn))*&
               &pdphidz(i1:i2,pj1,klmn)
       enddo
    enddo 

    PERFOFF_I

  end subroutine add_dchidz_df_adptv
  
end module dchidz_term
