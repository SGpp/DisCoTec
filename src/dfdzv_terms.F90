#include "redef.h"
#include "intrinsic_sizes.h"

!>Computes the parallel advection and trapping terms (i.e. z and v_parallel derivatives of f)
module dfdzv_terms
  use par_mod
  use discretization
  use discretization_adptv_module
  use blockindex
  use prefactors
  use boundaries
  use axpy
  use geometry, only: geom, C_xy, C_y
  use numerical_damping
  use par_geom
  use equilibrium_fields, only: with_comoving_other, phi0, dens_co

  implicit none
  public:: mem_est_dfdzv, initialize_dfdzv, equ_dfdzv, finalize_dfdzv
  public:: add_dfdz, add_dfdv, add_dfdz_block, add_dfdv_block
  public:: add_hypz, add_hypv, add_hypz_block, add_hypv_block
  ! functions necessary for adaptivity
  public:: equ_dfdzv_adptv

  private

  integer :: init_status = 0

  integer :: sten_bound

  logical:: replace_rhs

  ! Module arrays, which are shared on all threads
  integer,dimension(:),allocatable:: map_to_f, f_shifts
  integer,dimension(:,:),allocatable:: stenbound
  integer(kind=8),dimension(:,:), allocatable :: shifted_f_pos
  real,dimension(:,:,:,:),allocatable:: df_pref, trp_co
  complex,dimension(:,:,:,:),allocatable:: df_pref_shm

  logical, dimension(:,:,:,:), allocatable:: zbound1, zbound2, zbound3, zbound4
  complex, dimension(:,:), allocatable:: dfdz_pref3_shm
  real, dimension(:,:), allocatable:: dfdz_pref3
  real, dimension(:,:), allocatable:: dfdv_pref3
  real, dimension(:), allocatable:: v_sten, dfdv_pr
  real,dimension(-2:2):: comb_hyp_v,comb_hyp_z
  real, dimension(:,:,:,:,:),allocatable :: pdf1dz_4, trp_4


  INTERFACE add_dfdz
     MODULE PROCEDURE add_dfdz_block, add_dfdz_indices
  END INTERFACE

  INTERFACE add_dfdv
     MODULE PROCEDURE add_dfdv_block, add_dfdv_indices
  END INTERFACE

  INTERFACE add_hypz
     MODULE PROCEDURE add_hypz_block, add_hypz_indices
  END INTERFACE

  INTERFACE add_hypv
     MODULE PROCEDURE add_hypv_block, add_hypv_indices
  END INTERFACE

contains

  function mem_est_dfdzv(mem_req_in)
    real:: mem_req_in, mem_est_dfdzv
    real:: mem_loc
    
    select case(perf_vec(3))
    case(1)
       mem_loc=mem_est_dfdzv_1()
    case(2)
       mem_loc=mem_est_dfdzv_2()
    case(3)
       mem_loc=mem_est_dfdzv_3()
    case(4)
       mem_loc=mem_est_dfdzv_4()
    end select
    
    mem_est_dfdzv=mem_req_in+mem_loc
    
  end function mem_est_dfdzv


  subroutine initialize_dfdzv(dfdzv_replace_rhs)
    logical:: dfdzv_replace_rhs
    replace_rhs=dfdzv_replace_rhs

    sten_bound = MAX(par_sten_bound,2)

    if ((init_status.ne.perf_vec(3)).and.(init_status.gt.0)) &
         & call finalize_dfdzv

    if (init_status.eq.0) then
       if (with_comoving_other) call initialize_trp_co
       select case(perf_vec(3))
       case(1)
          call initialize_dfdzv_1
       case(2)
          call initialize_dfdzv_2
       case(3)
          call initialize_dfdzv_3
       case(4)
          call initialize_dfdzv_4          
       end select

       init_status = perf_vec(3)
    endif

  end subroutine initialize_dfdzv

  !!\todo dfdzv freezes somehow for 2nd order centered differences
  subroutine equ_dfdzv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    integer:: lbpm, ubpm

    select case(perf_vec(3))
    case(1)
       lbpm=(map_to_f(lb1)-1)*lij0+1
       ubpm=map_to_f(lb2)*lij0
       call equ_dfdzv_1(p_f,p_rhs,lb1,lb2,lbpm,ubpm)
    case(2)
       if(xy_local) then
          call equ_dfdzv_ff_2(p_f,p_rhs,lb1,lb2,df_pref)
       else
          if (shifted_metric) then
             call equ_dfdzv_df_2_shm(p_f,p_rhs,lb1,lb2,df_pref_shm)
          else
             ! standard: no shifted metric
             call equ_dfdzv_df_2(p_f,p_rhs,lb1,lb2,df_pref)
          endif
       endif
    case(3)
       if(xy_local) then
          call equ_dfdzv_ff_3(p_f,p_rhs,lb1,lb2,df_pref)
       else
          if (shifted_metric) then
             call equ_dfdzv_df_3_shm(p_f,p_rhs,lb1,lb2,df_pref_shm)
          else
             call equ_dfdzv_df_3(p_f,p_rhs,lb1,lb2,df_pref)
          endif
       end if
    case(4)
       call equ_dfdzv_ff_4(p_f,p_rhs,lb1,lb2)
    end select
  end subroutine equ_dfdzv

  subroutine finalize_dfdzv

    if (with_comoving_other) call finalize_trp_co

    select case(init_status)
    case(1)
       call finalize_dfdzv_1
    case(2)
       call finalize_dfdzv_2
    case(3)
       call finalize_dfdzv_3
    case(4)
       call finalize_dfdzv_4
    end select

    init_status = 0

  end subroutine finalize_dfdzv

  !different implementations

!***********************************************************************************************

  function mem_est_dfdzv_1()
    real:: mem_est_dfdzv_1

    if (shifted_metric) then
       !pref3's    
       mem_est_dfdzv_1=lij0*lzvwn0*(2*sten_bound+1)*(SIZE_OF_REAL_MB+SIZE_OF_COMPLEX_MB)
    else
       !pref3's    
       mem_est_dfdzv_1=2*lij0*lzvwn0*(2*sten_bound+1)*(SIZE_OF_REAL_MB)
    endif
    !map_to_f
    mem_est_dfdzv_1=mem_est_dfdzv_1+lklmn0*SIZE_OF_INTEGER_MB
    !locrhsg
    mem_est_dfdzv_1=mem_est_dfdzv_1+lij0*lzvwn0*SIZE_OF_COMPLEX_MB

#ifndef oldparbc
    if (xy_local) then
       mem_est_dfdzv_1=mem_est_dfdzv_1+4*lijkl0*SIZE_OF_LOGICAL_MB
    endif
#endif
    
  end function mem_est_dfdzv_1

  !>Initialize additional contribution to trapping term in comoving frame
  subroutine initialize_trp_co
    integer :: i,j,k,n
    real :: dphi0dz, lhs


    allocate(trp_co(pi1:pi2,pj1:pj2,lk1:lk2,ln1:ln2))
 
    do k=lk1,lk2
       do j=pj1,pj2
          do i=pi1,pi2
             dphi0dz = 0.0
             lhs = 0.0
             do n=0, n_spec-1
                dphi0dz = dphi0dz + spec(n)%charge*spec(n)%dens*dens_co(i,k,n)*&
                     spec(n)%mass/spec(n)%temp*Omega0_tor**2*geom%R_hat(i,k)*&
                !dRdz
                (-geom%jacobian(i,j,k)*C_y(i)/geom%R_hat(i,k)*geom%dxdZ(i,k))
                lhs = lhs + spec(n)%charge**2*spec(n)%dens*dens_co(i,k,n)/spec(n)%temp
             enddo
             dphi0dz = dphi0dz / lhs
             do n=ln1,ln2
                trp_co(i,j,k,n) = C_xy(i)/(geom%jacobian(i,j,k)*&
                     &geom%Bfield(i,j,k)*sqrt(2.0*spec(n)%temp/spec(n)%mass))*&
                     &(spec(n)%charge/spec(n)%mass*dphi0dz+Omega0_tor**2*&
                     geom%R_hat(i,k)*(-geom%jacobian(i,j,k)*C_y(i)/&
                     geom%R_hat(i,k)*geom%dxdZ(i,k)))
             enddo
          enddo
       enddo
    end do
    
  end subroutine initialize_trp_co

  subroutine finalize_trp_co
    deallocate(trp_co)
  end subroutine finalize_trp_co

  subroutine initialize_dfdzv_1
    integer:: i,j,k,l,m,n,sten,ijklmn_f,klmn, pii
    real, dimension(pi1:pi2) :: pdf1dz_1
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2) :: trp
    
    if (shifted_metric) then 
       allocate(dfdz_pref3_shm(lij0*lzvwn0,-sten_bound:sten_bound))
       dfdz_pref3_shm=0.
    else
       allocate(dfdz_pref3(lij0*lzvwn0,-sten_bound:sten_bound))
       dfdz_pref3=0.
    endif
    allocate(dfdv_pref3(lij0*lzvwn0,-sten_bound:sten_bound))
    
    allocate(map_to_f(lklmn0))

    dfdv_pref3=0.

#ifndef oldparbc
    if (xy_local) call ini_zbound
#endif
    
    pii = pi1
    do n=ln1,ln2
       do m=lm1,lm2 
          do k=lk1,lk2
             trp(:,pj1,k) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                  &geom%dBdz(pi1:pi2,pj1,k)/(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
             if (with_comoving_other) trp(:,:,k) = trp(:,:,k) + trp_co(:,:,k,n)
          enddo
          do l=ll1,ll2
             do k=lk1,lk2
                klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1

                pdf1dz_1(:) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                        &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))

                do j=lj1,lj2
                   do i=li1,li2
                      if (pmi0.gt.1) pii=i
                      !note: including jfac here will not be possible once we have bpar
                      !parallel dynamics
                      ijklmn_f=(n-ln1)*lij0*lzvw0 + (nwb+m-lm1)*lij0*lz0*lv0 + (nvb+l-ll1)*lij0*lz0 &
                           + (nzb+k-lk1)*lij0 + (j-lj1)*li0 + (i-li1)+1
                      if (xy_local) then
#ifndef oldparbc
                         if (parscheme.eq.'c4th') then
                            if (zbound1(i,j,k,l)) then
                               do sten=0,2
                                  dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*np1coeff(sten)
                               end do
                            elseif (zbound2(i,j,k,l)) then
                               do sten=-1,2
                                  dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*np2coeff(sten)
                               end do
                            elseif (zbound3(i,j,k,l)) then
                               do sten=-2,1
                                  dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*np3coeff(sten)
                               end do
                            elseif (zbound4(i,j,k,l)) then
                               do sten=-2,0
                                  dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*np4coeff(sten)
                               end do
                            else
                               do sten=-par_sten_bound,par_sten_bound
                                  dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*par_sten(sten)
                               end do
                            endif
                         else
#endif
                            do sten=-par_sten_bound,par_sten_bound
                               dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*par_sten(sten)
                            end do
#ifndef oldparbc
                         endif
#endif
                         !add parallel damping
#ifdef oldparbc
                         do sten=-hyp_z_order/2,hyp_z_order/2
                            dfdz_pref3(ijklmn_f,sten)=dfdz_pref3(ijklmn_f,sten)+hyp_z_scl*hyp_z_sten(sten)
                         enddo
#else
                         if (zbound1(i,j,k,l)) then
                         elseif (zbound2(i,j,k,l)) then             
                            do sten=-1,1
                               dfdz_pref3(ijklmn_f,sten)=dfdz_pref3(ijklmn_f,&
                                    &sten)+hyp_z_scl*dc2coeff(sten)
                            end do
                         elseif (zbound3(i,j,k,l)) then
                            do sten=-1,1 
                               dfdz_pref3(ijklmn_f,sten)=dfdz_pref3(ijklmn_f,&
                                    &sten)+hyp_z_scl*dc2coeff(sten)
                            end do
                         elseif (zbound4(i,j,k,l)) then
                         else
                            do sten=-hyp_z_order/2,hyp_z_order/2
                               dfdz_pref3(ijklmn_f,sten)=dfdz_pref3(ijklmn_f,&
                                    &sten)+hyp_z_scl*dc4coeff(sten)
                            end do
                         endif
#endif
                      else
                         if (shifted_metric) then
                            !parallel dynamics
                            do sten=-par_sten_bound,par_sten_bound
                               dfdz_pref3_shm(ijklmn_f,sten)=pdf1dz_1(pii)*&
                                    &par_sten(sten)*geom%phasefac(pii,j,k,sten)
                            enddo
                            do sten=-hyp_z_order/2,hyp_z_order/2
                               dfdz_pref3_shm(ijklmn_f,sten)=&
                                    &dfdz_pref3_shm(ijklmn_f,sten)+hyp_z_scl*&
                                    &hyp_z_sten(sten)*geom%phasefac(pii,j,k,sten)
                            enddo
                         else
                            !parallel dynamics
                            do sten=-par_sten_bound,par_sten_bound
                               dfdz_pref3(ijklmn_f,sten)=pdf1dz_1(pii)*&
                                    &par_sten(sten)
                            enddo
                            do sten=-hyp_z_order/2,hyp_z_order/2
                               dfdz_pref3(ijklmn_f,sten)=dfdz_pref3(ijklmn_f,&
                                    &sten)+hyp_z_scl*hyp_z_sten(sten)
                            enddo
                         endif
                      endif

                      !trapping term
                      if (Erad_acc.and.Erad.ne.0..and..not.y_local) then 
                         do sten=-2,2
                            dfdv_pref3(ijklmn_f,sten)=(trp(pii,pj1,k) - vp(l) * rhostar * &
                                 & Erad / geom%Bfield(pii,pj1, k) * geom%K_j(pii,pj1,k) ) * &
                                 & vderivative(sten,pii,k,l)
                         enddo
                      else
                         do sten=-2,2
                            dfdv_pref3(ijklmn_f,sten)=trp(pii,pj1,k)*&
                                 &vderivative(sten,pii,k,l)
                         enddo
                      endif
                      
                      !v_parallel damping
                      do sten=-2,2
                         dfdv_pref3(ijklmn_f,sten)=dfdv_pref3(ijklmn_f,sten)+&
                              &vdamping(sten,l)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

#ifndef oldparbc
    if (xy_local) deallocate(zbound1, zbound2, zbound3, zbound4)
#endif

    !include v_\parallel boundaries
    do sten=-2,2
       if ((sten.lt.0).and.(my_pev.eq.0)) then
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1,ll1+abs(sten)-1
                   do k=lk1,lk2
                      do j=lj1,lj2
                         do i=li1,li2
                            ijklmn_f=(n-ln1)*lij0*lzvw0 + (nwb+m-lm1)*lij0*&
                                 & lz0*lv0 + (nvb+l-ll1)*lij0*lz0 &
                                 & + (nzb+k-lk1)*lij0 + (j-lj1)*li0 + (i-li1)+1
                            dfdv_pref3(ijklmn_f,sten)=0.
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       elseif((sten.gt.0).and.(my_pev.eq.(n_procs_v-1))) then
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll2-abs(sten)+1,ll2
                   do k=lk1,lk2
                      do j=lj1,lj2
                         do i=li1,li2
                            ijklmn_f=(n-ln1)*lij0*lzvw0 + (nwb+m-lm1)*lij0*&
                                 & lz0*lv0 + (nvb+l-ll1)*lij0*lz0 &
                                 & + (nzb+k-lk1)*lij0 + (j-lj1)*li0 + (i-li1)+1
                            dfdv_pref3(ijklmn_f,sten)=0.
                         enddo
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif
    enddo

  end subroutine initialize_dfdzv_1


  subroutine equ_dfdzv_1(p_f,p_rhs,lb1,lb2,lbp,ubp)
    integer, intent(in):: lb1, lb2, lbp, ubp
    complex, dimension(lij0*lzvwn0), intent(in):: p_f
    complex, dimension(lij0,lb1:lb2), intent(inout):: p_rhs
    complex, dimension(lbp:ubp):: locrhsg

    integer,dimension(-sten_bound:sten_bound):: constr_lb, constr_ub
    integer:: sten, constr
    
    PERFON_I('dfdzv1')
    !z direction 

    if (shifted_metric) then
       sten = -par_sten_bound
       locrhsg(lbp:ubp)= dfdz_pref3_shm(lbp:ubp,sten)*p_f(lbp+sten*lij0:ubp+sten*lij0)
       
       do sten=-par_sten_bound+1,par_sten_bound
          locrhsg(lbp:ubp)= locrhsg(lbp:ubp) + &
               &dfdz_pref3_shm(lbp:ubp,sten)*p_f(lbp+sten*lij0:ubp+sten*lij0)
       enddo
    else
       sten = -par_sten_bound
       locrhsg(lbp:ubp)= dfdz_pref3(lbp:ubp,sten)*p_f(lbp+sten*lij0:ubp+sten*lij0)
       
       do sten=-par_sten_bound+1,par_sten_bound
          locrhsg(lbp:ubp)= locrhsg(lbp:ubp) + &
               &dfdz_pref3(lbp:ubp,sten)*p_f(lbp+sten*lij0:ubp+sten*lij0)
       enddo
    endif

    !boundary condition in v_\parallel
    constr_lb=lbp
    constr_ub=ubp

    if (n_procs_v.le.2) then
       constr=lbp-2*lijz0-1
       if(constr.lt.0) constr_lb(-2)=lbp-constr
       constr=lbp-lijz0-1
       if(constr.lt.0) constr_lb(-1)=lbp-constr
       constr=ubp+lijz0-lij0*lzvwn0
       if(constr.gt.0) constr_ub(1)=ubp-constr
       constr=ubp+2*lijz0-lij0*lzvwn0
       if(constr.gt.0) constr_ub(2)=ubp-constr
    endif

    do sten=-2,2
       locrhsg(constr_lb(sten):constr_ub(sten))= locrhsg(constr_lb(sten):constr_ub(sten))+&
            dfdv_pref3(constr_lb(sten):constr_ub(sten),sten)&
            *p_f(constr_lb(sten)+sten*lijz0:constr_ub(sten)+sten*lijz0)
    enddo

    call rm_ghostcells(locrhsg,p_rhs,lb1,lb2,map_to_f(lb1),map_to_f(lb2))

    PERFOFF_I

  end subroutine equ_dfdzv_1
 
  subroutine rm_ghostcells(loc_rhs,p_rhs,lb1,lb2,flb1,flb2)
    integer,intent(in):: lb1, lb2, flb1, flb2
    complex, dimension(lij0,flb1:flb2), intent(in):: loc_rhs
    complex, dimension(lij0,lb1:lb2), intent(inout):: p_rhs
    integer:: klmn

    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,klmn)=loc_rhs(:,map_to_f(klmn))
       enddo
    else
       do klmn=lb1,lb2
          p_rhs(:,klmn)=p_rhs(:,klmn)+loc_rhs(:,map_to_f(klmn))
       enddo
    endif

  end subroutine rm_ghostcells


  subroutine finalize_dfdzv_1  

    if (shifted_metric) then 
       deallocate(dfdz_pref3_shm)
    else
       deallocate(dfdz_pref3)
    endif
    deallocate(dfdv_pref3, map_to_f)  
    
  end subroutine finalize_dfdzv_1

!***********************************************************************************************
  
  function mem_est_dfdzv_2()
    real:: mem_est_dfdzv_2

    if (shifted_metric) then
       !df_pref
       mem_est_dfdzv_2=pmi0*lj0*(2*par_sten_bound+5)*lklmn0*SIZE_OF_COMPLEX_MB
    else
       !df_pref
       mem_est_dfdzv_2=pmi0*pj0*(2*par_sten_bound+5)*lklmn0*SIZE_OF_REAL_MB
    endif

    !map_to_f, stenbound, f_shifts
    mem_est_dfdzv_2=mem_est_dfdzv_2+1.0*(3*lklmn0+2*&
         &par_sten_bound+5)*SIZE_OF_INTEGER_MB

  end function mem_est_dfdzv_2

  subroutine initialize_dfdzv_2
    integer:: i,j,k,l,m,n,sten,klmn
    real, dimension(pi1:pi2) :: pdf1dz_1
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2) :: trp
    real:: Erad_term 

    if (shifted_metric) then
       allocate(df_pref_shm(pi1:pi2,lj1:lj2,-par_sten_bound-2:par_sten_bound+2,lklmn0))
       df_pref_shm=0.
    else
       allocate(df_pref(pi1:pi2,pj1:pj2,-par_sten_bound-2:par_sten_bound+2,lklmn0))
       df_pref=0.
    endif
    
    allocate(map_to_f(lklmn0),stenbound(2,lklmn0),f_shifts(-par_sten_bound-2:par_sten_bound+2))
    
    do sten=-par_sten_bound,par_sten_bound
       f_shifts(sten)=sten
    enddo
    f_shifts(-par_sten_bound-1)=-lz0
    f_shifts(-par_sten_bound-2)=-2*lz0
    f_shifts(par_sten_bound+1)=lz0
    f_shifts(par_sten_bound+2)=2*lz0
    stenbound(1,:)=-par_sten_bound-2
    stenbound(2,:)=par_sten_bound+2
    
    if (shifted_metric) then
       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                trp(:,pj1,k) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                     &geom%dBdz(pi1:pi2,pj1,k)/(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
                if (with_comoving_other) trp(:,:,k) = trp(:,:,k) + trp_co(:,:,k,n)
             enddo
             do l=ll1,ll2 
                do k=lk1,lk2
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + &
                        & (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1

                   pdf1dz_1(:) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                        &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
                   do j=lj1,lj2
                      do i=pi1,pi2
                         if((Erad.ne.0.).and.(.not.y_local)) then 
                            Erad_term =  - vp(l) * rhostar * Erad / geom%Bfield(i, pj1, k) * geom%K_j(i, pj1, k)
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref_shm(i,j,sten-par_sten_bound,klmn)= &
                                    &(trp(i,pj1,k) + Erad_term)*vderivative(sten,i,k,l)
                               df_pref_shm(i,j,sten-par_sten_bound,klmn)=&
                                    & df_pref_shm(i,j,sten-par_sten_bound,klmn)+&
                                    & hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref_shm(i,j,sten,klmn)= (trp(i,pj1,k) + Erad_term)*vderivative(sten,i,k,l)
                            df_pref_shm(i,j,sten,klmn)=df_pref_shm(i,j,sten,klmn)+&
                                 & hyp_v*dc4coeff(sten)
                            
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref_shm(i,j,sten+par_sten_bound,klmn)=&
                                    & (trp(i,pj1,k) + Erad)*vderivative(sten,i,k,l)
                               df_pref_shm(i,j,sten+par_sten_bound,klmn)=&
                                    & df_pref_shm(i,j,sten+par_sten_bound,klmn)+&
                                    & hyp_v*dc4coeff(sten)
                            enddo
                         else
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref_shm(i,j,sten-par_sten_bound,klmn)= trp(i,pj1,k)*vderivative(sten,i,k,l)
                               
                               df_pref_shm(i,j,sten-par_sten_bound,klmn)=&
                                    & df_pref_shm(i,j,sten-par_sten_bound,klmn)+&
                                    & hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref_shm(i,j,sten,klmn)= trp(i,pj1,k)*vderivative(sten,i,k,l)
                            df_pref_shm(i,j,sten,klmn)=df_pref_shm(i,j,sten,klmn)+&
                                 & hyp_v*dc4coeff(sten)
                         
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref_shm(i,j,sten+par_sten_bound,klmn)=&
                                    & trp(i,pj1,k) *vderivative(sten,i,k,l)
                               df_pref_shm(i,j,sten+par_sten_bound,klmn)=&
                                    & df_pref_shm(i,j,sten+par_sten_bound,klmn)+&
                                    & hyp_v*dc4coeff(sten)
                            enddo
                         endif

                         !parallel dynamics
                         do sten=-par_sten_bound,par_sten_bound
                            df_pref_shm(i,j,sten,klmn)=df_pref_shm(i,j,sten,klmn)+&
                                 &pdf1dz_1(i)*par_sten(sten)*geom%phasefac(i,j,k,sten)
                         enddo
                         !add parallel damping
                         do sten=-hyp_z_order/2,hyp_z_order/2
                            df_pref_shm(i,j,sten,klmn)=df_pref_shm(i,j,sten,klmn)+&
                                 &hyp_z_scl*hyp_z_sten(sten)*geom%phasefac(i,j,k,sten)
                         enddo
                         
                         if((l.eq.ll1).and.(my_pev.eq.0)) stenbound(1,klmn)=-par_sten_bound
                         if((l.eq.ll1+1).and.(my_pev.eq.0)) stenbound(1,klmn)=-par_sten_bound-1
                         if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=par_sten_bound+1
                         if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=par_sten_bound
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                trp(:,pj1,k) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                     &geom%dBdz(pi1:pi2,pj1,k)/(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
                if (with_comoving_other) trp(:,:,k) = trp(:,:,k) + trp_co(:,:,k,n)
             enddo
             do l=ll1,ll2 
                do k=lk1,lk2
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1

                   pdf1dz_1(:) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                        &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))

                   do j=pj1,pj2
                      do i=pi1,pi2
                         if((Erad.ne.0.).and.(.not.y_local)) then 
                            Erad_term = - vp(l) * rhostar * Erad / geom%Bfield(i, pj1, k) * geom%K_j(i, pj1, k)
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref(i,j,sten-par_sten_bound,klmn)=(trp(i,j,k)+Erad_term)*&
                                    &vderivative(sten,i,k,l)
                               df_pref(i,j,sten-par_sten_bound,klmn)=df_pref(i,j,sten-par_sten_bound,klmn)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref(i,j,sten,klmn)=(trp(i,j,k)+Erad_term)*vderivative(sten,i,k,l)
                            df_pref(i,j,sten,klmn)=df_pref(i,j,sten,klmn)+hyp_v*dc4coeff(sten)
                            
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref(i,j,sten+par_sten_bound,klmn)=(trp(i,j,k)+Erad_term)*&
                                    &vderivative(sten,i,k,l)
                               df_pref(i,j,sten+par_sten_bound,klmn)=df_pref(i,j,sten+par_sten_bound,klmn)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                            
                         else
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref(i,j,sten-par_sten_bound,klmn)=trp(i,j,k)*&
                                    &vderivative(sten,i,k,l)
                               df_pref(i,j,sten-par_sten_bound,klmn)=df_pref(i,j,sten-par_sten_bound,klmn)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref(i,j,sten,klmn)=trp(i,j,k)*vderivative(sten,i,k,l)
                            df_pref(i,j,sten,klmn)=df_pref(i,j,sten,klmn)+hyp_v*dc4coeff(sten)
                            
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref(i,j,sten+par_sten_bound,klmn)=trp(i,j,k)*vderivative(sten,i,k,l)
                               df_pref(i,j,sten+par_sten_bound,klmn)=df_pref(i,j,sten+par_sten_bound,klmn)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                         endif

                         !parallel dynamics
                         do sten=-par_sten_bound,par_sten_bound
                            df_pref(i,j,sten,klmn)=df_pref(i,j,sten,klmn)+&
                                 & pdf1dz_1(i)*par_sten(sten)
                         enddo
                         !add parallel damping
                         do sten=-hyp_z_order/2,hyp_z_order/2
                            df_pref(i,j,sten,klmn)=df_pref(i,j,sten,klmn)+hyp_z_scl*hyp_z_sten(sten)
                         enddo
                         
                         if((l.eq.ll1).and.(my_pev.eq.0)) &
                              & stenbound(1,klmn)=-par_sten_bound
                         if((l.eq.ll1+1).and.(my_pev.eq.0)) &
                              & stenbound(1,klmn)=-par_sten_bound-1
                         if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) &
                              & stenbound(2,klmn)=par_sten_bound+1
                         if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) &
                              & stenbound(2,klmn)=par_sten_bound
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif

  end subroutine initialize_dfdzv_2
  
  subroutine equ_dfdzv_ff_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(1:lij0,lzvwn0), intent(in):: p_f
    complex, dimension(1:lij0, lb1:lb2), intent(inout):: p_rhs
    real, dimension(-par_sten_bound-2:par_sten_bound+2,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind

    PERFON_I('dfdzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          p_rhs(:,klmn)=prefac(stenbound(1,klmn),klmn)*p_f(:,f_ind)
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             call axpy_ij(lij0,prefac(sten,klmn),p_f(:,f_ind),p_rhs(:,klmn))
          enddo
       enddo
    else
       do klmn=lb1,lb2
          do sten=stenbound(1,klmn),stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             call axpy_ij(lij0,prefac(sten,klmn),p_f(:,f_ind),p_rhs(:,klmn))
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dfdzv_ff_2

  subroutine equ_dfdzv_df_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-par_sten_bound-2:par_sten_bound+2,lklmn0), intent(in) ::  prefac
    integer:: j, klmn, sten, f_ind

    PERFON_I('dfdzv2')
    ! standard case is replace_rhs, only with collisions and some other 
    ! conditions, this is false
    if(replace_rhs) then
       do klmn=lb1,lb2
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          do j=lj1,lj2
             p_rhs(:,j,klmn)=prefac(:,stenbound(1,klmn),klmn)*p_f(:,j,f_ind)
          enddo
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,sten,klmn)*p_f(:,j,f_ind)
             enddo
          enddo
       enddo
    else
       do klmn=lb1,lb2
          do sten=stenbound(1,klmn),stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,sten,klmn)*p_f(:,j,f_ind)
             enddo
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dfdzv_df_2
  
  subroutine equ_dfdzv_df_2_shm(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-par_sten_bound-2:par_sten_bound+2,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind

    PERFON_I('dfdzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          p_rhs(:,:,klmn)=prefac(:,:,stenbound(1,klmn),klmn)*p_f(:,:,f_ind)
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,sten,klmn)*p_f(:,:,f_ind)
          enddo
       enddo
    else
       do klmn=lb1,lb2
          do sten=stenbound(1,klmn),stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,sten,klmn)*p_f(:,:,f_ind)
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dfdzv_df_2_shm

  subroutine finalize_dfdzv_2
    if (shifted_metric) then
       deallocate(df_pref_shm,map_to_f,stenbound,f_shifts)
    else
       deallocate(df_pref,map_to_f,stenbound,f_shifts)
    endif
  end subroutine finalize_dfdzv_2

!***********************************************************************************************
  
  function mem_est_dfdzv_3()
    real:: mem_est_dfdzv_3
    
    if (shifted_metric) then
       !df_pref_shm
       mem_est_dfdzv_3=pmi0*lj0*(2*par_sten_bound+5)*lklmn0*SIZE_OF_COMPLEX_MB
    else
       !df_pref
       mem_est_dfdzv_3=pmi0*pj0*(2*par_sten_bound+5)*lklmn0*SIZE_OF_REAL_MB
    endif

    !shifted_f_pos
    mem_est_dfdzv_3=mem_est_dfdzv_3+&
         &(2*par_sten_bound+5)*lklmn0*2*SIZE_OF_INTEGER_MB

  end function mem_est_dfdzv_3

  subroutine initialize_dfdzv_3
    integer,dimension(-par_sten_bound-2:par_sten_bound+2):: f_shift
    integer:: i,j,k,l,m,n,sten,klmn
    integer(kind=8) :: f_index
    real, dimension(pi1:pi2) :: pdf1dz_1
    real, dimension(pi1:pi2,pj1:pj2,lk1:lk2) :: trp
     real:: Erad_term 

    if (shifted_metric) then
       allocate(df_pref_shm(pi1:pi2,lj1:lj2,lklmn0,-par_sten_bound-2:par_sten_bound+2))
       df_pref_shm=0.
    else
       allocate(df_pref(pi1:pi2,pj1:pj2,lklmn0,-par_sten_bound-2:par_sten_bound+2))
       df_pref=0.
    endif
    
    allocate(shifted_f_pos(lklmn0,-par_sten_bound-2:par_sten_bound+2))
    shifted_f_pos=0

    do sten=-par_sten_bound,par_sten_bound
       f_shift(sten)=sten
    enddo
    f_shift(-par_sten_bound-1)=-lz0
    f_shift(-par_sten_bound-2)=-2*lz0
    f_shift(par_sten_bound+1)=lz0
    f_shift(par_sten_bound+2)=2*lz0
    
    if (shifted_metric) then
       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                trp(:,pj1,k) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                     &geom%dBdz(pi1:pi2,pj1,k)/(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
                if (with_comoving_other) trp(:,:,k) = trp(:,:,k) + trp_co(:,:,k,n)
             enddo
             do l=ll1,ll2 
                do k=lk1,lk2 
                   pdf1dz_1(:) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                        &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k)) 

                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   f_index= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   do sten=-par_sten_bound-2,par_sten_bound+2
                      shifted_f_pos(klmn,sten)=f_index+f_shift(sten)
                      if(shifted_f_pos(klmn,sten).lt.1) shifted_f_pos(klmn,sten)=1
                      if(shifted_f_pos(klmn,sten).gt.lzvwn0) shifted_f_pos(klmn,sten)=lzvwn0
                   end do
                   do j=lj1,lj2
                      do i=pi1,pi2
                         if((Erad.ne.0.).and.(.not.y_local)) then 
                            Erad_term =  - vp(l) * rhostar * Erad / geom%Bfield(i, pj1, k) * geom%K_j(i, pj1, k)
                            !trapping term + damping (negative shifts)
                               do sten=-2,-1
                                  df_pref_shm(i,j,klmn,sten-par_sten_bound)=&
                                       & (trp(i,pj1,k) + Erad_term) *vderivative(sten,i,k,l)
                                  df_pref_shm(i,j,klmn,sten-par_sten_bound)=&
                                       & df_pref_shm(i,j,klmn,sten-par_sten_bound)+&
                                       & hyp_v*dc4coeff(sten)
                               enddo
                               !trapping term + damping (zero shift)
                               sten=0
                               df_pref_shm(i,j,klmn,sten)=(trp(i,pj1,k)+Erad_term)*vderivative(sten,i,k,l)
                               df_pref_shm(i,j,klmn,sten)=df_pref_shm(i,j,klmn,sten)+&
                                    & hyp_v*dc4coeff(sten)
                               
                               !trapping term + damping (positive shifts)
                               do sten=1,2
                                  df_pref_shm(i,j,klmn,sten+par_sten_bound)= &
                                       & (trp(i,pj1,k)+Erad_term)*vderivative(sten,i,k,l)
                                  df_pref_shm(i,j,klmn,sten+par_sten_bound)= &
                                       & df_pref_shm(i,j,klmn,sten+par_sten_bound)+ &
                                       & hyp_v*dc4coeff(sten)
                               enddo
                            else

                               !trapping term + damping (negative shifts)
                               do sten=-2,-1
                                  df_pref_shm(i,j,klmn,sten-par_sten_bound)=&
                                       & trp(i,pj1,k)*vderivative(sten,i,k,l)
                                  df_pref_shm(i,j,klmn,sten-par_sten_bound)=&
                                       & df_pref_shm(i,j,klmn,sten-par_sten_bound)+&
                                       & hyp_v*dc4coeff(sten)
                               enddo
                               !trapping term + damping (zero shift)
                               sten=0
                               df_pref_shm(i,j,klmn,sten)=trp(i,pj1,k)*vderivative(sten,i,k,l)
                               df_pref_shm(i,j,klmn,sten)=df_pref_shm(i,j,klmn,sten)+&
                                    & hyp_v*dc4coeff(sten)
                               
                               !trapping term + damping (positive shifts)
                               do sten=1,2
                                  df_pref_shm(i,j,klmn,sten+par_sten_bound)= &
                                       & trp(i,pj1,k)*vderivative(sten,i,k,l)
                                  df_pref_shm(i,j,klmn,sten+par_sten_bound)= &
                                       & df_pref_shm(i,j,klmn,sten+par_sten_bound)+ &
                                       & hyp_v*dc4coeff(sten)
                               enddo
                         endif

                         !parallel dynamics
                         do sten=-par_sten_bound,par_sten_bound
                            df_pref_shm(i,j,klmn,sten)= & 
                                 & df_pref_shm(i,j,klmn,sten)+pdf1dz_1(i)*&
                                 & par_sten(sten)*geom%phasefac(i,j,k,sten)
                         enddo
                         !add parallel damping
                         do sten=-hyp_z_order/2,hyp_z_order/2
                            df_pref_shm(i,j,klmn,sten)=df_pref_shm(i,j,klmn,sten)+&
                                 & hyp_z_scl*hyp_z_sten(sten)*geom%phasefac(i,j,k,sten)
                         enddo
                         
                         !boundary conditions in v_parallel
                         if((l.eq.ll1).and.(my_pev.eq.0)) df_pref_shm(i,j,klmn,-par_sten_bound-2:-par_sten_bound-1)=0.
                         if((l.eq.ll1+1).and.(my_pev.eq.0)) df_pref_shm(i,j,klmn,-par_sten_bound-2)=0.
                         if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) df_pref_shm(i,j,klmn,par_sten_bound+2)=0.
                         if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) df_pref_shm(i,j,klmn,par_sten_bound+1:par_sten_bound+2)=0.
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do n=ln1,ln2
          do m=lm1,lm2
             do k=lk1,lk2
                trp(:,pj1,k) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                     &geom%dBdz(pi1:pi2,pj1,k)/(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
                if (with_comoving_other) trp(:,:,k) = trp(:,:,k) + trp_co(:,:,k,n)
             enddo
             do l=ll1,ll2 
                do k=lk1,lk2                    
                   pdf1dz_1(:) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                        &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k)) 

                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   f_index= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   do sten=-par_sten_bound-2,par_sten_bound+2
                      shifted_f_pos(klmn,sten)=f_index+f_shift(sten)
                      if(shifted_f_pos(klmn,sten).lt.1) shifted_f_pos(klmn,sten)=1
                      if(shifted_f_pos(klmn,sten).gt.lzvwn0) shifted_f_pos(klmn,sten)=lzvwn0
                   end do
                   do j=pj1,pj2
                      do i=pi1,pi2
                         if((Erad.ne.0.).and.(.not.y_local)) then 
                            Erad_term =  - vp(l) * rhostar * Erad / geom%Bfield(i, pj1, k) * geom%K_j(i, pj1, k)
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref(i,j,klmn,sten-par_sten_bound)=(trp(i,j,k)+Erad_term)*&
                                    &vderivative(sten,i,k,l)
                               df_pref(i,j,klmn,sten-par_sten_bound)=df_pref(i,j,klmn,sten-par_sten_bound)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref(i,j,klmn,sten)=(trp(i,j,k)+Erad_term)*vderivative(sten,i,k,l)
                            df_pref(i,j,klmn,sten)=df_pref(i,j,klmn,sten)+hyp_v*dc4coeff(sten)
                            
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref(i,j,klmn,sten+par_sten_bound)=(trp(i,j,k)+Erad_term)*&
                                    &vderivative(sten,i,k,l)
                               df_pref(i,j,klmn,sten+par_sten_bound)=df_pref(i,j,klmn,sten+par_sten_bound)+&
                                    &hyp_v*dc4coeff(sten)
                            enddo
                         else
                            
                            !trapping term + damping (negative shifts)
                            do sten=-2,-1
                               df_pref(i,j,klmn,sten-par_sten_bound)=trp(i,j,k)*vderivative(sten,i,k,l)
                               df_pref(i,j,klmn,sten-par_sten_bound)=df_pref(i,j,klmn,sten-par_sten_bound)+hyp_v*dc4coeff(sten)
                            enddo
                            !trapping term + damping (zero shift)
                            sten=0
                            df_pref(i,j,klmn,sten)=trp(i,j,k)*vderivative(sten,i,k,l)
                            df_pref(i,j,klmn,sten)=df_pref(i,j,klmn,sten)+hyp_v*dc4coeff(sten)
                            
                            !trapping term + damping (positive shifts)
                            do sten=1,2
                               df_pref(i,j,klmn,sten+par_sten_bound)=trp(i,j,k)*vderivative(sten,i,k,l)
                               df_pref(i,j,klmn,sten+par_sten_bound)=df_pref(i,j,klmn,sten+par_sten_bound)+hyp_v*dc4coeff(sten)
                            enddo
                         
                      endif

                         !parallel dynamics
                         do sten=-par_sten_bound,par_sten_bound
                            df_pref(i,j,klmn,sten)=df_pref(i,j,klmn,sten)+pdf1dz_1(i)*par_sten(sten)
                         enddo
                         !add parallel damping
                         do sten=-hyp_z_order/2,hyp_z_order/2
                            df_pref(i,j,klmn,sten)=df_pref(i,j,klmn,sten)+hyp_z_scl*hyp_z_sten(sten)
                         enddo
                         
                         !boundary conditions in v_parallel
                         if((l.eq.ll1).and.(my_pev.eq.0)) df_pref(i,j,klmn,-par_sten_bound-2:-par_sten_bound-1)=0.
                         if((l.eq.ll1+1).and.(my_pev.eq.0)) df_pref(i,j,klmn,-par_sten_bound-2)=0.
                         if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) df_pref(i,j,klmn,par_sten_bound+2)=0.
                         if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) df_pref(i,j,klmn,par_sten_bound+1:par_sten_bound+2)=0.
                         
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif

  end subroutine initialize_dfdzv_3

  subroutine equ_dfdzv_ff_3(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lij0,lzvwn0), intent(in):: p_f
    complex, dimension(lij0, lb1:lb2), intent(inout):: p_rhs
    real, dimension(lklmn0,-par_sten_bound-2:par_sten_bound+2), intent(in) ::  prefac

    integer:: klmn,sten

    PERFON_I('dfdzv3')
    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,klmn)=prefac(klmn,-par_sten_bound-2)*&
               &p_f(:,shifted_f_pos(klmn,-par_sten_bound-2))
       enddo
       do sten=-par_sten_bound-1,par_sten_bound+2
          do klmn=lb1,lb2
             call axpy_ij(lij0,prefac(klmn,sten),p_f(:,&
                  &shifted_f_pos(klmn,sten)),p_rhs(:,klmn))
          enddo
       enddo
    else
       do sten=-par_sten_bound-2,par_sten_bound+2
          do klmn=lb1,lb2
             call axpy_ij(lij0,prefac(klmn,sten),p_f(:,shifted_f_pos(klmn,sten)),p_rhs(:,klmn))
          enddo
       enddo
    endif
    PERFOFF_I
  end subroutine equ_dfdzv_ff_3

  subroutine equ_dfdzv_df_3(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lklmn0,-par_sten_bound-2:par_sten_bound+2), intent(in) ::  prefac

    integer:: j,klmn,sten

    PERFON_I('dfdzv3')
    if(replace_rhs) then
       do klmn=lb1,lb2
          do j=lj1,lj2
             p_rhs(:,j,klmn)=prefac(:,klmn,-par_sten_bound-2)*&
                  &p_f(:,j,shifted_f_pos(klmn,-par_sten_bound-2))
          enddo
       enddo
       do sten=-par_sten_bound-1,par_sten_bound+2
          do klmn=lb1,lb2
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,klmn,sten)*p_f(:,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    else
       do sten=-par_sten_bound-2,par_sten_bound+2
          do klmn=lb1,lb2
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,klmn,sten)*p_f(:,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dfdzv_df_3

  subroutine equ_dfdzv_df_3_shm(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,lklmn0,-par_sten_bound-2:par_sten_bound+2), intent(in) ::  prefac 

    integer:: klmn,sten

    PERFON_I('dfdzv3')
    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,:,klmn)=prefac(:,:,klmn,-par_sten_bound-2)*&
               &p_f(:,:,shifted_f_pos(klmn,-par_sten_bound-2))
       enddo
       do sten=-par_sten_bound-1,par_sten_bound+2
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    else
       do sten=-par_sten_bound-2,par_sten_bound+2
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dfdzv_df_3_shm

  subroutine finalize_dfdzv_3
    if (shifted_metric) then
       deallocate(df_pref_shm,shifted_f_pos)
    else
       deallocate(df_pref,shifted_f_pos)
    endif
  end subroutine finalize_dfdzv_3
 

!***********************************************************************************************

  subroutine ini_zbound
    integer:: i,j,k,l
    real :: pdf1dz_m1n1
    allocate(zbound1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),zbound2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),&
         zbound3(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),zbound4(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2))
    zbound1=.false.
    zbound2=.false.
    zbound3=.false.
    zbound4=.false.

    do l=ll1,ll2
       do k=lk1,lk2
          pdf1dz_m1n1 = - C_xy(pi1)*vTvpar(l,ln1)/&
               &(geom%jacobian(pi1,pj1,k)*geom%Bfield(pi1,pj1,k))
          do j=lj1,lj2
             do i=li1,li2              
                if ((k.eq.lk1).and.(my_pez.eq.0).and.delete_entry_lower(i,j)) then
                   if(pdf1dz_m1n1.gt.0) then
                      zbound1(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.(lk1+1)).and.(my_pez.eq.0).and.delete_entry_lower(i,j)) then             
                   if(pdf1dz_m1n1.gt.0) then
                      zbound2(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.(lk2-1)).and.(my_pez.eq.(n_procs_z-1)).and.delete_entry_upper(i,j)) then
                   if(pdf1dz_m1n1.lt.0) then
                      zbound3(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.lk2).and.(my_pez.eq.(n_procs_z-1)).and.delete_entry_upper(i,j)) then
                   if(pdf1dz_m1n1.lt.0) then
                      zbound4(i,j,k,l)=.true.
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine ini_zbound

!***********************************************************************************************
!The remaining part of this module contains the 4 different terms separated in subroutines 
!(currently only considering the "local" approximation)
!It should not be used in the time scheme (bad performance) but may prove useful, 
!e.g., for energy diagnostics or for checking against parallel upwind schemes which are only
!availabe here
!For testing purposes, i.e. comparison with the other implementations, it can be activated as 
!perf_vec(3) = 4
  function mem_est_dfdzv_4()
    real:: mem_est_dfdzv_4

    !pdf1dz_4
    mem_est_dfdzv_4=pmi0*pj0*lk0*ll0*ln0*SIZE_OF_REAL_MB

    !trp_4
    mem_est_dfdzv_4=mem_est_dfdzv_4 + pmi0*pj0*lk0*lm0*ln0*SIZE_OF_REAL_MB
  end function mem_est_dfdzv_4

  subroutine initialize_dfdzv_4
    integer :: j,k,l,m,n
    allocate(pdf1dz_4(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, ln1:ln2))
    allocate(trp_4(pi1:pi2, pj1:pj2, lk1:lk2, lm1:lm2, ln1:ln2))

    if(Erad.ne.0.) then 
       if(mype.eq.0) WRITE(*, "(A)") 'No Erad allowed in perf_vec(3)=4, exiting the code ...'
       stop
    endif

    do n=ln1,ln2
       do l=ll1,ll2
          do k=lk1,lk2
             do j=pj1,pj2
                pdf1dz_4(:,j,k,l,n) = - C_xy(pi1:pi2)*vTvpar(l,n)/&
                     &(geom%jacobian(pi1:pi2,pj1,k)*geom%Bfield(pi1:pi2,pj1,k))
             enddo
          enddo
       enddo
    enddo

    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             do j=pj1,pj2
                trp_4(:,j,k,m,n) = C_xy(pi1:pi2)*mu(m)/sqrt(2.*spec(n)%mass/spec(n)%temp)*&
                     &geom%dBdz(pi1:pi2,j,k)/(geom%jacobian(pi1:pi2,j,k)*&
                     &geom%Bfield(pi1:pi2,j,k))
                if (with_comoving_other) trp_4(:,:,k,m,n) = trp_4(:,:,k,m,n) + trp_co(:,:,k,n)
             enddo
          enddo
       enddo
    enddo
  end subroutine initialize_dfdzv_4

  subroutine add_dfdz_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k,l,m,n
    integer :: j, sten, stensign
#ifndef oldparbc
    integer :: i,pii
#endif

    select case(parscheme)
#ifndef oldparbc
    case('c4th') !newparbc only implemented for c4th
       if(.not.allocated(zbound1)) call ini_zbound
       pii = pi1
       do j=lj1,lj2
          do i=li1,li2
             if (pi0.gt.1) pii=i
             if (zbound1(i,j,k,l)) then
                do sten=0,2
                   localrhs(i,j)=localrhs(i,j)+np1coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound2(i,j,k,l)) then
                do sten=-1,2
                   localrhs(i,j)=localrhs(i,j)+np2coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound3(i,j,k,l)) then
                do sten=-2,1
                   localrhs(i,j)=localrhs(i,j)+np3coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound4(i,j,k,l)) then
                do sten=-2,0
                   localrhs(i,j)=localrhs(i,j)+np4coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             else
                do sten=-par_sten_bound,par_sten_bound
                   localrhs(i,j)=localrhs(i,j)+par_sten(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             endif
          enddo
       end do
#endif
    case('u3rd')
       ! Upwind in parallel direction considering sign of parallel velocity
       stensign = sign(1.0,vp(l))
       if (xy_local) then
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             call axpy_ij(lij0,pdf1dz_4(pi1,pj1,k,l,n)*stensign*par_sten(sten),&
                  &p_f_(:,:,k+stensign*sten,l,m,n),localrhs)
          End Do
       else
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             do j=lj1,lj2
                localrhs(li1:li2,j) = localrhs(li1:li2,j)+stensign*par_sten(sten)*&
                     &pdf1dz_4(pi1:pi2,pj1,k,l,n)*p_f_(li1:li2,j,k+stensign*sten,l,m,n)
             enddo
          End Do
       endif
    case default
       if (xy_local) then
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             call axpy_ij(lij0,pdf1dz_4(pi1,pj1,k,l,n)*par_sten(sten),p_f_(:,:,k+sten,l,m,n),&
                  localrhs)
          End Do
       else
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             do j=lj1,lj2
                localrhs(li1:li2,j) = localrhs(li1:li2,j)+&
                     &pdf1dz_4(pi1:pi2,pj1,k,l,n)*par_sten(sten)*p_f_(li1:li2,j,k+sten,l,m,n)
             enddo
          End Do
       endif
    end select

#ifndef oldparbc
    if(allocated(zbound1)) deallocate(zbound1,zbound2,zbound3,zbound4)
#endif

  end subroutine add_dfdz_indices

  subroutine add_dfdz_block(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_dfdz_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_dfdz_block

  
  !> add hyper diffusion in z direction
  !! by adding hyp_z_order'th derivative of f1 (not g1)  
  subroutine add_hypz_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: p_f_
    complex, dimension(li1:li2,lj1:lj2),intent(inout) :: localrhs
    integer,intent(in):: k, l, m, n
    
    integer :: sten,i,j, pii
    
#ifdef oldparbc
    if (shifted_metric) then
       pii=pi1
       do sten=-hyp_z_order/2,hyp_z_order/2
          do j=lj1,lj2
             do i=li1,li2
                if (pmi0.gt.1) pii=i
                localrhs(i,j)=localrhs(i,j)+&
                     & hyp_z_scl*hyp_z_sten(sten)*geom%phasefac(pii,j,k,sten)*p_f_(i,j,k+sten,l,m,n)
             enddo
          enddo
       enddo
    else
       Do sten=-hyp_z_order/2,hyp_z_order/2
          call axpy_ij(lij0,hyp_z_scl*hyp_z_sten(sten),p_f_(:,:,k+sten,l,m,n),localrhs)
       End Do
    endif
#else
    if(.not.allocated(zbound1)) call ini_zbound
    do j=lj1,lj2
       do i=li1,li2
          if (zbound1(i,j,k,l)) then
          elseif (zbound2(i,j,k,l)) then
             do sten=-1,1
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc2coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          elseif (zbound3(i,j,k,l)) then
             do sten=-1,1
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc2coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          elseif (zbound4(i,j,k,l)) then
          else
             do sten=-2,2
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc4coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          endif
       enddo
    enddo
    if(allocated(zbound1)) deallocate(zbound1,zbound2,zbound3,zbound4)
#endif    
    
  end subroutine add_hypz_indices

  subroutine add_hypz_block(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypz_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypz_block


  subroutine add_dfdv_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n
    integer :: sten, j, kb, lb, ub
    !4th order derivatives with boundaries set to 0
    
    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)
    if (xy_local) then
       do sten=lb,ub
          call axpy_ij(lij0,trp_4(pi1,pj1,k,m,n) * vderivative(sten,pi1,k,l),&
               &p_f_(:,:,k,l+sten,m,n),localrhs)
       enddo
    else
       do sten=lb,ub
          do j=lj1,lj2
             localrhs(:,j) = localrhs(:,j) + &
                  trp_4(:,pj1,k,m,n) * vderivative(sten,:,k,l)*p_f_(:,j,k,l+sten,m,n)
          enddo
       enddo
    endif
  end subroutine add_dfdv_indices

  subroutine add_dfdv_block(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_dfdv_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_dfdv_block  


  subroutine add_hypv_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n

    integer :: sten, kb, lb, ub

    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)
    do sten=lb,ub
       call axpy_ij(lij0,vdamping(sten,l),p_f_(:,:,k,l+sten,m,n),localrhs)
    enddo
    
  end subroutine add_hypv_indices

  subroutine add_hypv_block(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypv_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypv_block


  subroutine equ_dfdzv_ff_4(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    if (replace_rhs) p_rhs = 0.0

    do klmn=lb1,lb2
       call add_dfdz_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_hypz_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_dfdv_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_hypv_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine equ_dfdzv_ff_4

  subroutine finalize_dfdzv_4
    deallocate(pdf1dz_4)
    deallocate(trp_4)
  end subroutine finalize_dfdzv_4


! ------ next follow mirrors necessary for adaptive grid routines ---

subroutine equ_dfdzv_adptv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs

    select case(perf_vec(3))
    case(1)
       ! stop "this scheme perf_vec(3) doesn't work yet!"
       if(shifted_metric) then
          call equ_dfdzv_1_adptv(p_f,p_rhs,lb1,lb2,dfdv_pref3,dfdz_pref3_shm=dfdz_pref3_shm)
       else
          call equ_dfdzv_1_adptv(p_f,p_rhs,lb1,lb2,dfdv_pref3,dfdz_pref3=dfdz_pref3)
       endif
    case(2)
       if(xy_local) then
          stop "no adaptivity for local xy version!"
       else
          if (shifted_metric) then
             call equ_dfdzv_df_2_shm_adptv(p_f,p_rhs,lb1,lb2,df_pref_shm)
          else
             ! standard: no shifted metric
             call equ_dfdzv_df_2_adptv(p_f,p_rhs,lb1,lb2,df_pref)
          endif
       endif
    case(3)
       if(xy_local) then
          stop "no adaptivity for local xy version!"
       else
          if (shifted_metric) then
             call equ_dfdzv_df_3_shm_adptv(p_f,p_rhs,lb1,lb2,df_pref_shm)
          else
             call equ_dfdzv_df_3_adptv(p_f,p_rhs,lb1,lb2,df_pref)
          endif
       end if
    case(4)
       call equ_dfdzv_ff_4_adptv(p_f,p_rhs,lb1,lb2)
    end select
  end subroutine equ_dfdzv_adptv


  subroutine equ_dfdzv_1_adptv(p_f,p_rhs,lb1,lb2,dfdv_pref3,dfdz_pref3_shm,dfdz_pref3)
    integer, intent(in):: lb1, lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lj1:lj2,lzvwn0,-sten_bound:sten_bound), intent(in) :: dfdv_pref3
    complex, dimension(li1:li2,lj1:lj2,lzvwn0,-sten_bound:sten_bound), optional, intent(in) :: dfdz_pref3_shm
    real, dimension(li1:li2,lj1:lj2,lzvwn0,-sten_bound:sten_bound), optional, intent(in) :: dfdz_pref3
    integer,dimension(-sten_bound:sten_bound):: constr_lb1, constr_lb2
    integer:: sten, constr, klmn, f_ind, l, m, i1, i2
    
    PERFON_I('dfdzv1')
    !z direction

    if(shifted_metric) then
       do klmn = lb1, lb2
          !>todo ask if it is correct
          f_ind = map_to_f(klmn)
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          sten = -par_sten_bound
          p_rhs(i1:i2,:,klmn) = &
               & dfdz_pref3_shm(i1:i2,:,f_ind,sten)* &
               & p_f(i1:i2,:,f_ind+sten)

          do sten=-par_sten_bound+1,par_sten_bound
             p_rhs(i1:i2,:,klmn) = &
                  & p_rhs(i1:i2,:,klmn) + &
                  & dfdz_pref3_shm(i1:i2,:,f_ind,sten)* &
                  & p_f(i1:i2,:,f_ind+sten)
          enddo
       enddo
    else
       do klmn = lb1, lb2
          !>todo ask if it is correct
          f_ind = map_to_f(klmn)
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          sten = -par_sten_bound
          p_rhs(i1:i2,:,klmn) = &
               & dfdz_pref3(i1:i2,:,f_ind,sten)* &
               & p_f(i1:i2,:,f_ind+sten)

          do sten=-par_sten_bound+1,par_sten_bound
             p_rhs(i1:i2,:,klmn) = &
                  & p_rhs(i1:i2,:,klmn) + &
                  & dfdz_pref3(i1:i2,:,f_ind,sten)* &
                  & p_f(i1:i2,:,f_ind+sten)
          enddo
       enddo
    endif
       

    !boundary condition in v_\parallel
    constr_lb1=lb1
    constr_lb2=lb2

    if (n_procs_v.le.2) then
       constr = lb1 - 2*lk0 - 1
       if(constr.lt.0) constr_lb1(-2) = lb1-constr
       constr = lb1 - lk0-1
       if(constr.lt.0) constr_lb1(-1) = lb1-constr
       constr = lb2 + lk0 - lklmn0
       if(constr.gt.0) constr_lb2(1) = lb2-constr
       constr = lb2 + 2*lk0 - lklmn0
       if(constr.gt.0) constr_lb2(2) = lb2-constr
    endif

    do sten=-2,2
       do klmn = constr_lb1(sten), constr_lb2(sten)
          !>todo ask if it is correct
          f_ind = map_to_f(klmn)
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p_rhs(i1:i2,:,klmn) = &
               & p_rhs(i1:i2,:,klmn) + &
               & dfdv_pref3(i1:i2,:,f_ind,sten)* &
               & p_f(i1:i2,:,f_ind+sten*lz0)
       enddo
    enddo

    PERFOFF_I

  end subroutine equ_dfdzv_1_adptv

  subroutine equ_dfdzv_df_2_shm_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-par_sten_bound-2:par_sten_bound+2,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind, l, m, i1, i2

    PERFON_I('dfdzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          
          p_rhs(i1:i2,:,klmn) = &
               & prefac(i1:i2,:,stenbound(1,klmn),klmn)* &
               & p_f(i1:i2,:,f_ind)
          
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             p_rhs(i1:i2,:,klmn) = &
                  & p_rhs(i1:i2,:,klmn) + &
                  & prefac(i1:i2,:,sten,klmn)*p_f(i1:i2,:,f_ind)
          enddo
       enddo
    else
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do sten=stenbound(1,klmn),stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             p_rhs(i1:i2,:,klmn) = &
                  & p_rhs(i1:i2,:,klmn) + &
                  & prefac(i1:i2,:,sten,klmn)*p_f(i1:i2,:,f_ind)
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dfdzv_df_2_shm_adptv

  subroutine equ_dfdzv_df_2_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-par_sten_bound-2:par_sten_bound+2,lklmn0), intent(in) ::  prefac
    integer:: j, klmn, sten, f_ind, l, m, i1, i2

    PERFON_I('dfdzv2')
    ! standard case is replace_rhs, only with collisions and some other 
    ! conditions, this is false
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          do j=lj1,lj2
             p_rhs(i1:i2,j,klmn)=prefac(i1:i2,stenbound(1,klmn),klmn)* &
             & p_f(i1:i2,j,f_ind)
          enddo
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn)=p_rhs(i1:i2,j,klmn) + &
                & prefac(i1:i2,sten,klmn)*p_f(i1:i2,j,f_ind)
             enddo
          enddo
       enddo
    else
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do sten=stenbound(1,klmn),stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                     & prefac(i1:i2,sten,klmn)*p_f(i1:i2,j,f_ind)
             enddo
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dfdzv_df_2_adptv

    subroutine equ_dfdzv_df_3_shm_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,lklmn0,-par_sten_bound-2:par_sten_bound+2), intent(in) ::  prefac 

    integer:: klmn,sten, l, m, i1, i2

    PERFON_I('dfdzv3')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p_rhs(i1:i2,:,klmn)=prefac(i1:i2,:,klmn,-par_sten_bound-2)*&
               &p_f(i1:i2,:,shifted_f_pos(klmn,-par_sten_bound-2))
       enddo
       do sten=-par_sten_bound-1,par_sten_bound+2
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p_rhs(i1:i2,:,klmn)=p_rhs(i1:i2,:,klmn) + &
                  & prefac(i1:i2,:,klmn,sten)* &
                  & p_f(i1:i2,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    else
       do sten=-par_sten_bound-2,par_sten_bound+2
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p_rhs(i1:i2,:,klmn)=p_rhs(i1:i2,:,klmn) + &
                  & prefac(i1:i2,:,klmn,sten)* &
                  & p_f(i1:i2,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    endif    
    PERFOFF_I

  end subroutine equ_dfdzv_df_3_shm_adptv

  subroutine equ_dfdzv_df_3_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lklmn0,-par_sten_bound-2:par_sten_bound+2), intent(in) ::  prefac

    integer:: j, klmn, sten, l, m, i1, i2

    PERFON_I('dfdzv3')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do j=lj1,lj2
             p_rhs(i1:i2,j,klmn)=prefac(i1:i2,klmn,-par_sten_bound-2)*&
                  &p_f(i1:i2,j,shifted_f_pos(klmn,-par_sten_bound-2))
          enddo
       enddo
       do sten=-par_sten_bound-1,par_sten_bound+2
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn)=p_rhs(i1:i2,j,klmn) + &
                     & prefac(i1:i2,klmn,sten)* &
                     & p_f(i1:i2,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    else
       do sten=-par_sten_bound-2,par_sten_bound+2
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn)=p_rhs(i1:i2,j,klmn) + &
                     & prefac(i1:i2,klmn,sten)* &
                     & p_f(i1:i2,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dfdzv_df_3_adptv

  subroutine equ_dfdzv_ff_4_adptv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    if (replace_rhs) p_rhs = 0.0

    do klmn=lb1,lb2
       call add_dfdz_indices_adptv(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_hypz_indices_adptv(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_dfdv_indices_adptv(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       call add_hypv_indices_adptv(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine equ_dfdzv_ff_4_adptv

  ! --- equ_dfdzv_ff_4_adptv subroutines ---

  subroutine add_dfdz_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k,l,m,n
    integer :: j, sten, stensign, i1, i2, lij0_vwadp
#ifndef oldparbc
    integer :: i, pii
#endif

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    lij0_vwadp = li0_vwadp(l,m)*lj0

    select case(parscheme)
#ifndef oldparbc
    case('c4th')
       if(.not.allocated(zbound1)) call ini_zbound_adptv
       pii = pi1
       do j=lj1,lj2
          do i=i1,i2
             if (pi0.gt.1) pii = i
             if (zbound1(i,j,k,l)) then
                do sten=0,2
                   localrhs(i,j)=localrhs(i,j)+np1coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound2(i,j,k,l)) then
                do sten=-1,2
                   localrhs(i,j)=localrhs(i,j)+np2coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound3(i,j,k,l)) then
                do sten=-2,1
                   localrhs(i,j)=localrhs(i,j)+np3coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             elseif (zbound4(i,j,k,l)) then
                do sten=-2,0
                   localrhs(i,j)=localrhs(i,j)+np4coeff(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             else
                do sten=-par_sten_bound,par_sten_bound
                   localrhs(i,j)=localrhs(i,j)+par_sten(sten)*pdf1dz_4(pii,pj1,k,l,n)*p_f_(i,j,k+sten,l,m,n)
                end do
             endif
          enddo
       end do
#endif
    case('u3rd')
       ! Upwind in parallel direction considering sign of parallel velocity
       stensign = sign(1.0,vp(l))          
       if (xy_local) then
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             call axpy_ij(lij0_vwadp,pdf1dz_4(pi1,pj1,k,l,n)*stensign*par_sten(sten),&
                  &p_f_(:,:,k+stensign*sten,l,m,n),localrhs)
          End Do
       else
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             do j=lj1,lj2
                localrhs(i1:i2,j) = localrhs(i1:i2,j)+stensign*par_sten(sten)*&
                     &pdf1dz_4(i1:i2,pj1,k,l,n)*p_f_(i1:i2,j,k+stensign*sten,l,m,n)
             enddo
          End Do
       endif
    case default
       if (xy_local) then
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             call axpy_ij(lij0_vwadp,pdf1dz_4(pi1,pj1,k,l,n)*par_sten(sten),&
                  &p_f_(:,:,k+sten,l,m,n),localrhs)
          End Do
       else
          Do sten=-par_sten_bound,par_sten_bound
             if (par_sten(sten).eq.0.0) cycle
             do j=lj1,lj2
                localrhs(i1:i2,j) = localrhs(i1:i2,j)+par_sten(sten)*&
                     &pdf1dz_4(i1:i2,pj1,k,l,n)*p_f_(i1:i2,j,k+sten,l,m,n)
             enddo
          End Do
       endif
    end select

#ifndef oldparbc
    if(allocated(zbound1)) deallocate(zbound1,zbound2,zbound3,zbound4)
#endif

  end subroutine add_dfdz_indices_adptv

  !> add hyper diffusion in z direction
  !! by adding hyp_z_order'th derivative of f1 (not g1)  
  subroutine add_hypz_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: p_f_
    complex, dimension(li1:li2,lj1:lj2),intent(inout) :: localrhs
    integer,intent(in):: k, l, m, n
    
    integer :: sten,i,j, pii, i1, i2, lij0_vwadp

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    lij0_vwadp = li0_vwadp(l,m)*lj0
    
#ifdef oldparbc
    if (shifted_metric) then
       pii=pi1
       do sten=-hyp_z_order/2,hyp_z_order/2
          do j=lj1,lj2
             do i=i1,i2
                if (pmi0.gt.1) pii=i
                localrhs(i,j)=localrhs(i,j)+&
                     & hyp_z_scl*hyp_z_sten(sten)*geom%phasefac(pii,j,k,sten)*p_f_(i,j,k+sten,l,m,n)
             enddo
          enddo
       enddo
    else
       Do sten=-hyp_z_order/2,hyp_z_order/2
          call axpy_ij(lij0_vwadp,hyp_z_scl*hyp_z_sten(sten), &
               & p_f_(i1:i2,:,k+sten,l,m,n),localrhs(i1:i2,:))
       End Do
    endif
#else
    if(.not.allocated(zbound1)) call ini_zbound_adptv
    do j=lj1,lj2
       do i=i1,i2
          if (zbound1(i,j,k,l)) then
          elseif (zbound2(i,j,k,l)) then
             do sten=-1,1
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc2coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          elseif (zbound3(i,j,k,l)) then
             do sten=-1,1
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc2coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          elseif (zbound4(i,j,k,l)) then
          else
             do sten=-2,2
                localrhs(i,j)=localrhs(i,j)+hyp_z_scl*dc4coeff(sten)*p_f_(i,j,k+sten,l,m,n)
             end do
          endif
       enddo
    enddo
    if(allocated(zbound1)) deallocate(zbound1,zbound2,zbound3,zbound4)
#endif    
    
  end subroutine add_hypz_indices_adptv

  subroutine add_dfdv_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n
    integer :: sten, kb, lb, ub, i1, i2, j, lij0_vwadp

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    lij0_vwadp = li0_vwadp(l,m)*lj0
    
    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)

    !finite differences with boundaries set to 0 (Dirichlet b.c.)
    if (xy_local) then
       do sten=lb,ub
          call axpy_ij(lij0_vwadp,trp_4(pi1,pj1,k,m,n)*vderivative(sten,pi1,k,l),&
               &p_f_(i1:i2,:,k,l+sten,m,n),localrhs(i1:i2,:))
       enddo
    else
       do sten=lb,ub
          do j=lj1,lj2
             localrhs(i1:i2,j) = localrhs(i1:i2,j) + trp_4(i1:i2,pj1,k,m,n) * &
                  &vderivative(sten,i1:i2,k,l)*p_f_(i1:i2,j,k,l+sten,m,n)
          enddo
       enddo
    endif
  end subroutine add_dfdv_indices_adptv

  subroutine add_hypv_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n

    integer :: sten, kb, lb, ub, i1, i2, lij0_vwadp

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    lij0_vwadp = li0_vwadp(l,m)*lj0

    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)
    do sten=lb,ub
       call axpy_ij(lij0_vwadp,vdamping(sten,l), &
            & p_f_(i1:i2,:,k,l+sten,m,n),localrhs(i1:i2,:))
    enddo
    
  end subroutine add_hypv_indices_adptv

  subroutine ini_zbound_adptv
    integer:: i,j,k,l, i1,i2
    real :: pdf1dz_m1n1
    allocate(zbound1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),zbound2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),&
         zbound3(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2),zbound4(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2))
    zbound1=.false.
    zbound2=.false.
    zbound3=.false.
    zbound4=.false.

    do l=ll1,ll2
       do k=lk1,lk2
          pdf1dz_m1n1 = - C_xy(pi1)*vTvpar(l,ln1)/&
               &(geom%jacobian(pi1,pj1,k)*geom%Bfield(pi1,pj1,k))
          do j=lj1,lj2
             i1 = li1_vadp(l)
             i2 = li2_vadp(l)
             do i= i1,i2              
                if ((k.eq.lk1).and.(my_pez.eq.0).and.delete_entry_lower(i,j)) then
                   if(pdf1dz_m1n1.gt.0) then
                      zbound1(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.(lk1+1)).and.(my_pez.eq.0).and.delete_entry_lower(i,j)) then             
                   if(pdf1dz_m1n1.gt.0) then
                      zbound2(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.(lk2-1)).and.(my_pez.eq.(n_procs_z-1)).and.delete_entry_upper(i,j)) then
                   if(pdf1dz_m1n1.lt.0) then
                      zbound3(i,j,k,l)=.true.
                   endif
                elseif ((k.eq.lk2).and.(my_pez.eq.(n_procs_z-1)).and.delete_entry_upper(i,j)) then
                   if(pdf1dz_m1n1.lt.0) then
                      zbound4(i,j,k,l)=.true.
                   endif
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine ini_zbound_adptv

end module dfdzv_terms
