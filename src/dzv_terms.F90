#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the parallel advection and trapping terms 
!\todo add a criterion for maxclimit 3 and 4 
module dzv_terms
  use par_mod
  use discretization
  use discretization_adptv_module ! adaptivity addons to discretization
  use blockindex
  use prefactors
  use boundaries
  use boundary_exchange_z, only: exchange_z_3d_equil
  use geometry
  use vel_space
  use axpy
  use numerical_damping
  use lagrange_interpolation
  use fourier
  use gyro_average_ff_mod, only: jfac, I1_factor
  use equilibrium_fields, only: dens_co, with_comoving_other

  implicit none
  public:: mem_est_dzv, initialize_dzv, equ_dzv, finalize_dzv, add_hypv_ak, add_hypz_ak
  public:: arakawa_cons_bc, arakawa_zv_order, equ_comp_hypz, set_dzv_defaults
  public:: add_Erad_acc, initialize_Erad_acc, finalize_Erad_acc 
  ! adaptivity functions
  public:: equ_dzv_adptv, add_hypv_ak_adptv, add_hypz_ak_adptv, equ_comp_hypz_adptv

  ! make some subroutines and variables public for unit tests
#ifdef GIUTESTS
  public:: df_pref_r, df_pref_c
  public:: equ_dzv_xydep_c_1, equ_dzv_xydep_c_1_adptv
  public:: equ_dzv_xydep_c_2, equ_dzv_xydep_c_2_adptv
  public:: equ_dzv_xdep_r_1, equ_dzv_xdep_r_1_adptv
  public:: equ_dzv_xdep_r_2, equ_dzv_xdep_r_2_adptv
  public:: replace_rhs
#endif
  
  private

  logical:: arakawa_cons_bc
  integer :: init_status = 0, arakawa_zv_order

  logical:: replace_rhs, dissipative_bc=.false.
  real, dimension(:,:,:,:,:,:), allocatable:: pdzv,poisson_h
  real,dimension(:,:,:,:),allocatable:: df_pref_r
  complex,dimension(:,:,:,:),allocatable:: df_pref_c
  integer,dimension(:),allocatable:: map_to_f, f_shifts
  integer,dimension(:,:),allocatable:: stenbound, shifted_f_pos
  real, dimension(:),allocatable:: vplocal
  real, dimension(:,:,:),allocatable:: Blocal
  integer:: j, k, l, m, n, pni
  integer:: nsten
  real, dimension(:,:,:,:), allocatable ::fac_co
  real,dimension(:,:,:,:,:,:),allocatable:: hv_pref
  real,dimension(:,:,:,:,:,:),allocatable:: hz_pref
  real,dimension(:,:,:,:,:,:,:),allocatable:: hz_comp_pref
  complex,dimension(:,:,:,:,:,:,:),allocatable:: hz_comp_pref_c
  real:: signchange

  integer, dimension(:), allocatable::ll_bound, ul_bound 
  real, dimension(:,:,:,:), allocatable:: Erad_array 

contains

  function mem_est_dzv(mem_req_in)
    real:: mem_req_in, mem_est_dzv, mem_loc
    real:: me
    !pdzv
    mem_loc=SIZE_OF_REAL_MB*pmi0*pj0*lklmn0
    !poisson_h
    mem_loc=mem_loc+SIZE_OF_REAL_MB*pmi0*pj0*lz0*lv0*lm0*ln0
    
    if (hyp_on_h) then
       !use simple criterion for hypz_compensation: kymin*rho_s<0.05
       me=1
       do n=0,n_spec-1
          if (spec(n)%charge==-1) me=spec(n)%mass
       enddo
       if (kymin.lt.3./sqrt(me)) then
          if (shifted_metric) then
             mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*lj0*5*lk0*ll0*lm0*ln0
          else
             mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*5*lk0*ll0*lm0*ln0
          endif
       endif
    endif

    select case(perf_vec(3))
    case(1)
       mem_loc=mem_est_dzv_1()
    case(2)
       mem_loc=mem_est_dzv_2()
    end select
    
    mem_est_dzv=mem_req_in+mem_loc
    
  end function mem_est_dzv


  subroutine set_dzv_defaults
    arakawa_cons_bc = .false.
    arakawa_zv_order = 4
    dissipative_bc = .false.
  end subroutine set_dzv_defaults


  subroutine initialize_dzv(dzv_replace_rhs)
    logical:: dzv_replace_rhs
    integer:: i 
  
    replace_rhs=dzv_replace_rhs

    if(xy_local.and.(.not.arakawa_cons_bc).and.arakawa_zv_order.eq.2) then
       dissipative_bc=.true.
    endif

    if(init_status.eq.0) then    !initialize from scratch
       if (with_comoving_other) allocate(fac_co(pi1:pi2,pj1:pj2,lbz:ubz, ln1:ln2))
       allocate(vplocal(ll1-2:ll2+2), Blocal(pi1:pi2,pj1:pj2,lbz:ubz))
       Blocal(pi1:pi2,pj1:pj2,lk1:lk2)=geom%Bfield(pi1:pi2,pj1:pj2,lk1:lk2)
       call exchange_z_3d_equil(Blocal,n0_global*q0)

       do l=ll1-2,ll2+2
#ifdef COMBI
          vplocal(l)=-lv+l*dv + SHIFT
#else
          vplocal(l)=-lv+l*dv
#endif
       enddo
       
       allocate(pdzv(pi1:pi2, pj1:pj2, lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
       if (.not.x_local.and..not.lilo) then
          pni = pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=pj1,pj2
                         pdzv(:,j,k,l,m,n) = C_xy(pi1:pi2)/sqrt(2.*spec(n)%mass/spec(n)%temp)&
                              &/(geom%jacobian(pi1:pi2,j,k)*geom%Bfield(pi1:pi2,j,k))
                      end do
                   end do
                end do
             end do
          end do
       else
          pni = pn1
          do n=ln1,ln2
             if (pn0.gt.1) pni=n
             do m=lm1,lm2
                do l=ll1,ll2
                   do k=lk1,lk2
                      do j=pj1,pj2
                         pdzv(:,j,k,l,m,n) = C_xy(pi1:pi2)/sqrt(2.*spec(n)%mass/spec(n)%temp)&
                              &/(geom%jacobian(pi1:pi2,j,k)*geom%Bfield(pi1:pi2,j,k))*fm(:,j,k,l,m,pni)*spec(n)%temp_prof(pi1:pi2)
                      end do
                   end do
                end do
             end do
          end do
       endif
       allocate(poisson_h(pi1:pi2, pj1:pj2, lbz:ubz, ll1-2:ll2+2, lm1:lm2, ln1:ln2))
       poisson_h=0.      
       if (.not.x_local.and..not.lilo) then
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1-2,ll2+2
                   do k=lbz,ubz
                      do j=pj1,pj2
                         poisson_h(:,j,k,l,m,n)=mu(m)*Blocal(pi1:pi2,j,k)+vplocal(l)**2
                      enddo
                   enddo
                enddo
             enddo
          enddo
       else
          if (with_comoving_other) then
             do n=ln1,ln2
                do k=lk1,lk2
                   do j=pj1,pj2
                      do i=pi1,pi2
                         fac_co(i,j,k,n)=dens_co(i,k,n)
                      enddo
                   enddo
                enddo

                call exchange_z_3d_equil(fac_co(:,:,:,n),n0_global*q0)

             enddo
          endif
          do n=ln1,ln2
             do m=lm1,lm2
                do l=ll1-2,ll2+2
                   do k=lbz,ubz
                      do j=pj1,pj2
                         poisson_h(:,j,k,l,m,n)=1./spec(n)%dens_prof(pi1:pi2)*(pi*spec(n)%temp_prof(pi1:pi2))**(1.5) *&
                              & exp((mu(m)*Blocal(pi1:pi2,j,k)+vplocal(l)**2)/spec(n)%temp_prof(pi1:pi2))
                         if (with_comoving_other) poisson_h(:,j,k,l,m,n)=poisson_h(:,j,k,l,m,n)/fac_co(:,j,k,n)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       endif
       signchange=0.5*(nv0-1)

       if (hyp_on_h) then
          call initialize_hv_pref
          call initialize_hz_pref
          if (hypz_compensation) call initialize_hypz_compensation
       endif
       !perf initialization
       select case(perf_vec(3))
       case(1)
          call initialize_dzv_1
       case(2)
          call initialize_dzv_2
       end select
       init_status = perf_vec(3)
    elseif (init_status.ne.perf_vec(3)) then
       !only delete old perf initialization
       select case(init_status)
       case(1)
          call finalize_dzv_1
       case(2)
          call finalize_dzv_2
       end select
       !initialization new perf
       select case(perf_vec(3))
       case(1)
          call initialize_dzv_1
       case(2)
          call initialize_dzv_2
       end select
       init_status = perf_vec(3)
    endif

    if (Erad_acc.and.Erad.ne.0..and..not.y_local) then 
       call initialize_Erad_acc 
    endif

  end subroutine initialize_dzv
  
  subroutine equ_dzv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs

    if(xy_local) then
       if (dissipative_bc) then
          if (perf_vec(3).eq.1) then
             call equ_dzv_xydep_r_1(p_f,p_rhs,lb1,lb2,df_pref_r)
          else
             call equ_dzv_xydep_r_2(p_f,p_rhs,lb1,lb2,df_pref_r)
          end if
       else
          if (perf_vec(3).eq.1) then
             call equ_dzv_r_1(p_f,p_rhs,lb1,lb2,df_pref_r)
          else
             call equ_dzv_r_2(p_f,p_rhs,lb1,lb2,df_pref_r)
          end if
       end if
    else
       if (shifted_metric) then
          if (perf_vec(3).eq.1) then
             call equ_dzv_xydep_c_1(p_f,p_rhs,lb1,lb2,df_pref_c)
          else
             call equ_dzv_xydep_c_2(p_f,p_rhs,lb1,lb2,df_pref_c)
          end if
       else
          if (perf_vec(3).eq.1) then
             call equ_dzv_xdep_r_1(p_f,p_rhs,lb1,lb2,df_pref_r)
          else
             call equ_dzv_xdep_r_2(p_f,p_rhs,lb1,lb2,df_pref_r)
          end if
       end if
    end if

    if (Erad_acc.and.Erad.ne.0..and..not.y_local) then 
       call add_Erad_acc(p_f, lb1, lb2, p_rhs)
    end if
  end subroutine equ_dzv

  subroutine finalize_dzv
    
    select case(init_status)
    case(1)
       call finalize_dzv_1
    case(2)
       call finalize_dzv_2
    end select

    if (hyp_on_h) then
       deallocate(hv_pref, hz_pref)
       if (hypz_compensation) then
          if (shifted_metric) then
             deallocate(hz_comp_pref_c)
          else
             deallocate(hz_comp_pref)
          endif
       endif
    endif

    deallocate(vplocal, Blocal)
    deallocate(pdzv,poisson_h)
    if (with_comoving_other) deallocate(fac_co)
    init_status = 0
    if (Erad_acc.and.Erad.ne.0..and..not.y_local) then 
       call finalize_Erad_acc 
    end if 
  end subroutine finalize_dzv


  !------------------------------------------
!***********************************************************************************************
  
  function mem_est_dzv_1()
    real:: mem_est_dzv_1
    
    !df_pref
    if (shifted_metric) then
       mem_est_dzv_1=13*pmi0*lj0*lklmn0*SIZE_OF_COMPLEX_MB
    elseif(xy_local.and.(.not.arakawa_cons_bc)) then
       mem_est_dzv_1=13*li0*lj0*lklmn0*SIZE_OF_REAL_MB
    else
       mem_est_dzv_1=13*pmi0*pj0*lklmn0*SIZE_OF_REAL_MB
    endif

    mem_est_dzv_1=mem_est_dzv_1+3*lklmn0*SIZE_OF_INTEGER_MB
    
  end function mem_est_dzv_1

  subroutine initialize_dzv_1
    integer:: i,j,k,l,m,n,sten,klmn, pii, pji
    integer:: ilb, iub, jlb, jub

    if(dissipative_bc.or.shifted_metric) then
       ilb=li1
       iub=li2
       jlb=lj1
       jub=lj2
    else
       ilb=pi1
       iub=pi2
       jlb=pj1
       jub=pj2
    end if

    pii=li1
    pji=lj1

    if (shifted_metric) then
       allocate(df_pref_c(ilb:iub,jlb:jub,-6:6,lklmn0))
       df_pref_c=0.
    else
       allocate(df_pref_r(ilb:iub,jlb:jub,-6:6,lklmn0))
       df_pref_r=0.
    endif

    allocate(map_to_f(lklmn0),stenbound(2,lklmn0),f_shifts(-6:6))
 
    f_shifts(-6)=-2*lz0
    f_shifts(-5)=-lz0-1
    f_shifts(-4)=-lz0
    f_shifts(-3)=-lz0+1
    do sten=-2,2
       f_shifts(sten)=sten
    enddo
    f_shifts(3)=lz0-1
    f_shifts(4)=lz0
    f_shifts(5)=lz0+1
    f_shifts(6)=2*lz0

    stenbound(1,:)=-6
    stenbound(2,:)=6


    if(shifted_metric) then
       !prefactor is complex
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2 
                do k=lk1,lk2
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   if((l.eq.ll1).and.(my_pev.eq.0)) stenbound(1,klmn)=-2
                   if((l.eq.ll1+1).and.(my_pev.eq.0)) stenbound(1,klmn)=-5
                   if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=5
                   if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=2
                   do j=jlb,jub
                      if(pj0.gt.1) pji=j
                      do i=ilb,iub
                         if(pmi0.gt.1) pii=i
                         df_pref_c(i,j,-6,klmn)=0.
                         
                         df_pref_c(i,j,-5,klmn)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                         
                         df_pref_c(i,j,-4,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                              +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order
                         
                         df_pref_c(i,j,-3,klmn)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,1)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n))  !second order
                         
                         df_pref_c(i,j,-2,klmn)=0.
                         
                         df_pref_c(i,j,-1,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                              +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                         
                         df_pref_c(i,j,0,klmn)=0.
                         
                         df_pref_c(i,j,1,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,1)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                              +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                         
                         df_pref_c(i,j,2,klmn)=0.
                         
                         df_pref_c(i,j,3,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n))  !second order
                         
                         df_pref_c(i,j,4,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)+&  !second order
                              &poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))  !second order
                         
                         df_pref_c(i,j,5,klmn)= -1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,1) &
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n))  !second order
                         
                         df_pref_c(i,j,6,klmn)=0.

                         if (arakawa_zv_order==4) then
                            df_pref_c(i,j,-6,klmn)=df_pref_c(i,j,-6,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_c(i,j,-5,klmn)=2*df_pref_c(i,j,-5,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (-poisson_h(pii,pji,k-2,l,m,n)+poisson_h(pii,pji,k,l-2,m,n)&
                                 -poisson_h(pii,pji,k-1,l+1,m,n)+poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,-4,klmn)=2*df_pref_c(i,j,-4,klmn)
                            df_pref_c(i,j,-3,klmn)=2*df_pref_c(i,j,-3,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (poisson_h(pii,pji,k+2,l,m,n)-poisson_h(pii,pji,k,l-2,m,n)&
                                 +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,-2,klmn)=df_pref_c(i,j,-2,klmn)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))&
                              *geom%phasefac(i,j,k,-2)
                            df_pref_c(i,j,-1,klmn)=2*df_pref_c(i,j,-1,klmn)
                            df_pref_c(i,j,1,klmn)=2*df_pref_c(i,j,1,klmn)
                            df_pref_c(i,j,2,klmn)=df_pref_c(i,j,2,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,2)
                            df_pref_c(i,j,3,klmn)=2*df_pref_c(i,j,3,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k-2,l,m,n)-poisson_h(pii,pji,k,l+2,m,n)&
                                 -poisson_h(pii,pji,k+1,l+1,m,n)+poisson_h(pii,pji,k-1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,4,klmn)=2*df_pref_c(i,j,4,klmn)
                            df_pref_c(i,j,5,klmn)=2*df_pref_c(i,j,5,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(-poisson_h(pii,pji,k+2,l,m,n)+poisson_h(pii,pji,k,l+2,m,n)&
                                 +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,6,klmn)=df_pref_c(i,j,6,klmn)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))
                         endif

                         
                         if (hyp_on_h) then
                            df_pref_c(i,j,-6,klmn)=df_pref_c(i,j,-6,klmn) + hv_pref(pii,pji,k,-2,l,n)
                            df_pref_c(i,j,-4,klmn)=df_pref_c(i,j,-4,klmn) + hv_pref(pii,pji,k,-1,l,n)
                            df_pref_c(i,j,-2,klmn)=df_pref_c(i,j,-2,klmn) + hz_pref(pii,pji,-2,k,m,n)*geom%phasefac(i,j,k,-2)
                            df_pref_c(i,j,-1,klmn)=df_pref_c(i,j,-1,klmn) + hz_pref(pii,pji,-1,k,m,n)*geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,0,klmn)=df_pref_c(i,j,0,klmn) + hz_pref(pii,pji,0,k,m,n) + hv_pref(pii,pji,k,0,l,n)
                            df_pref_c(i,j,1,klmn)=df_pref_c(i,j,1,klmn) + hz_pref(pii,pji,1,k,m,n)*geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,2,klmn)=df_pref_c(i,j,2,klmn) + hz_pref(pii,pji,2,k,m,n)*geom%phasefac(i,j,k,2)
                            df_pref_c(i,j,4,klmn)=df_pref_c(i,j,4,klmn) + hv_pref(pii,pji,k,1,l,n)
                            df_pref_c(i,j,6,klmn)=df_pref_c(i,j,6,klmn) + hv_pref(pii,pji,k,2,l,n)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo       
    else
       !prefactor is real
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2 
                do k=lk1,lk2
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   map_to_f(klmn)= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   if((l.eq.ll1).and.(my_pev.eq.0)) stenbound(1,klmn)=-2
                   if((l.eq.ll1+1).and.(my_pev.eq.0)) stenbound(1,klmn)=-5
                   if((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=5
                   if((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) stenbound(2,klmn)=2
                   do j=jlb,jub
                      if(pj0.gt.1) pji=j
                      do i=ilb,iub
                         if(pmi0.gt.1) pii=i
                         df_pref_r(i,j,-6,klmn)=0.
                         
                         df_pref_r(i,j,-5,klmn)=pdzv(pii,pji,k,l,m,n)/dv/dz*(-1/12.&
                              *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n))) !second order
                         
                         df_pref_r(i,j,-4,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                              +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order
                         
                         df_pref_r(i,j,-3,klmn)=pdzv(pii,pji,k,l,m,n)/dv/dz*(-1/12.&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n)))  !second order
                         
                         df_pref_r(i,j,-2,klmn)=0.
                         
                         df_pref_r(i,j,-1,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                              +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                         
                         df_pref_r(i,j,0,klmn)=0.
                         
                         df_pref_r(i,j,1,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                              +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                         
                         df_pref_r(i,j,2,klmn)=0.
                         
                         df_pref_r(i,j,3,klmn)=pdzv(pii,pji,k,l,m,n)/dv/dz*(-1/12.&
                              *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n)))  !second order

                         
                         df_pref_r(i,j,4,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)+&  !second order
                              &poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))  !second order
                         
                         df_pref_r(i,j,5,klmn)= pdzv(pii,pji,k,l,m,n)/dv/dz*(-1/12.&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n)))  !second order

                         
                         df_pref_r(i,j,6,klmn)=0.
                         
                         if (arakawa_zv_order==4) then
                            df_pref_r(i,j,-6,klmn)=df_pref_r(i,j,-6,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,-5,klmn)=2*df_pref_r(i,j,-5,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (-poisson_h(pii,pji,k-2,l,m,n)+poisson_h(pii,pji,k,l-2,m,n)&
                                 -poisson_h(pii,pji,k-1,l+1,m,n)+poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,-4,klmn)=2*df_pref_r(i,j,-4,klmn)
                            df_pref_r(i,j,-3,klmn)=2*df_pref_r(i,j,-3,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (poisson_h(pii,pji,k+2,l,m,n)-poisson_h(pii,pji,k,l-2,m,n)&
                                 +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,-2,klmn)=df_pref_r(i,j,-2,klmn)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,-1,klmn)=2*df_pref_r(i,j,-1,klmn)
                            df_pref_r(i,j,1,klmn)=2*df_pref_r(i,j,1,klmn)
                            df_pref_r(i,j,2,klmn)=df_pref_r(i,j,2,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,3,klmn)=2*df_pref_r(i,j,3,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k-2,l,m,n)-poisson_h(pii,pji,k,l+2,m,n)&
                                 -poisson_h(pii,pji,k+1,l+1,m,n)+poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,4,klmn)=2*df_pref_r(i,j,4,klmn)
                            df_pref_r(i,j,5,klmn)=2*df_pref_r(i,j,5,klmn)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(-poisson_h(pii,pji,k+2,l,m,n)+poisson_h(pii,pji,k,l+2,m,n)&
                                 +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,6,klmn)=df_pref_r(i,j,6,klmn)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))
                         endif
                         
                         if (dissipative_bc) then
                            if ((k.eq.0).and.delete_entry_lower(i,j).and.(l.lt.signchange)) then    
                               df_pref_r(i,j,-6,klmn)=0.
                               df_pref_r(i,j,-5,klmn)=0.
                               
                               df_pref_r(i,j,-4,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                                    +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)&
                                    +poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                               
                               df_pref_r(i,j,-3,klmn)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,-2,klmn)=0.
                               df_pref_r(i,j,-1,klmn)=0.
                               
                               df_pref_r(i,j,0,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                                    +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                               
                               df_pref_r(i,j,1,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,2,klmn)=0.
                               df_pref_r(i,j,3,klmn)=0.
                               
                               if(l.ne.floor(signchange)) then     
                                  df_pref_r(i,j,4,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                       *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&  !second order
                                       +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n)&
                                       -poisson_h(pii,pji,k-1,l,m,n)+poisson_h(pii,pji,k,l+1,m,n))  !second order
                               end if
                               
                               df_pref_r(i,j,5,klmn)= -1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n))  !second order
                               
                               df_pref_r(i,j,6,klmn)=0.
                            elseif ((k.eq.nz0-1).and.delete_entry_upper(i,j).and.(l.gt.signchange)) then
                               df_pref_r(i,j,-6,klmn)=0.
                               
                               df_pref_r(i,j,-5,klmn)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                               
                               if(l.ne.ceiling(signchange)) then     
                                  df_pref_r(i,j,-4,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                       *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                                       +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)&
                                       +poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n)) !second order
                               end if
                               
                               df_pref_r(i,j,-3,klmn)=0.
                               df_pref_r(i,j,-2,klmn)=0.
                               
                               df_pref_r(i,j,-1,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                                    +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                               
                               df_pref_r(i,j,0,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,1,klmn)=0.
                               df_pref_r(i,j,2,klmn)=0.
                               
                               df_pref_r(i,j,3,klmn)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n))  !second order
                               
                               df_pref_r(i,j,4,klmn)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n)&
                                    -poisson_h(pii,pji,k,l+1,m,n)+poisson_h(pii,pji,k+1,l,m,n))  !second order
                               
                               df_pref_r(i,j,5,klmn)=0.                            
                               df_pref_r(i,j,6,klmn)=0.
                            end if
                         end if
                         
                         if (hyp_on_h) then
                            df_pref_r(i,j,-6,klmn)=df_pref_r(i,j,-6,klmn) + hv_pref(pii,pji,k,-2,l,n)
                            df_pref_r(i,j,-4,klmn)=df_pref_r(i,j,-4,klmn) + hv_pref(pii,pji,k,-1,l,n)
                            df_pref_r(i,j,-2,klmn)=df_pref_r(i,j,-2,klmn) + hz_pref(pii,pji,-2,k,m,n)
                            df_pref_r(i,j,-1,klmn)=df_pref_r(i,j,-1,klmn) + hz_pref(pii,pji,-1,k,m,n)
                            df_pref_r(i,j,0,klmn)=df_pref_r(i,j,0,klmn) + hz_pref(pii,pji,0,k,m,n) + hv_pref(pii,pji,k,0,l,n)
                            df_pref_r(i,j,1,klmn)=df_pref_r(i,j,1,klmn) + hz_pref(pii,pji,1,k,m,n)
                            df_pref_r(i,j,2,klmn)=df_pref_r(i,j,2,klmn) + hz_pref(pii,pji,2,k,m,n)
                            df_pref_r(i,j,4,klmn)=df_pref_r(i,j,4,klmn) + hv_pref(pii,pji,k,1,l,n)
                            df_pref_r(i,j,6,klmn)=df_pref_r(i,j,6,klmn) + hv_pref(pii,pji,k,2,l,n)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    end if
  end subroutine initialize_dzv_1
  
  subroutine equ_dzv_r_1(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(1:lij0,lzvwn0), intent(in):: p_f
    complex, dimension(1:lij0, lb1:lb2), intent(inout):: p_rhs
    real, dimension(-6:6,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind

    PERFON_I('dzv1')
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
    
  end subroutine equ_dzv_r_1

  subroutine equ_dzv_xdep_r_1(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-6:6,lklmn0), intent(in) ::  prefac
    integer:: j, klmn, sten, f_ind

    PERFON_I('dzv1')
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
    
  end subroutine equ_dzv_xdep_r_1

  subroutine equ_dzv_xydep_r_1(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lj1:lj2,-6:6,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind

    PERFON_I('dzv1')
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
    
  end subroutine equ_dzv_xydep_r_1

  subroutine equ_dzv_xydep_c_1(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-6:6,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind

    PERFON_I('dzv1')
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
    
  end subroutine equ_dzv_xydep_c_1
  
  subroutine finalize_dzv_1
    if (shifted_metric) then
       deallocate(df_pref_c)
    else
       deallocate(df_pref_r)
    endif
    deallocate(map_to_f,stenbound,f_shifts)
  end subroutine finalize_dzv_1

!***********************************************************************************************


  function mem_est_dzv_2()
    real:: mem_est_dzv_2
    
    !df_pref
    if (shifted_metric) then
       mem_est_dzv_2=13*pmi0*lj0*lklmn0*SIZE_OF_COMPLEX_MB
    elseif(xy_local.and.(.not.arakawa_cons_bc)) then
       mem_est_dzv_2=13*li0*lj0*lklmn0*SIZE_OF_REAL_MB
    else
       mem_est_dzv_2=13*pmi0*pj0*lklmn0*SIZE_OF_REAL_MB
    endif

    !shifted_f_pos
    mem_est_dzv_2=mem_est_dzv_2+&
         13*lklmn0*SIZE_OF_INTEGER_MB
    
  end function mem_est_dzv_2


  subroutine initialize_dzv_2
    integer,dimension(-6:6):: f_shift
    integer:: i, j, k, l, m, n, ilb, iub, jlb, jub, pii, pji
    integer:: sten, klmn, f_index

    nsten=6

    if(dissipative_bc.or.shifted_metric) then
       ilb=li1
       iub=li2
       jlb=lj1
       jub=lj2
    else
       ilb=pi1
       iub=pi2
       jlb=pj1
       jub=pj2
    end if

    pii=li1
    pji=lj1

    if (shifted_metric) then
       allocate(df_pref_c(ilb:iub,jlb:jub,lklmn0,-6:6))
       df_pref_c=0.
    else
       allocate(df_pref_r(ilb:iub,jlb:jub,lklmn0,-6:6))
       df_pref_r=0.
    endif

    allocate(shifted_f_pos(lklmn0,-6:6))
    shifted_f_pos=0

    f_shift(-6)=-2*lz0
    f_shift(-5)=-lz0-1
    f_shift(-4)=-lz0
    f_shift(-3)=-lz0+1
    do sten=-2,2
       f_shift(sten)=sten
    enddo
    f_shift(3)=lz0-1
    f_shift(4)=lz0
    f_shift(5)=lz0+1
    f_shift(6)=2*lz0

    if (shifted_metric) then
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2 
                do k=lk1,lk2 
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   f_index= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   do sten=-6,6
                      shifted_f_pos(klmn,sten)=f_index+f_shift(sten)
                      if(shifted_f_pos(klmn,sten).lt.1) then
                         shifted_f_pos(klmn,sten)=1
                      end if
                      if(shifted_f_pos(klmn,sten).gt.lzvwn0) then
                         shifted_f_pos(klmn,sten)=lzvwn0
                      endif
                   end do                  
                   do j=jlb,jub
                      do i=ilb,iub
                         if(pmi0.gt.1) pii=i
                         if(pj0.gt.1) pji=j
                         df_pref_c(i,j,klmn,-6)=0.
                         
                         df_pref_c(i,j,klmn,-5)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                         
                         df_pref_c(i,j,klmn,-4)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                              +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order
                         
                         df_pref_c(i,j,klmn,-3)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,1)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n))  !second order
                         
                         df_pref_c(i,j,klmn,-2)=0.
                         
                         df_pref_c(i,j,klmn,-1)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                              +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                         
                         df_pref_c(i,j,klmn,0)=0.
                         
                         df_pref_c(i,j,klmn,1)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,1)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                              +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                         
                         df_pref_c(i,j,klmn,2)=0.
                         
                         df_pref_c(i,j,klmn,3)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)*geom%phasefac(i,j,k,-1)&
                              *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n))  !second order
                         
                         df_pref_c(i,j,klmn,4)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)+&  !second order
                              &poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))  !second order
                         
                         df_pref_c(i,j,klmn,5)= -1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz*geom%phasefac(i,j,k,1) &
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n))  !second order
                         
                         df_pref_c(i,j,klmn,6)=0.

                         if (arakawa_zv_order==4) then
                            df_pref_c(i,j,klmn,-6)=df_pref_c(i,j,klmn,-6)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_c(i,j,klmn,-5)=2*df_pref_c(i,j,klmn,-5)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (-poisson_h(pii,pji,k-2,l,m,n)+poisson_h(pii,pji,k,l-2,m,n)&
                                 -poisson_h(pii,pji,k-1,l+1,m,n)+poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,klmn,-4)=2*df_pref_c(i,j,klmn,-4)
                            df_pref_c(i,j,klmn,-3)=2*df_pref_c(i,j,klmn,-3)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (poisson_h(pii,pji,k+2,l,m,n)-poisson_h(pii,pji,k,l-2,m,n)&
                                 +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,klmn,-2)=df_pref_c(i,j,klmn,-2)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))&
                              *geom%phasefac(i,j,k,-2)
                            df_pref_c(i,j,klmn,-1)=2*df_pref_c(i,j,klmn,-1)
                            df_pref_c(i,j,klmn,1)=2*df_pref_c(i,j,klmn,1)
                            df_pref_c(i,j,klmn,2)=df_pref_c(i,j,klmn,2)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,2)
                            df_pref_c(i,j,klmn,3)=2*df_pref_c(i,j,klmn,3)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k-2,l,m,n)-poisson_h(pii,pji,k,l+2,m,n)&
                                 -poisson_h(pii,pji,k+1,l+1,m,n)+poisson_h(pii,pji,k-1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,klmn,4)=2*df_pref_c(i,j,klmn,4)
                            df_pref_c(i,j,klmn,5)=2*df_pref_c(i,j,klmn,5)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(-poisson_h(pii,pji,k+2,l,m,n)+poisson_h(pii,pji,k,l+2,m,n)&
                                 +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))&
                                 *geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,klmn,6)=df_pref_c(i,j,klmn,6)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))
                         endif

                         if (hyp_on_h) then
                            df_pref_c(i,j,klmn,-6)=df_pref_c(i,j,klmn,-6) + hv_pref(pii,pji,k,-2,l,n)
                            df_pref_c(i,j,klmn,-4)=df_pref_c(i,j,klmn,-4) + hv_pref(pii,pji,k,-1,l,n)
                            df_pref_c(i,j,klmn,-2)=df_pref_c(i,j,klmn,-2) + hz_pref(pii,pji,-2,k,m,n)*geom%phasefac(i,j,k,-2)
                            df_pref_c(i,j,klmn,-1)=df_pref_c(i,j,klmn,-1) + hz_pref(pii,pji,-1,k,m,n)*geom%phasefac(i,j,k,-1)
                            df_pref_c(i,j,klmn,0)=df_pref_c(i,j,klmn,0) + hz_pref(pii,pji,0,k,m,n) + hv_pref(pii,pji,k,0,l,n)
                            df_pref_c(i,j,klmn,1)=df_pref_c(i,j,klmn,1) + hz_pref(pii,pji,1,k,m,n)*geom%phasefac(i,j,k,1)
                            df_pref_c(i,j,klmn,2)=df_pref_c(i,j,klmn,2) + hz_pref(pii,pji,2,k,m,n)*geom%phasefac(i,j,k,2)
                            df_pref_c(i,j,klmn,4)=df_pref_c(i,j,klmn,4) + hv_pref(pii,pji,k,1,l,n)
                            df_pref_c(i,j,klmn,6)=df_pref_c(i,j,klmn,6) + hv_pref(pii,pji,k,2,l,n)
                         endif                         
                         !boundary conditions in v_parallel
                         if((l.eq.0)) df_pref_c(i,j,klmn,-6:-3)=0.
                         if((l.eq.1)) df_pref_c(i,j,klmn,-6)=0.
                         if((l.eq.(nv0-2))) df_pref_c(i,j,klmn,6)=0.
                         if((l.eq.(nv0-1))) df_pref_c(i,j,klmn,3:6)=0.
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    else
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2 
                do k=lk1,lk2 
                   klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                   f_index= (n-ln1)*lzvw0 + (nwb+m-lm1)*lz0*lv0 + (nvb+l-ll1)*lz0 + (nzb+k-lk1)+1
                   do sten=-6,6
                      shifted_f_pos(klmn,sten)=f_index+f_shift(sten)
                      if(shifted_f_pos(klmn,sten).lt.1) then
                         shifted_f_pos(klmn,sten)=1
                      end if
                      if(shifted_f_pos(klmn,sten).gt.lzvwn0) then
                         shifted_f_pos(klmn,sten)=lzvwn0
                      endif
                   end do
                 do j=jlb,jub
                      do i=ilb,iub
                         if(pmi0.gt.1) pii=i
                         if(pj0.gt.1) pji=j
                         
                         df_pref_r(i,j,klmn,-6)=0.
                         
                         df_pref_r(i,j,klmn,-5)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                              *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                         
                         df_pref_r(i,j,klmn,-4)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                              +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order
                         
                         df_pref_r(i,j,klmn,-3)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n))  !second order
                         
                         df_pref_r(i,j,klmn,-2)=0.
                         
                         df_pref_r(i,j,klmn,-1)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                              +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                         
                         df_pref_r(i,j,klmn,0)=0.
                         
                         df_pref_r(i,j,klmn,1)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                              +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                         
                         df_pref_r(i,j,klmn,2)=0.
                         
                         df_pref_r(i,j,klmn,3)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n))  !second order
                         
                         df_pref_r(i,j,klmn,4)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)+&  !second order
                              &poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))  !second order
                         
                         df_pref_r(i,j,klmn,5)= -1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                              *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n))  !second order
                         
                         df_pref_r(i,j,klmn,6)=0.
                         
                         if (arakawa_zv_order==4) then
                            df_pref_r(i,j,klmn,-6)=df_pref_r(i,j,klmn,-6)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,klmn,-5)=2*df_pref_r(i,j,klmn,-5)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (-poisson_h(pii,pji,k-2,l,m,n)+poisson_h(pii,pji,k,l-2,m,n)&
                                 -poisson_h(pii,pji,k-1,l+1,m,n)+poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,klmn,-4)=2*df_pref_r(i,j,klmn,-4)
                            df_pref_r(i,j,klmn,-3)=2*df_pref_r(i,j,klmn,-3)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)*&
                                 (poisson_h(pii,pji,k+2,l,m,n)-poisson_h(pii,pji,k,l-2,m,n)&
                                 +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,klmn,-2)=df_pref_r(i,j,klmn,-2)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                              *(poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,klmn,-1)=2*df_pref_r(i,j,klmn,-1)
                            df_pref_r(i,j,klmn,1)=2*df_pref_r(i,j,klmn,1)
                            df_pref_r(i,j,klmn,2)=df_pref_r(i,j,klmn,2)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,klmn,3)=2*df_pref_r(i,j,klmn,3)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k-2,l,m,n)-poisson_h(pii,pji,k,l+2,m,n)&
                                 -poisson_h(pii,pji,k+1,l+1,m,n)+poisson_h(pii,pji,k-1,l-1,m,n))
                            df_pref_r(i,j,klmn,4)=2*df_pref_r(i,j,klmn,4)
                            df_pref_r(i,j,klmn,5)=2*df_pref_r(i,j,klmn,5)+1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(-poisson_h(pii,pji,k+2,l,m,n)+poisson_h(pii,pji,k,l+2,m,n)&
                                 +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))
                            df_pref_r(i,j,klmn,6)=df_pref_r(i,j,klmn,6)-1/24./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                 *(poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n))
                         endif

                         
                         if (dissipative_bc) then
                            if ((k.eq.0).and.delete_entry_lower(i,j).and.(l.lt.signchange)) then    
                               df_pref_r(i,j,klmn,-6)=0.
                               df_pref_r(i,j,klmn,-5)=0.
                               
                               df_pref_r(i,j,klmn,-4)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                                    +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)&
                                    +poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                               
                               df_pref_r(i,j,klmn,-3)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,-2)=0.
                               df_pref_r(i,j,klmn,-1)=0.
                               
                               df_pref_r(i,j,klmn,0)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                                    +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                               
                               df_pref_r(i,j,klmn,1)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,2)=0.
                               df_pref_r(i,j,klmn,3)=0.
                               
                               if(l.ne.floor(signchange)) then     
                                  df_pref_r(i,j,klmn,4)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                       *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&  !second order
                                       +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n)&
                                       -poisson_h(pii,pji,k-1,l,m,n)+poisson_h(pii,pji,k,l+1,m,n))  !second order
                               end if
                               
                               df_pref_r(i,j,klmn,5)= -1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k+1,l,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,6)=0.
                            elseif ((k.eq.nz0-1).and.delete_entry_upper(i,j).and.(l.gt.signchange)) then
                               df_pref_r(i,j,klmn,-6)=0.
                               
                               df_pref_r(i,j,klmn,-5)=-1/12.*pdzv(pii,pji,k,l,m,n)/dv/dz&
                                    *(poisson_h(pii,pji,k,l-1,m,n)-poisson_h(pii,pji,k-1,l,m,n)) !second order
                               
                               if(l.ne.ceiling(signchange)) then     
                                  df_pref_r(i,j,klmn,-4)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                       *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&
                                       +poisson_h(pii,pji,k+1,l-1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)&
                                       +poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k,l-1,m,n)) !second order
                               end if
                               
                               df_pref_r(i,j,klmn,-3)=0.
                               df_pref_r(i,j,klmn,-2)=0.
                               
                               df_pref_r(i,j,klmn,-1)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)& !second order
                                    +poisson_h(pii,pji,k-1,l+1,m,n)-poisson_h(pii,pji,k-1,l-1,m,n)) !second order 
                               
                               df_pref_r(i,j,klmn,0)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k,l+1,m,n)-poisson_h(pii,pji,k,l-1,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k+1,l-1,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,1)=0.
                               df_pref_r(i,j,klmn,2)=0.
                               
                               df_pref_r(i,j,klmn,3)=-1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k-1,l,m,n)-poisson_h(pii,pji,k,l+1,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,4)=1/12./dv/dz*pdzv(pii,pji,k,l,m,n)&
                                    *(poisson_h(pii,pji,k+1,l,m,n)-poisson_h(pii,pji,k-1,l,m,n)&  !second order
                                    +poisson_h(pii,pji,k+1,l+1,m,n)-poisson_h(pii,pji,k-1,l+1,m,n)&
                                    -poisson_h(pii,pji,k,l+1,m,n)+poisson_h(pii,pji,k+1,l,m,n))  !second order
                               
                               df_pref_r(i,j,klmn,5)=0.                            
                               df_pref_r(i,j,klmn,6)=0.
                            end if
                         end if
                         
                         if (hyp_on_h) then
                            df_pref_r(i,j,klmn,-6)=df_pref_r(i,j,klmn,-6) + hv_pref(pii,pji,k,-2,l,n)
                            df_pref_r(i,j,klmn,-4)=df_pref_r(i,j,klmn,-4) + hv_pref(pii,pji,k,-1,l,n)
                            df_pref_r(i,j,klmn,-2)=df_pref_r(i,j,klmn,-2) + hz_pref(pii,pji,-2,k,m,n)
                            df_pref_r(i,j,klmn,-1)=df_pref_r(i,j,klmn,-1) + hz_pref(pii,pji,-1,k,m,n)
                            df_pref_r(i,j,klmn,0)=df_pref_r(i,j,klmn,0) + hz_pref(pii,pji,0,k,m,n) + hv_pref(pii,pji,k,0,l,n)
                            df_pref_r(i,j,klmn,1)=df_pref_r(i,j,klmn,1) + hz_pref(pii,pji,1,k,m,n)
                            df_pref_r(i,j,klmn,2)=df_pref_r(i,j,klmn,2) + hz_pref(pii,pji,2,k,m,n)
                            df_pref_r(i,j,klmn,4)=df_pref_r(i,j,klmn,4) + hv_pref(pii,pji,k,1,l,n)
                            df_pref_r(i,j,klmn,6)=df_pref_r(i,j,klmn,6) + hv_pref(pii,pji,k,2,l,n)
                         endif
                         !boundary conditions in v_parallel
                         if((l.eq.0)) df_pref_r(i,j,klmn,-6:-3)=0.
                         if((l.eq.1)) df_pref_r(i,j,klmn,-6)=0.
                         if((l.eq.(nv0-2))) df_pref_r(i,j,klmn,6)=0.
                         if((l.eq.(nv0-1))) df_pref_r(i,j,klmn,3:6)=0.
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
  end subroutine initialize_dzv_2
  
  subroutine equ_dzv_r_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lij0,lzvwn0), intent(in):: p_f
    complex, dimension(lij0, lb1:lb2), intent(inout):: p_rhs
    real, dimension(lklmn0,-nsten:nsten)::  prefac
   
    integer:: klmn,sten

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,klmn)=prefac(klmn,-nsten)*p_f(:,shifted_f_pos(klmn,-nsten))
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             call axpy_ij(lij0,prefac(klmn,sten),p_f(:,shifted_f_pos(klmn,sten)),&
                  &p_rhs(:,klmn))
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             call axpy_ij(lij0,prefac(klmn,sten),p_f(:,shifted_f_pos(klmn,sten)),&
                  &p_rhs(:,klmn))
          enddo
       enddo
    endif
    PERFOFF_I
  end subroutine equ_dzv_r_2

  subroutine equ_dzv_xdep_r_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lklmn0,-nsten:nsten)::  prefac

    integer:: j,klmn,sten

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          do j=lj1,lj2
             p_rhs(:,j,klmn)=prefac(:,klmn,-nsten)*p_f(:,j,shifted_f_pos(klmn,-nsten))
          enddo
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,klmn,sten)*p_f(:,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             do j=lj1,lj2
                p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,klmn,sten)*p_f(:,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xdep_r_2

  subroutine equ_dzv_xydep_r_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lj1:lj2,lklmn0,-nsten:nsten)::  prefac

    integer:: klmn,sten

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,:,klmn)=prefac(:,:,klmn,-nsten)*p_f(:,:,shifted_f_pos(klmn,-nsten))
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xydep_r_2

  subroutine equ_dzv_xydep_c_2(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,lklmn0,-nsten:nsten)::  prefac

    integer:: klmn,sten

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          p_rhs(:,:,klmn)=prefac(:,:,klmn,-nsten)*p_f(:,:,shifted_f_pos(klmn,-nsten))
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,klmn,sten)*p_f(:,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xydep_c_2

  subroutine finalize_dzv_2
    if (shifted_metric) then
       deallocate(df_pref_c,shifted_f_pos)
    else
       deallocate(df_pref_r,shifted_f_pos)
    endif
  end subroutine finalize_dzv_2

  !> add hyper diffusion in z direction
  !! by adding hyp_z_order'th derivative of f1 (not g1)  
  subroutine add_hypz_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: p_f_
    complex, dimension(li1:li2,lj1:lj2),intent(inout) :: localrhs
    integer,intent(in):: k, l, m, n
    
    integer :: j, sten
    if (.not.shifted_metric) then
       if (hyp_on_h) then
          Do sten=-hyp_z_order/2,hyp_z_order/2
             if (xy_local) then
                call axpy_ij(lij0,hz_pref(li1,lj1,sten,k,m,n),&
                     p_f_(:,:,k+sten,l,m,n),localrhs)
             else
                do j=lj1,lj2
                   localrhs(:,j)=localrhs(:,j)+hz_pref(:,lj1,sten,k,m,n)*p_f_(:,j,k+sten,l,m,n)
                end do
             end if
          End Do
       else
          Do sten=-hyp_z_order/2,hyp_z_order/2
             call axpy_ij(lij0,hyp_z_scl*hyp_z_sten(sten),p_f_(:,:,k+sten,l,m,n),&
                  localrhs)
          End Do
       endif
    else
       if (hyp_on_h) then
          do sten=-hyp_z_order/2,hyp_z_order/2
             do j=lj1,lj2
                localrhs(:,j)=localrhs(:,j)+hz_pref(:,lj1,sten,k,m,n)*geom%phasefac(pi1:pi2,j,k,sten)&
                     *p_f_(:,j,k+sten,l,m,n)
             end do
          enddo
       else
          do sten=-hyp_z_order/2,hyp_z_order/2
             localrhs(:,:)=localrhs(:,:)+hyp_z_scl*hyp_z_sten(sten)*&
                  &geom%phasefac(pi1:pi2,:,k,sten)*p_f_(:,:,k+sten,l,m,n)
          enddo
       endif
    endif
  end subroutine add_hypz_indices

  subroutine add_hypz_ak(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypz_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypz_ak

  subroutine add_hypv_indices(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n

    integer :: sten, kb, lb, ub, j

    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)
    if (hyp_on_h) then
       do sten=lb,ub
          if (xy_local) then
             call axpy_ij(lij0,hv_pref(li1,lj1,k,sten,l,n),p_f_(li1,lj1,k,l+sten,m,n),localrhs(li1,lj1))
          else
             do j=lj1,lj2
                localrhs(:,j)=localrhs(:,j)+hv_pref(:,lj1,k,sten,l,n)*p_f_(:,j,k,l+sten,m,n)
             end do
          end if
       enddo
    else
       do sten=lb,ub
          call axpy_ij(lij0,vdamping(sten,l),p_f_(li1,lj1,k,l+sten,m,n),localrhs(li1,lj1))
       enddo
    endif
    
  end subroutine add_hypv_indices

  subroutine add_hypv_ak(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypv_indices(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypv_ak

  subroutine initialize_hv_pref
    integer:: i,j,k,l,n
    
    allocate(hv_pref(pi1:pi2,pj1:pj2,lk1:lk2,-2:2,ll1:ll2,ln1:ln2))
    do n=ln1,ln2
       do l=ll1,ll2
          do k=lk1,lk2
             do j=pj1,pj2
                do i=pi1,pi2
#ifdef hyp_over_j
                   hv_pref(i,j,k,-2,l,n) = -hyp_v/geom%jacobian(i,j,k)*dc2coeff(-1)**2&
                        *exp((vplocal(l-2)**2-vplocal(l-1)**2)/spec(n)%temp_prof(i))
                   
                   hv_pref(i,j,k,-1,l,n) = -hyp_v/geom%jacobian(i,j,k)*dc2coeff(-1)*dc2coeff(0)&
                        *(1+exp((vplocal(l-1)**2-vplocal(l)**2)/spec(n)%temp_prof(i)))
                   
                   hv_pref(i,j,k,0,l,n) = -hyp_v/geom%jacobian(i,j,k)*(dc2coeff(-1)*dc2coeff(1)&
                        *(exp((vplocal(l)**2-vplocal(l+1)**2)/spec(n)%temp_prof(i))&
                        +exp((vplocal(l)**2-vplocal(l-1)**2)/spec(n)%temp_prof(i)))&
                        +dc2coeff(0)**2)
                   
                   hv_pref(i,j,k,1,l,n) = -hyp_v/geom%jacobian(i,j,k)*dc2coeff(1)*dc2coeff(0)&
                        *(1+exp((vplocal(l+1)**2-vplocal(l)**2)/spec(n)%temp_prof(i)))
                   
                   hv_pref(i,j,k,2,l,n) = -hyp_v/geom%jacobian(i,j,k)*dc2coeff(1)**2&
                        *exp((vplocal(l+2)**2-vplocal(l+1)**2)/spec(n)%temp_prof(i))
#else
                   hv_pref(i,j,k,-2,l,n) = -hyp_v*dc2coeff(-1)**2&
                        *exp((vplocal(l-2)**2-vplocal(l-1)**2)/spec(n)%temp_prof(i))
                   
                   hv_pref(i,j,k,-1,l,n) = -hyp_v*dc2coeff(-1)*dc2coeff(0)&
                        *(1+exp((vplocal(l-1)**2-vplocal(l)**2)/spec(n)%temp_prof(i)))
                   
                   hv_pref(i,j,k,0,l,n) = -hyp_v*(dc2coeff(-1)*dc2coeff(1)&
                        *(exp((vplocal(l)**2-vplocal(l+1)**2)/spec(n)%temp_prof(i))&
                        +exp((vplocal(l)**2-vplocal(l-1)**2)/spec(n)%temp_prof(i)))&
                        +dc2coeff(0)**2)
                   
                   hv_pref(i,j,k,1,l,n) = -hyp_v*dc2coeff(1)*dc2coeff(0)&
                        *(1+exp((vplocal(l+1)**2-vplocal(l)**2)/spec(n)%temp_prof(i)))
                
                   hv_pref(i,j,k,2,l,n) = -hyp_v*dc2coeff(1)**2&
                        *exp((vplocal(l+2)**2-vplocal(l+1)**2)/spec(n)%temp_prof(i))
#endif
                end do
             end do
          end do
       end do
    end do

  end subroutine initialize_hv_pref

  subroutine initialize_hz_pref
    integer:: i,j,k,m,n

    allocate(hz_pref(pi1:pi2,pj1:pj2,-2:2,lk1:lk2,lm1:lm2,ln1:ln2))
    do n=ln1,ln2
       do m=lm1,lm2
          do k=lk1,lk2
             do j=pj1,pj2
                do i=pi1,pi2
#ifdef hyp_over_j
                   hz_pref(i,j,-2,k,m,n)=-hyp_z_scl/geom%jacobian(i,j,k)*dc2coeff(-1)**2&
                        *exp(mu(m)*(Blocal(i,j,k-2)-Blocal(i,j,k-1))/spec(n)%temp_prof(i))
                   
                   hz_pref(i,j,-1,k,m,n)=-hyp_z_scl/geom%jacobian(i,j,k)*dc2coeff(-1)*dc2coeff(0)&
                        *(1+exp(mu(m)*(Blocal(i,j,k-1)-Blocal(i,j,k))/spec(n)%temp_prof(i)))
                   
                   hz_pref(i,j,0,k,m,n)=-hyp_z_scl/geom%jacobian(i,j,k)*(dc2coeff(-1)*dc2coeff(1)&
                        *(exp(mu(m)*(Blocal(i,j,k)-Blocal(i,j,k+1))/spec(n)%temp_prof(i))&
                        +exp(mu(m)*(Blocal(i,j,k)-Blocal(i,j,k-1))/spec(n)%temp_prof(i)))&
                        +dc2coeff(0)**2)
                   
                   hz_pref(i,j,1,k,m,n)=-hyp_z_scl/geom%jacobian(i,j,k)*dc2coeff(1)*dc2coeff(0)&
                        *(1+exp(mu(m)*(Blocal(i,j,k+1)-Blocal(i,j,k))/spec(n)%temp_prof(i)))
                   
                   hz_pref(i,j,2,k,m,n)=-hyp_z_scl/geom%jacobian(i,j,k)*dc2coeff(1)**2&
                        *exp(mu(m)*(Blocal(i,j,k+2)-Blocal(i,j,k+1))/spec(n)%temp_prof(i))
#else
                   hz_pref(i,j,-2,k,m,n)=-hyp_z_scl*dc2coeff(-1)**2&
                        *exp(mu(m)*(Blocal(i,j,k-2)-Blocal(i,j,k-1))/spec(n)%temp_prof(i))
                   
                   hz_pref(i,j,-1,k,m,n)=-hyp_z_scl*dc2coeff(-1)*dc2coeff(0)&
                        *(1+exp(mu(m)*(Blocal(i,j,k-1)-Blocal(i,j,k))/spec(n)%temp_prof(i)))
                   
                   hz_pref(i,j,0,k,m,n)=-hyp_z_scl*(dc2coeff(-1)*dc2coeff(1)&
                        *(exp(mu(m)*(Blocal(i,j,k)-Blocal(i,j,k+1))/spec(n)%temp_prof(i))&
                        +exp(mu(m)*(Blocal(i,j,k)-Blocal(i,j,k-1))/spec(n)%temp_prof(i)))&
                        +dc2coeff(0)**2)
                   
                   hz_pref(i,j,1,k,m,n)=-hyp_z_scl*dc2coeff(1)*dc2coeff(0)&
                        *(1+exp(mu(m)*(Blocal(i,j,k+1)-Blocal(i,j,k))/spec(n)%temp_prof(i)))
                   
                   hz_pref(i,j,2,k,m,n)=-hyp_z_scl*dc2coeff(1)**2&
                        *exp(mu(m)*(Blocal(i,j,k+2)-Blocal(i,j,k+1))/spec(n)%temp_prof(i))
#endif
                end do
             end do
          end do
       end do
    end do

  end subroutine initialize_hz_pref

  subroutine initialize_hypz_compensation
    integer:: j,k,m,n, pni

    if (shifted_metric) then
       allocate(hz_comp_pref_c(pi1:pi2,lj1:lj2,-2:2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
       pni=ln1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      hz_comp_pref_c(:,j,-2,k,l,m,n)=-hz_pref(:,lj1,-2,k,m,n)*fm(:,lj1,k-2,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)&
                           *geom%phasefac(pi1:pi2,j,k,-2)
                      hz_comp_pref_c(:,j,-1,k,l,m,n)=-hz_pref(:,lj1,-1,k,m,n)*fm(:,lj1,k-1,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)&
                           *geom%phasefac(pi1:pi2,j,k,-1)
                      hz_comp_pref_c(:,j,0,k,l,m,n)=-hz_pref(:,lj1,0,k,m,n)*fm(:,lj1,k,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                      hz_comp_pref_c(:,j,1,k,l,m,n)=-hz_pref(:,lj1,1,k,m,n)*fm(:,lj1,k+1,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)&
                           *geom%phasefac(pi1:pi2,j,k,1)
                      hz_comp_pref_c(:,j,2,k,l,m,n)=-hz_pref(:,lj1,2,k,m,n)*fm(:,lj1,k+2,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)&
                           *geom%phasefac(pi1:pi2,j,k,2)
                   end do
                end do
             end do
          end do
       enddo
    else
       allocate(hz_comp_pref(pi1:pi2,pj1:pj2,-2:2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
       pni=ln1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=pj1,pj2
                      hz_comp_pref(:,j,-2,k,l,m,n)=-hz_pref(:,j,-2,k,m,n)*fm(:,j,k-2,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                      hz_comp_pref(:,j,-1,k,l,m,n)=-hz_pref(:,j,-1,k,m,n)*fm(:,j,k-1,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                      hz_comp_pref(:,j,0,k,l,m,n)=-hz_pref(:,j,0,k,m,n)*fm(:,j,k,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                      hz_comp_pref(:,j,1,k,l,m,n)=-hz_pref(:,j,1,k,m,n)*fm(:,j,k+1,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                      hz_comp_pref(:,j,2,k,l,m,n)=-hz_pref(:,j,2,k,m,n)*fm(:,j,k+2,l,m,pni)&
                           *spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(pi1:pi2)
                   end do
                end do
             end do
          end do
       enddo
    endif
    
  end subroutine initialize_hypz_compensation

  subroutine equ_comp_hypz(p_emfields,p_bar_emfields,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs

    if(xy_local) then
       call equ_compensate_hypz(p_emfields,p_rhs,lb1,lb2,hz_comp_pref)
    else
       if (shifted_metric) then
          call equ_compensate_hypz_xydep_c(p_bar_emfields,p_rhs,lb1,lb2,hz_comp_pref_c)
       else
          call equ_compensate_hypz_xdep(p_bar_emfields,p_rhs,lb1,lb2,hz_comp_pref)
       endif
    end if

  end subroutine equ_comp_hypz

  subroutine equ_compensate_hypz(p_emfields,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz):: gy_xi
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(-2:2,lklmn0)::  prefac
   
    integer:: klmn,sten

    PERFON_I('cmphypz')

       do klmn=lb1,lb2
          if ((compute_gy_av(klmn)).or.(klmn.eq.lb1)) then
             !gyroaverage
             gy_xi = p_emfields(:,:,:,1) * jfac(:,:,:,sm(klmn),sn(klmn))
             if (n_fields .gt. 2) gy_xi = gy_xi + mu_tjqj(sm(klmn),sn(klmn)) * &
                  p_emfields(:,:,:,3) * I1_factor(:,:,:,sm(klmn),sn(klmn))
          endif
          do sten=-2,2
             call axpy_ij(lij0,prefac(sten,klmn),gy_xi(:,:,sk(klmn)+sten),&
                  &p_rhs(:,:,klmn))
          enddo
       enddo

    PERFOFF_I
  end subroutine equ_compensate_hypz

  subroutine equ_compensate_hypz_xdep(p_bar_emfields,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-2:2,lklmn0)::  prefac
   
    integer:: j,klmn,sten

    PERFON_I('cmphypz')

    do klmn=lb1,lb2
       do sten=-2,2
          do j=lj1,lj2
             p_rhs(:,j,klmn)=p_rhs(:,j,klmn)+prefac(:,sten,klmn)*p_bar_emfields(:,j,sk(klmn)+sten,sm(klmn),sn(klmn),1)
          enddo
       enddo
    enddo
    
    PERFOFF_I
  end subroutine equ_compensate_hypz_xdep

  subroutine equ_compensate_hypz_xydep_c(p_bar_emfields,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-2:2,lklmn0)::  prefac
    integer:: klmn,sten

    PERFON_I('cmphypz')

       do klmn=lb1,lb2
          do sten=-2,2
             p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+prefac(:,:,sten,klmn)*p_bar_emfields(:,:,sk(klmn)+sten,sm(klmn),sn(klmn),1)
          enddo
       enddo

    PERFOFF_I
  end subroutine equ_compensate_hypz_xydep_c
  
  subroutine equ_dzv_adptv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs

    if(xy_local) then
       stop 'no adaptive block structured grid for local simulations!'
    else
       if (shifted_metric) then
          if (perf_vec(3).eq.1) then
             call equ_dzv_xydep_c_1_adptv(p_f,p_rhs,lb1,lb2,df_pref_c)
          else
             call equ_dzv_xydep_c_2_adptv(p_f,p_rhs,lb1,lb2,df_pref_c)
          end if
       else
          if (perf_vec(3).eq.1) then
             call equ_dzv_xdep_r_1_adptv(p_f,p_rhs,lb1,lb2,df_pref_r)
          else
             call equ_dzv_xdep_r_2_adptv(p_f,p_rhs,lb1,lb2,df_pref_r)
          end if
       end if
    end if

  end subroutine equ_dzv_adptv

  subroutine equ_dzv_xdep_r_1_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-6:6,lklmn0), intent(in) ::  prefac
    integer:: j, klmn, sten, f_ind, l, m, i1, i2

    PERFON_I('dzv1')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          do j=lj1,lj2
             p_rhs(i1:i2,j,klmn)= &
                  & prefac(i1:i2,stenbound(1,klmn),klmn)* &
                  & p_f(i1:i2,j,f_ind)
          enddo
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn)= &
                     & p_rhs(i1:i2,j,klmn)+ &
                     & prefac(i1:i2,sten,klmn)* &
                     & p_f(i1:i2,j,f_ind)
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
                p_rhs(i1:i2,j,klmn)= &
                     & p_rhs(i1:i2,j,klmn)+ &
                     & prefac(i1:i2,sten,klmn)* &
                     & p_f(i1:i2,j,f_ind)
             enddo
          enddo
       enddo
    endif
    PERFOFF_I
    
  end subroutine equ_dzv_xdep_r_1_adptv

  subroutine equ_dzv_xydep_c_1_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-6:6,lklmn0), intent(in) ::  prefac
    integer:: klmn, sten, f_ind, l, m, i1, i2
    
    PERFON_I('dzv1')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          f_ind=map_to_f(klmn)+f_shifts(stenbound(1,klmn))
          p_rhs(i1:i2,:,klmn)= &
               prefac(i1:i2,:,stenbound(1,klmn),klmn)* &
               p_f(i1:i2,:,f_ind)
          do sten=stenbound(1,klmn)+1,stenbound(2,klmn)
             f_ind=map_to_f(klmn)+f_shifts(sten)
             p_rhs(i1:i2,:,klmn)= &
                  p_rhs(i1:i2,:,klmn)+ &
                  prefac(i1:i2,:,sten,klmn)* &
                  p_f(i1:i2,:,f_ind)
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
             p_rhs(i1:i2,:,klmn)= &
                  p_rhs(i1:i2 ,:,klmn)+ &
                  prefac(i1:i2,:,sten,klmn)* &
                  p_f(i1:i2,:,f_ind)
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xydep_c_1_adptv

  subroutine equ_dzv_xdep_r_2_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,lklmn0,-nsten:nsten)::  prefac

    integer:: j,klmn,sten, l, m, i1, i2

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          do j=lj1,lj2
             p_rhs(i1:i2,j,klmn)= &
                  & prefac(i1:i2,klmn,-nsten)* &
                  & p_f(i1:i2,j,shifted_f_pos(klmn,-nsten))
          enddo
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn) = &
                     & p_rhs(i1:i2,j,klmn)+ &
                     & prefac(i1:i2,klmn,sten)* &
                     & p_f(i1:i2,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn)= &
                     & p_rhs(i1:i2,j,klmn)+ &
                     & prefac(i1:i2,klmn,sten)* &
                     & p_f(i1:i2,j,shifted_f_pos(klmn,sten))
             enddo
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xdep_r_2_adptv
  
  subroutine equ_dzv_xydep_c_2_adptv(p_f,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,lklmn0,-nsten:nsten)::  prefac

    integer:: klmn, sten, l, m, i1, i2

    PERFON_I('dzv2')
    if(replace_rhs) then
       do klmn=lb1,lb2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p_rhs(i1:i2,:,klmn)= &
               & prefac(i1:i2,:,klmn,-nsten)* &
               & p_f(i1:i2,:,shifted_f_pos(klmn,-nsten))
       enddo
       do sten=-nsten+1,nsten
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p_rhs(i1:i2,:,klmn)= &
                  & p_rhs(i1:i2,:,klmn)+ &
                  & prefac(i1:i2,:,klmn,sten)* &
                  & p_f(i1:i2,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    else
       do sten=-nsten,nsten
          do klmn=lb1,lb2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p_rhs(i1:i2,:,klmn)= &
                  & p_rhs(i1:i2,:,klmn)+ &
                  & prefac(i1:i2,:,klmn,sten)* &
                  & p_f(i1:i2,:,shifted_f_pos(klmn,sten))
          enddo
       enddo
    endif
    PERFOFF_I

  end subroutine equ_dzv_xydep_c_2_adptv

  !> add hyper diffusion in z direction
  !! by adding hyp_z_order'th derivative of f1 (not g1)  
  subroutine add_hypz_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in) :: p_f_
    complex, dimension(li1:li2,lj1:lj2),intent(inout) :: localrhs
    integer,intent(in):: k, l, m, n
    integer :: j, sten, i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)
    
    if (.not.shifted_metric) then
       if (hyp_on_h) then
          Do sten=-hyp_z_order/2,hyp_z_order/2
             if (xy_local) then
                stop "no adaptivity for local xy version!"
             else
                do j=lj1,lj2
                   localrhs(i1:i2,j)=localrhs(i1:i2,j)+&
                        &hz_pref(i1:i2,lj1,sten,k,m,n)*&
                        &p_f_(i1:i2,j,k+sten,l,m,n)
                end do
             end if
          End Do
       else
          Do sten=-hyp_z_order/2,hyp_z_order/2
             call axpy_ij(lij0,hyp_z_scl*hyp_z_sten(sten),p_f_(i1:i2,:,k+sten,l,m,n),&
                  localrhs)
          End Do
       endif
    else
       if (hyp_on_h) then
          do sten=-hyp_z_order/2,hyp_z_order/2
             do j=lj1,lj2
                localrhs(i1:i2,j)=localrhs(i1:i2,j)+&
                     &hz_pref(i1:i2,lj1,sten,k,m,n)*&
                     &geom%phasefac(pi1_vwadp(l,m):pi2_vwadp(l,m),j,k,sten)*&
                     &p_f_(i1:i2,j,k+sten,l,m,n)
             end do
          enddo
       else
          do sten=-hyp_z_order/2,hyp_z_order/2
             localrhs(i1:i2,:)=localrhs(i1:i2,:)+&
                  &hyp_z_scl*hyp_z_sten(sten)*geom%phasefac(pi1_vwadp(l,m):pi2_vwadp(l,m),:,k,sten)*&
                  &p_f_(i1:i2,:,k+sten,l,m,n)
          enddo
       endif
    endif
  end subroutine add_hypz_indices_adptv

  subroutine add_hypz_ak_adptv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypz_indices_adptv(p_f, p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypz_ak_adptv

  subroutine add_hypv_indices_adptv(p_f_,localrhs,k,l,m,n)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f_
    complex, dimension(li1:li2,lj1:lj2), intent(inout):: localrhs
    integer, intent(in) :: k, l, m, n

    integer :: sten, kb, lb, ub, j, i1, i2

    i1 = li1_vwadp(l,m)
    i2 = li2_vwadp(l,m)

    kb=2
    lb=max(-l,-kb)
    ub=min(nv0-1-l,kb)
    if (hyp_on_h) then
       do sten=lb,ub
          if (xy_local) then
             call axpy_ij(lij0,hv_pref(i1,lj1,k,sten,l,n),&
                  &p_f_(i1,lj1,k,l+sten,m,n),localrhs(i1,lj1))
          else
             do j=lj1,lj2
                localrhs(i1:i2,j)=localrhs(i1:i2,j)+&
                     &hv_pref(i1:i2,lj1,k,sten,l,n)*p_f_(i1:i2,j,k,l+sten,m,n)
             end do
          end if
       enddo
    else
       do sten=lb,ub
          call axpy_ij(lij0,vdamping(sten,l),p_f_(i1,lj1,k,l+sten,m,n),localrhs(i1,lj1))
       enddo
    endif
    
  end subroutine add_hypv_indices_adptv
  
  subroutine add_hypv_ak_adptv(p_f,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lzvwn0), intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs
    
    integer :: klmn

    do klmn=lb1,lb2
       call add_hypv_indices_adptv(p_f,p_rhs(:,:,klmn),sk(klmn),sl(klmn),sm(klmn),sn(klmn))
    enddo

  end subroutine add_hypv_ak_adptv

  subroutine equ_comp_hypz_adptv(p_emfields,p_bar_emfields,p_rhs,lb1,lb2)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2,lb1:lb2), intent(inout):: p_rhs

    if(xy_local) then
       stop "no adaptivity for local grid version"
    else
       if (shifted_metric) then
          call equ_compensate_hypz_xydep_c_adptv(p_bar_emfields,p_rhs,lb1,lb2,hz_comp_pref_c)
       else
          call equ_compensate_hypz_xdep_adptv(p_bar_emfields,p_rhs,lb1,lb2,hz_comp_pref)
       endif
    end if

  end subroutine equ_comp_hypz_adptv

  subroutine equ_compensate_hypz_xdep_adptv(p_bar_emfields,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    real, dimension(li1:li2,-2:2,lklmn0)::  prefac
   
    integer:: j,klmn,sten

    ! block-structured grid stuff
    integer:: l, m, i1, i2

    PERFON_I('cmphypz')

    do klmn=lb1,lb2
       do sten=-2,2
          do j=lj1,lj2
             l = sl(klmn)
             m = sm(klmn)
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             p_rhs(i1:i2,j,klmn)=p_rhs(i1:i2,j,klmn)+prefac(i1:i1,sten,klmn)&
                  &*p_bar_emfields(i1:i2,j,sk(klmn)+sten,sm(klmn),sn(klmn),1)
          enddo
       enddo
    enddo
    
    PERFOFF_I
  end subroutine equ_compensate_hypz_xdep_adptv

  subroutine equ_compensate_hypz_xydep_c_adptv(p_bar_emfields,p_rhs,lb1,lb2,prefac)
    integer, intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer :: p_bar_emfields
    complex, dimension(li1:li2,lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(li1:li2,lj1:lj2,-2:2,lklmn0)::  prefac
    integer:: klmn,sten
    ! block-structured grid stuff
    integer:: l, m, i1, i2

    PERFON_I('cmphypz')
    
    do klmn=lb1,lb2
       do sten=-2,2
          l = sl(klmn)
          m = sm(klmn)
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p_rhs(i1:i2,:,klmn)=p_rhs(i1:i2,:,klmn)+prefac(i1:i2,:,sten,klmn)&
               &*p_bar_emfields(i1:i2,:,sk(klmn)+sten,sm(klmn),sn(klmn),1)
       enddo
    enddo
    
    PERFOFF_I
  end subroutine equ_compensate_hypz_xydep_c_adptv
 

  subroutine add_Erad_acc(p_f, lb1, lb2, p_rhs)
!Parallel acceleration due to Erad 
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2), intent(in):: p_f
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    integer:: sten, klmn, n, m,l, k, j
   
     do n=ln1,ln2
        do m=lm1,lm2 
           do l=ll1,ll2
              do k=lk1,lk2
                 klmn = (n-ln1)*lklm0 + (m-lm1)*lkl0 + (l-ll1)*lk0 + (k-lk1)+1
                 do j =lj1, lj2
                    do sten=ll_bound(l),-1
                       p_rhs(li1:li2,j,klmn) = p_rhs(li1:li2,j,klmn) + Erad_array(li1:li2,pj1,k,l) * & 
                            & vderivative(sten,:,k,l) * p_f(li1:li2,j,k,l+sten,m,n)
                    enddo
                    do sten=1,ul_bound(l)
                       p_rhs(li1:li2,j,klmn) = p_rhs(li1:li2,j,klmn) + Erad_array(li1:li2,pj1,k,l) * & 
                            & vderivative(sten,:,k,l) * p_f(li1:li2,j,k,l+sten,m,n)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
   end subroutine add_Erad_acc

   subroutine initialize_Erad_acc 
     integer:: j,k,l
     if (.not.allocated(ll_bound))  allocate(ll_bound(li1:li2))
     if (.not.allocated(ul_bound))  allocate(ul_bound(li1:li2))
     allocate(Erad_array(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2))


     !initialize time-independent part of Erad_acc
     do l=ll1,ll2
        do k=lk1,lk2
           do j =lj1, lj2
              Erad_array(li1:li2, pj1, k, l) = vp(l) * rhostar * &
                   & Erad / geom%Bfield(li1:li2,pj1, k) * geom%K_j(li1:li2,pj1,k)
           enddo
        enddo
     enddo

     !initialize stencil boundaries for vpar derivatives
    ll_bound(:) = -2
    ul_bound(:) = +2
  
    if (li1.eq.0) ll_bound(li1) = 0
    if (li1.eq.1) ll_bound(li1) = -1

    if(li2.eq.(nv0-1)) ul_bound(li2) = 0
    if(li2.eq.(nv0-2)) ul_bound(li2) = 1

   end subroutine initialize_Erad_acc

   subroutine finalize_Erad_acc 
     
     deallocate(Erad_array)
     deallocate(ll_bound)
     deallocate(ul_bound)

   end subroutine finalize_Erad_acc


end module dzv_terms
