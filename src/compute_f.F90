#include "redef.h"
#include "intrinsic_sizes.h"

module compute_f
  
  use par_mod
  use communications
  use vel_space
  use gyro_average_df_mod, only: gyro_average_df,&
       &gyro_average_df_wb
  use gyro_average_ff_mod, only: gyro_average_ff,jfac, I1_factor
  use boundaries
  use antenna, only: antenna_type, Apar_pre_antenna
  use discretization_adptv_module

  implicit none
  public:: initialize_calc_f, calc_f, finalize_calc_f, mem_est_compute_f
  public:: calc_h_from_f
  public:: f_,h_
  ! routines for adptive grids
  public:: calc_f_adptv, calc_h_from_f_adptv

  private

  integer :: init_status = 0

  complex, dimension(:,:,:,:,:,:), allocatable::  f_
  complex, dimension(:,:,:,:,:,:), allocatable::  h_
  real, dimension(:,:,:,:,:,:), allocatable :: papbar, phi_prefac, bpar_prefac
  real, dimension(:,:,:,:,:,:), allocatable :: gy_pre   
  complex,dimension(:,:,:,:), allocatable :: gy_field_4d !used as Apar or Phi
  integer:: lbijk, ubijk

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_compute_f(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    !f_
    mem_loc=SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lv0*lw0*ln0
    if (n_fields.gt.1) then
       !papbar
       mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*pj0*lklmn0
       if (perf_vec(6)==2) then
          !gy_pre
          mem_loc=mem_loc+SIZE_OF_REAL_MB*lijklmn0
       else
          !gy_field_4d
          if (.not.xy_local) mem_loc=&
               &mem_loc+SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lm0
       endif
    endif
    if (arakawa_zv) then
       !phi_prefac
       mem_loc=mem_loc+SIZE_OF_REAL_MB*li0*lj0*lz0*lv0*lm0*ln0
       !bpar_prefac
       if(n_fields.gt.2) &
            mem_loc=mem_loc+SIZE_OF_REAL_MB*li0*lj0*lz0*lv0*lm0*ln0
       !gy_field_4d(in phiterm)
       if (.not.xy_local) mem_loc=&
           mem_loc+SIZE_OF_COMPLEX_MB*li0*lj0*lz0*lm0
    endif

    mem_est_compute_f=mem_req_in+mem_loc
  End Function mem_est_compute_f 


  subroutine initialize_calc_f
    
    if ((init_status.ne.perf_vec(6)).and.(init_status.gt.0)) &
         & call finalize_calc_f

    call initialize_exchange_z5d 
    call initialize_exchange_v

    if (init_status==0) then
       allocate(f_(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2))
       !initialize f_ (some boundary points in vp are not initialized otherwise)
       f_=0.
       if (arakawa_zv) then
          allocate(h_(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2))
          !initialize h_ (some boundary points in vp are not initialized otherwise)
          h_=0.
       else
          allocate(h_(1:1,1:1,1:1,1:1,1:1,1:1))
       endif

       if (arakawa_zv.or.(istep_fe_transfer.gt.0)) then
          call initialize_h_from_f
       endif
       
       if(n_fields.gt.1) then
          call initialize_papbar
          if(perf_vec(6).eq.1) call initialize_calc_f_em_1
          if(perf_vec(6).eq.2) call initialize_calc_f_em_2
       endif

       init_status = perf_vec(6)
    endif

    
  end subroutine initialize_calc_f

  subroutine calc_f(p_g,p_emfields,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: p_emfields
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz):: p_apar
    PERFON('calc_f')
    if (n_fields.eq.1) then
       if (perf_vec(6).eq.1) then
          call calc_f_es_1(p_g,p_f)
       else
          call calc_f_es_2(p_g,p_f)
       endif
    else
       if (antenna_type.ne.3) then
          p_apar=p_emfields(:,:,:,2)
       else
          p_apar=Apar_pre_antenna
       endif
       if (perf_vec(6).eq.1) then
          call calc_f_em_1(p_g,p_apar,p_f)
       else
          call calc_f_em_2(p_g,gy_pre,p_apar,p_f)
       endif
    endif

    call exchange_z(p_f)
    !boundaries of h_ are exchanged later on

    if (n_procs_v.gt.1) call exchange_v(p_f)
    PERFOFF
  end subroutine calc_f

  subroutine finalize_calc_f
    deallocate(f_)
    deallocate(h_)
    if(n_fields.gt.1) then
       call finalize_papbar
       if(init_status.eq.1) call finalize_calc_f_em_1
       if(init_status.eq.2) call finalize_calc_f_em_2
    endif

    if (arakawa_zv.or.(istep_fe_transfer.gt.0)) then 
       call finalize_h_from_f
    endif

    call finalize_exchange_z5d 
    call finalize_exchange_v
    init_status = 0

  end subroutine finalize_calc_f

!---------

  !>Initializes the prefactor for computation of $q/T\cdot F_M\phi$, the difference
  !! between f and h.
  subroutine initialize_h_from_f
    integer:: pni, j, k, l, m, n

    allocate(phi_prefac(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lm1:lm2,ln1:ln2))
    if (n_fields.gt.2) allocate(bpar_prefac(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lm1:lm2,ln1:ln2))
    if (.not.xy_local.and.n_fields.eq.1) &
         & allocate(gy_field_4d(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2))
    pni = pn1

    !phi_prefac and bpar_prefac are initialized with correct v and z boundaries 
    if (xy_local) then
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          Do m=lm1,lm2
             do l=lbv,ubv
                Do k=lbz,ubz
                   phi_prefac(:,:,k,l,m,n)=spec(n)%charge/spec(n)%temp&
                        *fm(pi1,pj1,k,l,m,pni)*jfac(:,:,k,m,n)
                   if (n_fields.gt.2) &
                        &bpar_prefac(:,:,k,l,m,n)=mu(m)*&
                        &fm(pi1,pj1,k,l,m,pni)*I1_factor(:,:,k,m,n)
                End Do
             End Do
          enddo
       enddo
    else
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             do l=lbv,ubv
                Do k=lbz,ubz
                   Do j=lj1,lj2
                      phi_prefac(:,j,k,l,m,n)=spec(n)%charge/spec(n)%temp/spec(n)%temp_prof(li1:li2)&
                           *fm(pi1:pi2,pj1,k,l,m,pni)
                   End Do
                End Do
             End Do
          enddo
       enddo
    endif


  end subroutine initialize_h_from_f

  !> Computes an h-distribution from an f-distribution by adding a term proportional to the
  !! electrostatic potential and b_parallel (if pesent).
  subroutine calc_h_from_f(p_f,p_emfields,p_h)
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: p_emfields
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_h
    integer:: l, m, n, pni

    PERFON('h_from_f')

    if (xy_local) then
       if (n_fields.gt.2) then
          do n=ln1,ln2
#if WITHOMP_COMPUTE_F
          !$omp parallel do default(none)&
          !$omp shared(phi_prefac,bpar_prefac,p_f,p_h)&
          !$omp shared(n,lm1,lm2,lbv,ubv,p_emfields)&
          !$omp private(m,l)
#endif
             do m=lm1,lm2
                do l=lbv,ubv
                   p_h(:,:,:,l,m,n) =p_f(:,:,:,l,m,n)+phi_prefac(:,:,:,l,m,n)*p_emfields(:,:,:,1)&
                        &+bpar_prefac(:,:,:,l,m,n)*p_emfields(:,:,:,3)
                end do
             end do
#if WITHOMP_COMPUTE_F
          !$omp end parallel do
#endif
          end do
       else
          do n=ln1,ln2
#if WITHOMP_COMPUTE_F
          !$omp parallel do default(none)&
          !$omp shared(phi_prefac,p_f,p_h)&
          !$omp shared(n,lm1,lm2,lbv,ubv,p_emfields)&
          !$omp private(m,l)
#endif
             do m=lm1,lm2
                do l=lbv,ubv
                   p_h(:,:,:,l,m,n) =p_f(:,:,:,l,m,n)+phi_prefac(:,:,:,l,m,n)*p_emfields(:,:,:,1)
                end do
             end do
#if WITHOMP_COMPUTE_F
          !$omp end parallel do
#endif
          end do
       endif
    else
       pni = pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             CALL gyro_average_df_wb(p_emfields(:,:,:,1),gy_field_4d(:,:,:,m),m,n)          
             !in nonlocal mode we require the z exchange of phibar during runtime
             call exchange_z(gy_field_4d(:,:,:,m))
          enddo
          do m=lm1,lm2
             do l=lbv,ubv
                p_h(:,:,:,l,m,n)=p_f(:,:,:,l,m,n)+phi_prefac(:,:,:,l,m,n)*gy_field_4d(:,:,:,m)
             End Do
          enddo
       enddo
    endif

    PERFOFF

  end subroutine calc_h_from_f


  subroutine finalize_h_from_f

    deallocate(phi_prefac)
    if (n_fields.gt.2) deallocate(bpar_prefac)
    if (.not.xy_local.and.n_fields.eq.1) deallocate(gy_field_4d)

  end subroutine finalize_h_from_f

  subroutine initialize_calc_f_em_1
    if (.not.xy_local) &
         & allocate(gy_field_4d(li1:li2,lj1:lj2,lbz:ubz,lm1:lm2))
  end subroutine initialize_calc_f_em_1

  subroutine calc_f_em_1(p_g,p_Apar,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz), intent(in) :: p_Apar
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f

    integer:: j,k,l,m,n

    PERFON_I('cfem1')

    if (xy_local) then
       do n=ln1,ln2
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   p_f(:,:,k,l,m,n) =&
                        &p_g(:,:,k,l,m,n) + papbar(pi1,pj1,k,l,m,n) * jfac(:,:,k,m,n) * p_Apar(:,:,k)
                end do
             end do
          end do
       end do
    else
       do n=ln1,ln2
          do m=lm1,lm2
             CALL gyro_average_df_wb(p_Apar(:,:,:),gy_field_4d(:,:,:,m),m,n)          
          enddo
          do m=lm1,lm2
             do l=ll1,ll2
                do k=lk1,lk2
                   do j=lj1,lj2
                      p_f(li1:li2,j,k,l,m,n) = p_g(li1:li2,j,k,l,m,n) + &
                           & papbar(li1:li2,pj1,k,l,m,n) * gy_field_4d(li1:li2,j,k,m)
                   end do
                end do
             end do
          end do
       end do
    end if

    PERFOFF_I
  end subroutine calc_f_em_1

  subroutine finalize_calc_f_em_1
    if (.not.xy_local) deallocate(gy_field_4d)
  end subroutine finalize_calc_f_em_1

  subroutine initialize_calc_f_em_2
    integer:: k,l,m,n

   allocate(gy_pre(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                gy_pre(:,:,k,l,m,n)=jfac(:,:,k,m,n)*papbar(pi1,pj1,k,l,m,n)
             enddo
          enddo
       enddo
    enddo

    lbijk=li0*lj0*nzb
    ubijk=li0*lj0*(lk0+nzb)-1
  end subroutine initialize_calc_f_em_2

  subroutine calc_f_em_2(g_s,pre_s,psi_s,f_s)
    complex,dimension(0:li0*lj0*lk0-1,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: g_s
    complex,dimension(0:li0*lj0*lz0-1,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: f_s
    complex,dimension(0:li0*lj0*lz0-1), intent(in) :: psi_s
    real,dimension(0:li0*lj0*lk0-1,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: pre_s

    integer:: l,m,n

    PERFON_I('cfem2')
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             f_s(lbijk:ubijk,l,m,n) =&
                  &g_s(:,l,m,n) + pre_s(:,l,m,n)* psi_s(lbijk:ubijk)
          End Do
       enddo
    enddo
    PERFOFF_I

  end subroutine calc_f_em_2

  subroutine finalize_calc_f_em_2
    deallocate(gy_pre)
  end subroutine finalize_calc_f_em_2

  subroutine calc_f_es_1(p_g,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f
    integer:: k,l,m,n

    PERFON_I('cfes1')

    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                call ccopy(lij0,p_g(li1,lj1,k,l,m,n),1,p_f(li1,lj1,k,l,m,n),1)
             end do
          end do
       end do
    end do

    PERFOFF_I

  end subroutine calc_f_es_1
  
  subroutine calc_f_es_2(p_g,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f

    PERFON_I('cfes2')
    p_f(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)=p_g
    PERFOFF_I

  end subroutine calc_f_es_2


  subroutine initialize_papbar
    integer:: j,k,l,m,n, pni
    
    allocate(papbar(pi1:pi2, pj1:pj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))

    pni = pn1
    do n=ln1,ln2
       if (pn0.gt.1) pni=n
       do m=lm1,lm2
          do l=ll1,ll2
             do k=lk1,lk2
                do j=pj1,pj2
                   papbar(:,j,k,l,m,n) = - spec(n)%charge*sqrt(2./(spec(n)%temp*spec(n)%mass))/ &
                        & spec(n)%temp_prof(pi1:pi2) * vp(l) * fm(:,j,k,l,m,pni)
                enddo
             enddo
          enddo
       enddo
    enddo

  end subroutine initialize_papbar

  subroutine finalize_papbar
    deallocate(papbar)
  end subroutine finalize_papbar

  ! modified routines for adaptive grids

  subroutine calc_f_adptv(p_g,p_emfields,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: p_emfields
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f

    PERFON('calc_f')
    if (n_fields.eq.1) then
       if (perf_vec(6).eq.1) then
          call calc_f_es_1_adptv(p_g,p_f)
       else
          call calc_f_es_2_adptv(p_g,p_f)
       endif
    else
       if (perf_vec(6).eq.1) then
          call calc_f_em_1_adptv(p_g,p_emfields(:,:,:,2),p_f)
       else
          call calc_f_em_2_adptv(p_g,gy_pre,p_emfields(:,:,:,2),p_f)
       endif
    endif

    call exchange_z(p_f)
    !boundaries of h_ are exchanged later on

    if (n_procs_v.gt.1) call exchange_v(p_f)
    PERFOFF
  end subroutine calc_f_adptv

  subroutine calc_f_es_1_adptv(p_g,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f
    integer:: k,l,m,n, i1, i0

    PERFON_I('cfes1')

    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             i0 = li0_vwadp(l,m)
             i1 = li1_vwadp(l,m)
             do k=lk1,lk2
                call ccopy(lj0*i0,p_g(i1,lj1,k,l,m,n),1,p_f(i1,lj1,k,l,m,n),1)
             end do
          end do
       end do
    end do

    PERFOFF_I

  end subroutine calc_f_es_1_adptv

  subroutine calc_f_es_2_adptv(p_g,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f
    integer:: l, m, i1, i2

    PERFON_I('cfes2')
    do m = lm1, lm2
       do l = ll1, ll2
          i1 = li1_vwadp(l,m)
          i2 = li2_vwadp(l,m)
          p_f(i1:i2,lj1:lj2,lk1:lk2,l,m,ln1:ln2)=p_g(i1:i2,:,:,l,m,:)
       end do
    end do
    PERFOFF_I

  end subroutine calc_f_es_2_adptv

  subroutine calc_f_em_1_adptv(p_g,p_Apar,p_f)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: p_g
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz), intent(in) :: p_Apar
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_f

    integer:: j,k,l,m,n, i1, i2

    PERFON_I('cfem1')

    if (xy_local) then
       stop "no adaptive grids for the xy local version"
    else
       ! TODO: return to block-structured grid and gyro-average computation after consultation!
       do n=ln1,ln2
          do m=lm1,lm2
             CALL gyro_average_df_wb(p_Apar(:,:,:),gy_field_4d(:,:,:,m),m,n)
          enddo
          do m=lm1,lm2
             do l=ll1,ll2
                i1 = li1_vwadp(l,m)
                i2 = li2_vwadp(l,m)
                do k=lk1,lk2
                   do j=lj1,lj2
                      p_f(i1:i2,j,k,l,m,n) = p_g(i1:i2,j,k,l,m,n) + &
                           & papbar(i1:i2,pj1,k,l,m,n) * gy_field_4d(i1:i2,j,k,m)
                   end do
                end do
             end do
          end do
       end do
    end if

    PERFOFF_I
  end subroutine calc_f_em_1_adptv

  subroutine calc_f_em_2_adptv(g_s,pre_s,psi_s,f_s)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: g_s
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: f_s
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz), intent(in) :: psi_s
    real,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: pre_s

    integer:: l,m,n, i1, i2

    PERFON_I('cfem2')
    do n=ln1,ln2
       do m=lm1,lm2
          do l=ll1,ll2
             i1 = li1_vwadp(l,m)
             i2 = li2_vwadp(l,m)
             f_s(i1:i2,:,:,l,m,n) =&
                  &g_s(i1:i2,:,:,l,m,n) + pre_s(i1:i2,:,:,l,m,n)*psi_s(i1:i2,:,:)
          End Do
       enddo
    enddo
    PERFOFF_I

  end subroutine calc_f_em_2_adptv

  !> Computes an h-distribution from an f-distribution by adding a term proportional to the
  !! electrostatic potential and b_parallel (if pesent).
  subroutine calc_h_from_f_adptv(p_f,p_emfields,p_h)
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(in):: p_f
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(in) :: p_emfields
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: p_h
    integer:: l, m, n, pni, i1, i2

    PERFON('h_from_f')

    if (xy_local) then
       stop "no adaptive grids for the xy local version"
    else
       pni = pn1
       do n=ln1,ln2
          if (pn0.gt.1) pni=n
          do m=lm1,lm2
             CALL gyro_average_df_wb(p_emfields(:,:,:,1),gy_field_4d(:,:,:,m),m,n)          
             !in nonlocal mode we require the z exchange of phibar during runtime
             call exchange_z(gy_field_4d(:,:,:,m))
          enddo
          do m=lm1,lm2
             do l=lbv,ubv
                i1 = li1_b_vwadp(l,m)
                i2 = li2_b_vwadp(l,m)
                p_h(i1:i2,:,:,l,m,n)=p_f(i1:i2,:,:,l,m,n)+phi_prefac(i1:i2,:,:,l,m,n)&
                     &*gy_field_4d(i1:i2,:,:,m)
             end do
          end do          
       enddo
    endif

    PERFOFF

  end subroutine calc_h_from_f_adptv

end module compute_f
