#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"

!>Computes the terms containing x and y derivatives of chi and g_1 (i.e. the drive and curvature terms)
module dchidxy_terms
  use par_mod
  use discretization
  use discretization_adptv_module
  use blockindex
  use coordinates
  use prefactors
  use external_contr
  use gyro_average_ff_mod, only: jfac, I1_factor
  use gyro_average_df_mod, only: gyro_average_df
  use x_derivatives, only: x_deriv_exc, x_deriv, x_deriv_1d_adptv

  implicit none
  public:: initialize_dchidxy, add_dchidxy, finalize_dchidxy, &
       &mem_est_dchidxy
  public :: add_dchidxy_orig
  ! routines for adaptive grids
  public :: add_dchidxy_adptv
  private

  integer :: init_status = 0
  complex,dimension(:,:),allocatable:: iki, ikj

  complex, dimension(:,:,:),allocatable:: chi_pref, p_bar_chi_3d
  complex, dimension(:,:,:,:),allocatable:: mp_dbar_chi_dxy

  !auxiliary array for add_dchidxy_2
  complex, dimension(:,:,:),allocatable :: chi

#ifdef WITHOMP_BLOCKLOOP
  !$OMP THREADPRIVATE(p_bar_chi_3d,mp_dbar_chi_dxy,chi)
#endif

contains

  function mem_est_dchidxy(mem_req_in)
    real:: mem_req_in, mem_est_dchidxy
    real:: mem_loc

    if (xy_local) then
       if (perf_vec(2).eq.1) then
          mem_loc=mem_est_dchidxy_1()
       else 
          mem_loc=mem_est_dchidxy_2()
          if (istep_energy.gt.0.or.istep_energy3d.gt.0) &
               &mem_loc = mem_loc + mem_est_dchidxy_1()
       end if
    else
       mem_loc=mem_est_dchidxy_df_1()
    endif

    mem_est_dchidxy=mem_req_in+mem_loc

  end function mem_est_dchidxy


  subroutine initialize_dchidxy
    if ((init_status.ne.perf_vec(2)).and.(init_status.gt.0)) &
         call finalize_dchidxy

    if (init_status.eq.0) then
       if (xy_local) then
          !if (perf_vec(2).eq.1.or.istep_energy.gt.0) then
          if (perf_vec(2).eq.1) then
             call initialize_dchidxy_1
          else
             call initialize_dchidxy_2
             if (istep_energy.gt.0.or.istep_energy3d.gt.0.or.&
                  &diag_GyroLES) call initialize_dchidxy_1
          end if
       else
          call initialize_dchidxy_df_1
       endif
       
       init_status = perf_vec(2)
    endif

  end subroutine initialize_dchidxy

  subroutine add_dchidxy(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lbi:ubi,lj1:lj2,1:lbg0), intent(IN) :: chi_block
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer :: p_dbarchidxy
    logical, intent(in) :: with_curv

    ! this routine is only called for .not.xy_local
    if (xy_local) then
       if ((perf_vec(2).eq.1).or.(.not.with_curv)) then
          !call add_dchidxy_1(p_emfields,p_rhs,p_barchi,lb1,lb2,with_curv)
       else
          !call add_dchidxy_2(p_emfields,p_rhs,lb1,lb2,chi)
!          p_barchi=>NULL()
       end if
    else
       call add_dchidxy_df_1(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    endif
    
  end subroutine add_dchidxy

  subroutine add_dchidxy_orig(p_emfields,p_bar_emfields,p_rhs,p_barchi,p_dbarchidxy,lb1,lb2,with_curv)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields), intent(in):: p_emfields
    complex, dimension(:,:,:,:,:,:), pointer:: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(:,:,:), pointer:: p_barchi
    complex, dimension(:,:,:,:), pointer:: p_dbarchidxy
    logical, intent(in) :: with_curv

    if (xy_local) then
       if ((perf_vec(2).eq.1).or.(.not.with_curv)) then
          call add_dchidxy_1(p_emfields,p_rhs,p_barchi,lb1,lb2,with_curv)
       else
          call add_dchidxy_2(p_emfields,p_rhs,lb1,lb2,chi)
       end if
    else
       call add_dchidxy_df_1_orig(p_bar_emfields,p_rhs,p_barchi,&
            &p_dbarchidxy,lb1,lb2,with_curv)
    endif
    
  end subroutine add_dchidxy_orig


  subroutine finalize_dchidxy

    if (xy_local) then
       !if (init_status.eq.1.or.istep_energy.gt.0) then
       if (init_status.eq.1) then
          call finalize_dchidxy_1
       else
          call finalize_dchidxy_2
          if (istep_energy.gt.0.or.istep_energy3d.gt.0.or.&
               &diag_GyroLES) call finalize_dchidxy_1
       end if
    else
       call finalize_dchidxy_df_1
    endif

    init_status = 0

  end subroutine finalize_dchidxy

!-------------------------------------------------------------------------------------------------

  function mem_est_dchidxy_1()
    real:: mem_est_dchidxy_1

    !iki, ikj
    mem_est_dchidxy_1=2.*lij0*SIZE_OF_COMPLEX_MB

    !p_bar_chi_3d
    mem_est_dchidxy_1 = mem_est_dchidxy_1 + lijk0*SIZE_OF_COMPLEX_MB
    
  end function mem_est_dchidxy_1

  subroutine initialize_dchidxy_1
    Integer :: i,j

    allocate(iki(li1:li2,lj1:lj2),ikj(li1:li2,lj1:lj2))

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) shared(li1,li2,lj1,lj2,lk1,lk2)
    !$OMP CRITICAL
#endif
    allocate(p_bar_chi_3d(li1:li2,lj1:lj2,lk1:lk2))
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif

    do j=lj1,lj2
       do i=li1,li2
          iki(i,j)=imag*ki(i)
          ikj(i,j)=imag*kj(j)
       enddo
    enddo

  end subroutine initialize_dchidxy_1

  subroutine add_dchidxy_1(p_emfields,p_rhs,p_barchi,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields), intent(in):: p_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(:,:,:), pointer :: p_barchi
    logical, intent(in) :: with_curv
    complex, dimension(li1:li2, lj1:lj2) :: localpref
    integer:: klmn, klmn1
    
    PERFON_I('dcdxy1')

    do klmn=lb1,lb2
       if (with_curv) then
          localpref=( pdchibardj(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * ikj +&
               &pdchibardi(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * iki)
       else
          if (xy_local.and.yx_order) then
             !drive term ~ky=ki
             localpref=edr_for_energy(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn))* iki
          else
             localpref=edr_for_energy(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn))* ikj
          endif
       endif

       if (.not.associated(p_barchi)) then
          !for linear, single field simulations, we can skip the "l" loop
          !and we do not save bar_chi
          if (n_fields.gt.1) then
             p_bar_chi_3d(:,:,sk(klmn))=(p_emfields(li1:li2,lj1:lj2,sk(klmn),1) - vTvpar(sl(klmn),sn(klmn)) * &
                  p_emfields(li1:li2,lj1:lj2,sk(klmn),2) )*jfac(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn))
             if (n_fields.gt.2) p_bar_chi_3d(:,:,sk(klmn)) = p_bar_chi_3d(:,:,sk(klmn)) + &
                  & mu_tjqj(sm(klmn),sn(klmn)) * p_emfields(li1:li2,lj1:lj2,sk(klmn),3)*&
                  I1_factor(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn))
          else
#ifdef WITHOMP_BLOCKLOOP
             if (compute_gy_av(klmn) .or. (klmn.eq.lb1)) then
#else
             if (compute_gy_av(klmn)) then  !ignore vpar loop
#endif
                p_bar_chi_3d(:,:,lk1:lk2)=p_emfields(li1:li2,lj1:lj2,lk1:lk2,1)*&
                     &jfac(li1:li2,lj1:lj2,lk1:lk2,sm(klmn),sn(klmn))
             end if
          endif

          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + localpref * p_bar_chi_3d(:,:,sk(klmn))

          if ((with_omt_ext).or.(with_omn_ext)) then  
             if (with_omt_ext) call add_external_contr(p_rhs(:,:,klmn),&
                  &1,p_bar_chi_3d(:,:,sk(klmn)),kxind_omt_ext,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
             if (with_omn_ext) call add_external_contr(p_rhs(:,:,klmn),&
                  &2,p_bar_chi_3d(:,:,sk(klmn)), kxind_omn_ext, sk(klmn),sl(klmn),sm(klmn),sn(klmn))
          endif
       else
          klmn1 = klmn-lb1+1

          if (n_fields.gt.1) then
             p_barchi(:,:,klmn1)=(p_emfields(li1:li2,lj1:lj2,sk(klmn),1) - &
                  & vTvpar(sl(klmn),sn(klmn)) * p_emfields(li1:li2,lj1:lj2,sk(klmn),2) )*&
                  & jfac(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn))
             if (n_fields.gt.2) p_barchi(:,:,klmn1) = p_barchi(:,:,klmn1) + &
                  & mu_tjqj(sm(klmn),sn(klmn)) * p_emfields(li1:li2,lj1:lj2,sk(klmn),3)*&
                  I1_factor(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn))
          else
             p_barchi(:,:,klmn1)=p_emfields(li1:li2,lj1:lj2,sk(klmn),1)*&
                  &jfac(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn))
          endif

          p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + localpref * p_barchi(:,:,klmn1)

          if ((with_omt_ext).or.(with_omn_ext)) then  
             if (with_omt_ext) call add_external_contr(p_rhs(:,:,klmn),&
                  &1,p_barchi(:,:,klmn1),kxind_omt_ext,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
             if (with_omn_ext) call add_external_contr(p_rhs(:,:,klmn),&
                  &2,p_barchi(:,:,klmn1), kxind_omn_ext, sk(klmn),sl(klmn),sm(klmn),sn(klmn))
          endif
       endif
    enddo

    PERFOFF_I

  end subroutine add_dchidxy_1

  subroutine finalize_dchidxy_1
    deallocate(iki,ikj)

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL
#endif
    deallocate(p_bar_chi_3d)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END PARALLEL
#endif
  end subroutine finalize_dchidxy_1


!------------------------------------------------------------------------------------------------------


  function mem_est_dchidxy_2()
    real:: mem_est_dchidxy_2

    !chi_pref
    mem_est_dchidxy_2=lijklmn0*SIZE_OF_COMPLEX_MB
    
    !chi
    mem_est_dchidxy_2=mem_est_dchidxy_2+lij0*lklmn0/nblocks*SIZE_OF_COMPLEX_MB

  end function mem_est_dchidxy_2

  subroutine initialize_dchidxy_2
    integer:: i,j,klmn
    
    allocate(chi_pref(li1:li2,lj1:lj2,1:lklmn0))
    
    IF (n_fields .GT. 2) THEN
      do klmn=1,lklmn0
        do j=lj1,lj2
          do i=li1,li2
            chi_pref(i,j,klmn)= imag*( pdchibardj(pi1,pj1,sk(klmn),sl(klmn),&
              & sm(klmn),sn(klmn)) * kj(j) &
              & + pdchibardi(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * ki(i))
          enddo
        enddo
      enddo
    ELSE
      do klmn=1,lklmn0
         do j=lj1,lj2
            do i=li1,li2
               chi_pref(i,j,klmn)= imag*( pdchibardj(pi1,pj1,sk(klmn),sl(klmn),&
                    & sm(klmn),sn(klmn)) * kj(j) &
                    & + pdchibardi(pi1,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * ki(i)) &
                    & *jfac(i,j,sk(klmn),sm(klmn),sn(klmn))
            enddo
         enddo
      enddo
    END IF

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL shared(li1,li2,lj1,lj2,lbg0)
    !$OMP CRITICAL
#endif
    allocate(chi(li1:li2,lj1:lj2,1:lbg0))
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif

  end subroutine initialize_dchidxy_2

  subroutine add_dchidxy_2(p_emfields,p_rhs,lb1,lb2,p_chi)
    integer,intent(in) :: lb1,lb2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz,1:n_fields), intent(in):: p_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs, p_chi

    complex, dimension(li1:li2, lj1:lj2) :: tmparr

    integer:: klmn

    PERFON_I('dcdxy2')    
   
    if (n_fields .LT. 2) then
       do klmn=lb1,lb2
          p_chi(:,:,klmn)=p_emfields(:,:,sk(klmn),1)
          !drive,curvature and pressure
          p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+chi_pref(:,:,klmn)*p_chi(:,:,klmn)
       enddo
    else IF (n_fields .EQ. 2) THEN
       do klmn=lb1,lb2
          p_chi(:,:,klmn)=p_emfields(:,:,sk(klmn),1)-&
               &vTvpar(sl(klmn),sn(klmn))*p_emfields(:,:,sk(klmn),2)
          !drive,curvature and pressure
          p_rhs(:,:,klmn)=p_rhs(:,:,klmn)+chi_pref(:,:,klmn)*p_chi(:,:,klmn)
       enddo
    ELSE
      DO klmn = lb1, lb2
        p_chi(:,:,klmn) = (p_emfields(:,:,sk(klmn),1) - &
          vTvpar(sl(klmn),sn(klmn)) * p_emfields(:,:,sk(klmn),2)) * &
          jfac(:,:,sk(klmn),sm(klmn),sn(klmn)) + &
          mu_tjqj(sm(klmn),sn(klmn)) * p_emfields(:,:,sk(klmn),3) * &
          I1_factor(:,:,sk(klmn),sm(klmn),sn(klmn))
        p_rhs(:,:,klmn) = p_rhs(:,:,klmn) + chi_pref(:,:,klmn) * p_chi(:,:,klmn)
      END DO
    END IF

    ! (radial derivative) contributions of external equilibrium 
    ! temperature and/or density variations
    ! (needed for tertiary instability analysis)
    if ((with_omt_ext).or.(with_omn_ext)) then     
       do klmn=lb1,lb2
          tmparr = p_chi(li1:li2,lj1:lj2,klmn)
          if (n_fields.lt.2) tmparr = tmparr*jfac(li1:li2,lj1:lj2,sk(klmn),&
                  &sm(klmn),sn(klmn))
             
          if (with_omt_ext) call add_external_contr(&
               &p_rhs(li1:li2,lj1:lj2,klmn),1,tmparr,&
               &kxind_omt_ext,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
          if (with_omn_ext) call add_external_contr(&
               &p_rhs(li1:li2,lj1:lj2,klmn),2,tmparr,&
               &kxind_omn_ext,sk(klmn),sl(klmn),sm(klmn),sn(klmn))
       enddo
    endif
  
    PERFOFF_I

  end subroutine add_dchidxy_2

  subroutine finalize_dchidxy_2
    deallocate(chi_pref)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL
    !$OMP CRITICAL
#endif
    deallocate(chi)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif
  end subroutine finalize_dchidxy_2

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

  function mem_est_dchidxy_df_1()
    real:: mem_est_dchidxy_df_1

    mem_est_dchidxy_df_1=(3.*lijk0+2.*lij0+lx0*lj0)*SIZE_OF_COMPLEX_MB

  end function mem_est_dchidxy_df_1
  
  subroutine initialize_dchidxy_df_1
    Integer :: i,j

    allocate(ikj(li1:li2,lj1:lj2))

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) shared(li1,li2,lj1,lj2,lk1,lk2,n_fields,lbg0)
    !$OMP CRITICAL
#endif
    !allocate(p_bar_chi_3d(li1:li2,lj1:lj2,lk1:lk2))
    if (n_fields.gt.1) then
       allocate(mp_dbar_chi_dxy(li1:li2,lj1:lj2,1:lbg0,2))
    else
       allocate(mp_dbar_chi_dxy(li1:li2,lj1:lj2,lk1:lk2,2))
    end if
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif

    do j=lj1,lj2
       do i=li1,li2
          ikj(i,j)=imag*kj(j)
       enddo
    enddo

  end subroutine initialize_dchidxy_df_1

  subroutine add_dchidxy_df_1(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lbi:ubi,lj1:lj2,1:lbg0), intent(IN) :: chi_block
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer :: p_dbarchidxy
    logical, intent(in) :: with_curv

    !complex, dimension(lbi:ubi, lj1:lj2) :: bar_chi_2d_wb

    integer:: j,klmn
    
    PERFON_I('dcdxyg1')

    !ve_x,ve_y are not needed for linear runs

    do klmn=1,lbg0
       !chi terms
       if (n_fields.gt.1) then
          ! this is the case when chi depends on x,y,z,v_par,mu and species
          call x_deriv(chi_block(:,:,klmn), mp_dbar_chi_dxy(:,:,klmn,1))
          mp_dbar_chi_dxy(:,:,klmn,2) = ikj(:,:) * chi_block(li1:li2,lj1:lj2,klmn)

       else
!if (compute_gy_av(klmn+lb1-1)&
#ifdef WITHOMP_BLOCKLOOP
!            &.or.(klmn.eq.1)&
#endif
!            &) then
          ! this is the case when chi is independent of v_par (it is equal to phi_bar)
          !write(*,"(3I3,ES17.10)") mype,lb1+klmn-1,sk(klmn-1+lb1),&
          !     &real(sum(conjg(chi_block(li1:li2,lj1:lj2,klmn))*chi_block(li1:li2,lj1:lj2,klmn)))
          call x_deriv(chi_block(:,:,klmn),mp_dbar_chi_dxy(:,:,sk(klmn-1+lb1),1))
          mp_dbar_chi_dxy(:,:,sk(klmn-1+lb1),2) = ikj(:,:) * chi_block(li1:li2,lj1:lj2,klmn)
          !write(*,"(3I3,2ES17.10)") mype,lb1+klmn-1,sk(klmn-1+lb1),&
          !     &real(sum(conjg(mp_dbar_chi_dxy(li1:li2,lj1:lj2,sk(klmn-1+lb1),1))&
          !     &*mp_dbar_chi_dxy(li1:li2,lj1:lj2,sk(klmn-1+lb1),1))),&
          !     &real(sum(conjg(mp_dbar_chi_dxy(li1:li2,lj1:lj2,sk(klmn-1+lb1),2))&
          !     &*mp_dbar_chi_dxy(li1:li2,lj1:lj2,sk(klmn-1+lb1),2)))
       endif

       if (n_fields.gt.1) then
          if (associated(p_dbarchidxy)) then
             !save bar_chi and dbar_chi_dxy
             !p_barchi(:,:,klmn-lb1+1) = p_bar_chi_3d(:,:,sk(klmn))
             p_dbarchidxy(:,:,1,klmn) = mp_dbar_chi_dxy(:,:,klmn,1)
             p_dbarchidxy(:,:,2,klmn) = mp_dbar_chi_dxy(:,:,klmn,2)
          end if

          if (with_curv) then
             do j=lj1,lj2
                p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                     &pdchibardj(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(:,j,klmn,2) +&
                     &pdchibardi(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(:,j,klmn,1)
             enddo
          else
             if (yx_order) then
                do j=lj1,lj2
                   p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                        &edr_for_energy(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(:,j,klmn,1)
                enddo
             else
                do j=lj1,lj2
                   p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                        &edr_for_energy(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(:,j,klmn,2)
                enddo
             endif
          endif
       else
          if (associated(p_dbarchidxy)) then
             !p_barchi(:,:,klmn-lb1+1) = p_bar_chi_3d(:,:,sk(klmn))
             p_dbarchidxy(:,:,1,klmn) = mp_dbar_chi_dxy(:,:,sk(klmn-1+lb1),1)
             p_dbarchidxy(:,:,2,klmn) = mp_dbar_chi_dxy(:,:,sk(klmn-1+lb1),2)
          end if

          if (with_curv) then
             do j=lj1,lj2
                p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                     &pdchibardj(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(:,j,sk(klmn-1+lb1),2) +&
                     &pdchibardi(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(:,j,sk(klmn-1+lb1),1)
             enddo
          else
             if (yx_order) then
                do j=lj1,lj2
                   p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                        &edr_for_energy(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(:,j,sk(klmn-1+lb1),1)
                enddo
             else
                do j=lj1,lj2
                   p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                        &edr_for_energy(:,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(:,j,sk(klmn-1+lb1),2)
                enddo
             endif
          endif
       endif
    enddo

    PERFOFF_I

  end subroutine add_dchidxy_df_1

  subroutine add_dchidxy_df_1_orig(p_bar_emfields,p_rhs,p_barchi,p_dbarchidxy,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(:,:,:,:,:,:), pointer:: p_bar_emfields
    complex, dimension(li1:li2, lj1:lj2, lb1:lb2), intent(inout):: p_rhs
    complex, dimension(:,:,:), pointer :: p_barchi
    complex, dimension(:,:,:,:), pointer :: p_dbarchidxy
    logical, intent(in) :: with_curv

    complex, dimension(lbi:ubi, lj1:lj2) :: bar_chi_2d_wb

    integer:: j, k, klmn
    
    PERFON_I('dcdxyg1')

    !ve_x,ve_y are not needed for linear runs

    do klmn=lb1,lb2
       !chi terms
       if (n_fields.gt.1) then
          p_bar_chi_3d(:,:,sk(klmn)) = &
               &p_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),1)&
               &- vTvpar(sl(klmn),sn(klmn)) * &
               & p_bar_emfields(li1:li2,lj1:lj2,sk(klmn),sm(klmn),sn(klmn),2) 
          do j=lj1,lj2
             bar_chi_2d_wb(li1:li2,j) = p_bar_chi_3d(li1:li2,j,sk(klmn))
          enddo
          call x_deriv_exc(bar_chi_2d_wb,mp_dbar_chi_dxy(:,:,sk(klmn),1))
          mp_dbar_chi_dxy(:,:,sk(klmn),2) = ikj(:,:) * p_bar_chi_3d(:,:,sk(klmn))

       elseif (compute_gy_av(klmn)&
#ifdef WITHOMP_BLOCKLOOP
            &.or.(klmn.eq.lb1)&
#endif
            &) then
          !write(*,"(A,I4,A,4I3)") "compute_gy_av(",klmn,") is true. klmn = ",sk(klmn),sl(klmn),sm(klmn),sn(klmn)
          !this part is independent of 'l' and thus only called in the 
          !m and n loops (which requires an additional k loop)
          do k=lk1,lk2
             p_bar_chi_3d(:,:,k) = &
                  &p_bar_emfields(li1:li2,lj1:lj2,k,sm(klmn),sn(klmn),1)            
             do j=lj1,lj2
                bar_chi_2d_wb(li1:li2,j) = p_bar_chi_3d(li1:li2,j,k)
             enddo
             call x_deriv_exc(bar_chi_2d_wb,mp_dbar_chi_dxy(:,:,k,1))
             mp_dbar_chi_dxy(:,:,k,2) = ikj(:,:) * p_bar_chi_3d(:,:,k)
          enddo
       endif

       if (associated(p_barchi)) then
          !save bar_chi and dbar_chi_dxy
          p_barchi(:,:,klmn-lb1+1) = p_bar_chi_3d(:,:,sk(klmn))
          p_dbarchidxy(:,:,1,klmn-lb1+1) = mp_dbar_chi_dxy(:,:,sk(klmn),1)
          p_dbarchidxy(:,:,2,klmn-lb1+1) = mp_dbar_chi_dxy(:,:,sk(klmn),2)
       endif
       
       if (with_curv) then
          do j=lj1,lj2
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                  &pdchibardj(:,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                  &mp_dbar_chi_dxy(:,j,sk(klmn),2) +&
                  &pdchibardi(:,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn)) * &
                  &mp_dbar_chi_dxy(:,j,sk(klmn),1)
          enddo
       else
          do j=lj1,lj2
             p_rhs(:,j,klmn) = p_rhs(:,j,klmn) + &
                  &edr_for_energy(:,pj1,sk(klmn),sl(klmn),sm(klmn),sn(klmn))*&
                  &mp_dbar_chi_dxy(:,j,sk(klmn),2)
          enddo
       endif
    enddo

    PERFOFF_I

  end subroutine add_dchidxy_df_1_orig


  subroutine finalize_dchidxy_df_1
    deallocate(ikj)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL 
    !$OMP CRITICAL
#endif
    deallocate(mp_dbar_chi_dxy)
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END CRITICAL
    !$OMP END PARALLEL
#endif
    
  end subroutine finalize_dchidxy_df_1

  ! ----- modifidied routines for adaptive grids ------

  subroutine add_dchidxy_adptv(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lbi:ubi,lj1:lj2,1:lbg0), intent(IN) :: chi_block
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer :: p_dbarchidxy
    logical, intent(in) :: with_curv

    ! this routine is only called for .not.xy_local
    if (xy_local) then
       stop "no adaptive grids for xy local version!"
    else
       call add_dchidxy_df_1_adptv(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    endif
    
  end subroutine add_dchidxy_adptv
  
  subroutine add_dchidxy_df_1_adptv(chi_block,p_rhs,p_dbarchidxy,lb1,lb2,with_curv)
    integer, intent(in) :: lb1,lb2
    complex, dimension(lbi:ubi,lj1:lj2,1:lbg0), intent(IN) :: chi_block
    complex, dimension(li1:li2, lj1:lj2, 1:lbg0), intent(inout):: p_rhs
    complex, dimension(:,:,:,:), pointer :: p_dbarchidxy
    logical, intent(in) :: with_curv

    integer:: j, klmn, l, m, i1, i2
    
    PERFON_I('dcdxyg1')

    !ve_x,ve_y are not needed for linear runs

    do klmn=1,lbg0
       l = sl(lb1+klmn-1)
       m = sm(lb1+klmn-1)
       i1 = li1_vwadp(l,m)
       i2 = li2_vwadp(l,m)
       !chi terms
       if (n_fields.gt.1) then
          ! this is the case when chi depends on x,y,z,v_par,mu and species
          call x_deriv_1d_adptv(chi_block(:,:,klmn), mp_dbar_chi_dxy(:,:,klmn,1),i1, i2)
          mp_dbar_chi_dxy(i1:i2,:,klmn,2) = ikj(i1:i2,:) * chi_block(i1:i2,lj1:lj2,klmn)

       else
          call x_deriv_1d_adptv(chi_block(:,:,klmn),mp_dbar_chi_dxy(:,:,sk(klmn-1+lb1),1), i1, i2)
          mp_dbar_chi_dxy(i1:i2,:,sk(klmn-1+lb1),2) = ikj(i1:i2,:) * chi_block(i1:i2,lj1:lj2,klmn)
       endif

       if (n_fields.gt.1) then
          if (associated(p_dbarchidxy)) then
             p_dbarchidxy(i1:i2,:,1,klmn) = mp_dbar_chi_dxy(i1:i2,:,klmn,1)
             p_dbarchidxy(i1:i2,:,2,klmn) = mp_dbar_chi_dxy(i1:i2,:,klmn,2)
          end if

          if (with_curv) then
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                     &pdchibardj(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(i1:i2,j,klmn,2) +&
                     &pdchibardi(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(i1:i2,j,klmn,1)
             enddo
          else
             if (yx_order) then
                do j=lj1,lj2
                   p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                        &edr_for_energy(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(i1:i2,j,klmn,1)
                enddo
             else
                do j=lj1,lj2
                   p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                        &edr_for_energy(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(i1:i2,j,klmn,2)
                enddo
             endif
          endif
       else
          if (associated(p_dbarchidxy)) then
             p_dbarchidxy(i1:i2,:,1,klmn) = mp_dbar_chi_dxy(i1:i2,:,sk(klmn-1+lb1),1)
             p_dbarchidxy(i1:i2,:,2,klmn) = mp_dbar_chi_dxy(i1:i2,:,sk(klmn-1+lb1),2)
          end if

          if (with_curv) then
             do j=lj1,lj2
                p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                     &pdchibardj(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(i1:i2,j,sk(klmn-1+lb1),2) +&
                     &pdchibardi(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1)) * &
                     &mp_dbar_chi_dxy(i1:i2,j,sk(klmn-1+lb1),1)
             enddo
          else
             if (yx_order) then
                do j=lj1,lj2
                   p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                        &edr_for_energy(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(i1:i2,j,sk(klmn-1+lb1),1)
                enddo
             else
                do j=lj1,lj2
                   p_rhs(i1:i2,j,klmn) = p_rhs(i1:i2,j,klmn) + &
                        &edr_for_energy(i1:i2,pj1,sk(klmn-1+lb1),sl(klmn-1+lb1),sm(klmn-1+lb1),sn(klmn-1+lb1))*&
                        &mp_dbar_chi_dxy(i1:i2,j,sk(klmn-1+lb1),2)
                enddo
             endif
          endif
       endif
    enddo

    PERFOFF_I

  end subroutine add_dchidxy_df_1_adptv

end module dchidxy_terms
