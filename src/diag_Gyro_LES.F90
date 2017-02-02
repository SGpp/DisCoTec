#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! Fre energy versus time
! Free energy spectral analysis (kx,ky,kperp ...)
! Transfer function
Module diag_Gyro_LES
  
  Use par_mod 
  use diag_Gyro_LES_common
  use diag_Gyro_LES_fetime
  use diag_Gyro_LES_spec2d
  use diag_Gyro_LES_transfer
  use diag_Gyro_LES_transfer3d

  Implicit None
    
Contains

!!!******************************************************************!!!
!!!******************************************************************!!!
  !>Give an estimate of the memory requirements of diag_GyroLES
  Real Function mem_est_diag_GyroLES(mem_req_in)
    real, parameter:: mem_cn=SIZE_OF_COMPLEX/(1024.)**2,&
         & mem_rn=SIZE_OF_REAL/(1024.)**2
    real:: mem_req_in
    real:: mem_loc=0.
    real:: mem_loc_max=0 !one diag a time, arrays are deallocated

    mem_loc = 0.
    if (istep_fe_time.gt.0) then
       !g_last,g_rhs,diss_rhs,temp_rhs,temp_fe6d
       !cfgammas
       mem_loc = 9*lijklmn0*mem_cn
       !f_last, h_last
       mem_loc = mem_loc + lij0*lzvwn0*mem_cn
       !if (arakawa_zv) mem_loc = mem_loc + lij0*lzvwn0*mem_cn
       !emfields_last  
       mem_loc = mem_loc + lij0*lz0*n_fields*mem_rn
       if (fracy.gt.0.and.nonlinear) then
          mem_loc = mem_loc + 5*lijklmn0*mem_cn
          mem_loc = mem_loc + lx0*lj0*lzvwn0*mem_cn
          mem_loc = mem_loc + lij0*lz0*n_fields*mem_rn
       endif
    endif
    mem_loc_max=max(mem_loc,mem_loc_max)

    mem_loc = 0.
    if (istep_fe_twoD.gt.0) then
       !v7d,v7d_spec,g_rhs,gc_rhs,d_v_rhs,d_z_rhs,d_kperp_rhs,coll_rhs,p_rhs,c_rhs...
       !gfgamma012
       mem_loc = 6*lijklmn0*mem_cn
       !phi_term
       mem_loc = mem_loc + lij0*lzvwn0*mem_cn
       !cfgamma3
       mem_loc = mem_loc + lij0*lz0*ln0*mem_cn
    endif
    mem_loc_max=max(mem_loc,mem_loc_max)

    mem_loc = 0.
    if (istep_fe_transfer.gt.0) then
       !g_ls1,g_ls2,nlt_ls2_g,nlt_ls2_f, calc_transfer arrays
       mem_loc = 6*lijklmn0*mem_cn
       !f_ls2
       mem_loc = mem_loc + lij0*lzvwn0*mem_cn
       !fields_ls2
       mem_loc = mem_loc + lij0*lz0*ln0*mem_rn
       !to be considered?
       !allocate(fe_tf_tot(n_shells,n_shells,ln1:ln2,n_transfer_models))
       !allocate(fe_tot(n_shells,ln1:ln2,3,8))
       !allocate(k_value(1:n_shells))
    endif
    mem_loc_max=max(mem_loc,mem_loc_max)
    
    mem_est_diag_GyroLES = mem_loc_max+mem_req_in

  End Function mem_est_diag_GyroLES


  Subroutine initialize_all_diags_GyroLES

#ifdef with_extended_diags
    IF (istep_fe_time.gt.0) call initialize_diag_fe_time
    IF (istep_fe_twoD.gt.0) call initialize_diag_spectral_twoD
    IF (istep_fe_transfer.gt.0) call initialize_diag_transfer_general
#else
    if (mype==0) Write(*,'(a)') 'ERROR: with_extended_diags must be set in switches.h for using diag_GyroLES'
    stop
#endif

  End Subroutine initialize_all_diags_GyroLES

  Subroutine exec_all_diags_GyroLES(itime,time)

    Integer, intent(in) :: itime
    Real, intent(in) :: time
    logical :: no_aux_yet
    no_aux_yet = .true.

    if (istep_fe_time.gt.0) then
       if (itime.gt.0) then
          if (MODULO(itime,istep_fe_time) == 0) then
             if(arakawa_zv) then 
                call calc_aux_fields(g_1,emfields,f_,.true.,h_)  !Calculate f and h
             else
                call calc_aux_fields(g_1,emfields,f_,.false.)    !Calculate f
             endif
             no_aux_yet = .false.
             CALL diagnostic_fe_time(2)
          else if (MODULO(itime + 1,istep_fe_time) == 0) then
             if(arakawa_zv) then 
                call calc_aux_fields(g_1,emfields,f_,.true.,h_)  !Calculate f and h
             else
                call calc_aux_fields(g_1,emfields,f_,.false.)    !Calculate f
             endif
             no_aux_yet = .false.
             CALL diagnostic_fe_time(1)
          endif
       endif
    endif

    IF (istep_fe_twoD.gt.0.and.itime.gt.0) THEN
       if (modulo(itime,istep_fe_twoD) .eq. 0) then
          if(arakawa_zv) then 
             call calc_aux_fields(g_1,emfields,f_,.true.,h_)  !Calculate f and h
          else
             call calc_aux_fields(g_1,emfields,f_,.false.)    !Calculate f
          endif
          no_aux_yet = .false.
          call diag_spectral_twoD
       endif
    endif

    if (istep_fe_transfer.gt.0) then
         if (itime.gt.0) then    
            if (modulo(itime,istep_fe_transfer) .eq. 0) then
                !?call calc_aux_fields(g_1,emfields,f_,.true.,h_)  !Calculate h
                !?no_aux_yet = .false.
                call diagnostic_transfer_general
            endif
        endif
    endif

  End Subroutine exec_all_diags_GyroLES

  Subroutine finalize_all_diags_GyroLES

    IF (istep_fe_time.gt.0) call finalize_diag_fe_time
    IF (istep_fe_twoD.gt.0) call finalize_diag_spectral_twoD
    IF (istep_fe_transfer.gt.0) call finalize_diag_transfer_general

  End Subroutine finalize_all_diags_GyroLES

  Subroutine initialize_diag_transfer_general

    if (any((/tky,tkx,tkperp,t_log_kperp,t_log_ky,t_log_kperp_1,t_log_kperp_fb,t_log_kperp_h/))) &
         call initialize_diag_transfer_std
    if (any((/tk3_ky,tk3_kx,tk3_log_kperp,tk3_log_kperp_1,tk3_log_kperp_h/))) &
         call initialize_diag_transfer_3d

  End Subroutine initialize_diag_transfer_general

  Subroutine diagnostic_transfer_general
    
    if (any((/tky,tkx,tkperp,t_log_kperp,t_log_ky,t_log_kperp_1,t_log_kperp_fb,t_log_kperp_h/)))&
         call diagnostic_transfer_std
    if (any((/tk3_ky,tk3_kx,tk3_log_kperp,tk3_log_kperp_1,tk3_log_kperp_h/))) &
         call diagnostic_transfer_3d

  End Subroutine diagnostic_transfer_general

  Subroutine finalize_diag_transfer_general

    if (any((/tky,tkx,tkperp,t_log_kperp,t_log_ky,t_log_kperp_1,t_log_kperp_fb,t_log_kperp_h/))) &
         call finalize_diag_transfer_std
    if (any((/tk3_ky,tk3_kx,tk3_log_kperp,tk3_log_kperp_1,tk3_log_kperp_h/))) &
         call finalize_diag_transfer_3d
    
  end Subroutine finalize_diag_transfer_general



End Module diag_Gyro_LES

