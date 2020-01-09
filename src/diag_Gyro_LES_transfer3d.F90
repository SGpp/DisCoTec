#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! Transfer function
Module diag_Gyro_LES_transfer3d
  
  Use par_mod 
  Use file_io, only: get_unit_nr
  use communications
  use diagnostics_energy, only: get_cfgamma_f,get_cfgamma_h,get_cfgamma_ant,&
       &energy_integral
  use spatial_averages, only: sum_int_3d, sum_int_z
  use prefactors
  USE calc_rhs, only: this_nonlinear_term
  use diag_Gyro_LES_common

  Implicit None

  public:: diagnostic_transfer_3d, initialize_diag_transfer_3d, finalize_diag_transfer_3d

  private
  Character(Len=8) :: filestat='replace', filepos='rewind'
  Integer, dimension(:), allocatable:: transfer_counter_3d
  integer, dimension(:,:,:),allocatable :: FE_TF_SPEC_3D_FILE
  real, dimension(:,:,:,:,:), allocatable :: tf_cumul_3d_1
  real, dimension(:,:,:,:,:), allocatable :: tf_cumul_3d_2
  real, dimension(:,:,:,:,:), allocatable :: tf_cumul_3d_4
  real, dimension(:,:,:,:,:), allocatable :: tf_cumul_3d_6


contains


  Subroutine initialize_diag_transfer_3d

#ifdef with_extended_diags
    allocate(FE_TF_SPEC_3D_FILE(ln1:ln2,3,n_transfer_models))
    allocate(transfer_counter_3d(1:6))
    call init_kperp
    
    IF (tk3_ky) call calc_initialize_diag_transfer_3d(n_shells_ky,1)
    IF (tk3_kx) call calc_initialize_diag_transfer_3d(n_shells_kx,2)
    IF (tk3_log_kperp) call calc_initialize_diag_transfer_3d(n_shells_log_kperp,4)
    IF (tk3_log_kperp_1) call calc_initialize_diag_transfer_3d(n_shells_log_kperp_1,6)
    IF (tk3_log_kperp_h) call calc_initialize_diag_transfer_3d_h(n_shells_log_kperp,4)
#endif

  End Subroutine initialize_diag_transfer_3d

  Subroutine diagnostic_transfer_3d
    
#ifdef with_extended_diags
    IF (tk3_ky) call calc_diagnostic_transfer_3d(n_shells_ky,1)
    IF (tk3_kx) call calc_diagnostic_transfer_3d(n_shells_kx,2)
    IF (tk3_log_kperp) call calc_diagnostic_transfer_3d(n_shells_log_kperp,4)
    IF (tk3_log_kperp_1) call calc_diagnostic_transfer_3d(n_shells_log_kperp_1,6)
    IF (tk3_log_kperp_h) call calc_diagnostic_transfer_3d_h(n_shells_log_kperp,4)
#endif

  End Subroutine diagnostic_transfer_3d

  Subroutine finalize_diag_transfer_3d

#ifdef with_extended_diags
    IF (tk3_ky) call calc_finalize_diag_transfer_3d(1)
    IF (tk3_kx) call calc_finalize_diag_transfer_3d(2)
    IF (tk3_log_kperp) call calc_finalize_diag_transfer_3d(4)
    IF (tk3_log_kperp_1) call calc_finalize_diag_transfer_3d(6)
    IF (tk3_log_kperp_h) call calc_finalize_diag_transfer_3d_h(4)

    if (allocated(kperp)) deallocate(kperp)
    deallocate(FE_TF_SPEC_3D_FILE)
    Deallocate(transfer_counter_3d)
#endif

  End Subroutine finalize_diag_transfer_3d


#ifdef with_extended_diags
!!!******************************************************************!!!
!!!******************************** Locality transfer functions *****
!********************************************************************
Subroutine calc_initialize_diag_transfer_3d(n_shells,which_model)

    Implicit None
    Integer, intent(in) :: n_shells, which_model
    Integer :: n,o 
    Character(len=4) :: xstring
    Character(len=80) :: label, header
    Character(len=40), dimension(3) :: tf_spec_3d_file_label

    tf_spec_3d_file_label = (/'/w_transfer_3d_ ',&
         &'/e_transfer_3d_ ','/fe_transfer_3d_'/)

    select case(which_model)
    case(1)
       allocate(tf_cumul_3d_1(n_shells,n_shells,n_shells,ln1:ln2,3))
       tf_cumul_3d_1 = 0
       transfer_counter_3d(1) = 0
       label ='ky_'
       Write(header,*) '#', kymax
    case(2)
       allocate(tf_cumul_3d_2(n_shells,n_shells,n_shells,ln1:ln2,3))
       tf_cumul_3d_2 = 0
       transfer_counter_3d(2) = 0
       label ='kx_'
       Write(header,*) '#'
    case(4)
       allocate(tf_cumul_3d_4(n_shells,n_shells,n_shells,ln1:ln2,3))
       tf_cumul_3d_4 = 0
       transfer_counter_3d(4) = 0
       label ='log_kperp_'
       Write(header,*) '#', k0, kperp_max
    case(6)
       allocate(tf_cumul_3d_6(n_shells,n_shells,n_shells,ln1:ln2,3))
       tf_cumul_3d_6 = 0
       transfer_counter_3d(6) = 0
       label ='log_kperp1_'
       Write(header,*) '#', lambda1, kperp_max
    case default
       return
    end select

    Write(xstring, '(i4)') n_shells

    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       Do o=1,3
          Do n=ln1,ln2
             call get_unit_nr(FE_TF_SPEC_3D_FILE(n,o,which_model))
             OPEN(FE_TF_SPEC_3D_FILE(n,o,which_model), file=trim(diagdir)//trim(tf_spec_3d_file_label(o))&
                  //trim(adjustl(label))//trim(spec(n)%name)//trim(file_extension),form= 'FORMATTED', &
                  status=filestat, position=filepos)
             Write(FE_TF_SPEC_3D_FILE(n,o,which_model),*) trim(header)
          enddo
       enddo
    endif

End subroutine calc_initialize_diag_transfer_3d

Subroutine calc_initialize_diag_transfer_3d_h(n_shells,which_model)

    Implicit None
    Integer, intent(in) :: n_shells, which_model
    Integer :: n
    Character(len=4) :: xstring
    Character(len=80) :: label, header

    select case(which_model)
    case(1)
       allocate(tf_cumul_3d_1(n_shells,n_shells,n_shells,ln1:ln2,1))
       tf_cumul_3d_1 = 0
       transfer_counter_3d(1) = 0
       label ='ky_'
       Write(header,*) '#', kymax
    case(2)
       allocate(tf_cumul_3d_2(n_shells,n_shells,n_shells,ln1:ln2,1))
       tf_cumul_3d_2 = 0
       transfer_counter_3d(2) = 0
       label ='kx_'
       Write(header,*) '#'
    case(4)
        allocate(tf_cumul_3d_4(n_shells,n_shells,n_shells,ln1:ln2,1))
        tf_cumul_3d_4 = 0
        transfer_counter_3d(4) = 0
        label ='log_kperp_'
        Write(header,*) '#', k0, kperp_max
     case(6)
        allocate(tf_cumul_3d_6(n_shells,n_shells,n_shells,ln1:ln2,1))
        tf_cumul_3d_6 = 0
        transfer_counter_3d(6) = 0
        label ='log_kperp1_'
        Write(header,*) '#', lambda1, kperp_max
     case default
        return
     end select
     
     Write(xstring, '(i4)') n_shells
     
     IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) THEN
        Do n=ln1,ln2
           call get_unit_nr(FE_TF_SPEC_3D_FILE(n,1,which_model))
           OPEN(FE_TF_SPEC_3D_FILE(n,1,which_model), file=trim(diagdir)//'/fe_transfer_3d_'&
                &//trim(adjustl(label))//trim(spec(n)%name)//trim(file_extension),&
                &form= 'FORMATTED', status=filestat, position=filepos)
           Write(FE_TF_SPEC_3D_FILE(n,1,which_model),*) trim(header)
        enddo
    Endif
End Subroutine calc_initialize_diag_transfer_3D_h

Subroutine calc_diagnostic_transfer_3d(n_shells,which_model)

    integer, intent(in) :: n_shells, which_model
    integer  :: ls1, ls2, ls3, o
    complex, dimension(:,:,:),pointer :: ptr2
    integer :: lbg1,lbg2,stage,n    
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    real, dimension(n_shells,n_shells,n_shells,ln1:ln2,3):: fe_tf_tot_3d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2):: e_n_temp
    
    allocate(g_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(g_ls2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(g_ls3(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(nlt_ls2_ls3_g(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(nlt_ls2_ls3_f(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tf_term(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    allocate(f_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(f_ls2(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(f_ls3(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(fields_ls1(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    allocate(fields_ls2(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    allocate(fields_ls3(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1
    lbg2 = lklmn0
    stage = 0

    ! This initialization that should be here and cant be inside the ls1,ls2 loop
    fe_tf_tot_3d = 0.0

    do ls3=1,n_shells
       select case(which_model)
          case(1) 
             call calc_n_shells_ky(ls3,g_1,g_ls3)
          case(2) 
             call calc_n_shells_kx(ls3,g_1,g_ls3)
          case(4) 
             call calc_n_shells_log_kperp(k0,ls3,g_1,g_ls3)
          case(6) 
             call calc_n_shells_log_kperp_1(ls3,g_1,g_ls3)
          case default
          end select
        call calc_aux_fields(g_ls3, fields_ls3, f_ls3,.false.)    
        do ls2=1,n_shells
            if (which_model.eq.1) call calc_n_shells_ky(ls2,g_1,g_ls2)
            if (which_model.eq.2) call calc_n_shells_kx(ls2,g_1,g_ls2)
            if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls2)
            if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls2,g_1,g_ls2)
            call calc_aux_fields(g_ls2, fields_ls2, f_ls2,.false.)    
            nlt_ls2_ls3_g = cmplx(0.0,0.0) 
            nlt_ls2_ls3_f = cmplx(0.0,0.0)

            CALL this_nonlinear_term%add(g_ls2, ptr1, fields_ls3,&
                 &ptr2, ptr3, nlt_ls2_ls3_g, lbg1, lbg2, stage)
            CALL this_nonlinear_term%add(g_ls3, ptr1, fields_ls2,&
                 &ptr2, ptr3, nlt_ls2_ls3_f, lbg1, lbg2, stage)
            do ls1=1,n_shells
                if (which_model.eq.1) call calc_n_shells_ky(ls1,g_1,g_ls1)
                if (which_model.eq.2) call calc_n_shells_kx(ls1,g_1,g_ls1)
                if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
                if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls1,g_1,g_ls1)
                call calc_aux_fields(g_ls1, fields_ls1, f_ls1,.false.)    
                call get_cfgamma(f_ls1,fields_ls1,cfgamma1,cfgamma2)
                fe_tf_term = cmplx(0.0,0.0)
                ! Calculte the tranfers of E, W and FE different for ions than for electrons
                do n =ln1,ln2
                    e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
                    fe_tf_term(:,:,:,:,:,n,1) = cfgamma1(:,:,:,:,:,n)*nlt_ls2_ls3_g(:,:,:,:,:,n)
                    call integral_vpw(fe_tf_term(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
                    if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
                    call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot_3d(ls3,ls2,ls1,n,1))
                enddo
                do n=ln1,ln2
                    e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
                    fe_tf_term(:,:,:,:,:,n,2) = cfgamma2(:,:,:,:,:,n)*nlt_ls2_ls3_f(:,:,:,:,:,n)
                    call integral_vpw(fe_tf_term(:,:,:,:,:,n,2),e_n_temp(:,:,:,n))
                    if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
                    call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot_3d(ls3,ls2,ls1,n,2))
                enddo
                do n=ln1,ln2
                    fe_tf_tot_3d(ls3,ls2,ls1,n,3) = fe_tf_tot_3d(ls3,ls2,ls1,n,1) + fe_tf_tot_3d(ls3,ls2,ls1,n,2)
                enddo
            enddo
        enddo
    enddo        
    
    deallocate(g_ls1) 
    deallocate(g_ls2) 
    deallocate(g_ls3) 
    deallocate(nlt_ls2_ls3_g) 
    deallocate(nlt_ls2_ls3_f) 
    deallocate(cfgamma1) 
    deallocate(cfgamma2) 
    deallocate(fe_tf_term) 
    deallocate(f_ls1) 
    deallocate(f_ls2) 
    deallocate(f_ls3) 
    deallocate(fields_ls1) 
    deallocate(fields_ls2) 
    deallocate(fields_ls3) 
    
    if ( .not.( mod(itime,istep_fe_transfer*avg_window) .eq. 0) ) then
       ! We cumulate in the first case
        if (which_model.eq.1) then
            tf_cumul_3d_1 = fe_tf_tot_3d + tf_cumul_3d_1
            transfer_counter_3d(1) = transfer_counter_3d(1) + 1 
       endif
        if (which_model.eq.2) then
            tf_cumul_3d_2 = fe_tf_tot_3d + tf_cumul_3d_2
            transfer_counter_3d(2) = transfer_counter_3d(2) + 1 
       endif
       if (which_model.eq.4) then
            tf_cumul_3d_4 = fe_tf_tot_3d + tf_cumul_3d_4
            transfer_counter_3d(4) = transfer_counter_3d(4) + 1 
       endif
       if (which_model.eq.6) then
            tf_cumul_3d_6 = fe_tf_tot_3d + tf_cumul_3d_6
            transfer_counter_3d(6) = transfer_counter_3d(6) + 1 
       endif
   endif

    ! Then we do output files
    if ( itime .gt. 0 .and. mod(itime,istep_fe_transfer*avg_window) .eq. 0 ) then
       ! Here we use v2d as a buffer for write routine and average it

       select case(which_model)
       case(1)
          tf_cumul_3d_1 = fe_tf_tot_3d + tf_cumul_3d_1
          transfer_counter_3d(1) = transfer_counter_3d(1) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_1 /real(transfer_counter_3d(1))
       case(2)
          tf_cumul_3d_2 = fe_tf_tot_3d + tf_cumul_3d_2
          transfer_counter_3d(2) = transfer_counter_3d(2) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_2 /real(transfer_counter_3d(2))
       case(4)
          tf_cumul_3d_4 = fe_tf_tot_3d + tf_cumul_3d_4
          transfer_counter_3d(4) = transfer_counter_3d(4) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_4 /real(transfer_counter_3d(4))
       case(6)
          tf_cumul_3d_6 = fe_tf_tot_3d + tf_cumul_3d_6
          transfer_counter_3d(6) = transfer_counter_3d(6) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_6 /real(transfer_counter_3d(6))
       case default
       end select
       
       IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
          Do o=1,3
             Do n=ln1,ln2
                Write(FE_TF_SPEC_3D_FILE(n,o,which_model),*) '#', time
                Do ls3=1,n_shells
                   Do ls2=1,n_shells
                      Write(FE_TF_SPEC_3D_FILE(n,o,which_model),*) 'Transfer spectra T(L1) for shells L2, L3' 
                      Write(FE_TF_SPEC_3D_FILE(n,o,which_model),"(2I4)") ls2, ls3
                      Do ls1=1,n_shells
                         Write(FE_TF_SPEC_3D_FILE(n,o,which_model),"(I4,ES12.4)") ls1, fe_tf_tot_3d(ls3,ls2,ls1,n,o)
                      Enddo
                   Enddo
                Enddo
                call flush(FE_TF_SPEC_3D_FILE(n,o,which_model))
             enddo
          enddo
        Endif
 
        select case(which_model)
        case(1) 
           tf_cumul_3d_1 = 0
           transfer_counter_3d(1) = 0
        case(2) 
           tf_cumul_3d_2 = 0
           transfer_counter_3d(2) = 0
        case(4)
           tf_cumul_3d_4 = 0
           transfer_counter_3d(4) = 0
        case(6)
           tf_cumul_3d_6 = 0
           transfer_counter_3d(6) = 0
        case default
        end select
    endif

  End subroutine calc_diagnostic_transfer_3d



Subroutine calc_diagnostic_transfer_3d_h(n_shells,which_model)

    integer, intent(in) :: n_shells, which_model
    integer  :: ls1, ls2,ls3
    complex, dimension(:,:,:),pointer :: ptr2
    integer :: lbg1,lbg2,stage,n    
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    real, dimension(n_shells,n_shells,n_shells,ln1:ln2,1):: fe_tf_tot_3d
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2):: e_n_temp
    
    allocate(g_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(g_ls2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(g_ls3(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
!    allocate(h_mod_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
!    allocate(h_mod_ls2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
!    allocate(h_mod_ls3(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(nlt_ls2_ls3_h(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(cfgamma(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tf_term(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    allocate(f_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(f_ls2(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(f_ls3(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(h_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(h_ls2(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(h_ls3(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(fields_ls1(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    allocate(fields_ls2(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    allocate(fields_ls3(li1:li2, lj1:lj2, lbz:ubz,n_fields)) 
    
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1
    lbg2 = lklmn0
    stage = 0

    ! This initialization that should be here and cant be inside the ls1,ls2 loop
    fe_tf_tot_3d = 0.0

    do ls3=1,n_shells
        if (which_model.eq.1) call calc_n_shells_ky(ls3,g_1,g_ls3)
        if (which_model.eq.2) call calc_n_shells_kx(ls3,g_1,g_ls3)
        if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls3,g_1,g_ls3)
        if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls3,g_1,g_ls3)
        call calc_aux_fields(g_ls3, fields_ls3, f_ls3,.false.)!comp_h=.true.,h_out=h_ls3)    
!        h_mod_ls3 = h_ls3(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) 
        do ls2=1,n_shells
            if (which_model.eq.1) call calc_n_shells_ky(ls2,g_1,g_ls2)
            if (which_model.eq.2) call calc_n_shells_kx(ls2,g_1,g_ls2)
            if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls2)
            if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls2,g_1,g_ls2)
            call calc_aux_fields(g_ls2, fields_ls2, f_ls2,comp_h=.true.,h_out=h_ls2)    
            nlt_ls2_ls3_h = cmplx(0.0,0.0) 
            if (.not.nonlin_h) then
               allocate(h_mod_ls2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
               h_mod_ls2 = h_ls2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) 
               call this_nonlinear_term%add(h_mod_ls2, ptr1, fields_ls3,&
                    &ptr2, ptr3, nlt_ls2_ls3_h, lbg1, lbg2, stage)
               deallocate(h_mod_ls2) 
            else
               call this_nonlinear_term%add(h_ls2, ptr1, fields_ls3,&
                    &ptr2, ptr3, nlt_ls2_ls3_h, lbg1, lbg2, stage)          
            endif
            do ls1=1,n_shells
                if (which_model.eq.1) call calc_n_shells_ky(ls1,g_1,g_ls1)
                if (which_model.eq.2) call calc_n_shells_kx(ls1,g_1,g_ls1)
                if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
                if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls1,g_1,g_ls1)
                call calc_aux_fields(g_ls1, fields_ls1, f_ls1,comp_h=.true.,h_out=h_ls1)    
                call get_cfgamma_h(h_ls1,cfgamma)
                if (any(antenna_type.eq.(/2,3/))) then
                   allocate(cfgamma_ant(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),&
                        cfgamma_ant_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
                   call get_cfgamma_ant(cfgamma_ant)
                   if (which_model.eq.1) call calc_n_shells_ky(ls1,cfgamma_ant,cfgamma_ant_ls1)
                   if (which_model.eq.2) call calc_n_shells_kx(ls1,cfgamma_ant,cfgamma_ant_ls1)
                   if (which_model.eq.4) call calc_n_shells_log_kperp(k0,ls1,cfgamma_ant,cfgamma_ant_ls1)
                   if (which_model.eq.6) call calc_n_shells_log_kperp_1(ls1,cfgamma_ant,cfgamma_ant_ls1)
                   cfgamma=cfgamma+cfgamma_ant_ls1
                   deallocate(cfgamma_ant,cfgamma_ant_ls1)
                endif
                   
!                h_mod_ls1 = h_ls1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) 

                fe_tf_term = cmplx(0.0,0.0)
                ! Calculate the transfers FE 
                do n =ln1,ln2
                    e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
                    fe_tf_term(:,:,:,:,:,n,1) = cfgamma(:,:,:,:,:,n)*nlt_ls2_ls3_h(:,:,:,:,:,n)
                    call integral_vpw(fe_tf_term(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
                    if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
                    call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot_3d(ls3,ls2,ls1,n,1))
                enddo
            enddo
        enddo
    enddo        
    
    deallocate(g_ls1) 
    deallocate(g_ls2) 
    deallocate(g_ls3) 
!    deallocate(h_mod_ls1) 
!    deallocate(h_mod_ls2) 
!    deallocate(h_mod_ls3) 
    deallocate(nlt_ls2_ls3_h) 
    deallocate(cfgamma) 
    deallocate(fe_tf_term) 
    deallocate(f_ls1) 
    deallocate(f_ls2) 
    deallocate(f_ls3) 
    deallocate(fields_ls1) 
    deallocate(fields_ls2) 
    deallocate(fields_ls3) 
    deallocate(h_ls1) 
    deallocate(h_ls2) 
    deallocate(h_ls3) 
    
    if ( .not.( mod(itime,istep_fe_transfer*avg_window) .eq. 0) ) then
       ! We cumulate in the first case
       select case(which_model)
       case(1)
          tf_cumul_3d_1 = fe_tf_tot_3d + tf_cumul_3d_1
          transfer_counter_3d(1) = transfer_counter_3d(1) + 1 
       case(2)
          tf_cumul_3d_2 = fe_tf_tot_3d + tf_cumul_3d_2
          transfer_counter_3d(2) = transfer_counter_3d(2) + 1 
       case(4)
          tf_cumul_3d_4 = fe_tf_tot_3d + tf_cumul_3d_4
          transfer_counter_3d(4) = transfer_counter_3d(4) + 1 
       case(6)
          tf_cumul_3d_6 = fe_tf_tot_3d + tf_cumul_3d_6
          transfer_counter_3d(6) = transfer_counter_3d(6) + 1 
       case default
       end select
   endif

    ! Then we do output files
    if ( itime .gt. 0 .and. mod(itime,istep_fe_transfer*avg_window) .eq. 0 ) then
       ! Here we use v2d as a buffer for write routine and average it

       select case(which_model)
       case(1)
          tf_cumul_3d_1 = fe_tf_tot_3d + tf_cumul_3d_1
          transfer_counter_3d(1) = transfer_counter_3d(1) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_1 /real(transfer_counter_3d(1))
       case(2)
          tf_cumul_3d_2 = fe_tf_tot_3d + tf_cumul_3d_2
          transfer_counter_3d(2) = transfer_counter_3d(2) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_2 /real(transfer_counter_3d(2))
       case(4)
          tf_cumul_3d_4 = fe_tf_tot_3d + tf_cumul_3d_4
          transfer_counter_3d(4) = transfer_counter_3d(4) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_4 /real(transfer_counter_3d(4))
       case(6)
          tf_cumul_3d_6 = fe_tf_tot_3d + tf_cumul_3d_6
          transfer_counter_3d(6) = transfer_counter_3d(6) + 1 
          fe_tf_tot_3d =  tf_cumul_3d_6 /real(transfer_counter_3d(6))
       case default
       end select
       
        IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
            Do n=ln1,ln2
                Write(FE_TF_SPEC_3D_FILE(n,1,which_model),*) '#', time
                Do ls3=1,n_shells
                    Do ls2=1,n_shells
                        Write(FE_TF_SPEC_3D_FILE(n,1,which_model),*) 'Transfer spectra T(L1) for shells L2, L3' 
                        Write(FE_TF_SPEC_3D_FILE(n,1,which_model),"(2I4)") ls2, ls3
                        Do ls1=1,n_shells
                            Write(FE_TF_SPEC_3D_FILE(n,1,which_model),"(I4,ES12.4)") ls1, fe_tf_tot_3d(ls3,ls2,ls1,n,1)
                        Enddo
                   Enddo
                Enddo
               call flush(FE_TF_SPEC_3D_FILE(n,1,which_model))
           enddo
        Endif
 
       select case(which_model)
       case(1)
          tf_cumul_3d_1 = 0 
          transfer_counter_3d(1) = 0
       case(2)
          tf_cumul_3d_2 = 0 
          transfer_counter_3d(2) = 0
       case(4)
          tf_cumul_3d_4 = 0 
          transfer_counter_3d(4) = 0
       case(6)
          tf_cumul_3d_6 = 0 
          transfer_counter_3d(6) = 0
       case default
       end select
    endif

  End subroutine calc_diagnostic_transfer_3d_h

  Subroutine calc_finalize_diag_transfer_3d(which_model)

    integer,intent(in) :: which_model
    Integer :: n 

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       Do n=ln1,ln2
          close(FE_TF_SPEC_3D_FILE(n,1,which_model))
          close(FE_TF_SPEC_3D_FILE(n,2,which_model))
          close(FE_TF_SPEC_3D_FILE(n,3,which_model))
       Enddo
    Endif

    if (which_model.eq.1) deallocate(tf_cumul_3d_1)
    if (which_model.eq.2) deallocate(tf_cumul_3d_2)
    if (which_model.eq.4) deallocate(tf_cumul_3d_4)
    if (which_model.eq.6) deallocate(tf_cumul_3d_6)

  End subroutine calc_finalize_diag_transfer_3d

  Subroutine calc_finalize_diag_transfer_3d_h(which_model)

    integer,intent(in) :: which_model
    Integer :: n 

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       Do n=ln1,ln2
          close(FE_TF_SPEC_3D_FILE(n,1,which_model))
       Enddo
    Endif

    if (which_model.eq.1) deallocate(tf_cumul_3d_1)
    if (which_model.eq.2) deallocate(tf_cumul_3d_2)
    if (which_model.eq.4) deallocate(tf_cumul_3d_4)
    if (which_model.eq.6) deallocate(tf_cumul_3d_6)

  End subroutine calc_finalize_diag_transfer_3d_h
#endif

end Module diag_Gyro_LES_transfer3d
