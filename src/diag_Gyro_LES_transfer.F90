#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! Transfer function
Module diag_Gyro_LES_transfer
  
  Use par_mod 
  Use file_io, only: get_unit_nr
  use communications
  use collisions, only: equ_collisions
  use diagnostics_energy, only: get_cfgamma_f,get_cfgamma_h,&
       &energy_integral, get_cfgamma_ant
  use dfdzv_terms 
  use dgdxy_terms
  use dzv_terms
  use dfdxy_terms, only: add_dfdxy
  use dchidz_term, only: add_dchidz
  use dchidxy_terms, only: add_dchidxy, add_dchidxy_orig
  use spatial_averages, only: sum_int_3d, sum_int_z
  use prefactors
  use numerical_damping
  USE calc_rhs, only: this_nonlinear_term
  use diag_Gyro_LES_common

  Implicit None

  public:: initialize_diag_transfer_std, finalize_diag_transfer_std, diagnostic_transfer_std

  private
  !****** Parameters for transfer function  
  Character(Len=8) :: filestat='replace', filepos='rewind'
  Real, dimension(:), allocatable :: K_VALUE
  Integer,dimension(:,:,:),allocatable  :: FE_TF_FILE
  Integer,dimension(:,:,:,:),allocatable  :: FE_SP_FILE
  real, dimension(:,:,:,:), allocatable :: fe_tf_tot
  real, dimension(:,:,:,:), allocatable :: fe_tot
  integer, dimension(:), allocatable :: transfer_counter ! Stores the number of windows for average
  real, dimension(:,:,:,:), allocatable :: tf_cumul_1
  real, dimension(:,:,:,:), allocatable :: tf_cumul_2
  real, dimension(:,:,:,:), allocatable :: tf_cumul_3
  real, dimension(:,:,:,:), allocatable :: tf_cumul_4
  real, dimension(:,:,:,:), allocatable :: tf_cumul_5
  real, dimension(:,:,:,:), allocatable :: tf_cumul_6
  real, dimension(:,:,:,:), allocatable :: tf_cumul_7
  real, dimension(:,:,:,:), allocatable :: spec_cumul_1
  real, dimension(:,:,:,:), allocatable :: spec_cumul_2
  real, dimension(:,:,:,:), allocatable :: spec_cumul_3
  real, dimension(:,:,:,:), allocatable :: spec_cumul_4
  real, dimension(:,:,:,:), allocatable :: spec_cumul_5
  real, dimension(:,:,:,:), allocatable :: spec_cumul_6
  real, dimension(:,:,:,:), allocatable :: spec_cumul_7
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_ls_g
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_ls_f



  contains

!!!**************************************************************************!!!
!!!********************* transfer functions analysis ************************!!! 
!!!**************************************************************************!!!

  Subroutine initialize_diag_transfer_std

#ifdef with_extended_diags

    Allocate(FE_TF_FILE(ln1:ln2,n_transfer_models,n_transfer_models)) 
    Allocate(FE_SP_FILE(ln1:ln2,3,8,n_transfer_models)) 
    allocate(transfer_counter(n_transfer_models))

    call init_kperp

    IF (tky) call calc_initialize_diag_transfer(n_shells_ky,1)
    IF (tkx) call calc_initialize_diag_transfer(n_shells_kx,2)
    IF (tkperp) call calc_initialize_diag_transfer(n_shells_kperp,3)
    IF (t_log_kperp) call calc_initialize_diag_transfer(n_shells_log_kperp,4)
    IF (t_log_ky) call calc_initialize_diag_transfer(n_shells_log_ky,5)
    IF (t_log_kperp_1) call calc_initialize_diag_transfer(n_shells_log_kperp_1,6)
    IF (t_log_kperp_fb) call calc_initialize_diag_transfer(n_shells_log_kperp,7)
    IF (t_log_kperp_h) call calc_initialize_diag_transfer_h(n_shells_log_kperp,4)

#endif
   
  End Subroutine initialize_diag_transfer_std


  Subroutine diagnostic_transfer_std

#ifdef with_extended_diags
    IF (.not.with_nvp) THEN
        IF (tky) call calc_diagnostic_transfer(n_shells_ky,1)
        IF (tkx) call calc_diagnostic_transfer(n_shells_kx,2)
        IF (tkperp) call calc_diagnostic_transfer(n_shells_kperp,3)
        IF (t_log_kperp) call calc_diagnostic_transfer(n_shells_log_kperp,4)
        IF (t_log_ky) call calc_diagnostic_transfer(n_shells_log_ky,5)
        IF (t_log_kperp_1) call calc_diagnostic_transfer(n_shells_log_kperp_1,6)
        IF (t_log_kperp_fb) call calc_diagnostic_transfer(n_shells_log_kperp,7)
        IF (t_log_kperp_h) call calc_diagnostic_transfer_h(n_shells_log_kperp,4)
    ELSE 
        IF (tky) call calc_diagnostic_transfer_mod(n_shells_ky,1)
        IF (tkx) call calc_diagnostic_transfer_mod(n_shells_kx,2)
        IF (tkperp) call calc_diagnostic_transfer_mod(n_shells_kperp,3)
        IF (t_log_kperp) call calc_diagnostic_transfer_mod(n_shells_log_kperp,4)
        IF (t_log_ky) call calc_diagnostic_transfer_mod(n_shells_log_ky,5)
        IF (t_log_kperp_1) call calc_diagnostic_transfer_mod(n_shells_log_kperp_1,6)
        IF (t_log_kperp_fb) call calc_diagnostic_transfer_mod(n_shells_log_kperp,7)
    ENDIF
#endif

  End Subroutine diagnostic_transfer_std


  Subroutine finalize_diag_transfer_std

#ifdef with_extended_diags
    IF (tky) call calc_finalize_diag_transfer(1)
    IF (tkx) call calc_finalize_diag_transfer(2)
    IF (tkperp) call calc_finalize_diag_transfer(3)
    IF (t_log_kperp) call calc_finalize_diag_transfer(4)
    IF (t_log_ky) call calc_finalize_diag_transfer(5)
    IF (t_log_kperp_1) call calc_finalize_diag_transfer(6)
    IF (t_log_kperp_fb) call calc_finalize_diag_transfer(7)
    IF (t_log_kperp_h) call calc_finalize_diag_transfer_h(4)
    
    deallocate(kperp)
    Deallocate(FE_TF_FILE)
    Deallocate(FE_SP_FILE)
    Deallocate(transfer_counter)

#endif

  End Subroutine finalize_diag_transfer_std

!!!! ************************************************************************ !!!
#ifdef with_extended_diags
  Subroutine calc_initialize_diag_transfer(n_shells,which_model)

    Implicit None
    Integer, intent(in) :: n_shells, which_model
    Integer :: n,o,t
    Character(len=4) :: xstring
    Character(len=20) :: label
    Character(len=20), dimension(8) :: label1
    Character(len=20), dimension(3) :: tf_file_label, spectral_file_label
    Write(xstring, '(i4)') n_shells
    
    label1 = (/'fe       ','drive    ','hyp_z    ',&
              &'hyp_v    ','hyp_kperp','coll     ',&
              &'curvature','parallel '/)

    tf_file_label = (/'/w_transfer_ ','/e_transfer_ ','/fe_transfer_'/)
    spectral_file_label = (/'/Spectral_W_ ','/Spectral_E_ ','/Spectral_FE_'/)

    select case(which_model)
    case(1)
       allocate(tf_cumul_1(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_1(n_shells,ln1:ln2,3,8))
       tf_cumul_1 = 0
       spec_cumul_1 = 0 
       label ='ky_'
    case(2)
       allocate(tf_cumul_2(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_2(n_shells,ln1:ln2,3,8))
       tf_cumul_2 = 0
       spec_cumul_2 = 0 
       label ='kx_'
    case(3)
       allocate(tf_cumul_3(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_3(n_shells,ln1:ln2,3,8))
       tf_cumul_3 = 0 
       spec_cumul_3 = 0 
       label ='kperp_'
    case(4)
       allocate(tf_cumul_4(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_4(n_shells,ln1:ln2,3,8))
       tf_cumul_4 = 0 
       spec_cumul_4 = 0 
       label ='log_kperp_'
    case(5)
       allocate(tf_cumul_5(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_5(n_shells,ln1:ln2,3,8))
       tf_cumul_5 = 0 
       spec_cumul_5 = 0 
       label ='log_ky_'
    case(6)
       allocate(tf_cumul_6(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_6(n_shells,ln1:ln2,3,8))
       tf_cumul_6 = 0 
       spec_cumul_6 = 0 
       label ='log_kperp1_'
     case(7)
       allocate(tf_cumul_7(n_shells,n_shells,ln1:ln2,n_transfer_models))
       allocate(spec_cumul_7(n_shells,ln1:ln2,3,8))
       tf_cumul_7 = 0 
       spec_cumul_7 = 0 
       label ='log_kperp_fb_'
   case default
   end select
    transfer_counter(which_model)=0

    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
    IF (.not.with_nvp) THEN
       Do o=1,3
          Do n=ln1,ln2
             call get_unit_nr(FE_TF_FILE(n,o,which_model))
             OPEN(FE_TF_FILE(n,o,which_model), file=trim(diagdir)//trim(tf_file_label(o))&
                  //trim(adjustl(label))//trim(spec(n)%name)//trim(file_extension),form= 'FORMATTED', &
                  status=filestat, position=filepos)
          enddo
       enddo
    ENDIF
       do t=1,8
          do o=1,3
             do n=ln1,ln2
                call get_unit_nr(FE_SP_FILE(n,o,t,which_model))
                OPEN(FE_SP_FILE(n,o,t,which_model), file=trim(diagdir)//trim(spectral_file_label(o))&
                     //trim(adjustl(label1(t)))//'_'//trim(adjustl(label))//&
                     trim(spec(n)%name)//trim(file_extension),form = 'FORMATTED', &
                     status=filestat, position=filepos)         
             enddo
          enddo
       enddo
    endif

  End Subroutine calc_initialize_diag_transfer

  Subroutine calc_initialize_diag_transfer_h(n_shells,which_model)

    Implicit None
    Integer, intent(in) :: n_shells, which_model
    Integer :: n,o
    Character(len=4) :: xstring
    Character(len=20) :: label
    Character(len=20), dimension(8) :: label1
    Character(len=20), dimension(3) :: tf_file_label
    Write(xstring, '(i4)') n_shells
    
    label1 = (/'fe       ','drive    ','hyp_z    ',&
              &'hyp_v    ','hyp_kperp','coll     ',&
              &'curvature','parallel '/)

    tf_file_label = (/'/w_transfer_ ','/e_transfer_ ','/fe_transfer_'/)

    select case(which_model)
    case(1)
       allocate(tf_cumul_1(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_1 = 0
       label ='ky_'
    case(2)
       allocate(tf_cumul_2(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_2 = 0
       label ='kx_'
    case(3)
       allocate(tf_cumul_3(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_3 = 0 
       label ='kperp_'
    case(4)
       allocate(tf_cumul_4(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_4 = 0 
       label ='log_kperp_'
    case(5)
       allocate(tf_cumul_5(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_5 = 0 
       label ='log_ky_'
    case(6)
       allocate(tf_cumul_6(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_6 = 0 
       label ='log_kperp1_'
     case(7)
       allocate(tf_cumul_7(n_shells,n_shells,ln1:ln2,n_transfer_models))
       tf_cumul_7 = 0 
       label ='log_kperp_fb_'
   case default
   end select
   transfer_counter(which_model)=0
   
   IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
      IF (.not.with_nvp) THEN
         o=1
         Do n=ln1,ln2
            call get_unit_nr(FE_TF_FILE(n,o,which_model))
            OPEN(FE_TF_FILE(n,o,which_model), file=trim(diagdir)//trim(tf_file_label(o))&
                 //trim(adjustl(label))//trim(spec(n)%name)//trim(file_extension),form= 'FORMATTED', &
                 status=filestat, position=filepos)
         enddo
      ENDIF
   endif
   
  End Subroutine calc_initialize_diag_transfer_h


!!!*****************************************************************************************!!!

  Subroutine calc_diagnostic_transfer(n_shells,which_model)


    integer, intent(in) :: which_model,n_shells
    !g,f,h and nonlinearities and fields arrays

    !Other factors used
    complex, dimension(:,:,:),pointer :: ptr2
    integer :: lbg1,lbg2,stage, ls1,ls2    
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    !transfer array and flux

    allocate(g_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(nlt_ls_g(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(nlt_ls_f(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(f_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(fields_ls1(li1:li2, lj1:lj2, lbz:ubz,1:n_fields)) 
    allocate(fe_tf_tot(n_shells,n_shells,ln1:ln2,n_transfer_models))
    allocate(fe_tot(n_shells,ln1:ln2,3,8))
    allocate(k_value(1:n_shells))

    !abn 8mar2010  some initialization values
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1
    lbg2 = lklmn0
    stage = 0

    ! This initialization that should be here and cant be inside the ls1,ls2 loop
    fe_tf_tot = 0.0
    fe_tot = 0.0

    do ls2=1,n_shells
       select case(which_model)
       case(1) ! ky shells           
          call calc_n_shells_ky(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_ky 
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:))
       case(2) ! kx shells
          call calc_n_shells_kx(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_kx
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
       case(3) ! kperp shells
          call calc_n_shells_kperp(ls2,g_1,g_ls1)
          k_value(ls2) = ls2*delta_kperp
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
       case(4) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
      case(5) ! log ky shells
          call calc_n_shells_log_ky(k0_y,ls2,g_1,g_ls1)
          k_value(ls2) = k0_y*2**(real((ls2-1)/6.))
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
      case(6) ! perpendicular shells modification 
          call calc_n_shells_log_kperp_1(ls2,g_1,g_ls1)
          If (ls2.eq.1) k_value(ls2) = lambda1
          If (ls2.eq.2) k_value(ls2) = lambda1 +lambda2
          IF (ls2.ge.3 .and. ls2.lt.n_shells) k_value(ls2) = &
            (lambda1+lambda2)*2**(real((ls2-2)/5.))
          If (ls2.eq.n_shells) k_value(ls2) = kperp_max
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
       case(7) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
          call calc_spectral(g_ls1,fe_tot(ls2,:,:,:)) 
       end select

       call calc_aux_fields(g_ls1, fields_ls1, f_ls1,.false.)    
       nlt_ls_g = cmplx(0.0,0.0) 
       nlt_ls_f = cmplx(0.0,0.0)
       CALL this_nonlinear_term%add(g_ls1, ptr1, emfields,&
            &ptr2, ptr3, nlt_ls_g, lbg1, lbg2, stage)

       CALL this_nonlinear_term%add(g_1, ptr1, fields_ls1,&
            &ptr2, ptr3, nlt_ls_f, lbg1, lbg2, stage)

       do ls1=1,n_shells
          select case(which_model)
          case(1) ! ky shells
             call calc_n_shells_ky(ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(2) ! kx shells
             call calc_n_shells_kx(ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(3) ! kperp1 shells
             call calc_n_shells_kperp(ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(4) ! log kperp shells
             call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(5) ! log ky shells
             call calc_n_shells_log_ky(k0_y,ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(6) ! log ky shells
             call calc_n_shells_log_kperp_1(ls1,g_1,g_ls1)
             call calc_transfer(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
          case(7) ! log kperp shells for back scattering
             call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
             call calc_transfer_fb(ls2,ls1,g_ls1,nlt_ls_f,nlt_ls_g,fe_tf_tot(ls2,ls1,:,:))
         end select
       enddo
    enddo
    
    deallocate(g_ls1) 
    deallocate(nlt_ls_g) 
    deallocate(nlt_ls_f) 
    deallocate(f_ls1) 
    deallocate(fields_ls1) 

    if ( .not.( mod(itime,istep_fe_transfer*avg_window) .eq. 0) ) then
       ! We cumulate in the first case
       select case(which_model)
       case(1)
          tf_cumul_1 = fe_tf_tot + tf_cumul_1
          spec_cumul_1 = fe_tot + spec_cumul_1
       case(2)
          tf_cumul_2 = fe_tf_tot + tf_cumul_2
          spec_cumul_2 = fe_tot + spec_cumul_2
       case(3)
          tf_cumul_3 = fe_tf_tot + tf_cumul_3
          spec_cumul_3 = fe_tot + spec_cumul_3
       case(4)
          tf_cumul_4 = fe_tf_tot + tf_cumul_4
          spec_cumul_4 = fe_tot + spec_cumul_4
       case(5)
          tf_cumul_5 = fe_tf_tot + tf_cumul_5
          spec_cumul_5 = fe_tot + spec_cumul_5
       case(6)
          tf_cumul_6 = fe_tf_tot + tf_cumul_6
          spec_cumul_6 = fe_tot + spec_cumul_6
       case(7)
          tf_cumul_7 = fe_tf_tot + tf_cumul_7
          spec_cumul_7 = fe_tot + spec_cumul_7
       case default
       end select

       transfer_counter(which_model) = transfer_counter(which_model) + 1 
    endif

    ! Then we do output files
    if ( itime .gt. 0 .and. mod(itime,istep_fe_transfer*avg_window) .eq. 0 ) then
       ! Here we use v2d as a buffer for write routine and average it

       select case(which_model)
       case(1)
          tf_cumul_1 = fe_tf_tot + tf_cumul_1
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_1 /real(transfer_counter(which_model))

          spec_cumul_1 = fe_tot + spec_cumul_1
          fe_tot =  spec_cumul_1/real(transfer_counter(which_model))
       case(2)
          tf_cumul_2 = fe_tf_tot + tf_cumul_2
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_2/real(transfer_counter(which_model))
       
          spec_cumul_2 = fe_tot + spec_cumul_2
          fe_tot =  spec_cumul_2/real(transfer_counter(which_model))
       case(3)
          tf_cumul_3 = fe_tf_tot + tf_cumul_3
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_3/real(transfer_counter(which_model))
       
          spec_cumul_3 = fe_tot + spec_cumul_3
          fe_tot =  spec_cumul_3/real(transfer_counter(which_model))
       case(4)
          tf_cumul_4 = fe_tf_tot + tf_cumul_4
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_4/real(transfer_counter(which_model))
          
          spec_cumul_4 = fe_tot + spec_cumul_4
          fe_tot =  spec_cumul_4/real(transfer_counter(which_model))
       case(5)
          tf_cumul_5 = fe_tf_tot + tf_cumul_5
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_5/real(transfer_counter(which_model))
      
          spec_cumul_5 = fe_tot + spec_cumul_5
          fe_tot =  spec_cumul_5/real(transfer_counter(which_model))
       case(6)
          tf_cumul_6 = fe_tf_tot + tf_cumul_6
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_6/real(transfer_counter(which_model))
      
          spec_cumul_6 = fe_tot + spec_cumul_6
          fe_tot =  spec_cumul_6/real(transfer_counter(which_model))
        case(7)
          tf_cumul_7 = fe_tf_tot + tf_cumul_7
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_7/real(transfer_counter(which_model))
          
          spec_cumul_7 = fe_tot + spec_cumul_7
          fe_tot =  spec_cumul_7/real(transfer_counter(which_model))
       case default
       end select

       Call write_transfer(n_shells,which_model)    

       select case(which_model)
       case(1)
           tf_cumul_1 = 0
           Call write_spectral(n_shells,which_model)
           spec_cumul_1 = 0 
        case(2)
           tf_cumul_2 = 0
           Call write_spectral(n_shells,which_model)
           spec_cumul_2 = 0 
        case(3)
           tf_cumul_3 = 0 
           Call write_spectral(n_shells,which_model)
           spec_cumul_3 = 0 
        case(4)
           tf_cumul_4 = 0 
           Call write_spectral(n_shells,which_model)
           spec_cumul_4 = 0 
        case(5)
           tf_cumul_5 = 0 
           Call write_spectral(n_shells,which_model)
           spec_cumul_5 = 0 
        case(6)
           tf_cumul_6 = 0 
           Call write_spectral(n_shells,which_model)
           spec_cumul_6 = 0 
        case(7)
           tf_cumul_7 = 0 
           Call write_spectral(n_shells,which_model)
           spec_cumul_7 = 0 
        case default
        end select
 
       transfer_counter(which_model)=0
    endif
    deallocate(fe_tf_tot,fe_tot,k_value) 

  End subroutine calc_diagnostic_transfer

  Subroutine calc_diagnostic_transfer_h(n_shells,which_model)

    integer, intent(in) :: which_model,n_shells
    !g,f,h and nonlinearities and fields arrays

    !Other factors used
    complex, dimension(:,:,:),pointer :: ptr2
    integer :: lbg1,lbg2,stage, ls1,ls2,n    
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2):: e_n_temp
    
   !transfer array and flux
    allocate(g_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(h_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(cfgamma(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(f_ls1(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    allocate(fields_ls1(li1:li2, lj1:lj2, lbz:ubz,1:n_fields)) 
    allocate(nlt_ls_g(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tf_tot(n_shells,n_shells,ln1:ln2,n_transfer_models))
    allocate(fe_tot(n_shells,ln1:ln2,3,8))
    allocate(k_value(1:n_shells))
    allocate(fe_tf_term(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    
    !abn 8mar2010  some initialization values
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1
    lbg2 = lklmn0
    stage = 0

    ! This initialization that should be here and cant be inside the ls1,ls2 loop
    fe_tf_tot = 0.0
    fe_tot = 0.0

    !
    do ls2=1,n_shells
       select case(which_model)
       case(1) ! ky shells           
          call calc_n_shells_ky(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_ky 
       case(2) ! kx shells
          call calc_n_shells_kx(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_kx
       case(3) ! kperp shells
          call calc_n_shells_kperp(ls2,g_1,g_ls1)
          k_value(ls2) = ls2*delta_kperp
       case(4) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
      case(5) ! log ky shells
          call calc_n_shells_log_ky(k0_y,ls2,g_1,g_ls1)
          k_value(ls2) = k0_y*2**(real((ls2-1)/6.))
      case(6) ! perpendicular shells modification 
          call calc_n_shells_log_kperp_1(ls2,g_1,g_ls1)
          If (ls2.eq.1) k_value(ls2) = lambda1
          If (ls2.eq.2) k_value(ls2) = lambda1 +lambda2
          IF (ls2.ge.3 .and. ls2.lt.n_shells) k_value(ls2) = &
            (lambda1+lambda2)*2**(real((ls2-2)/5.))
          If (ls2.eq.n_shells) k_value(ls2) = kperp_max
       case(7) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
       !   call calc_spectral_h(g_ls2,fe_tot(ls2,:,:,:)) 
       end select
       call calc_aux_fields(g_ls1, fields_ls1, f_ls1,comp_h=.true.,h_out=h_ls1)    
     
       nlt_ls_g = cmplx(0.0,0.0) 

       !if nonlin_h, the nonlinearity accepts h-shaped arrays, otherwise an intermediate
       !step is needed
       if (.not.nonlin_h) then
          allocate(h_mod_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
          h_mod_ls1 = h_ls1(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2) 
          call this_nonlinear_term%add(h_mod_ls1, ptr1, emfields,&
               &ptr2, ptr3, nlt_ls_g, lbg1, lbg2, stage)
          deallocate(h_mod_ls1) 
       else
          call this_nonlinear_term%add(h_ls1, ptr1, emfields,&
               &ptr2, ptr3, nlt_ls_g, lbg1, lbg2, stage)          
       endif

       !from here, we can re-use the arrays f_ls1, g_ls1, and h_ls1 for the second shell
       !only nlt_ls_g has to be kept from the outer loop
       do ls1=1,n_shells
          if (any(antenna_type.eq.(/2,3/))) then
             allocate(cfgamma_ant(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),&
                  cfgamma_ant_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
             call get_cfgamma_ant(cfgamma_ant)
          endif
          select case(which_model)
          case(1) ! ky shells
             call calc_n_shells_ky(ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_ky(ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(2) ! kx shells
             call calc_n_shells_kx(ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_kx(ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(3) ! kperp1 shells
             call calc_n_shells_kperp(ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_kperp(ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(4) ! log kperp shells
             call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_log_kperp(k0,ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(5) ! log ky shells
             call calc_n_shells_log_ky(k0_y,ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_log_ky(k0_y,ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(6) ! log ky shells
             call calc_n_shells_log_kperp_1(ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_log_kperp_1(ls1,cfgamma_ant,cfgamma_ant_ls1)
          case(7) ! log kperp shells for back scattering
             call calc_n_shells_log_kperp(k0,ls1,g_1,g_ls1)
             if (any(antenna_type.eq.(/2,3/))) call calc_n_shells_log_kperp(k0,ls1,cfgamma_ant,cfgamma_ant_ls1)
         end select
        
         call calc_aux_fields(g_ls1, fields_ls1, f_ls1,comp_h=.true.,h_out=h_ls1)    
         call get_cfgamma_h(h_ls1,cfgamma)
         if (any(antenna_type.eq.(/2,3/))) then
            cfgamma=cfgamma+cfgamma_ant_ls1
            deallocate(cfgamma_ant,cfgamma_ant_ls1)
         endif

         fe_tf_term = cmplx(0.0,0.0)
         ! Calculate the transfers FE 
         do n =ln1,ln2
             e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
             fe_tf_term(:,:,:,:,:,n,1) = cfgamma(:,:,:,:,:,n)*nlt_ls_g(:,:,:,:,:,n)
             call integral_vpw(fe_tf_term(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot(ls2,ls1,n,1))
         enddo

       enddo
    enddo
   
    deallocate(fe_tf_term) 
    deallocate(g_ls1) 
    deallocate(h_ls1) 
    deallocate(f_ls1) 
    deallocate(fields_ls1)
    deallocate(nlt_ls_g) 
    deallocate(cfgamma) 

    if ( .not.( mod(itime,istep_fe_transfer*avg_window) .eq. 0) ) then
       ! We cumulate in the first case
       select case(which_model)
       case(1)
          tf_cumul_1 = fe_tf_tot + tf_cumul_1
          !spec_cumul_1 = fe_tot + spec_cumul_1
       case(2)
          tf_cumul_2 = fe_tf_tot + tf_cumul_2
          !spec_cumul_2 = fe_tot + spec_cumul_2
       case(3)
          tf_cumul_3 = fe_tf_tot + tf_cumul_3
          !spec_cumul_3 = fe_tot + spec_cumul_3
       case(4)
          tf_cumul_4 = fe_tf_tot + tf_cumul_4
          !spec_cumul_4 = fe_tot + spec_cumul_4
       case(5)
          tf_cumul_5 = fe_tf_tot + tf_cumul_5
          !spec_cumul_5 = fe_tot + spec_cumul_5
       case(6)
          tf_cumul_6 = fe_tf_tot + tf_cumul_6
          !spec_cumul_6 = fe_tot + spec_cumul_6
       case(7)
          tf_cumul_7 = fe_tf_tot + tf_cumul_7
          !spec_cumul_7 = fe_tot + spec_cumul_7
       case default
       end select

       transfer_counter(which_model) = transfer_counter(which_model) + 1 
    endif

    ! Then we do output files
    if ( itime .gt. 0 .and. mod(itime,istep_fe_transfer*avg_window) .eq. 0 ) then
       ! Here we use v2d as a buffer for write routine and average it

       select case(which_model)
       case(1)
          tf_cumul_1 = fe_tf_tot + tf_cumul_1
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_1 /real(transfer_counter(which_model))

          !spec_cumul_1 = fe_tot + spec_cumul_1
          !fe_tot =  spec_cumul_1/real(transfer_counter(which_model))
       case(2)
          tf_cumul_2 = fe_tf_tot + tf_cumul_2
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_2/real(transfer_counter(which_model))
       
          !spec_cumul_2 = fe_tot + spec_cumul_2
          !fe_tot =  spec_cumul_2/real(transfer_counter(which_model))
       case(3)
          tf_cumul_3 = fe_tf_tot + tf_cumul_3
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_3/real(transfer_counter(which_model))
       
          !spec_cumul_3 = fe_tot + spec_cumul_3
          !fe_tot =  spec_cumul_3/real(transfer_counter(which_model))
       case(4)
          tf_cumul_4 = fe_tf_tot + tf_cumul_4
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_4/real(transfer_counter(which_model))
          
          !spec_cumul_4 = fe_tot + spec_cumul_4
          !fe_tot =  spec_cumul_4/real(transfer_counter(which_model))
       case(5)
          tf_cumul_5 = fe_tf_tot + tf_cumul_5
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_5/real(transfer_counter(which_model))
      
          !spec_cumul_5 = fe_tot + spec_cumul_5
          !fe_tot =  spec_cumul_5/real(transfer_counter(which_model))
       case(6)
          tf_cumul_6 = fe_tf_tot + tf_cumul_6
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_6/real(transfer_counter(which_model))
      
          !spec_cumul_6 = fe_tot + spec_cumul_6
          !fe_tot =  spec_cumul_6/real(transfer_counter(which_model))
        case(7)
          tf_cumul_7 = fe_tf_tot + tf_cumul_7
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          fe_tf_tot =  tf_cumul_7/real(transfer_counter(which_model))
          
          !spec_cumul_7 = fe_tot + spec_cumul_7
          !fe_tot =  spec_cumul_7/real(transfer_counter(which_model))
       case default
       end select

       Call write_transfer_h(n_shells,which_model)    

       select case(which_model)
       case(1)
           tf_cumul_1 = 0
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_1 = 0 
        case(2)
           tf_cumul_2 = 0
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_2 = 0 
        case(3)
           tf_cumul_3 = 0 
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_3 = 0 
        case(4)
           tf_cumul_4 = 0 
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_4 = 0 
        case(5)
           tf_cumul_5 = 0 
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_5 = 0 
        case(6)
          tf_cumul_6 = 0 
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_6 = 0 
        case(7)
           tf_cumul_7 = 0 
          ! Call write_spectral(n_shells,which_model)
          ! spec_cumul_7 = 0 
        case default
        end select
 
       transfer_counter(which_model)=0
    endif
    deallocate(fe_tf_tot,fe_tot,k_value) 

  End subroutine calc_diagnostic_transfer_h

!!!*****************************************************************************************!!!

  Subroutine calc_diagnostic_transfer_mod(n_shells,which_model)


    integer, intent(in) :: which_model,n_shells
    integer :: ls2    
    allocate(g_ls1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tot(n_shells,ln1:ln2,3,8))
    allocate(k_value(1:n_shells))

    ! This initialization that should be here and cant be inside the ls1,ls2 loop
    fe_tot = 0.0

    do ls2=1,n_shells
       select case(which_model)
       case(1) ! ky shells           
          call calc_n_shells_ky(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_ky 
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
       case(2) ! kx shells
          call calc_n_shells_kx(ls2,g_1,g_ls1)
          k_value(ls2) = (ls2-1)*delta_kx
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
       case(3) ! kperp shells
          call calc_n_shells_kperp(ls2,g_1,g_ls1)
          k_value(ls2) = ls2*delta_kperp
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
       case(4) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
      case(5) ! log ky shells
          call calc_n_shells_log_ky(k0_y,ls2,g_1,g_ls1)
          k_value(ls2) = k0_y*2**(real((ls2-1)/6.))
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
      case(6) ! perpendicular shells modification 
          call calc_n_shells_log_kperp_1(ls2,g_1,g_ls1)
          If (ls2.eq.1) k_value(ls2) = lambda1
          If (ls2.eq.2) k_value(ls2) = lambda1 +lambda2
          IF (ls2.ge.3 .and. ls2.lt.n_shells) k_value(ls2) = &
            (lambda1+lambda2)*2**(real((ls2-2)/5.))
          If (ls2.eq.n_shells) k_value(ls2) = kperp_max
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
       case(7) ! log kperp shells
          call calc_n_shells_log_kperp(k0,ls2,g_1,g_ls1)
          k_value(ls2) = k0*2**(real((ls2-1)/5.))
          call calc_spectral_mod(g_ls1,fe_tot(ls2,:,:,:)) 
       end select
    enddo
    
    deallocate(g_ls1) 

    if ( .not.( mod(itime,istep_fe_transfer*avg_window) .eq. 0) ) then
       ! We cumulate in the first case
       select case(which_model)
       case(1)
          spec_cumul_1 = fe_tot + spec_cumul_1
       case(2)
          spec_cumul_2 = fe_tot + spec_cumul_2
       case(3)
          spec_cumul_3 = fe_tot + spec_cumul_3
       case(4)
          spec_cumul_4 = fe_tot + spec_cumul_4
       case(5)
          spec_cumul_5 = fe_tot + spec_cumul_5
       case(6)
          spec_cumul_6 = fe_tot + spec_cumul_6
       case(7)
          spec_cumul_7 = fe_tot + spec_cumul_7
       case default
       end select

       transfer_counter(which_model) = transfer_counter(which_model) + 1 
    endif

    ! Then we do output files
    if ( itime .gt. 0 .and. mod(itime,istep_fe_transfer*avg_window) .eq. 0 ) then
       ! Here we use v2d as a buffer for write routine and average it

       select case(which_model)
       case(1)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_1 = fe_tot + spec_cumul_1
          fe_tot =  spec_cumul_1/real(transfer_counter(which_model))
       case(2)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_2 = fe_tot + spec_cumul_2
          fe_tot =  spec_cumul_2/real(transfer_counter(which_model))
       case(3)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_3 = fe_tot + spec_cumul_3
          fe_tot =  spec_cumul_3/real(transfer_counter(which_model))
       case(4)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_4 = fe_tot + spec_cumul_4
          fe_tot =  spec_cumul_4/real(transfer_counter(which_model))
       case(5)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_5 = fe_tot + spec_cumul_5
          fe_tot =  spec_cumul_5/real(transfer_counter(which_model))
       case(6)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_6 = fe_tot + spec_cumul_6
          fe_tot =  spec_cumul_6/real(transfer_counter(which_model))
        case(7)
          transfer_counter(which_model) = transfer_counter(which_model) + 1 
          spec_cumul_7 = fe_tot + spec_cumul_7
          fe_tot =  spec_cumul_7/real(transfer_counter(which_model))
       case default
       end select

       select case(which_model)
       case(1)
           Call write_spectral(n_shells,which_model)
           spec_cumul_1 = 0 
        case(2)
           Call write_spectral(n_shells,which_model)
           spec_cumul_2 = 0 
        case(3)
           Call write_spectral(n_shells,which_model)
           spec_cumul_3 = 0 
        case(4)
           Call write_spectral(n_shells,which_model)
           spec_cumul_4 = 0 
        case(5)
           Call write_spectral(n_shells,which_model)
           spec_cumul_5 = 0 
        case(6)
           Call write_spectral(n_shells,which_model)
           spec_cumul_6 = 0 
        case(7)
           Call write_spectral(n_shells,which_model)
           spec_cumul_7 = 0 
        case default
        end select
 
       transfer_counter(which_model)=0
    endif
    deallocate(fe_tot,k_value) 

  End subroutine calc_diagnostic_transfer_mod

!!!*****************************************************************************************!!!

  Subroutine calc_transfer(ls2,ls1,g_ls1,nlt_ls2_f,nlt_ls2_g,fe_tf_tot)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) ::g_ls1, nlt_ls2_f, nlt_ls2_g
    integer, intent(in) :: ls1, ls2
    real, dimension(ln1:ln2,3),intent(inout) :: fe_tf_tot
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1: n_fields) :: fields_ls1
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2) :: f_ls1
    !temporary array to integrate over vp , mu and species
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2)::e_n_temp
    !Arrays of terms to be analyzed for transfer analysis
    integer :: n
    
    allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tf_term(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    
    call calc_aux_fields(g_ls1, fields_ls1, f_ls1,.false.)    
    call get_cfgamma(f_ls1,fields_ls1,cfgamma1,cfgamma2)

    ! Calculte the tranfers of E, W and FE different for ions than for electrons

    do n =ln1,ln2
       e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
       fe_tf_term(:,:,:,:,:,n,1) = cfgamma1(:,:,:,:,:,n)*nlt_ls2_g(:,:,:,:,:,n)
       call integral_vpw(fe_tf_term(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
       if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
       call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot(n,1))
    enddo

    do n=ln1,ln2
       e_n_temp(:,:,:,n) = cmplx(0.0,0.0)
       fe_tf_term(:,:,:,:,:,n,2) = cfgamma2(:,:,:,:,:,n)*nlt_ls2_f(:,:,:,:,:,n)
       call integral_vpw(fe_tf_term(:,:,:,:,:,n,2),e_n_temp(:,:,:,n))
       if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
       call sum_int_3d(e_n_temp(:,:,:,n),fe_tf_tot(n,2))
    enddo


    do n=ln1,ln2
       fe_tf_tot(n,3) = fe_tf_tot(n,1) + fe_tf_tot(n,2)
    enddo

    deallocate(cfgamma1) 
    deallocate(cfgamma2) 
    deallocate(fe_tf_term) 
  

  End subroutine calc_transfer

  Subroutine calc_transfer_fb(ls2,ls1,g_ls1,nlt_ls2_f,nlt_ls2_g,fe_tf_tot)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_ls1, nlt_ls2_f, nlt_ls2_g
    integer, intent(in) :: ls1, ls2
    real, dimension(ln1:ln2,6),intent(inout) :: fe_tf_tot
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1: n_fields) :: fields_ls1
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2) :: f_ls1
    !temporary array to integrate over vp , mu and species
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2,2)::e_n_temp
    !Arrays of terms to be analyzed for transfer analysis
    integer :: n
    
    allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
    allocate(fe_tf_term(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    allocate(fe_tf_term_p(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    allocate(fe_tf_term_n(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2,2)) 
    
    call calc_aux_fields(g_ls1, fields_ls1, f_ls1,.false.)    
    call get_cfgamma(f_,emfields,cfgamma1,cfgamma2)

    ! Calculte the tranfers of E, W and FE different for ions than for electrons

    do n =ln1,ln2
       e_n_temp(:,:,:,n,:) = cmplx(0.0,0.0)
       fe_tf_term(:,:,:,:,:,n,1) = cfgamma1(:,:,:,:,:,n)*nlt_ls2_g(:,:,:,:,:,n)
       call calc_n_shells_log_kperp_pn(k0,ls1,fe_tf_term(:,:,:,:,:,n,1),fe_tf_term_p(:,:,:,:,:,n,1),&
                                        fe_tf_term_n(:,:,:,:,:,n,1))
       call integral_vpw(fe_tf_term_p(:,:,:,:,:,n,1),e_n_temp(:,:,:,n,1))
       call integral_vpw(fe_tf_term_n(:,:,:,:,:,n,1),e_n_temp(:,:,:,n,2))
       if(p_has_00_mode) e_n_temp(0,0,:,n,1) =cmplx(0.0,0.0)
       if(p_has_00_mode) e_n_temp(0,0,:,n,2) =cmplx(0.0,0.0)
       call sum_int_3d(e_n_temp(:,:,:,n,1),fe_tf_tot(n,1))
       call sum_int_3d(e_n_temp(:,:,:,n,2),fe_tf_tot(n,4))
    enddo

    do n =ln1,ln2
       e_n_temp(:,:,:,n,:) = cmplx(0.0,0.0)
       fe_tf_term(:,:,:,:,:,n,2) = cfgamma2(:,:,:,:,:,n)*nlt_ls2_f(:,:,:,:,:,n)
       call calc_n_shells_log_kperp_pn(k0,ls1,fe_tf_term(:,:,:,:,:,n,2),fe_tf_term_p(:,:,:,:,:,n,2),&
                                        fe_tf_term_n(:,:,:,:,:,n,2))
       call integral_vpw(fe_tf_term_p(:,:,:,:,:,n,2),e_n_temp(:,:,:,n,1))
       call integral_vpw(fe_tf_term_n(:,:,:,:,:,n,2),e_n_temp(:,:,:,n,2))
       if(p_has_00_mode) e_n_temp(0,0,:,n,1) =cmplx(0.0,0.0)
       if(p_has_00_mode) e_n_temp(0,0,:,n,2) =cmplx(0.0,0.0)
       call sum_int_3d(e_n_temp(:,:,:,n,1),fe_tf_tot(n,2))
       call sum_int_3d(e_n_temp(:,:,:,n,2),fe_tf_tot(n,5))
    enddo

    do n=ln1,ln2
       fe_tf_tot(n,3) = fe_tf_tot(n,1) + fe_tf_tot(n,2)
       fe_tf_tot(n,6) = fe_tf_tot(n,4) + fe_tf_tot(n,5)
    enddo

    deallocate(cfgamma1) 
    deallocate(cfgamma2) 
    deallocate(fe_tf_term) 
    deallocate(fe_tf_term_p) 
    deallocate(fe_tf_term_n) 
  
  End subroutine calc_transfer_fb
!!*****************************************************************************************!!!
Subroutine calc_spectral(g_ls2,fe_tot)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_ls2
    real, dimension(ln1:ln2,3,8),intent(inout) :: fe_tot
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1: n_fields) :: fields_ls2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2) :: f_ls2
    !temporary array to integrate over vp , mu and species
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2)::e_n_temp
    !Other factors used
    complex, dimension(:,:,:),pointer :: ptr2
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    integer ::n,lbg1,lbg2

    Allocate(v7d_spec(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2,ln1:ln2,1:14)) 
    Allocate(g_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(gc_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_v_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_z_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_kperp_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(coll_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(p_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(c_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma3(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
    
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1 
    lbg2 = lklmn0
    v7d_spec=cmplx(0.0,0.0)

    call calc_aux_fields(g_ls2, fields_ls2, f_ls2,.false.) !compute f and h 
    call get_cfgamma(f_ls2,fields_ls2,cfgamma1,cfgamma2)

    d_v_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then 
       if (hyp_on_h) then 
          call add_hypv_ak(h_,d_v_rhs,lbg1,lbg2)
       else
          call add_hypv_ak(f_,d_v_rhs,lbg1,lbg2)
       endif
    else
        call add_hypv_block(f_,d_v_rhs,lbg1,lbg2)
    endif

    d_z_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
       if (hyp_on_h) then 
          call add_hypz_ak(h_,d_z_rhs,lbg1,lbg2)
          if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,d_z_rhs,lbg1,lbg2)
       else  
          call add_hypz_ak(f_,d_z_rhs,lbg1,lbg2)
       endif  
    else
        call add_hypz_block(f_,d_z_rhs,lbg1,lbg2)
    endif

    d_kperp_rhs = cmplx(0.0,0.0)
    if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
       !/todo:check if hyp_on_h should affect perpendicular hyperdiffusions!
       call add_dfdxy(f_,d_kperp_rhs,lbg1,lbg2)
    endif

    !Collisions works with f
    coll_rhs = cmplx(0.0,0.0)                 
    if (collision_op.ne.'none') then
       call equ_collisions(f_,coll_rhs,.true.)
    endif

    !Gradient and curvature
    gc_rhs = cmplx(0.0,0.0)
    call add_dgdxy(g_1, gc_rhs, ptr1, pdg1di, pdg1dj,lbg1,lbg2)
    call add_dchidxy_orig(emfields,ptr4, gc_rhs,ptr2,ptr3,lbg1,lbg2,.true.)
    !call add_dchidxy_1(emfields,gc_rhs,ptr2,lbg1,lbg2,.true.)

    g_rhs  = cmplx(0.0,0.0)
    call add_dchidxy_orig(emfields,ptr4, g_rhs,ptr2,ptr3,lbg1,lbg2,.false.)
    
    c_rhs = cmplx(0.0,0.0)
    c_rhs = gc_rhs - g_rhs   

    !Parallel term
    p_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
       call equ_dzv(h_,p_rhs,lbg1,lbg2) !f_ is h
       if (hyp_on_h) then
          if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,p_rhs,lbg1,lbg2)
       endif
    else
        call equ_dfdzv(f_,p_rhs,lbg1,lbg2)
        call add_dchidz(emfields, ptr4, p_rhs, lbg1, lbg2)
    end if

    !TERM 1  FREE ENERGY      
    !W
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,1) = cfgamma1(:,:,:,:,:,n)*g_1(:,:,:,:,:,n)/2.
        call integral_vpw(v7d_spec(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,1))
    enddo                     
    
    
    if (n_spec.eq.1) then
    ! Electrostatic
        Do n=ln1,ln2
            e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
            v7d_spec(:,:,:,:,:,n,2) = cfgamma2(:,:,:,:,:,n)*g_1(:,:,:,:,:,n)/2.
            Call integral_vpw(v7d_spec(:,:,:,:,:,n,2),e_n_temp(:,:,:,n))
            if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
            call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,1))
        Enddo
    else
        ! Electrostatic with modification
        call get_electrostatic(fields_ls2,cfgamma3)
        Do n=ln1,ln2
            if(p_has_00_mode) cfgamma3(0,0,:,n)=cmplx(0.0,0.0)
            call sum_int_3d(cfgamma3(:,:,:,n),fe_tot(n,2,1))
        Enddo
    endif 
    
    
    ! FE
    do n =ln1,ln2
        fe_tot(n,3,1) = fe_tot(n,1,1) + fe_tot(n,2,1)
    enddo
 
    ! TERM 2  DRIVE
    ! W
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,3)= cfgamma1(:,:,:,:,:,n)*g_rhs(:,:,:,:,:,n)    
        call integral_vpw(v7d_spec(:,:,:,:,:,n,3),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,2))
    enddo    
    ! E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,4) = cfgamma2(:,:,:,:,:,n)*g_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,4),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,2))
    enddo       
    ! FE
    do n =ln1,ln2
        fe_tot(n,3,2) = fe_tot(n,1,2) + fe_tot(n,2,2)
    enddo
    
        
    !TERM 3 Hyp Z       
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,5) = cfgamma1(:,:,:,:,:,n)*d_z_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,5),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,3))
    enddo
    !E 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,6) = cfgamma2(:,:,:,:,:,n)*d_z_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,6),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,3))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,3) = fe_tot(n,1,3) + fe_tot(n,2,3)
    enddo 

    !TERM 4 Hyp V 
    !W
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,7) = cfgamma1(:,:,:,:,:,n)*d_v_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,7),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,4))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,8) = cfgamma2(:,:,:,:,:,n)*d_v_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,8),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,4))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,4) = fe_tot(n,1,4) + fe_tot(n,2,4)
    enddo    

    !TERM 5 hyp perp 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,9) = cfgamma1(:,:,:,:,:,n)*d_kperp_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,9),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,5))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,10) = cfgamma2(:,:,:,:,:,n)*d_kperp_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,10),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,5))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,5) = fe_tot(n,1,5) + fe_tot(n,2,5)
    enddo    
 
    !TERM 6 Collisions  
    if (collision_op.ne.'none') then
    ! W collisions
        do n =ln1,ln2
             e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
             v7d_spec(:,:,:,:,:,n,11) = cfgamma1(:,:,:,:,:,n)*coll_rhs(:,:,:,:,:,n)
             call integral_vpw(v7d_spec(:,:,:,:,:,n,11),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,6))
        enddo            
    ! E collisions
        do n =ln1,ln2
            e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
            v7d_spec(:,:,:,:,:,n,12) = cfgamma2(:,:,:,:,:,n)*coll_rhs(:,:,:,:,:,n)
            call integral_vpw(v7d_spec(:,:,:,:,:,n,12),e_n_temp(:,:,:,n))
            if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
            call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,6))
        enddo
       ! FE collisions
        do n =ln1,ln2
            fe_tot(n,3,6) = fe_tot(n,1,6) + fe_tot(n,2,6)
        enddo
     else
        do n =ln1,ln2
             fe_tot(n,1,6) = 0.0
             fe_tot(n,2,6) = 0.0
             fe_tot(n,3,6) = 0.0
        enddo
     endif

    !TERM 7 curvature 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,13) = cfgamma1(:,:,:,:,:,n)*c_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,13),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,7))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,13) = cfgamma2(:,:,:,:,:,n)*c_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,13),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,7))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,7) = fe_tot(n,1,7) + fe_tot(n,2,7)
    enddo    
 
    !TERM 8 Parallel 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,14) = cfgamma1(:,:,:,:,:,n)*p_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,14),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,8))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,14) = cfgamma2(:,:,:,:,:,n)*p_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,14),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,8))
    enddo    
  
    if (hyp_on_h) then
       do n =ln1,ln2
          fe_tot(n,1,8) = fe_tot(n,1,8) - fe_tot(n,1,4) - fe_tot(n,1,3)
          fe_tot(n,2,8) = fe_tot(n,2,8) - fe_tot(n,2,4) - fe_tot(n,2,3)
       enddo
    else
       do n =ln1,ln2
          if (.not.arakawa_zv) fe_tot(n,1,8) = fe_tot(n,1,8)-fe_tot(n,1,4)-fe_tot(n,1,3) !eliminated hyperdiffusions
          if (.not.arakawa_zv) fe_tot(n,2,8) = fe_tot(n,2,8)-fe_tot(n,2,4)-fe_tot(n,2,3) !eliminated hyperdiffusions
       enddo
    endif
    !FE
    do n =ln1,ln2
        fe_tot(n,3,8) = fe_tot(n,1,8) + fe_tot(n,2,8)
    enddo    
 
    deallocate(v7d_spec)
    deallocate(g_rhs)
    deallocate(gc_rhs)
    deallocate(p_rhs)
    deallocate(c_rhs)
    deallocate(d_v_rhs)
    deallocate(d_z_rhs)
    deallocate(d_kperp_rhs)
    deallocate(coll_rhs)
    deallocate(cfgamma1)
    deallocate(cfgamma2)
    deallocate(cfgamma3)
 
End Subroutine calc_spectral

Subroutine calc_spectral_mod(g_ls2,fe_tot)

    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2),intent(in) :: g_ls2
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2) :: cfgamma1m
    real, dimension(ln1:ln2,3,8),intent(inout) :: fe_tot
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, lm1:lm2,ln1:ln2,1:16) :: v6d_spec
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, lm1:lm2, ln1:ln2) :: cfgamma1_nvp
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, lm1:lm2, ln1:ln2) ::  g_1_nvp, g_rhs_nvp,&
                                                     &d_v_rhs_nvp, d_z_rhs_nvp, d_kperp_rhs_nvp, h_nvp
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, lm1:lm2, ln1:ln2) ::  coll_rhs_nvp, p_rhs_nvp, c_rhs_nvp
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1: n_fields) :: fields_ls2
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2) :: f_ls2
    !temporary array to integrate over vp , mu and species
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2)::e_n_temp
    !Other factors used
    complex, dimension(:,:,:),pointer :: ptr2
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    integer ::n,lbg1,lbg2

    Allocate(v7d_spec(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2,ln1:ln2,8)) 
    allocate(h_ls2(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)) 
    Allocate(g_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(gc_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_v_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_z_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(d_kperp_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(coll_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(p_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(c_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
    Allocate(cfgamma3(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
    
    ptr1 => null()
    ptr2 => null()
    ptr3 => null()
    ptr4 => null()
    lbg1 = 1 
    lbg2 = lklmn0
    
    v7d_spec=cmplx(0.0,0.0)
    v6d_spec=cmplx(0.0,0.0)

    call calc_aux_fields(g_ls2, fields_ls2, f_ls2,comp_h=.true.,h_out=h_ls2)    
    call get_cfgamma_mod(f_ls2,fields_ls2,cfgamma1m)
    call get_cfgamma(f_ls2,fields_ls2,cfgamma1,cfgamma2)
    
    d_v_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then 
       if (hyp_on_h) then 
          call add_hypv_ak(h_,d_v_rhs,lbg1,lbg2)
       else
          call add_hypv_ak(f_,d_v_rhs,lbg1,lbg2)
       endif
    else
        call add_hypv_block(f_,d_v_rhs,lbg1,lbg2)
    endif

    d_z_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
       if (hyp_on_h) then 
          call add_hypz_ak(h_,d_z_rhs,lbg1,lbg2)
       else  
          call add_hypz_ak(f_,d_z_rhs,lbg1,lbg2)
       endif  
    else
        call add_hypz_block(f_,d_z_rhs,lbg1,lbg2)
    endif

    d_kperp_rhs = cmplx(0.0,0.0)
    if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then
       !/todo:check if hyp_on_h should affect perpendicular hyperdiffusions!
       call add_dfdxy(f_,d_kperp_rhs,lbg1,lbg2)
    endif

    !Collisions works with f
    coll_rhs = cmplx(0.0,0.0)                 
    if (collision_op.ne.'none') then
       call equ_collisions(f_,coll_rhs,.true.)
    endif

    !Gradient and curvature
    gc_rhs = cmplx(0.0,0.0)
    call add_dgdxy(g_1, gc_rhs, ptr1, pdg1di, pdg1dj,lbg1,lbg2)
    call add_dchidxy_orig(emfields,ptr4, gc_rhs,ptr2,ptr3,lbg1,lbg2,.true.)

    g_rhs  = cmplx(0.0,0.0)
    call add_dchidxy_orig(emfields,ptr4, g_rhs,ptr2,ptr3,lbg1,lbg2,.false.)
    
    c_rhs = cmplx(0.0,0.0)
    c_rhs = gc_rhs - g_rhs   

    !Parallel term
    p_rhs = cmplx(0.0,0.0)
    if (arakawa_zv) then
       call equ_dzv(h_,p_rhs,lbg1,lbg2) !f_ is h
    else
        call equ_dfdzv(f_,p_rhs,lbg1,lbg2)
        call add_dchidz(emfields, ptr4, p_rhs, lbg1, lbg2)
    end if

    do n=ln1,ln2
        call integral_vp(cfgamma1m(:,:,:,:,:,n),cfgamma1_nvp(:,:,:,:,n))
        call integral_vp(g_ls2(:,:,:,:,:,n),g_1_nvp(:,:,:,:,n))
        call integral_vp(h_ls2(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,n),h_nvp(:,:,:,:,n))
        call integral_vp(g_rhs(:,:,:,:,:,n),g_rhs_nvp(:,:,:,:,n))
        call integral_vp(d_v_rhs(:,:,:,:,:,n),d_v_rhs_nvp(:,:,:,:,n))
        call integral_vp(d_z_rhs(:,:,:,:,:,n),d_z_rhs_nvp(:,:,:,:,n))
        call integral_vp(d_kperp_rhs(:,:,:,:,:,n),d_kperp_rhs_nvp(:,:,:,:,n))
        call integral_vp(coll_rhs(:,:,:,:,:,n),coll_rhs_nvp(:,:,:,:,n))
        call integral_vp(p_rhs(:,:,:,:,:,n),p_rhs_nvp(:,:,:,:,n))
        call integral_vp(c_rhs(:,:,:,:,:,n),c_rhs_nvp(:,:,:,:,n))
    enddo

    !TERM 1  FREE ENERGY      
    
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,1) = (conjg(g_1_nvp(:,:,:,:,n)))*g_1_nvp(:,:,:,:,n)/(2.*cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,1),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,1))
    enddo                     
    
    
    if (n_spec.eq.1) then
    ! Electrostatic
        Do n=ln1,ln2
            e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
            v7d_spec(:,:,:,:,:,n,1) = cfgamma2(:,:,:,:,:,n)*g_1(:,:,:,:,:,n)/2.
            Call integral_vpw(v7d_spec(:,:,:,:,:,n,1),e_n_temp(:,:,:,n))
            if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
            call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,1))
        Enddo
    else
        ! Electrostatic with modification
        call get_electrostatic(fields_ls2,cfgamma3)
        Do n=ln1,ln2
            if(p_has_00_mode) cfgamma3(0,0,:,n)=cmplx(0.0,0.0)
            call sum_int_3d(cfgamma3(:,:,:,n),fe_tot(n,2,1))
        Enddo
    endif 
    
   !FE
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,2) = (conjg(h_nvp(:,:,:,:,n)))*g_1_nvp(:,:,:,:,n)/(2.*(cfgamma1_nvp(:,:,:,:,n)))
        call integral_w(v6d_spec(:,:,:,:,n,2),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,3,1))
    enddo                     
    
    ! TERM 2  DRIVE
    ! W
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,3) = (conjg(g_1_nvp(:,:,:,:,n)))*g_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,3),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,2))
    enddo    
    ! E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,2) = cfgamma2(:,:,:,:,:,n)*g_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,2),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,2))
    enddo       
    ! FE
    do n =ln1,ln2
        fe_tot(n,3,2) = fe_tot(n,1,2) + fe_tot(n,2,2)
    enddo
    
        
    !TERM 3 Hyp Z       
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,5) = (conjg(g_1_nvp(:,:,:,:,n)))*d_z_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,5),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,3))
    enddo
    !E 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,3) = cfgamma2(:,:,:,:,:,n)*d_z_rhs(:,:,:,:,:,n)
        call integral_w(v7d_spec(:,:,:,:,:,n,3),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,3))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,3) = fe_tot(n,1,3) + fe_tot(n,2,3)
    enddo 

    !TERM 4 Hyp V 
    !W
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,7) = (conjg(g_1_nvp(:,:,:,:,n)))*d_v_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,7),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,4))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,4) = cfgamma2(:,:,:,:,:,n)*d_v_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,4),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,4))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,4) = fe_tot(n,1,4) + fe_tot(n,2,4)
    enddo    

    !TERM 5 hyp perp 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,9) = (conjg(g_1_nvp(:,:,:,:,n)))*d_kperp_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,9),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,5))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,5) = cfgamma2(:,:,:,:,:,n)*d_kperp_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,5),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,5))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,5) = fe_tot(n,1,5) + fe_tot(n,2,5)
    enddo    
 
    !TERM 6 Collisions  
    if (collision_op.ne.'none') then
    ! W collisions
        do n =ln1,ln2
             e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
             v6d_spec(:,:,:,:,n,11) = (conjg(g_1_nvp(:,:,:,:,n)))*coll_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
             call integral_w(v6d_spec(:,:,:,:,n,11),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,6))
        enddo            
    ! E collisions
        do n =ln1,ln2
            e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
            v7d_spec(:,:,:,:,:,n,6) = cfgamma2(:,:,:,:,:,n)*coll_rhs(:,:,:,:,:,n)
            call integral_vpw(v7d_spec(:,:,:,:,:,n,6),e_n_temp(:,:,:,n))
            if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
            call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,6))
        enddo
       ! FE collisions
        do n =ln1,ln2
            fe_tot(n,3,6) = fe_tot(n,1,6) + fe_tot(n,2,6)
        enddo
     else
        do n =ln1,ln2
             fe_tot(n,1,6) = 0.0
             fe_tot(n,2,6) = 0.0
             fe_tot(n,3,6) = 0.0
        enddo
     endif

    !TERM 7 curvature 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,13) = (conjg(g_1_nvp(:,:,:,:,n)))*c_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,13),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,7))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,7) = cfgamma2(:,:,:,:,:,n)*c_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,7),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,7))
    enddo    
    !FE
    do n =ln1,ln2
        fe_tot(n,3,7) = fe_tot(n,1,7) + fe_tot(n,2,7)
    enddo    
 
    !TERM 8 Parallel 
    !W 
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v6d_spec(:,:,:,:,n,15) = (conjg(g_1_nvp(:,:,:,:,n)))*p_rhs_nvp(:,:,:,:,n)/(cfgamma1_nvp(:,:,:,:,n))
        call integral_w(v6d_spec(:,:,:,:,n,15),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,1,8))
    enddo
    !E
    do n =ln1,ln2
        e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
        v7d_spec(:,:,:,:,:,n,8) = cfgamma2(:,:,:,:,:,n)*p_rhs(:,:,:,:,:,n)
        call integral_vpw(v7d_spec(:,:,:,:,:,n,8),e_n_temp(:,:,:,n))
        if(p_has_00_mode) e_n_temp(0,0,:,n) =cmplx(0.0,0.0)
        call sum_int_3d(e_n_temp(:,:,:,n),fe_tot(n,2,8))
    enddo    
  
    if (hyp_on_h) then
       do n =ln1,ln2
          fe_tot(n,1,8) = fe_tot(n,1,8) - fe_tot(n,1,4) - fe_tot(n,1,3)
          fe_tot(n,2,8) = fe_tot(n,2,8) - fe_tot(n,2,4) - fe_tot(n,2,3)
       enddo
    else
       do n =ln1,ln2
          if (.not.arakawa_zv) fe_tot(n,1,8) = fe_tot(n,1,8)-fe_tot(n,1,4)-fe_tot(n,1,3) !eliminated hyperdiffusions
          if (.not.arakawa_zv) fe_tot(n,2,8) = fe_tot(n,2,8)-fe_tot(n,2,4)-fe_tot(n,2,3) !eliminated hyperdiffusions
       enddo
    endif
    !FE
    do n =ln1,ln2
        fe_tot(n,3,8) = fe_tot(n,1,8) + fe_tot(n,2,8)
    enddo    
 
    deallocate(v7d_spec)
    deallocate(g_rhs)
    deallocate(h_ls2)
    deallocate(gc_rhs)
    deallocate(p_rhs)
    deallocate(c_rhs)
    deallocate(d_v_rhs)
    deallocate(d_z_rhs)
    deallocate(d_kperp_rhs)
    deallocate(coll_rhs)
    deallocate(cfgamma1)
    deallocate(cfgamma2)
    deallocate(cfgamma3)
 
End Subroutine calc_spectral_mod


!!!************************************************************************************!!!

  Subroutine calc_finalize_diag_transfer(which_model)

    Integer, intent(in) :: which_model
    Integer :: n,o,t

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       IF (.not.with_nvp) THEN
            Do o=1,3
                 Do n=ln1,ln2
                    close(FE_TF_FILE(n,o,which_model))
                Enddo
            Enddo
       ENDIF

       Do t=1,8
          Do o=1,3
             Do n=ln1,ln2
                close(FE_SP_FILE(n,o,t,which_model))
             Enddo
          Enddo
       enddo
    Endif

    select case(which_model)
    case(1)
       IF (.not.with_nvp) deallocate(tf_cumul_1)
       deallocate(spec_cumul_1)
    case(2)
       IF (.not.with_nvp) deallocate(tf_cumul_2)
       deallocate(spec_cumul_2)
    case(3)
       IF (.not.with_nvp) deallocate(tf_cumul_3)
       deallocate(spec_cumul_3)
    case(4)
       IF (.not.with_nvp) deallocate(tf_cumul_4)
       deallocate(spec_cumul_4)
    case(5)
       IF (.not.with_nvp) deallocate(tf_cumul_5)
       deallocate(spec_cumul_5)
    case(6)
       IF (.not.with_nvp) deallocate(tf_cumul_6)
       deallocate(spec_cumul_6)
    case(7)
       IF (.not.with_nvp) deallocate(tf_cumul_7)
       deallocate(spec_cumul_7)
    case default
    end select

  End subroutine calc_finalize_diag_transfer

  Subroutine calc_finalize_diag_transfer_h(which_model)

    Integer, intent(in) :: which_model
    Integer :: n,o

    If (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       IF (.not.with_nvp) THEN
          o=1
          Do n=ln1,ln2
             close(FE_TF_FILE(n,o,which_model))
          Enddo
       ENDIF
    Endif

    select case(which_model)
    case(1)
       IF (.not.with_nvp) deallocate(tf_cumul_1)
    case(2)
       IF (.not.with_nvp) deallocate(tf_cumul_2)
    case(3)
       IF (.not.with_nvp) deallocate(tf_cumul_3)
    case(4)
       IF (.not.with_nvp) deallocate(tf_cumul_4)
    case(5)
       IF (.not.with_nvp) deallocate(tf_cumul_5)
    case(6)
       IF (.not.with_nvp) deallocate(tf_cumul_6)
    case(7)
       IF (.not.with_nvp) deallocate(tf_cumul_7)
    case default
    end select

  End subroutine calc_finalize_diag_transfer_h

!!!*****************************************************************************************!!!

  Subroutine write_transfer(n_shells,which_model)

    integer, intent(in) :: which_model, n_shells
    integer :: ls1,ls2,n,o
  

    select case(which_model)
    case(1,2,3,4,5,6) 
    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       Do o=1,3
          Do n=ln1,ln2
             Write(FE_TF_FILE(n,o,which_model),*) '#', time
             Do ls2=1,n_shells
                Do ls1=1,n_shells
                   Write(FE_TF_FILE(n,o,which_model),"(2I4,2ES12.4,ES12.4)") ls2,ls1,k_value(ls2),&
                        k_value(ls1),fe_tf_tot(ls2,ls1,n,o)
                Enddo
                Write(FE_TF_FILE(n,o,which_model),*)
             Enddo
             Write(FE_TF_FILE(n,o,which_model),*)
             call flush(FE_TF_FILE(n,o,which_model))
          enddo
       Enddo
    Endif 

   case(7) 
    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        do o=1,3    
            Do n=ln1,ln2
                Write(FE_TF_FILE(n,o,which_model),*) '#', time
                Do ls2=1,n_shells
                    Do ls1=1,n_shells
                        Write(FE_TF_FILE(n,o,which_model),"(2I4,2ES12.4,3ES12.4)") ls2,ls1,k_value(ls2),&
                            k_value(ls1),fe_tf_tot(ls2,ls1,n,o),fe_tf_tot(ls2,ls1,n,o+3),&
                                     fe_tf_tot(ls2,ls1,n,o) + fe_tf_tot(ls2,ls1,n,o+3)
                    Enddo
                    Write(FE_TF_FILE(n,o,which_model),*)
                Enddo
                Write(FE_TF_FILE(n,o,which_model),*)
                call flush(FE_TF_FILE(n,o,which_model))
            Enddo
        Enddo
   Endif
   case default
   End Select

  End subroutine write_transfer


  Subroutine write_transfer_h(n_shells,which_model)

    integer, intent(in) :: which_model, n_shells
    integer :: ls1,ls2,n,o
  

    select case(which_model)
    case(1,2,3,4,5,6) 
    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       o=1
          Do n=ln1,ln2
             Write(FE_TF_FILE(n,o,which_model),*) '#', time
             Do ls2=1,n_shells
                Do ls1=1,n_shells
                   Write(FE_TF_FILE(n,o,which_model),"(2I4,2ES12.4,ES12.4)") ls2,ls1,k_value(ls2),&
                        k_value(ls1),fe_tf_tot(ls2,ls1,n,o)
                Enddo
                Write(FE_TF_FILE(n,o,which_model),*)
             Enddo
             Write(FE_TF_FILE(n,o,which_model),*)
             call flush(FE_TF_FILE(n,o,which_model))
          enddo
    Endif 

   case(7) 
    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
        o=1    
            Do n=ln1,ln2
                Write(FE_TF_FILE(n,o,which_model),*) '#', time
                Do ls2=1,n_shells
                    Do ls1=1,n_shells
                        Write(FE_TF_FILE(n,o,which_model),"(2I4,2ES12.4,3ES12.4)") ls2,ls1,k_value(ls2),&
                            k_value(ls1),fe_tf_tot(ls2,ls1,n,o),fe_tf_tot(ls2,ls1,n,o+3),&
                                     fe_tf_tot(ls2,ls1,n,o) + fe_tf_tot(ls2,ls1,n,o+3)
                    Enddo
                    Write(FE_TF_FILE(n,o,which_model),*)
                Enddo
                Write(FE_TF_FILE(n,o,which_model),*)
                call flush(FE_TF_FILE(n,o,which_model))
            Enddo
   Endif
   case default
   End Select

  End subroutine write_transfer_h

!!!*****************************************************************************************!!!

  Subroutine write_spectral(n_shells,which_model)

    integer, intent(in) :: n_shells, which_model
    integer :: ls1,n,o,t
    real, dimension(n_shells) :: w_k_value,delta_k_value
   
    select case(which_model)
    
    case(1,2)
        w_k_value(:) = k_value(:)
        delta_k_value(:) = k_value(2)
    
    case(3,4,5,6)
        do ls1=1,n_shells
            if (ls1.eq.1) then
                w_k_value(ls1) = k_value(ls1)/2.
                delta_k_value = k_value(ls1)
            else 
                w_k_value(ls1) = (k_value(ls1) + k_value(ls1-1))/2.
                delta_k_value(ls1) = (k_value(ls1) - k_value(ls1-1))
            endif
        enddo
    end select

    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then
       do t = 1,8
          Do o=1,3
             Do n=ln1,ln2
                Do ls1=1,n_shells
              !  Write(FE_SP_FILE(n,o,t,which_model),"(I4,4ES12.4)") ls1, k_value(ls1), w_k_value(ls1), delta_k_value(ls1), fe_tot(ls1,n,o,t)
                    Write(FE_SP_FILE(n,o,t,which_model),"(I4,2ES12.4)") ls1, k_value(ls1), fe_tot(ls1,n,o,t)
                Enddo
                Write(FE_SP_FILE(n,o,t,which_model),*)
                call flush(FE_SP_FILE(n,o,t,which_model))
             enddo
          enddo
!         write(*,*) "curvature and parallel w", sum(fe_tot(:,n,1,7)),sum(fe_tot(:,n,1,8))
!         write(*,*) "curvature and parallel e", sum(fe_tot(:,n,2,7)),sum(fe_tot(:,n,2,8))
!         write(*,*) "hyp_z and hyp_v", sum(fe_tot(:,n,3,3)),sum(fe_tot(:,n,3,4))
!         write(*,*) "collisoins", sum(fe_tot(:,n,3,6))
!         write(*,*) "gradient", sum(fe_tot(:,n,3,2))
!         write(*,*) "free energy", sum(fe_tot(:,n,3,1))
       enddo
    Endif

  End subroutine write_spectral

#endif

end Module diag_Gyro_LES_transfer
