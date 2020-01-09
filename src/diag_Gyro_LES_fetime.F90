#include "redef.h"
#include "intrinsic_sizes.h"
#include "switches.h"
! Fre energy versus time
Module diag_Gyro_LES_fetime
  
  Use par_mod 
  Use file_io, only: get_unit_nr
  use communications
  use antenna, only: antenna_type
  use collisions, only: equ_collisions
  use diagnostics_energy, only: get_cfgamma_f,get_cfgamma_h,get_cfgamma_ant,&
       &energy_integral
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

  PUBLIC :: initialize_diag_fe_time, diagnostic_fe_time, finalize_diag_fe_time

  !quantities for free energy conservation vs time
  integer, public :: istep_fe_time = 0

  PRIVATE

  !***** Parameters for free energy conservation vs time *****
  Character(Len=8) :: filestat='replace', filepos='rewind'
  Integer :: FE_TIME_FILE
  Integer, dimension(:,:),allocatable :: FE_TIME_FILE_SPEC
  Complex, dimension(:,:,:,:,:,:),allocatable  :: temp_rhs ,diss_rhs, temp_fe6d
  Complex, dimension(:,:,:,:,:,:),allocatable  :: cfgamma_last
  ! Conditional terms if filter defined
  real, dimension(:), allocatable :: time_storage
  Complex, Dimension(:,:,:,:,:,:), Allocatable  :: g_last
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: f_last
  !Complex, Dimension(:,:,:,:,:,:), Allocatable :: h_last
  Complex, Dimension(:,:,:,:), Allocatable::  emfields_last
  Complex, Dimension(:,:,:,:,:,:), Allocatable  :: filter_g_1
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: filter_f_ 
  Complex, Dimension(:,:,:,:), Allocatable::  filter_emfields
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: filter_nlt    
  Complex, Dimension(:,:,:,:,:,:), Allocatable :: nlt_filter_g_1 
  Complex, Dimension(:,:,:,:,:,:), Allocatable:: filter_nlt_filter_g_1  
!  Complex, Dimension(:,:,:,:,:,:), Allocatable  :: t_fgs

contains

  !!!**************************************************************************!!!
!!!******************** free energy balance (vs time)************************!!! 
!!!**************************************************************************!!!

  Subroutine initialize_diag_fe_time

    Implicit None

    Integer ::  n, o
    Character(len=20),dimension(3) :: time_file_spec_label

#ifdef with_extended_diags

    time_file_spec_label = (/'/w_time_ ','/e_time_ ','/fe_time_'/)

    allocate(FE_TIME_FILE_SPEC(ln1:ln2,3))

    IF (mype.eq.0) then    
       call get_unit_nr(FE_TIME_FILE)
       OPEN(FE_TIME_FILE, file=trim(diagdir)//&
            &'/fe_time'//trim(file_extension),form = 'FORMATTED', &
            status=filestat, position=filepos)

       WRITE(FE_TIME_FILE,*) "#  1. time"
       WRITE(FE_TIME_FILE,*) "#  2. Free energy"
       WRITE(FE_TIME_FILE,*) "#  3. Drive"
       WRITE(FE_TIME_FILE,*) "#  4. Dissipation v"
       WRITE(FE_TIME_FILE,*) "#  5. Dissipation z"
       WRITE(FE_TIME_FILE,*) "#  6. Dissipation kperp"
       WRITE(FE_TIME_FILE,*) "#  7. Collision"
       WRITE(FE_TIME_FILE,*) "#  8. Parallel "        
       WRITE(FE_TIME_FILE,*) "#  9. Curvature "        
       WRITE(FE_TIME_FILE,*) "#  10. D FE/ dt (RHS)" 
       WRITE(FE_TIME_FILE,*) "#  11. D Fe /dt (LHS) "        
       WRITE(FE_TIME_FILE,*) "#  12. LHS - RHS"
       if (fracy.gt.0.and.nonlinear) then
100       format ('#',"        1",15I12)
          WRITE(FE_TIME_FILE,*) "#  13. subgrid"
          WRITE(FE_TIME_FILE,*) "#  14. filter fe drive"
          WRITE(FE_TIME_FILE,*) "#  15. filter fe dissipation"
          write(FE_TIME_FILE,100) 2,3,4,5,6,7,8,9,10,11,12,13,14,15
       else
110       format ('#',"        1",12I12)
          write(FE_TIME_FILE,110) 2,3,4,5,6,7,8,9,10,11,12
       endif
    Endif


    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then    
       Do o=1,3
          DO n=ln1,ln2                
             call get_unit_nr(FE_TIME_FILE_SPEC(n,o))
             OPEN(FE_TIME_FILE_SPEC(n,o), file=trim(diagdir)//&
                  &trim(time_file_spec_label(o))//trim(spec(n)%name)//&
                  &trim(file_extension),form = 'FORMATTED', &
                  status=filestat, position=filepos)

             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  1. time"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  2. Free energy"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  3. Drive"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  4. Dissipation v"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  5. Dissipation z"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  6. Dissipation kperp"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  7. Collision"
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  8. Parallel "        
             WRITE(FE_TIME_FILE_SPEC(n,o),*) "#  9. Curvature "        
190          format ('#',"        1",9I12)
             write(FE_TIME_FILE_SPEC(n,o),190) 2,3,4,5,6,7,8,9
          enddo
       enddo
    Endif
#endif

  End Subroutine initialize_diag_fe_time

!!!*****************************************************************************************!!!
  Subroutine diagnostic_fe_time(call_number)

    integer, intent(in) :: call_number
#ifdef with_extended_diags
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2):: e_temp
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2):: e_n_temp
    complex, dimension(:,:,:),pointer :: ptr2
    complex, dimension(:,:,:,:),pointer :: ptr1, ptr3
    complex, dimension(:,:,:,:,:,:), pointer :: ptr4
    real :: rhs, lhs, fe_1, fe_2
    real,dimension(12) :: e_tot
    real,dimension(ln1:ln2,9,3) :: e_tot_n
    real,dimension(3) :: e_opt
    integer ::n,lbg1,lbg2,stage

    if(.not.(allocated(time_storage))) allocate(time_storage(2)) 
    if (.not.allocated(g_last))  allocate(g_last(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2))
    if (.not.allocated(f_last))  allocate(f_last(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    !if (.not.allocated(h_last).and.arakawa_zv)  allocate(h_last(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2))
    if (.not.allocated(emfields_last))  allocate(emfields_last(li1:li2,lj1:lj2,lbz:ubz,1:n_fields))

    if(call_number.eq.1) then
       g_last = g_1
       f_last = f_
       !if (arakawa_zv) h_last = h_
       emfields_last = emfields
       time_storage(1) = time
    endif

    if (call_number.eq.2) then

       Allocate(g_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(diss_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(temp_rhs(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(temp_fe6d(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(cfgamma(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(cfgamma_last(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(cfgamma1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(cfgamma2(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
       Allocate(cfgamma3(li1:li2, lj1:lj2, lk1:lk2, ln1:ln2))
          
       if (fracy.gt.0.and.nonlinear) then
          Allocate(filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
          Allocate(filter_f_(lbi:ubi, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2))   
          Allocate(filter_emfields(li1:li2,lj1:lj2,lbz:ubz,1:n_fields))
          Allocate(nlt(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2)) 
          Allocate(filter_nlt (li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))     
          Allocate(nlt_filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))  
          Allocate(filter_nlt_filter_g_1(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))   
       endif

       time_storage(2) = time
       ptr1 => null()
       ptr2 => null()
       ptr3 => null()
       ptr4 => null()
       lbg1 = 1 
       lbg2 = lklmn0  
       stage = 0

       e_tot = 0.0
       e_tot_n = 0.0

       call get_cfgamma_f(f_,emfields,cfgamma)  !compute f
       call get_cfgamma(f_,emfields,cfgamma1,cfgamma2)
       !with antenna, the free energy as well as the free energy per species will contain the 
       !antenna contribution -- however, with electromagnetics the entropy and electrostatic parts
       !are not fully separated -- the entropy contains picks up some Bperp terms, the electrostatic part 
       !picks up some Bpar terms!
       if (any(antenna_type.eq.(/2,3/))) then
          allocate(cfgamma_ant(li1:li2, lj1:lj2, lk1:lk2, ll1:ll2, lm1:lm2, ln1:ln2))
          call get_cfgamma_ant(cfgamma_ant)
          !add antenna term to cfgamma for calculation of total energy
          cfgamma=cfgamma+cfgamma_ant
          !add antenna term to cfgamma1 (as it already contains some electromagnetic contributions)
          !to have it also in the species-dependent energy
          cfgamma1=cfgamma1+cfgamma_ant
          deallocate(cfgamma_ant)
       endif
          
       !--------------------------------------------------------------
       !Term 4 Dissipation  v  
       temp_rhs = cmplx(0.0,0.0)
       diss_rhs = cmplx(0.0,0.0)
       temp_fe6d = cmplx(0.0,0.0)
          
       if (arakawa_zv) then 
          if (hyp_on_h) then 
             call add_hypv_ak(h_,temp_rhs,lbg1,lbg2)
          else
             call add_hypv_ak(f_,temp_rhs,lbg1,lbg2)
          endif
       else
           call add_hypv_block(f_,temp_rhs,lbg1,lbg2)
       endif

       diss_rhs=diss_rhs+temp_rhs
       e_temp=cmplx(0.0,0.0)
       temp_fe6d  = cfgamma*temp_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(4))


       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,4,1))
       Enddo

       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,4,2))
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,4,3) = e_tot_n(n,4,1) + e_tot_n(n,4,2) 
       Enddo

       !--------------------------------------------------------------------------
       !Term 5 Dissipation  z  
       temp_rhs = cmplx(0.0,0.0)
       temp_fe6d = cmplx(0.0,0.0)
       
       if (arakawa_zv) then
          if (hyp_on_h) then 
             call add_hypz_ak(h_,temp_rhs,lbg1,lbg2)
             if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,temp_rhs,lbg1,lbg2)
          else  
             call add_hypz_ak(f_,temp_rhs,lbg1,lbg2)
          endif  
       else
           call add_hypz_block(f_,temp_rhs,lbg1,lbg2)
       endif


       diss_rhs=diss_rhs+temp_rhs
       e_temp=cmplx(0.0,0.0)
       temp_fe6d  = cfgamma*temp_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(5))

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,5,1))
       Enddo

       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,5,2))
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,5,3) = e_tot_n(n,5,1) + e_tot_n(n,5,2) 
       Enddo

       !------------------------------------------------------------------------------------
       !Term 6 Dissipation  kperp  
       temp_rhs = cmplx(0.0,0.0)
       temp_fe6d = cmplx(0.0,0.0)
       if ((hyp_x.gt.0.0).or.(hyp_y.gt.0.0).or.(hyp_perp.gt.0.0).or.(GyroLES)) then !same here
          !/todo:check if hyp_on_h should affect perpendicular hyperdiffusions!
          !if (hyp_on_h) then
          !   call add_dfdxy(h_,temp_rhs,lbg1,lbg2) 
          !else
             call add_dfdxy(f_,temp_rhs,lbg1,lbg2) 
          !endif
       endif

       diss_rhs=diss_rhs+temp_rhs
       e_temp = cmplx(0.0,0.0)
       temp_fe6d  = cfgamma*temp_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(6))

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,6,1))
       Enddo

       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,6,2))
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,6,3) = e_tot_n(n,6,1) + e_tot_n(n,6,2) 
       Enddo

       !------------------------------------------------------------------------------
       !Curvature and drive term
       temp_rhs = cmplx(0.0,0.0)
       call add_dgdxy(g_1, temp_rhs, ptr1, pdg1di, pdg1dj,lbg1,lbg2)
       call add_dchidxy_orig(emfields,ptr4, temp_rhs,ptr2,ptr3,lbg1,lbg2,.true.)

       g_rhs  = cmplx(0.0,0.0)
       temp_fe6d  = cmplx(0.0,0.0)
       call add_dchidxy_orig(emfields,ptr4, g_rhs,ptr2,ptr3,lbg1,lbg2,.false.)
       temp_rhs = temp_rhs - g_rhs   


       !Term 9 curvature  
       e_temp =cmplx(0.0,0.0)
       temp_fe6d  = cfgamma*temp_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(9))
       if (INDEX(magn_geometry,'slab').eq.1) e_tot(9) = 0.0

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,9,1))
       Enddo

       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,9,2))
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,9,3) = e_tot_n(n,9,1) + e_tot_n(n,9,2) 
       Enddo

       !Term 3 Drive
       e_temp=cmplx(0.0,0.0)
       temp_fe6d = cfgamma*g_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(3))

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*g_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,3,1))
       Enddo

       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*g_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,3,2))
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,3,3) = e_tot_n(n,3,1) + e_tot_n(n,3,2) 
       Enddo
        
       !----------------------------------------------------------------------
       !Term 8 parallel  
       temp_rhs = cmplx(0.0,0.0)
       temp_fe6d = cmplx(0.0,0.0)

       !Parallel term goes always with h
       if (arakawa_zv) then
           call equ_dzv(h_,temp_rhs,lbg1,lbg2)
           if (hyp_on_h) then
              if (hypz_compensation) call equ_comp_hypz(emfields,ptr4,temp_rhs,lbg1,lbg2)
           endif
       else
           call equ_dfdzv(f_,temp_rhs,lbg1,lbg2)
           call add_dchidz(emfields, ptr4, temp_rhs, lbg1, lbg2)
       end if

       e_temp=cmplx(0.0,0.0)
       temp_fe6d = cfgamma*temp_rhs
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(8))
       !Substract hyper diffusions
       if (hyp_on_h) then
          e_tot(8) = e_tot(8) - e_tot(4) - e_tot(5)
       else
          if (.not.arakawa_zv) e_tot(8)= e_tot(8)-e_tot(4)-e_tot(5) !eliminated hyperdiffusions
       endif

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,8,1))
          !Substract hyper diffusions
          if (hyp_on_h) then
             e_tot_n(n,8,1) = e_tot_n(n,8,1) - e_tot_n(n,4,1) - e_tot_n(n,5,1)
          else
             if (.not.arakawa_zv) e_tot_n(n,8,1)= e_tot_n(n,8,1)-e_tot_n(n,4,1)-e_tot_n(n,5,1) !eliminated hyperdiffusions
          endif
       Enddo
       ! Electrostatic
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,8,2))
          !Substract hyper diffusions
          if (hyp_on_h) then
             e_tot_n(n,8,2) = e_tot_n(n,8,2) - e_tot_n(n,4,2) - e_tot_n(n,5,2)
          else
             if (.not.arakawa_zv) e_tot_n(n,8,2)= e_tot_n(n,8,2)-e_tot_n(n,4,2)-e_tot_n(n,5,2) !eliminated hyperdiffusions
          endif
       Enddo

       ! Total
       Do n =ln1,ln2
          e_tot_n(n,8,3) = e_tot_n(n,8,1) + e_tot_n(n,8,2) 
       Enddo

       !----------------------------------------------------------------------------------
       !Term 7  collision     
       temp_fe6d = cmplx(0.0,0.0)                 
       if (collision_op.ne.'none') then
          !collisions on f_
          call equ_collisions(f_,temp_rhs,replace_rhs=.true.)
       endif

       diss_rhs=diss_rhs+temp_rhs

       if (collision_op.ne.'none') then 
          e_temp=cmplx(0.0,0.0)
          temp_fe6d = cfgamma*temp_rhs
          Call energy_integral(temp_fe6d,e_temp)
          if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
          call sum_int_3d(e_temp,e_tot(7))
       else
          e_tot(7) = 0.0
       endif

       if (collision_op.ne.'none') then 
          ! Entropy
          Do n=ln1,ln2
             e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
             temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
             Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,7,1))
          Enddo
          !    
          ! Electrostatic
          Do n=ln1,ln2
             e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
             temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*temp_rhs(:,:,:,:,:,n)
             Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,7,2))
          Enddo

          ! Total
          Do n =ln1,ln2
             e_tot_n(n,7,3) = e_tot_n(n,7,1) + e_tot_n(n,7,2) 
          Enddo

       else
          Do n =ln1,ln2
             e_tot_n(n,7,1) = 0  
             e_tot_n(n,7,2) = 0  
             e_tot_n(n,7,3) = 0 
          enddo
       endif

       !--------------------------------------------------------------------------
       !Term 1 Time
       e_tot(1) = time
       Do n=ln1,ln2
          e_tot_n(n,1,1) = time
          e_tot_n(n,1,2) = time
          e_tot_n(n,1,3) = time
       Enddo

       !--------------------------------------------------------------------------
       !Term 2 Free energy  
       temp_fe6d = cmplx(0.0,0.0)    
       e_temp=cmplx(0.0,0.0)
       temp_fe6d = cfgamma*g_1/2.
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,e_tot(2))

       ! Entropy
       Do n=ln1,ln2
          e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
          temp_fe6d(:,:,:,:,:,n) = cfgamma1(:,:,:,:,:,n)*g_1(:,:,:,:,:,n)/2.
          Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
          if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
          call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,2,1))
       Enddo


       if (n_spec.eq.1) then
          ! Electrostatic
          Do n=ln1,ln2
             e_n_temp(:,:,:,n)=cmplx(0.0,0.0)
             temp_fe6d(:,:,:,:,:,n) = cfgamma2(:,:,:,:,:,n)*g_1(:,:,:,:,:,n)/2.
             Call integral_vpw(temp_fe6d(:,:,:,:,:,n),e_n_temp(:,:,:,n))
             if(p_has_00_mode) e_n_temp(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(e_n_temp(:,:,:,n),e_tot_n(n,2,2))
          Enddo

       else
          ! Electrostatic with modification
          call get_electrostatic(emfields,cfgamma3)
          Do n=ln1,ln2
             if(p_has_00_mode) cfgamma3(0,0,:,n)=cmplx(0.0,0.0)
             call sum_int_3d(cfgamma3(:,:,:,n),e_tot_n(n,2,2))
          Enddo
       endif 
          
       ! Total
       Do n =ln1,ln2
          e_tot_n(n,2,3) = e_tot_n(n,2,1) + e_tot_n(n,2,2) 
       Enddo



       !-------------------------------------------------------------------------
       ! Term 10: RHS, 11: LHS and 12: LHS - RHS 
       !if (arakawa_zv) then
       !   call get_cfgamma_h(h_last,cfgamma_last)
       !else
          call get_cfgamma_f(f_last,emfields_last,cfgamma_last)
       !endif


       temp_fe6d = cmplx (0.0,0.0)
       e_temp  = cmplx(0.0)
       fe_1 = 0.0
       temp_fe6d =cfgamma_last*g_last/2.
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,fe_1)

       temp_fe6d = cmplx (0.0,0.0)
       e_temp = cmplx(0.0)
       fe_2 = 0.0
       temp_fe6d = cfgamma*g_1/2.
       Call energy_integral(temp_fe6d,e_temp)
       if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
       call sum_int_3d(e_temp,fe_2)

       lhs = 0.0
       rhs = 0.0
       lhs  = (fe_2 - fe_1)/(time_storage(2) - time_storage(1))
       rhs  = e_tot(3) + e_tot(4) + e_tot(5) + e_tot(6) + e_tot(7) + e_tot(8) + e_tot(9)
       e_tot(10) = rhs
       e_tot(11) = lhs
       e_tot(12) = lhs - rhs

       !--------------------------------------------------------------------------------
       !Terms 13,14,15 sub-grids contribution (optional)
       if  (fracy.gt.0.and.nonlinear)  then
          e_opt = 0.0  
          e_temp = cmplx(0.0,0.0)
          Call filter(g_1, filter_g_1)
          call calc_aux_fields(filter_g_1,filter_emfields,filter_f_,comp_h=.false.)

          nlt = cmplx(0.0,0.0)
          CALL this_nonlinear_term%add(g_1, ptr1, emfields,&
               &ptr2, ptr3, nlt, lbg1, lbg2, stage)

          Call filter(nlt, filter_nlt)

          nlt_filter_g_1 = cmplx(0.0,0.0)
          CALL this_nonlinear_term%add(filter_g_1, ptr1, filter_emfields,&
               &ptr2, ptr3, nlt_filter_g_1, lbg1, lbg2, stage)
          Call filter(nlt_filter_g_1, filter_nlt_filter_g_1)     

          temp_rhs = filter_nlt - filter_nlt_filter_g_1
          call filter(temp_rhs,temp_rhs)
          call filter(cfgamma,cfgamma)

          temp_fe6d = cmplx(0.0,0.0)
          temp_fe6d  = cfgamma*temp_rhs
          Call energy_integral(temp_fe6d,e_temp)
          call sum_int_3d(e_temp,e_opt(1))

          temp_fe6d  = cmplx(0.0,0.0)
          temp_fe6d  = cfgamma*g_rhs
          Call energy_integral(temp_fe6d,e_temp)
          if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
          call sum_int_3d(e_temp,e_opt(2))

          temp_fe6d = cmplx(0.0,0.0)
          temp_fe6d   = cfgamma*diss_rhs
          Call energy_integral(temp_fe6d,e_temp)
          if(p_has_00_mode) e_temp(0,0,:)=cmplx(0.0,0.0)
          call sum_int_3d(e_temp,e_opt(3))

       endif

       IF (mype.eq.0) then    
          if (fracy.gt.0.and.nonlinear) then
             write(FE_TIME_FILE,"(15ES12.4)") e_tot, e_opt
          else
             WRITE(FE_TIME_FILE,"(12ES12.4)") e_tot
          endif
         call flush(FE_TIME_FILE)
       Endif


       IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then    
          DO n=ln1,ln2
             write(FE_TIME_FILE_SPEC(n,1),"(9ES12.4)") e_tot_n(n,:,1)
             write(FE_TIME_FILE_SPEC(n,2),"(9ES12.4)") e_tot_n(n,:,2) 
             write(FE_TIME_FILE_SPEC(n,3),"(9ES12.4)") e_tot_n(n,:,3) 
            call flush(FE_TIME_FILE_SPEC(n,1))
            call flush(FE_TIME_FILE_SPEC(n,2))
            call flush(FE_TIME_FILE_SPEC(n,3))
          enddo
       Endif


       Deallocate(time_storage)    
       Deallocate(g_last)
       Deallocate(f_last)
       !Deallocate(h_last)
       Deallocate(emfields_last)
       Deallocate(g_rhs,diss_rhs,temp_rhs,temp_fe6d)
       Deallocate(cfgamma,cfgamma_last, cfgamma1, cfgamma2, cfgamma3)
       
       if (fracy.gt.0.and.nonlinear) then
          Deallocate(filter_g_1)
          Deallocate(filter_f_)
          Deallocate(filter_emfields)
          Deallocate(nlt)
          Deallocate(filter_nlt)
          Deallocate(nlt_filter_g_1)
          Deallocate(filter_nlt_filter_g_1)
       endif

    endif
#endif
  End Subroutine diagnostic_fe_time

!!!************************************************************************************!!!

  Subroutine finalize_diag_fe_time  
#ifdef with_extended_diags

    Integer :: n

    IF (mype.eq.0) then    
       Close (FE_TIME_FILE)
    endif

    IF (mype.eq.pexyzvwspec(0,0,0,0,0,my_pespec)) then    
       do n=ln1,ln2 
          Close (FE_TIME_FILE_SPEC(n,1))
          Close (FE_TIME_FILE_SPEC(n,2))
          Close (FE_TIME_FILE_SPEC(n,3))
       enddo
    endif

    Deallocate(FE_TIME_FILE_SPEC)
#endif
  End Subroutine finalize_diag_fe_time

end Module diag_Gyro_LES_fetime
