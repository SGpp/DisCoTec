#include "switches.h"
#include "redef.h"
#include "intrinsic_sizes.h"

module field_solver

  use discretization
  use par_in
  use par_other, only: print_ini_msg, n_fields
  use charge_curr_dens
  use adiabatic_response  
  use gyro_average_df_mod
  use gyro_average_ff_mod
  use fieldsolver_ff 
  use fieldsolver_df
  use antenna, only: check_antenna,initialize_antenna,antenna_type

  implicit none

  public::  mem_est_field_solver, check_field_solver, initialize_field_solver, compute_emfields, finalize_field_solver
  ! routines for adaptive grids
  public:: compute_emfields_adptv
  private

  integer :: init_status = 0

  complex, dimension(:,:,:,:), allocatable :: moments   

contains

  Real Function mem_est_field_solver(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

   !moments
    mem_loc=SIZE_OF_COMPLEX_MB*lijk0*n_fields

    mem_loc=mem_est_adiabatic_response(mem_loc)
    mem_loc=mem_est_charge_curr_dens(mem_loc)
    if (xy_local) then
       mem_loc=mem_est_gyro_average_ff(mem_loc)
       mem_loc=mem_est_field_solve_ff(mem_loc)
    else
       mem_loc=mem_est_gyro_average_df(mem_loc)
       mem_loc=mem_loc+estimate_memory_fieldsolver_df(nx0)
    endif

    mem_est_field_solver=mem_req_in+mem_loc
  end Function mem_est_field_solver

  subroutine check_field_solver
    logical :: write_pe
    integer:: n

    write_pe = ((mype.le.0).AND.(print_ini_msg))

    !check for adiabatic electrons
    do n=0,n_spec-1
       if (.not.spec(n)%passive.and.spec(n)%charge.eq.-1) adiabatic_electrons=.false.
    enddo

    !set n_fields  (default: 1)
    if (beta.ne.0.) then
       if (adiabatic_electrons.or.trap_pass) then
          stop "electromagnetic simulation with adiabatic electrons does not make sense"
       else
          n_fields = 2
          if (bpar) n_fields = n_fields + 1
       end if
    end if

    !field modifications
    if (write_pe) then 
       if (delzonal) write(*,"(A)") "The zonal component of phi will be deleted."
       if (delzonal_fields) then
          if (n_fields.ge.2) then
             write(*,"(A)") "The zonal component of A_par will be deleted."
          else
             delzonal_fields = .false.
             write(*,"(A)") "n_fields < 2: set parameter delzonal_fields=.false."
          endif
       endif

       if (only_Er) write(*,"(A)") "With purely radial electric field (phi=flux_surface_average(phi))"
       if (del_phi) write(*,"(A)") "Without electric field: phi=0"
       if (del_fields) write(*,"(A)") "All electromagnetic fields set to zero."
       if (abs(add_zonal_phi).gt.epsilon(0.0)) write(*,"(A,F12.6)") &
            & "Adding an zonal component to phi with amp. ", add_zonal_phi
    endif
    
    IF(trap_pass.and.n_fields.gt.1) THEN
       STOP "Hybrid model only available for electrostatic cases" 
    END IF
    if (xy_local) call check_antenna

    if (adiabatic_electrons.and.(debye2.gt.0)) stop "Debye2 does not make sense with adiabatic electrons"
    if (no_electron_response) then 
       debye2=0.
       adiabatic_electrons=.false.
       if (.not.xy_local) stop 'No-response electrons are currently only implemented for the local code.'
    endif
       

  end subroutine check_field_solver

  subroutine initialize_field_solver


    !no alternative routines for the following modules
    if(init_status.eq.0) then
       if (xy_local) then
          call initialize_gyro_average_ff
       else
          call initialize_gyro_average_df
       end if
       call check_for_adiabatic_species
       if (adiabatic_electrons.or.trap_pass) call initialize_adiabatic_response

       if(xy_local) then
          call initialize_field_solve_ff
       else
          call initialize_field_solve_df
       end if
       allocate(moments(li1:li2,lj1:lj2,lk1:lk2,1:n_fields))
    end if

    !charge_curr_dens has alternative routines and therefore has to be re-initialized
    call initialize_charge_curr_dens
    !Langevin antenna amplitudes must be reset to avoid errors during perf_opt
    if (antenna_type.eq.2.or.antenna_type.eq.3) call initialize_antenna
    init_status = 1

  end subroutine initialize_field_solver

  subroutine compute_emfields(g_1_in,emfields_out)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in):: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(out):: emfields_out
    logical,parameter :: OUTPUT=.false.
    real :: local_sum, global_sum
   
    if (.not.del_fields) then 
       call calc_charge_curr_dens(g_1_in,moments) 
       IF (OUTPUT) THEN
          CALL calculate_test_sum(moments(:,:,:,1),local_sum, global_sum)
          IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "phi = ",global_sum
          CALL calculate_test_sum(moments(:,:,:,2),local_sum, global_sum)
          IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "apar = ",global_sum
       END IF


       if (adiabatic_electrons) call add_adiabatic(moments(:,:,:,1))
       if (trap_pass) call add_adiabatic_passing_electrons(moments(:,:,:,1))

       if(xy_local) then
          call field_solve_ff(moments,emfields_out)
       else 
          call field_solve_df(moments,emfields_out)
       end if
    else
       !delete all electromagnetic fields
       emfields_out = 0.0
    end if

  end subroutine compute_emfields

  subroutine finalize_field_solver

    call finalize_charge_curr_dens
    deallocate(moments)
    if (adiabatic_electrons.or.trap_pass) call finalize_adiabatic_response
    if(xy_local) then
       call finalize_gyro_average_ff
       call finalize_field_solve_ff
    else
       call finalize_gyro_average_df
       call finalize_field_solve_df
    end if

    init_status = 0

  end subroutine finalize_field_solver

  ! routines for adaptive grids
  subroutine compute_emfields_adptv(g_1_in,emfields_out)
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2), intent(in):: g_1_in
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields), intent(out):: emfields_out
    logical,parameter :: OUTPUT=.false.
    real :: local_sum, global_sum
    
    call calc_charge_curr_dens_adptv(g_1_in,moments) 
    IF (OUTPUT) THEN
       CALL calculate_test_sum(moments(:,:,:,1),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "phi = ",global_sum
       CALL calculate_test_sum(moments(:,:,:,2),local_sum, global_sum)
       IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") "apar = ",global_sum
    END IF


    if (adiabatic_electrons) call add_adiabatic(moments(:,:,:,1))
    if (trap_pass) call add_adiabatic_passing_electrons(moments(:,:,:,1))

    if(xy_local) then
       stop "no adaptive grids for xy local version"
    else 
       call field_solve_df(moments,emfields_out)
    end if

  end subroutine compute_emfields_adptv

end module field_solver
