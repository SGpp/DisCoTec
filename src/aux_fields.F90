#include "redef.h"
#include "switches.h"
#include "intrinsic_sizes.h"
module aux_fields
  
  use par_mod
  use fieldsolver_ff
  use fieldsolver_df
  use field_solver
  use gyro_average_df_mod, only: gyro_average_df_wb
  use gyro_average_ff_mod, only: gyro_average_ff, I1_factor
  USE prefactors, only: mu_tjqj
  use boundaries
  use charge_curr_dens
  use compute_f
  use geometry, only: geom
  use fourier

  implicit none
  public:: initialize_calc_aux_fields, calc_aux_fields, finalize_calc_aux_fields,  calc_dxi_dz, &
       mem_est_aux_fields, calc_aux_fields_from_f
  PUBLIC:: f_, h_, emfields
  ! routines for adaptive grids
  public:: calc_aux_fields_adptv, antenna_contrib
  private
 
  COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE,TARGET :: emfields
  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: phi0
  Integer:: init_status = 0
 
contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_aux_fields(mem_req_in)
    real:: mem_req_in
    real:: mem_loc

    mem_loc=0.
    !em_fields
    mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*lz0*n_fields
   
    mem_loc=mem_est_field_solver(mem_loc)
    mem_loc=mem_est_compute_f(mem_loc)
    !gy_xi
    mem_loc=mem_loc+SIZE_OF_COMPLEX_MB*lij0*lz0
    mem_est_aux_fields=mem_req_in+mem_loc
  End Function mem_est_aux_fields


  subroutine initialize_calc_aux_fields

    ! allocate electromagnetic fields
    if (init_status.eq.0) then
       allocate(emfields(li1:li2,lj1:lj2,lbz:ubz,1:n_fields))
       if (.not.y_local) call initialize_fourier_boundary
    end if

    emfields=0.

    call initialize_field_solver   
    
    !initialize parallel exchange routine for f_
    !(also sets f_ and h_)
    call initialize_calc_f

    init_status = 1

  end subroutine initialize_calc_aux_fields

  !>computes the electromagnetic fields and auxiliary distribution functions
  !! \param  g_in input g_1 type array
  !! \param  emfields_out output electrostatic(magnetic) fields
  !! \param  f_out output f_1 type distribution function
  !! \param  comp_h optional switch to compute h_1 type distribution function
  !! \param  h_out optional output h_1 type distribution function
  !!               (only needed for arakawa_zv parallel scheme)
  subroutine calc_aux_fields(g_in,emfields_out,f_out,comp_h,h_out)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: g_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(out) :: emfields_out
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: f_out
    logical,intent(in):: comp_h
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),optional,intent(out):: h_out

    PERFON('calcaux')

    !switch off electromagnetic fields if only the L_g part of the linear 
    !operator is to be computed
    if(.not.precond_approx) then
       call compute_emfields(g_in,emfields_out)
    endif

    if (xy_local) then
       call exchange_z(emfields_out(:,:,:,1))
       if (n_fields .gt. 2) call exchange_z(emfields_out(:,:,:,3))
    endif

    call calc_f(g_in,emfields_out,f_out)

    !in case of comp_h=arakawa_zv, h_ is allocated.
    if (comp_h) call calc_h_from_f(f_out,emfields_out,h_out)
        
    PERFOFF

  end subroutine calc_aux_fields

  !Warning: this routine expects f without boundaries and gives g with boundaries, the inverse
  !of the usual setup! Currently only used for current_sheet initial condition, see init_cond.F90 
  subroutine calc_aux_fields_from_f(f_in,emfields_out,g_1_out)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: f_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: g_1_out
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(out) :: emfields_out

    PERFON('calcrest')

    if (xy_local) then
 
       !switch off electromagnetic fields if only the L_g part of the linear 
       !operator is to be computed
       if(.not.precond_approx) then
          call field_solve_ff_f(f_in,emfields_out)
       endif

       CALL exchange_z(emfields_out(:,:,:,1))
       IF (n_fields .GT. 2) CALL exchange_z(emfields_out(:,:,:,3))
    else
       stop "field_solve from f is not implemented for global simulations"
    endif
   
    call calc_f(f_in,-emfields_out,g_1_out)
    PERFOFF

  end subroutine calc_aux_fields_from_f


  subroutine finalize_calc_aux_fields

    deallocate(emfields)
    call finalize_field_solver
    call finalize_calc_f
    if(.not.y_local) call finalize_fourier_boundary
    init_status = 0

  end subroutine finalize_calc_aux_fields


  subroutine calc_dxi_dz(p_emfields,p_dxi_dz)
    complex, dimension(li1:li2, lj1:lj2, lbz:ubz, 1:n_fields),intent(inout):: p_emfields
    complex, dimension(li1:li2, lj1:lj2, lk1:lk2, lm1:lm2,ln1:ln2), intent(out):: p_dxi_dz

    complex, dimension(li1:li2,lj1:lj2,lbz:ubz):: gy_xi
    integer:: n, m, sten

    PERFON('dxi_dz')

    if (xy_local) then
       call exchange_z(p_emfields(:,:,:,1))
       IF (n_fields .GT. 2) CALL exchange_z(p_emfields(:,:,:,3))
    endif


    do n=ln1,ln2
       do m=lm1,lm2 
          if (xy_local) then
             call gyro_average_ff(p_emfields(:,:,:,1),gy_xi(li1:li2,:,:),m,n)
             IF (n_fields .GT. 2) CALL add_bparcontr(p_emfields(:,:,:,3),gy_xi,m,n)
          else
             call gyro_average_df_wb(p_emfields(:,:,:,1),gy_xi(:,:,:),m,n)
             !bpar cannot be used in the global version at the moment
             !IF (n_fields .GT. 2) CALL add_bparcontr(p_emfields(:,:,:,3),gy_xi,m,n)
             call exchange_z(gy_xi)
          endif

          sten = -par_sten_bound  
          if (shifted_metric) then
             p_dxi_dz(:,:,:,m,n) = par_sten(sten)*gy_xi(:,:,lk1+sten:lk2+sten)*&
                  &geom%phasefac(pi1:pi2,:,:,sten)
             do sten=-par_sten_bound+1,par_sten_bound
                p_dxi_dz(:,:,:,m,n) = p_dxi_dz(:,:,:,m,n)+par_sten(sten)*&
                     &gy_xi(:,:,lk1+sten:lk2+sten)*&
                     &geom%phasefac(pi1:pi2,:,:,sten)
             enddo
          else
             p_dxi_dz(:,:,:,m,n) = par_sten(sten)*gy_xi(:,:,lk1+sten:lk2+sten)
             do sten=-par_sten_bound+1,par_sten_bound
                call saxpy(2*lij0*lk0,par_sten(sten),gy_xi(li1,lj1,lk1+sten),1,p_dxi_dz(li1,lj1,lk1,m,n),1)
             enddo
          endif
       enddo
    enddo
    PERFOFF

  end subroutine calc_dxi_dz

  subroutine add_bparcontr(p_bpar,p_gy_xi,m,n)
    INTEGER, INTENT(in) :: m, n
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz), intent(in) :: p_bpar
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz), intent(inout) :: p_gy_xi

    p_gy_xi = p_gy_xi + mu_tjqj(m,n) * p_bpar * I1_factor(:,:,:,m,n)
  end subroutine add_bparcontr

  subroutine calc_aux_fields_adptv(g_in,emfields_out,f_out,comp_h,h_out)
    complex,dimension(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2),intent(in):: g_in
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,1:n_fields),intent(out) :: emfields_out
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(out):: f_out
    logical,intent(in):: comp_h
    complex,dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),optional,intent(out):: h_out

    PERFON('calcaux')

    !switch off electromagnetic fields if only the L_g part of the linear 
    !operator is to be computed
    if(.not.precond_approx) then
       call compute_emfields_adptv(g_in,emfields_out)
    endif

    if (xy_local) then
       stop "no adaptive grids for the xy local version"
    endif

    call calc_f_adptv(g_in,emfields_out,f_out)

    !in case of comp_h=arakawa_zv, h_ is allocated.
    if (comp_h) call calc_h_from_f_adptv(f_out,emfields_out,h_out)
        
    PERFOFF

  end subroutine calc_aux_fields_adptv

end module aux_fields
