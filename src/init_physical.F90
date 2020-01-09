#include "redef.h"
#include "intrinsic_sizes.h"
module phys_ini
  use par_mod
  use par_in, only: arakawa_zv
  use prefactors
  use checkpoint, only: initialize_checkpoint_module, finalize_checkpoint_module
  use communications
  use discretization_adptv_module
  use discretization
  use aux_func
  use geometry
  use collisions
  use vel_space
  use spectype_mod
  use boundaries
  use numerical_damping
  USE x_derivatives, only: initialize_radscheme, finalize_radscheme,&
       &CENTERED4TH
  use external_contr
  use equilibrium_fields
  use arrays
  use aux_fields
  use box_data_module
  use Grid1DModule
  use blockindex
  use dzv_terms, only: arakawa_zv_order
  use profiles_mod, only: set_spec_profiles, unset_spec_profiles, &
       &calc_derived_ref_vals, calc_dpdx_pm
  use lagrange_interpolation


  Implicit None
  public:: init_physical, finalize_physical, initialize_current_parall, finalize_current_parall,&
       mem_est_init_physical, initialize_boundaries, finalize_boundaries
  private

  Logical :: write_pe
  TYPE(box_data_type),save :: simulation_box

contains

  !>Give an estimate of the memory requirements of this module
  Real Function mem_est_init_physical(mem_req_in)
    real:: mem_req_in
    real:: mem_loc

    mem_loc=0.
    mem_loc=mem_est_geometry(mem_loc)
    mem_loc=mem_est_equilibrium_fields(mem_loc)
    mem_loc=mem_est_external_contr(mem_loc)
    mem_loc=mem_est_vel_space(mem_loc)
    mem_loc=mem_est_prefactors(mem_loc)
    if (collision_op.ne.'none') mem_loc=mem_est_collisions(mem_loc)

    mem_est_init_physical=mem_req_in+mem_loc
  End Function mem_est_init_physical
  

  Subroutine init_physical 

    !Local variables
    INTEGER :: i, j, k, l, n, pes, myspecrank
    Real :: c2coeff(-1:1),c4coeff(-2:2),c6coeff(-3:3),u3coeff(-2:2)


    !**************************************************************************!
    !**************************** Initialisation  *****************************!
    !**************************************************************************!
    
    write_pe = ((mype==0).AND.(print_ini_msg))

    call initialize_geometry
    call initialize_boundaries

    !note: deli needs to be initialized
    if (.not.xy_local) call initialize_radscheme(CENTERED4TH,deli,rad_bc_type)

    p_has_00_mode=.false.
    p_has_0_mode=.false.
   
    !the "0_mode" is the 0 Fourier mode of the spatial direction (x or y) in which 
    !only positive modes are considered, due to the reality constraint
    do j=lj1,lj2
       do i=li1,li2
          !Nyquist modes don't count
          if (evenx.eq.1) then
             if (yx_order) then
                if (j.eq.hkx+1) cycle
             else
                if (i.eq.hkx+1) cycle
             end if
          end if
          if (yx_order.and.xy_local) then
             !in this case the relevant mode is ky=ki=0
             if ((abs(ki(i)).lt.epsilon(ki(i)))) then
                p_has_0_mode=.true.
                if (abs(kj(j)).lt.epsilon(kj(j))) &
                     p_has_00_mode=.true.
             end if
          else
             !in this case the relevant mode is kj=0
             if ((abs(kj(j)).lt.epsilon(kj(j)))) then
                p_has_0_mode=.true.
                if (xy_local.and.(abs(ki(i)).lt.epsilon(ki(i)))) &
                     p_has_00_mode=.true.
             end if
          endif
       enddo
    end do

    ! Set T, n, grad T, grad n profiles (after lx adaption!)
    call set_spec_profiles(spec,xval_a,n_spec,pi1gl,pi2gl,write_pe)
    call calc_derived_ref_vals(beta,coll,debye2,rhostar,Omega0_tor,write_pe)
    !check parameters now?

    !calculate the negative pressure gradient from n and T gradients 
    !if requested (dpdx_pm = -1)
    if (dpdx_pm.eq.-1) then 
       call calc_dpdx_pm(dpdx_pm_arr)
    elseif (dpdx_pm.gt.0) then
       dpdx_pm_arr = dpdx_pm
    endif

    if (x_local) dpdx_pm = dpdx_pm_arr((pi1gl+pi2gl)/2)
   
    !Initialize external (radial) fields 
    !(must be called after geometry, before prefactor!)
    CALL init_external_contributions
    CALL initialize_equilibrium_fields

    CALL initialize_vel_space ! initialize velocity space quantities
    CALL initialize_prefactors
    
    !stencils for finite differences    
    if (write_pe) then
       if (prec.EQ.'DOUBLE') then
          write(*,"(A)") "double precision computation"
       else
          write(*,"(A)") "single precision computation"
       endif
    endif

    !coefficients for 1st order derivatives
    c2coeff=(/-1,0,1/)/(2.0)
    u3coeff=(/1,-6,3,2,0/)/(6.0)
    c4coeff=(/1,-8,0,8,-1/)/(12.0)
    c6coeff=(/-1,9,-45,0,45,-9,1/)/(60.0)

    np1coeff=(/0, 0, -3, 4, -1/)/(2.0*dz)
    np2coeff=(/0, -2, -3, 6, -1/)/(6.0*dz)
    np3coeff=(/1, -6, 3, 2, 0/)/(6.0*dz)
    np4coeff=(/1, -4, 3, 0, 0/)/(2.0*dz)

    allocate(vderivative(-2:2,pi1:pi2,lk1:lk2,ll1:ll2))
    do l=ll1,ll2
       do k=lk1,lk2
          do i=pi1,pi2
             select case(vparscheme)
             case('c4th')
                vderivative(:,i,k,l)=c4coeff/dv
#ifndef oldvpbc
                if((l.eq.ll1).and.(my_pev.eq.0)) then
                   if (geom%dBdz(i,lj1,k).gt.0) then
                      vderivative(:,i,k,l)=(/0, 0, -3, 4, -1/)/(2.0*dv)
                   else
                      !vderivative(:,i,k,l)=(/0, -1, 1, 0, 0/)/dv !assumes zero outside vp box
                   endif
                elseif((l.eq.ll1+1).and.(my_pev.eq.0)) then
                   if (geom%dBdz(i,lj1,k).gt.0) then
                      vderivative(:,i,k,l)=(/0, -2, -3, 6, -1/)/(6.0*dv)
                   else
                      !vderivative(:,i,k,l)=(/0, -1, 0, 1, 0/)/(2.0*dv)
                   endif
                elseif((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) then
                   if (geom%dBdz(i,lj1,k).ge.0) then
                      !vderivative(:,i,k,l)=(/0, -1, 0, 1, 0/)/(2.0*dv)
                   else
                      vderivative(:,k,l)=u3coeff/dv
                   endif
                elseif((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) then
                   if (geom%dBdz(i,lj1,k).ge.0) then
                      !vderivative(:,i,k,l)=(/0, 0, -1, 1, 0/)/dv !assumes zero outside vp box
                   else
                      vderivative(:,i,k,l)=(/1, -4, 3, 0, 0/)/(2.0*dv)
                   endif
                endif
#endif
             case('u3rd') 
                if (geom%dBdz(i,lj1,k).ge.0) then
                   !reverse stencil due to opp. sign of prefactor
                   vderivative(:,i,k,l)=-u3coeff(5:1:-1)/dv 
                else
                   vderivative(:,i,k,l)=-u3coeff/dv 
                endif
             case default
                stop 'chosen vparscheme not implemented'
             end select
          end do
       end do
    enddo

    if ((.not.arakawa_zv).or.(parallel_nl)) then
       select case(parscheme)
       case('c0th')
          if (write_pe) write(*,"(A)") "parallel direction: derivatives zeroed"
          par_sten_bound=0
          !       allocate(par_sten(-par_sten_bound:par_sten_bound))
          allocate(par_sten(0:0))
          par_sten=0.0 ! c2coeff/dz
       case('c2nd')
          if (write_pe) write(*,"(A)") "parallel direction: 2nd order centered differences"
          par_sten_bound=1
          allocate(par_sten(-par_sten_bound:par_sten_bound))
          par_sten=c2coeff/dz
       case('c4th')
          if (write_pe) write(*,"(A)") "parallel direction: 4th order centered differences"
          par_sten_bound=2
          allocate(par_sten(-par_sten_bound:par_sten_bound))
          par_sten=c4coeff/dz
       case('u3rd')
          if (write_pe) write(*,"(A)") "parallel direction: 3rd order upwind scheme" 
          par_sten_bound=2
          allocate(par_sten(-par_sten_bound:par_sten_bound))
          par_sten=u3coeff/dz
       case('c6th')
          if (write_pe) write(*,"(A)") "parallel direction: 6th order centered differences"
          par_sten_bound=3
          allocate(par_sten(-par_sten_bound:par_sten_bound))
          par_sten=c6coeff/dz
       case default
          if (write_pe) write(*,"(A)") "no valid choice for parscheme"
          stop
       end select
    endif
    if (arakawa_zv) then
       if (arakawa_zv_order==2.and.write_pe) then
          write(*,"(A)") "parallel direction: 2nd order Arakawa scheme"
       elseif (arakawa_zv_order==4.and.write_pe) then
          write(*,"(A)") "parallel direction: 4th order Arakawa scheme"
       endif
    endif

    !calculate hyper diffusion prefactors
    call initialize_numerical_damping

    !recurrence time
    myspecrank = modulo(mype, n_procs_x*n_procs_y*n_procs_z*n_procs_v*n_procs_w)
    if((mype.eq.0).and.(print_ini_msg)) write(*,"(A)") 'recurrence times for the different species:'
    do pes=0,n_procs_s-1
       if (my_pespec.eq.pes) then
          do n = ln1,ln2
             if ((myspecrank == 0).and.(print_ini_msg)) then
                write(*,"(A,A,F10.5)") spec(n)%name,' : ', &
                     & 2.0*pi/(sqrt(2.*spec(n)%temp/spec(n)%mass)*dv)
             endif
          enddo
       endif
       call my_barrier
    enddo

    Call my_barrier()

    ! Initialisation of collision operator
    If (collision_op.ne.'none') call initialize_collisions

    !correct location for initialize_checkpoint_module?
    Call initialize_checkpoint_module

    Call my_barrier()

  end Subroutine init_physical


  subroutine initialize_boundaries
    TYPE(grid1d) :: igrid, kjgrid,zgrid,vgrid,wgrid
    !initialize boundary conditions (requires lx/kymin to be adapted to
    !boundary conditions!)

    if (xy_local) then
       call initialize_boundary_exchange(q0, shat)
    else
       CALL initialize(simulation_box, q0, shat, mag_prof, &
            &rhostar) !mapped to box_data_initialize

       CALL initialize(igrid,ni0,.FALSE.)
!       IF (rad_bc_type.eq.0) THEN
          ! periodic boundary condition in x implies periodic grid in x
          CALL set_boundaries(igrid,-lp1/2,lp1/2,(rad_bc_type.eq.0))
!       ELSE
          ! all other boundary conditions imply a nonperiodic grid in x
!          CALL set_boundaries(igrid,-lp1/2,lp1/2,.FALSE.)
!       END IF
       CALL initialize(kjgrid,nj0,.true.)
       CALL set_fourier_boundaries(kjgrid,kj0_ind*kjmin,kjmin)
       CALL initialize(zgrid,nz0,.false.)
       CALL set_boundaries(zgrid,-pi,pi,.true.)
       CALL initialize(vgrid,nv0,.false.)
#ifdef COMBI
       CALL set_boundaries(vgrid,-lv+SHIFT,lv+SHIFT,.false.)
#else
       CALL set_boundaries(vgrid,-lv,lv,.false.)
#endif
       CALL initialize(wgrid,nw0,.false.)
       CALL set_boundaries(wgrid,0.0,lw,.false.)
       CALL set_grids(simulation_box,igrid,kjgrid,zgrid,vgrid,wgrid)

       !IF (mype.EQ.0) CALL print_box(simulation_box)

       CALL initialize_boundary_exchange(q_prof,rad_bc_type,simulation_box)
       if (shifted_metric) CALL init_shifted_metric
    endif
  end subroutine initialize_boundaries

  subroutine finalize_boundaries
    call finalize_boundary_exchange
  end subroutine finalize_boundaries

  !> Initializes the phase factor that is applied to parallel derivatives with shifted_metric=.t.
  Subroutine init_shifted_metric
    integer:: i, j, k, sten
    integer:: nx_fine
    real,dimension(:,:),allocatable:: shift_int_fine, shift_int
    real,dimension(:),allocatable:: xval_a_fine, alpha, tmparr, tmparr_fine
    real:: ctr_shift, dx_fine

    call exchange_z_nopb(geom%shift)
    !apply shift of 2*pi*shat at the parallel ends of the domain
    if (my_pez==n_procs_z-1.or.n_procs_z==1) then
       do k=1,nzb
          if (lilo) then
             geom%shift(:,:,lk2+k)=geom%shift(:,:,lk2+k)+2*pi*shat
          else
             do j=pj1,pj2
                geom%shift(:,j,lk2+k)=geom%shift(:,j,lk2+k)+2*pi*x0/q0*dqdx_prof*minor_r
             enddo
          endif
       enddo
    endif
    if (my_pez==0.or.n_procs_z==1) then
       do k=1,nzb
          if (lilo) then
             geom%shift(:,:,lk1-k)=geom%shift(:,:,lk1-k)-2*pi*shat
          else
             do j=pj1,pj2
                geom%shift(:,j,lk1-k)=geom%shift(:,j,lk1-k)-2*pi*x0/q0*dqdx_prof*minor_r
             enddo
          endif
       enddo
    endif
 
    nx_fine=nx0*16
    allocate(shift_int_fine(1:nx_fine,lbz:ubz),xval_a_fine(1:nx_fine),&
         alpha(1:nx_fine),shift_int(1:nx0,lbz:ubz),tmparr(1:nx0),&
         tmparr_fine(1:nx_fine))

    dx_fine=(maxval(xval_a)-minval(xval_a))/(nx_fine-1)
    do i=1,nx_fine
       xval_a_fine(i)=minval(xval_a)+dx_fine*(i-1)
    enddo
    do k=lbz,ubz
       tmparr=real(geom%shift(:,pj1,k))
       call lag3interp(tmparr,xval_a,nx0,alpha,xval_a_fine,nx_fine)
       do i=1,nx_fine
          shift_int_fine(i,k)=sum(alpha(1:i))*dx_fine/rhostar
       enddo
       tmparr_fine=shift_int_fine(:,k)
       call lag3interp(tmparr_fine,xval_a_fine,nx_fine,ctr_shift,x0)
       tmparr_fine=tmparr_fine-ctr_shift
       call lag3interp(tmparr_fine,xval_a_fine,nx_fine,tmparr,xval_a,nx0)
       shift_int(:,k)=tmparr
    enddo
    do sten=-2,2
       do k=lk1,lk2
          do j=lj1,lj2
             geom%phasefac(pi1gl:pi2gl,j,k,sten)=exp(imag*ky(j)*(shift_int(1:nx0,k)-shift_int(1:nx0,k+sten)))
          enddo
       enddo
    enddo

    deallocate(shift_int_fine,shift_int,xval_a_fine,alpha,tmparr,tmparr_fine)

  End Subroutine init_shifted_metric


  Subroutine finalize_physical 
    IMPLICIT NONE

    deallocate(vderivative)

    call finalize_external_contributions
    call finalize_equilibrium_fields
    call finalize_geometry

    Call unset_spec_profiles(spec,n_spec)

    call finalize_boundaries

    CALL finalize_prefactors
    CALL finalize_vel_space

    if ((.not.arakawa_zv).or.(parallel_nl)) &
         &deallocate(par_sten)

    if (.not.xy_local) then
       call finalize_radscheme
    endif

    call finalize_numerical_damping

    If (collision_op.ne.'none') then
       call finalize_collisions
    endif

    !correct location for finalize_checkpoint_module?
    Call finalize_checkpoint_module

  end Subroutine finalize_physical


  subroutine initialize_current_parall

    call split_comm
    call initialize_discretization(print_ini_msg)
    call allocate_arrays
    call initialize_blockindex
    call init_physical

  end subroutine initialize_current_parall


  subroutine finalize_current_parall

    Call finalize_physical
    Call deallocate_arrays
    Call free_comm
    
  end subroutine finalize_current_parall


End module phys_ini
