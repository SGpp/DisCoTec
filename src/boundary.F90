#include "redef.h"
!> All boundary exchange functions for all dimensions.
!!
!! This module contains all functions used for implementing
!! the boundary conditions and furthermore the closely
!! related exchange functions for inner and outer
!! boundary points. Also the ghost cell exchange for parallelized
!! dimensions are implemented here
!!
Module boundaries
  Use par_mod
  Use mpi
  USE communications, ONLY: communicators, COMM_X, COMM_Z, COMM_V, COMM_W,&
       & N_PE_DIMS, comm_cart, mpi_comm_z, mpi_comm_w, mpi_comm_v
 
  !USE lagrange_interpolation, ONLY: lag_interp, lag_filtering
  USE boundary_exchange_general
  USE boundary_exchange_x
  USE boundary_exchange_z
  use boundary_exchange_vw
  USE BoundaryDescriptionModule
  use box_data_module
  use exchange_z_ff
  use blockindex, only: lbg0

  Implicit None
  
  PUBLIC :: exchange_x, exchange_z, exchange_z_nopb, exchange_v, exchange_mu
  PUBLIC :: exchange_z_noinner
  PUBLIC :: initialize_boundary_exchange, finalize_boundary_exchange

  PUBLIC :: initialize_exchange_z5d, finalize_exchange_z5d
  PUBLIC :: initialize_exchange_v, finalize_exchange_v
  PUBLIC :: ubexc, lx0_boundary, lx0_boundary_1D, lx0_boundary_2D, lx0_boundary_2D_real, pb_phase
  public :: x_boundary_block
  PUBLIC :: pb_xshift, delete_entry_lower, delete_entry_upper
  PRIVATE 

  INTERFACE exchange_z
     MODULE PROCEDURE exchange_3d, exchange_5df
  END INTERFACE

  INTERFACE exchange_z_nopb
     MODULE PROCEDURE exchange_3d_nopb
  END INTERFACE

  INTERFACE exchange_v
     MODULE PROCEDURE bnd_exchange_v
  END INTERFACE

  INTERFACE exchange_mu
     MODULE PROCEDURE bnd_exchange_mu
  END INTERFACE

  INTERFACE initialize_boundary_exchange
     MODULE PROCEDURE initialize_boundary_exchange_ff, initialize_boundary_exchange_df
  END INTERFACE


  INTEGER :: vectype_vexc

  TYPE(BoundaryDescription) :: bdesc_z_5df_pv1,bdesc_z_5df_pv2
  TYPE(BoundaryDescription) :: bdesc_z_3D !< Boundary description for 3D z exchange, needed for the fields
  TYPE(BoundaryDescription) :: bdesc_z_3D_metric
  TYPE(BoundaryDescription) :: bdesc_mu, bdesc_v_1, bdesc_v_2
  TYPE(BoundaryDescription) :: lx0_boundary, lx0_boundary_1D, lx0_boundary_2D, lx0_boundary_2D_real
  TYPE(BoundaryDescription) :: x_boundary_block

#ifdef WITHOMP_BLOCKLOOP
  !$OMP THREADPRIVATE(lx0_boundary,lx0_boundary_1D,lx0_boundary_2D,lx0_boundary_2D_real,x_boundary_block)
#endif
  Integer :: init_status_exchange_z5d = 0
  Integer :: init_status_exchange_v = 0

Contains



  SUBROUTINE initialize_boundary_exchange_ff(ai_q0, ai_shat)
    REAL, INTENT(IN) :: ai_q0, ai_shat

    ! Local variables
    INTEGER, DIMENSION(N_PE_DIMS) :: pe_dims, my_pe_coords
    LOGICAL, dimension(N_PE_DIMS) :: periods
    integer :: ierr

    CALL mpi_cart_get(comm_cart,N_PE_DIMS, pe_dims, periods, my_pe_coords, ierr)
    
    n_procs_w = pe_dims(2)
    n_procs_v = pe_dims(3)
    n_procs_z = pe_dims(4)
    n_procs_x = pe_dims(6)
    my_pew    = my_pe_coords(2)
    my_pev    = my_pe_coords(3)
    my_pez    = my_pe_coords(4)
    my_pex    = my_pe_coords(6)

    ! call the specialized initialization routines for the different
    ! directions
    CALL initialize_boundary_exchange_general(comm_cart)
    call init_parallel_boundary_ff(ai_q0,ai_shat)

  END SUBROUTINE initialize_boundary_exchange_ff
  
  SUBROUTINE initialize_boundary_exchange_df(ai_qprof, ai_rad_bc_type, &
       &ai_simulation_box)
    REAL, DIMENSION(pi1gl:pi2gl), INTENT(IN) :: ai_qprof
    INTEGER,INTENT(IN) :: ai_rad_bc_type
    TYPE(box_data_type), INTENT(IN) :: ai_simulation_box

    ! Local variables
    INTEGER, DIMENSION(N_PE_DIMS) :: pe_dims, my_pe_coords
    LOGICAL, dimension(N_PE_DIMS) :: periods
    integer :: ierr

    CALL mpi_cart_get(comm_cart,N_PE_DIMS, pe_dims, periods, my_pe_coords, ierr)
    
    n_procs_w = pe_dims(2)
    n_procs_v = pe_dims(3)
    n_procs_z = pe_dims(4)
    n_procs_x = pe_dims(6)
    my_pew    = my_pe_coords(2)
    my_pev    = my_pe_coords(3)
    my_pez    = my_pe_coords(4)
    my_pex    = my_pe_coords(6)

    ! call the specialized initialization routines for the different
    ! directions
    CALL initialize_boundary_exchange_general(comm_cart)
    CALL initialize_boundary_exchange_x(communicators(COMM_X), ai_rad_bc_type)

#ifdef WITHOMP_BLOCKLOOP
    !$OMP PARALLEL default(none) &
    !$OMP shared(lx0,nib,lj0, li0da,ly0da,n_procs_x,lbg0)
#endif
    !CALL initialize_type(lx0_boundary_1D_y,1,li0da/n_procs_y+2*nib,nib,nib,1)
    CALL initialize_type(lx0_boundary,1,lx0,nib,nib,lj0)
    CALL initialize_type(lx0_boundary_2D,1,li0da+2*nib,nib,nib,ly0da)
    CALL initialize_type(lx0_boundary_2D_real,1,li0da+2*nib,nib,nib,ly0da)
    CALL initialize_type(lx0_boundary_1D,1,lx0,nib,nib,1)
    call initialize_type(x_boundary_block,1,lx0,nib,nib,lj0*lbg0)
    
    IF (n_procs_x > 1) THEN
       CALL set_mpi_type(lx0_boundary)
       CALL set_mpi_type(lx0_boundary_1D)
       CALL set_mpi_type(lx0_boundary_2D)
       CALL set_mpi_type_real(lx0_boundary_2D_real)
       call set_mpi_type(x_boundary_block)
    ENDIF
#ifdef WITHOMP_BLOCKLOOP
    !$OMP END PARALLEL
#endif

    IF (y_local) THEN
       CALL initialize_boundary_exchange_z(communicators(COMM_Z),ai_simulation_box, ai_qprof)

       ! set the boundary descriptions
       CALL initialize_type(bdesc_z_5df_pv1,lij0,lz0,nzb,nzb,ll0)
       CALL initialize_type(bdesc_z_5df_pv2,lij0,lz0,nzb,nzb,lv0*lw0*ln0)
       CALL initialize_type(bdesc_z_3D,     lij0,lz0,nzb,nzb,1)
       CALL initialize_type(bdesc_z_3D_metric,px0*pj0,lz0,nzb,nzb,1)
       !hkd: add besc for checkpoint_read here
       !don't know if I have to set mpi type, if I want each z proc exchanging its own array.
       
       IF (n_procs_z.GT.1) THEN
          CALL set_mpi_type(bdesc_z_5df_pv1)
          CALL set_mpi_type(bdesc_z_5df_pv2)
          CALL set_mpi_type(bdesc_z_3D)
          CALL set_mpi_type(bdesc_z_3D_metric)
       END IF
    ELSE
       call init_parallel_boundary_ff(get_q0(ai_simulation_box),get_shat(ai_simulation_box))       
    ENDIF

  END SUBROUTINE initialize_boundary_exchange_df

  SUBROUTINE finalize_boundary_exchange
    
    IF (x_local) THEN
       call finalize_parallel_boundary_ff
    ELSE
#ifdef WITHOMP_BLOCKLOOP
       !$OMP PARALLEL
#endif
       CALL finalize_type(lx0_boundary)
       CALL finalize_type(lx0_boundary_1D)
       CALL finalize_type(lx0_boundary_2D)
       CALL finalize_type(lx0_boundary_2D_real)
       call finalize_type(x_boundary_block)
#ifdef WITHOMP_BLOCKLOOP
       !$OMP END PARALLEL
#endif
       CALL finalize_boundary_exchange_x
       CALL finalize_type(bdesc_z_5df_pv1)
       call finalize_type(bdesc_z_5df_pv2)
       call finalize_type(bdesc_z_3D)
       call finalize_type(bdesc_z_3D_metric)
       call finalize_boundary_exchange_z
    END IF  
    CALL finalize_boundary_exchange_general
      
    !CALL finalize_boundary_exchange_vw
    !CALL finalize_boundary_exchange_w
  END SUBROUTINE finalize_boundary_exchange

  !**************************************************!
  !************ Exchange in velocity direction ******!
  !**************************************************!
  SUBROUTINE initialize_exchange_v

    if (init_status_exchange_v==1) return

    CALL initialize_boundary_exchange_vw(communicators(COMM_V),communicators(COMM_W))
    CALL initialize_type(bdesc_mu,lij0*lz0*lv0,lw0,nwb,nwb,ln0)
    if (n_procs_w.gt.1) call set_mpi_type(bdesc_mu)

    CALL initialize_type(bdesc_v_1,lij0*lz0,lv0,nvb,nvb,1)
    CALL initialize_type(bdesc_v_2,lij0*lz0,lv0,nvb,nvb,lw0*ln0)
    IF (n_procs_v.GT.1) THEN
       CALL set_mpi_type(bdesc_v_1)
       CALL set_mpi_type(bdesc_v_2)
    END IF

    init_status_exchange_v = 1

  END SUBROUTINE initialize_exchange_v

  SUBROUTINE finalize_exchange_v
    call finalize_type(bdesc_mu)
    CALL finalize_type(bdesc_v_1)
    CALL finalize_type(bdesc_v_2)
    call finalize_boundary_exchange_vw
    init_status_exchange_v = 0
  END SUBROUTINE finalize_exchange_v

  SUBROUTINE bnd_exchange_v(p_f)
    complex, dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),intent(inout):: p_f

    if (perf_vec(8).eq.1) then
       CALL exchange_v(bdesc_v_1,p_f)
    else
       CALL exchange_v(bdesc_v_2,p_f)
    endif
  END SUBROUTINE bnd_exchange_v

  SUBROUTINE bnd_exchange_mu(p_f )
    complex,Intent(INOUT):: p_f(li1:li2, lj1:lj2, lbz:ubz, lbv:ubv, lbw:ubw, ln1:ln2)

    CALL exchange_mu(bdesc_mu,p_f)
  END SUBROUTINE bnd_exchange_mu
  !**************************************************!
  !************ Exchange in parallel direction ******!
  !**************************************************!

  !this routine assumes that the input array u is not distributed in z
  !no mpi exchanges needed
  !(for example in read-in of chpt with changing z resolution)
  subroutine exchange_z_noinner(u)
    complex, dimension(li1:li2,lj1:lj2, lbz:ubz), intent(inout):: u
    IF (x_local) THEN
       !supposed to use n_procs_z = 1 branch of the xy_local routine
       CALL parallel_boundary_ff(u,1,0)
    ELSE
       !extra routine for global version
       CALL exchange_z_df_noinner(u)
    ENDIF
  end subroutine exchange_z_noinner

  Subroutine exchange_3d(u)
    complex, Dimension(li1:li2, lj1:lj2, lbz:ubz), Intent(Inout):: u
    PERFON_I('z_ex_3d')

    IF (x_local) THEN
       CALL parallel_boundary_ff(u,1,0)
    ELSE
       CALL exchange_z_df(bdesc_z_3d,u)
    END IF

    PERFOFF_I
  End Subroutine exchange_3d

  Subroutine exchange_3d_nopb(u)
    complex, Dimension(pi1gl:pi2gl,pj1:pj2,lbz:ubz), Intent(Inout):: u
    !PERFON('z_ex_3d_metric')

    CALL exchange_z_nopb(bdesc_z_3D_metric,u)
    
    !PERFOFF
  End Subroutine exchange_3d_nopb
  
  subroutine initialize_exchange_z5d

    if ((init_status_exchange_z5d.ne.perf_vec(7)).and. &
         & (init_status_exchange_z5d.gt.0)) &
         & call finalize_exchange_z5d

    if (perf_vec(7).eq.1) then
       !nothing
!       call initialize_exchange_5df_ff_1
    else
       if (n_procs_z.gt.1) call initialize_exch_z_cyc
    endif

    init_status_exchange_z5d = perf_vec(7)

  end subroutine initialize_exchange_z5d

  subroutine exchange_5df(u)
    complex, Dimension(li1:li2,lj1:lj2,lbz:ubz,lbv:ubv,lbw:ubw,ln1:ln2),Intent(Inout):: u
    integer:: m,n
    LOGICAL :: OUTPUT=.false.

    PERFON_I('z_ex_5df')
 
    if (x_local) then
       if (perf_vec(7).eq.1) then
          call exchange_5df_ff_1(u)
       else
          call exchange_5df_ff_2(u)
       endif
    else 
       if (perf_vec(7).eq.1) then
          IF (OUTPUT) WRITE(*,"(I3,A,2ES20.10)") mype,": -deb- before1 u = ",&
               &DBLE(SUM(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)))),&
               &DBLE(SUM(u(:,:,:,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(:,:,:,ll1:ll2,lm1:lm2,ln1:ln2))))
          do n=ln1,ln2
             do m=lm1,lm2
                CALL exchange_z_df(bdesc_z_5df_pv1,u(:,:,:,ll1:ll2,m,n))
             enddo
          enddo
          IF (OUTPUT) WRITE(*,"(I3,A,2ES20.10)") mype,": -deb- after1  u = ",&
               &DBLE(SUM(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)))),&
               &DBLE(SUM(u(li1:li2,:,:,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,:,:,ll1:ll2,lm1:lm2,ln1:ln2))))
       ELSE
          IF (OUTPUT) WRITE(*,"(I3,A,2ES20.10)") mype,": -deb- before2 u = ",&
               &DBLE(SUM(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)))),&
               &DBLE(SUM(u(:,:,:,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(:,:,:,ll1:ll2,lm1:lm2,ln1:ln2))))
          CALL exchange_z_df(bdesc_z_5df_pv2,u)
          IF (OUTPUT) WRITE(*,"(I3,A,2ES20.10)") mype,": -deb- after2  u = ",&
               &DBLE(SUM(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,lj1:lj2,lk1:lk2,ll1:ll2,lm1:lm2,ln1:ln2)))),&
               &DBLE(SUM(u(li1:li2,:,:,ll1:ll2,lm1:lm2,ln1:ln2)&
               &*CONJG(u(li1:li2,:,:,ll1:ll2,lm1:lm2,ln1:ln2))))
       endif
    endif
    
    PERFOFF_I

  end subroutine exchange_5df

  subroutine finalize_exchange_z5d

    if (init_status_exchange_z5d.eq.1) then
!       call finalize_exchange_5df_ff_1
    else
       if (n_procs_z.gt.1) call finalize_exch_z_cyc
    endif

    init_status_exchange_z5d = 0
  end subroutine finalize_exchange_z5d

  
End Module boundaries
