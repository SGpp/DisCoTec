#include "redef.h"
#include "intrinsic_sizes.h"
module numerical_damping
  use discretization
  use coordinates
  use geometry
  use par_other, only: print_ini_msg
  use communications
  use par_in

  implicit none
  public:: initialize_numerical_damping, finalize_numerical_damping, set_nd_defaults
  public:: hyp_x, hyp_y, hyp_perp, hyp_z, hyp_v, hyp_x_order, hyp_y_order, hyp_perp_order, hyp_z_order
  public:: damp_i, damp_j, damp_perp, didiff, djdiff, dperpdiff
  public:: didiff_max, djdiff_max, dperpdiff_max
  public:: dc2coeff, dc4coeff, hyp_z_sten, vdamping
  public:: hyp_i, hyp_j, hyp_i_order, hyp_j_order, hyp_z_scl
  real, dimension(:), allocatable, public :: hyp_x_spec, hyp_y_spec
  private

  real:: hyp_x, hyp_y, hyp_perp, hyp_z, hyp_v
  integer:: hyp_x_order, hyp_y_order, hyp_perp_order, hyp_z_order
  real:: hyp_i, hyp_j, hyp_z_scl
  integer:: hyp_i_order, hyp_j_order

  logical:: damp_i, damp_j, damp_perp
  real :: dc2coeff(-1:1), dc4coeff(-2:2), dc6coeff(-3:3)
  real, dimension(:), allocatable:: didiff, djdiff, hyp_z_sten
  real, dimension(:,:,:), allocatable:: dperpdiff
  real, dimension(:,:),allocatable:: vdamping
  real :: didiff_max, djdiff_max, dperpdiff_max
  logical:: write_pe

contains

  subroutine set_nd_defaults

    hyp_x=0.0 ; hyp_y=0.0 ; hyp_perp=0.0 ; hyp_z=0.0 ; hyp_v=0.0
    hyp_x_order=4 ; hyp_y_order=4 ; hyp_perp_order=2 ; hyp_z_order=4
  end subroutine set_nd_defaults

  subroutine initialize_numerical_damping
    real :: k_perp2
    integer:: i,j,k,l,ierr

    write_pe = ((mype==0).AND.(print_ini_msg))
    if (write_pe) write(*,"(A)") "=========== hyper diffusion ============"

    !coefficients for 2nd, 4th and 6th derivative
    dc2coeff=(/1,-2,1/)*0.25 
    dc4coeff=(/-1,4,-6,4,-1/)* 0.0625
    dc6coeff=(/1,-6,15,-20,15,-6,1/)* 0.015625


    !x/y direction
    if (xy_local.and..not.nonlinear.and.nky0.gt.1) then
       if (write_pe) write(*,"(A)") "Switching off radial hyperdiffusion for linear run with nky0>1"
       hyp_x=0.
    end if

    if (hyp_x.lt.0.0) hyp_x = - hyp_x
    if (hyp_y.lt.0.0) hyp_y = - hyp_y
    if ((.not.nonlinear).and.(hyp_y.gt.0).and.(y_local)) &
         &stop "ERROR: Using hyp_y in linear runs not allowed"
    if (hyp_perp.lt.0.0) hyp_perp = - hyp_perp
    if (write_pe) then
       if (hyp_x.gt.0.0) write(*,"(A,F6.3)") "hyp_x =     ",hyp_x
       if (hyp_y.gt.0.0) write(*,"(A,F6.3)") "hyp_y =     ",hyp_y
       if (hyp_perp.gt.0.0) write(*,"(A,F6.3)") "hyp_perp =  ",hyp_perp
    end if

    if (yx_order) then
       hyp_i=hyp_y; hyp_j=hyp_x
       hyp_i_order=hyp_y_order; hyp_j_order=hyp_x_order
       GyroLES=.false. !GyroLES doest work with yx order
       diag_GyroLES=.false. !diagnostics for GyroLES doest work with yx order
    else
       hyp_i=hyp_x; hyp_j=hyp_y
       hyp_i_order=hyp_x_order; hyp_j_order=hyp_y_order
    end if

    if (hyp_i.gt.0.) then
       damp_i=.true.
       if (xy_local) then
          !i is in Fourier space
          allocate(didiff(li1:li2))
          didiff = -(1/2.*deli*ki(:))**hyp_i_order
       else
          allocate(didiff(-hyp_i_order/2:hyp_i_order/2))
          select case(hyp_i_order)
          case(6)
             didiff(:) = dc6coeff
          case(4)
             didiff(:) = dc4coeff
          case(2)
             didiff(:) = dc2coeff
          case default
             hyp_i=0.0
          end select
       end if
       Call mpi_allreduce (maxval(abs(didiff)),didiff_max,1,&
            MPI_REAL_TYPE, MPI_MAX,mpi_comm_x, ierr)         
    else
       damp_i=.false.
    end if

    if (hyp_j.gt.0.0) then
       damp_j=.true.
       allocate(djdiff(lj1:lj2))
       djdiff = -(1/2.*delj*kj(lj1:lj2))**hyp_j_order
       Call mpi_allreduce (maxval(abs(djdiff)),djdiff_max,1,&
            MPI_REAL_TYPE, MPI_MAX,mpi_comm_y, ierr)        
    end if

    !perpendicular hyperdiffusion   
    if (hyp_perp.gt.0) then
       if (xy_local) then
          damp_perp=.true.
          allocate(dperpdiff(li1:li2,lj1:lj2,lk1:lk2))
          if (hyp_perp.gt.0.0) damp_perp=.true.
          do k=lk1,lk2
             do j=lj1,lj2
                do i=li1,li2
                   k_perp2=geom%gii(pi1,pj1,k)*ki(i)**2+&
                        2*geom%gij(pi1,pj1,k)*ki(i)*kj(j)+&
                        geom%gjj(pi1,pj1,k)*kj(j)**2
                   dperpdiff(i,j,k)=k_perp2**(hyp_perp_order/2.)
                end do
             enddo
          enddo
          Call mpi_allreduce (maxval(abs(dperpdiff)),dperpdiff_max,1,&
               MPI_REAL_TYPE, MPI_MAX,mpi_comm_xyz, ierr)   
       else
          stop 'hyp_perp not implemented for nonlocal code versions'
       endif
    endif

    !parallel direction
    allocate(hyp_z_sten(-hyp_z_order/2:hyp_z_order/2))
    !if (hyp_z .lt. 0.0) hyp_z = - hyp_z
    if (hyp_z .lt. 0.0) then
       if (write_pe) write(*,"(A)") "prefactor for z-hyperdiffusion is scaled with dz"
       !introduce scaling factor to hyperdiffusion stencil:
       !if prefactor to parallel advection is (-)unity, c4th+hyp_z yields the upwind stencil for hyp_z = -1 
       !the sign compensates for hyp_z being negative
       hyp_z_scl = -hyp_z*4.0/(3.0*dz)
    else
       !hyperdiffusion stencil is not scaled with dz
       hyp_z_scl = hyp_z
    endif

    select case(hyp_z_order)
    case(6)
       hyp_z_sten = dc6coeff
    case(4)
       hyp_z_sten = dc4coeff
    case(2)
       hyp_z_sten = dc2coeff
    case default
       hyp_z=0.0
       if (write_pe) write(*,"(A)") "no explicit damping in z-direction"
    end select
    if ((write_pe).and.(hyp_z.ne.0.0)) &
         & write(*,"(A,F6.3)") "hyp_z =     ",hyp_z


    !vparallel direction
    allocate(vdamping(-2:2,ll1:ll2))
    if (write_pe) then
       if (abs(hyp_v).gt.1e-5) then
          write(*,"(A,F6.3)") "hyp_v =     ", hyp_v
       else
          write(*,"(A)") "no explicit parallel velocity damping"
       end if
       write(*,*)
    end if
    do l=ll1,ll2
       vdamping(:,l)=hyp_v*dc4coeff
#ifndef oldvpbc
       if((l.eq.ll1).and.(my_pev.eq.0)) then
          vdamping(:,l)=0.
       elseif((l.eq.ll1+1).and.(my_pev.eq.0)) then
          vdamping(:,l)=hyp_v*(/0,1,-2,1,0/)* 0.0625
       elseif((l.eq.ll2-1).and.(my_pev.eq.n_procs_v-1)) then
          vdamping(:,l)=hyp_v*(/0,1,-2,1,0/)* 0.0625
       elseif((l.eq.ll2).and.(my_pev.eq.n_procs_v-1)) then
          vdamping(:,l)=0.
       endif
#endif
    enddo

    !GyroLES hypperdiffusion
    if (GyroLES) then
        allocate(hyp_x_spec(ln1:ln2))
        hyp_x_spec = 0.0001
        allocate(hyp_y_spec(ln1:ln2))
        hyp_y_spec = 0.0001
        if (.not.damp_i) then
            allocate(didiff(li1:li2))
            didiff = -(1/2.*deli*ki(:))**hyp_i_order
            Call mpi_allreduce (maxval(abs(didiff)),didiff_max,1,&
            MPI_REAL_TYPE, MPI_MAX,mpi_comm_x, ierr) 
        endif 
        if (.not.damp_j) then
            allocate(djdiff(lj1:lj2))
            djdiff = -(1/2.*delj*kj(lj1:lj2))**hyp_j_order
            Call mpi_allreduce (maxval(abs(djdiff)),djdiff_max,1,&
            MPI_REAL_TYPE, MPI_MAX,mpi_comm_y, ierr)        
        endif
    end if

  end subroutine initialize_numerical_damping

  subroutine finalize_numerical_damping
    if (hyp_i.gt.0.) deallocate(didiff)
    if (hyp_j.gt.0.) deallocate(djdiff)
    if (hyp_perp.gt.0) deallocate(dperpdiff)
    deallocate(hyp_z_sten)    
    deallocate(vdamping)
    if (GyroLES) then
        deallocate(hyp_x_spec,hyp_y_spec)
        if(.not.damp_i) deallocate(didiff)
        if(.not.damp_j) deallocate(djdiff)
    endif

    end subroutine finalize_numerical_damping

end module numerical_damping
