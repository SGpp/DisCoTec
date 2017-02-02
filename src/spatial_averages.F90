#include "redef.h"

module spatial_averages
  use discretization
  use geometry
  use communications
  use par_other

  implicit none

  public:: flux_surface_average, volume_average_abssq, volume_average
  public:: sum_int_3d, sum_int_z
  private

  interface flux_surface_average
     module procedure flux_surface_average_c, flux_surface_average_r, flux_surface_average_c_spec
  end interface

  interface volume_average_abssq
     module procedure volume_average_abssq_c
  end interface  
  
  interface volume_average
     module procedure volume_average_r
  end interface  

contains

  !> computes flux surface average of a 4D field
  !> (extension to species dependent input/output fields)
  !! \param dat_4D is complex input field without ghost cells in z
  !! \param avg is real output field (flux surface average, species dependent)
  subroutine flux_surface_average_c_spec(dat_4D, avg)
    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2,ln1:ln2), intent(in):: dat_4D
    real, dimension(lg1:lg2,ln1:ln2), intent(out):: avg
    !local variables
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2) :: dat_3D
    complex, dimension(lg1:lg2) :: tmp
    integer::n
    do n=ln1,ln2
       dat_3D(:,:,:) = dat_4D(:,:,:,n)
       call flux_surface_average_c(dat_3D, .false., tmp)
       avg(:,n)=real(tmp)
    enddo
  end subroutine flux_surface_average_c_spec

  !> computes flux surface average of a 3D field (complex)
  !! \param dat_3D is the field to average
  !! \param withbounds indicates presence of ghost cells in 3rd index
  !! \param avg is the result (flux surface average of dat_3D)
  subroutine flux_surface_average_c(dat_3D, withbounds, avg)
    complex, dimension(:,:,:), intent(in):: dat_3D
    logical, intent(in):: withbounds
    complex, dimension(lg1:lg2), intent(out):: avg
    integer:: j,k,ioff,joff,koff
    
    ! the passed array dat_3D is running locally from 1 to the end
    ! depending on if it has boundary points in the third index,
    ! the inner point indices are calculated with different offsets 
    ! the first two dimensions are always without boundaries. they are
    ! running over 1:li0 1:lj0
    ! 
    if (withbounds) then
       koff = 1-lk1+nzb
    else
       koff = 1-lk1
    endif
    ioff = 1-li1
    joff = 1-lj1

    avg=0.
    if(yx_order) then
       if (y_local) then
          if (p_has_0_mode) then
             do k=lk1,lk2
                avg = avg + dat_3D(li1+ioff,:,k+koff)*geom%jacobian(pi1,pj1,k)
             end do
             avg = avg/(geom%avg_jaco*nz0)
          end if
       else
          do k=lk1,lk2
             do j=lj1,lj2
                avg(j) = avg(j) + sum(dat_3D(:,j+joff,k+koff)*geom%jacobian(li1:li2,pj1,k))
             end do
          end do
          avg = avg/(geom%avg_jaco*nz0*nky0)
          !sum over y
          call sum_to_all_complex(avg, lg0, mpi_comm_x)
       end if
    else
       if (p_has_0_mode) then
          if (x_local) then
             do k=lk1,lk2
                avg = avg + dat_3D(:,lj1+joff,k+koff)*geom%jacobian(pi1,pj1,k)
             end do
             avg = avg/(geom%avg_jaco*nz0)
          else
             do k=lk1,lk2
                avg = avg + dat_3D(:,lj1+joff,k+koff)*geom%jacobian(li1:li2,pj1,k)
             end do
             avg = avg/(geom%avg_jaco_yz(li1:li2)*nz0)
          end if
       end if
    end if
    call sum_to_all_complex(avg, lg0, mpi_comm_z)
 
  end subroutine flux_surface_average_c

  !> computes flux surface average of a 3D field (real)
  !! \param dat_3D is the field to average
  !! \param withbounds indicates presence of ghost cells
  !! \param avg is the radial dependent result (flux surface average of dat_3D)
  subroutine flux_surface_average_r(dat_3D, withbounds, avg)
    real, dimension(:,:,:), intent(in):: dat_3D
    logical,intent(in)::withbounds
    real, dimension(lg1:lg2), intent(out):: avg
    integer:: j,k,ioff,joff,koff
   
    ! the passed array dat_3D is running locally from 1 to the end
    ! dependent on if it has boundary points in the third index,the calculation of
    ! the inner point indices is different
    ! the first two dimensions are always without boundaries, that means
    ! running over 1:li0 1:lj0
    ! 
    if (withbounds) then
       koff = 1-lk1+nzb
    else
       koff = 1-lk1
    endif
    ioff = 1-li1
    joff = 1-lj1

    avg=0.
    if(yx_order) then
       if (y_local) then
          if (p_has_0_mode) then
             do k=lk1,lk2
                !li1+ioff has the ky=0 mode that is selected
                avg = avg + dat_3D(li1+ioff,:,k+koff)*geom%jacobian(pi1,pj1,k)
             end do
             avg = avg/(geom%avg_jaco*nz0) 
          end if
       else
          do k=lk1,lk2
             do j=lj1,lj2
                avg(j) = avg(j) + sum(dat_3D(:,j+joff,k+koff)*geom%jacobian(li1:li2,pj1,k))
             end do
          end do
          avg = avg/(geom%avg_jaco*nz0*nky0)
          !sum over y
          call my_sum_to_all_real(avg, lg0, mpi_comm_x)
       end if
    else
       !lj1+joff has the ky=0 mode that is selected
       if (p_has_0_mode) then
          if(x_local) then
             do k=lk1,lk2
                avg = avg + dat_3D(:,lj1+joff,k+koff)*geom%jacobian(pi1,pj1,k)
             end do
             avg = avg/(geom%avg_jaco*nz0)
          else
             do k=lk1,lk2
                avg = avg + dat_3D(:,lj1+joff,k+koff)*geom%jacobian(li1:li2,pj1,k)
             end do
             avg = avg/(geom%avg_jaco_yz(li1:li2)*nz0)
          end if
       end if
    end if
    call my_sum_to_all_real(avg, lg0, mpi_comm_z)
    
  end subroutine flux_surface_average_r

  !> computes volume average of the absolute square value of 
  !! a 3D field (complex)
  !! \param dat_3D is the field to average
  !! \param withbounds indicates presence of ghost cells in 3rd index
  !! \param avg is the result (volume average of dat_3D)
  subroutine volume_average_abssq_c(dat_3D, withbounds, avg)
    complex, dimension(:,:,:), intent(in):: dat_3D
    logical, intent(in):: withbounds
    real, intent(out):: avg
    integer:: j,k,ioff,joff,koff
    
    ! the passed array dat_3D is running locally from 1 to the end
    ! depending on if it has boundary points in the third index,
    ! the inner point indices are calculated with differen offsets 
    ! the first two dimensions are always without boundaries. they are
    ! running over 1:li0 1:lj0
    ! 
    if (withbounds) then
       koff = 1-lk1+nzb
    else
       koff = 1-lk1
    endif
    ioff = 1-li1
    joff = 1-lj1

    avg=0.
    if(yx_order) then
       if (y_local) then
          do k = lk1,lk2
             avg = avg + 2.0*Sum(abs(dat_3D(:,:,k+koff))**2)*geom%jacobian(pi1,pj1,k)
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(abs(dat_3D(li1+ioff,:,k+koff))**2)*&
                  &geom%jacobian(pi1,pj1,k)
          enddo
          avg = avg / (Real(nz0)*geom%avg_jaco)
       else
          do k = lk1,lk2
             do j=lj1,lj2
                avg = avg + 2.0*Sum(abs(dat_3D(:,j+joff,k+koff))**2*geom%jacobian(li1:li2,pj1,k))
             enddo
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(abs(dat_3D(:,lj1+joff,k+koff))**2*&
                  &geom%jacobian(li1:li2,pj1,k))
          enddo
          avg = avg / (Real(nz0*nky0)*geom%avg_jaco)
       end if
    else
       if (x_local) then
          do k = lk1,lk2
             avg = avg + 2.0*Sum(abs(dat_3D(:,:,k+koff))**2)*geom%jacobian(pi1,pj1,k)
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(abs(dat_3D(:,lj1+joff,k+koff))**2)*&
                  &geom%jacobian(pi1,pj1,k)
          enddo
          avg = avg / (Real(nz0)*geom%avg_jaco)
       else
          do k = lk1,lk2
             do j=lj1,lj2
                avg = avg + 2.0*Sum(abs(dat_3D(:,j+joff,k+koff))**2*geom%jacobian(li1:li2,pj1,k))
             enddo
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(abs(dat_3D(:,lj1+joff,k+koff))**2*&
                  &geom%jacobian(li1:li2,pj1,k))
          enddo
          avg = avg / (Real(nz0*nx0)*geom%avg_jaco)
       end if
    end if

    call my_sum_to_all_real_0d(avg, mpi_comm_xyz)
 
  end subroutine volume_average_abssq_c

!> computes volume average of 
  !! a 3D field (real)
  !! \param dat_3D is the field to average
  !! \param withbounds indicates presence of ghost cells in 3rd index
  !! \param avg is the result (volume average of dat_3D)
  subroutine volume_average_r(dat_3D, withbounds, avg)
    real, dimension(:,:,:), intent(in):: dat_3D
    logical, intent(in):: withbounds
    real, intent(out):: avg
    integer:: j,k,ioff,joff,koff
    
    ! the passed array dat_3D is running locally from 1 to the end
    ! depending on if it has boundary points in the third index,
    ! the inner point indices are calculated with different offsets 
    ! the first two dimensions are always without boundaries. they are
    ! running over 1:li0 1:lj0
    ! 
    if (withbounds) then
       koff = 1-lk1+nzb
    else
       koff = 1-lk1
    endif
    ioff = 1-li1
    joff = 1-lj1

    avg=0.
    if(yx_order) then
       if (y_local) then
          do k = lk1,lk2
             avg = avg + 2.0*Sum(dat_3D(:,:,k+koff)*geom%jacobian(pi1,pj1,k))
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(dat_3D(li1+ioff,:,k+koff))*&
                  &geom%jacobian(pi1,pj1,k)
          enddo
          avg = avg / (Real(nz0)*geom%avg_jaco)
       else
          do k = lk1,lk2
             do j=lj1,lj2
                avg = avg + 2.0*Sum(dat_3D(:,j+joff,k+koff)*geom%jacobian(li1:li2,pj1,k))
             enddo
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(dat_3D(:,lj1+joff,k+koff)*&
                  &geom%jacobian(li1:li2,pj1,k))
          enddo
          avg = avg / (Real(nz0*nky0)*geom%avg_jaco)
       end if
    else
       if (x_local) then
          do k = lk1,lk2
             avg = avg + 2.0*Sum(dat_3D(:,:,k+koff))*geom%jacobian(pi1,pj1,k)
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(dat_3D(:,lj1+joff,k+koff))*&
                  &geom%jacobian(pi1,pj1,k)
          enddo
          avg = avg / (Real(nz0)*geom%avg_jaco)
       else
          do k = lk1,lk2
             do j=lj1,lj2
                avg = avg + 2.0*Sum(dat_3D(:,j+joff,k+koff)*geom%jacobian(li1:li2,pj1,k))
             enddo
             !0 mode is not mirrored to negative k
             if (p_has_0_mode) avg = avg - Sum(dat_3D(:,lj1+joff,k+koff)*&
                  &geom%jacobian(li1:li2,pj1,k))
          enddo
          avg = avg / (Real(nz0*nx0)*geom%avg_jaco)
       end if
    end if

    call my_sum_to_all_real_0d(avg, mpi_comm_xyz)
 
  end subroutine volume_average_r

  !!*****************************************************************************************

  !>Sums / integrates over kx (x),ky,z
  subroutine sum_int_3d(v3d,v1d)
    !This subroutine performs the z integral and sum over kx,ky

    implicit none
    complex, dimension(li1:li2,lj1:lj2,lk1:lk2), intent(in) :: v3d
    real, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d_tmp

    real, intent(out) :: v1d
#if 0
    real :: v1d_loc
    integer :: ierr
#endif
      
    v3d_tmp = real(v3d)
    if (yx_order) then
       if ((xy_local).and.(evenx.eq.1)) v3d_tmp(:,hkx+1,:) = 0.0
    else
       if ((xy_local).and.(evenx.eq.1)) v3d_tmp(hkx+1,:,:) = 0.0
    endif
#if 0
    if (p_has_0_mode) v3d_tmp(:,lj1,:) = 0.5*v3d_tmp(:,lj1,:)

    if (xy_local) then
       v1d_loc = 2.0*sum(sum(sum(v3d_tmp,1),1)*geom%jacobian(pi1,pj1,:))/&
            (real(nz0)*geom%avg_jaco)
    else
       if (y_local) then
          v1d_loc = 2.0*sum(sum(v3d_tmp,2)*geom%jacobian(pi1:pi2,pj1,:))/&
               (real(nx0*nz0)*geom%avg_jaco)
       else
          v1d_loc = 2.0*sum(sum(v3d_tmp,2)*geom%jacobian(pi1:pi2,pj1,:))/&
               (real(nky0*nz0)*geom%avg_jaco)
       endif
    endif
    call mpi_allreduce(v1d_loc,v1d,1,MPI_REAL_TYPE,MPI_SUM, mpi_comm_xyz, ierr)
#else
    !there is now an interface to the spatial_average.F90 module, 
    !that yields the same result in diag_energy
    call volume_average(v3d_tmp,.false.,v1d)
#endif

  end subroutine sum_int_3d




  !*****************************************************************************************
  Subroutine sum_int_z(v3d,v2d)
    !This subroutine performpe the z integration

    implicit none
    real, dimension(li1:li2,lj1:lj2,lk1:lk2) :: v3d
    Real, dimension(li1:li2,lj1:lj2),intent(out) :: v2d
    Real, dimension(li1:li2,lj1:lj2) :: v2d_loc
    integer :: k,ierr

    v2d_loc=0.0
    v2d = 0.0
    do k =lk1,lk2
       v2d_loc(:,:) = v2d_loc(:,:) + v3d(:,:,k)*geom%jacobian(pi1,pj1,k)
    enddo

    v2d_loc=v2d_loc/(real(nz0)*geom%avg_jaco)

    call mpi_allreduce(v2d_loc,v2d,size(v2d_loc),MPI_REAL_TYPE,MPI_SUM,mpi_comm_z, ierr)

  end subroutine sum_int_z


end module spatial_averages
