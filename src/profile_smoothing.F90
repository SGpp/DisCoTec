#include "redef.h"
!>Module containing routines for smoothing background 
!profiles for F0 adaptation.
Module profile_smoothing
  Use par_mod
  use geometry, only: minor_r, rhostar
  use lagrange_interpolation

  Implicit None

  PUBLIC :: smooth_profiles, smooth_reset_type, smooth_reset_width, smooth_reset
  PUBLIC :: set_profile_smoothing_defaults
 


  character(len=20):: smooth_reset_type = 'gauss'
  real:: smooth_reset_width
  logical :: smooth_reset = .false.

  PRIVATE
  
Contains

  subroutine set_profile_smoothing_defaults
    smooth_reset=.false.
    smooth_reset_type='gauss'
    smooth_reset_width= -1.0 ! Cannot be set to default 5*rhostar yet
  end subroutine set_profile_smoothing_defaults

  subroutine smooth_profiles(rawprof, smoothprof, width)
    real, intent(inout) :: width ! Width of region used for filtering 
    real,dimension(0:nx0-1,1:4), intent(in) :: rawprof
    real,dimension(0:nx0-1,1:4), intent(inout) :: smoothprof
    real,dimension(:),allocatable:: exp_rawprof !used for the boundary points
    real,dimension(:),allocatable:: srt_tmp, weight ! used for selection algorithm, weights for filters
    real,dimension(0:nx0-1):: temparr
    integer :: lim ! number of points for smoothing 2*lim+1
    integer :: minIndex
    real :: minValue, val_newcenter, val_oldcenter
    integer :: i,o,p,q
    
    lim=width*nx0
    lim=ceiling(real((lim-1))/2)
    allocate(exp_rawprof(-lim:nx0-1+lim),srt_tmp(-lim:lim),weight(-lim:lim))
    
    select case (smooth_reset_type)
    !smooth profiles by taking a running median of the gradients
    case('median')
       !smooth the gradient profiles
       do q=3,4
          exp_rawprof(0:nx0-1) = rawprof(:,q)
          do i=1,lim !mirror boundary points
             exp_rawprof(-i) = rawprof(i,q)
             exp_rawprof(nx0-1+i) = rawprof(nx0-1-i,q)
          enddo

          ! These loops replace the value of a profile point with the median of the
          ! 2*lim+1 points around it (including the point itself)
          do i=0,nx0-1
             srt_tmp = exp_rawprof(i-lim:i+lim)
             do o=-lim,0
                minIndex=o
                minValue=srt_tmp(o)
                do p=o+1,lim
                   if (srt_tmp(p) < minValue) then
                      minIndex = p
                      minValue = srt_tmp(p)
                   endif
                enddo
                srt_tmp(minIndex) = srt_tmp(o)
                srt_tmp(o) = minValue
             enddo
             smoothprof(i,q) = srt_tmp(0)
          enddo
       enddo

      ! integrate the logarithmic gradients (trapezoidal rule) to generate n(x) and T(x)
       do q=1,2
          smoothprof(:,q) = 0.0
          do i=0,nx0-2
             do o=0,i
                smoothprof(i,q)= smoothprof(i,q) + (smoothprof(o+1,q+2)+smoothprof(o,q+2))/2/nx0*lx*rhostar*minor_r
             enddo
          enddo
          !outward boundary point can be set to second last point (the boundary has
          !enough distortions already)
          smoothprof(nx0-1,q) = smoothprof(nx0-2,q)
          smoothprof(:,q) = exp(-smoothprof(:,q))
          call calc_center_value(rawprof(:,q),val_oldcenter)
          call calc_center_value(smoothprof(:,q),val_newcenter)
          smoothprof(:,q) = smoothprof(:,q)*val_oldcenter/val_newcenter
       enddo
    !smooth profiles by taking a running average of the profiles
    case('average')
       !boundary conditions: zero gradient
       do q=1,2
          exp_rawprof(0:nx0-1)=rawprof(:,q)
          do i=1,lim
             exp_rawprof(-i)=rawprof(i,q)
             exp_rawprof(nx0-1+i)=rawprof(nx0-1-i,q)
          enddo
          do i=0,nx0-1
             smoothprof(i,q)=sum(exp_rawprof(i-lim:i+lim))/(2*lim+1)
          enddo
          temparr = -log(smoothprof(:,q))
          call lag3deriv(temparr, xval, nx0,smoothprof(:,q+2), xval, nx0)
          smoothprof(:,q+2) = smoothprof(:,q+2)/(rhostar*minor_r)
       enddo
    !smooth profiles with a gaussian filter
    case('gauss')
       do q=1,2
          exp_rawprof(0:nx0-1)=rawprof(:,q)
          do i=1,lim
             exp_rawprof(-i)=rawprof(i,q)
             exp_rawprof(nx0-1+i)=rawprof(nx0-1-i,q)
          enddo
          do o=-lim,lim
             weight(o) = exp(-(real(o)/real(lim))**2/0.3)
          enddo
          weight(-lim:lim) = weight(-lim:lim)/sum(weight(-lim:lim))
          !if(mype == 0) print*, weight(-lim:lim)
          smoothprof(:,q) = 0.0
          do i=0,nx0-1
             do o=-lim,lim
                smoothprof(i,q)= smoothprof(i,q) + (exp_rawprof(i+o))*weight(o)
             enddo
          enddo
          temparr = -log(smoothprof(:,q))
          call lag3deriv(temparr, xval, nx0,smoothprof(:,q+2), xval, nx0)
          smoothprof(:,q+2) = smoothprof(:,q+2)/(rhostar*minor_r)
       enddo
    case default
       stop 'smooth_reset_type specified incorrectly!'
    end select

    deallocate(exp_rawprof,srt_tmp,weight)

  end subroutine smooth_profiles

end Module profile_smoothing
