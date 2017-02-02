#include "redef.h"
#include "intrinsic_sizes.h"
  
  !> Computes the variables for the hybrid model
module hybrid
  use communications
  use geometry
  use discretization
  use coordinates, only: mu,lv,dv,vp

  implicit none
  public:: mem_est_hybrid
  public:: vp_star,trap,Bfieldmax,initialize_hyb_vars,finalize_hyb_vars
  public:: vp_weight_trap

  private

  real, dimension(:), allocatable :: Bfieldmax
  real, dimension(:,:,:), allocatable :: vp_star
  real, dimension(:,:), allocatable :: trap
  real, dimension(:), allocatable :: temparr_loc
  real, dimension(:,:,:,:), allocatable :: vp_weight_trap

contains

  !>Give an estimate of the memory requirements of this module
  Real function mem_est_hybrid(mem_req_in)
    real:: mem_req_in
    real:: mem_loc=0

    ! vp_star
    mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*lk0*lm0
    ! Bfieldmax and trap
    mem_loc=mem_loc+2*SIZE_OF_REAL_MB*px0
    mem_loc=mem_loc+SIZE_OF_REAL_MB*px0*lk0
    ! vp_weight_trap
    mem_loc=mem_loc+SIZE_OF_REAL_MB*pi0*lk0*ll0*lm0

    mem_est_hybrid=mem_req_in+mem_loc

  end function mem_est_hybrid

  subroutine initialize_hyb_vars

    call calc_Bfieldmax
    call calc_trap_fraction
    call calc_vpstar
    call calc_vp_weight_trap

  end subroutine initialize_hyb_vars

  subroutine calc_Bfieldmax

    integer :: i,ierr
    real :: temparr_loc(pi1gl:pi2gl)

    allocate(Bfieldmax(pi1gl:pi2gl))

    temparr_loc=0.D0
    do i=pi1gl,pi2gl
       temparr_loc(i) = maxval(geom%Bfield(i,:,:))
    end do
    Call mpi_allreduce(temparr_loc, Bfieldmax, px0,&
         MPI_REAL_TYPE, MPI_MAX, mpi_comm_z, ierr)

  end subroutine calc_Bfieldmax

  subroutine calc_trap_fraction

    integer :: i,k
    allocate(trap(pi1gl:pi2gl,lk1:lk2))

    ! Fraction of trapped particles
    Do k=lk1,lk2
       Do i=pi1gl,pi2gl
          trap(i,k)=sqrt(1.0-geom%Bfield(i,pj1,k)/Bfieldmax(i))
       End Do
    End Do
  end subroutine calc_trap_fraction

  subroutine calc_vpstar

    integer :: i,k,m
    allocate(vp_star(pi1:pi2,lk1:lk2,lm1:lm2))
    ! Defining the trapped model v*=sqrt[(2/m)*(B_max-B_o)*mu]
    ! v*=sqrt[(2/m)]*alpha_t*sqrt(B_max)*sqrt(mu) here normalized.

    do m = lm1, lm2
       do k = lk1,lk2
          do i = pi1, pi2
             vp_star(i,k,m) =trap(i,k)*sqrt(Bfieldmax(i)*mu(m))
             vp_star(i,k,m) = min(lv,vp_star(i,k,m))
          end do
       end do
    end do
  end subroutine calc_vpstar

  subroutine calc_vp_weight_trap

    integer :: i, k, m
    integer :: ll, lu, diff_pts
    real :: zp, w0, w1, w2, w3, w4
    real :: weight(0:nv0-1)

    allocate(vp_weight_trap(pi1:pi2,lk1:lk2,ll1:ll2,lm1:lm2))


    do m = lm1, lm2
       do k = lk1, lk2                                                       
          do i = pi1, pi2         
             ! upper index of the interval where +v* is located 
             lu = ceiling((vp_star(i,k,m)+lv)/dv)
             lu = min(lu,nv0-1)
             ! lower index of the interval where -v* is located 
             ll = nv0-1-lu
             ! distance to be used in the weights
             zp = (vp_star(i,k,m)-vp(lu-1))/dv
             ! number of grid points to be used for the integration 
             diff_pts = lu-ll+1

             ! Assigning weights for v_par integration within trapping region
             ! initialize 	     
             weight = 0.0
             if(diff_pts.ge.9) then ! Alternative extended Simpson's rule with cut cell 

                w0=17.0/48.0+(1.0/48.0)*(zp-1.0)*(zp**3+9.0*zp**2+19.0*zp+19.0)
                w1=59.0/48.0-(1.0/48.0)*(zp-1.0)**2*(3.0*zp**2+22.0*zp+35.0)         
                w2=43.0/48.0+(1.0/48.0)*(zp-1.0)**2*(3.0*zp**2+14.0*zp+7.0)          
                w3=49.0/48.0-(1.0/48.0)*(zp**2-1.0)**2                                       
                w4=1.0

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1
                weight(ll+2)=w2 ; weight(lu-2)=w2
                weight(ll+3)=w3 ; weight(lu-3)=w3

                weight(ll+4:lu-4) = w4

             elseif(diff_pts.eq.8) then

                w0=1.0/3.0+1.0/12.0*(zp-1.0)*(2.0*zp**2+5.0*zp+5.0)    
                w1=4.0/3.0-1.0/3.0 *(zp-1.0)**2*(zp+2.0)                
                w2=17.0/24.0+1.0/6.0 *(zp-1.0)**2*(zp+1.0/2.0)          
                w3=9.0/8.0                                             
                
                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1
                weight(ll+2)=w2 ; weight(lu-2)=w2
                weight(ll+3)=w3 ; weight(lu-3)=w3

             elseif(diff_pts.eq.7) then

                w0=17.0/48.0+(1.0/48.0)*(zp-1.0)*(zp**3+9.0*zp**2+19.0*zp+19.0)
                w1=59.0/48.0-(1.0/48.0)*(zp-1.0)**2*(3.0*zp**2+22.0*zp+35.0)
                w2=43.0/48.0+(1.0/48.0)*(zp-1.0)**2*(3.0*zp**2+14.0*zp+7.0)
                w3=25.0/24.0-(1.0/24.0)*(zp**2-1.0)**2

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1
                weight(ll+2)=w2 ; weight(lu-2)=w2
                weight(ll+3)=w3 

             elseif(diff_pts.eq.6) then

                w0=17.0/48.0+(1.0/48.0)*(zp-1.0)*(zp**3+9.0*zp**2+19.0*zp+19.0)
                w1=59.0/48.0-(1.0/48.0)*(zp-1.0)**2*(3.0*zp**2+22.0*zp+35.0)
                w2=11.0/12.0+(1.0/24.0)*(zp-1.0)**2*(zp**2+6.0*zp+3.0)

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1
                weight(ll+2)=w2 ; weight(lu-2)=w2

             elseif(diff_pts.eq.5) then

                w0=1.0/3.0+(1.0/12.0)*(zp-1.0)*(2.0*zp**2+5.0*zp+5.0)
                w1=4.0/3.0-(1.0/3.0)*(zp-1.0)**2*(zp+2.0)
                w2=2.0/3.0+(1.0/3.0)*(zp-1.0)**2*(zp+1.0/2.0) 

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1
                weight(ll+2)=w2 ;

             elseif(diff_pts.eq.4) then

                w0=3.0/8.0+1.0/12.0*(zp-1.0)*(2.0*zp**2+5.0*zp+5.0)
                w1=9.0/8.0-1.0/12.0*(zp-1.0)**2*(2.0*zp+7.0)

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1

             elseif(diff_pts.eq.3) then

                w0=zp**3/3.0
                w1=2.0*zp*(1.0-zp**2/3.0)

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ;

             elseif(diff_pts.eq.2) then

                zp=zp-1.0
                ll=ll-1
                lu=lu+1

                w0=3.0/8.0+1.0/12.0*(zp-1.0)*(2.0*zp**2+5.0*zp+5.0)
                w1=9.0/8.0-1.0/12.0*(zp-1.0)**2*(2.0*zp+7.0)

                weight(ll)=w0 ; weight(lu)=w0
                weight(ll+1)=w1 ; weight(lu-1)=w1

             elseif(diff_pts.eq.1) then
                ! do nothing
             end if
             ! only keeping weights relative to mesh points of domain
             vp_weight_trap(i, k, ll1:ll2,m) = weight(ll1:ll2)
          end do
       end do
    end do

    vp_weight_trap = vp_weight_trap*dv

  end subroutine calc_vp_weight_trap

  subroutine finalize_hyb_vars
    deallocate(Bfieldmax, trap, vp_star, vp_weight_trap)
  end subroutine finalize_hyb_vars

end module hybrid
