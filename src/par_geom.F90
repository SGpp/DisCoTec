#include "redef.h"
#include "intrinsic_sizes.h"
!Module providing parameters and derived data types
!for magnetic geometry related quantities
!used by geometry.F90 and geometry interfaces
module par_geom

  implicit none

  public :: geomtype, initialize_geomtype, finalize_geomtype
  public :: Bref, Lref, magn_geometry, magn_geo_geofile
  public :: geomdir, geomfile, flux_label, x_def, flux_pos
  public :: q0, major_R, minor_r, trpeps, shat, amhd, rhostar, Bprof_coeffs

  private

  ! magnetic geometry
  Character(len=80)  :: magn_geometry=""
  Character(len=100) :: geomdir="./"
  Character(len=FILENAME_MAX)  :: geomfile="", magn_geo_geofile=""
  Character(len=80)  :: flux_label= " "

  ! ... used in chease module
  Character(len=80)  :: x_def = "arho_t"
  Real:: flux_pos=0.0

  Real:: q0=-1111.0, major_R=1.,minor_r=1.
  Real:: trpeps=0.0, rhostar=0.0

  Real:: shat=0.0, amhd=0.0

  real, dimension(0:8) :: Bprof_coeffs


  TYPE geomtype
     character(len=10) :: name  ! geometry type
 
     !metric coefficients
     real, dimension(:,:,:), allocatable :: gii, gij, giz
     real, dimension(:,:,:), allocatable :: gjj, gjz, gzz
     complex, dimension(:,:,:), allocatable :: shift !defined complex for the exchange
     complex, dimension(:,:,:,:), allocatable:: phasefac !for shifted metric

     !Bfield and jacobian
     real, dimension(:,:,:), allocatable :: Bfield, jacobian
     real, dimension(:), allocatable :: avg_jaco_yz
     real :: avg_jaco

     !Bfield derivatives
     real, dimension(:,:,:), allocatable :: dBdi, dBdj, dBdz

     !curvature 
     real, dimension(:,:,:), allocatable :: K_i, K_j

     !magnetic eq. coordinates (only required for output to diagnostics)
     real, dimension(:,:),allocatable:: R, Z, PHI, dxdR, dxdZ
     real, dimension(:,:),allocatable:: R_hat,Z_hat
  END TYPE geomtype

  real:: Bref=0, Lref=0

contains
  
  !>Allocates the arrays in variable geom of type geomtype
  !!Note: the arrays bounds are input variables since some
  !!interfaces might return "full" arrays while GENE should
  !!make use of distributed arrays
  subroutine initialize_geomtype(geom,pi1gl,pi2gl,pg1gl,pg2gl,pi1,pi2,&
       lj1,lj2,pj1,pj2,lk1,lk2,nzb)
    implicit none

    TYPE(geomtype) :: geom
    integer, intent(in) :: pi1gl,pi2gl,pg1gl,pg2gl,pi1,pi2,&
       lj1,lj2,pj1,pj2,lk1,lk2,nzb

    allocate(geom%gii(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%gij(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%giz(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%gjz(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%gjj(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%gzz(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%shift(pi1gl:pi2gl,pj1:pj2,lk1-nzb:lk2+nzb))
    allocate(geom%Bfield(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%dBdi(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%dBdj(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%dBdz(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%jacobian(pi1gl:pi2gl,pj1:pj2,lk1:lk2))
    allocate(geom%avg_jaco_yz(pg1gl:pg2gl))
    Allocate(geom%K_i(pi1:pi2,pj1:pj2,lk1:lk2))
    Allocate(geom%K_j(pi1:pi2,pj1:pj2,lk1:lk2))
    allocate(geom%phasefac(pi1gl:pi2gl,lj1:lj2,lk1:lk2,-nzb:nzb))
    allocate(geom%R(pi1gl:pi2gl,lk1:lk2), geom%Z(pi1gl:pi2gl,lk1:lk2), &
         & geom%PHI(pi1gl:pi2gl,lk1:lk2), geom%dxdR(pi1gl:pi2gl,lk1:lk2), &
         & geom%dxdZ(pi1gl:pi2gl,lk1:lk2),&
         geom%R_hat(pi1gl:pi2gl,lk1:lk2), geom%Z_hat(pi1gl:pi2gl,lk1:lk2))

    geom%giz=0.; geom%gjz=0.; geom%gzz=0.; geom%dBdi=0.; geom%dBdj=0.; 
    geom%phasefac=1.; geom%shift=0.
    geom%R=0.; geom%z=0.; geom%phi=0.; geom%dxdR=0.; geom%dxdZ=0.  

  end subroutine initialize_geomtype

  subroutine finalize_geomtype(geom)
    TYPE(geomtype) :: geom

    deallocate(geom%gii,geom%gij,geom%giz,geom%gjj)
    deallocate(geom%gjz,geom%dBdz,geom%Bfield,geom%jacobian)
    deallocate(geom%avg_jaco_yz,geom%phasefac)
    deallocate(geom%K_i,geom%K_j)
    deallocate(geom%gzz,geom%shift)
    deallocate(geom%dBdi, geom%dBdj)
    deallocate(geom%R, geom%Z, geom%PHI, geom%dxdR, geom%dxdZ)
    deallocate(geom%R_hat, geom%Z_hat)

  end subroutine finalize_geomtype
    

end module par_geom
