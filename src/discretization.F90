#include "redef.h"
!>Contains the GENE parameters related to the discretization (index ranges, number of points in the
!!various directions, points per process)
module discretization
  USE par_in, ONLY: turbdeal, y_local, x_local, &
       &parscheme, nonlinear, collision_op, arakawa, nblocks, Omega0_tor
  use par_other, only: xy_local
  use, intrinsic :: iso_c_binding

  implicit none

 ! Parameters defining the grid
  Integer:: nx0         ! Number of points in x direction
  Integer:: nky0        ! Number of positive modes in y direction
  Integer:: ni0
  Integer(C_INT),bind(c):: nj0
  Integer:: nz0         ! Number of points in z direction
  Integer:: nv0         ! Number of v-parallel points
  Integer:: nw0         ! Number of points for mu-direction
  Integer:: n_spec         ! Number of species

  ! Number of boundary cells for finite differences in x, z, v, w direction
  Integer:: nxb
  Integer:: nyb
  Integer:: nib
  Integer:: nzb
  Integer:: nvb
  Integer:: nwb

  ! Parameters for parallelization, default: no parallelization
  Integer:: n_procs_s = -1     ! Number of processes for species
  Integer:: n_procs_w = -1     ! Number of processes for mu direction
  Integer:: n_procs_v = -1     ! Number of processes for v-parallel direction
  Integer:: n_procs_y = -1     ! Number of processes for y direction
  Integer:: n_procs_z = -1     ! Number of processes for z direction
  Integer:: n_procs_x = -1     ! Number of processes for x direction
  Integer:: n_procs_sim = -1   ! Total number of MPI processes for one gene run
  Integer:: n_parallel_sims=-1 ! Number of parallel gene simulations
  ! min and max number of processes for each direction, which is tested int
  ! the performance optimization
  integer :: min_npx=0, max_npx = 10000
  integer :: min_npy=0, max_npy = 10000
  integer :: min_npz=0, max_npz = 10000
  integer :: min_npv=0, max_npv = 10000
  integer :: min_npw=0, max_npw = 10000
  integer :: min_nps=0, max_nps = 10000

  !processor coordinates
  INTEGER:: mype=-1, my_pespec, my_pew, my_pez, my_pey, my_pev, my_pex, my_sim, mype_gl

  !number of points/process without boundaries
  Integer(C_INT), bind(c) :: lg0, lh0, li0, lj0, lk0, ll0, lm0, ln0

  ! total number of profile points
  integer:: px0, pmx0, pmg0

  !number of profile points/process
  Integer:: pi0, pmi0, pj0, pn0

  !number of points/process including boundaries
  Integer:: lx0, ly0, lz0, lv0, lw0

  ! Lower and upper bounds of arrays and inner Points:
  Integer::&
       & lbi, li1, li2, ubi,&
       &      pi1, pi2,     &
       &    pi1gl, pi2gl,   &
       &    pg1gl, pg2gl,   &
       &     pmi1, pmi2,    &
       &   pmi1gl, pmi2gl,  &
       &      pj1, pj2,     &
       &      lj1, lj2,     &
       & lbz, lk1, lk2, ubz,&
       & lbv, ll1, ll2, ubv,&
       & lbw, lm1, lm2, ubw,&
       &      ln1, ln2,     &
       &      pn1, pn2,     &
       &      lg1, lg2,     &
       &      lh1, lh2

  !dealiasing bounds
  integer(C_INT),bind(c):: li1da, li2da, li0da, ly0da , lbida, ubida, lbida2, ubida2
  integer:: nibda

  Integer:: lkx, hkx, hky, hki
  integer:: lij0, lijk0, llm0, lkl0, lklm0, lklmn0, lijz0
  integer(KIND=8) :: ljklmn0, lijklmn0, lzvwn0, lzvw0, lijkl0

  integer:: evenx, even_ni0
  Integer:: vlen, vlen_gl
  integer:: kx_offset

  !mode shift for linear runs
  integer:: ky0_ind=1, kj0_ind

  !auxiliary variable (set to input nky0)
  integer:: nky0_in

  !order of the x and y coordinate
  logical:: yx_order
contains

  subroutine set_discretization_defaults
    ky0_ind = 1
  end subroutine set_discretization_defaults

  !>Check parallelization relevant parameters
  !!\param print_stop if true check_discretization will stop and print an error message
  !!       whenever an invalid parallelization is detected; otherwise a logical value will be returned
  !!\return {status of parallelization}
  function check_parallelization(print_stop) result(valid_parall)
    Logical, intent(in) :: print_stop
    Logical :: valid_parall

    if (mod(n_spec,n_procs_s).ne.0) then
       valid_parall=.false.
       if (print_stop) Stop &
            "species-parallelization failed: n_spec has to be a multiple of n_procs_s!"
    elseif ((mod(nv0,n_procs_v).ne.0).or.(nv0/n_procs_v.lt.2)) then
       valid_parall=.false.
       if (print_stop) Stop &
            "v-parallelization failed: nv0 has to be a multiple of n_procs_v and nv0/n_procs_v has to be greater 1!"
    elseif (mod(nw0, n_procs_w).ne.0) then
       valid_parall=.false.
       if (print_stop) Stop&
         "w-parallelization failed: nw0 has to be a multiple of n_procs_w!"
    elseif ((n_procs_x .ne. 1) .and. x_local) then
       valid_parall = .false.
       If (print_stop) Stop &
            "parallelization failed: n_procs_x must be 1 in x_local mode!"
    elseif ((n_procs_y .ne. 1) .and. (.not.y_local)) Then
       valid_parall = .false.
       If (print_stop) Stop &
            "parallelization failed: n_procs_y must be 1 for y_local=.false.!"
    elseif (mod(nky0,n_procs_y).ne.0) then
       valid_parall=.false.
       if (print_stop) Stop &
         "y-parallelization failed: nky0 has to be a multiple of n_procs_y!"
    elseif (mod(nx0,n_procs_x).ne.0) then
       valid_parall=.false.
       if (print_stop) Stop &
         "x-parallelization failed: nx0 has to be a multiple of n_procs_x!"
    elseif (n_procs_x*n_procs_y > nx0) then
       valid_parall=.false.
       if (print_stop) Stop &
         "parallelization failed: n_procs_x*n_procs_y must not exceed nx0!"
    elseif (nxb>(nx0/n_procs_x)) then
       valid_parall = .false.
       If (print_stop) Stop &
         "x-parallelization failed: number of x grid points is smaller than required number of ghost cells"
    elseif ((mod(nz0, n_procs_z).ne.0).or.((nz0/n_procs_z.lt.2) .AND. (nz0 .GT. 1))) then
       valid_parall=.false.
       if (print_stop) Stop&
         "z-parallelization failed: nz0 has to be a multiple of n_procs_z and nz0/n_procs_z has to be greater 1!"
         !note: the second requirement need not be fulfilled for nz0 = 1
    elseif (nonlinear.and.(mod(nx0+1-evenx,2*n_procs_y).ne.0)) then
       valid_parall=.false.
       if (print_stop) Stop &
         "nonlinear y-parallelization failed: mod(nx0,2*n_procs_y) has to be 0 (nonlinear run)!"
    ELSEIF (arakawa.AND.(nx0.LT.2*nxb*n_procs_x*n_procs_y)) THEN
       valid_parall=.false.
       if (print_stop) stop &
            "Arakawa and x,y parallelization failed: nx0 has to be >= 2*nxb*n_procs_x*n_procs_y!"
    elseif (.not.xy_local.and.mod(nx0/n_procs_x,n_procs_y).ne.0) then
       valid_parall=.false.
       if (print_stop) stop &
            "x,y parallelization failed: nx0/n_procs_x has to be a multiple of n_procs_y!"
    else
       valid_parall=.true.
    endif
    
  end function check_parallelization

  !>Check parameters
  subroutine check_discretization
    if (nonlinear.and.(nx0.eq.1)) Stop &
         "nx0 must be greater than 1 for nonlinear run"
    if (nonlinear.and.(nky0.eq.1)) Stop &
         "nky0 must be greater than 1 for nonlinear run"

    if (ky0_ind.lt.0) then
       if (mype==0) write(*,"(a)") "only positive values of ky0_ind allowed, ky0_ind is set to zero"
       ky0_ind=0
    endif

    !save total number of ky (nky0 is set to 1 for time step computations)
    nky0_in=nky0

  end subroutine check_discretization


  !>Sets all index ranges etc., except for the indices that depend on nvb
  subroutine initialize_discretization(print_ini_msg)
    logical, intent(in) :: print_ini_msg
    logical :: write_pe

    write_pe = ((mype==0).and.(print_ini_msg))

    if (x_local) then
       nxb=0
    else
       nxb=2 ! centered4th is default
    endif

    if (y_local) then
       nyb=0
    else
       nyb=2
    endif

    if (parscheme.eq.'c2nd') then
       nzb = 1
    elseif (parscheme.eq.'c6th') then
       nzb = 3
    elseif (parscheme.eq.'c0th') then
       nzb = 0
    else
       nzb = 2
    endif

    if (n_procs_v.gt.1) then
       nvb=2
    else
       nvb=0
    endif    

    If ((collision_op.ne.'none').and.(n_procs_w.gt.1)) then
       nwb = 1
    else
       nwb = 0
    endif
    
    !number of physical gridpoints per process
    lg0 = nx0/n_procs_x
    lh0 = nky0/n_procs_y 
    lk0 = nz0/n_procs_z
    ll0 = (nv0)/n_procs_v 
    lm0 = nw0/n_procs_w
    ln0 = n_spec/n_procs_s

    If (lh0 < 1) Stop "lh0 must not be < 1"
    If (lk0 < nzb) Stop&
         "Number of local points must be more or equal to number of boundary points!!"

    !index limits

    ! x: radial direction
    lg1 = lg0*my_pex
    lg2 = lg1 + lg0 - 1
    hkx = (nx0-1)/2         !position of highest kx (without Nyquist)
    lkx = nx0-(nx0-1)/2     !position of lowest kx (without Nyquist)

    ! y: toroidal direction
    lh1 = my_pey*lh0
    lh2 = lh1 + lh0-1
    if (y_local) then
       hky=nky0-1
    else
       hky = nky0/2
    end if

    if (yx_order) then
       li0=lh0; li1=lh1; li2=lh2; ni0=nky0;
       lj0=lg0; lj1=lg1; lj2=lg2; nj0=nx0;
       nib=nyb
       kj0_ind=0
       hki = hky
    else
       li0=lg0; li1=lg1; li2=lg2; ni0=nx0;
       lj0=lh0; lj1=lh1; lj2=lh2; nj0=nky0;
       nib=nxb
       kj0_ind=ky0_ind
       hki = hkx
    end if
    even_ni0 = mod(ni0+1,2)
    lbi=li1-nib
    ubi=li2+nib

    ! profile information
    if (xy_local) then
       pi0  = 1
       pmi0 = 1
       px0  = 1
       pmx0 = 1
       pmg0 = 1 !number of radial magnetic profile points
    elseif (.not.x_local) then
       pi0  = lg0        !li0 if we really have a profile else 1
       pmi0 = lg0
       px0  = nx0
       pmx0 = nx0
       pmg0 = nx0
    elseif(.not.y_local) then
       pi0  = 1
       pmi0 = lh0
       px0  = 1
       pmx0 = nky0
       pmg0 = 1
       !introduce magnetic profiles
    endif
    if(pmi0.gt.1) then 
       pmi1 = li1   !local start and end of magnetic profile index
    else
       pmi1=pmi1gl  !start with pmi1gl (for globally initialized geometry arrays)
    endif
    pmi2 = pmi1 + pmi0 - 1
    pmi1gl = 0               !global start and end of magnetic profile index
    pmi2gl = pmx0 - 1
!    pi1 = li1                   !local start and end of x-profile index
!    pi2 = pi1 + pi0 - 1
!    pi1gl = 0               !global start and end of x-profile index
!    pi2gl = px0 - 1
    pi1 = pmi1
    pi2 = pmi2
    pg1gl = 0
    pg2gl = pmg0 - 1
    pi1gl = pmi1gl
    pi2gl = pmi2gl
    

    !x and y dependent geometry??
    pj0 = 1
    pj1 = lj1
    pj2 = pj1 + pj0 - 1


    ! z: parallel direction
    lk1 = lk0*my_pez
    lk2 = lk1 + lk0 - 1 
    lbz = lk1 - nzb
    ubz = lk2 + nzb

    ! v: v-parallel direction
    ll1 = ll0*my_pev
    ll2 = ll1 + ll0 - 1
    lbv = ll1 - nvb
    ubv = ll2 + nvb
 
    ! w: magnetic moment direction
    lm1 = lm0*my_pew
    lm2 = lm1 + lm0 - 1
    lbw = lm1 - nwb
    ubw = lm2 + nwb

    ! s: species direction
    ln1= ln0*my_pespec
    ln2= ln1 + ln0 - 1

    ! profile information
    If ((xy_local).and.(Omega0_tor.eq.0.0)) then
       pn0 = 1
    else
       pn0 = ln0
    endif
    pn1 = ln1
    pn2 = pn1 + pn0 - 1

    !number of physical gridpoints per process incl. boundaries
    lx0 = li0 + 2*nib
    ly0 = lj0
    lz0 = lk0 + 2*nzb
    lv0 = ll0 + 2*nvb
    lw0 = lm0 + 2*nwb

    !multidimensional grid sizes
    lij0 = li0*lj0
    lijk0 = li0*lj0*lk0
    lijkl0 = li0*lj0*lk0*ll0
    lijklmn0 = li0*lj0*lk0*ll0*lm0*ln0
    llm0 = ll0*lm0
    lkl0 = lk0*ll0
    lklm0 = lk0*llm0
    lklmn0 = lklm0*ln0
    ljklmn0 = lj0*lk0*ll0*lm0*ln0

    lzvw0 = lz0*lv0*lw0
    lzvwn0 = lz0*lv0*lw0*ln0
    lijz0 = li0*lj0*lz0

    !total length of state vector on each processor
    if ((evenx.eq.1).and.(xy_local)) then
       vlen=(lg0-evenx)*lh0*lk0*ll0*lm0*ln0
    else
       vlen=li0*lj0*lk0*ll0*lm0*ln0
    endif
    vlen_gl=vlen*n_procs_sim
    
    call initialize_da_bounds

  end subroutine initialize_discretization
  
  subroutine initialize_da_bounds
    ! boundaries for dealiasing
    If (xy_local) then
       if ((.not.nonlinear).or.turbdeal) then
          ! no additional kx modes for dealiasing
          li0da= li0
          li1da= li1
          li2da= li1da+li0da-1
          kx_offset= 0
          ly0da= 2*nky0
       else
          ! boundaries for dealiasing with 3/2-rule
          li0da= 3*(nx0+1-evenx)/2
          li1da= li1
          li2da= li1da+li0da-1
          kx_offset=nx0/2 +2*(1-evenx)  !number of modes added in kx direction for dealiasing
          ly0da= 3*nky0
       endif
    Else
       ! Always y (j) dealiasing
       if (turbdeal) then
          ly0da= 2*nj0
       else
          ly0da= 3*nj0
       endif
       ! no explicit dealiasing in x (i) direction
       li0da  = li0;  li1da  = li1
       li2da  = li1+li0da-1
       kx_offset = 0;
    Endif
  end subroutine initialize_da_bounds

end module discretization
