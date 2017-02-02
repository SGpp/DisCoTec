! 30-08-2012

module discretization_adptv_module

  use discretization, only: nv0, nx0, nw0, nky0, nz0, n_spec, nvb, nib, nwb, &
       & lbv, ubv, lbw, ubw, &
       & n_procs_v, n_procs_x, n_procs_w, my_pev, mype_gl
  !  use blockindex
  use coordinates, only: dv, deli, mu, lv
  use geometry, only: rhostar

  implicit none
  private

  integer :: n_vx_blks, n_wx_blks ! number of blocks for different block structured grids
  real :: prob_vx, prob_wx ! probabilities to find a 'particle' in the block-str. subgrid domains
  character(len=20) :: opt_mthd_vx, opt_mthd_wx ! optimization method names
  real, dimension(:), allocatable :: blk_mks_r_v, blk_mks_v, blk_mks_r_w, blk_mks_w

  ! boolean flags to trigger adaptive grid usage in the computational routines
  logical :: is_grid_adptv = .false.
  ! flag to track if grid is initialized or not
  logical :: is_block_grid_vx_initialized_flag = .false., &
       & is_block_grid_wx_initialized_flag = .false.

  ! -------------- mirrored adaptive versions of variables ---------------------

  ! number of parallel velocity points for block-structured grid
  integer, dimension(:), allocatable :: nv0_adp

  ! number of magnetic moment points for block-structured grid
  integer, dimension(:), allocatable :: nw0_adp

  ! number of x points for block-structured grid
  integer, dimension(:), allocatable :: nx0_vadp, nx0_wadp
  integer, dimension(:,:), allocatable :: nx0_vwadp

  ! number of points/process without boundaries
  integer, dimension(:), allocatable :: ll0_adp, lm0_adp, &
       & li0_vadp, li0_wadp
  integer, dimension(:,:), allocatable :: li0_vwadp

  ! number of points/process including boundaries
  integer, dimension(:), allocatable :: lv0_adp, lw0_adp, &
       & lx0_vadp, lx0_wadp
  integer, dimension(:,:), allocatable :: lx0_vwadp

  ! lower and upper bounds of arrays and inner points:
  integer, dimension(:), allocatable:: &
       & lbv_adp, ll1_adp, ll2_adp, ubv_adp, &
       & lbw_adp, lm1_adp, lm2_adp, ubw_adp, &
       & lbi_vadp, li1_vadp, li2_vadp, ubi_vadp, pi1_vadp, pi2_vadp, &
       & lbi_wadp, li1_wadp, li2_wadp, ubi_wadp, pi1_wadp, pi2_wadp
  integer, dimension(:,:), allocatable:: &
       lbi_vwadp, li1_vwadp, li2_vwadp, ubi_vwadp, pi1_vwadp, pi2_vwadp

  ! special case taking parallel velocity and magnetic moment
  ! boundaries into consideration
  integer, dimension(:), allocatable :: li1_b_vadp, li2_b_vadp
  integer, dimension(:), allocatable :: li1_b_wadp, li2_b_wadp
  integer, dimension(:,:), allocatable :: li1_b_vwadp, li2_b_vwadp

  ! some development parameters
  logical, parameter :: DOUTPUT=.true.

  ! second (2) version adaptivity variables
  ! each block has the same number of parallel velocity points
  logical :: is_numofvpoints_const = .false.
  integer :: nv0_adp2, nv0_par                      ! recomputed number of points in each block
  integer, dimension(:), allocatable :: rv_mks_adp2 ! block marks for "splitted" loops

  ! ----------------------------------------------------------------------------

  public :: nv0_adp, n_vx_blks, n_wx_blks, &
       & prob_vx, prob_wx, opt_mthd_vx, opt_mthd_wx, &
       & blk_mks_r_v, blk_mks_v, blk_mks_r_w, blk_mks_w, &
       & li0_vadp, li1_vadp, li2_vadp, pi1_vadp, pi2_vadp, &
       & li0_wadp, li1_wadp, li2_wadp, pi1_wadp, pi2_wadp, &
       & li0_vwadp, li1_vwadp, li2_vwadp, pi1_vwadp, pi2_vwadp, &
       & ll0_adp, ll1_adp, ll2_adp, &
       & lm0_adp, lm1_adp, lm2_adp, &
       & li1_b_vwadp, li2_b_vwadp
  public :: is_numofvpoints_const, nv0_adp2, rv_mks_adp2
  public :: initialize_discretization_adptv, &
       & set_discretization_adptv_defaults, &
       & extract_blk_mks, &
       & finalize_discretization_adptv, &
       & is_block_grid_initialized, &
       & set_number_of_v_points
  public :: is_grid_adptv

#ifdef GIUTESTS
  public :: initialize_vx_subgrid, initialize_wx_subgrid, estimate_number_of_points, &
       & initialize_vx_subgrid_2
#endif

contains

  logical function is_block_grid_initialized() result(res)
    res = is_block_grid_vx_initialized_flag .and. is_block_grid_wx_initialized_flag
  end function is_block_grid_initialized

  subroutine initialize_discretization_adptv
    ! should be called always after initialiaze_discretization
    integer :: l, m

    if (.not.is_grid_adptv) then
       ! omit initialization
       return
    end if

    ! construct x, mu subgrid
    call initialize_wx_subgrid

    if(is_numofvpoints_const) then
       ! construct x, v|| subgrid with constant v|| number of points in each block
       call initialize_vx_subgrid_2
       ! estimate number of points in the block structured grid
       call estimate_number_of_points_2
    else
       ! construct x, v|| subgrid with constant resolution in each block
       call initialize_vx_subgrid
       ! initialization of hybrid v, w, x ranges
       allocate(nx0_vwadp(0:nv0-1, 0:nw0-1), li0_vwadp(0:nv0-1, 0:nw0-1), &
            & lx0_vwadp(0:nv0-1, 0:nw0-1), lbi_vwadp(0:nv0-1, 0:nw0-1), &
            & li1_vwadp(0:nv0-1, 0:nw0-1), li2_vwadp(0:nv0-1, 0:nw0-1), &
            & ubi_vwadp(0:nv0-1, 0:nw0-1), pi1_vwadp(0:nv0-1, 0:nw0-1), &
            & pi2_vwadp(0:nv0-1, 0:nw0-1))
       allocate(li1_b_vwadp(-nvb:nv0-1+nvb, -nwb:nw0-1+nwb), &
            & li2_b_vwadp(-nvb:nv0-1+nvb, -nwb:nw0-1+nwb))

       do m = 0, nw0-1
          do l = 0, nv0-1
             nx0_vwadp(l,m) = min(nx0_vadp(l), nx0_wadp(m))
             ! temperal solution only valid without parallelization
             ! TODO: reconsider a version with paralelization
             lx0_vwadp(l,m) = min(lx0_vadp(l), lx0_wadp(m))
             li0_vwadp(l,m) = min(li0_vadp(l), li0_wadp(m))
             lbi_vwadp(l,m) = max(lbi_vadp(l), lbi_wadp(m))
             li1_vwadp(l,m) = max(li1_vadp(l), li1_wadp(m))
             li2_vwadp(l,m) = min(li2_vadp(l), li2_wadp(m))
             li1_b_vwadp(l,m) = max(li1_b_vadp(l), li1_b_wadp(m))
             li2_b_vwadp(l,m) = min(li2_b_vadp(l), li2_b_wadp(m))
             ubi_vwadp(l,m) = min(ubi_vadp(l), ubi_wadp(m))
          end do
       end do
       pi1_vwadp = li1_vwadp
       pi2_vwadp = li2_vwadp

       if (DOUTPUT .and. (mype_gl.eq.0)) call estimate_number_of_points
       
    end if
    
  end subroutine initialize_discretization_adptv

  subroutine estimate_number_of_points
    integer :: np_rg, np_vx, np_wx, np_vwx
    integer :: l, m
    ! full (regular) grid
    np_rg = nx0*nky0*nz0*nv0*nw0*n_spec
    ! x, v|| block-structured grid
    np_vx = 0
    do l = 0, nv0-1
       np_vx = np_vx + nx0_vadp(l)*nky0*nz0*nw0*n_spec
    end do
    ! x, mu block-structured grid
    np_wx = 0
    do m = 0, nw0-1
       np_wx = np_wx + nx0_wadp(m)*nky0*nz0*nv0*n_spec
    end do
    ! x, v||, mu block-structured grid
    np_vwx = 0
    do m = 0, nw0-1
       do l = 0, nv0-1
          np_vwx = np_vwx + nx0_vwadp(l,m)*nky0*nz0*n_spec
       end do
    end do
    ! output of the information
    print *, "~~~~ Constant v|| resolution in each block grid version ~~~~"
    print *, "number of blocks ->"
    print *, "in v|| direction:  ", n_vx_blks
    print *, "in mu direction:   ", n_wx_blks
    print *, "number of points ->"
    print *, "regular full grid:    ", np_rg
    print *, "x, v|| b.s. grid:     ", np_vx
    print *, "x, mu b.s. grid:      ", np_wx
    print *, "x, v||, mu b.s. grid: ", np_vwx
  end subroutine estimate_number_of_points

  subroutine estimate_number_of_points_2
    integer :: np_rg, np_vx, np_wx, np_vwx
    integer :: m
    ! full (regular) grid
    np_rg = nx0*nky0*nz0*nv0_par*nw0*n_spec
    ! x, v|| block-structured grid
    np_vx = nx0*nky0*nz0*nv0_adp2*nw0*n_spec
    ! x, mu block-structured grid
    np_wx = 0
    do m = 0, nw0-1
       np_wx = np_wx + nx0_wadp(m)*nky0*nz0*nv0_par*n_spec
    end do
    ! x, v||, mu block-structured grid
    np_vwx = 0
    do m = 0, nw0-1
       np_vwx = np_vwx + nx0_wadp(m)*nky0*nv0_adp2*nz0*n_spec
    end do
    ! output of the information
    print *, "~~~~ Constant number of v|| points in each block grid version 2.0 ~~~~"
    print *, "number of v|| points in initial regular grid:", nv0_par
    print *, "number of v|| points in block structured grid:", nv0_adp2
    print *, "number of blocks ->"
    print *, "in v|| direction:  ", n_vx_blks
    print *, "in mu direction:   ", n_wx_blks
    print *, "number of points ->"
    print *, "regular full grid:    ", np_rg
    print *, "x, v|| b.s. grid:     ", np_vx
    print *, "x, mu b.s. grid:      ", np_wx
    print *, "x, v||, mu b.s. grid: ", np_vwx
  end subroutine estimate_number_of_points_2

  subroutine initialize_vx_subgrid
    integer :: indx
    integer, dimension(n_vx_blks+1) :: iblk_mks_r_v, iblk_mks_v
    integer, dimension(0:nx0-1) :: ll_shift
    integer :: r1, r2, nv0mod2, dr12

    ! some temperory arrays
    integer, dimension(:), allocatable :: nx0_b_vadp

    if (is_block_grid_vx_initialized_flag) then
       print *, "block structured x, v|| subgrid is already initialized!"
       return
    end if

    if (is_numofvpoints_const) then
       print *, "wrong choice of vx grid initialization (constant resolution)!"
       return
    end if

    iblk_mks_r_v = nint((blk_mks_r_v-blk_mks_r_v(1))/(deli*rhostar)) + 1
    iblk_mks_r_v(1) = 0

    iblk_mks_v = nint(2*blk_mks_v/dv) + 1
    nv0mod2 = mod(nv0,2)
    iblk_mks_v = iblk_mks_v + abs(nv0mod2 - mod(iblk_mks_v,2))
    iblk_mks_v(n_vx_blks+1) = 0

    allocate(nv0_adp(0:nx0-1), ll0_adp(0:nx0-1), lv0_adp(0:nx0-1), &
         & lbv_adp(0:nx0-1), ll1_adp(0:nx0-1), ll2_adp(0:nx0-1), &
         & ubv_adp(0:nx0-1))
    allocate(nx0_vadp(0:nv0-1), li0_vadp(0:nv0-1), lx0_vadp(0:nv0-1), &
         & lbi_vadp(0:nv0-1), li1_vadp(0:nv0-1), li2_vadp(0:nv0-1), &
         & ubi_vadp(0:nv0-1), pi1_vadp(0:nv0-1), pi2_vadp(0:nv0-1))
    allocate( nx0_b_vadp(-nvb:nv0-1+nvb), &
         & li1_b_vadp(-nvb:nv0-1+nvb), li2_b_vadp(-nvb:nv0-1+nvb))

    r1 = 0
    r2 = nv0-1
    do indx = 1, n_vx_blks
       nv0_adp(iblk_mks_r_v(indx) : iblk_mks_r_v(indx+1) - 1) = &
            & iblk_mks_v(indx)
       ll_shift(iblk_mks_r_v(indx) : iblk_mks_r_v(indx+1) - 1) = &
            & (nv0 - iblk_mks_v(indx))/2
       dr12 = (iblk_mks_v(indx)-iblk_mks_v(indx+1))/2
       nx0_vadp(r1 : r1 + dr12) = &
            & iblk_mks_r_v(indx+1)
       nx0_vadp(r2 - dr12 : r2) = &
            & iblk_mks_r_v(indx+1)
       nx0_b_vadp(r1 - nvb : r1 + dr12 - nvb) = iblk_mks_r_v(indx+1)
       nx0_b_vadp(r2 - dr12 + nvb : r2 + nvb) = iblk_mks_r_v(indx+1)
       r1 = r1 + dr12
       r2 = r2 - dr12
    end do

    ! print *, "nv0_adp(0) : ", nv0_adp(0), ", nv0 : ", nv0
    if (nv0_adp(0) .ne. nv0) stop "wrong maximum number of grid points"

    ! print *, "second version of adaptive array nv0_adp: ", nv0_adp
    ! print *, "maximum number of points in nx0 direction: ", nx0
    ! print *, "adaptive array nx0_adp: ", nx0_adp

    !number of physical gridpoints per process
    ll0_adp = nv0_adp/n_procs_v
    li0_vadp = nx0_vadp/n_procs_x

    !index limits

    ! v: v-parallel direction
    ll1_adp = ll_shift + ll0_adp*my_pev
    ll2_adp = ll1_adp + ll0_adp - 1
    lbv_adp = ll1_adp - nvb
    ubv_adp = ll2_adp + nvb
    lv0_adp = ll0_adp + 2*nvb

    ! x: flux surface direction, no "twisted" lg version!!!
    ! \todo consider "twisted" x <-> y variant if the whole thing
    !  makes sense

    li1_vadp = li0_vadp*my_pev
    li2_vadp = li1_vadp + li0_vadp - 1
    li1_b_vadp = nx0_b_vadp*my_pev/n_procs_x
    li2_b_vadp = li1_b_vadp + nx0_b_vadp/n_procs_x - 1
    lbi_vadp = li1_vadp - nib
    ubi_vadp = li2_vadp + nib
    lx0_vadp = li0_vadp + 2*nib

    ! print *, "nv0: ", nv0, ", ll0: ", ll0, ", ll1: ", ll1, ", ll2: ", ll2
    ! print *, "iblock_marks_v: ", iblock_marks_v
    ! print *, "nv0_adp: ", nv0_adp
    ! print *, "ll1_adp: ", ll1_adp
    ! print *, "ll2_adp: ", ll2_adp
    ! print *, "ll_shift: ", ll_shift
    ! print *, " nx0: ", nx0, ", li0: ", li0, ", li1: ", li1, ", li2: ", li2
    ! print *, "iblock_marks_r: ", iblock_marks_r
    ! print *, "nx0_adp: ", nx0_adp
    ! print *, "li1_vadp: ", li1_vadp
    ! print *, "li2_vadp: ", li2_vadp
    !stop "testing"

    ! cheating with the profile indices
    ! \todo figure out how should pi* indices should be treated
    pi1_vadp = li1_vadp
    pi2_vadp = li2_vadp

    is_block_grid_vx_initialized_flag = .true.

  end subroutine initialize_vx_subgrid

  subroutine initialize_vx_subgrid_2

    if (is_block_grid_vx_initialized_flag) then
       print *, "block structured x, v|| subgrid is already initialized!"
       return
    end if

    if (.not.is_numofvpoints_const) then
       print *, "wrong choice of vx grid initialization (constant number of v|| points)!"
       return
    end if

    allocate(rv_mks_adp2(0:n_vx_blks-1))

    rv_mks_adp2 = nint((blk_mks_r_v(2:)-blk_mks_r_v(1))/(deli*rhostar)) + 1

    is_block_grid_vx_initialized_flag = .true.

  end subroutine initialize_vx_subgrid_2

  subroutine initialize_wx_subgrid
    integer :: indx, indxb
    integer, dimension(n_wx_blks+1) :: iblk_mks_r_w, iblk_mks_w

    ! some temperory arrays
    integer, dimension(:), allocatable :: nx0_b_wadp

    if (is_block_grid_wx_initialized_flag) then
       print *, "block structured x, mu subgrid is already initialized!"
       return
    end if

    iblk_mks_r_w = nint((blk_mks_r_w-blk_mks_r_w(1))/(deli*rhostar)) + 1
    iblk_mks_r_w(1) = 0

    do indxb = 1, n_wx_blks
       do indx = 0, nw0-1
          if ((mu(indx).ge.blk_mks_w(indxb)).or.(indx.eq.nw0-1)) then
             iblk_mks_w(indxb) = indx + 1
             exit
          end if
       end do
    end do
    iblk_mks_w(n_wx_blks+1) = 0

    allocate(nw0_adp(0:nx0-1), lm0_adp(0:nx0-1), lw0_adp(0:nx0-1), &
         & lbw_adp(0:nx0-1), lm1_adp(0:nx0-1), lm2_adp(0:nx0-1), &
         & ubw_adp(0:nx0-1))
    allocate(nx0_wadp(0:nw0-1), li0_wadp(0:nw0-1), lx0_wadp(0:nw0-1), &
         & lbi_wadp(0:nw0-1), li1_wadp(0:nw0-1), li2_wadp(0:nw0-1), &
         & ubi_wadp(0:nw0-1), pi1_wadp(0:nw0-1), pi2_wadp(0:nw0-1))
    allocate(nx0_b_wadp(-nwb:nw0-1+nwb), &
         & li1_b_wadp(-nwb:nw0-1+nwb), li2_b_wadp(-nwb:nw0-1+nwb))

    do indx = 1, n_wx_blks
       nw0_adp(iblk_mks_r_w(indx) : iblk_mks_r_w(indx+1) - 1) = &
            & iblk_mks_w(indx)
       nx0_wadp(iblk_mks_w(n_wx_blks-indx+2):iblk_mks_w(n_wx_blks-indx+1)-1) = &
            & iblk_mks_r_w(n_wx_blks-indx+2)
       nx0_b_wadp(iblk_mks_w(n_wx_blks-indx+2)+nwb:iblk_mks_w(n_wx_blks-indx+1)-1+nwb) = &
            & iblk_mks_r_w(n_wx_blks-indx+2)
    end do

    if (nw0_adp(0) .ne. nw0) stop "wrong maximum number of mu grid points"

    !number of physical gridpoints per process
    lm0_adp = nw0_adp/n_procs_w
    li0_wadp = nx0_wadp/n_procs_x

    !index limits
    ! magnetic moment direction
    lm1_adp = lm0_adp*my_pev
    lm2_adp = lm1_adp + lm0_adp - 1
    lbw_adp = lm1_adp - nwb
    ubw_adp = lm2_adp + nwb
    lw0_adp = lm0_adp + 2*nwb

    ! x: flux surface direction, no "twisted" lg version!!!
    ! \todo consider "twisted" x <-> y variant if the whole thing
    !  makes sense

    li1_wadp = li0_wadp*my_pev
    li2_wadp = li1_wadp + li0_wadp - 1
    li1_b_wadp = nx0_b_wadp*my_pev/n_procs_x
    li2_b_wadp = li1_b_wadp + nx0_b_wadp/n_procs_x - 1
    lbi_wadp = li1_wadp - nib
    ubi_wadp = li2_wadp + nib
    lx0_wadp = li0_wadp + 2*nib

    ! cheating with the profile indices
    ! \todo figure out how should pi* indices should be treated
    pi1_wadp = li1_wadp
    pi2_wadp = li2_wadp

    is_block_grid_wx_initialized_flag = .true.

  end subroutine initialize_wx_subgrid

  subroutine set_discretization_adptv_defaults
    n_vx_blks = 1
    n_wx_blks = 1
  end subroutine set_discretization_adptv_defaults

  subroutine extract_blk_mks( &
       & blk_mks_r_v_big, blk_mks_v_big, &
       & blk_mks_r_w_big, blk_mks_w_big )

    real, dimension(:), intent(in) :: blk_mks_r_v_big, blk_mks_v_big
    real, dimension(:), intent(in) :: blk_mks_r_w_big, blk_mks_w_big

    ! allocate blk_mks arrays for r, v|| subgrid
    allocate(blk_mks_r_v(n_vx_blks+1))
    allocate(blk_mks_v(n_vx_blks+1))

    ! allocate blk_mks arrays for r, mu subgrid
    allocate(blk_mks_r_w(n_wx_blks+1))
    allocate(blk_mks_w(n_wx_blks+1))

    ! extract blk_mks array values
    blk_mks_r_v = blk_mks_r_v_big(1:n_vx_blks+1)
    blk_mks_v = blk_mks_v_big(1:n_vx_blks+1)

    blk_mks_r_w = blk_mks_r_w_big(1:n_wx_blks+1)
    blk_mks_w = blk_mks_w_big(1:n_wx_blks+1)
  end subroutine extract_blk_mks

  subroutine set_number_of_v_points
    real :: dv
    dv = 2.0*lv/(nv0-1.0)
    if(is_numofvpoints_const) then
       nv0_adp2 = 2*blk_mks_v(n_vx_blks)/dv + 1
       ! nv0_adp2 should be even
       nv0_adp2 = nv0_adp2 + mod(nv0_adp2,2)
       ! nv0_adp2 should be divisible by n_procs_v
       if(mod(nv0_adp2,n_procs_v).ne. 0) then
          print *, "nv0_adp2 should be divisible by n_procs_v"
          nv0_adp2 = nv0_adp2 + mod(nv0_adp2,n_procs_v)
       end if

       ! TODO: discuss about better solution
       ! workaround to safe number of points in data structures
       nv0_par = nv0
       nv0 = nv0_adp2
    else
       ! to be on the safe side
       nv0_par  = nv0
       nv0_adp2 = nv0
    end if
  end subroutine set_number_of_v_points
     
  subroutine finalize_discretization_adptv

    ! deallocation of v, w, x block-structured grid arrays
    if(is_block_grid_vx_initialized_flag.and. &
         & is_block_grid_wx_initialized_flag) then
       deallocate(nx0_vwadp, li0_vwadp, &
            & lx0_vwadp, lbi_vwadp, &
            & li1_vwadp, li2_vwadp, &
            & li1_b_vwadp, li2_b_vwadp, &
            & ubi_vwadp, pi1_vwadp, &
            & pi2_vwadp)
    end if
    if(is_block_grid_vx_initialized_flag) then
       ! deallocation of x, v|| block-structured grid arrays
       if(is_numofvpoints_const) then
          deallocate(rv_mks_adp2)
          is_numofvpoints_const = .false.
       else
          deallocate(nv0_adp, nx0_vadp)
          deallocate(ll0_adp, li0_vadp)
          deallocate(lv0_adp, lx0_vadp)
          deallocate(lbv_adp, ll1_adp, ll2_adp, ubv_adp)
          deallocate(lbi_vadp, li1_vadp, li2_vadp, ubi_vadp)
          deallocate(pi1_vadp, pi2_vadp)
          deallocate(li1_b_vadp, li2_b_vadp)
       end if
       is_block_grid_vx_initialized_flag = .false.
    end if
    if(is_block_grid_wx_initialized_flag) then
       ! deallocation of x, mu block-structured grid arrays
       deallocate(nw0_adp, nx0_wadp)
       deallocate(lm0_adp, li0_wadp)
       deallocate(lw0_adp, lx0_wadp)
       deallocate(lbw_adp, lm1_adp, lm2_adp, ubw_adp)
       deallocate(lbi_wadp, li1_wadp, li2_wadp, ubi_wadp)
       deallocate(pi1_wadp, pi2_wadp)
       deallocate(li1_b_wadp, li2_b_wadp)
       is_block_grid_wx_initialized_flag = .false.
    end if

    if (allocated(blk_mks_r_v)) deallocate(blk_mks_r_v)
    if (allocated(blk_mks_v)) deallocate(blk_mks_v)
    if (allocated(blk_mks_r_w)) deallocate(blk_mks_r_w)
    if (allocated(blk_mks_w)) deallocate(blk_mks_w)

  end subroutine finalize_discretization_adptv

end module discretization_adptv_module
