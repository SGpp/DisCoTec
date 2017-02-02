#include "redef.h"

module df_nonlinear_c_mod
  use,intrinsic :: iso_c_binding
  use df_nonlinear_term_mod
  use discretization, only: pi1,pi2, li1,li2,lj1,lj2
  implicit none

  private

  type,public,extends(df_nonlinear_term_t) :: df_nonlinear_c_t
   contains
     procedure :: initialize => c_initialize_nonlinearity
     procedure :: finalize => c_finalize_nonlinearity
     !procedure :: add => add_nonlinearity
     procedure,private :: calc => c_calc_nonlinearity
     !procedure :: mem_est => mem_est_nonlinearity
     !procedure :: construct => construct_df_nonlinear_term
     !final :: destruct
     !procedure :: destruct => destruct_df_nonlinear_term
     procedure :: getType => getThisType
  end type df_nonlinear_c_t

  interface
     subroutine c_calc_nonlinearity_df(gy_chi,g2d,vexy,dgdxy,localrhs,cptr_pnl_1d,first) bind(c)
       import
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,1:*),Intent(in) :: gy_chi, g2d
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
       Complex(C_DOUBLE_COMPLEX), Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
       real(C_DOUBLE),dimension(pi1:pi2) :: cptr_pnl_1d
       Logical(C_BOOL), value, Intent(in):: first
     end subroutine c_calc_nonlinearity_df

     !subroutine initialize_fourier_hp(a_li0da,a_ly0da) bind(c)
     !  import
     !  integer(C_INT), value :: a_li0da,a_ly0da
     !end subroutine initialize_fourier_hp

     !subroutine finalize_fourier_hp() bind(c)
     !end subroutine finalize_fourier_hp

     subroutine c_initialize_nonlinearity_df(a_lbg0) bind(c)
       import
       integer(C_INT),value :: a_lbg0
     end subroutine c_initialize_nonlinearity_df

     subroutine c_finalize_nonlinearity_df() bind(c)
     end subroutine c_finalize_nonlinearity_df
  end interface
contains
  function getThisType(this)
    class(df_nonlinear_c_t) :: this
    character(MAX_TYPENAME_LENGTH) :: getThisType

    getThisType = "df_nonlinear_c_t"
  end function getThisType

  subroutine c_initialize_nonlinearity(this)
    class(df_nonlinear_c_t) :: this

    if (this%init_status.eq.1) return
    call this%prefactor%initialize()
    call c_initialize_nonlinearity_df(this%lbg0)
    this%init_status = 1
  end subroutine c_initialize_nonlinearity

  subroutine c_finalize_nonlinearity(this)
    class(df_nonlinear_c_t) :: this

    call this%prefactor%finalize()
    call c_finalize_nonlinearity_df()
    this%init_status = 0
  end subroutine c_finalize_nonlinearity

  subroutine c_calc_nonlinearity(this,gy_chi,g2d,vexy,dgdxy,localrhs,first)
    class(df_nonlinear_c_t) :: this
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(in) :: gy_chi, g2d
    Complex, Dimension(li1:li2,lj1:lj2,2,1:*),Intent(inout) :: vexy, dgdxy
    Complex, Dimension(li1:li2,lj1:lj2,1:*),Intent(inout) :: localrhs  
    Logical, Intent(in) :: first

    !print*,"prefactor status : ",allocated(this%prefactor),allocated(this%prefactor%pnl_1d)

    call c_calc_nonlinearity_df(gy_chi,&
         &g2d,vexy,dgdxy,localrhs,this%prefactor%pnl_1d,logical(first,kind=C_BOOL))
  end subroutine c_calc_nonlinearity

end module df_nonlinear_c_mod
