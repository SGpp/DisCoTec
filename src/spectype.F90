#include "intrinsic_sizes.h"
#include "redef.h"
!>Definition of a species type and corresponding routines
MODULE spectype_mod
  IMPLICIT NONE

  ! constants
  REAL, PRIVATE,parameter :: pi= 3.141592653589793239d0

  TYPE spectype
     character(len=10) :: name  ! identifier for each species
     logical :: passive         ! back reaction on fields or tracer

     integer :: charge   ! signed charge in units of |e|
     real :: mass
     real :: temp
     real :: dens
     real :: omn ! L_ref/L_n
     real :: omt ! L_ref/L_T

     !temperature and density profile information for global code
     real :: kappa_T
     real :: LT_center
     real :: LT_width
     real :: kappa_n
     real :: Ln_center
     real :: Ln_width
     real :: delta_x_T
     real :: delta_x_n
     integer :: prof_type

     !source term information for global code
     real, DIMENSION(0:4) :: src_amp       !source amplitude
     real, DIMENSION(0:4) :: src_width     !source width
     real, DIMENSION(0:4) :: src_x0        !source center
     real, DIMENSION(0:4) :: src_prof_type !source profile type

     ! profile information for global code
     real, dimension(:),allocatable :: temp_prof,dens_prof
     real, dimension(:),allocatable :: omt_prof,omn_prof

     ! profile input file
     character(len=FILENAME_MAX) :: prof_file

  END TYPE spectype

  INTERFACE ASSIGNMENT (=)
     module procedure spec_to_spec
!     SUBROUTINE spec_to_spec(left,right)
!       TYPE(spectype), INTENT(OUT) :: left
!       TYPE(spectype), INTENT(IN)  :: right
!     END SUBROUTINE spec_to_spec
  END INTERFACE
     
CONTAINS
  FUNCTION static_size_of_spectype() RESULT(memory_need)
    integer :: memory_need

    memory_need = 10*SIZE_OF_CHAR &
         & + SIZE_OF_LOGICAL &
         & + 2*SIZE_OF_INTEGER &
         & + 34*SIZE_OF_REAL
  END FUNCTION static_size_of_spectype
  
  SUBROUTINE spec_to_spec(left,right)
    TYPE(spectype), INTENT(OUT) :: left
    TYPE(spectype), INTENT(IN)  :: right
    left%name = right%name
    left%passive = right%passive
    left%charge = right%charge
    left%mass = right%mass
    left%temp = right%temp
    left%dens = right%dens
    left%omn = right%omn
    left%omt = right%omt


    left%kappa_T   = right%kappa_T
    left%LT_center = right%LT_center
    left%LT_width  = right%LT_width
    left%delta_x_T = right%delta_x_T

    left%kappa_n   = right%kappa_n
    left%Ln_center = right%Ln_center
    left%Ln_width  = right%Ln_width
    left%delta_x_n = right%delta_x_n

    left%prof_type = right%prof_type

    left%src_amp   = right%src_amp
    left%src_width = right%src_width
    left%src_x0    = right%src_x0
    left%src_prof_type = right%src_prof_type
    left%prof_file = right%prof_file

    IF (.NOT.ALLOCATED(left%temp_prof)) THEN
       ALLOCATE(left%temp_prof(LBOUND(right%temp_prof,1):UBOUND(right%temp_prof,1)))
    END IF
    left%temp_prof = right%temp_prof
    IF (.NOT.ALLOCATED(left%dens_prof)) THEN
       ALLOCATE(left%dens_prof(LBOUND(right%dens_prof,1):UBOUND(right%dens_prof,1)))
    END IF
    left%dens_prof = right%dens_prof
    IF (.NOT.ALLOCATED(left%omt_prof)) THEN
       ALLOCATE(left%omt_prof(LBOUND(right%omt_prof,1):UBOUND(right%omt_prof,1)))
    END IF
    left%omt_prof = right%omt_prof
    IF (.NOT.ALLOCATED(left%omn_prof)) THEN
       ALLOCATE(left%omn_prof(LBOUND(right%omn_prof,1):UBOUND(right%omn_prof,1)))
    END IF
    left%omn_prof = right%omn_prof

  END SUBROUTINE spec_to_spec

  SUBROUTINE set_spec_nonprof(sp,a_name,a_passive,a_charge,a_mass, a_temp, &
       & a_dens, a_omn, a_omt, a_kappa_T,a_LT_center, a_LT_width, &
       & a_kappa_n, a_Ln_center, a_Ln_width, a_delta_x_T, a_delta_x_n,&
       & a_prof_type, a_src_amp, a_src_width,a_src_x0, a_src_prof_type,&
       & a_prof_file)
    type(spectype) :: sp
    CHARACTER(*) :: a_name, a_prof_file
    logical :: a_passive    
    integer :: a_charge
    REAL :: a_mass, a_dens, a_omn,a_omt,a_temp
    !temperature and density profile information for global code
    REAL :: a_kappa_T,a_LT_center, a_LT_width, a_kappa_n
    REAL :: a_Ln_center, a_Ln_width, a_delta_x_T, a_delta_x_n
    integer :: a_prof_type
    !source term information for global code
    REAL, DIMENSION(0:4):: a_src_amp, a_src_width,a_src_x0, a_src_prof_type

    sp%name=a_name
    sp%passive=a_passive
    sp%omn    = a_omn
    sp%omt    = a_omt
    sp%mass   = a_mass
    sp%charge = a_charge
    sp%temp   = a_temp
    sp%dens   = a_dens

    sp%kappa_T = a_kappa_T
    sp%LT_center= a_LT_center
    sp%LT_width= a_LT_width
    sp%kappa_n = a_kappa_n
    sp%Ln_center= a_Ln_center
    sp%Ln_width= a_Ln_width
    sp%Ln_width= a_Ln_width
    sp%delta_x_T= a_delta_x_T
    sp%delta_x_n= a_delta_x_n
    sp%prof_type= a_prof_type

    sp%src_amp   = a_src_amp
    sp%src_width = a_src_width
    sp%src_x0    = a_src_x0
    sp%src_prof_type = a_src_prof_type
    sp%prof_file = a_prof_file

  END SUBROUTINE set_spec_nonprof


  SUBROUTINE PRINT(sp)
    type(spectype) :: sp

    WRITE(*,"(3A)") "Species ",TRIM(sp%name)," has the following properties:"
    IF (sp%passive) THEN
       WRITE(*,"(A)") "passive"
    ELSE
       WRITE(*,"(A)") "active"
    END IF
    WRITE(*,"(A,I3)")   "charge = ", sp%charge
    write(*,"(A,F9.4)") "mass   = ", sp%mass
    write(*,"(A,F5.2)") "temp   = ", sp%temp
    WRITE(*,"(A,F5.2)") "dens   = ", sp%dens
    write(*,"(A,F5.2)") "omn    = ", sp%omn
    write(*,"(A,F5.2)") "omt    = ", sp%omt
       if ((sp%kappa_T.ne.-1.0).and.(sp%kappa_n.ne.-1.0)) Then
          Write(*,"(A,F9.4)") "kappa_T  = ", sp%kappa_T
          Write(*,"(A,F9.4)") "LT_center  = ", sp%LT_center
          Write(*,"(A,F9.4)") "LT_width  = ", sp%LT_width
          Write(*,"(A,F9.4)") "kappa_n  = ", sp%kappa_n
          Write(*,"(A,F9.4)") "Ln_center  = ", sp%Ln_center
          Write(*,"(A,F9.4)") "Ln_width  = ", sp%Ln_width
       end if
  END SUBROUTINE PRINT

  
END MODULE spectype_mod
  
