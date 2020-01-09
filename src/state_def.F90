module state_def_mod
  implicit none
  
  type, abstract, public :: state_t
     ! some flags showing what is consistent with g
     logical :: fields_are_consistent = .false.
     logical :: gyro_fields_are_consistent = .false.
  end type state_t
  
end module state_def_mod

module complex_state
  use state_def_mod
  implicit none

  type, public, extends(state_t) :: complex_state_t
     ! distribution function
     complex, dimension(:,:,:,:,:,:), allocatable :: g
     ! electromagnetic fields
     complex, dimension(:,:,:,:), allocatable :: emfields
     ! gyro-averaged fields
  end type complex_state_t
end module complex_state

module real_state
  use state_def_mod
  implicit none

  type, public, extends(state_t) :: real_state_t
     ! distribution function
     real, dimension(:,:,:,:,:,:), allocatable :: g
     ! electromagnetic fields
     real, dimension(:,:,:,:), allocatable :: emfields
     ! gyro-averaged fields
  end type real_state_t
end module real_state
