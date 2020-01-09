#include "intrinsic_sizes.h"
module page_locked_memory_mod
  use, intrinsic :: iso_c_binding
  implicit none

  public :: cuda_register_array,cuda_unregister_array,lock_memory, unlock_memory

  interface
     function mlock(ptr,number_of_bytes) bind(c)
       import
       type(C_PTR), value :: ptr
       integer(C_INT), value :: number_of_bytes
       integer(C_INT) :: mlock
     end function mlock

     function munlock(ptr,number_of_bytes) bind(c)
       import
       type(C_PTR), value :: ptr
       integer(C_INT), value :: number_of_bytes
       integer(C_INT) :: mlock
     end function munlock

     subroutine cuda_register_array(arr,arr_size) bind(c)
       import
       integer(C_INT), value :: arr_size
       Complex(C_DOUBLE_COMPLEX), Dimension(arr_size),Intent(inout) :: arr
     end subroutine cuda_register_array

     subroutine cuda_unregister_array(arr) bind(c)
       import
       Complex(C_DOUBLE_COMPLEX), Dimension(*),Intent(inout) :: arr
     end subroutine cuda_unregister_array
  end interface

  interface lock_memory
     module procedure lock_memory_6D,lock_memory_4D
  end interface

  interface unlock_memory
     module procedure unlock_memory_6D, unlock_memory_4D
  end interface unlock_memory
  private
contains
  subroutine lock_memory_6D(arr)
    complex, dimension(:,:,:,:,:,:) :: arr

    integer(C_INT) :: number_of_elements

    number_of_elements = size(arr)
    call cuda_register_array(arr,number_of_elements)
  end subroutine lock_memory_6D

  subroutine lock_memory_4D(arr)
    complex, dimension(:,:,:,:) :: arr

    integer(C_INT) :: number_of_elements

    number_of_elements = size(arr)
    call cuda_register_array(arr,number_of_elements)
  end subroutine lock_memory_4D

  subroutine unlock_memory_6D(arr)
    complex, dimension(:,:,:,:,:,:) :: arr

    call cuda_unregister_array(arr)
  end subroutine unlock_memory_6D

  subroutine unlock_memory_4D(arr)
    complex, dimension(:,:,:,:) :: arr

    call cuda_unregister_array(arr)
  end subroutine unlock_memory_4D
end module page_locked_memory_mod
    
