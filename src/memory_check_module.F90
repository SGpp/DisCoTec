module memory_checker
  use, intrinsic :: iso_c_binding
  implicit none

  interface 
     function get_free_memory() bind(c)
       use, intrinsic :: iso_c_binding
       integer( kind = c_long ) :: get_free_memory
     end function get_free_memory
     
     subroutine get_malloc_stat(max_all,unfreed_blocks,still_all) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=c_long) :: max_all
       integer(kind=c_int) :: unfreed_blocks
       integer(kind=c_long) :: still_all
     end subroutine get_malloc_stat

     subroutine show_actual_blocks() bind(c)
     end subroutine show_actual_blocks

     function get_allocated_memory() bind(c)
       use,intrinsic :: iso_c_binding
       integer(kind=c_long) :: get_allocated_memory
     end function get_allocated_memory
  end interface
end module memory_checker

