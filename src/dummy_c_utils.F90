FUNCTION fortran_mallinfo(arena, ordblks, uordblks, fordblks, hblkhd) RESULT(res)
  IMPLICIT NONE
  INTEGER :: arena,ordblks, uordblks, fordblks,hblkhd,res
  arena = 0
  ordblks = 0
  uordblks = 0
  fordblks = 0
  hblkhd = 0
  res = 0
END FUNCTION fortran_mallinfo

function get_max_used_memory(max_used_mem) RESULT(res)
  implicit none
  integer(8) :: max_used_mem
  integer :: res
  
  res = 0
end function get_max_used_memory


function get_datasegment_size() RESULT(res)
  implicit none
  integer :: res

  res=0

end function get_datasegment_size

function short_wait() result(res)
  implicit none
  integer:: res

#if defined(AIX) || defined(BGP)
  call sleep_(1)
#else
  call sleep(1) 
#endif
  res=0
end function short_wait

#if defined(WITHSLEPC)
#define SLEPC_VERSION_DEFAULT -1
integer function get_slepc_version_major()
  get_slepc_version_major = SLEPC_VERSION_DEFAULT
end function get_slepc_version_major
integer function get_slepc_version_minor()
  get_slepc_version_minor = SLEPC_VERSION_DEFAULT
end function get_slepc_version_minor
integer function get_slepc_version_subminor()
  get_slepc_version_subminor = SLEPC_VERSION_DEFAULT
end function get_slepc_version_subminor
integer function get_slepc_version_patch()
  get_slepc_version_patch = SLEPC_VERSION_DEFAULT
end function get_slepc_version_patch
#endif
