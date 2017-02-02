/*
int arena
    This is the total size of memory allocated with sbrk by malloc, in bytes.
int ordblks
    This is the number of chunks not in use. (The memory allocator internally 
    gets chunks of memory from the operating system, and then carves them up 
    to satisfy individual malloc requests; see Efficiency and Malloc.)
int smblks
    This field is unused.
int hblks
    This is the total number of chunks allocated with mmap.
int hblkhd
    This is the total size of memory allocated with mmap, in bytes.
int usmblks
    This field is unused.
int fsmblks
    This field is unused.
int uordblks
    This is the total size of memory occupied by chunks handed out by malloc.
int fordblks
    This is the total size of memory occupied by free (not in use) chunks.
int keepcost
    This is the size of the top-most releasable chunk that normally borders the end of the heap (i.e., the high end of the virtual address space's data segment). 
*/

#include <malloc.h>
#include <unistd.h>
#if defined(WITHSLEPC)
#include "slepcversion.h"
#endif

/*int fortran_mallinfo(int *arena, int *ordblks, int *smblks, int *hblks, int *hblkhd, int *usmblks, int *fsmblks, int *uordblks, int *fordblks, int *keepcost) {*/

int fortran_mallinfo(int *arena, int *ordblks, int *uordblks, int *fordblks, int *hblkhd) {
  struct mallinfo mi;

  mi = mallinfo();
  *arena = mi.arena;
  *ordblks = mi.ordblks;
  /* *smblks = mi.smblks;
   *hblks =  mi.hblks; */
  *hblkhd = mi.hblkhd;
  /*
  *usmblks = mi.usmblks;
  *fsmblks = mi.fsmblks; */
  *uordblks = mi.uordblks;
  *fordblks = mi.fordblks;
  /* *keepcost = mi.keepcost;*/
  return 0;
}

int fortran_long_mallinfo(long long *arena, long long *ordblks, long long *uordblks, long long *fordblks, long long *hblkhd) {
  struct mallinfo mi;

  unsigned int u_arena,u_ordblks, u_hblkhd, u_uordblks, u_fordblks;


  mi = mallinfo();
  u_arena = mi.arena;
  u_ordblks = mi.ordblks;
  /* *smblks = mi.smblks;
   *hblks =  mi.hblks; */
  u_hblkhd = mi.hblkhd;
  /*
  *usmblks = mi.usmblks;
  *fsmblks = mi.fsmblks; */
  u_uordblks = mi.uordblks;
  u_fordblks = mi.fordblks;

  *arena = u_arena;
  *ordblks = u_ordblks;
  /* *smblks = mi.smblks;
   *hblks =  mi.hblks; */
  *hblkhd = u_hblkhd;
  /*
  *usmblks = mi.usmblks;
  *fsmblks = mi.fsmblks; */
  *uordblks = u_uordblks;
  *fordblks = u_fordblks;
  /* *keepcost = mi.keepcost;*/
  return 0;
}

int fortran_mallinfo_(int *arena, int *ordblks, int *uordblks, int *fordblks, int *hblkhd) {
  return fortran_mallinfo(arena, ordblks, uordblks, fordblks, hblkhd);
}

int fortran_long_mallinfo_(long long *arena, long long *ordblks, long long *uordblks, long long *fordblks, long long *hblkhd) {
  return fortran_long_mallinfo(arena, ordblks, uordblks, fordblks, hblkhd);
}

#if defined(POWER)
unsigned long long get_datasegment_size()
{
   extern char     _data[];
   return sbrk(0);/* - (unsigned long long) _data;*/
}
#else
unsigned long long get_datasegment_size() {
   extern char     _start[], _end[];
   struct mallinfo mi;
   mi = mallinfo();

   return (_end - _start) + mi.arena + mi.hblkhd;
}
#endif

unsigned long long get_datasegment_size_() {
   return get_datasegment_size();
}

int get_max_used_memory(long long *max_used_mem) {
  *max_used_mem = get_datasegment_size();
  return 0;
}

int get_max_used_memory_(long long *max_used_mem) {
  return get_max_used_memory(max_used_mem);
}

int short_wait() {
  unsigned int millisec;
  millisec=1000;
  return usleep(millisec);
}

int short_wait_() {
  unsigned int millisec;
  millisec=1000;
  return usleep(millisec);
}

#if defined(WITHSLEPC)
int get_slepc_version_major() {
  return SLEPC_VERSION_MAJOR;
}
int get_slepc_version_major_() {
  return SLEPC_VERSION_MAJOR;
}
int get_slepc_version_minor() {
  return SLEPC_VERSION_MINOR;
}
int get_slepc_version_minor_() {
  return SLEPC_VERSION_MINOR;
}
int get_slepc_version_subminor() {
  return SLEPC_VERSION_SUBMINOR;
}
int get_slepc_version_subminor_() {
  return SLEPC_VERSION_SUBMINOR;
}
int get_slepc_version_patch() {
  return SLEPC_VERSION_PATCH;
}
int get_slepc_version_patch_() {
  return SLEPC_VERSION_PATCH;
}
#endif
