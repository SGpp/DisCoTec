!>Wrapper for 'use mpi' being required on some machines
!!
!!Required if the MPI library is not directly provided;
!!Typically, this wrapper should be linked to comm.F90
!!(see e.g. BlueGene.mk)
module mpi
#include "mpif.h"
end module mpi
