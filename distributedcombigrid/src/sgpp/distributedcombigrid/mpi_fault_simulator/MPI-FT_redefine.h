/*
 *  MPI-FT_redefine.h
 *
 *  Created on: 31.07.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

#ifndef MPI_FT_REDEFINE_H_
#define MPI_FT_REDEFINE_H_

static_assert( false, "if you include this header file, it is very probable that you are doing something wrong. onyl remove this assert if you know what you are doing");

#ifdef SIMFTMPI
#undef MPI_COMM_WORLD
#define MPI_COMM_WORLD simft::Sim_FT_MPI_COMM_WORLD

#undef MPI_COMM_NULL
#define MPI_COMM_NULL simft::Sim_FT_MPI_COMM_NULL

#undef MPI_STATUS_IGNORE
//#define MPI_STATUS_IGNORE simft::Sim_FT_MPI_STATUS_IGNORE

#undef MPI_STATUSES_IGNORE
//#define MPI_STATUSES_IGNORE simft::Sim_FT_MPI_STATUSES_IGNORE

#define MPI_STATUSES_IGNORE (simft::Sim_FT_MPI_Status *)1
#define MPI_STATUS_IGNORE (simft::Sim_FT_MPI_Status *)1

#define MPI_Abort simft::Sim_FT_MPI_Abort
#define MPI_Allreduce simft::Sim_FT_MPI_Allreduce
#define MPI_Barrier simft::Sim_FT_MPI_Barrier
#define MPI_Bcast simft::Sim_FT_MPI_Bcast
#define MPI_Cart_coords simft::Sim_FT_MPI_Cart_coords
#define MPI_Cart_create simft::Sim_FT_MPI_Cart_create
#define MPI_Comm_create simft::Sim_FT_MPI_Comm_create
#define MPI_Comm_free simft::Sim_FT_MPI_Comm_free
#define MPI_Cart_rank simft::Sim_FT_MPI_Cart_rank
#define MPI_Comm simft::Sim_FT_MPI_Comm
//#define MPI_Comm_create_errhandler simft::Sim_FT_MPI_Comm_create_errhandler
#define MPI_Comm_group simft::Sim_FT_MPI_Comm_group
#define MPI_Comm_rank simft::Sim_FT_MPI_Comm_rank
#define MPI_Comm_set_errhandler simft::Sim_FT_MPI_Comm_set_errhandler
#define MPI_Comm_size simft::Sim_FT_MPI_Comm_size
#define MPI_Comm_split simft::Sim_FT_MPI_Comm_split
#define MPI_Finalize simft::Sim_FT_MPI_Finalize
#define MPI_Get_count simft::Sim_FT_MPI_Get_count
#ifndef DISABLE_NBC
#define MPI_Iallreduce simft::Sim_FT_MPI_Iallreduce
#define MPI_Ibarrier simft::Sim_FT_MPI_Ibarrier
#define MPI_Ibcast simft::Sim_FT_MPI_Ibcast
#endif
#define MPI_Irecv simft::Sim_FT_MPI_Irecv
#define MPI_Isend simft::Sim_FT_MPI_Isend
#define MPI_Probe simft::Sim_FT_MPI_Probe
#define MPI_Reduce simft::Sim_FT_MPI_Reduce
#define MPI_Request simft::Sim_FT_MPI_Request
#define MPI_Request_free simft::Sim_FT_MPI_Request_free
#define MPI_Status simft::Sim_FT_MPI_Status
#define MPI_Test simft::Sim_FT_MPI_Test
#define MPI_Wait simft::Sim_FT_MPI_Wait
#define MPI_Waitall simft::Sim_FT_MPI_Waitall

#define MPI_Init simft::Sim_FT_MPI_Init
#define MPI_Recv simft::Sim_FT_MPI_Recv
#define MPI_Send simft::Sim_FT_MPI_Send

#endif

#endif /* MPI_FT_REDEFINE_H_ */
