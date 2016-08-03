/*
 *  Sim_FT_allreduce.cpp
 *
 *  Created on: 11.10.2015
 *      Author: Johannes Walter
 */

#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>


int simft::Sim_FT_MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, simft::Sim_FT_MPI_Comm f_comm){
	simft::Sim_FT_decide_kill();
	simft::Sim_FT_Perform_background_operations();

	int ErrLayer = simft::Sim_FT_Check_dead_processes(f_comm);
	if( MPI_ERR_PROC_FAILED == ErrLayer ){
		return MPI_ERR_PROC_FAILED;
	}else if( MPI_ERR_REVOKED == ErrLayer ){
		return MPI_ERR_REVOKED;
	}

	return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, f_comm->c_comm);
}
