/*
 *  Sim_FT_barrier.cpp
 *
 *  Created on: 27.07.2015
 *      Author: Johannes Walter
 */

#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

int simft::Sim_FT_MPI_Barrier(simft::Sim_FT_MPI_Comm f_comm)
{
	simft::Sim_FT_decide_kill();
	simft::Sim_FT_Perform_background_operations();

	int ErrLayer = simft::Sim_FT_Check_dead_processes(f_comm);

	if( MPI_ERR_PROC_FAILED == ErrLayer ){
		return MPI_ERR_PROC_FAILED;
	}else if( MPI_ERR_REVOKED == ErrLayer ){
		return MPI_ERR_REVOKED;
	}

    int err = MPI_Barrier(f_comm->c_comm);

    return err;
}
