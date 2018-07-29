/*
 *  MPI-FT_ULFM.cpp
 *
 *  Created on: 18.10.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE

#include <iostream>



int MPI_Comm_revoke(simft::Sim_FT_MPI_Comm f_comm){
	simft::Sim_FT_Send_revoke_message(f_comm);
	return MPI_SUCCESS;
}

int MPI_Comm_shrink(simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Comm* newcomm){
	simft::Sim_FT_Check_dead_processes(f_comm, simft::NextOp::Default); //in case of revoked comm make sure all processes switched sync tags

	simft::Sim_FT_Check_dead_processes(f_comm, simft::NextOp::Shrink);

	(*newcomm) = new simft::Sim_FT_MPI_Comm_struct;
	int Ret = MPI_Comm_split(f_comm->c_comm, 1, f_comm->comm_rank, &(*newcomm)->c_comm);

	simft::Sim_FT_Initialize_new_comm(newcomm, true);

	return Ret;
}

int MPI_Comm_failure_ack(simft::Sim_FT_MPI_Comm f_comm){
	MPI_Group tempgroup;
	MPI_Comm_group(f_comm->c_comm, &tempgroup);

	if(f_comm->dead_set.size() > 0){

		if(f_comm->failedgrp != MPI_GROUP_NULL){
			//if failedgrp already exists, free its resources
			MPI_Group_free(&f_comm->failedgrp);
		}

		//create a group containing the current dead processes
		std::vector<int> deadP (f_comm->dead_set.begin(),f_comm->dead_set.end());
		MPI_Group_incl(tempgroup, deadP.size(), &deadP[0], &f_comm->failedgrp);
	}

	f_comm->Ack_failed_processes = true;
	return MPI_SUCCESS;
}

int MPI_Comm_failure_get_acked(simft::Sim_FT_MPI_Comm f_comm, MPI_Group* failedgrp){
	*failedgrp = f_comm->failedgrp;
	return MPI_SUCCESS;
}

/*
 * following function is used to agree on the boolean flag, if flag is false at any
 * other(alive) process in the communicator, the value is set to false at ALL processes.
 * This function is blocking collective. That means all alive processes have to participate.
 */
int MPI_Comm_agree(simft::Sim_FT_MPI_Comm f_comm, int* flag){
	//important: if comm is revoked, an initial Default sync ensures, all processes have switched to the revoked sync tags
	simft::Sim_FT_Check_dead_processes(f_comm,simft::NextOp::Default);

	if (*flag == 0){
		simft::Sim_FT_Check_dead_processes(f_comm,simft::NextOp::AFalse);
	}else{
		simft::Sim_FT_Check_dead_processes(f_comm,simft::NextOp::ATrue);
	}

	if(f_comm->Last_Send_Vector[0] == 4){ //This value will be the same at all processes. 4 => false; 5 => true
		*flag = 0;
	}else{
		*flag = 1;
	}

	if(f_comm->dead_nodes.size() == 0){
		return MPI_SUCCESS;
	}else{
		return MPI_ERR_PROC_FAILED;
	}
}
