/*
 *  Sim_FT_finalize.cpp
 *
 *  Created on: 27.08.2015
 *      Author: Johannes Walter
 */


#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

int simft::Sim_FT_MPI_Finalize(void){
	simft::Sim_FT_Check_dead_processes(simft::Sim_FT_MPI_COMM_WORLD, simft::NextOp::Default); //in case of revoked comm make sure all processes switched sync tags
	simft::Sim_FT_Check_dead_processes(simft::Sim_FT_MPI_COMM_WORLD, simft::NextOp::Finalize);

	simft::Sim_FT_Perform_background_operations();


	#ifdef USE_NON_BLOCKING_COLLECTIVES
	MPI_Op_free(&customOp);
	#endif

	/* Delete the dynamically created custom communicators. If necessary,
	 * (e.g. the finalize exits with an error), they can additionally be freed using MPI_Comm_free
	 * together with our protocol to free outstanding communication
	 *  */
	for( auto Actives : simft::Sim_FT_Current_Active_Communicators){
		delete Actives;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return MPI_Finalize();
}
