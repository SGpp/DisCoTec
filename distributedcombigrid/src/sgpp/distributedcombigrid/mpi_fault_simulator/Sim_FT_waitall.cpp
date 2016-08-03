/*
 * Sim_FT_waitall.cpp
 *
 *  Created on: 18.10.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>
#include <unordered_set>


int simft::Sim_FT_MPI_Waitall(int count, simft::Sim_FT_MPI_Request array_of_requests[],
		simft::Sim_FT_MPI_Status array_of_statuses[]){
	simft::Sim_FT_decide_kill();
	simft::Sim_FT_Perform_background_operations();

	std::unordered_set<int> ActiveWait;
	ActiveWait.reserve(count);

	for(int i = 0; i < count; i++){
		ActiveWait.insert(i);
	}

	bool errFlag = false;

	while(ActiveWait.size() > 0){
		for (const auto& id: ActiveWait) {
		    int Flag = 0;
		    int Err;
		    if(array_of_statuses != SIM_FT_MPI_STATUSES_IGNORE){
		    	Err = simft::Sim_FT_MPI_Test(&array_of_requests[id], &Flag, &array_of_statuses[id]);
		    }else{
		    	Err = simft::Sim_FT_MPI_Test(&array_of_requests[id], &Flag, SIM_FT_MPI_STATUS_IGNORE);
		    }
		    if(Flag == 1){
		    	if(MPI_SUCCESS != Err){
		    		errFlag = true;
		    	}

		    	//Wait complete => remove id from wait set
		    	ActiveWait.erase(id);
		    }
		}
	}

	if(errFlag == true){
		return MPI_ERR_IN_STATUS;
	}else{
		return MPI_SUCCESS;
	}

}
