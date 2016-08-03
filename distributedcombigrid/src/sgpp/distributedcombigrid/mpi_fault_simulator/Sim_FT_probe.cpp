/*
 *  Sim_FT_probe.cpp
 *
 *  Created on: 31.07.2015
 *      Author: Johannes Walter
 */
#include "MPI-FT.h"
#include REAL_MPI_INCLUDE
#include <iostream>

//note this currently only works to probe messages from blocking MPI_Send
int simft::Sim_FT_MPI_Probe(int source, int tag, simft::Sim_FT_MPI_Comm f_comm, simft::Sim_FT_MPI_Status *status){
	simft::Sim_FT_decide_kill();
	simft::Sim_FT_Perform_background_operations();

	//send probe request
	MPI_Request Send_Request;
	char Probe_Req = PROBE_REQUEST;
	MPI_Isend(&Probe_Req, 1, MPI_CHAR, source, tag, f_comm->c_comm_copy_p2p, &Send_Request);//, &Send_Request);
	MPI_Request_free(&Send_Request);
	MPI_Status ProbeStatus, probeRecvStatus;
	int Flag = 0;
	while(Flag == 0){
		simft::Sim_FT_Perform_background_operations();
		MPI_Iprobe(source, tag, f_comm->c_comm_copy_probe, &Flag, &ProbeStatus);
	}
	//the sending process sends the information, how big the sending message is via a special probe communicator
	int ProbeCount = 0;
	MPI_Recv(&ProbeCount, 1, MPI_INT, ProbeStatus.MPI_SOURCE, ProbeStatus.MPI_TAG, f_comm->c_comm_copy_probe, &probeRecvStatus);

	/* Normally we would store the count inside the status object, but on the HLRS HazelHen implementation of MPI, count does not exist,
	 * so we just use our own Status object and store the count there.
	 */
	status->count = ProbeCount;
	status->MPI_ERROR = probeRecvStatus.MPI_ERROR;
	status->MPI_SOURCE = probeRecvStatus.MPI_SOURCE;
	status->MPI_TAG = probeRecvStatus.MPI_TAG;


	return 0;
}
