/*
 * ProcTaskAssoc.hpp
 *
 *  Created on: 29.05.2018
 *      Author: simon
 */

#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_PROCTASKASSOC_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_PROCTASKASSOC_HPP_

#include <cassert>
#include <vector>
#include <algorithm>
#include <mpi.h>

class ProcTaskAssoc{
public:

	ProcTaskAssoc() : procRank {}, taskID {}{

	}

	ProcTaskAssoc(std::vector<int> procRank, std::vector<int> taskID) : procRank {std::move(procRank)}, taskID {std::move(taskID)}{
		assert(procRank.size() == taskID.size());
	}


	int getProc(int task){
		size_t index = std::find(std::begin(taskID), std::end(taskID), task) - std::begin(taskID);
		return procRank.at(index);
	}

	void broadcast(int sourceRank, MPI_Comm comm){
		size_t size = procRank.size();
		MPI_Bcast(&size, 1, MPI_INT, sourceRank, comm);
		procRank.resize(size);
		taskID.resize(size);
		MPI_Bcast(procRank.data(), size, MPI_INT, sourceRank, comm);
		MPI_Bcast(taskID.data(), size, MPI_INT, sourceRank, comm);
	}

	void send(int dest, MPI_Comm comm){
		size_t size = procRank.size();
		MPI_Send(&size, 1, MPI_INT, dest, messageID, comm);
		procRank.resize(size);
		taskID.resize(size);
		MPI_Send(&procRank.data(), size, MPI_INT, dest, messageID+1, comm);
		MPI_Send(&taskID.data(), size, MPI_INT, dest, messageID+2, comm);
	}

	void receive(int src, MPI_Comm comm){
		size_t size;
		MPI_Recv(&size, 1, MPI_INT, src, messageID, comm, MPI_STATUS_IGNORE);
		procRank.resize(size);
		taskID.resize(size);
		MPI_Recv(&procRank.data(), size, MPI_INT, src, messageID+1, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&taskID.data(), size, MPI_INT, src, messageID+2, comm, MPI_STATUS_IGNORE);
	}

private:
	static constexpr int messageID = 5930;


	std::vector<int> procRank;
	std::vector<int> taskID;
};


#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_PROCTASKASSOC_HPP_ */
