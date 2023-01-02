/*
 * HelperFunctions.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: sccs
 */
#ifndef HELPERFUNCTIONS_HPP_
#define HELPERFUNCTIONS_HPP_

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include "mpi.h"
#include <vector>

namespace combigrid {

void createCommunicators( size_t ngroup, size_t nprocs, int globalID, int globalSize,
    int& managerIDgcomm, int& grank, int& lrank, MPI_Comm& gcomm, MPI_Comm& lcomm){
  /* determine global rank of each process
   * the manager process always has the highest rank
   * all other processes are worker processes */


  /* create a local communicator for each process group
   * lcomm is the local communicator of its own process group for each worker process
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  int color = globalID / nprocs;
  int key = globalID - color * nprocs;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);
  MPI_Comm_rank(lcomm, &lrank);
  const int managerIDworld = globalSize - 1;

  /* create global communicator which contains only the manager and the master
   * process of each process group
   * the master processes of the process groups are the processes which have
   * rank 0 in lcomm
   * this communicator is used for communication between master processes of the
   * process groups and the manager and the master processes to each other
   */
  MPI_Group worldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

  std::vector<int> ranks(ngroup + 1);
  for (size_t i = 0; i < ngroup; i++) {
    ranks[i] = i * nprocs;
  }
  ranks.back() = managerIDworld;

  MPI_Group rootGroup;
  MPI_Group_incl(worldGroup, (int) ranks.size(), &ranks[0], &rootGroup);

  MPI_Comm_create(MPI_COMM_WORLD, rootGroup, &gcomm);

  managerIDgcomm = MPI_PROC_NULL;
  if (gcomm != MPI_COMM_NULL) {
    int gcommSize;
    MPI_Comm_size(gcomm, &gcommSize);
    managerIDgcomm = gcommSize - 1;
    MPI_Comm_rank(gcomm, &grank);
  }
}

void readParameterFile(const std::string &fileName, size_t &ngroup, size_t &nprocs, DimType &dim,
                       LevelVector &lmin, LevelVector &lmax, LevelVector &leval,
                       std::vector<int> &p, size_t &ncombi, combigrid::real &dt, size_t &nsteps,
                       FaultsInfo &faultsInfo) {
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(fileName, cfg);

  ngroup = cfg.get<size_t>("manager.ngroup");
  nprocs = cfg.get<size_t>("manager.nprocs");

  dim = cfg.get<DimType>("ct.dim");

  lmin.resize(dim);
  lmax.resize(dim);
  leval.resize(dim);
  p.resize(dim);

  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  cfg.get<std::string>("ct.leval") >> leval;
  cfg.get<std::string>("ct.p") >> p;
  ncombi = cfg.get<size_t>("ct.ncombi");

  dt = cfg.get<combigrid::real>("application.dt");
  nsteps = cfg.get<size_t>("application.nsteps");

  faultsInfo.numFaults_ = cfg.get<int>("faults.num_faults");

  faultsInfo.iterationFaults_.resize(faultsInfo.numFaults_);
  faultsInfo.globalRankFaults_.resize(faultsInfo.numFaults_);

if( faultsInfo.numFaults_ > 0 ){
  cfg.get<std::string>("faults.iteration_faults") >> faultsInfo.iterationFaults_;
  cfg.get<std::string>("faults.global_rank_faults") >> faultsInfo.globalRankFaults_;
}
}

void writeSolutionToFile(std::ofstream& outFile, const FullGrid<CombiDataType>& fg_eval){

  DimType dim = fg_eval.getDimension();

  std::vector<double> coords(dim, 0.0);

  for (int i = 0; i < fg_eval.getNrElements(); i++) {
    if (i % fg_eval.length(0) == 0 && i > 0) {
      outFile << std::endl;
    }

    fg_eval.getCoords(i, coords);
    outFile << coords[0] << "\t" << coords[1] << "\t"
           << fg_eval.getElementVector()[i] << std::endl;
  }

  outFile << std::endl << std::endl;
}

} // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
