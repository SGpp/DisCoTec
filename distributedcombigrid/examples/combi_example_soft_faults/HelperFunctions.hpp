/*
 * HelperFunctions.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: sccs
 */
#ifndef HELPERFUNCTIONS_HPP_
#define HELPERFUNCTIONS_HPP_

#include "mpi.h"
#include <vector>

namespace combigrid {

void createCommunicators( size_t ngroup, size_t nprocs, int grank, int gsize,
    int& key, int& managerID, MPI_Comm& gcomm, MPI_Comm& lcomm){
  /* determine global rank of each process
   * the manager process always has the highest rank
   * all other processes are worker processes */


  /* create a local communicator for each process group
   * lcomm is the local communicator of its own process group for each worker process
   * for manager, lcomm is a group which contains only manager process and can be ignored
   */
  int color = grank / nprocs;
  key = grank - color*nprocs;
  MPI_Comm_split(MPI_COMM_WORLD, color, key, &lcomm);

  // create global communicator containing manager and pgroup roots
  MPI_Group worldGroup;
  MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);

  std::vector<int> ranks(ngroup+1);
  for (size_t i = 0; i < ngroup; ++i) {
    ranks[i] = i*nprocs;
  }
  ranks.back() = managerID;

  MPI_Group rootGroup;
  MPI_Group_incl(worldGroup, (int)ranks.size(), &ranks[0], &rootGroup);

  MPI_Comm_create(MPI_COMM_WORLD, rootGroup, &gcomm);
}

void readParameterFile(const std::string& fileName, size_t &ngroup, size_t &nprocs, DimType &dim,
    LevelVector &lmin, LevelVector &lmax, LevelVector &leval,
    IndexVector &p, double &time_step, double &time_start,
    double &time_end, size_t &ncombi, FaultsInfo &faultsInfo, bool &plot ){

  // parser for the ini parameter file

  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(fileName, cfg);
  // there are ngroup*nprocs+1 processes needed
  ngroup = cfg.get<size_t>("manager.ngroup");
  nprocs = cfg.get<size_t>("manager.nprocs");

  time_step = cfg.get<double>("simulation.time_step");
  time_start = cfg.get<double>("simulation.time_start");
  time_end = cfg.get<double>("simulation.time_end");
  ncombi = cfg.get<double>("simulation.ncombi");
  plot = cfg.get<bool>("simulation.plot");

  dim = cfg.get<DimType>("ct.dim");

  lmin.resize(dim);
  lmax.resize(dim);
  leval.resize(dim);
  p.resize(dim);

  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  cfg.get<std::string>("ct.leval") >> leval;
  cfg.get<std::string>("ct.p") >> p;

  faultsInfo.numFaults_ = cfg.get<int>("faults.num_faults");

  if ( faultsInfo.numFaults_ == 0 ){
    faultsInfo.iterationFaults_.resize(1);
    faultsInfo.taskFaults_.resize(1);
  } else {
    faultsInfo.iterationFaults_.resize(faultsInfo.numFaults_);
    faultsInfo.taskFaults_.resize(faultsInfo.numFaults_);
  }

  faultsInfo.sdcIndex_.resize(dim);

  cfg.get<std::string>("faults.iteration_faults") >> faultsInfo.iterationFaults_;
  cfg.get<std::string>("faults.task_faults") >> faultsInfo.taskFaults_;
  cfg.get<std::string>("faults.sdc_index") >> faultsInfo.sdcIndex_;

  int sdcMag;
  sdcMag = cfg.get<int>("faults.sdc_mag");
  switch( sdcMag ){
    case 1 : faultsInfo.sdcMag_ = 1e-300; break;
    case 2 : faultsInfo.sdcMag_ = pow(10,-0.5); break;
    case 3 : faultsInfo.sdcMag_ = 1e150; break;
  }

  faultsInfo.sdcMethod_ = cfg.get<int>("method.sdc_method");

  // parameter for gnuplot
  // TODO: use VTK
  std::ofstream paramfile;
  paramfile.open("param.plt");
  paramfile << "time_step = "  << time_step << std::endl;
  paramfile << "time_start = " << time_start << std::endl;
  paramfile << "time_end = "   << time_end << std::endl;
  paramfile << "dim = "        << dim << std::endl;
  paramfile.close();
}

// TODO: redundant definition of exact solution,
//       use solution from Problem.hh
double exact(std::vector<double>& x, double t)
{
  const int dim = x.size();

  double exponent = 0;
  for (int d = 0; d < dim; d++) {
    x[d] -= 0.5 * t;
    exponent -= std::pow(x[d]-0.5,2);
  }

  return std::exp(exponent*100.0);
}

void writeSolutionToFile(std::ofstream& outFile, const FullGrid<CombiDataType>& fg){

  DimType dim = fg.getDimension();
  if ( dim == 2 ) {
    std::vector<double> coords(dim, 0.0);

    for (int i = 0; i < fg.getNrElements(); ++i) {
      if (i % fg.length(0) == 0 && i > 0) {
        outFile << std::endl;
      }
      fg.getCoords(i, coords);
      outFile << coords[0] << "\t"
          << coords[1] << "\t"
          << fg.getElementVector()[i] << std::endl;
    }
    outFile << std::endl << std::endl;
  }
}

void writeErrorToFile(std::ofstream& outFile, FullGrid<CombiDataType>& fg,
    const double &time_step, const double &step){

  // calculate error
  std::vector<double> coords(fg.getDimension(), 0.0);
  for (int i = 0; i < fg.getNrElements(); ++i) {
    fg.getCoords(i, coords);
    fg.getElementVector()[i] -= exact(coords, time_step*step);
  }

  // output for approximation error in gnuplot
  outFile << time_step*step << "\t"
      << fg.getlpNorm(0) << "\t"
      << fg.getlpNorm(2) << std::endl;
}

void checkProcs(const size_t &nprocs, const IndexVector &p){
  IndexType check = 1;
  for (auto k : p)
    check *= k;
  assert(check == IndexType(nprocs));
}

} // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
