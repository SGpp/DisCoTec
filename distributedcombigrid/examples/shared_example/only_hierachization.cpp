#include <mpi.h>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

//double * chache_flush_pointer;

double fillfunction(std::vector<double>& coords) {
  double result=1;
  for (double d:coords) {
    result*=d*d;
  }
  return result;
}


/**
 * Acts as main function <br>
 * Has to be split up from main to run destructors before MPI_FINALIZE
 * @param comm:
 */
void main2(CommunicatorType comm)
{

  
  int rank;
  int size_comm;
  MPI_Comm_size (comm, &size_comm);
  MPI_Comm_rank (comm, &rank);

  DimType dim=4;

  LevelVector levels = {6, 6, 6,6};//Todo read in
  IndexVector procs = {1,1,1,1};//! product has to be equal to processes used

  std::vector<bool> boundary(4, false);

  DistributedFullGrid<double> dfg(dim, levels, comm, boundary, procs,false);
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);
    dfg.getData()[li] = fillfunction(coords);
  }
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    auto start=std::chrono::high_resolution_clock::now();
    //DistributedHierarchization::hierarchize(dfg);

    auto middle1=std::chrono::high_resolution_clock::now();
    // ? flush cache?
    auto middle2=std::chrono::high_resolution_clock::now();
    DistributedHierarchization::dehierarchize(dfg);
    auto end=std::chrono::high_resolution_clock::now();
    std::cout << "hierachize:   " << std::chrono::duration_cast<std::chrono::milliseconds>(middle1-start).count() <<"ms\n";
    std::cout << "dehierachize: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-middle2).count() <<"ms\n";

  }
  else
  {
    DistributedHierarchization::hierarchize(dfg);
    // ? flush cache?
    DistributedHierarchization::dehierarchize(dfg);
  }
  
  /* 
  WORLD_MANAGER_EXCLUSIVE_SECTION {
    Stats::stopEvent("total time");
  }

  Stats::finalize();
  Stats::write("timeoh.json");*/
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  CommunicatorType comm =MPI_COMM_WORLD;
  main2(comm);

  MPI_Finalize();
}

