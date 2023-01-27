// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <filesystem>
#include <string>
#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "combischeme/CombiThirdLevelScheme.hpp"
#include "io/H5InputOutput.hpp"
#include "manager/ProcessGroupWorker.hpp"

using namespace combigrid;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  combigrid::Stats::initialize();

  // read in parameter file -- use the same one as for the simulation, but add the other ct scheme
  // (and adapt p/nprocs and spread if you are using fewer processes than in your target scenario)
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

  // only need one process group here
  size_t ngroup = 1;
  size_t nprocs = cfg.get<size_t>("manager.nprocs");

  theMPISystem()->initWorldReusable(MPI_COMM_WORLD, ngroup, nprocs, false);

  /* read other parameters from ctparam */
  DimType dim = cfg.get<DimType>("ct.dim");
  LevelVector lmin(dim), lmax(dim);
  std::vector<int> p(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;
  cfg.get<std::string>("ct.p") >> p;
  std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme");
  std::string otherCtschemeFile = cfg.get<std::string>("subspace.otherctscheme", "");
  auto divide = cfg.get<DimType>("subspace.divide", 1);

  // periodic boundary conditions
  std::vector<BoundaryType> boundary(dim, 1);
  auto forwardDecomposition = false;

  // check whether parallelization vector p agrees with nprocs
  int checkProcs = 1;
  for (auto k : p) checkProcs *= k;
  if (checkProcs != IndexType(nprocs)) {
    throw std::invalid_argument("process group size and parallelization do not match");
  }

  // create combiparameters
  auto reduceCombinationDimsLmax = LevelVector(dim, 1);
  CombiParameters params(dim, lmin, lmax, boundary, 2, 1, p, LevelVector(dim, 0),
                         reduceCombinationDimsLmax, forwardDecomposition);
  IndexVector minNumPoints(dim), maxNumPoints(dim);
  for (DimType d = 0; d < dim; ++d) {
    minNumPoints[d] = combigrid::getNumDofNodal(lmin[d], boundary[d]);
    maxNumPoints[d] = combigrid::getNumDofNodal(lmax[d], boundary[d]);
  }
  // first, test if decomposition possible for small resolution
  auto decomposition = combigrid::getDefaultDecomposition(minNumPoints, p, forwardDecomposition);
  // then assign the actual used one
  decomposition = combigrid::getDefaultDecomposition(maxNumPoints, p, forwardDecomposition);
  // default decomposition works only for powers of 2!
  params.setDecomposition(decomposition);
  MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "generated parameters" << std::endl;

  // make local communicator cartesian
  std::vector<int> periods(dim);
  int reorder = false;
  MPI_Comm new_communicator;
  MPI_Cart_create(theMPISystem()->getLocalComm(), dim, p.data(), periods.data(), reorder,
                  &new_communicator);
  theMPISystem()->storeLocalComm(new_communicator);
  MPI_Barrier(MPI_COMM_WORLD);
  MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "stored cartesian communicator" << std::endl;
  std::string firstSubspaceFileName =
      ctschemeFile.substr(0, ctschemeFile.length() - std::string(".json").length()) + ".sizes";
  std::string conjointSubspaceFileName =
      ctschemeFile.substr(0, ctschemeFile.length() - std::string("split1_40groups.json").length()) +
      "conjoint.sizes";
  {
    // read in first CT scheme
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
        new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));
    const auto& allLevels = scheme->getCombiSpaces();

    // generate distributed sparse grid
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(dim, lmax, lmin,
                                                        theMPISystem()->getLocalComm()));
    MIDDLE_PROCESS_EXCLUSIVE_SECTION {
      std::cout << "sparse grid contains " << uniDSG->getNumSubspaces() << " subspaces."
                << std::endl;
    }

    // register all component grid levels in this sparse grid
    for (const auto& l : allLevels) {
      auto dfgDecomposition = combigrid::downsampleDecomposition(decomposition, lmax, l, boundary);
      auto uniDFG = std::unique_ptr<DistributedFullGrid<CombiDataType>>(
          new DistributedFullGrid<CombiDataType>(dim, l, theMPISystem()->getLocalComm(), boundary,
                                                 p, false, dfgDecomposition));
      // MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "registering " << l << std::endl;
      uniDSG->registerDistributedFullGrid(*uniDFG);
    }
    auto numDOF = std::accumulate(uniDSG->getSubspaceDataSizes().begin(),
                                  uniDSG->getSubspaceDataSizes().end(), 0);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "sparse grid has " << numDOF << " DOF per rank"
                                               << std::endl;
    // write the resulting sparse grid sizes to file
    uniDSG->writeSubspaceSizesToFile(firstSubspaceFileName);
  }
  if (otherCtschemeFile == "") {
    // we are done already
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "wrote subspace file" << std::endl;
    return 0;
  }
  {
    // read in second CT scheme
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
        new CombiMinMaxSchemeFromFile(dim, lmin, lmax, otherCtschemeFile));
    const auto& allLevels = scheme->getCombiSpaces();

    // another sparse grid
    auto uniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(dim, lmax, lmin,
                                                        theMPISystem()->getLocalComm()));

    // register levels from other CT scheme
    for (const auto& l : allLevels) {
      auto dfgDecomposition = combigrid::downsampleDecomposition(decomposition, lmax, l, boundary);
      auto uniDFG = std::unique_ptr<DistributedFullGrid<CombiDataType>>(
          new DistributedFullGrid<CombiDataType>(dim, l, theMPISystem()->getLocalComm(), boundary,
                                                 p, false, dfgDecomposition));
      // MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "registering " << l << std::endl;
      uniDSG->registerDistributedFullGrid(*uniDFG);
    }
    auto numDOF = std::accumulate(uniDSG->getSubspaceDataSizes().begin(),
                                  uniDSG->getSubspaceDataSizes().end(), 0);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "other sparse grid has " << numDOF
                                               << " DOF per rank" << std::endl;
    // write the resulting sparse grid sizes to file
    uniDSG->writeSubspaceSizesToFile(
        otherCtschemeFile.substr(0, otherCtschemeFile.length() - std::string(".json").length()) +
        ".sizes");

    // read written sparse grid sizes from file
    // for extra sparse grid / conjoint subspaces, min-reduce the sizes
    auto minFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
      return std::min(a, b);
    };
    // use extra sparse grid
    uniDSG->readReduceSubspaceSizesFromFile(firstSubspaceFileName, minFunctionInstantiation);

    auto numDOFconjoint = std::accumulate(uniDSG->getSubspaceDataSizes().begin(),
                                          uniDSG->getSubspaceDataSizes().end(), 0);
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "conjoint sparse grid has " << numDOFconjoint
                                               << " DOF per rank" << std::endl;
    // write final sizes to file
    uniDSG->writeSubspaceSizesToFile(conjointSubspaceFileName);

    // output first rank's sizes
  }
  MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "wrote conjoint subspace file" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  combigrid::Stats::finalize();
  MPI_Finalize();

  return 0;
}
