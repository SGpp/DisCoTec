// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

#include <boost/serialization/export.hpp>
#include <filesystem>
#include <string>
#include <vector>

#include "discotec/combischeme/CombiMinMaxScheme.hpp"
#include "discotec/combischeme/CombiThirdLevelScheme.hpp"
#include "discotec/io/BroadcastParameters.hpp"
#include "discotec/manager/ProcessGroupWorker.hpp"
#include "discotec/sparsegrid/DistributedSparseGridIO.hpp"

using namespace combigrid;

DistributedSparseGridUniform<CombiDataType>* schemeFileToSparseGridAndSizesFile(
    const std::string& schemeFileName, DimType dim, const LevelVector& lmin,
    const LevelVector& lmax, const LevelVector& reducedLmax,
    const std::vector<BoundaryType>& boundary, const std::vector<int>& p,
    const LevelVectorList& decomposition) {
  // read in the scheme
  std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
      new CombiMinMaxSchemeFromFile(dim, lmin, lmax, schemeFileName));
  const auto& allLevels = scheme->getCombiSpaces();

  // create sparse grid
  auto uniDSG = new DistributedSparseGridUniform<CombiDataType>(dim, reducedLmax, lmin,
                                                                theMPISystem()->getLocalComm());

  // register levels from other CT scheme
  for (const auto& l : allLevels) {
    auto dfgDecomposition = combigrid::downsampleDecomposition(decomposition, lmax, l, boundary);
    auto uniDFG = std::unique_ptr<OwningDistributedFullGrid<CombiDataType>>(
        new OwningDistributedFullGrid<CombiDataType>(dim, l, theMPISystem()->getLocalComm(),
                                                     boundary, p, false, dfgDecomposition));
    // MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "registering " << l << std::endl;
    uniDSG->registerDistributedFullGrid(*uniDFG);
  }
  auto numDOF =
      std::accumulate(uniDSG->getSubspaceDataSizes().begin(), uniDSG->getSubspaceDataSizes().end(),
                      static_cast<size_t>(0), std::plus<size_t>());
  MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "other sparse grid has " << numDOF
                                             << " DOF per rank, which is " << numDOF * 8. / 1e9
                                             << " GB" << std::endl;
  // write the resulting sparse grid sizes to file
  DistributedSparseGridIO::writeSubspaceSizesToFile(
      *uniDSG, schemeFileName.substr(
                   0, schemeFileName.length() - std::string("_00008groups.json").length()) +
                   ".sizes");
  return uniDSG;
}

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);
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
  std::string thirdCtschemeFile = cfg.get<std::string>("subspace.thirdctscheme", "");
  DimType numSystems = cfg.get<DimType>("thirdLevel.numSystems", 2);

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
  CombiParameters params(dim, lmin, lmax, boundary, 2, 1, CombinationVariant::sparseGridReduce, p,
                         LevelVector(dim, 0), reduceCombinationDimsLmax, 1, forwardDecomposition);
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

  LevelVector reducedLmax = lmax;
  for (DimType d = 0; d < dim; ++d) reducedLmax[d] -= reduceCombinationDimsLmax[d];

  // make local communicator cartesian
  std::vector<int> periods(dim);
  int reorder = false;
  MPI_Comm new_communicator;
  MPI_Cart_create(theMPISystem()->getLocalComm(), dim, p.data(), periods.data(), reorder,
                  &new_communicator);
  theMPISystem()->storeLocalComm(new_communicator);

  std::string conjointSubspaceFileName =
      ctschemeFile.substr(0,
                          ctschemeFile.length() - std::string("_part0_00008groups.json").length()) +
      "conjoint.sizes";

  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> uniDSGs;
  {
    // read in first CT scheme
    std::unique_ptr<CombiMinMaxSchemeFromFile> scheme(
        new CombiMinMaxSchemeFromFile(dim, lmin, lmax, ctschemeFile));

    // generate distributed sparse grid
    uniDSGs.emplace_back(schemeFileToSparseGridAndSizesFile(
        ctschemeFile, dim, lmin, lmax, reducedLmax, boundary, p, decomposition));

    MIDDLE_PROCESS_EXCLUSIVE_SECTION {
      std::cout << "sparse grid contains " << uniDSGs[0]->getNumSubspaces() << " subspaces."
                << std::endl;
    }
  }
  if (otherCtschemeFile == "") {
    // we are done already
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "wrote subspace file" << std::endl;
    return 0;
  }
  {
    // read in second CT scheme
    uniDSGs.emplace_back(schemeFileToSparseGridAndSizesFile(
        otherCtschemeFile, dim, lmin, lmax, reducedLmax, boundary, p, decomposition));

    // read in third CT scheme if applicable -- this could be a loop for > 3 systems
    if (numSystems > 2) {
      uniDSGs.emplace_back(schemeFileToSparseGridAndSizesFile(
          thirdCtschemeFile, dim, lmin, lmax, reducedLmax, boundary, p, decomposition));
    }

    // create conjoint sparse grid
    std::unique_ptr<DistributedSparseGridUniform<CombiDataType>> conjointDSG(
        new DistributedSparseGridUniform<CombiDataType>(dim, reducedLmax, lmin,
                                                        theMPISystem()->getLocalComm()));

    // iterate all uniDSGs and compare their subspace sizes -- if more than one has > 0, add it to
    // conjointDSG
    for (auto i = 0; i < conjointDSG->getNumSubspaces(); ++i) {
      SubspaceSizeType size = 0;
      int numUniDSGsWithSubspaceSet = 0;
      for (const auto& uniDSG : uniDSGs) {
        auto thisUniDSGsSize = uniDSG->getSubspaceDataSizes()[i];
        assert(thisUniDSGsSize == 0 || size == 0 || thisUniDSGsSize == size);
        if (thisUniDSGsSize > 0) {
          ++numUniDSGsWithSubspaceSet;
          size = thisUniDSGsSize;
        }
      }
      if (numUniDSGsWithSubspaceSet > 1) {
        conjointDSG->setDataSize(i, size);
      }
    }

    auto numDOFconjoint = std::accumulate(conjointDSG->getSubspaceDataSizes().begin(),
                                          conjointDSG->getSubspaceDataSizes().end(),
                                          static_cast<size_t>(0), std::plus<size_t>());
    MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "conjoint sparse grid has " << numDOFconjoint
                                               << " DOF per rank" << std::endl;
    // write final sizes to file
    DistributedSparseGridIO::writeSubspaceSizesToFile(*conjointDSG, conjointSubspaceFileName);

    // output first rank's sizes
  }
  MIDDLE_PROCESS_EXCLUSIVE_SECTION std::cout << "wrote conjoint subspace file" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  combigrid::Stats::finalize();

  return 0;
}
