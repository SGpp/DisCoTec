#pragma once
/**
 * @file ParameterIO.hpp
 * @brief Read, broadcast, and interpret combination technique parameters.
 *
 * Provides:
 *   - Low-level I/O: read files and broadcast from rank 0
 *   - High-level: readParameterFile() + buildCombiParametersFromConfig()
 */

#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "discotec/combischeme/CombiMinMaxScheme.hpp"
#include "discotec/io/BroadcastParameters.hpp"
#include "discotec/manager/CombiParameters.hpp"
#include "discotec/mpi/MPISystem.hpp"
#include "discotec/utils/Types.hpp"

namespace combigrid {

struct ParameterFileContent {
  size_t ngroup;
  size_t nprocs;
  boost::property_tree::ptree cfg;
};

/**
 * @brief Read a ctparam file and extract ngroup/nprocs.
 *
 * If MPI is initialized, rank 0 reads and broadcasts.
 * Otherwise reads locally (useful for pre-MPI-init extraction).
 *
 * @param paramfile  Path to ctparam INI file
 * @param comm       Communicator for broadcast (default: MPI_COMM_WORLD)
 */
inline ParameterFileContent readParameterFile(const std::string& paramfile,
                                              MPI_Comm comm = MPI_COMM_WORLD) {
  ParameterFileContent result;
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized) {
    result.cfg = broadcastParameters::getParametersFromRankZero(paramfile, comm);
  } else {
    boost::property_tree::ini_parser::read_ini(paramfile, result.cfg);
  }
  result.ngroup = result.cfg.get<size_t>("manager.ngroup");
  result.nprocs = result.cfg.get<size_t>("manager.nprocs");
  return result;
}

/**
 * @brief Per-group task assignment from the combination scheme.
 */
struct LocalTasks {
  std::vector<LevelVector> levels;
  std::vector<real> coeffs;
  std::vector<size_t> taskNumbers;
  size_t totalNumTasks;
};

/**
 * @brief Build CombiParameters and distribute tasks from a property tree.
 *
 * Call after MPI groups are set up.  Reads [ct] section, builds the
 * combination scheme (from lmin/lmax or a JSON ctscheme file), distributes
 * tasks across process groups, and returns ready-to-use CombiParameters.
 *
 * Supported ctparam fields (all in [ct] section unless noted):
 *   dim, lmin, lmax, p, ncombi, boundary, basis, chunkSize,
 *   hierarchization_dims, ctscheme, combinationVariant,
 *   forwardDecomposition, reduceCombinationDimsLmin/Lmax,
 *   ncombiLocal, [thirdLevel] host/port
 *
 * @param cfg  Property tree from readParameterFile()
 * @return     {CombiParameters, LocalTasks}
 */
inline std::pair<CombiParameters, LocalTasks> buildCombiParametersFromConfig(
    const boost::property_tree::ptree& cfg) {
  DimType dim = cfg.get<DimType>("ct.dim");

  LevelVector lmin(dim), lmax(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;
  cfg.get<std::string>("ct.lmax") >> lmax;

  std::vector<int> p(dim, 1);
  if (cfg.get_child_optional("ct.p")) {
    cfg.get<std::string>("ct.p") >> p;
  }

  size_t ncombi = cfg.get<size_t>("ct.ncombi");
  uint32_t chunkSize = cfg.get<uint32_t>("ct.chunkSize", 128);
  std::string basis = cfg.get<std::string>("ct.basis", "hat_periodic");

  std::vector<BoundaryType> boundary(dim, 1);
  if (cfg.get_child_optional("ct.boundary")) {
    cfg.get<std::string>("ct.boundary") >> boundary;
  }

  std::vector<bool> hierarchizationDims(dim, true);
  if (cfg.get_child_optional("ct.hierarchization_dims")) {
    cfg.get<std::string>("ct.hierarchization_dims") >> hierarchizationDims;
  }

  // Validate parallelization
  if (cfg.get_child_optional("ct.p")) {
    IndexType checkProcs = 1;
    for (auto k : p) checkProcs *= k;
    size_t nprocs = cfg.get<size_t>("manager.nprocs");
    if (checkProcs != static_cast<IndexType>(nprocs)) {
      throw std::invalid_argument("product of p (" + std::to_string(checkProcs) + ") != nprocs (" +
                                  std::to_string(nprocs) + ")");
    }
  }

  // Build combination scheme and distribute tasks
  LocalTasks tasks;
  std::string ctschemeFile = cfg.get<std::string>("ct.ctscheme", "");
  auto pgroupNumber = theMPISystem()->getProcessGroupNumber();

  if (!ctschemeFile.empty()) {
    CombiMinMaxSchemeFromFile scheme(dim, lmin, lmax, ctschemeFile);
    if (scheme.getProcessGroupNumbers().size() > 0) {
      tasks.totalNumTasks =
          getAssignedLevels(scheme, pgroupNumber, tasks.levels, tasks.coeffs, tasks.taskNumbers);
    } else {
      tasks.totalNumTasks =
          getLoadBalancedLevels(scheme, pgroupNumber, theMPISystem()->getNumGroups(), boundary,
                                tasks.levels, tasks.coeffs, tasks.taskNumbers);
    }
  } else {
    CombiMinMaxScheme scheme(dim, lmin, lmax);
    scheme.createClassicalCombischeme();
    tasks.totalNumTasks =
        getLoadBalancedLevels(scheme, pgroupNumber, theMPISystem()->getNumGroups(), boundary,
                              tasks.levels, tasks.coeffs, tasks.taskNumbers);
  }

  // Combination variant
  std::string variantStr = cfg.get<std::string>("ct.combinationVariant", "subspaceReduce");
  CombinationVariant variant = CombinationVariant::subspaceReduce;
  if (variantStr == "sparseGridReduce")
    variant = CombinationVariant::sparseGridReduce;
  else if (variantStr == "outgroupSparseGridReduce")
    variant = CombinationVariant::outgroupSparseGridReduce;
  else if (variantStr == "chunkedOutgroupSparseGridReduce")
    variant = CombinationVariant::chunkedOutgroupSparseGridReduce;

  // forward decomposition: false at least for periodic setups
  bool forwardDecomposition = cfg.get<bool>("ct.forwardDecomposition", false);

  LevelVector reduceLmin(dim, 0), reduceLmax(dim, 0);
  if (cfg.get_child_optional("ct.reduceCombinationDimsLmin")) {
    cfg.get<std::string>("ct.reduceCombinationDimsLmin") >> reduceLmin;
  }
  if (cfg.get_child_optional("ct.reduceCombinationDimsLmax")) {
    cfg.get<std::string>("ct.reduceCombinationDimsLmax") >> reduceLmax;
  }

  std::string thirdLevelHost = cfg.get<std::string>("thirdLevel.host", "");
  unsigned short thirdLevelPort = cfg.get<unsigned short>("thirdLevel.port", 0);
  size_t ncombiLocal = cfg.get<size_t>("ct.ncombiLocal", ncombi);

  CombiParameters params(dim, lmin, lmax, boundary, tasks.levels, tasks.coeffs, hierarchizationDims,
                         tasks.taskNumbers, ncombiLocal, 1, variant, reduceLmin, reduceLmax,
                         chunkSize, forwardDecomposition, thirdLevelHost, thirdLevelPort, 0);
  setCombiParametersHierarchicalBasesUniform(params, basis);
  params.setParallelization(p);

  // Hierarchization backend (optional, default: discotec)
  std::string backendStr = cfg.get<std::string>("ct.hierarchizationBackend", "discotec");
  if (backendStr == "paliwa") {
    params.setHierarchizationBackend(HierarchizationBackend::PALIWA);
  } else {
    params.setHierarchizationBackend(HierarchizationBackend::DISCOTEC);
  }

  return {std::move(params), std::move(tasks)};
}
}  // namespace combigrid
