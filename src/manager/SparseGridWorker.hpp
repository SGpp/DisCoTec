#pragma once

#include <thread>

#include "combicom/CombiCom.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "manager/TaskWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "mpi/MPIUtils.hpp"
#include "sparsegrid/DistributedSparseGridIO.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
namespace combigrid {

// template <typename CombiDataType, DimType NumDimensions>
class SparseGridWorker {
 public:
  explicit SparseGridWorker(TaskWorker& taskWorkerToReference);

  SparseGridWorker(SparseGridWorker const&) = delete;
  SparseGridWorker& operator=(SparseGridWorker const&) = delete;

  ~SparseGridWorker() = default;

  template <bool keepValuesInExtraSparseGrid = false>
  inline void collectReduceDistribute(CombinationVariant combinationVariant,
                                      uint32_t maxMiBToSendPerThread,
                                      bool collectMinMaxCoefficients = false);

  inline void copyFromExtraDsgToPartialDSG(size_t gridNumber = 0);

  inline void copyFromPartialDsgToExtraDSG(size_t gridNumber = 0);

  /* free DSG memory as intermediate step */
  inline void deleteDsgsData();

  inline void distributeChunkedBroadcasts(uint32_t maxMiBToSendPerThread);

  inline void distributeCombinedSolutionToTasks();

  inline void fillDFGFromDSGU(Task& t, const std::vector<bool>& hierarchizationDims,
                              const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                              const LevelVector& lmin) const;

  inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedUniDSGVector();

  inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedUniDSGVector() const;

  inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getExtraUniDSGVector();

  inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getExtraUniDSGVector() const;

  inline size_t getNumberOfGrids() const;

  inline std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>&
  getSparseGridToUseForThirdLevel(bool thirdLevelExtraSparseGrid);

  inline void initCombinedUniDSGVector(const LevelVector& lmin, LevelVector lmax,
                                       const LevelVector& reduceLmaxByVector, size_t numGrids,
                                       CombinationVariant combinationVariant,
                                       bool clearLevels = false);

  inline void interpolateAndPlotOnLevel(
      const std::string& filename, const LevelVector& levelToEvaluate,
      const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
      const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin,
      const std::vector<int>& parallelization, const std::vector<LevelVector>& decomposition) const;

  inline void maxReduceSubspaceSizesInOutputGroup();

  inline int readDSGsFromDisk(const std::string& filenamePrefixToRead,
                              bool addSpeciesNumberToFileName, bool alwaysReadFullDSG = false);

  inline int readDSGsFromDiskAndReduce(const std::string& filenamePrefixToRead,
                                       uint32_t maxMiBToReadPerThread,
                                       bool addSpeciesNumberToFileName,
                                       bool alwaysReadFullDSG = false);

  inline int readReduce(const std::vector<std::string>& filenamePrefixesToRead,
                        const std::vector<std::string>& startReadingTokenFileNames,
                        uint32_t maxMiBToReadPerThread, bool overwrite);

  inline int reduceExtraSubspaceSizes(const std::vector<std::string>& filenamesToRead,
                                      CombinationVariant combinationVariant, bool overwrite);

  /* reduction within and between process groups */
  inline void reduceLocalAndGlobal(CombinationVariant combinationVariant,
                                   uint32_t maxMiBToSendPerThread,
                                   RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  inline void reduceSubspaceSizesBetweenGroups(CombinationVariant combinationVariant);

  inline void removeReadingFiles(const std::vector<std::string>& filenamePrefixToRead,
                                 const std::vector<std::string>& startReadingTokenFileName,
                                 bool keepSparseGridFiles) const;

  inline void setExtraSparseGrid(bool initializeSizes = true);

  inline void startSingleBroadcastDSGs(CombinationVariant combinationVariant,
                                       RankType broadcastSender, MPI_Request* request);

  inline int writeDSGsToDisk(const std::string& filenamePrefix,
                             CombinationVariant combinationVariant);

  inline int writeExtraSubspaceSizesToFile(const std::string& filenamePrefixToWrite) const;

  inline void writeMinMaxCoefficients(const std::string& fileNamePrefix) const;

  inline void writeTokenFiles(const std::string& writeCompleteTokenFileName) const;

  inline int writeSubspaceSizesToFile(const std::string& filenamePrefixToWrite) const;

  inline void zeroDsgsData(CombinationVariant combinationVariant);

 private:
  TaskWorker& taskWorkerRef_;

  /**
   * Vector containing the combined distributed sparse grids
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> combinedUniDSGVector_;

  /**
   * Vector containing the third level extra distributed sparse grids (conjoint grids)
   */
  std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>> extraUniDSGVector_;

  /**
   * @brief copy the sparse grid data into the full grid and dehierarchize
   *
   * @param dfg the distributed full grid to fill
   * @param g the dimension index (in the case that there are multiple different full grids per
   * task)
   */
  inline void fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg, size_t g,
                              const std::vector<BoundaryType>& boundary,
                              const std::vector<bool>& hierarchizationDims,
                              const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                              const LevelVector& lmin) const;
};

inline SparseGridWorker::SparseGridWorker(TaskWorker& taskWorkerToReference)
    : taskWorkerRef_(taskWorkerToReference) {}

template <bool keepValuesInExtraSparseGrid>
inline void SparseGridWorker::collectReduceDistribute(CombinationVariant combinationVariant,
                                                      uint32_t maxMiBToSendPerThread,
                                                      bool collectMinMaxCoefficients) {
  auto numGrids = this->getNumberOfGrids();
  assert(combinationVariant == CombinationVariant::chunkedOutgroupSparseGridReduce);
  assert(numGrids == 1 && "Initialize dsgu first with initCombinedUniDSGVector()");
  assert(this->getCombinedUniDSGVector()[0]->getSubspacesByCommunicator().size() < 2 &&
         "Initialize dsgu's outgroup communicator");

  // allow only up to the specified MiB per reduction
  auto chunkSize =
      combigrid::CombiCom::getGlobalReduceChunkSize<CombiDataType>(maxMiBToSendPerThread);

  for (size_t g = 0; g < numGrids; ++g) {
    auto& dsg = this->getCombinedUniDSGVector()[g];
    if (!dsg->getSubspacesByCommunicator().empty()) {
      auto& chunkedSubspaces = combigrid::CombiCom::getChunkedSubspaces(
          *dsg, dsg->getSubspacesByCommunicator()[0].second, chunkSize);
      for (auto& subspaceChunk : chunkedSubspaces) {
        // allocate new subspace vector
        dsg->allocateDifferentSubspaces(std::move(subspaceChunk));
        assert(dsg->getRawDataSize() <= static_cast<size_t>(chunkSize));
#ifndef NDEBUG
        auto myRawDataSize = dsg->getRawDataSize();
        decltype(myRawDataSize) maxRawDataSize = 0;
        MPI_Allreduce(&myRawDataSize, &maxRawDataSize, 1,
                      getMPIDatatype(abstraction::getabstractionDataType<size_t>()), MPI_MAX,
                      theMPISystem()->getGlobalReduceComm());
        if (myRawDataSize != maxRawDataSize) {
          throw std::runtime_error(
              "collectReduceDistribute: Raw data size is not the same across all ranks; my size: " +
              std::to_string(myRawDataSize) + ", max size: " + std::to_string(maxRawDataSize));
        }
#endif

        // local reduce (fg -> sg, within rank)
        for (const auto& t : this->taskWorkerRef_.getTasks()) {
          const DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(g);
          dsg->addDistributedFullGrid<false>(dfg, t->getCoefficient());
        }
        // global reduce (across process groups)
        CombiCom::distributedGlobalSubspaceReduce<DistributedSparseGridUniform<CombiDataType>,
                                                  true>(*dsg, maxMiBToSendPerThread);
        // assert(CombiCom::sumAndCheckSubspaceSizes(*dsg)); // todo adapt for allocated spaces
        if constexpr (keepValuesInExtraSparseGrid) {
          this->copyFromPartialDsgToExtraDSG(g);
        }
        if (collectMinMaxCoefficients) {
          dsg->accumulateMinMaxCoefficients();
        }

        // distribute (sg -> fg, within rank)
        for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
          // fill dfg with hierarchical coefficients from distributed sparse grid
          taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG<false>(*dsg);
        }
      }
    }
    // update our ingroup subspaces, ie the ones that are not communicated
    // (given by those that have a data size set but no communicator)
    auto& chunkedSubspaces =
        combigrid::CombiCom::getChunkedSubspaces(*dsg, dsg->getIngroupSubspaces(), chunkSize);
    for (auto& subspaceChunk : chunkedSubspaces) {
      // allocate new subspace vector
      dsg->allocateDifferentSubspaces(std::move(subspaceChunk));

      // local reduce (fg -> sg, within rank)
      for (const auto& t : this->taskWorkerRef_.getTasks()) {
        const DistributedFullGrid<CombiDataType>& dfg =
            t->getDistributedFullGrid(static_cast<int>(g));
        dsg->addDistributedFullGrid<false>(dfg, t->getCoefficient());
      }
      if constexpr (keepValuesInExtraSparseGrid) {
        this->copyFromPartialDsgToExtraDSG(g);
      }
      if (collectMinMaxCoefficients) {
        dsg->accumulateMinMaxCoefficients();
      }

      // distribute (sg -> fg, within rank)
      for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
        // fill dfg with hierarchical coefficients from distributed sparse grid
        taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG<false>(*dsg);
      }
    }
    dsg->deleteSubspaceData();
  }
}

inline void SparseGridWorker::copyFromPartialDsgToExtraDSG(size_t gridNumber) {
  assert(gridNumber == 0);
  assert(this->getCombinedUniDSGVector().size() == 1);
  assert(this->getExtraUniDSGVector().size() == 1);
  if (this->getExtraUniDSGVector()[0]->getRawDataSize() == 0) {
    throw std::runtime_error("Extra DSG is not initialized");
  }
  auto& dsg = this->getCombinedUniDSGVector()[gridNumber];
  auto& extraDSG = this->getExtraUniDSGVector()[gridNumber];
  auto subspacesToCopy = std::vector(dsg->getCurrentlyAllocatedSubspaces().cbegin(),
                                     dsg->getCurrentlyAllocatedSubspaces().cend());
  subspacesToCopy.erase(
      std::remove_if(subspacesToCopy.begin(), subspacesToCopy.end(),
                     [&extraDSG](const auto& s) { return extraDSG->getDataSize(s) == 0; }),
      subspacesToCopy.end());
  for ([[maybe_unused]] const auto& s : subspacesToCopy) {
    assert(extraDSG->getDataSize(s) == extraDSG->getAllocatedDataSize(s));
    assert(dsg->getDataSize(s) == dsg->getAllocatedDataSize(s));
  }
  extraDSG->copyDataFrom(*dsg, subspacesToCopy);
}

inline void SparseGridWorker::copyFromExtraDsgToPartialDSG(size_t gridNumber) {
  assert(gridNumber == 0);
  assert(this->getCombinedUniDSGVector().size() == 1);
  assert(this->getExtraUniDSGVector().size() == 1);
  auto& dsg = this->getCombinedUniDSGVector()[gridNumber];
  auto& extraDSG = this->getExtraUniDSGVector()[gridNumber];
  auto subspacesToCopy = std::vector(dsg->getCurrentlyAllocatedSubspaces().cbegin(),
                                     dsg->getCurrentlyAllocatedSubspaces().cend());
  subspacesToCopy.erase(
      std::remove_if(subspacesToCopy.begin(), subspacesToCopy.end(),
                     [&extraDSG](const auto& s) { return extraDSG->getDataSize(s) == 0; }),
      subspacesToCopy.end());
  for ([[maybe_unused]] const auto& s : subspacesToCopy) {
    assert(extraDSG->getDataSize(s) == extraDSG->getAllocatedDataSize(s));
    assert(dsg->getDataSize(s) == dsg->getAllocatedDataSize(s));
  }
  dsg->copyDataFrom(*extraDSG, subspacesToCopy);
}

inline void SparseGridWorker::deleteDsgsData() {
  for (auto& dsg : this->getCombinedUniDSGVector()) dsg->deleteSubspaceData();
  for (auto& dsg : this->getExtraUniDSGVector()) dsg->deleteSubspaceData();
}

inline void SparseGridWorker::distributeChunkedBroadcasts(uint32_t maxMiBToSendPerThread) {
  RankType broadcastSender = theMPISystem()->getOutputRankInGlobalReduceComm();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();

  // communicate the subspaces that will be broadcast
  std::vector<typename AnyDistributedSparseGrid::SubspaceIndexType> subspaces;
  assert(this->getCombinedUniDSGVector().size() == 1);
  auto& uniDSG = this->getCombinedUniDSGVector()[0];
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    assert(broadcastSender == theMPISystem()->getGlobalReduceRank());
    auto& extraDSG = this->getExtraUniDSGVector()[0];
    subspaces.reserve(extraDSG->getNumSubspaces());
    std::vector<typename AnyDistributedSparseGrid::SubspaceIndexType> emptySubspaces = {};
    const auto* pointerToCandidateSubspaces = &emptySubspaces;
    if (!uniDSG->getSubspacesByCommunicator().empty()) {
      pointerToCandidateSubspaces = &uniDSG->getSubspacesByCommunicator()[0].second;
      assert(uniDSG->getSubspacesByCommunicator().size() == 1);
    }
    for (typename AnyDistributedSparseGrid::SubspaceIndexType i = 0;
         i < extraDSG->getNumSubspaces(); ++i) {
      if (extraDSG->getDataSize(i) > 0 &&  // TODO find something smarter than this
          std::find(pointerToCandidateSubspaces->cbegin(), pointerToCandidateSubspaces->cend(),
                    i) != pointerToCandidateSubspaces->cend()) {
        subspaces.push_back(i);
      }
    }
  }
  else {
    assert(this->getExtraUniDSGVector().empty());
  }
  MPIUtils::broadcastContainer(subspaces, broadcastSender, globalReduceComm);

  auto chunkSize =
      combigrid::CombiCom::getGlobalReduceChunkSize<CombiDataType>(maxMiBToSendPerThread);
  assert(this->getCombinedUniDSGVector().size() == 1);
  for (auto& dsg : this->getCombinedUniDSGVector()) {
    // make sure data sizes are set
    for ([[maybe_unused]] auto& subspace : subspaces) {
      assert(dsg->getDataSize(subspace) > 0);
    }

    auto& chunkedSubspaces = combigrid::CombiCom::getChunkedSubspaces(*dsg, subspaces, chunkSize);
    size_t roundNumber = 0;
    for (auto& subspaceChunk : chunkedSubspaces) {
      dsg->allocateDifferentSubspaces(std::move(subspaceChunk));

      OUTPUT_GROUP_EXCLUSIVE_SECTION {
        // I need to initialize my dsg from the extra dsg
        this->copyFromExtraDsgToPartialDSG(0);
      }

      MPI_Request request;
      CombiCom::asyncBcastOutgroupDsgData<decltype(*dsg), true>(*dsg, broadcastSender,
                                                                globalReduceComm, &request);

      OUTPUT_GROUP_EXCLUSIVE_SECTION {}
      else {
        // non-output ranks need to wait before they can extract
        if (roundNumber == 0) Stats::startEvent("wait 1st bcast");
        [[maybe_unused]] auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
        if (roundNumber == 0) Stats::stopEvent("wait 1st bcast");
        assert(returnedValue == MPI_SUCCESS);
      }
      // distribute (sg -> fg, within rank)
      for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
        // fill dfg with hierarchical coefficients from distributed sparse grid
        taskToUpdate->getDistributedFullGrid(0).extractFromUniformSG<false>(*dsg);
      }
      OUTPUT_GROUP_EXCLUSIVE_SECTION {
        // output ranks can wait later
        [[maybe_unused]] auto returnedValue = MPI_Wait(&request, MPI_STATUS_IGNORE);
        assert(returnedValue == MPI_SUCCESS);
      }
      ++roundNumber;
    }
  }
  subspaces.clear();

  // and the output group may also need to update the subspaces from extra dsg that it has alone
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    auto& extraDSG = this->getExtraUniDSGVector()[0];
    const auto& ingroupSubpspaces = uniDSG->getIngroupSubspaces();
    for (typename AnyDistributedSparseGrid::SubspaceIndexType i = 0; i < uniDSG->getNumSubspaces();
         ++i) {
      if (ingroupSubpspaces.find(i) != ingroupSubpspaces.end() && extraDSG->getDataSize(i) > 0) {
        subspaces.push_back(i);
      }
    }
    for (auto& dsg : this->getCombinedUniDSGVector()) {
      // make sure data sizes are set
      for ([[maybe_unused]] auto& subspace : subspaces) {
        assert(dsg->getDataSize(subspace) > 0);
      }

      auto& chunkedSubspaces = combigrid::CombiCom::getChunkedSubspaces(*dsg, subspaces, chunkSize);
      size_t roundNumber = 0;
      for (auto& subspaceChunk : chunkedSubspaces) {
        dsg->allocateDifferentSubspaces(std::move(subspaceChunk));
        this->copyFromExtraDsgToPartialDSG(0);
        // distribute (sg -> fg, within rank)
        for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
          taskToUpdate->getDistributedFullGrid(0).extractFromUniformSG<false>(*dsg);
        }
        ++roundNumber;
      }
    }
  }
}

inline void SparseGridWorker::distributeCombinedSolutionToTasks() {
  // not better than the parallel loop within extractFromUniformSG
  // #pragma omp parallel for collapse(2) default(none) schedule(static)
  for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
    for (size_t g = 0; g < this->getNumberOfGrids(); ++g) {
      // fill dfg with hierarchical coefficients from distributed sparse grid
      taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG(
          *this->getCombinedUniDSGVector()[g]);
    }
  }
}

inline void SparseGridWorker::fillDFGFromDSGU(
    DistributedFullGrid<CombiDataType>& dfg, size_t g, const std::vector<BoundaryType>& boundary,
    const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) const {
  // fill dfg with hierarchical coefficients from distributed sparse grid
  dfg.extractFromUniformSG(*this->getCombinedUniDSGVector()[g]);

  bool anyNotBoundary =
      std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });

  if (anyNotBoundary) {
    std::remove_reference_t<decltype(lmin)> zeroLMin(lmin.size(), 0);
    DistributedHierarchization::dehierarchizeDFG(dfg, hierarchizationDims, hierarchicalBases,
                                                 zeroLMin);
  } else {
    DistributedHierarchization::dehierarchizeDFG(dfg, hierarchizationDims, hierarchicalBases, lmin);
  }
}

inline void SparseGridWorker::fillDFGFromDSGU(
    Task& t, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) const {
  for (size_t g = 0; g < this->getNumberOfGrids(); g++) {
    assert(this->getCombinedUniDSGVector()[g] != nullptr);
    this->fillDFGFromDSGU(t.getDistributedFullGrid(g), g, t.getBoundary(), hierarchizationDims,
                          hierarchicalBases, lmin);
  }
}

inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
SparseGridWorker::getCombinedUniDSGVector() {
  return combinedUniDSGVector_;
}

inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
SparseGridWorker::getCombinedUniDSGVector() const {
  return combinedUniDSGVector_;
}

inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
SparseGridWorker::getExtraUniDSGVector() {
  return extraUniDSGVector_;
}

inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
SparseGridWorker::getExtraUniDSGVector() const {
  return extraUniDSGVector_;
}

inline size_t SparseGridWorker::getNumberOfGrids() const {
  return this->combinedUniDSGVector_.size();
}

inline std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>&
SparseGridWorker::getSparseGridToUseForThirdLevel(bool thirdLevelExtraSparseGrid) {
  if (this->getNumberOfGrids() != 1) {
    throw std::runtime_error("number of grids > 1 -- not implemented on pg manager's side");
  }
  if (thirdLevelExtraSparseGrid) {
    this->setExtraSparseGrid(false);  // don't initialize, would be overwritten
    return this->getExtraUniDSGVector()[0];
  }
  return this->getCombinedUniDSGVector()[0];
}

/** Initializes the dsgu for each species by setting the subspace sizes of all
 * dfgs in the global reduce comm. After calling, all workers which share the
 * same spatial distribution of the dsgu (those who combine during global
 * reduce) have the same sized subspaces and thus share the same sized dsg.
 *
 * Attention: No data is created here, only subspace sizes are shared.
 */
inline void SparseGridWorker::initCombinedUniDSGVector(const LevelVector& lmin, LevelVector lmax,
                                                       const LevelVector& reduceLmaxByVector,
                                                       size_t numGrids,
                                                       CombinationVariant combinationVariant,
                                                       bool clearLevels) {
  if (this->taskWorkerRef_.getTasks().size() == 0) {
    std::cout << "Possible error: task size is 0! \n";
  }
  // the dsg can be smaller than lmax because the highest subspaces do not have
  // to be exchanged: then reduceLmaxByVector is > 0

  for (size_t i = 0; i < reduceLmaxByVector.size(); ++i) {
    lmax[i] = std::max(lmin[i], static_cast<LevelType>(lmax[i] - reduceLmaxByVector[i]));
  }

  // get all subspaces in the (optimized) combischeme, create dsgs
  combinedUniDSGVector_.reserve(static_cast<size_t>(numGrids));
  for (size_t g = 0; g < numGrids; ++g) {
    combinedUniDSGVector_.emplace_back(std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(static_cast<DimType>(lmax.size()), lmax,
                                                        lmin, theMPISystem()->getLocalComm())));
  }
  for (size_t g = 0; g < combinedUniDSGVector_.size(); ++g) {
    for (const auto& t : this->taskWorkerRef_.getTasks()) {
      DistributedFullGrid<CombiDataType>& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      // set subspace sizes locally
      combinedUniDSGVector_[g]->registerDistributedFullGrid(dfg);
    }
    // we may clear the levels_ member of the sparse grids here to save memory
    // but only if we need no new full grids initialized from the sparse grids!
    // ...such as for rescheduling or interpolation (parallelEval/ evalNorm / ...)
    if (clearLevels) {
      combinedUniDSGVector_[(size_t)g]->resetLevels();
    }

    // create the kahan buffer now, so it has only the subspaces present on the grids in this
    // process group
    if (combinationVariant != CombinationVariant::chunkedOutgroupSparseGridReduce) {
      combinedUniDSGVector_[g]->createKahanBuffer();
    }
  }

  this->reduceSubspaceSizesBetweenGroups(combinationVariant);
}

// cf https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
static bool endsWith(const std::string& str, const std::string& suffix) {
  return str.size() >= suffix.size() &&
         0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

inline void SparseGridWorker::interpolateAndPlotOnLevel(
    const std::string& filename, const LevelVector& levelToEvaluate,
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin,
    const std::vector<int>& parallelization, const std::vector<IndexVector>& decomposition) const {
  assert(levelToEvaluate.size() == parallelization.size());
  // create dfg
  OwningDistributedFullGrid<CombiDataType> dfg(static_cast<DimType>(levelToEvaluate.size()),
                                               levelToEvaluate, theMPISystem()->getLocalComm(),
                                               boundary, parallelization, false, decomposition);
  for (size_t g = 0; g < this->getNumberOfGrids(); ++g) {  // loop over all grids and plot them
    this->fillDFGFromDSGU(dfg, g, boundary, hierarchizationDims, hierarchicalBases, lmin);
    // save dfg to file with MPI-IO
    if (endsWith(filename, ".vtk")) {
      dfg.writePlotFileVTK(filename.c_str());
    } else {
      std::string fn = filename;
      auto pos = fn.find(".");
      if (pos != std::string::npos) {
        // if filename contains ".", insert grid number before that
        fn.insert(pos, "_" + std::to_string(g));
      }
      dfg.writePlotFile(fn.c_str());
    }
  }
}

inline void SparseGridWorker::maxReduceSubspaceSizesInOutputGroup() {
  RankType globalReduceRankThatCollects = theMPISystem()->getOutputRankInGlobalReduceComm();
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    auto dsgToUse = this->getCombinedUniDSGVector()[0].get();
    if (this->getExtraUniDSGVector().empty() == false) {
      dsgToUse = this->getExtraUniDSGVector()[0].get();
    }
    CombiCom::maxReduceSubspaceSizesAcrossGroups(*dsgToUse, globalReduceRankThatCollects,
                                                 globalReduceComm);
  }
  else {
    assert(this->getExtraUniDSGVector().empty());
    CombiCom::maxReduceSubspaceSizesAcrossGroups(*this->getCombinedUniDSGVector()[0],
                                                 globalReduceRankThatCollects, globalReduceComm);
  }
}

inline int SparseGridWorker::readDSGsFromDisk(const std::string& filenamePrefixToRead,
                                              bool addSpeciesNumberToFileName,
                                              bool alwaysReadFullDSG) {
  int numRead = 0;
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    auto filename = filenamePrefixToRead;
    if (addSpeciesNumberToFileName) {
      filename += "_" + std::to_string(i);
    }
    if (this->getExtraUniDSGVector().size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = this->getExtraUniDSGVector()[i].get();
      numRead += DistributedSparseGridIO::readSomeFiles(*dsgToUse, filename);
    } else {
      numRead += DistributedSparseGridIO::readOneFile(*dsgToUse, filename);
    }
    if (this->getExtraUniDSGVector().size() > 0 && uniDsg->isSubspaceDataCreated()) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }  // else: assert this only happens for chunkedOutgroupSparseGridReduce //TODO
  }
  return numRead;
}

inline int SparseGridWorker::readDSGsFromDiskAndReduce(const std::string& filenamePrefixToRead,
                                                       uint32_t maxMiBToReadPerThread,
                                                       bool addSpeciesNumberToFileName,
                                                       bool alwaysReadFullDSG) {
  int numReduced = 0;
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto filename = filenamePrefixToRead;
    if (addSpeciesNumberToFileName) {
      filename += "_" + std::to_string(i);
    }
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getExtraUniDSGVector().size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = this->getExtraUniDSGVector()[i].get();
      numReduced += DistributedSparseGridIO::readSomeFilesAndReduce(*dsgToUse, filename,
                                                                    maxMiBToReadPerThread);
    } else {
      numReduced +=
          DistributedSparseGridIO::readOneFileAndReduce(*dsgToUse, filename, maxMiBToReadPerThread);
    }
    if (this->getExtraUniDSGVector().size() > 0 && uniDsg->isSubspaceDataCreated()) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }
  }
  return numReduced;
}

inline int SparseGridWorker::readReduce(const std::vector<std::string>& filenamePrefixesToRead,
                                        const std::vector<std::string>& startReadingTokenFileNames,
                                        uint32_t maxMiBToReadPerThread, bool overwrite) {
  int numRead = 0;
  assert(this->getNumberOfGrids() == 1);  // for more species, todo loop

  // wait until we can start to read any of the files
  std::set<size_t> indicesStillToReadReduce;
  for (size_t i = 0; i < filenamePrefixesToRead.size(); ++i) {
    indicesStillToReadReduce.insert(i);
  }
  auto filePart = theMPISystem()->getFilePartNumber();
  bool firstTry = true;
  while (!indicesStillToReadReduce.empty()) {
    std::string filePartTokenToRead;
    std::string filePartNameToRead;
    for (auto it = indicesStillToReadReduce.begin(); it != indicesStillToReadReduce.end(); ++it) {
      if (combigrid::theMPISystem()->getOutputComm() != MPI_COMM_NULL) {
        filePartTokenToRead = startReadingTokenFileNames[*it] + ".part" + std::to_string(filePart);
      } else {
        filePartTokenToRead = startReadingTokenFileNames[*it];
      }
      if (firstTry) {
        std::cout << "trying to read " << filePartTokenToRead << std::endl;
        firstTry = false;
      }
      if (combigrid::getFileExistsRootOnly(filePartTokenToRead, theMPISystem()->getOutputComm())) {
        if (combigrid::theMPISystem()->getOutputComm() != MPI_COMM_NULL) {
          filePartNameToRead =
              filenamePrefixesToRead[*it] + "_0" + ".part" + std::to_string(filePart);
        } else {
          filePartNameToRead = filenamePrefixesToRead[*it] + "_0";
        }
        std::cout << " read from " << filePartNameToRead << std::endl;
        overwrite ? Stats::startEvent("read SG") : Stats::startEvent("read/reduce SG");
        if (overwrite) {
          numRead += this->readDSGsFromDisk(filePartNameToRead, false);
        } else {
          numRead +=
              this->readDSGsFromDiskAndReduce(filePartNameToRead, maxMiBToReadPerThread, false);
        }
        assert(numRead > 0);
        overwrite ? Stats::stopEvent("read SG") : Stats::stopEvent("read/reduce SG");
        indicesStillToReadReduce.erase(it);
        break;  // because iterators may be invalidated
      }
    }
    // wait for 200ms
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }
  if (this->getNumberOfGrids() != 1) {
    throw std::runtime_error("Combining more than one DSG is not implemented yet");
  }
  return numRead;
}

inline int SparseGridWorker::reduceExtraSubspaceSizes(
    const std::vector<std::string>& filenamesToRead, CombinationVariant combinationVariant,
    bool overwrite) {
  int numSubspacesReduced = 0;
  if (!overwrite) {
    // make output group's dsg "larger" to be able to reduce
    this->maxReduceSubspaceSizesInOutputGroup();
  }
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    if (this->getExtraUniDSGVector().empty()) {
      this->setExtraSparseGrid(true);
    }
  }
  OUTPUT_GROUP_EXCLUSIVE_SECTION {
    auto dsgToUse = this->getExtraUniDSGVector()[0].get();
#ifndef NDEBUG
    // duplicate subspace sizes to validate later
    std::vector<SubspaceSizeType> subspaceSizesToValidate = dsgToUse->getSubspaceDataSizes();
#endif
    if (overwrite) {
      assert(filenamesToRead.size() == 1);
      numSubspacesReduced =
          DistributedSparseGridIO::readSubspaceSizesFromFile(*dsgToUse, filenamesToRead[0], false);
#ifndef NDEBUG
      for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
        if (!(dsgToUse->getSubspaceDataSizes()[i] == 0 || subspaceSizesToValidate[i] == 0 ||
              dsgToUse->getSubspaceDataSizes()[i] == subspaceSizesToValidate[i])) {
          std::cout << "i " << i << " dsgToUse->getSubspaceDataSizes()[i] "
                    << dsgToUse->getSubspaceDataSizes()[i] << " subspaceSizesToValidate[i] "
                    << subspaceSizesToValidate[i] << std::endl;
        }
        assert(dsgToUse->getSubspaceDataSizes()[i] == 0 || subspaceSizesToValidate[i] == 0 ||
               dsgToUse->getSubspaceDataSizes()[i] == subspaceSizesToValidate[i]);
      }
#endif
    } else {
      if (filenamesToRead.size() == 1) {
        auto minFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
          return std::min(a, b);
        };
        numSubspacesReduced = DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
            *dsgToUse, filenamesToRead[0], minFunctionInstantiation, 0, false);
#ifndef NDEBUG
        for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
          if (!(dsgToUse->getSubspaceDataSizes()[i] == 0 ||
                dsgToUse->getSubspaceDataSizes()[i] == subspaceSizesToValidate[i])) {
            std::cout << "i " << i << " dsgToUse->getSubspaceDataSizes()[i] "
                      << dsgToUse->getSubspaceDataSizes()[i] << " subspaceSizesToValidate[i] "
                      << subspaceSizesToValidate[i] << std::endl;
          }
          assert(dsgToUse->getSubspaceDataSizes()[i] == 0 ||
                 dsgToUse->getSubspaceDataSizes()[i] == subspaceSizesToValidate[i]);
        }
#endif
      } else if (filenamesToRead.size() == 2) {
        auto moreThanOneFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b,
                                                   SubspaceSizeType c) {
          auto maxValue = std::max(a, std::max(b, c));
          return ((a + b + c) > maxValue ? maxValue : 0);
        };
        numSubspacesReduced = DistributedSparseGridIO::readReduceSubspaceSizesFromFiles(
            *dsgToUse, filenamesToRead, moreThanOneFunctionInstantiation, 0, false);
      } else {
        throw std::runtime_error("More than 3 systems not implemented yet");
      }
    }
#ifndef NDEBUG
    assert(subspaceSizesToValidate.size() == dsgToUse->getSubspaceDataSizes().size());
    auto numDOFtoValidate =
        std::accumulate(subspaceSizesToValidate.begin(), subspaceSizesToValidate.end(), 0);
    auto numDOFnow = std::accumulate(dsgToUse->getSubspaceDataSizes().begin(),
                                     dsgToUse->getSubspaceDataSizes().end(), 0);
    assert(numDOFtoValidate >= numDOFnow || filenamesToRead.size() > 1);
#endif
    // may need to re-size original spaces if original sparse grid was too small
    this->getCombinedUniDSGVector()[0]->maxReduceSubspaceSizes(*dsgToUse);
  }
  this->reduceSubspaceSizesBetweenGroups(combinationVariant);
  return numSubspacesReduced;
}

inline void SparseGridWorker::reduceLocalAndGlobal(CombinationVariant combinationVariant,
                                                   uint32_t maxMiBToSendPerThread,
                                                   RankType globalReduceRankThatCollects) {
  assert(this->getNumberOfGrids() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  auto numGrids = this->getNumberOfGrids();
  this->zeroDsgsData(combinationVariant);
  // local reduce (within rank)
  for (const auto& t : this->taskWorkerRef_.getTasks()) {
    for (size_t g = 0; g < numGrids; ++g) {
      const DistributedFullGrid<CombiDataType>& dfg =
          t->getDistributedFullGrid(static_cast<int>(g));
      this->getCombinedUniDSGVector()[g]->addDistributedFullGrid(dfg, t->getCoefficient());
    }
  }
  // global reduce (across process groups)
  for (size_t g = 0; g < numGrids; ++g) {
    if (combinationVariant == CombinationVariant::sparseGridReduce) {
      CombiCom::distributedGlobalSparseGridReduce(
          *this->getCombinedUniDSGVector()[g], maxMiBToSendPerThread, globalReduceRankThatCollects);
    } else if (combinationVariant == CombinationVariant::subspaceReduce ||
               combinationVariant == CombinationVariant::outgroupSparseGridReduce) {
      CombiCom::distributedGlobalSubspaceReduce(
          *this->getCombinedUniDSGVector()[g], maxMiBToSendPerThread, globalReduceRankThatCollects);
    } else {
      throw std::runtime_error("Combination variant not implemented");
    }
    assert(CombiCom::sumAndCheckSubspaceSizes(*this->getCombinedUniDSGVector()[g]));
    if (std::isnan(this->getCombinedUniDSGVector()[g]->getData(0)[0])) {
      throw std::runtime_error("aborting due to nans in combined sparse grid");
    }
  }
}

inline void SparseGridWorker::reduceSubspaceSizesBetweenGroups(
    CombinationVariant combinationVariant) {
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  if (combinationVariant == CombinationVariant::sparseGridReduce) {
    for (auto& uniDSG : combinedUniDSGVector_) {
      CombiCom::reduceSubspaceSizes(*uniDSG, globalReduceComm);
    }
  } else if (combinationVariant == CombinationVariant::subspaceReduce) {
    for (auto& uniDSG : combinedUniDSGVector_) {
      uniDSG->setSubspaceCommunicators(globalReduceComm, theMPISystem()->getGlobalReduceRank());
      uniDSG->deleteSubspaceData();
    }
  } else if (combinationVariant == CombinationVariant::outgroupSparseGridReduce ||
             combinationVariant == CombinationVariant::chunkedOutgroupSparseGridReduce) {
    for (auto& uniDSG : combinedUniDSGVector_) {
      uniDSG->setOutgroupCommunicator(globalReduceComm, theMPISystem()->getGlobalReduceRank());
      uniDSG->deleteSubspaceData();
    }
    assert(this->getCombinedUniDSGVector()[0]->getSubspacesByCommunicator().size() < 2);
  } else {
    throw std::runtime_error("Combination variant not implemented");
  }
}

void SparseGridWorker::removeReadingFiles(
    const std::vector<std::string>& filenamePrefixesToRead,
    const std::vector<std::string>& startReadingTokenFileNames, bool keepSparseGridFiles) const {
  assert(filenamePrefixesToRead.size() == startReadingTokenFileNames.size());
  OUTPUT_ROOT_EXCLUSIVE_SECTION {
    for (size_t i = 0; i < filenamePrefixesToRead.size(); ++i) {
      auto filePart = theMPISystem()->getFilePartNumber();
      // remove reading token
      auto partTokenString = startReadingTokenFileNames[i];
      if (combigrid::theMPISystem()->getOutputComm() != MPI_COMM_NULL) {
        partTokenString += ".part" + std::to_string(filePart);
      }
      std::filesystem::remove(partTokenString);
      // remove sparse grid file(s)
      if (!keepSparseGridFiles) {
        std::filesystem::remove(filenamePrefixesToRead[i] + "_" + std::to_string(0));
        for (const auto& entry : std::filesystem::directory_iterator(".")) {
          if (entry.path().string().find(filenamePrefixesToRead[i] + "_" + std::to_string(0) +
                                         ".part" + std::to_string(filePart)) != std::string::npos) {
            assert(entry.is_regular_file());
            std::filesystem::remove(entry.path());
          }
        }
      }
    }
  }
}

inline void SparseGridWorker::setExtraSparseGrid(bool initializeSizes) {
  if (this->getNumberOfGrids() != 1) {
    throw std::runtime_error("this->getCombinedUniDSGVector() is empty");
  }
  if (!this->getExtraUniDSGVector().empty()) {
    throw std::runtime_error(
        "this->getExtraUniDSGVector() is not empty-- if you think this is ok, try to remove the "
        "if-else "
        "here");
  }

  // create new vector for extra sparse grids (that will be only on this process group)
  this->getExtraUniDSGVector().resize(this->getNumberOfGrids());
  for (auto& extraUniDSG : this->getExtraUniDSGVector()) {
    extraUniDSG = std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(
            this->getCombinedUniDSGVector()[0]->getDim(),
            this->getCombinedUniDSGVector()[0]->getAllLevelVectors(),
            theMPISystem()->getOutputGroupComm()));
    // create Kahan buffer now (at zero size), because summation is not needed on this sparse grid
    extraUniDSG->createKahanBuffer();
    if (initializeSizes) {
      for (typename AnyDistributedSparseGrid::SubspaceIndexType i = 0;
           i < extraUniDSG->getNumSubspaces(); ++i) {
        extraUniDSG->setDataSize(i, this->getCombinedUniDSGVector()[0]->getDataSize(i));
      }
    }
    // level vectors are not required; read from the initial sparse grid if needed
    extraUniDSG->resetLevels();
  }
}

inline void SparseGridWorker::startSingleBroadcastDSGs(CombinationVariant combinationVariant,
                                                       RankType broadcastSender,
                                                       MPI_Request* request) {
  if (this->getNumberOfGrids() != 1) {
    throw std::runtime_error("Combining more than one DSG is not implemented yet");
  }
  if (combinationVariant == CombinationVariant::sparseGridReduce) {
    // distribute solution in globalReduceComm to other pgs
    return CombiCom::asyncBcastDsgData(*this->getCombinedUniDSGVector()[0], broadcastSender,
                                       theMPISystem()->getGlobalReduceComm(), request);
  } else if (combinationVariant == CombinationVariant::outgroupSparseGridReduce) {
    return CombiCom::asyncBcastOutgroupDsgData(*this->getCombinedUniDSGVector()[0], broadcastSender,
                                               theMPISystem()->getGlobalReduceComm(), request);
  } else {
    throw std::runtime_error("Combination variant not implemented");
  }
}

inline int SparseGridWorker::writeDSGsToDisk(const std::string& filenamePrefix,
                                             CombinationVariant combinationVariant) {
  int numWritten = 0;
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto filename = filenamePrefix + "_" + std::to_string(i);
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getExtraUniDSGVector().size() > 0) {
      assert(theMPISystem()->getOutputComm() != MPI_COMM_NULL);
      dsgToUse = this->getExtraUniDSGVector()[i].get();
      if (combinationVariant == CombinationVariant::sparseGridReduce ||
          combinationVariant == CombinationVariant::outgroupSparseGridReduce) {
        dsgToUse->copyDataFrom(*uniDsg);
      }
      assert(dsgToUse->isSubspaceDataCreated());
      numWritten += DistributedSparseGridIO::writeSomeFiles(*dsgToUse, filename);

    } else {
      assert(dsgToUse->isSubspaceDataCreated());
      numWritten += DistributedSparseGridIO::writeOneFile(*dsgToUse, filename);
    }
  }
  return numWritten;
}

inline int SparseGridWorker::writeExtraSubspaceSizesToFile(
    const std::string& filenamePrefixToWrite) const {
  assert(this->getExtraUniDSGVector().size() == 1);
  return DistributedSparseGridIO::writeSubspaceSizesToFile(*this->getExtraUniDSGVector()[0],
                                                           filenamePrefixToWrite);
}

inline void SparseGridWorker::writeMinMaxCoefficients(const std::string& fileNamePrefix) const {
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    DistributedSparseGridIO::writeMinMaxCoefficents(
        *(this->getCombinedUniDSGVector()[i]), fileNamePrefix, i, theMPISystem()->getLocalComm(),
        theMPISystem()->getGlobalReduceComm());
  }
}

inline void SparseGridWorker::writeTokenFiles(const std::string& writeCompleteTokenFileName) const {
  OUTPUT_ROOT_EXCLUSIVE_SECTION {
    if (combigrid::theMPISystem()->getOutputComm() == MPI_COMM_NULL) {
      // write single token
      std::ofstream tokenFile(writeCompleteTokenFileName);
    } else {
      auto filePart = theMPISystem()->getFilePartNumber();
      std::ofstream tokenFile(writeCompleteTokenFileName + ".part" + std::to_string(filePart));
    }
  }
}

inline int SparseGridWorker::writeSubspaceSizesToFile(
    const std::string& filenamePrefixToWrite) const {
  return DistributedSparseGridIO::writeSubspaceSizesToFile(*this->getCombinedUniDSGVector()[0],
                                                           filenamePrefixToWrite);
}

inline void SparseGridWorker::zeroDsgsData(CombinationVariant combinationVariant) {
  for (auto& dsg : this->getCombinedUniDSGVector()) {
    if (combinationVariant != chunkedOutgroupSparseGridReduce && !dsg->isSubspaceDataCreated()) {
      dsg->createSubspaceData();
    }
    dsg->setZero();
  }
  // always create the extra dsgs' subspaces if they don't exist yet
  for (auto& dsg : this->getExtraUniDSGVector()) {
    if (!dsg->isSubspaceDataCreated()) {
      dsg->createSubspaceData();
    }
    dsg->setZero();
  }
}
} /* namespace combigrid */
