#pragma once

#include "combicom/CombiCom.hpp"
#include "manager/TaskWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"

namespace combigrid {

class SparseGridWorker {
 public:
  explicit SparseGridWorker(TaskWorker& taskWorkerToReference);

  SparseGridWorker(SparseGridWorker const&) = delete;
  SparseGridWorker& operator=(SparseGridWorker const&) = delete;

  ~SparseGridWorker() = default;

  inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedUniDSGVector();

  inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getCombinedUniDSGVector() const;

  inline std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>& getExtraUniDSGVector();

  inline const std::vector<std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>>&
  getExtraUniDSGVector() const;

  inline int getNumberOfGrids() const;

  inline void initCombinedUniDSGVector(const LevelVector& lmin, LevelVector lmax,
                                const LevelVector& reduceLmaxByVector, int numGrids,
                                bool clearLevels = false);

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
};

inline SparseGridWorker::SparseGridWorker(TaskWorker& taskWorkerToReference)
    : taskWorkerRef_(taskWorkerToReference) {}

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

inline int SparseGridWorker::getNumberOfGrids() const { return this->combinedUniDSGVector_.size(); }

/** Initializes the dsgu for each species by setting the subspace sizes of all
 * dfgs in the global reduce comm. After calling, all workers which share the
 * same spatial distribution of the dsgu (those who combine during global
 * reduce) have the same sized subspaces and thus share the same sized dsg.
 *
 * Attention: No data is created here, only subspace sizes are shared.
 */
inline void SparseGridWorker::initCombinedUniDSGVector(const LevelVector& lmin, LevelVector lmax,
                                                const LevelVector& reduceLmaxByVector, int numGrids,
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
  for (int g = 0; g < numGrids; ++g) {
    combinedUniDSGVector_.emplace_back(std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>(
        new DistributedSparseGridUniform<CombiDataType>(static_cast<DimType>(lmax.size()), lmax,
                                                        lmin, theMPISystem()->getLocalComm())));
  }

  // register dsgs in all dfgs
  Stats::startEvent("register dsgus");
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
    combinedUniDSGVector_[g]->createKahanBuffer();
  }
  Stats::stopEvent("register dsgus");

  // global reduce of subspace sizes
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  for (auto& uniDSG : combinedUniDSGVector_) {
    CombiCom::reduceSubspaceSizes(*uniDSG, globalReduceComm);
  }
}

} /* namespace combigrid */
