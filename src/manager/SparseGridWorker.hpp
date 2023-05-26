#pragma once

#include "combicom/CombiCom.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "manager/TaskWorker.hpp"
#include "mpi/MPISystem.hpp"
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

  /* local reduce */
  inline void addFullGridsToUniformSG();

  /* free DSG memory as intermediate step */
  inline void deleteDsgsData();

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

  inline int getNumberOfGrids() const;

  inline std::unique_ptr<DistributedSparseGridUniform<CombiDataType>>&
  getSparseGridToUseForThirdLevel(bool thirdLevelExtraSparseGrid);

  inline void initCombinedUniDSGVector(const LevelVector& lmin, LevelVector lmax,
                                       const LevelVector& reduceLmaxByVector, int numGrids,
                                       bool clearLevels = false);

  inline void integrateCombinedSolutionToTasks(
      const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
      const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin);

  inline void interpolateAndPlotOnLevel(
      const std::string& filename, const LevelVector& levelToEvaluate,
      const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
      const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin,
      const std::vector<int>& parallelization, const std::vector<LevelVector>& decomposition) const;

  inline void readDSGsFromDisk(const std::string& filenamePrefix, bool alwaysReadFullDSG = false);

  inline void readDSGsFromDiskAndReduce(const std::string& filenamePrefixToRead,
                                        bool alwaysReadFullDSG = false);

  inline MPI_Request readReduceStartBroadcastDSGs(const std::string& filenamePrefixToRead,
                                                  bool overwrite);

  inline void reduceSubspaceSizes(const std::string& filenameToRead, bool extraSparseGrid,
                                  bool overwrite);
  /* global reduction between process groups */
  inline void reduceUniformSG(RankType globalReduceRankThatCollects = MPI_PROC_NULL);

  inline void setExtraSparseGrid(bool initializeSizes = true);

  inline void writeDSGsToDisk(std::string filenamePrefix);

  inline void writeMinMaxCoefficients(std::string fileNamePrefix) const;

  inline void writeSubspaceSizesToFile(const std::string& filenamePrefixToWrite) const;

  inline void zeroDsgsData();

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
  inline void fillDFGFromDSGU(DistributedFullGrid<CombiDataType>& dfg, int g,
                              const std::vector<BoundaryType>& boundary,
                              const std::vector<bool>& hierarchizationDims,
                              const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                              const LevelVector& lmin) const;
};

inline SparseGridWorker::SparseGridWorker(TaskWorker& taskWorkerToReference)
    : taskWorkerRef_(taskWorkerToReference) {}

inline void SparseGridWorker::addFullGridsToUniformSG() {
  assert(this->getNumberOfGrids() > 0 &&
         "Initialize dsgu first with "
         "initCombinedUniDSGVector()");
  auto numGrids = this->getNumberOfGrids();
  for (const auto& t : this->taskWorkerRef_.getTasks()) {
    for (int g = 0; g < numGrids; ++g) {
      const DistributedFullGrid<CombiDataType>& dfg =
          t->getDistributedFullGrid(static_cast<int>(g));

      // lokales reduce auf sg
      this->getCombinedUniDSGVector()[g]->addDistributedFullGrid(dfg, t->getCoefficient());
    }
  }
}

inline void SparseGridWorker::deleteDsgsData() {
  for (auto& dsg : this->getCombinedUniDSGVector()) dsg->deleteSubspaceData();
  for (auto& dsg : this->getExtraUniDSGVector()) dsg->deleteSubspaceData();
}

inline void SparseGridWorker::fillDFGFromDSGU(
    DistributedFullGrid<CombiDataType>& dfg, int g, const std::vector<BoundaryType>& boundary,
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
  for (int g = 0; g < this->getNumberOfGrids(); g++) {
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

inline int SparseGridWorker::getNumberOfGrids() const { return this->combinedUniDSGVector_.size(); }

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
                                                       int numGrids, bool clearLevels) {
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

  // global reduce of subspace sizes
  CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
  for (auto& uniDSG : combinedUniDSGVector_) {
    CombiCom::reduceSubspaceSizes(*uniDSG, globalReduceComm);
  }
}

inline void SparseGridWorker::integrateCombinedSolutionToTasks(
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) {
  for (auto& taskToUpdate : this->taskWorkerRef_.getTasks()) {
    for (int g = 0; g < this->getNumberOfGrids(); ++g) {
      // fill dfg with hierarchical coefficients from distributed sparse grid
      taskToUpdate->getDistributedFullGrid(g).extractFromUniformSG(
          *this->getCombinedUniDSGVector()[g]);
    }
  }
  this->taskWorkerRef_.dehierarchizeFullGrids(boundary, hierarchizationDims, hierarchicalBases,
                                              lmin);
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
  for (IndexType g = 0; g < this->getNumberOfGrids(); ++g) {  // loop over all grids and plot them
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

inline void SparseGridWorker::readDSGsFromDisk(const std::string& filenamePrefix,
                                               bool alwaysReadFullDSG) {
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getExtraUniDSGVector().size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = this->getExtraUniDSGVector()[i].get();
    }
    DistributedSparseGridIO::readOneFile(*dsgToUse, filenamePrefix + "_" + std::to_string(i));
    if (this->getExtraUniDSGVector().size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }
  }
}

inline void SparseGridWorker::readDSGsFromDiskAndReduce(const std::string& filenamePrefixToRead,
                                                        bool alwaysReadFullDSG) {
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getExtraUniDSGVector().size() > 0 && !alwaysReadFullDSG) {
      dsgToUse = this->getExtraUniDSGVector()[i].get();
    }
    // assume that at least for four process groups, we should have enough spare RAM
    // to read all of the sparse grid at once
    // if fewer, chunk the read/reduce
    int numberReduceChunks = 1;
    if (theMPISystem()->getNumGroups() == 1) {
      numberReduceChunks = 4;
    } else if (theMPISystem()->getNumGroups() < 4) {
      numberReduceChunks = 2;
    }
    DistributedSparseGridIO::readOneFileAndReduce(
        *dsgToUse, filenamePrefixToRead + "_" + std::to_string(i), numberReduceChunks);
    if (this->getExtraUniDSGVector().size() > 0) {
      // copy partial data from extraDSG back to uniDSG
      uniDsg->copyDataFrom(*dsgToUse);
    }
  }
}

inline MPI_Request SparseGridWorker::readReduceStartBroadcastDSGs(
    const std::string& filenamePrefixToRead, bool overwrite) {
  if (overwrite) {
    this->readDSGsFromDisk(filenamePrefixToRead);
  } else {
    this->readDSGsFromDiskAndReduce(filenamePrefixToRead);
  }
  if (this->getNumberOfGrids() != 1) {
    throw std::runtime_error("Combining more than one DSG is not implemented yet");
  }
  // distribute solution in globalReduceComm to other pgs
  return CombiCom::asyncBcastDsgData(*this->getCombinedUniDSGVector()[0],
                                     theMPISystem()->getGlobalReduceRank(),
                                     theMPISystem()->getGlobalReduceComm());
}

inline void SparseGridWorker::reduceSubspaceSizes(const std::string& filenameToRead,
                                                  bool extraSparseGrid, bool overwrite) {
  if (extraSparseGrid) {
    OUTPUT_GROUP_EXCLUSIVE_SECTION {
      this->setExtraSparseGrid(true);
#ifndef NDEBUG
      // duplicate subspace sizes to validate later
      std::vector<SubspaceSizeType> subspaceSizesToValidate =
          this->getExtraUniDSGVector()[0]->getSubspaceDataSizes();
#endif
      // use extra sparse grid
      if (overwrite) {
        DistributedSparseGridIO::readSubspaceSizesFromFile(*this->getExtraUniDSGVector()[0],
                                                           filenameToRead, false);
      } else {
        auto minFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
          return std::min(a, b);
        };
        DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
            *this->getExtraUniDSGVector()[0], filenameToRead, minFunctionInstantiation, 0, false);
      }
#ifndef NDEBUG
      assert(subspaceSizesToValidate.size() ==
             this->getExtraUniDSGVector()[0]->getSubspaceDataSizes().size());
      for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
        assert(this->getExtraUniDSGVector()[0]->getSubspaceDataSizes()[i] == 0 ||
               this->getExtraUniDSGVector()[0]->getSubspaceDataSizes()[i] ==
                   subspaceSizesToValidate[i]);
      }
      auto numDOFtoValidate =
          std::accumulate(subspaceSizesToValidate.begin(), subspaceSizesToValidate.end(), 0);
      auto numDOFnow =
          std::accumulate(this->getExtraUniDSGVector()[0]->getSubspaceDataSizes().begin(),
                          this->getExtraUniDSGVector()[0]->getSubspaceDataSizes().end(), 0);
      assert(numDOFtoValidate >= numDOFnow);
#endif
    }
  } else {
    if (!this->getExtraUniDSGVector().empty()) {
      throw std::runtime_error(
          "this->getExtraUniDSGVector() not empty, but extraSparseGrid is false");
    }
#ifndef NDEBUG
    std::vector<SubspaceSizeType> subspaceSizesToValidate =
        this->getCombinedUniDSGVector()[0]->getSubspaceDataSizes();
#endif
    FIRST_GROUP_EXCLUSIVE_SECTION {
      if (overwrite) {
        DistributedSparseGridIO::readSubspaceSizesFromFile(*this->getCombinedUniDSGVector()[0],
                                                           filenameToRead, true);
      } else {
        // if no extra sparse grid, max-reduce the normal one
        auto maxFunctionInstantiation = [](SubspaceSizeType a, SubspaceSizeType b) {
          return std::max(a, b);
        };
        DistributedSparseGridIO::readReduceSubspaceSizesFromFile(
            *this->getCombinedUniDSGVector()[0], filenameToRead, maxFunctionInstantiation, 0, true);
      }
      if (theMPISystem()->getGlobalReduceRank() != 0) {
        throw std::runtime_error("read rank is not the global reduce rank");
      }
    }
    else {
      if (theMPISystem()->getGlobalReduceRank() == 0) {
        throw std::runtime_error("read rank IS the global reduce rank");
      }
    }
    // reduce to all other process groups
    CommunicatorType globalReduceComm = theMPISystem()->getGlobalReduceComm();
    RankType senderRank = 0;
    CombiCom::broadcastSubspaceSizes(*this->getCombinedUniDSGVector()[0], globalReduceComm,
                                     senderRank);
#ifndef NDEBUG
    assert(subspaceSizesToValidate.size() ==
           this->getCombinedUniDSGVector()[0]->getSubspaceDataSizes().size());
    for (size_t i = 0; i < subspaceSizesToValidate.size(); ++i) {
      assert(subspaceSizesToValidate[i] == 0 ||
             subspaceSizesToValidate[i] ==
                 this->getCombinedUniDSGVector()[0]->getSubspaceDataSizes()[i]);
    }
    auto numDOFtoValidate =
        std::accumulate(subspaceSizesToValidate.begin(), subspaceSizesToValidate.end(), 0);
    auto numDOFnow =
        std::accumulate(this->getCombinedUniDSGVector()[0]->getSubspaceDataSizes().begin(),
                        this->getCombinedUniDSGVector()[0]->getSubspaceDataSizes().end(), 0);
    assert(numDOFtoValidate <= numDOFnow);
#endif
  }
}

inline void SparseGridWorker::reduceUniformSG(RankType globalReduceRankThatCollects) {
  auto numGrids = this->getNumberOfGrids();
  for (int g = 0; g < numGrids; ++g) {
    CombiCom::distributedGlobalReduce(*this->getCombinedUniDSGVector()[g],
                                      globalReduceRankThatCollects);
    assert(CombiCom::sumAndCheckSubspaceSizes(*this->getCombinedUniDSGVector()[g]));
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
      for (size_t i = 0; i < extraUniDSG->getNumSubspaces(); ++i) {
        extraUniDSG->setDataSize(i, this->getCombinedUniDSGVector()[0]->getDataSize(i));
      }
    }
    // level vectors are not required; read from the initial sparse grid if needed
    extraUniDSG->resetLevels();
  }
}

inline void SparseGridWorker::writeDSGsToDisk(std::string filenamePrefix) {
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    auto filename = filenamePrefix + "_" + std::to_string(i);
    auto uniDsg = this->getCombinedUniDSGVector()[i].get();
    auto dsgToUse = uniDsg;
    if (this->getExtraUniDSGVector().size() > 0) {
      dsgToUse = this->getExtraUniDSGVector()[i].get();
      dsgToUse->copyDataFrom(*uniDsg);
    }
    DistributedSparseGridIO::writeOneFile(*dsgToUse, filename);
  }
}

inline void SparseGridWorker::writeMinMaxCoefficients(std::string fileNamePrefix) const {
  for (size_t i = 0; i < this->getNumberOfGrids(); ++i) {
    DistributedSparseGridIO::writeMinMaxCoefficents(*(this->getCombinedUniDSGVector()[i]),
                                                    fileNamePrefix, i);
  }
}

inline void SparseGridWorker::writeSubspaceSizesToFile(
    const std::string& filenamePrefixToWrite) const {
  DistributedSparseGridIO::writeSubspaceSizesToFile(*this->getCombinedUniDSGVector()[0],
                                                    filenamePrefixToWrite);
}

inline void SparseGridWorker::zeroDsgsData() {
  for (auto& dsg : this->getCombinedUniDSGVector()) dsg->setZero();
  for (auto& dsg : this->getExtraUniDSGVector()) dsg->setZero();
}
} /* namespace combigrid */
