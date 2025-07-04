#ifndef SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_
#define SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_

#include <boost/serialization/map.hpp>

#include "hierarchization/CombiLinearBasisFunction.hpp"
#include "../mpi/MPISystem.hpp"
#include "../utils/LevelSetUtils.hpp"
#include "../utils/LevelVector.hpp"
#include "../utils/PowerOfTwo.hpp"
#include "../utils/Types.hpp"
namespace combigrid {

/**
 * @brief Class for the parameters of the combination technique
 *
 * This class contains all the parameters that describe a combination technique setup.
 * The class is serializable and can be sent over MPI.
 * It is used, for instance, by ProcessManager to distribute the combination technique information
 * to the ProcessGroupWorkers.
 */
class CombiParameters {
 public:
  /**
   * @brief default constructor for serialization
   */
  CombiParameters() = default;

  /**
   * @brief Constructor for CombiParameters
   *
   * the constructor variant with levels, coeffs and taskIDs specified, such that Tasks can be
   * created and tracked suitable for manager-worker setups
   *
   * @param dim dimensionality of the PDE solver domain
   * @param lmin minimum level vector, size \p dim
   * @param lmax maximum level vector, size \p dim
   * @param boundary boundary conditions, size \p dim
   * @param levels level vectors for each task
   * @param coeffs coefficients for each task
   * @param taskIDs IDs for each task
   * @param numberOfCombinations total number of combinations, only used if
   * reduceCombinationDimsLmax is used
   * @param numGrids number of grids per task, \see Task::getDistributedFullGrid()
   * @param combinationVariant variant of the combination reduction
   * @param parallelization number of processes in (Cartesian) communicator
   * @param reduceCombinationDimsLmin unused
   * @param reduceCombinationDimsLmax reduce the maximum level vector for the distributed sparse
   * grid (to avoid unnecessary storage / communication use (1, 1, 1, ...))
   * @param sizeForChunkedCommunicationInMebibyte size of the chunks for chunked reductions in MiB
   * @param forwardDecomposition in case of ambiguity, decompose the full grid such that the middle
   * grid points are on the "lower" process
   * @param thirdLevelHost hostname of the third level manager (deprecated)
   * @param thirdLevelPort port to connect to third level manager (deprecated)
   * @param thirdLevelPG process group that should connect to the third level manager (deprecated)
   */
  CombiParameters(DimType dim, const LevelVector& lmin, const LevelVector& lmax,
                  const std::vector<BoundaryType>& boundary, const std::vector<LevelVector>& levels,
                  const std::vector<real>& coeffs, const std::vector<size_t>& taskIDs,
                  size_t numberOfCombinations, IndexType numGrids = 1,
                  CombinationVariant combinationVariant = CombinationVariant::sparseGridReduce,
                  const std::vector<int>& parallelization = {0},
                  const LevelVector& reduceCombinationDimsLmin = LevelVector(0),
                  const LevelVector& reduceCombinationDimsLmax = LevelVector(0),
                  uint32_t sizeForChunkedCommunicationInMebibyte = 64,
                  bool forwardDecomposition = false, const std::string& thirdLevelHost = "",
                  unsigned short thirdLevelPort = 0, size_t thirdLevelPG = 0)
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        forwardDecomposition_(forwardDecomposition),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        combinationVariant_{combinationVariant},
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax),
        sizeForChunkedCommunicationInMebibyte_{sizeForChunkedCommunicationInMebibyte},
        thirdLevelHost_(thirdLevelHost),
        thirdLevelPort_(thirdLevelPort),
        thirdLevelPG_(thirdLevelPG) {
    hierarchizationDims_ = std::vector<bool>(dim_, true);
    for (DimType d = 0; d < dim_; ++d) {
      if (boundary_[d] == 1) {
        hierarchicalBases_.push_back(new HierarchicalHatPeriodicBasisFunction());
      } else {
        hierarchicalBases_.push_back(new HierarchicalHatBasisFunction());
      }
    }
    setLevelsCoeffs(taskIDs, levels, coeffs);
    if (parallelization != std::vector<int>({0})) {
      this->setParallelization(parallelization);
    }
  }

  /**
   * @brief Constructor for CombiParameters
   *
   * constructor variant w/o combination scheme specified -- the workers have their partial list
   * imlpicitly as tasks vector, e.g. worker-only setups
   *
   * @param dim dimensionality of the PDE solver domain
   * @param lmin minimum level vector, size \p dim
   * @param lmax maximum level vector, size \p dim
   * @param boundary boundary conditions, size \p dim
   * @param numberOfCombinations total number of combinations, only used if
   * reduceCombinationDimsLmax is used
   * @param numGrids number of grids per task, \see Task::getDistributedFullGrid()
   * @param combinationVariant variant of the combination reduction
   * @param parallelization number of processes in (Cartesian) communicator
   * @param reduceCombinationDimsLmin unused
   * @param reduceCombinationDimsLmax reduce the maximum level vector for the distributed sparse
   * grid (to avoid unnecessary storage / communication use (1, 1, 1, ...))
   * @param sizeForChunkedCommunicationInMebibyte size of the chunks for chunked reductions in MiB
   * @param forwardDecomposition in case of ambiguity, decompose the full grid such that the middle
   * grid points are on the "lower" process
   * @param thirdLevelHost hostname of the third level manager (deprecated)
   * @param thirdLevelPort port to connect to third level manager (deprecated)
   * @param thirdLevelPG process group that should connect to the third level manager (deprecated)
   */
  CombiParameters(DimType dim, const LevelVector& lmin, const LevelVector& lmax,
                  const std::vector<BoundaryType>& boundary, size_t numberOfCombinations,
                  IndexType numGrids = 1,
                  CombinationVariant combinationVariant = CombinationVariant::sparseGridReduce,
                  const std::vector<int>& parallelization = {0},
                  const LevelVector& reduceCombinationDimsLmin = LevelVector(0),
                  const LevelVector& reduceCombinationDimsLmax = LevelVector(0),
                  uint32_t sizeForChunkedCommunicationInMebibyte = 64,
                  bool forwardDecomposition = false, const std::string& thirdLevelHost = "",
                  unsigned short thirdLevelPort = 0, size_t thirdLevelPG = 0)
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        forwardDecomposition_(forwardDecomposition),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        combinationVariant_{combinationVariant},
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax),
        sizeForChunkedCommunicationInMebibyte_{sizeForChunkedCommunicationInMebibyte},
        thirdLevelHost_(thirdLevelHost),
        thirdLevelPort_(thirdLevelPort),
        thirdLevelPG_(thirdLevelPG) {
    hierarchizationDims_ = std::vector<bool>(dim_, true);
    for (DimType d = 0; d < dim_; ++d) {
      if (boundary_[d] == 1) {
        hierarchicalBases_.push_back(new HierarchicalHatPeriodicBasisFunction());
      } else {
        hierarchicalBases_.push_back(new HierarchicalHatBasisFunction());
      }
    }
    if (parallelization != std::vector<int>({0})) {
      this->setParallelization(parallelization);
    }
  }

  /**
   * @brief Constructor for CombiParameters
   *
   * constructor variant with hierarchizationDims specified, to only hierarchize certain dimensions
   *
   * @param dim dimensionality of the PDE solver domain
   * @param lmin minimum level vector, size \p dim
   * @param lmax maximum level vector, size \p dim
   * @param boundary boundary conditions, size \p dim
   * @param levels level vectors for each task
   * @param coeffs coefficients for each task
   * @param hierarchizationDims vector of bools, whether to hierarchize the dimension
   * @param taskIDs IDs for each task
   * @param numberOfCombinations total number of combinations, only used if
   * reduceCombinationDimsLmax is used
   * @param numGrids number of grids per task, \see Task::getDistributedFullGrid()
   * @param combinationVariant variant of the combination reduction
   * @param reduceCombinationDimsLmin unused
   * @param reduceCombinationDimsLmax reduce the maximum level vector for the distributed sparse
   * grid (to avoid unnecessary storage / communication use (1, 1, 1, ...))
   * @param sizeForChunkedCommunicationInMebibyte size of the chunks for chunked reductions in MiB
   * @param forwardDecomposition in case of ambiguity, decompose the full grid such that the middle
   * grid points are on the "lower" process
   * @param thirdLevelHost hostname of the third level manager (deprecated)
   * @param thirdLevelPort port to connect to third level manager (deprecated)
   * @param thirdLevelPG process group that should connect to the third level manager (deprecated)
   */
  CombiParameters(DimType dim, const LevelVector& lmin, const LevelVector& lmax,
                  const std::vector<BoundaryType>& boundary, const std::vector<LevelVector>& levels,
                  const std::vector<real>& coeffs, const std::vector<bool>& hierarchizationDims,
                  const std::vector<size_t>& taskIDs, size_t numberOfCombinations,
                  IndexType numGrids = 1,
                  CombinationVariant combinationVariant = CombinationVariant::sparseGridReduce,
                  const LevelVector& reduceCombinationDimsLmin = LevelVector(0),
                  const LevelVector& reduceCombinationDimsLmax = LevelVector(0),
                  uint32_t sizeForChunkedCommunicationInMebibyte = 64,
                  bool forwardDecomposition = !isGENE, const std::string& thirdLevelHost = "",
                  unsigned short thirdLevelPort = 0, size_t thirdLevelPG = 0)
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        hierarchizationDims_(hierarchizationDims),
        forwardDecomposition_(forwardDecomposition),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        combinationVariant_{combinationVariant},
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax),
        sizeForChunkedCommunicationInMebibyte_{sizeForChunkedCommunicationInMebibyte},
        thirdLevelHost_(thirdLevelHost),
        thirdLevelPort_(thirdLevelPort),
        thirdLevelPG_(thirdLevelPG) {
    for (DimType d = 0; d < dim_; ++d) {
      if (boundary_[d] == 1) {
        hierarchicalBases_.push_back(new HierarchicalHatPeriodicBasisFunction());
      } else {
        hierarchicalBases_.push_back(new HierarchicalHatBasisFunction());
      }
    }
    setLevelsCoeffs(taskIDs, levels, coeffs);
  }

  ~CombiParameters() = default;
  // for (auto& b : hierarchicalBases_) {
  //   if (b!= nullptr) { delete b; }
  //   b = nullptr;
  // }
  // }

  inline const LevelVector& getLMin() const { return lmin_; }

  inline const LevelVector& getLMax() const { return lmax_; }

  inline const LevelVector& getLMaxReductionVector() const { return reduceCombinationDimsLmax_; }

  /**
   * @brief get the boundary flags
   */
  inline const std::vector<BoundaryType>& getBoundary() const { return boundary_; }

  /**
   * @brief get the combination coefficient for a specific task
   */
  inline real getCoeff(size_t taskID) const {
    if (coeffs_.find(taskID) == coeffs_.end()) {
      return 0;
    } else {
      return coeffs_.at(taskID);
    }
  }

  /**
   * @brief get the combination coefficients of all tasks
   */
  inline void getCoeffs(std::vector<size_t>& taskIDs, std::vector<real>& coeffs) const {
    for (const auto& it : coeffs_) {
      taskIDs.push_back(it.first);
      coeffs.push_back(it.second);
    }
  }

  /**
   * @brief get a mapping from levels to task IDs
   */
  inline std::map<LevelVector, size_t>& getLevelsToIDs() { return levelsToIDs_; }

  /**
   * @brief set a new combination coefficient for a task
   *
   * only used for fault tolerance
   */
  inline void setCoeff(size_t taskID, real coeff) {
    assert(coeffs_.find(taskID) != coeffs_.end());
    coeffs_[taskID] = coeff;
    combiDictionary_[levels_[taskID]] = coeff;
  }

  /**
   * @brief set the levels and coefficients for the tasks
   */
  inline void setLevelsCoeffs(const std::vector<size_t>& taskIDs,
                              const std::vector<LevelVector>& levels,
                              const std::vector<real>& coeffs) {
    if (taskIDs.empty() && ENABLE_FT) {
      throw std::runtime_error(
          "CombiParameters::setLevelsCoeffs: taskIDs is empty. Not sure if this should be possible "
          "but this is an error for now so I can see where it happens.");
    }
#ifndef NDEBUG
    assert(taskIDs.size() == coeffs.size() || taskIDs.empty());
    auto sorted = taskIDs;
    std::sort(sorted.begin(), sorted.end());
    assert(std::unique(sorted.begin(), sorted.end()) == sorted.end());
    assert(coeffs.size() == levels.size());
#endif  // NDEBUG

    if (taskIDs.empty()) {
      assert(false);
    } else {
      for (size_t i = 0; i < taskIDs.size(); ++i) {
        coeffs_[taskIDs[i]] = coeffs[i];
        levels_[taskIDs[i]] = levels[i];
        levelsToIDs_[levels[i]] = taskIDs[i];
        combiDictionary_[levels[i]] = coeffs[i];
      }
    }
  }

  /**
   * @brief get the level vector for a specific task by its \p taskID
   */
  inline const LevelVector& getLevel(size_t taskID) const {
    static thread_local LevelVector emptyLevelVector;
    emptyLevelVector.clear();
    if (levels_.find(taskID) == levels_.end()) {
      return emptyLevelVector;
    } else {
      return levels_.at(taskID);
    }
  }

  /**
   * @brief get the task ID for a specific level vector
   */
  inline size_t getID(LevelVector level) { return getLevelsToIDs()[level]; }

  /**
   * @brief get a reference to a mappping of task IDs to level vectors
   */
  inline std::map<size_t, LevelVector>& getLevelsDict() { return levels_; }

  /**
   * @brief get a reference to a mapping of level vectors to combination coefficients
   */
  inline std::map<LevelVector, real>& getCombiDict() { return combiDictionary_; }

  /**
   * @brief get the dimensionality of the PDE solver domain
   */
  inline DimType getDim() const { return dim_; }

  /**
   * @brief get the number of levels / component grids in the combination technique
   *
   * with fault tolerance, this number can change and may include grids with zero-coefficients
   */
  inline size_t getNumLevels() const { return levels_.size(); }

  /**
   * this method returns the number of grids a task contains
   * in case we have multiple grids in our simulation, for instance multiple particle species in a
   * plasma turbulence simulation
   */
  inline size_t getNumGrids() const { return numGridsPerTask_; }

  /**
   * @brief get the dimensions in which hierarchization should be performed, i.e., the combination
   * technique will be used
   */
  inline const std::vector<bool>& getHierarchizationDims() const { return hierarchizationDims_; }

  /**
   * @brief Set the Hierarchical Bases object
   *
   * set a vector of hierarchicas bases, one for each dimension
   * (not necessary if using hierarchical hats in all dimensions)
   * Takes over ownership of the contents of the bases object,
   * which are assumed to be on the heap
   */
  inline void setHierarchicalBases(std::vector<BasisFunctionBasis*>& bases) {
    assert(bases.size() == dim_);
    // delete old hierarchicalBases_
    for (auto& b : hierarchicalBases_) {
      if (b != nullptr) {
        delete b;
      }
      b = nullptr;
    }
    for (size_t i = 0; i < bases.size(); ++i) {
      if (bases[i] == nullptr) {
        assert(hierarchizationDims_[i] == false);
      }
      hierarchicalBases_[i] = bases[i];
      bases[i] = nullptr;
    }
  }

  /**
   * @brief Get the hierarchical bases, one for each dimension
   *        (assuming all the dfgs are using the same number of dimensions and the same bases)
   *
   * @return std::vector<BasisFunctionBasis*> pointers of the type of basis function
   *          may be nullptr or anything for a non-hierarchization dimension
   */
  inline const std::vector<BasisFunctionBasis*>& getHierarchicalBases() const {
    assert(hierarchicalBases_.size() == dim_);
    return hierarchicalBases_;
  }

  /**
   * get the common (Cartesian) parallelization, assuming all component grids share the same
   * parallelization
   *
   * @return const std::vector<int>& the number of processes in each dimension
   */
  inline const std::vector<int>& getParallelization() const {
    assert(uniformDecomposition);
    return procs_;
  }

  /**
   * @brief get the number of combinations
   */
  inline const size_t& getNumberOfCombinations() const { return numberOfCombinations_; }

  /**
   * @brief get the combination variant used
   */
  inline CombinationVariant getCombinationVariant() const { return combinationVariant_; }

  /**
   * @brief get the reduction chunk size in MiB per OpenMP thread
   *
   * for chunked reductions
   */
  inline uint32_t getChunkSizeInMebibybtePerThread() const {
    if (sizeForChunkedCommunicationInMebibyte_ == 0) {
      throw std::runtime_error("chunkSizeInMebibytePerThread_ is not set");
    }
    return sizeForChunkedCommunicationInMebibyte_;
  }

  /**
   * @brief get the URL where the third level manager is running
   */
  inline const std::string& getThirdLevelHost() { return thirdLevelHost_; }

  /**
   * @brief get the port where the third level manager is running
   */
  inline unsigned short getThirdLevelPort() { return thirdLevelPort_; }

  /**
   * @brief get the process group that should connect to the third level manager
   */
  inline size_t getThirdLevelPG() { return thirdLevelPG_; }

  /**
   * @brief set the parallelization
   */
  inline void setParallelization(const std::vector<int>& p) {
    assert(uniformDecomposition);

    procs_ = p;
  }

  inline bool isParallelizationSet() const { return !procs_.empty(); }

  /**
   * @brief Set the Decomposition, specified by 1d point indices at maximum resolution level
   *
   * @param decomposition a vector of index vectors, specifying for each dimension
   *        the lowest 1d index (on a full grid of level lmax)
   *        that will belong to each cartesian communicator slice
   */
  inline void setDecomposition(const std::vector<IndexVector>& decomposition) {
    decomposition_ = decomposition;
#ifndef NDEBUG
    assert(uniformDecomposition);
    for (DimType d = 0; d < dim_; ++d) {
      assert(decomposition[d][0] == 0);
      auto numPoints = combigrid::getNumDofNodal(lmax_[d], boundary_[d]);
      assert(decomposition[d].back() < numPoints);
      assert(static_cast<size_t>(procs_[d]) == decomposition[d].size());
    }
#endif  // not def NDEBUG
  }

  /**
   * @brief get the decomposition, specified by 1d point indices at maximum resolution level
   */
  inline const std::vector<IndexVector>& getDecomposition() const { return decomposition_; }

  inline bool getForwardDecomposition() const {
    if (isGENE) {
      assert(!forwardDecomposition_);
    }
    return forwardDecomposition_;
  }

 private:
  DimType dim_;

  LevelVector lmin_;

  LevelVector lmax_;

  std::vector<BoundaryType> boundary_;

  std::map<size_t, LevelVector> levels_;

  std::map<size_t, real> coeffs_;

  std::map<LevelVector, size_t> levelsToIDs_;

  std::map<LevelVector, real> combiDictionary_;

  std::vector<bool> hierarchizationDims_;

  std::vector<BasisFunctionBasis*> hierarchicalBases_;

  std::vector<int> procs_;

  std::vector<IndexVector> decomposition_;

  bool forwardDecomposition_ = false;

  friend class boost::serialization::access;
  size_t numberOfCombinations_;  // total number of combinations
  size_t numGridsPerTask_;       // number of grids per task

  CombinationVariant combinationVariant_;

  /**
   * This level vector indicates which dimension of lmin should be decreased by how many levels
   * for constructing the distributed sparse grid.
   * (0,0,0,0,0) would be the classical scheme were the sparse grid has the same minimum level as
   * the combi scheme
   * Higher values can decrease the communication volume but decrease the accuracy
   */
  LevelVector reduceCombinationDimsLmin_;
  /**
   * This level vector indicates which dimension of lmax should be decreased by how many levels
   * for constructing the distributed sparse grid.
   * (0,0,0,0,0,0) would be the classical scheme were the sparse grid has the same maximum level as
   * the combi scheme
   * (1,1,1,1,1,1) would be the classical optimization where we do not combine the highest
   * subspaces
   *               as they are only contained on the owning grid
   * Higher values can decrease the communication volume but decrease the accuracy
   * It is ensured that lmax >= lmin
   */
  LevelVector reduceCombinationDimsLmax_;

  uint32_t sizeForChunkedCommunicationInMebibyte_;

  std::string thirdLevelHost_;

  unsigned short thirdLevelPort_;

  size_t thirdLevelPG_;

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

template <class Archive>
void CombiParameters::serialize(Archive& ar, const unsigned int version) {
  ar & dim_;
  ar & lmin_;
  ar & lmax_;
  ar & boundary_;
  ar & levels_;
  ar & coeffs_;
  ar & hierarchizationDims_;
  ar & hierarchicalBases_;
  ar & procs_;
  ar & decomposition_;
  ar & forwardDecomposition_;
  ar & numberOfCombinations_;
  ar & numGridsPerTask_;
  ar & combinationVariant_;
  ar & sizeForChunkedCommunicationInMebibyte_;
  ar & reduceCombinationDimsLmin_;
  ar & reduceCombinationDimsLmax_;
  ar & thirdLevelHost_;
  ar & thirdLevelPort_;
  ar & thirdLevelPG_;
}

template <typename T>
static void setCombiParametersHierarchicalBasesUniform(CombiParameters& combiParameters) {
  std::vector<BasisFunctionBasis*> bases;
  for (DimType d = 0; d < combiParameters.getDim(); ++d) {
    bases.push_back(new T());
  }
  assert(bases.size() == combiParameters.getDim());
  combiParameters.setHierarchicalBases(bases);
}

inline static void setCombiParametersHierarchicalBasesUniform(CombiParameters& combiParameters,
                                                              const std::string& basisName) {
  if (basisName == "hat") {
    setCombiParametersHierarchicalBasesUniform<HierarchicalHatBasisFunction>(combiParameters);
  } else if (basisName == "hat_periodic") {
    setCombiParametersHierarchicalBasesUniform<HierarchicalHatPeriodicBasisFunction>(
        combiParameters);
  } else if (basisName == "fullweighting") {
    setCombiParametersHierarchicalBasesUniform<FullWeightingBasisFunction>(combiParameters);
  } else if (basisName == "fullweighting_periodic") {
    setCombiParametersHierarchicalBasesUniform<FullWeightingPeriodicBasisFunction>(combiParameters);
  } else if (basisName == "biorthogonal") {
    setCombiParametersHierarchicalBasesUniform<BiorthogonalBasisFunction>(combiParameters);
  } else if (basisName == "biorthogonal_periodic") {
    setCombiParametersHierarchicalBasesUniform<BiorthogonalPeriodicBasisFunction>(combiParameters);
  } else {
    throw std::invalid_argument("Hierarchical basis string not known.");
  }
}

}  // namespace combigrid

#endif /* SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_ */
