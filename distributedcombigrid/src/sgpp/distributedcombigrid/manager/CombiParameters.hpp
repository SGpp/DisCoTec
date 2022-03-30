#ifndef SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_
#define SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_

#include <boost/serialization/map.hpp>
#include "sgpp/distributedcombigrid/legacy/CombiLinearBasisFunction.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
namespace combigrid {

class CombiParameters {
 public:
  CombiParameters()
      : procsSet_(false) {}

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax, std::vector<bool>& boundary,
                  std::vector<LevelVector>& levels, std::vector<real>& coeffs,
                  std::vector<size_t>& taskIDs, IndexType numberOfCombinations, IndexType numGrids = 1,
                  const IndexVector parallelization = {0},
                  LevelVector reduceCombinationDimsLmin = std::vector<IndexType>(0),
                  LevelVector reduceCombinationDimsLmax = std::vector<IndexType>(0),
                  bool forwardDecomposition = !isGENE,
                  const std::string& thirdLevelHost = "",
                  unsigned short thirdLevelPort = 0,
                  size_t thirdLevelPG = 0
                  )
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        procsSet_(false),
        forwardDecomposition_(forwardDecomposition),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax),
        thirdLevelHost_(thirdLevelHost),
        thirdLevelPort_(thirdLevelPort),
        thirdLevelPG_(thirdLevelPG)
  {
    hierarchizationDims_ = std::vector<bool>(dim_, true);
    for (DimType d = 0; d < dim_; ++d) {
      hierarchicalBases_.push_back(new HierarchicalHatBasisFunction());
    }
    setLevelsCoeffs(taskIDs, levels, coeffs);
    numTasks_ = static_cast<long>(taskIDs.size());
    if (parallelization != IndexVector({0})){
      this->setParallelization(parallelization);
    }
  }

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax, std::vector<bool>& boundary,
                  std::vector<LevelVector>& levels, std::vector<real>& coeffs,
                  std::vector<bool>& hierarchizationDims, std::vector<size_t>& taskIDs,
                  IndexType numberOfCombinations, IndexType numGrids = 1,
                  LevelVector reduceCombinationDimsLmin = std::vector<IndexType>(0),
                  LevelVector reduceCombinationDimsLmax = std::vector<IndexType>(0),
                  bool forwardDecomposition = !isGENE, const std::string& thirdLevelHost = "",
                  unsigned short thirdLevelPort = 0, size_t thirdLevelPG = 0)
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        hierarchizationDims_(hierarchizationDims),
        procsSet_(false),
        forwardDecomposition_(forwardDecomposition),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax),
        thirdLevelHost_(thirdLevelHost),
        thirdLevelPort_(thirdLevelPort),
        thirdLevelPG_(thirdLevelPG) {
    for (DimType d = 0; d < dim_; ++d) {
      hierarchicalBases_.push_back(new HierarchicalHatBasisFunction());
    }
    setLevelsCoeffs(taskIDs, levels, coeffs);
    numTasks_ = taskIDs.size();
  }

  ~CombiParameters() {
    // for (auto& b : hierarchicalBases_) {
    //   if (b!= nullptr) { delete b; }
    //   b = nullptr;
    // }
  }

  inline const LevelVector& getLMin() { return lmin_; }

  inline const LevelVector& getLMax() { return lmax_; }

  inline const LevelVector& getLMinReductionVector() { return reduceCombinationDimsLmin_; }

  inline const LevelVector& getLMaxReductionVector() { return reduceCombinationDimsLmax_; }

  inline const std::vector<bool>& getBoundary() { return boundary_; }

  inline real getCoeff(size_t taskID) { return coeffs_[taskID]; }

  inline void getCoeffs(std::vector<size_t>& taskIDs, std::vector<real>& coeffs) {
    for (auto it : coeffs_) {
      taskIDs.push_back(it.first);
      coeffs.push_back(it.second);
    }
  }

  inline std::map<size_t, real>& getCoeffsDict() { return coeffs_; }

  inline std::map<LevelVector, size_t>& getLevelsToIDs() { return levelsToIDs_; }

  inline void setCoeff(size_t taskID, real coeff) {
    coeffs_[taskID] = coeff;
    combiDictionary_[levels_[taskID]] = coeff;
  }

  inline void setLevelsCoeffs(std::vector<size_t>& taskIDs, std::vector<LevelVector>& levels,
                              std::vector<real>& coeffs) {
    assert(taskIDs.size() == coeffs.size());
    assert(taskIDs.size() == levels.size());

    for (size_t i = 0; i < taskIDs.size(); ++i) {
      coeffs_[taskIDs[i]] = coeffs[i];
      levels_[taskIDs[i]] = levels[i];
      levelsToIDs_[levels[i]] = taskIDs[i];
      combiDictionary_[levels[i]] = coeffs[i];
    }
  }

  inline const LevelVector& getLevel(size_t taskID) { return levels_[taskID]; }

  inline size_t getID(LevelVector level) { return getLevelsToIDs()[level]; }

  inline void getLevels(std::vector<size_t>& taskIDs, std::vector<LevelVector>& levels) {
    for (auto it : levels_) {
      taskIDs.push_back(it.first);
      levels.push_back(it.second);
    }
  }

  inline std::map<size_t, LevelVector>& getLevelsDict() { return levels_; }

  inline std::map<LevelVector, real>& getCombiDict() { return combiDictionary_; }

  inline DimType getDim() { return dim_; }

  inline size_t getNumLevels() { return levels_.size(); }
  /**
   * this method returns the number of grids a task contains
   * in case we have multiple grids in our simulation
   */
  inline IndexType getNumGrids() const { return numGridsPerTask_; }
  /**
   * this method returns the number of tasks also referred to as component grids (one task might
   * contain multiple grids)
   */
  inline IndexType getNumTasks() { return numTasks_; }

  inline const std::vector<bool>& getHierarchizationDims() { return hierarchizationDims_; }

  /**
   * @brief Set the Hierarchical Bases object
   *        set a vector of hierarchicas bases, one for each dimension
   *        (not necessary if using hierarchical hats in all dimensions)
   *        Takes over ownership of the contents of the bases object,
   *        which are assumed to be on the heap
   */
  inline void setHierarchicalBases(std::vector<BasisFunctionBasis*>& bases) {
    assert(bases.size() == dim_);
    // delete old hierarchicalBases_
    for (auto& b : hierarchicalBases_) {
      if (b!= nullptr) { delete b; }
      b = nullptr;
    }
    for (size_t i = 0; i < bases.size(); ++i) {
      if (bases[i] == nullptr) {
        assert(hierarchizationDims_[i] == false);
      }
      hierarchicalBases_[i] =bases[i];
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
  inline const std::vector<BasisFunctionBasis*>& getHierarchicalBases() {
    assert(hierarchicalBases_.size() == dim_);
    return hierarchicalBases_;
  }

  /* get the common parallelization
   * this function can only be used in the uniform mode
   */
  inline const IndexVector getParallelization() const {
    assert(uniformDecomposition && procsSet_);
    return procs_;
  }

  inline const IndexType& getNumberOfCombinations() const { return numberOfCombinations_; }

  inline CommunicatorType getApplicationComm() const {
    assert(uniformDecomposition);
    return theMPISystem()->getLocalComm();
    // assert( uniformDecomposition && applicationCommSet_ );

    // return applicationComm_;
  }

  inline const std::string& getThirdLevelHost() {
    return thirdLevelHost_;
  }

  inline unsigned short getThirdLevelPort() {
    return thirdLevelPort_;
  }

  inline size_t getThirdLevelPG() {
    return thirdLevelPG_;
  }

  inline bool isApplicationCommSet() const {
    return false;
    // return applicationCommSet_;
  }

  inline void setApplicationComm(CommunicatorType comm) {
    assert(uniformDecomposition);
    return;  // outdated
    // make sure it is set only once
    if (applicationCommSet_ == true) return;

    MPI_Comm_dup(comm, &applicationComm_);
    applicationCommSet_ = true;
  }

  /* set the common parallelization
   * this function can only be used in the uniform mode
   */
  inline void setParallelization(const IndexVector p) {
    assert(uniformDecomposition);

    procs_ = p;
    procsSet_ = true;
  }

  inline bool isParallelizationSet() const {
    return procsSet_;
  }

  /**
   * @brief Set the Decomposition
   *
   * @param decomposition a vector of index vectors, specifying for each dimension
   *        the lowest 1d index (on a full grid of level lmax)
   *        that will belong to each cartesian communicator slice
   */
  inline void setDecomposition(const std::vector<IndexVector>& decomposition) {
    assert(uniformDecomposition);
    decomposition_ = decomposition;
    for (DimType d = 0; d < dim_; ++d) {
      assert(decomposition[d][0] == 0);
      auto numPoints = powerOfTwo[lmax_[d]] + (boundary_[d] ? 1 : -1);
      assert(decomposition[d].back() < numPoints);
      assert(procs_[d] == decomposition[d].size());
    }
  }

  inline const std::vector<IndexVector>& getDecomposition() const {
    return decomposition_;
  }

  inline bool getForwardDecomposition() const {
    if (isGENE){
      assert(!forwardDecomposition_);
    }
    return forwardDecomposition_;
  }

 private:
  DimType dim_;

  LevelVector lmin_;

  LevelVector lmax_;

  std::vector<bool> boundary_;

  std::map<size_t, LevelVector> levels_;

  std::map<size_t, real> coeffs_;

  std::map<LevelVector, size_t> levelsToIDs_;

  std::map<LevelVector, real> combiDictionary_;

  std::vector<bool> hierarchizationDims_;

  std::vector<BasisFunctionBasis*> hierarchicalBases_;

  IndexVector procs_;

  bool procsSet_;

  CommunicatorType applicationComm_;

  bool applicationCommSet_;

  std::vector<IndexVector> decomposition_;

  bool forwardDecomposition_;

  friend class boost::serialization::access;
  IndexType numberOfCombinations_;  // total number of combinations
  IndexType numGridsPerTask_;       // number of grids per task

  IndexType numTasks_;
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


  std::string thirdLevelHost_;

  unsigned short thirdLevelPort_;

  size_t thirdLevelPG_;

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

template <class Archive>
void CombiParameters::serialize(Archive& ar, const unsigned int version) {
  ar& dim_;
  ar& lmin_;
  ar& lmax_;
  ar& boundary_;
  ar& levels_;
  ar& coeffs_;
  ar& hierarchizationDims_;
  ar& hierarchicalBases_;
  ar& procs_;
  ar& procsSet_;
  ar& decomposition_;
  ar& forwardDecomposition_;
  ar& numberOfCombinations_;
  ar& numGridsPerTask_;
  ar& numTasks_;
  ar& reduceCombinationDimsLmin_;
  ar& reduceCombinationDimsLmax_;
  ar& thirdLevelHost_;
  ar& thirdLevelPort_;
  ar& thirdLevelPG_;
}


template<typename T>
static void setCombiParametersHierarchicalBasesUniform(CombiParameters& combiParameters) {
  std::vector<BasisFunctionBasis*> bases;
  for (DimType d = 0; d < combiParameters.getDim(); ++d) {
    bases.push_back(new T());
  }
  assert(bases.size() == combiParameters.getDim());
  combiParameters.setHierarchicalBases(bases);
}

inline static void setCombiParametersHierarchicalBasesUniform(CombiParameters& combiParameters,
                                                std::string basisName) {
  if (basisName == "hat") {
    setCombiParametersHierarchicalBasesUniform<HierarchicalHatBasisFunction>(combiParameters);
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
