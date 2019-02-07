/*
 * CombiParameters.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: heenemo
 */

#ifndef SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_
#define SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_

#include <boost/serialization/map.hpp>
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
namespace combigrid {

class CombiParameters {
 public:
  CombiParameters()
      : procsSet_(false), applicationComm_(MPI_COMM_NULL), applicationCommSet_(false) {}

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax, std::vector<bool>& boundary,
                  std::vector<LevelVector>& levels, std::vector<real>& coeffs,
                  std::vector<int>& taskIDs, IndexType numberOfCombinations, IndexType numGrids = 1,
                  LevelVector reduceCombinationDimsLmin = std::vector<IndexType>(0),
                  LevelVector reduceCombinationDimsLmax = std::vector<IndexType>(0))
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        procsSet_(false),
        applicationComm_(MPI_COMM_NULL),
        applicationCommSet_(false),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax) {
    hierarchizationDims_ = std::vector<bool>(dim_, true);
    setLevelsCoeffs(taskIDs, levels, coeffs);
    numTasks_ = taskIDs.size();
  }

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax, std::vector<bool>& boundary,
                  std::vector<LevelVector>& levels, std::vector<real>& coeffs,
                  std::vector<bool>& hierachizationDims, std::vector<int>& taskIDs,
                  IndexType numberOfCombinations, IndexType numGrids = 1,
                  LevelVector reduceCombinationDimsLmin = std::vector<IndexType>(0),
                  LevelVector reduceCombinationDimsLmax = std::vector<IndexType>(0))
      : dim_(dim),
        lmin_(lmin),
        lmax_(lmax),
        boundary_(boundary),
        hierarchizationDims_(hierachizationDims),
        procsSet_(false),
        applicationComm_(MPI_COMM_NULL),
        applicationCommSet_(false),
        numberOfCombinations_(numberOfCombinations),
        numGridsPerTask_(numGrids),
        reduceCombinationDimsLmin_(reduceCombinationDimsLmin),
        reduceCombinationDimsLmax_(reduceCombinationDimsLmax) {
    setLevelsCoeffs(taskIDs, levels, coeffs);
    numTasks_ = taskIDs.size();
  }

  ~CombiParameters() {}

  inline const LevelVector& getLMin() { return lmin_; }

  inline const LevelVector& getLMax() { return lmax_; }

  inline const LevelVector& getLMinReductionVector() { return reduceCombinationDimsLmin_; }

  inline const LevelVector& getLMaxReductionVector() { return reduceCombinationDimsLmax_; }

  inline const std::vector<bool>& getBoundary() { return boundary_; }

  inline real getCoeff(int taskID) { return coeffs_[taskID]; }

  inline void getCoeffs(std::vector<int>& taskIDs, std::vector<real>& coeffs) {
    for (auto it : coeffs_) {
      taskIDs.push_back(it.first);
      coeffs.push_back(it.second);
    }
  }

  inline std::map<int, real>& getCoeffsDict() { return coeffs_; }

  inline std::map<LevelVector, int>& getLevelsToIDs() { return levelsToIDs_; }

  inline void setCoeff(int taskID, real coeff) {
    coeffs_[taskID] = coeff;
    combiDictionary_[levels_[taskID]] = coeff;
  }

  inline void setLevelsCoeffs(std::vector<int>& taskIDs, std::vector<LevelVector>& levels,
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

  inline const LevelVector& getLevel(int taskID) { return levels_[taskID]; }

  inline int getID(LevelVector level) { return getLevelsToIDs()[level]; }

  inline void getLevels(std::vector<int>& taskIDs, std::vector<LevelVector>& levels) {
    for (auto it : levels_) {
      taskIDs.push_back(it.first);
      levels.push_back(it.second);
    }
  }

  inline std::map<int, LevelVector>& getLevelsDict() { return levels_; }

  inline std::map<LevelVector, real>& getCombiDict() { return combiDictionary_; }

  inline DimType getDim() { return dim_; }

  inline size_t getNumLevels() { return levels_.size(); }
  /**
   * this method returns the number of grids a task contains
   * in case we have multiple grids in our simulation
   */
  inline IndexType getNumGrids() { return numGridsPerTask_; }
  /**
   * this method returns the number of tasks also referred to as component grids (one task might
   * contain multiple grids)
   */
  inline IndexType getNumTasks() { return numTasks_; }

  inline const std::vector<bool>& getHierarchizationDims() { return hierarchizationDims_; }

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

 private:
  DimType dim_;

  LevelVector lmin_;

  LevelVector lmax_;

  std::vector<bool> boundary_;

  std::map<int, LevelVector> levels_;

  std::map<int, real> coeffs_;

  std::map<LevelVector, int> levelsToIDs_;

  std::map<LevelVector, real> combiDictionary_;

  std::vector<bool> hierarchizationDims_;

  IndexVector procs_;

  bool procsSet_;

  CommunicatorType applicationComm_;

  bool applicationCommSet_;

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
  ar& procs_;
  ar& procsSet_;
  ar& numberOfCombinations_;
  ar& numGridsPerTask_;
  ar& numTasks_;
  ar& reduceCombinationDimsLmin_;
  ar& reduceCombinationDimsLmax_;
}
}

#endif /* SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_ */
