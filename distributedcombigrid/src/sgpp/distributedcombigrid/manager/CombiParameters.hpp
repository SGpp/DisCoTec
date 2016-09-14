/*
 * CombiParameters.hpp
 *
 *  Created on: Dec 8, 2015
 *      Author: heenemo
 */

#ifndef SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_
#define SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_

#include <boost/serialization/map.hpp>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

class CombiParameters {
 public:
  CombiParameters() {
  }

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax,
                  std::vector<bool>& boundary, std::vector<LevelVector>& levels,
                  std::vector<real>& coeffs, std::vector<int>& taskIDs ) :
    dim_(dim), lmin_(lmin), lmax_(lmax), boundary_(boundary)
  {
    hierarchizationDims_ = std::vector<bool>(dim_,true);
    setLevelsCoeffs( taskIDs, levels, coeffs );
  }

  CombiParameters(DimType dim, LevelVector lmin, LevelVector lmax,
                  std::vector<bool>& boundary, std::vector<LevelVector>& levels,
                  std::vector<real>& coeffs, std::vector<bool>& hierachizationDims,
                  std::vector<int>& taskIDs ) :
    dim_(dim), lmin_(lmin), lmax_(lmax), boundary_(boundary),
    hierarchizationDims_(hierachizationDims)
  {
    setLevelsCoeffs( taskIDs, levels, coeffs );
  }

  ~CombiParameters() {
  }

  inline const LevelVector& getLMin() {
    return lmin_;
  }

  inline const LevelVector& getLMax() {
    return lmax_;
  }

  inline const std::vector<bool>& getBoundary() {
    return boundary_;
  }

  inline real getCoeff(int taskID) {
    return coeffs_[taskID];
  }

  inline void getCoeffs(std::vector<int>& taskIDs, std::vector<real>& coeffs) {
    for (auto it : coeffs_) {
      taskIDs.push_back( it.first );
      coeffs.push_back( it.second );
    }
  }

  inline std::map<int, real>& getCoeffsDict() {
    return coeffs_;
  }

  inline std::map<LevelVector, int>& getLevelsToIDs() {
    return levelsToIDs_;
  }

  inline void setCoeff(int taskID, real coeff) {
    coeffs_[taskID] = coeff;
    combiDictionary_[levels_[taskID]] = coeff;
  }

  inline void setLevelsCoeffs(std::vector<int>& taskIDs,
                              std::vector<LevelVector>& levels, std::vector<real>& coeffs) {
    assert(taskIDs.size() == coeffs.size());
    assert(taskIDs.size() == levels.size());

    for (size_t i = 0; i < taskIDs.size(); ++i) {
      coeffs_[taskIDs[i]] = coeffs[i];
      levels_[taskIDs[i]] = levels[i];
      levelsToIDs_[levels[i]] = taskIDs[i];
      combiDictionary_[levels[i]] = coeffs[i];
    }
  }

  inline const LevelVector& getLevel( int taskID ) {
    return levels_[ taskID ];
  }

  inline int getID( LevelVector level ) {
    return getLevelsToIDs()[level];
  }

  inline void getLevels(std::vector<int>& taskIDs, std::vector<LevelVector>& levels) {
    for (auto it : levels_) {
          taskIDs.push_back( it.first );
          levels.push_back( it.second );
        }
  }

  inline std::map<int, LevelVector>& getLevelsDict() {
    return levels_;
  }

  inline std::map<LevelVector, real>& getCombiDict() {
    return combiDictionary_;
  }

  inline DimType getDim() {
    return dim_;
  }

  inline size_t getNumLevels() {
    return levels_.size();
  }

  inline const std::vector<bool>& getHierarchizationDims(){
    return hierarchizationDims_;
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

  friend class boost::serialization::access;

  // serialize
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
};

template<class Archive>
void CombiParameters::serialize(Archive& ar, const unsigned int version) {
  ar& dim_;
  ar& lmin_;
  ar& lmax_;
  ar& boundary_;
  ar& levels_;
  ar& coeffs_;
  ar& hierarchizationDims_;
}

}

#endif /* SRC_SGPP_COMBIGRID_MANAGER_COMBIPARAMETERS_HPP_ */
