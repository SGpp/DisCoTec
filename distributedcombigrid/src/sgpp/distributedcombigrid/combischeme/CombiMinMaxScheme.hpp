#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <numeric>
#include "sgpp/distributedcombigrid/legacy/combigrid_utils.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelSetUtils.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

class CombiMinMaxScheme {
 public:
  CombiMinMaxScheme(DimType dim, const LevelVector& lmin, const LevelVector& lmax) {
    assert(dim > 0);

    assert(lmax.size() == dim);
    for (size_t i = 0; i < lmax.size(); ++i) assert(lmax[i] > 0);

    assert(lmin.size() == dim);

    for (size_t i = 0; i < lmin.size(); ++i) {
      assert(lmin[i] > 0);
      assert(lmax[i] >= lmin[i]);
    }

    n_ = 0;
    dim_ = dim;
    lmin_ = lmin;
    lmax_ = lmax;

    // Calculate the effective dimension
    effDim_ = dim_;
    LevelVector diff = lmax_ - lmin_;
    for (auto i : diff)
      if (i == 0) effDim_--;
  }

  CombiMinMaxScheme(const std::vector<LevelVector>& levels, const std::vector<real>& coefficients)
    : combiSpaces_(levels), coefficients_(coefficients) {
      assert(levels.size() > 0);
      assert(levels.size() == coefficients.size());
      n_ = 0;
      dim_ = static_cast<combigrid::DimType>(levels.back().size());
    }

  virtual ~CombiMinMaxScheme() = default;

  /* Generate the combischeme corresponding to the classical combination technique.
   * We need to ensure that lmax = lmin +c*ones(dim), and take special care
   * of dummy dimensions
   * */
  void createClassicalCombischeme();

  /* Generates the adaptive combination scheme (equivalent to CK's
   * Python code)
   * */
  void createAdaptiveCombischeme();

  /* Generates the fault tolerant combination technique with extra
   * grids used in case of faults
   * */
  void makeFaultTolerant();

  inline const std::vector<LevelVector>& getCombiSpaces() const { return combiSpaces_; }

  inline const std::vector<double>& getCoeffs() const { return coefficients_; }

  inline const std::vector<LevelVector>& getDownSet() {
    return levels_;
  }

  /**
   * @brief (re-)generate the downward closed set of the CombiScheme
   *
   * (this function is here because I don't understand the levels_ in the code so far,
   * so I re-implemented the downward closed set from the definition)
   *
   */
  void createDownSet() {
    std::set<LevelVector> subspaces;
    for (size_t i = 0; i < combiSpaces_.size(); ++ i) {
      // only check the active front of the scheme -- those grids with coefficient 1
      if (coefficients_[i] == 1.) {
        // insert all subspaces covered by grid i into the set
        auto thisGridsDownSet = combigrid::getDownSet(combiSpaces_[i]);
        subspaces.insert(thisGridsDownSet.begin(), thisGridsDownSet.end());
      }
    }
    levels_.assign(subspaces.begin(), subspaces.end());
  }

  inline void print(std::ostream& os) const;

 protected:
  /* L1 norm of combispaces on the highest diagonal */
  LevelType n_;

  /* Dimension of lmin_ and lmax_ */
  DimType dim_;

  /* Number of actual combination dimensions */
  DimType effDim_;

  /* Minimal resolution */
  LevelVector lmin_;

  /* Maximal resolution */
  LevelVector lmax_;

  /* Downset */
  std::vector<LevelVector> levels_;

  /* Component grid levels of the combination technique*/
  std::vector<LevelVector> combiSpaces_;

  /* Combination coefficients */
  std::vector<real> coefficients_;

  /* Calculate the coefficients of the classical CT (binomial coefficient)*/
  void computeCombiCoeffsClassical();

  /* Calculate the coefficients of the adaptive CT using the formula in Alfredo's
   * SDC paper (from Brendan Harding)*/
  void computeCombiCoeffsAdaptive();
};

inline std::ostream& operator<<(std::ostream& os, const combigrid::CombiMinMaxScheme& scheme) {
  scheme.print(os);
  return os;
}

inline void CombiMinMaxScheme::print(std::ostream& os) const {
  for (uint i = 0; i < combiSpaces_.size(); ++i)
    os << "\t" << i << ". " << combiSpaces_[i] << "\t" << coefficients_[i] << std::endl;

  os << std::endl;
}

class CombiMinMaxSchemeFromFile : public CombiMinMaxScheme {
 public:
  CombiMinMaxSchemeFromFile(DimType dim, LevelVector& lmin, LevelVector& lmax, std::string ctschemeFile)
      : CombiMinMaxScheme(dim, lmin, lmax) {
    // read in CT scheme, if applicable
    boost::property_tree::ptree pScheme;
    boost::property_tree::json_parser::read_json(ctschemeFile, pScheme);
    for (const auto& component : pScheme.get_child("")) {
      assert(component.first.empty());  // list elements have no names
      for (const auto& c : component.second) {
        if (c.first == "coeff") {
          coefficients_.push_back(c.second.get_value<real>());
        } else if (c.first == "level") {
          LevelVector lvl(dim);
          int i = 0;
          for (const auto& l : c.second) {
            lvl[i] = l.second.get_value<LevelType>();
            ++i;
          }
          assert(lvl <= lmax);
          assert(lmin <= lvl);
          combiSpaces_.push_back(lvl);
        } else if (c.first == "group_no") {
          processGroupNumbers_.push_back(c.second.get_value<size_t>());
        } else {
          assert(false);
        }
      }
      // process group numbers should either be always given or never
      assert(coefficients_.size() == combiSpaces_.size());
      assert(processGroupNumbers_.size() == combiSpaces_.size() || processGroupNumbers_.size() == 0);
    }
    assert(coefficients_.size() > 0);
    assert(coefficients_.size() == combiSpaces_.size());
  }

  inline const std::vector<size_t>& getProcessGroupNumbers() const { return processGroupNumbers_; }

 private:
  std::vector<size_t> processGroupNumbers_;
};

inline long long int getCombiDegreesOfFreedom(const LevelVector& level,
                                              const std::vector<BoundaryType>& boundary) {
  long long int numDOF = 1;
  assert(level.size() == boundary.size());
  auto dim = static_cast<DimType>(level.size());
  for (DimType d = 0; d < dim; ++d) {
    const auto& level_i = level[d];
    assert(level_i > -1);
    if (level_i > 0) {
      numDOF *= powerOfTwo[level_i] - 1 + boundary[d];
    } else {
      numDOF *= boundary[d];
    }
  }
  return numDOF;
}

inline long long int printCombiDegreesOfFreedom(const std::vector<LevelVector>& combiSpaces,
                                                const std::vector<BoundaryType>& boundary) {
  long long int numDOF = 0;
  for (const auto& space : combiSpaces) {
    numDOF += getCombiDegreesOfFreedom(space, boundary);
  }
  std::cout << "Combination scheme DOF : " << numDOF << " i.e. "
            << (static_cast<double>(numDOF * sizeof(CombiDataType)) / 1e9) << " GB " << std::endl;

  return numDOF;
}

inline long long int getSGDegreesOfFreedomFromDownSet(const std::vector<LevelVector>& downSet) {
  long long int numDOF = 0;
  for (const auto& subspaceLevel : downSet) {
    long long int numDOFSpace = 1;
    for (const auto& level_i : subspaceLevel) {
      assert(level_i > -1);
      // for the weird kind of boundary handling we currently have
      // (in the hierarchical spaces, the boundary points belong to level 1)
      assert(level_i > 0);
      if (level_i > 1) {
        numDOFSpace *= powerOfTwo[level_i - 1];
      } else {
        numDOFSpace *= 3;
      }
    }
    // std::cout << "Sparse grid subspace level : " << subspaceLevel << " has "
    //           << numDOFSpace << " DOF " << std::endl;
    numDOF += numDOFSpace;
  }
  return numDOF;
}

inline long long int printSGDegreesOfFreedomAdaptive(const LevelVector& lmin,
                                                     const LevelVector& lmax) {
  const auto dim = static_cast<DimType>(lmin.size());
  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  // combischeme.createDownSet();
  auto downSet = combischeme.getDownSet();
  auto numDOF = getSGDegreesOfFreedomFromDownSet(downSet);

  std::cout << "Sparse grid DOF : " << numDOF << " i.e. "
            << (static_cast<double>(numDOF * sizeof(CombiDataType)) / 1e9) << " GB " << std::endl;

  return numDOF;
}

inline IndexType getHighestIndexInHierarchicalSubspaceLowerThanNodalIndexOnLref(
    LevelType hierarchicalSubspaceLevel, IndexType nodalIndex, LevelType referenceLevel) {
  assert(hierarchicalSubspaceLevel <= referenceLevel);
  IndexType index = -1;
  if (hierarchicalSubspaceLevel == 0) {
    throw std::runtime_error(
        "this is basically untested, also if you use it pls remove next clause");
    if (nodalIndex > 0) {
      if (nodalIndex >= (powerOfTwo[referenceLevel] - 1)) {
        index = 1;
      } else {
        index = 0;
      }
    }
  } else if (hierarchicalSubspaceLevel == 1) {
    if (nodalIndex > 0) {
      if (nodalIndex >= (powerOfTwo[referenceLevel])) {
        index = 2;
      } else if (nodalIndex >= (powerOfTwo[referenceLevel] / 2)) {
        index = 1;
      } else {
        index = 0;
      }
    }
  } else {
    auto levelDiff = referenceLevel - hierarchicalSubspaceLevel;
    auto stepFactor = oneOverPowOfTwo[levelDiff];
    // if the nodal index' coordinate is larger than any of the subspace' indices' coordinates
    if (nodalIndex > powerOfTwo[levelDiff]) {
      index = static_cast<IndexType>(
          std::ceil(0.5 * (static_cast<double>(nodalIndex) * stepFactor - 1.)) - 1);
      assert(index < (nodalIndex * stepFactor - 1.) / 2.);
      assert((index + 1) >= (nodalIndex * stepFactor - 1.) / 2.);
    }
  }
  return index;
}

inline std::vector<long long int> getPartitionedNumDOFSG(
    std::vector<LevelVector> downSet, const LevelVector& referenceLevel,
    const std::vector<IndexVector>& decomposition) {
  if (downSet.size() == 0) {
    return {};
  }

  // this is only valid for with-boundary schemes!
  // cf downsampleDecomposition to extend to non-boundary
  auto dim = static_cast<DimType>(downSet[0].size());
  IndexVector decompositionOffsets;
  IndexType multiplier = 1;
  for (const auto& d : decomposition) {
    decompositionOffsets.push_back(multiplier);
    multiplier *= d.size();
  }
  auto numProcsPerGroup = multiplier;
  std::vector<long long int> numDOF(numProcsPerGroup, 0);

  // iterate subspaces
  for (const auto& subspaceLevel : downSet) {
    // assign decomposition just to set the right vector sizes
    std::vector<IndexVector> subspaceExtentsPerProcessPerDimension = decomposition;
    // iterate dimensions -> decomposition for subspace
    for (DimType d = 0; d < dim; ++d) {
      const auto& level_d = subspaceLevel[d];
      size_t numIndicesProcessedOnPoleMinusOne = -1;
      for (size_t k = 0; k < subspaceExtentsPerProcessPerDimension[d].size() - 1; ++k) {
        auto highestIndexOfSubspaceInThisDimensionInThisProcess =
            getHighestIndexInHierarchicalSubspaceLowerThanNodalIndexOnLref(
                level_d, decomposition[d][k + 1], referenceLevel[d]);
        if (highestIndexOfSubspaceInThisDimensionInThisProcess != -1) {
          subspaceExtentsPerProcessPerDimension[d][k] =
              highestIndexOfSubspaceInThisDimensionInThisProcess -
              numIndicesProcessedOnPoleMinusOne;
          numIndicesProcessedOnPoleMinusOne = highestIndexOfSubspaceInThisDimensionInThisProcess;
        } else {
          subspaceExtentsPerProcessPerDimension[d][k] = 0;
        }
      }
      // the rest belongs to the last partition
      subspaceExtentsPerProcessPerDimension[d][subspaceExtentsPerProcessPerDimension[d].size() -
                                               1] =
          powerOfTwo[level_d - 1] - numIndicesProcessedOnPoleMinusOne - 1;
      if (level_d == 1) {
        subspaceExtentsPerProcessPerDimension[d][subspaceExtentsPerProcessPerDimension[d].size() -
                                                 1] = 3 - numIndicesProcessedOnPoleMinusOne - 1;
        assert(std::accumulate(subspaceExtentsPerProcessPerDimension[d].begin(),
                               subspaceExtentsPerProcessPerDimension[d].end(), 0) == 3);
      } else {
        assert(std::accumulate(subspaceExtentsPerProcessPerDimension[d].begin(),
                               subspaceExtentsPerProcessPerDimension[d].end(),
                               0) == powerOfTwo[level_d - 1]);
      }
      // make sure last is not negative
      assert(subspaceExtentsPerProcessPerDimension[d].back() >= 0);
      // std::cout << subspaceLevel << " level_d " << level_d << " decomp " << decomposition[d] <<
      // subspaceExtentsPerProcessPerDimension[d] << std::endl;
    }
    // iterate all processes, add dof from subspace
    for (size_t i = 0; i < static_cast<size_t>(numProcsPerGroup); ++i) {
      size_t numDOFtoAdd = 1;
      // iterate the vector index entries belonging to linear index i
      auto tmp = i;
      for (int d = dim - 1; d >= 0; --d) {
        assert(d < subspaceExtentsPerProcessPerDimension.size());
        auto decompositionIndexInDimD = tmp / decompositionOffsets[d];
        assert(decompositionIndexInDimD < subspaceExtentsPerProcessPerDimension[d].size());
        numDOFtoAdd *= subspaceExtentsPerProcessPerDimension[d][decompositionIndexInDimD];
        tmp = tmp % decompositionOffsets[d];
      }
      numDOF[i] += numDOFtoAdd;
    }
  }
  // std::cout << "numDOF" << numDOF << std::endl;
  return numDOF;
}

inline std::vector<long long int> getPartitionedNumDOFSGAdaptive(
    LevelVector lmin, LevelVector lmax, const LevelVector& referenceLevel,
    const std::vector<IndexVector> decomposition) {
  assert((lmin.size() == lmax.size()) == (referenceLevel.size() == decomposition.size()));
  auto dim = static_cast<DimType>(lmin.size());
  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  // auto downSet = combischeme.getDownSet();
  // combischeme.createDownSet();
  auto downSet2 = combischeme.getDownSet();
  // for (const auto& s : downSet) {
  //   if (std::find(downSet2.begin(), downSet2.end(), s) == downSet2.end()) {
  //     std::cout << "in 1 but not 2 : " << s << std::endl;
  //   }
  // }
  // for (const auto& s : downSet2) {
  //   if (std::find(downSet.begin(), downSet.end(), s) == downSet.end()) {
  //     std::cout << "in 2 but not 1 : " << s << std::endl;
  //   }
  // }
  // std::cout << "Adaptive subspaces : " << downSet << std::endl;
  return getPartitionedNumDOFSG(downSet2, referenceLevel, decomposition);
}

inline std::vector<LevelVector> getConjointSet(const CombiMinMaxScheme& combischeme,
                                                      const LevelVector& lmin) {
  // we follow the idea in CombiMinMaxSchemeFromFile::createDownSet
  // but this time we count the occurrences of each grid in the individual downPoleSets
  // if a hierarchical level occurs d times, this means that this space is not required on the other
  // system (unless we have chosen a stupid decomposition -- this should maybe be checked)
  auto& combiSpaces = combischeme.getCombiSpaces();
  std::map<LevelVector, size_t> subspaces;
  // TODO actually, we need to check for the effective dimension
  //  get the length of level vector for now
  auto dim = combiSpaces[0].size();
  for (size_t i = 0; i < combiSpaces.size(); ++i) {
    // only check the active front of the scheme -- those grids with coefficient 1
    if (std::abs(combischeme.getCoeffs()[i] - 1.) < 1e-9) {
      // insert all subspaces "covered" by grid i into the map
      // iterate dimensions
      for (DimType d = 0; d < dim; ++d) {
        // iterate the downward pole in that dimension, leave out grid itself
        for (LevelType l_j = lmin[d]; l_j < combiSpaces[i][d]; ++l_j) {
          auto subspace = combiSpaces[i];
          subspace[d] = l_j;
          subspaces[subspace] += 1;
        }
      }
    }
  }
  // std::cout << "Conjoint subspace map : " << subspaces << std::endl;

  std::vector<LevelVector> conjointSet;
  for (const auto& item : subspaces) {
    if (item.second < dim) {
      conjointSet.push_back(item.first);
      // make the conjoint set downward closed
      auto thisGridsDownSet = combigrid::getDownSet(item.first);
      for (const auto& downSpace : thisGridsDownSet) {
        if (std::find(conjointSet.begin(), conjointSet.end(), downSpace) == conjointSet.end()) {
          conjointSet.push_back(downSpace);
        }
      }
    }
  }
  // std::cout << "Conjoint subspace level : " << conjointSet << std::endl;
  return conjointSet;
}

// for widely-distributed simulations, get the number of DOF that absolutely needs
// to be exchanged with the other system
inline long long int getNumDOFSGConjoint(
    const CombiMinMaxScheme& combischeme, const LevelVector& lmin) {
  auto conjointSet = getConjointSet(combischeme, lmin);
  return getSGDegreesOfFreedomFromDownSet(conjointSet);
}

// for widely-distributed simulations, get the number of DOF that absolutely needs
// to be exchanged with the other system -- partitioned
inline std::vector<long long int> getPartitionedNumDOFSGConjoint(
    const CombiMinMaxScheme& combischeme, const LevelVector& lmin, const LevelVector& referenceLevel,
    const std::vector<IndexVector> decomposition) {
  auto conjointSet = getConjointSet(combischeme, lmin);
  return getPartitionedNumDOFSG(conjointSet, referenceLevel, decomposition);
}

}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
