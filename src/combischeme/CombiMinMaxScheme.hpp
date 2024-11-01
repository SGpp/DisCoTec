#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <numeric>

#include "io/BroadcastParameters.hpp"
#include "utils/LevelSetUtils.hpp"
#include "utils/LevelVector.hpp"
#include "utils/PowerOfTwo.hpp"
#include "utils/Types.hpp"

namespace combigrid {
/**
 * @class CombiMinMaxScheme
 *
 * @brief Class to generate combination schemes based on a maximal and minimal level
 */
class CombiMinMaxScheme {
 public:
  /**
   * @brief Constructor
   *
   * @param dim dimensionality of the PDE problem domain
   * @param lmin minimum level of the combination technique, length dim
   * @param lmax maximum level of the combination technique, length dim
   */
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

  /**
   * @brief Constructor
   *
   * for externally-generated combination schemes
   *
   * @param levels level for each component grid
   * @param coefficients coefficient for each component grid
   */
  CombiMinMaxScheme(const std::vector<LevelVector>& levels, const std::vector<real>& coefficients)
      : combiSpaces_(levels), coefficients_(coefficients) {
    assert(levels.size() > 0);
    assert(levels.size() == coefficients.size());
    n_ = 0;
    dim_ = static_cast<combigrid::DimType>(levels.back().size());
  }

  virtual ~CombiMinMaxScheme() = default;

  /**
   * @brief Generate the combischeme corresponding to the classical combination technique.
   *
   * This sets the combiSpaces_ and coefficients_ member variables.
   */
  void createClassicalCombischeme();

  /**
   * @brief Generates an adaptive combination scheme (equivalent to CK's
   * Python code)
   *
   * This sets the combiSpaces_ and coefficients_ member variables.
   */
  void createAdaptiveCombischeme();

  /**
   * @brief Updates the combination scheme to be fault tolerant
   *
   * For fault tolerance, "cheap" grids with an initial coefficient of 0. are added to the scheme.
   * This extends the combiSpaces_ and coefficients_ member variables.
   */
  void makeFaultTolerant();

  /**
   * @brief Get the levels of the component grids
   */
  inline const std::vector<LevelVector>& getCombiSpaces() const { return combiSpaces_; }

  /**
   * @brief Get the coefficients of the component grids
   */
  inline const std::vector<double>& getCoeffs() const { return coefficients_; }

  /**
   * @brief Get the downward closed set of the CombiMinMaxScheme
   * 
   * i.e. all hierarchical subspaces / mixed resolutions that will be covered by the scheme
   */
  inline const std::vector<LevelVector>& getDownSet() { return levels_; }

  /**
   * @brief (re-)generate the downward closed set of the CombiScheme
   *
   * The downward closed set corresponds to the subspaces in the DistributedSparseGridUniform class.
   * (Re-)sets the levels_ data member.
   */
  void createDownSet() {
    std::set<LevelVector> subspaces;
    for (size_t i = 0; i < combiSpaces_.size(); ++i) {
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
  CombiMinMaxSchemeFromFile(DimType dim, const LevelVector& lmin, const LevelVector& lmax,
                            std::string ctschemeFile, const CommunicatorType& comm = MPI_COMM_WORLD)
      : CombiMinMaxScheme(dim, lmin, lmax) {
    // read in CT scheme, if applicable
    boost::property_tree::ptree pScheme =
        broadcastParameters::getJsonFromRankZero(ctschemeFile, comm);
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
          processGroupNumbers_.push_back(c.second.get_value<RankType>());
        } else {
          assert(false);
        }
      }
      // process group numbers should either be always given or never
      assert(coefficients_.size() == combiSpaces_.size());
      assert(processGroupNumbers_.size() == combiSpaces_.size() ||
             processGroupNumbers_.size() == 0);
    }
    assert(coefficients_.size() > 0);
    assert(coefficients_.size() == combiSpaces_.size());
  }

  inline const std::vector<RankType>& getProcessGroupNumbers() const {
    return processGroupNumbers_;
  }

 private:
  std::vector<RankType> processGroupNumbers_;
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
    multiplier *= static_cast<IndexType>(d.size());
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
              static_cast<IndexType>(highestIndexOfSubspaceInThisDimensionInThisProcess -
                                     numIndicesProcessedOnPoleMinusOne);
          numIndicesProcessedOnPoleMinusOne = highestIndexOfSubspaceInThisDimensionInThisProcess;
        } else {
          subspaceExtentsPerProcessPerDimension[d][k] = 0;
        }
      }
      // the rest belongs to the last partition
      subspaceExtentsPerProcessPerDimension[d][subspaceExtentsPerProcessPerDimension[d].size() -
                                               1] =
          powerOfTwo[level_d - 1] - static_cast<IndexType>(numIndicesProcessedOnPoleMinusOne) - 1;
      if (level_d == 1) {
        subspaceExtentsPerProcessPerDimension[d][subspaceExtentsPerProcessPerDimension[d].size() -
                                                 1] =
            static_cast<IndexType>(3 - numIndicesProcessedOnPoleMinusOne - 1);
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
      for (IndexType d = dim - 1; d >= 0; --d) {
        assert(static_cast<size_t>(d) < subspaceExtentsPerProcessPerDimension.size());
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
    const std::vector<IndexVector>& decomposition) {
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
inline long long int getNumDOFSGConjoint(const CombiMinMaxScheme& combischeme,
                                         const LevelVector& lmin) {
  auto conjointSet = getConjointSet(combischeme, lmin);
  return getSGDegreesOfFreedomFromDownSet(conjointSet);
}

// for widely-distributed simulations, get the number of DOF that absolutely needs
// to be exchanged with the other system -- partitioned
inline std::vector<long long int> getPartitionedNumDOFSGConjoint(
    const CombiMinMaxScheme& combischeme, const LevelVector& lmin,
    const LevelVector& referenceLevel, const std::vector<IndexVector>& decomposition) {
  auto conjointSet = getConjointSet(combischeme, lmin);
  return getPartitionedNumDOFSG(conjointSet, referenceLevel, decomposition);
}

inline size_t getAssignedLevels(const CombiMinMaxSchemeFromFile& combischeme,
                                RankType myProcessGroupNumber, std::vector<LevelVector>& levels,
                                std::vector<combigrid::real>& coeffs,
                                std::vector<size_t>& taskNumbers) {
  assert(levels.empty());
  assert(levels.size() == coeffs.size());
  assert(levels.size() == taskNumbers.size());

  const auto& pgNumbers = combischeme.getProcessGroupNumbers();
  const auto& allCoeffs = combischeme.getCoeffs();
  const auto& allLevels = combischeme.getCombiSpaces();
  assert(allLevels.size() == allCoeffs.size());
  assert(pgNumbers.size() == allCoeffs.size());
  if (pgNumbers.empty() || pgNumbers.size() != allCoeffs.size() ||
      pgNumbers.size() != allLevels.size()) {
    throw std::runtime_error(
        "CombiSchemeFromFile::getAssignedLevels : inconsistent scheme from file");
  }
  [[maybe_unused]] const auto [itMin, itMax] =
      std::minmax_element(pgNumbers.begin(), pgNumbers.end());
  assert(*itMin == 0);  // make sure it starts with 0
  // filter out only those tasks that belong to "our" process group
  for (size_t taskNo = 0; taskNo < pgNumbers.size(); ++taskNo) {
    if (pgNumbers[taskNo] == myProcessGroupNumber) {
      taskNumbers.push_back(taskNo);
      coeffs.push_back(allCoeffs[taskNo]);
      levels.push_back(allLevels[taskNo]);
    }
  }
  return allLevels.size();
}

inline size_t getRoundRobinLevels(const CombiMinMaxScheme& combischeme,
                                  RankType myProcessGroupNumber, size_t totalNumGroups,
                                  std::vector<LevelVector>& levels,
                                  std::vector<combigrid::real>& coeffs,
                                  std::vector<size_t>& taskNumbers) {
  assert(levels.empty());
  assert(levels.size() == coeffs.size());
  assert(levels.size() == taskNumbers.size());

  const auto& allCoeffs = combischeme.getCoeffs();
  const auto& allLevels = combischeme.getCombiSpaces();
  assert(allLevels.size() == allCoeffs.size());
  // make sure its enough tasks
  assert(allLevels.size() >= totalNumGroups);

  for (size_t taskNo = myProcessGroupNumber; taskNo < allLevels.size(); taskNo += totalNumGroups) {
    taskNumbers.push_back(taskNo);
    coeffs.push_back(allCoeffs[taskNo]);
    levels.push_back(allLevels[taskNo]);
  }
  return allLevels.size();
}

inline size_t getLoadBalancedLevels(const CombiMinMaxScheme& combischeme,
                                    RankType myProcessGroupNumber, size_t totalNumGroups,
                                    const std::vector<BoundaryType>& boundary,
                                    std::vector<LevelVector>& levels,
                                    std::vector<combigrid::real>& coeffs,
                                    std::vector<size_t>& taskNumbers) {
  assert(levels.empty());
  assert(levels.size() == coeffs.size());
  assert(levels.size() == taskNumbers.size());

  const auto& allCoeffs = combischeme.getCoeffs();
  const auto& allLevels = combischeme.getCombiSpaces();
  std::vector<size_t> allTaskIDs(allLevels.size());
  std::iota(allTaskIDs.begin(), allTaskIDs.end(), 0);
  assert(allLevels.size() == allCoeffs.size());
  // make sure its enough tasks
  assert(allLevels.size() >= totalNumGroups);

  // sort the tasks by their load
  // TODO use LoadModel class in this function
  std::sort(allTaskIDs.begin(), allTaskIDs.end(), [&allLevels, &boundary](size_t i1, size_t i2) {
    return getCombiDegreesOfFreedom(allLevels[i1], boundary) >
           getCombiDegreesOfFreedom(allLevels[i2], boundary);
  });

  std::map<RankType, long long> loadAssignedToGroup;
  for (size_t taskNo : allTaskIDs) {
    // find the group with the lowest load
    RankType minGroup = 0;
    long long minLoad = std::numeric_limits<long long>::max();
    for (RankType group = 0; group < static_cast<int>(totalNumGroups); ++group) {
      if (loadAssignedToGroup[group] < minLoad) {
        minGroup = group;
        minLoad = loadAssignedToGroup[group];
      }
    }
    loadAssignedToGroup[minGroup] += getCombiDegreesOfFreedom(allLevels[taskNo], boundary);
    if (minGroup == myProcessGroupNumber) {
      taskNumbers.push_back(taskNo);
      coeffs.push_back(allCoeffs[taskNo]);
      levels.push_back(allLevels[taskNo]);
    }
  }
  return allLevels.size();
}

}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
