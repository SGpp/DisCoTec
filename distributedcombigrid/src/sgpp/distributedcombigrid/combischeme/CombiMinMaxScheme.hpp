#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <numeric>
#include "sgpp/distributedcombigrid/legacy/combigrid_utils.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"

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

  /* Creates the downset recursively */
  void createLevelsRec(DimType dim, LevelType n, DimType d, LevelVector& l,
                       const LevelVector& lmax);

  /* Calculate the coefficients of the classical CT (binomial coefficient)*/
  void computeCombiCoeffsClassical();

  /* Calculate the coefficients of the adaptive CT using the formula in Alfredo's
   * SDC paper (from Brendan Harding)*/
  void computeCombiCoeffsAdaptive();

  LevelVector getLevelMinima();
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
            lvl[i] = l.second.get_value<int>();
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

static long long int printCombiDegreesOfFreedom(const std::vector<LevelVector>& combiSpaces) {
  long long int numDOF = 0;
  for (const auto& space : combiSpaces) {
    long int numDOFSpace = 1;
    for (const auto& level_i : space) {
      assert(level_i > -1);
      if (level_i > 0) {
        numDOFSpace *= powerOfTwo[level_i] + 1;
      } else {
        numDOFSpace *= 2;
      }
    }
    numDOF += numDOFSpace;
  }
  std::cout << "Combination scheme DOF : " << numDOF << " i.e. "
            << (static_cast<double>(numDOF * sizeof(CombiDataType)) / 1e9) << " GB " << std::endl;

  return numDOF;
}

static long long int printSGDegreesOfFreedomAdaptive(const LevelVector& lmin,
                                                     const LevelVector& lmax) {
  DimType dim = lmin.size();
  long long int numDOF = 0;
  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  auto downSet = combischeme.getDownSet();
  for (const auto& subspaceLevel : downSet) {
    long int numDOFSpace = 1;
    for (const auto& level_i : subspaceLevel) {
      assert(level_i > -1);
      // for the weird kind of boundary handling we currently have (the boundary points belong to
      // level 1)
      assert(level_i > 0);
      if (level_i > 1) {
        numDOFSpace *= powerOfTwo[level_i - 1];
      } else {
        numDOFSpace *= 3;
      }
    }
    numDOF += numDOFSpace;
  }

  std::cout << "Sparse grid DOF : " << numDOF << " i.e. "
            << (static_cast<double>(numDOF * sizeof(CombiDataType)) / 1e9) << " GB " << std::endl;

  return numDOF;
}

static IndexType getHighestIndexInHierarchicalSubspaceLowerThanNodalIndexOnLref(
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

// TODO add test for this -- to compare to actual sg assigned sizes
static std::vector<long long int> getPartitionedNumDOFSGAdaptive(
    LevelVector lmin, LevelVector lmax, const LevelVector& referenceLevel,
    const std::vector<IndexVector> decomposition) {
  // this is only valid for with-boundary schemes!
  // cf downsampleDecomposition to extend to non-boundary
  assert((lmin.size() == lmax.size()) == (referenceLevel.size() == decomposition.size()));
  DimType dim = lmin.size();
  IndexVector decompositionOffsets;
  IndexType multiplier = 1;
  for (const auto& d : decomposition) {
    decompositionOffsets.push_back(multiplier);
    multiplier *= d.size();
  }
  auto numProcsPerGroup = multiplier;
  std::vector<long long int> numDOF(numProcsPerGroup, 0);

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createAdaptiveCombischeme();
  auto downSet = combischeme.getDownSet();
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
    for (size_t i = 0; i < numProcsPerGroup; ++i) {
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

}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
