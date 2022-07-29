
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {
std::string toString(combigrid::LevelVector const& l) {
  std::stringstream ss;
  for (size_t i = 0; i < l.size(); ++i) {
    if (i != 0) ss << ",";
    ss << l[i];
  }
  return ss.str();
}

void getDownSetRecursively(combigrid::LevelVector const& l, combigrid::LevelVector fixedDimensions,
                           std::vector<LevelVector>& downSet) {
  if (fixedDimensions.size() == l.size()) {
    downSet.push_back(fixedDimensions);
  } else {
    // for the whole range of values in dimension d
    size_t d = fixedDimensions.size();
    for (LevelType i = 1; i <= l[d]; ++i) {  // levels start at 1 here, in compliance with our
                                         // definition of DistributedSparseGrid
      auto moreDimensionsFixed = fixedDimensions;
      moreDimensionsFixed.push_back(i);
      getDownSetRecursively(l, moreDimensionsFixed, downSet);
    }
  }
}

std::vector<LevelVector> getDownSet(combigrid::LevelVector const& l) {
  std::vector<LevelVector> downSet;
  getDownSetRecursively(l, LevelVector(), downSet);
  return downSet;
}

// cf.
// https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
std::vector<std::vector<DimType>> getAllKOutOfDDimensions(DimType k, DimType d) {
  std::vector<std::vector<DimType>> combinations;
  std::vector<DimType> selected;
  std::vector<DimType> selector(d);
  std::fill(selector.begin(), selector.begin() + k, 1);
  do {
    for (DimType i = 0; i < d; i++) {
      if (selector[i]) {
        selected.push_back(i);
      }
    }
    combinations.push_back(selected);
    selected.clear();
  } while (std::prev_permutation(selector.begin(), selector.end()));
  return combinations;
}

void createTruncatedHierarchicalLevelsRec(size_t dim, size_t n, LevelVector& l,
                                        const LevelVector& lmax, const LevelVector& lmin,
                                        std::vector<LevelVector>& created) {
  assert(lmax[dim-1] == lmin[dim-1] + n);

  // sum rightmost entries of level vector
  LevelType lsum(0);
  for (size_t i = dim; i < l.size(); ++i) {
    lsum += l[i];
  }

  // iterate everything below hyperplane
  for (LevelType ldim = 1; ldim <= LevelType(sum(lmin) + n) - lsum; ++ldim) {
    // smallereq than lmax in every dim
    if (ldim > lmax[dim - 1]) {
      continue;
    } else {
      l[dim - 1] = ldim;
      if (dim == 1) {
        // all mixed dimension sums
        bool pleaseAdd = true;
        DimType d = lmin.size();
        for (DimType k = 2; k <= d; ++k) {
          auto dimList = getAllKOutOfDDimensions(k, d);
          // for each subselection of dimensions, compute sum of l and lmin
          for (const auto& dimCombination : dimList) {
            LevelType partlsum(0), partlminsum(0);
            for (const auto& i : dimCombination) {
              partlsum += l[i];
              partlminsum += lmin[i];
            }
            if ((partlsum > static_cast<LevelType>(partlminsum + n))) {
              // std::cout << k << " k " << l << partlsum << " " << dimCombination << " other stop "
              // << partlminsum << " n " << n << std::endl;
              pleaseAdd = false;
              break;
            }
          }
        }
        if (pleaseAdd) {
          created.push_back(l);
        }
      } else {
        createTruncatedHierarchicalLevelsRec(dim - 1, n, l, lmax, lmin, created);
      }
    }
  }
}

void createTruncatedHierarchicalLevels(const LevelVector& lmax,
                                     const LevelVector& lmin, std::vector<LevelVector>& created) {
  auto dim = lmax.size();
  assert(lmin.size() == dim);

  LevelVector ldiff = lmax - lmin;
  LevelType minLevelDifference = *(std::min_element(ldiff.begin(), ldiff.end()));

  //TODO currently, we basically enlarge the sparse grid
  // (potentially more hierarchical subspaces than necessary for scheme!)
  // it is fine for evenly spaced level differences though
  LevelVector rlmin(dim);

  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = lmax[i] - minLevelDifference;
  }

  LevelType n = minLevelDifference;

  LevelVector l(dim);
  createTruncatedHierarchicalLevelsRec(dim, n, l, lmax, rlmin, created);
}

}  // namespace combigrid
