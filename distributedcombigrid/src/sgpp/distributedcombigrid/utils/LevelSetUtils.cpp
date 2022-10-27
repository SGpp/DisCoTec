#include "sgpp/distributedcombigrid/utils/LevelSetUtils.hpp"

#include <map>

namespace combigrid {
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
  assert(k <= d);
  assert(d > 0);
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

const std::vector<std::vector<DimType>>& AllKOutOfDDimensions::get(DimType k, DimType d) {
  auto key = std::make_pair(k, d);
  if (cache_.find(key) == cache_.end()) {
    cache_[key] = getAllKOutOfDDimensions(k, d);
  }
  return cache_[key];
}
std::map<std::pair<DimType, DimType>, std::vector<std::vector<DimType>>>
    AllKOutOfDDimensions::cache_;

void createTruncatedHierarchicalLevelsRec(DimType dim, size_t n, LevelVector& l,
                                          const LevelVector& lmax, const LevelVector& lmin,
                                          std::vector<LevelVector>& created) {
  assert(lmax.size() == lmin.size());
  auto dimensionality = static_cast<DimType>(lmax.size());
  assert(lmax[dim - 1] == lmin[dim - 1] + n);

  // sum rightmost entries of level vector
  LevelType lsum(0);
  for (size_t i = dim; i < l.size(); ++i) {
    lsum = static_cast<LevelType>(lsum + l[i]);
  }

  // iterate everything below hyperplane
  for (LevelType ldim = 1; ldim <= LevelType(combigrid::levelSum(lmin) + n) - lsum; ++ldim) {
    // smallereq than lmax in every dim
    if (ldim > lmax[dim - 1]) {
      continue;
    } else {
      l[dim - 1] = ldim;
      if (dim == 1) {
        // all mixed dimension sums
        bool pleaseAdd = true;
        for (DimType k = 2; k <= dimensionality; ++k) {
          auto dimList = AllKOutOfDDimensions::get(k, dimensionality);
          // auto dimList = getAllKOutOfDDimensions(k, dimensionality);
          // for each subselection of dimensions, compute sum of l and lmin
          for (const auto& dimCombination : dimList) {
            LevelType partlsum(0), partlminsum(0);
            for (const auto& i : dimCombination) {
              partlsum = static_cast<LevelType>(partlsum + l[i]);
              partlminsum = static_cast<LevelType>(partlminsum + lmin[i]);
            }
            if ((partlsum > static_cast<LevelType>(partlminsum + n))) {
              // std::cout << k << " k " << l << partlsum << " " << dimCombination << " other stop "
              // << partlminsum << " n " << n << std::endl;
              pleaseAdd = false;
              break;
            }
          }
        }
        if (pleaseAdd == true) {
          created.push_back(l);
        }
      } else {
        createTruncatedHierarchicalLevelsRec(dim - 1, n, l, lmax, lmin, created);
      }
    }
  }
}

void createTruncatedHierarchicalLevels(const LevelVector& lmax, const LevelVector& lmin,
                                       std::vector<LevelVector>& created) {
  assert(created.empty());
  auto dim = static_cast<DimType>(lmax.size());
  assert(lmin.size() == dim);

  for (DimType d = 0; d < dim; ++d) {
    assert(lmin[d] <= lmax[d]);
  }

  LevelVector ldiff = lmax - lmin;
  LevelType minLevelDifference = *(std::min_element(ldiff.begin(), ldiff.end()));

  // TODO currently, we basically enlarge the sparse grid
  //  (potentially more hierarchical subspaces than necessary for scheme!)
  //  it is fine for evenly spaced level differences though
  LevelVector rlmin(dim);

  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = static_cast<LevelType>(lmax[i] - minLevelDifference);
  }

  LevelType n = minLevelDifference;

  LevelVector l(dim);
  createTruncatedHierarchicalLevelsRec(dim, n, l, lmax, rlmin, created);
}

}  // namespace combigrid
