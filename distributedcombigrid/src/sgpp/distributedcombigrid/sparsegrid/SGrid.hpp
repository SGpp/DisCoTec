#ifndef SGRID_HPP_
#define SGRID_HPP_

#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>  // std::rand, std::srand
#include <ctime>    // std::time
#include <iostream>
#include <vector>

#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/LevelSetUtils.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

template <typename FG_ELEMENT>
class SGrid {
 public:
  /** create sparse grid of dimension d and level n with or without boundary
   *  in all dimensions
   */
  SGrid(DimType dim, LevelType n, bool boundary = false);

  /** create sparse grid of dimension d and specify for each dimension the
   *  the maximum discretization level and whether there is a boundary or not
   */
  SGrid(DimType dim, const LevelVector& nmax, const LevelVector& lmin,
        const std::vector<BoundaryType>& boundary);

  virtual ~SGrid();

  void print(std::ostream& os) const;

  inline size_t getSize() const;

  size_t getLevelIndex(const LevelVector& l) const;

  inline bool isContained(const LevelVector& l, size_t& index) const;

  inline bool isContained(const LevelVector& l) const;

  inline DimType getDim() const;

  inline const LevelVector& getNMax() const;

  inline const LevelVector& getNMin() const;

  inline const LevelVector& getLevelVector(size_t i) const;

  inline const std::vector<BoundaryType>& getBoundaryVector() const;

 private:
  void createLevels(DimType dim, const LevelVector& nmax, const LevelVector& lmin);

  DimType dim_;

  LevelVector nmax_;

  LevelVector lmin_;

  std::vector<LevelVector> levels_;

  std::vector<BoundaryType> boundary_;
};

}  // namespace

namespace combigrid {

template <typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const SGrid<FG_ELEMENT>& sg);

template <typename FG_ELEMENT>
inline size_t SGrid<FG_ELEMENT>::getSize() const {
  return levels_.size();
}


template <typename FG_ELEMENT>
inline DimType SGrid<FG_ELEMENT>::getDim() const {
  return dim_;
}

template <typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getNMax() const {
  return nmax_;
}

template <typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getNMin() const {
  return lmin_;
}

template <typename FG_ELEMENT>
inline const LevelVector& SGrid<FG_ELEMENT>::getLevelVector(size_t i) const {
  return levels_[i];
}

// at construction create only levels, no data
template <typename FG_ELEMENT>
SGrid<FG_ELEMENT>::SGrid(DimType dim, LevelType n, bool boundary)
    : dim_(dim) {
  assert(dim > 0);
  assert(n > 0);

  nmax_.resize(dim, n);
  lmin_.resize(dim, 1);
  boundary_.resize(dim, boundary);

  createLevels(dim, nmax_, lmin_);
}

// at construction create only levels, no data
template <typename FG_ELEMENT>
SGrid<FG_ELEMENT>::SGrid(DimType dim, const LevelVector& nmax, const LevelVector& lmin,
                         const std::vector<BoundaryType>& boundary)
    : dim_(dim) {
  assert(dim > 0);

  assert(nmax.size() == dim);

  for (size_t i = 0; i < nmax.size(); ++i) assert(nmax[i] > 0);

  assert(lmin.size() == dim);

  for (size_t i = 0; i < lmin.size(); ++i) assert(lmin[i] > 0);

  assert(boundary.size() == dim);

  nmax_ = nmax;
  lmin_ = lmin;
  boundary_ = boundary;

  createLevels(dim, nmax_, lmin_);
}

template <typename FG_ELEMENT>
std::ostream& operator<<(std::ostream& os, const SGrid<FG_ELEMENT>& sg) {
  sg.print(os);
  return os;
}


template <typename FG_ELEMENT>
SGrid<FG_ELEMENT>::~SGrid() = default;

template <typename FG_ELEMENT>
void SGrid<FG_ELEMENT>::createLevels(DimType dim, const LevelVector& nmax,
                                     const LevelVector& lmin) {
  combigrid::createTruncatedHierarchicalLevels(nmax, lmin, levels_);
}

template <typename FG_ELEMENT>
inline const std::vector<BoundaryType>& SGrid<FG_ELEMENT>::getBoundaryVector() const {
  return boundary_;
}


} /* namespace combigrid */

#endif /* SGRID_HPP_ */
