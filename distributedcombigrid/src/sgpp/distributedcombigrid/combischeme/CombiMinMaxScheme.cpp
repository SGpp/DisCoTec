#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

namespace combigrid {

void CombiMinMaxScheme::createClassicalCombischeme() {
  // Remove dummy dimensions (e.g. if lmin = (2,2,2) and lmax = (4,4,2), dimension 3 is dummy
  // and effDim_ = 2)
  LevelVector lmaxtmp = lmax_;
  LevelVector lmintmp = lmin_;
  LevelVector dummyDims;
  for (size_t i = 0; i < lmax_.size(); ++i) {
    if (lmax_[i] - lmin_[i] == 0) {
      lmaxtmp.erase(lmaxtmp.begin() + i - dummyDims.size());
      lmintmp.erase(lmintmp.begin() + i - dummyDims.size());
      dummyDims.push_back(i);
    }
  }

  // compute c which fulfills lmax - c*1  >= lmin
  LevelType c(0);
  LevelVector cv = lmaxtmp - lmintmp;

  // check if all elements in lmax-lmin are equal. If not, define c as the max of cv
  if (std::adjacent_find(cv.begin(), cv.end(), std::not_equal_to<int>()) == cv.end()) {
    c = cv[effDim_ - 1];
  } else {
    c = *std::min_element(cv.begin(), cv.end());
  }
  LevelVector rlmin(effDim_);
  for (size_t i = 0; i < rlmin.size(); ++i) {
    rlmin[i] = lmaxtmp[i] - c;
  }

  // Define new lmin (dummy dimensions will be fixed later)
  for (size_t i = 0; i < lmin_.size(); ++i) {
    lmin_[i] = lmax_[i] - c;
  }
  LevelType n = sum(rlmin) + c;

  LevelVector l(effDim_, 0);
  createLevelsRec(effDim_, n, effDim_, l, lmaxtmp);

  // re-insert dummy dimensions left out
  for (size_t i = 0; i < dummyDims.size(); ++i) {
    for (size_t j = 0; j < levels_.size(); ++j) {
      levels_[j].insert(levels_[j].begin() + dummyDims[i], lmax_[i]);
    }
    lmin_[dummyDims[i]] = lmax_[i];
  }
  n_ = sum(lmin_) + c;
  // create combi spaces
  for (size_t i = 0; i < levels_.size(); ++i) {
    LevelVector& l = levels_[i];
    for (LevelType p = 0; p < LevelType(effDim_); ++p) {
      if (l >= lmin_ && sum(l) == n_ - p) combiSpaces_.push_back(l);
    }
  }
  computeCombiCoeffsClassical();
}

void CombiMinMaxScheme::createAdaptiveCombischeme() {
  LevelVector lm = getLevelMinima();
  n_ = std::accumulate(lm.begin(), lm.end(), 0);
  LevelVector l(dim_);
  // TODO: Why levelSum-1 ?
  createLevelsRec(dim_, n_ - 1, dim_, l, lmax_);

  for (auto level : levels_) {
    int l1norm = std::accumulate(level.begin(), level.end(), 0);
    if (level >= lmin_ && l1norm <= n_) combiSpaces_.push_back(level);
  }

  computeCombiCoeffsAdaptive();
}

/*
 * Generates the fault tolerant combination technique with extra
 * grids used in case of faults
 * */
void CombiMinMaxScheme::makeFaultTolerant() {
  const int extraDiags = 2;
  if (lmin_ == lmax_) return;
  // Add extra combiSpaces to ensure fault tolerance
  for (size_t i = 0; i < levels_.size(); ++i) {
    LevelVector& l = levels_[i];
    for (LevelType p = LevelType(effDim_); p < LevelType(effDim_) + extraDiags; ++p) {
      if (l >= lmin_ && sum(l) == n_ - p) {
        combiSpaces_.push_back(l);
        coefficients_.push_back(0.0);
      }
    }
  }
}

LevelVector CombiMinMaxScheme::getLevelMinima() {
  LevelVector tmp = lmax_ - lmin_;
  std::vector<IndexType>::iterator result = std::max_element(tmp.begin(), tmp.end());
  LevelType indexMax = std::distance(tmp.begin(), result);
  LevelVector lm(lmin_);
  lm[indexMax] = lmax_[indexMax];
  return lm;
}

void CombiMinMaxScheme::createLevelsRec(DimType dim, LevelType n, DimType d, LevelVector& l,
                                        const LevelVector& lmax) {
  // sum rightmost entries of level vector
  LevelType lsum(0);
  for (size_t i = dim; i < l.size(); ++i) lsum += l[i];

  for (LevelType ldim = 1; ldim <= LevelType(n) + LevelType(d) - 1 - lsum; ++ldim) {
    l[dim - 1] = ldim;
    if (dim == 1) {
      if (l <= lmax) {
        levels_.push_back(l);
      }
    } else {
      createLevelsRec(dim - 1, n, d, l, lmax);
    }
  }
}

void CombiMinMaxScheme::computeCombiCoeffsAdaptive() {
  for (size_t i = 0; i < combiSpaces_.size(); i++) {
    real coeff = 0;
    LevelVector tmp(combiSpaces_[i]);
    tmp = tmp + LevelVector(dim_, 1);
    for (auto nbr : combiSpaces_) {
      if (nbr >= combiSpaces_[i] && nbr <= tmp) {
        LevelVector diff = nbr - combiSpaces_[i];
        coeff += std::pow(-1.0, std::accumulate(diff.begin(), diff.end(), 0));
      }
    }
    // Flag combiSpaces that will be erased later
    if (coeff == 0) combiSpaces_[i] = LevelVector(dim_, 0);
    coefficients_.push_back(coeff);
  }
  // Delete combiSpaces_ and coefficients_ for which the coeff is zero
  combiSpaces_.erase(std::remove(combiSpaces_.begin(), combiSpaces_.end(), LevelVector(dim_, 0)),
                     combiSpaces_.end());
  coefficients_.erase(std::remove(coefficients_.begin(), coefficients_.end(), 0),
                      coefficients_.end());
}

void CombiMinMaxScheme::computeCombiCoeffsClassical() {
  for (DimType i = 0; i < combiSpaces_.size(); ++i) {
    LevelType l1 = sum(combiSpaces_[i]);
    LevelType p = n_ - l1;
    // Classical combination coefficients
    coefficients_.push_back(std::pow(-1, p) *
                            boost::math::binomial_coefficient<real>(effDim_ - 1, p));
  }
}

}  // end namespace combigrid
