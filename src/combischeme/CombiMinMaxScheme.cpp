#include "combischeme/CombiMinMaxScheme.hpp"

namespace combigrid {

void CombiMinMaxScheme::createClassicalCombischeme() {
  // Remove dummy dimensions (e.g. if lmin = (2,2,2) and lmax = (4,4,2), dimension 3 is dummy
  // and effDim_ = 2)
  LevelVector lmaxtmp = lmax_;
  LevelVector lmintmp = lmin_;
  std::vector<DimType> dummyDims;
  for (DimType i = 0; i < static_cast<DimType>(lmax_.size()); ++i) {
    if (lmax_[i] - lmin_[i] == 0) {
      lmaxtmp.erase(lmaxtmp.begin() + i - dummyDims.size());
      lmintmp.erase(lmintmp.begin() + i - dummyDims.size());
      dummyDims.push_back(i);
    }
  }

  // compute c which fulfills lmax - c*1  >= lmin
  LevelType c(0);
  LevelVector cv = lmaxtmp - lmintmp;

  // check if all elements in lmax-lmin are equal. If not, define c as the min of cv
  if (std::adjacent_find(cv.begin(), cv.end(), std::not_equal_to<int>()) == cv.end()) {
    c = cv[effDim_ - 1];
  } else {
    c = *std::min_element(cv.begin(), cv.end());
  }
  for (DimType i = 0; i < static_cast<DimType>(lmin_.size()); ++i) {
    lmin_[i] = static_cast<LevelType>(lmax_[i] - c);
  }

  combigrid::createTruncatedHierarchicalLevels(lmaxtmp, lmintmp, levels_);

  // re-insert dummy dimensions left out
  for (size_t i = 0; i < dummyDims.size(); ++i) {
    for (size_t j = 0; j < levels_.size(); ++j) {
      levels_[j].insert(levels_[j].begin() + dummyDims[i], lmax_[i]);
    }
    lmin_[dummyDims[i]] = lmax_[i];
  }
  n_ = static_cast<LevelType>(combigrid::levelSum(lmin_) + c);
  // create combi spaces
  for (size_t i = 0; i < levels_.size(); ++i) {
    LevelVector& l = levels_[i];
    for (LevelType p = 0; p < LevelType(effDim_); ++p) {
      if (l >= lmin_ && combigrid::levelSum(l) == n_ - p) combiSpaces_.push_back(l);
    }
  }
  computeCombiCoeffsClassical();

  if (levels_.empty()) {
    throw std::runtime_error("No combination levels!");
  }
}

LevelVector getFurthestCorner(LevelVector& lmax, LevelVector& lmin) {
  LevelVector ldiff = lmax - lmin;
  LevelVector::iterator result = std::max_element(ldiff.begin(), ldiff.end());
  LevelType indexMax = static_cast<LevelType>(std::distance(ldiff.begin(), result));
  LevelVector lm(lmin);
  lm[indexMax] = lmax[indexMax];
  return lm;
}

void CombiMinMaxScheme::createAdaptiveCombischeme() {
  LevelVector lm = getFurthestCorner(lmax_, lmin_);
  n_ = std::accumulate(lm.begin(), lm.end(), static_cast<LevelType>(0)); // = sum(lmin_) + max(ldiff)
  LevelVector l(dim_);
  // // TODO: Why levelSum-1 ?
  // createLevelsRec(dim_, n_ - 1, dim_, l, lmax_);
  auto n = static_cast<LevelType>(n_ - combigrid::levelSum(lmin_));
  auto rlmin = lmax_;
  for (auto& rl: rlmin) {
    rl = static_cast<LevelType>(rl - n);
  }
  combigrid::createTruncatedHierarchicalLevels(lmax_, rlmin, levels_);

  for (auto level : levels_) {
    int l1norm = std::accumulate(level.begin(), level.end(), 0);
    if (level >= lmin_ && l1norm <= n_) combiSpaces_.push_back(level);
  }

  computeCombiCoeffsAdaptive();
  if (levels_.empty()) {
    throw std::runtime_error("No combination levels!");
  }
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
      if (l >= lmin_ && combigrid::levelSum(l) == n_ - p) {
        combiSpaces_.push_back(l);
        coefficients_.push_back(0.0);
      }
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
    LevelType l1 = combigrid::levelSum(combiSpaces_[i]);
    auto p = static_cast<LevelType>(n_ - l1);
    // Classical combination coefficients
    coefficients_.push_back(std::pow(-1, p) *
                            boost::math::binomial_coefficient<real>(static_cast<unsigned int>(effDim_ - 1), static_cast<unsigned int>(p)));
  }
}

}  // end namespace combigrid
