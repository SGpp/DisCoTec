/*
 * CombiMinMaxScheme.hpp
 *
 *  Created on: Oct 2, 2015
 *      Author: sccs
 */

#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <numeric>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/SGrid.hpp"

namespace combigrid {

class CombiMinMaxScheme {
 public:
  CombiMinMaxScheme(DimType dim, LevelVector& lmin, LevelVector& lmax) {
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

  ~CombiMinMaxScheme(){};

  /* Generate the combischeme corresponding to the classical combination technique.
   * We need to ensure that lmax = lmin +c*ones(dim), and take special care
   * of dummy dimensions
   * */
  void createClassicalCombischeme();

  /* Generates the adaptive combination scheme (equivalent to CK's
   * Python code)
   * */
  void createAdaptiveCombischeme();

  /*
   * Computes the distribution of a classical scheme to the systems of the
   * third level combination.
   */
  void createClassicalThirdLevelCombischeme(unsigned int systemNumber,
                                            std::vector<bool> boundary,
                                            std::vector<LevelVector>& commonSubspaces);

  /*
   * Computes the distribution of an adaptive scheme to the systems of the
   * third level combination.
   */
  void createAdaptiveThirdLevelCombischeme(unsigned int systemNumber,
                                           std::vector<bool> boundary,
                                           std::vector<LevelVector>& commonSubspaces);

  /* Generates the fault tolerant combination technique with extra
   * grids used in case of faults
   * */
  void makeFaultTolerant();

  inline const std::vector<LevelVector>& getCombiSpaces() const { return combiSpaces_; }

  inline const std::vector<double>& getCoeffs() const { return coefficients_; }

  inline void print(std::ostream& os) const;

 private:
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

  /* Subspaces of the combination technique*/
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

  /*
   * We assume only 2 systems participating.
   * For example purpose we just split the scheme in half, and assign each half
   * to a system.
   * TODO Implement for arbitrary number of systems
   */
  void splitScheme(std::vector<std::vector<LevelVector>>& splittedScheme,
                                      std::vector<std::vector<real>>& splittedCoeffs,
                                      size_t numSystems);
  /*
   * Computes the common subspaces for a given splitted
   */
  std::vector<LevelVector> computeCommonSubspaces(const std::vector<std::vector<LevelVector>>& splittedScheme, const std::vector<bool>& boundary);

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
}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
