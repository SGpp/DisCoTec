#ifndef COMBILINEARBASISFUNCTION_HPP_
#define COMBILINEARBASISFUNCTION_HPP_

#include <sgpp/distributedcombigrid/legacy/CombiBasisFunctionBasis.hpp>
#include <sgpp/distributedcombigrid/legacy/combigrid_utils.hpp>

namespace combigrid {

/** Linear basis function */
class LinearBasisFunction : public combigrid::BasisFunctionBasis {
 public:
  /** empty Ctror */
  LinearBasisFunction() = default;
  virtual ~LinearBasisFunction() = default;

  /** first method which returns the contribution of the first point in the 1D
   * cell
   * @param coord  1D coordonate idealy should be [0,1] but for extrapolation
   * could be different [-1,2]*/
  double functionEval1(double coord) const override { return (1.0 - coord); }

  /** second method which returns the contribution of the second point in the 1D
   * cell
   * @param coord  1D coordonate idealy should be [0,1] but for extrapolation
   * could be different [-1,2]*/
  double functionEval2(double coord) const override { return (coord); }

  // /** return the default basis function*/
  // static const BasisFunctionBasis* getDefaultBasis() { return defaultBasis_; }
 private:
};
}  // namespace combigrid

#endif /* COMBILINEARBASISFUNCTION_HPP_ */
