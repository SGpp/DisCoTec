#ifndef COMBIBASISFUNCTION_HPP_
#define COMBIBASISFUNCTION_HPP_

#include <boost/serialization/access.hpp>

namespace combigrid {

/** The basis function for 1D Cell, the generalization for ND is simply
 * the tensor product. This class contains two methods since in a 1D cell
 * we have two points at the end of the cell, and the two methods returns
 * the component of the first point and the second returns the component
 * of the second point to the evaluation point <br>
 * All the evaluations should be done on the reference cell [0,1], except
 * the extrapolation, that can be [-1,2]. */

class BasisFunctionBasis {
 public:
  /** empty Ctror */
  BasisFunctionBasis() = default;
  /** virtual Dtor */
  virtual ~BasisFunctionBasis() = default;

 private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {}
};
}  // namespace combigrid

#endif /* COMBIBASISFUNCTION_HPP_ */
