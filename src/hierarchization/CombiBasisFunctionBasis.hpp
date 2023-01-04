#ifndef COMBIBASISFUNCTION_HPP_
#define COMBIBASISFUNCTION_HPP_

#include <boost/serialization/access.hpp>

namespace combigrid {

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
