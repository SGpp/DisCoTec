#ifndef COMBILINEARBASISFUNCTION_HPP_
#define COMBILINEARBASISFUNCTION_HPP_

#include <boost/serialization/base_object.hpp>
#include <hierarchization/CombiBasisFunctionBasis.hpp>

namespace combigrid {

/** Linear basis function */
class LinearBasisFunction : public combigrid::BasisFunctionBasis {
 public:
  /** empty Ctror */
  LinearBasisFunction() = default;
  virtual ~LinearBasisFunction() = default;

 private:
   // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};

/** Linear basis functions (= based on B-splines of first order), but with different hierarchical
 * increments */
class HierarchicalHatBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};
class HierarchicalHatPeriodicBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};
class FullWeightingBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};
class FullWeightingPeriodicBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};
class BiorthogonalBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};
class BiorthogonalPeriodicBasisFunction : public combigrid::LinearBasisFunction {
  // serialize
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<BasisFunctionBasis>(*this);
  }
};

}  // namespace combigrid

#endif /* COMBILINEARBASISFUNCTION_HPP_ */
