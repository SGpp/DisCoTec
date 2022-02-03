#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <numeric>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

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

  virtual ~CombiMinMaxScheme() = default;

  /* Generate the combischeme corresponding to the classical combination technique.
   * We need to ensure that lmax = lmin +c*ones(dim), and take special care
   * of dummy dimensions
   * */
  void createClassicalCombischeme();

  /* Generates the adaptive combination scheme (equivalent to CK's
   * Python code)
   * */
  void createAdaptiveCombischeme();

  /* Generates the fault tolerant combination technique with extra
   * grids used in case of faults
   * */
  void makeFaultTolerant();

  inline const std::vector<LevelVector>& getCombiSpaces() const { return combiSpaces_; }

  inline const std::vector<double>& getCoeffs() const { return coefficients_; }

  inline void print(std::ostream& os) const;

 protected:
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

class CombiMinMaxSchemeFromFile : public CombiMinMaxScheme {
 public:
  CombiMinMaxSchemeFromFile(DimType dim, LevelVector& lmin, LevelVector& lmax, std::string ctschemeFile)
      : CombiMinMaxScheme(dim, lmin, lmax) {
    // read in CT scheme, if applicable
    boost::property_tree::ptree pScheme;
    boost::property_tree::json_parser::read_json(ctschemeFile, pScheme);
    for (const auto& component : pScheme.get_child("")) {
      assert(component.first.empty());  // list elements have no names
      for (const auto& c : component.second) {
        if (c.first == "coeff") {
          coefficients_.push_back(c.second.get_value<real>());
        } else if (c.first == "level") {
          LevelVector lvl(dim);
          int i = 0;
          for (const auto& l : c.second) {
            lvl[i] = l.second.get_value<int>();
            ++i;
          }
          assert(lvl <= lmax);
          assert(lmin <= lvl);
          combiSpaces_.push_back(lvl);
        } else if (c.first == "group_no") {
          processGroupNumbers_.push_back(c.second.get_value<size_t>());
        } else {
          assert(false);
        }
      }
      // process group numbers should either be always given or never
      assert(coefficients_.size() == combiSpaces_.size());
      assert(processGroupNumbers_.size() == combiSpaces_.size() || processGroupNumbers_.size() == 0);
    }
    assert(coefficients_.size() > 0);
    assert(coefficients_.size() == combiSpaces_.size());
  }

  inline const std::vector<size_t>& getProcessGroupNumbers() const { return processGroupNumbers_; }
 private:
  std::vector<size_t> processGroupNumbers_;
};
}  // namespace combigrid
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
