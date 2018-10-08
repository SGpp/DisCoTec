/*
 * CombiMinMaxScheme.hpp
 *
 *  Created on: Oct 2, 2015
 *      Author: sccs
 */

#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_STATICCOMBISCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_STATICCOMBISCHEME_HPP_

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include <boost/math/special_functions/binomial.hpp>
#include <numeric>
#include <vector>
#include <algorithm>

namespace combigrid {

/**
 * Used terminology: A full grid represented by its level in each dimension is !contained! in this combiScheme
 * if is directly used in the combination technique or if it contained a subset of one of these grids
 * (meaning there is a grid used in the combischeme for which every level in each dimension is bigger or equal).
 *
 * A dummy Dimension is a dimension where lmin[d] = lmax[d]. These dimensions are treaded specially.
 */
class StaticCombiScheme {
public:

	/**
	 * Creates a classical combi scheme. That means that for every non dummy dimension lmax[j] == lmin[j] + c, where c is a
	 * natural number.
	 * If this invariant doesn't hold for the given parameters a new lmin* vector is created where
	 * lmin*[j] = lmax[j] - k, where k is min_i(lmax[i] - lmin[i]). Both i and j are indices of non dummy dimensions.
	 */
	static StaticCombiScheme createClassicalScheme(LevelVector lmin, LevelVector lmax, bool faultTolerant = false){
		return StaticCombiScheme(makeClassicalMin(lmin, lmax),
				                 lmax, faultTolerant);
  }

	/**
	 * Creates an adaptive combi scheme.
	 */
	static StaticCombiScheme createAdaptiveScheme(LevelVector lmin, LevelVector lmax, bool faultTolerant = false){
		return StaticCombiScheme(lmin, lmax, faultTolerant);
  }

  /* Generate the combischeme corresponding to the classical combination technique.
   * We need to ensure that lmax = lmin +c*ones(dim), and take special care
   * of dummy dimensions
   * */
  void createCombiScheme();

  /**
   * Used in createClassicalScheme to force the lmin to the required format
   * (see createClassicalScheme for further info)
   */
  static LevelVector makeClassicalMin(LevelVector lmin, LevelVector lmax); 

  /**
   * Returns the fullgrids (represented by their level in every dimension)
   * for which a solution should be calculated
   */
   const std::vector<LevelVector>& getCombiSpaces() const noexcept{
    return combiSpaces_;
  }

  /**
   * Returns the coefficients needed to correctly combine the results combiSpaces.
   */
   const std::vector<double>& getCoeffs() const noexcept{
    return coefficients_;
  }

  /**
   * Returns the fullgrids (represented by their level in every dimension)
   * which are either in the comiSpaces or subgrids of these grid.
   *
   * This does not include grids which have levels smaller than lmin
   * since these are trivially in the grid.
   */
  const std::vector<LevelVector>& getLevels() const noexcept{
	 return levels_;
  };

  inline void print(std::ostream& os) const{
	  for (uint i = 0; i < combiSpaces_.size(); ++i)
	    os << "\t" << i << ". "<< combiSpaces_[i] << "\t" << coefficients_[i] << std::endl;

	  os << std::endl;
	}

  void printLevels(std::ostream& os) const {
	  for (uint i = 0; i < levels_.size(); ++i)
	    os << "\t" << i << ". "<< levels_[i] << std::endl;

	  os << std::endl;
  }

  /**
   * Returns the dimension of the grids used
   */
  DimType dim() const noexcept{
    return static_cast<LevelType>(lmin_.size());
  }

  /**
   * Returns the dimension of the grids used
   * minus the dummy dimensions
   */
  std::size_t effDim(){
	  std::size_t count = 0;
	  for(size_t i = 0; i < dim(); ++i){
		  if(!isDummyDim(i)){
			  ++count;
		  }
	  }
	  return count;
  }

  /**
   * Returns true if the given dimension is a dummy dimension
   */
  bool isDummyDim(DimType index) const{
    assert(index < dim());
    return lmax_.at(index) == lmin_.at(index);
  }

protected:
  /* Minimal resolution */
  LevelVector lmin_;

  /* Maximal resolution */
  LevelVector lmax_;

  /* Downset without grids smaller than lmin*/
  std::vector<LevelVector> levels_;

  /* Subspaces of the combination technique*/
  std::vector<LevelVector> combiSpaces_;

  /* Combination coefficients */
  std::vector<real> coefficients_;

  /**
   * Creates the adaptive combination technique from lmax and lmin.
   * If called with an apropriate lmin it creates a classical combi technique
   * (see also createClassicalScheme)
   */
  StaticCombiScheme(const LevelVector& lmin, const LevelVector& lmax, bool faultTolerant) :
  lmin_ {lmin}, lmax_{lmax}, levels_{}, combiSpaces_{}, coefficients_{}{
	assert(lmin_.size() > 0);
    assert(lmax_.size() == lmin_.size());
    std::cout << "test\n";

    for (size_t i = 0; i < dim(); ++i) {
      assert(lmax_[i] > 0);
      assert(lmin_[i] > 0);
      assert(lmax_[i] >= lmin_[i]);
    }

    createLevelsRec();
    printLevels(std::cout);
    computeCombiCoeffs();

    if(faultTolerant){
    	makeFaultTolerant();
    }


  }

  /**
   * Returns the maximum possible sum of the levels of the grids
   * that are contained in this scheme
   */
  LevelType getMaxLevelSum();

  /**
   * The following three functions are used to generate the
   * grids that make up the combi spaces and all of their
   * subgrids.
   *
   * To do that first the grid is transformed from (lmin, lmax)
   * to (1, lmax - lmin) (1 here represents the vector containing only ones).
   * Then the contained grids are generated and in the end lmin-1 is added to these grids again
   */
  LevelType getNameSum();

  void createLevelsRec(){
	  createLevelsRec(getMaxLevelSum(), dim() - 1, LevelVector(dim(), 0));
  }

  void createLevelsRec(LevelType currMaxSum,
          size_t currentIndex,
          LevelVector l);

  /**
   * TODO this method is simply copied from the old version and most likely won't
   * work yet
   */
  void makeFaultTolerant();


  /*
   * Calculate the coefficients for the grids and
   * stores the grid with a non zero coefficient in combiSpaces_
   */
  void computeCombiCoeffs();

  bool isBelowMin(const LevelVector& grid){
	  assert(grid.size() == lmin_.size());
	  for(size_t i = 0; i < grid.size(); ++i){
		  if(grid.at(i) < lmin_.at(i)){
			  return true;
		  }
	  }
	  return false;
  }

};


inline std::ostream& operator<<(std::ostream& os,
    const combigrid::StaticCombiScheme& scheme) {
  scheme.print(os);
  return os;
}

}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
