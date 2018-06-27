/*
 * CombiMinMaxScheme.hpp
 *
 *  Created on: Oct 2, 2015
 *      Author: sccs
 */

#ifndef SRC_SGPP_COMBIGRID_COMBISCHEME_DIMADAPTIVECOMBISCHEME_HPP_
#define SRC_SGPP_COMBIGRID_COMBISCHEME_DIMADAPTIVECOMBISCHEME_HPP_

#include <boost/math/special_functions/binomial.hpp>
#include <numeric>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/StaticCombiScheme.hpp"

namespace combigrid {

class DimAdaptiveCombiScheme : public StaticCombiScheme {
 public:
	/**
	 * Returns true if the scheme contains all backward neighbours of this grid
	 */
  bool containsAllBwdNeighbours(const LevelVector& grid) const;

  /**
   * Returns true if the scheme contains exactly one backward neighbour that is active
   */
  bool containsOneActiveBwdNeighbour(const LevelVector& grid) const;

  /**
   * Returns true if the scheme contains all forward neighbours of this grid
   */
  bool containsAllFwdNeighbours(const LevelVector& grid) const;

  bool hasActiveFwdNeighbour(const LevelVector& grid) const;
/**
 * Returns true if the scheme contains the given grid
 */
  bool contains(const LevelVector& level) const;

  bool containsActive(const LevelVector& level) const;

  void addNeighbouringExpansions(const LevelVector& grid);

  void generateActiveNodes();

  void generatePossibleExpansions();

  void expand(std::size_t index);

  bool addExpansion(const LevelVector& grid);

  const std::vector<LevelVector>& getActiveNodes() const noexcept{
	  return activeNodes;
  }

  const LevelVector getCoeffBwdNeighbour(const LevelVector& activeNode, DimType dimension) const noexcept{
	  assert(dimension < dim() && dimension >= 0);
	  assert(activeNode.size() == dim());
	  if(isDummyDim(dimension)){
		  return LevelVector {};
	  }
	  LevelVector potentialCoeffBwdNeigh = activeNode;
	  while(potentialCoeffBwdNeigh.at(dimension) > 1){
		  --potentialCoeffBwdNeigh.at(dimension);
		  if(contains(potentialCoeffBwdNeigh)){
			  return potentialCoeffBwdNeigh;
		  }

	  }

	  return LevelVector {};
  }

  bool isExpansion(const LevelVector& grid) const{
	  if(grid == LevelVector{11, 3}){
		  std::cout << containsAllBwdNeighbours(grid) << " bla " << containsOneActiveBwdNeighbour(grid) << "\n";
	  }
	  return containsAllBwdNeighbours(grid) && containsOneActiveBwdNeighbour(grid);
  }

  //void expandBest(){
  //  expand(findBestExpansion());
  //}

DimAdaptiveCombiScheme(LevelVector lmin, LevelVector lmax)
    : StaticCombiScheme{StaticCombiScheme::createClassicalScheme(lmin, lmax)}, possibleExpansions{} {
  generateActiveNodes();
}

    inline void printActive(std::ostream& os) const{
  	  for (uint i = 0; i < activeNodes.size(); ++i)
  	    os << "\t" << i << ". "<< activeNodes[i] << std::endl;

  	  os << std::endl;
  	}

private:
    std::vector<LevelVector> activeNodes;
    std::vector<LevelVector> possibleExpansions;
};

}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
