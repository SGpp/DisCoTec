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
#include <utility>
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/combischeme/StaticCombiScheme.hpp"

namespace combigrid {


class DimAdaptiveCombiScheme : public StaticCombiScheme {
 public:
  //returns true if and only if l1 and l2 are equal except for the value at position dim
  static bool equalsExceptDim(const LevelVector& l1, const LevelVector& l2, DimType dim){
	  assert(l1.size() == l2.size());
	  for(int i = 0; i < l1.size(); ++i){
		  if(l1.at(i) != l2.at(i) && i != dim){
			  return false;
		  }
	  }

	  return true;
  }

  std::vector<LevelVector> generateNonBoundaryBwdNeighbours(const LevelVector& grid) const;
	/**
	 * Returns true if the scheme contains all backward neighbours of this grid
	 */
  bool containsAllBwdNeighbours(const LevelVector& grid) const;

  /**
   * Returns true if the scheme contains exactly one backward neighbour that is active
   */
  bool containsOneActiveBwdNeighbour(const LevelVector& grid) const;

  bool containsOneCriticalBwdNeighbour(const LevelVector& grid) const;

  /**
   * Returns true if the scheme contains all forward neighbours of this grid
   */
  bool containsAllFwdNeighbours(const LevelVector& grid) const;

  bool containsOneFwdNeighbour(const LevelVector& grid) const;

  bool isCritical(const LevelVector& grid) const;

  bool hasActiveFwdNeighbour(const LevelVector& grid) const;

/**
 * Returns true if the scheme contains the given grid
 */
  bool contains(const LevelVector& level) const;

  bool containsActive(const LevelVector& level) const;

  void addNeighbouringExpansions(const LevelVector& grid);

  void generateActiveNodes();

  void expand(std::size_t index);

  bool addExpansion(const LevelVector& grid);

  bool isBorder(const LevelVector& grid, DimType dimension){
	  assert(grid.size() == dim());
	  assert(dimension < dim());
	  return grid.at(dimension) == lmin_.at(dimension);
  }

  std::pair<LevelVector, LevelVector> getPosNegPair(const LevelVector& grid, DimType dimension){
	  assert(dimension < dim());
	  assert(grid.size() == dim());
	  assert(grid.at(dimension) >= lmin_.at(dimension));
	  //If the grid doesn't have a fwd neighbour then it has to be critical so it
	  //has coefficient 1
	  DimType searchDim = dimension;
	  for(DimType i = 0; i < dim(); ++i){
		  LevelVector fwdNeigh {grid};
		  ++fwdNeigh.at(i);
		  if(contains(fwdNeigh)){
			  searchDim = i;
			  break;
		  }
	  }

	  int maxLevel = std::numeric_limits<int>::min();
	  size_t posIndex = std::numeric_limits<size_t>::max();
	  for(size_t i = 0; i < combiSpaces_.size(); ++i){
		  const LevelVector& combiSpace = combiSpaces_.at(i);
		  if(coefficients_.at(i) == 1 && equalsExceptDim(grid, combiSpace, searchDim) && combiSpace.at(searchDim) > maxLevel){
			  maxLevel = combiSpace.at(searchDim);
			  posIndex = i;
		  }
	  }
	  assert(!posIndex <= combiSpaces_.size());
	  LevelVector posGrid = combiSpaces_.at(posIndex);
	  maxLevel = std::numeric_limits<int>::min();
	  size_t negIndex = std::numeric_limits<size_t>::max();
	  for(size_t i = 0; i < combiSpaces_.size(); ++i){
		  const LevelVector& combiSpace = combiSpaces_.at(i);
		  if(coefficients_.at(i) <= -1 && equalsExceptDim(posGrid, combiSpace, dimension) && combiSpace.at(dimension) > maxLevel){
			  maxLevel = combiSpace.at(dimension);
			  negIndex = i;
		  }
	  }
	  assert(!negIndex <= combiSpaces_.size());

	  return std::make_pair(std::move(posGrid), combiSpaces_.at(negIndex));
  }

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
	  //a valid expansions has to fullfill the following constraints:
	  //1. it must not be contained in the scheme already
	  if(contains(grid)){
		  if(grid == LevelVector{11, 3}){
			  std::cout << "contained\n";
		  }
		  return false;
	  }

	  //2. all bwd Neighbours must be active (since they have this grid as possible expansion)
	  int criticalCount = 0;
	  for(const auto& bwdNeigh: generateNonBoundaryBwdNeighbours(grid)){
		  if(!containsActive(bwdNeigh)){
			  if(grid == LevelVector{11, 3}){
				  std::cout << "bwdNeigh: " << bwdNeigh <<  " not contained\n";
			  }
			  return false;
		  }
		  if(isCritical(bwdNeigh)){
			  if(grid == LevelVector{11, 3}){
				  std::cout << "critical: " << bwdNeigh <<  "\n";
			  }
			  ++criticalCount;
		  }
	  }
	  //3. have at most one critical bwd neighbour
	  assert(criticalCount >= 0);
	  return criticalCount <= 1;
  }

  //void expandBest(){
  //  expand(findBestExpansion());
  //}

DimAdaptiveCombiScheme(LevelVector lmin, LevelVector lmax)
    : StaticCombiScheme{StaticCombiScheme::createClassicalScheme(lmin, lmax)} {
  generateActiveNodes();
}

    inline void printActive(std::ostream& os) const{
  	  for (uint i = 0; i < activeNodes.size(); ++i)
  	    os << "\t" << i << ". "<< activeNodes[i] << std::endl;

  	  os << std::endl;
  	}

private:
    std::vector<LevelVector> activeNodes;
};

}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
