#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"

namespace combigrid{
bool DimAdaptiveCombiScheme::containsAllBwdNeighbours(const LevelVector& grid) const{
  for(std::size_t i = 0; i < dim(); ++i){
      LevelVector bwdNeigh{grid};
      --bwdNeigh.at(i);
      if(!contains(bwdNeigh)){
          return false;
      }
  }

  return true;
}
bool DimAdaptiveCombiScheme::containsOneActiveBwdNeighbour(const LevelVector& grid) const{
	bool foundOne = false;

	for(std::size_t i = 0; i < dim(); ++i){
	      LevelVector bwdNeigh{grid};
	      --bwdNeigh.at(i);
	      if(containsActive(bwdNeigh) && !foundOne){
	          foundOne = true;
	      } else {
	    	  //found at least two bwd neighbours that are active
	    	  return false;
	      }
	  }

	  return foundOne;
}

bool DimAdaptiveCombiScheme::containsAllFwdNeighbours(const LevelVector& grid) const{
  for(std::size_t i = 0; i < dim(); ++i){
      LevelVector fwdNeigh{grid};
      ++fwdNeigh.at(i);
      if(!contains(fwdNeigh)){
          return false;
      }
  }

  return true;
}

bool DimAdaptiveCombiScheme::contains(const LevelVector& level) const{
  return std::find(std::begin(levels_), std::end(levels_), level) != std::end(levels_);
}

bool DimAdaptiveCombiScheme::containsActive(const LevelVector& level) const{
  return std::find(std::begin(activeNodes), std::end(activeNodes), level) != std::end(activeNodes);
}

void DimAdaptiveCombiScheme::addNeighbouringExpansions(const LevelVector& grid){
    for(std::size_t i = 0; i < dim(); ++i){
      LevelVector fwdNeigh{grid};
      ++fwdNeigh.at(i);
      if (!isExpansion(fwdNeigh)) {
    	  assert(!contains(fwdNeigh));
          possibleExpansions.insert(std::upper_bound(std::begin(possibleExpansions), std::end(possibleExpansions), fwdNeigh ),
          fwdNeigh);
      }
    }
}

void DimAdaptiveCombiScheme::generateActiveNodes(){
	for(const LevelVector& grid: levels_){
		if(!containsAllFwdNeighbours(grid)){
			assert(containsAllBwdNeighbours(grid));
			activeNodes.push_back(grid);
		}
	}
}

void DimAdaptiveCombiScheme::generatePossibleExpansions() {
  // TODO implement more intelligent approach that doesn't generate possible expansions more than once
  for (const LevelVector& grid : activeNodes) {
      	addNeighbouringExpansions(grid);
  }

  possibleExpansions.erase(std::unique(possibleExpansions.begin(), possibleExpansions.end()), possibleExpansions.end());
}

void DimAdaptiveCombiScheme::expand(std::size_t index){
  assert(index < possibleExpansions.size());
  const LevelVector& expansionGrid = possibleExpansions.at(index);
  possibleExpansions.erase(std::begin(possibleExpansions) + index);
  addNeighbouringExpansions(expansionGrid);
  levels_.push_back(expansionGrid);
  computeCombiCoeffs();
}

bool DimAdaptiveCombiScheme::addExpansion(const LevelVector& grid){
  auto gridPos = std::find(std::begin(possibleExpansions), std::end(possibleExpansions), grid);
  if(gridPos == std::end(possibleExpansions)){
	  return false;
  } else {
	  expand(gridPos - std::begin(possibleExpansions));
	  return true;
  }
}

} // end namespace combigrid
