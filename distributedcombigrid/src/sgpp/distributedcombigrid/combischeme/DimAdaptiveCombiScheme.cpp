#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"

namespace combigrid{

std::vector<LevelVector> DimAdaptiveCombiScheme::generateNonBoundaryBwdNeighbours(const LevelVector& grid) const{
	std::vector<LevelVector> bwdNeighs {};
	bwdNeighs.reserve(grid.size());
	for(int i = 0; i < grid.size(); ++i){
		auto potBwdNeigh = grid;
		--potBwdNeigh.at(i);
		if(potBwdNeigh.at(i) >= lmin_.at(i)){
			bwdNeighs.push_back(potBwdNeigh);
		}
	}
	return bwdNeighs;
}

bool DimAdaptiveCombiScheme::containsAllBwdNeighbours(const LevelVector& grid) const{
  for(std::size_t i = 0; i < dim(); ++i){
      LevelVector bwdNeigh{grid};
      --bwdNeigh.at(i);
      if(bwdNeigh.at(i) >= lmin_.at(i) && !contains(bwdNeigh)){
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
	      if(isActive(bwdNeigh)){
	    	  if(!foundOne){
		          foundOne = true;
	    	  } else {
		    	  //found at least two bwd neighbours that are active
		    	  return false;
	    	  }
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

bool DimAdaptiveCombiScheme::containsFwdNeighbours(const LevelVector& grid) const{
  for(std::size_t i = 0; i < dim(); ++i){
      LevelVector fwdNeigh{grid};
      ++fwdNeigh.at(i);
      if(contains(fwdNeigh)){
          return true;
      }
  }

  return false;
}

bool DimAdaptiveCombiScheme::contains(const LevelVector& level) const{
  return std::find(std::begin(levels_), std::end(levels_), level) != std::end(levels_);
}

bool DimAdaptiveCombiScheme::isActive(const LevelVector& level) const{
  return std::find(std::begin(activeNodes), std::end(activeNodes), level) != std::end(activeNodes);
}

void DimAdaptiveCombiScheme::generateActiveNodes(){
	for(const LevelVector& grid: levels_){
		if(!containsFwdNeighbours(grid)){
			assert(grid >= lmin_);
			assert(containsAllBwdNeighbours(grid));
			activeNodes.push_back(grid);
		}
	}
}

void DimAdaptiveCombiScheme::addExpansion(const LevelVector& grid){
	addExpansionAllDirections(grid);
}

void DimAdaptiveCombiScheme::addBelowMinNeighbours(const LevelVector& grid){
	  for(int i = 0; i < grid.size(); ++i){
		  if(grid.at(i) == lmin_.at(i)){
			  for(int j = 1; j < lmin_.at(i); ++j){
				  LevelVector temp {grid};
				  temp.at(i) = j;
				  levels_.push_back(temp);
			  }
		  }
	  }
}

void DimAdaptiveCombiScheme::addExpansionAllDirections(const LevelVector& grid){
  assert(isActive(grid));
  assert(grid.size() == dim());

  std::vector<LevelVector> expansions;

  for(DimType i = 0; i < dim(); ++i){
	  LevelVector possibleExpansion = grid;
	  ++possibleExpansion.at(i);
	  if(isExpansion(possibleExpansion)){
		  expansions.push_back(possibleExpansion);
	  }
  }

  for(const LevelVector& exp : expansions){
	  addBelowMinNeighbours(exp);
	  activeNodes.push_back(exp);
	  levels_.push_back(exp);
  }

  //remove all bwd neighbours
  activeNodes.erase(std::remove_if(std::begin(activeNodes), std::end(activeNodes), [this](const LevelVector& vec){
	  return containsFwdNeighbours(vec);
  }), std::end(activeNodes));

  computeCombiCoeffs();
}

} // end namespace combigrid
