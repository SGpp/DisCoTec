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
	      if(containsActive(bwdNeigh)){
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

bool DimAdaptiveCombiScheme::containsOneCriticalBwdNeighbour(const LevelVector& grid) const{
	bool foundOne = false;

	for(std::size_t i = 0; i < dim(); ++i){
	      LevelVector bwdNeigh{grid};
	      --bwdNeigh.at(i);
	      if(isCritical(bwdNeigh)){
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

bool DimAdaptiveCombiScheme::containsOneFwdNeighbour(const LevelVector& grid) const{
  for(std::size_t i = 0; i < dim(); ++i){
      LevelVector fwdNeigh{grid};
      ++fwdNeigh.at(i);
      if(contains(fwdNeigh)){
          return true;
      }
  }

  return false;
}

bool DimAdaptiveCombiScheme::isCritical(const LevelVector& grid) const{
	return !containsOneFwdNeighbour(grid) && containsActive(grid);
}

bool DimAdaptiveCombiScheme::contains(const LevelVector& level) const{
  return std::find(std::begin(levels_), std::end(levels_), level) != std::end(levels_);
}

bool DimAdaptiveCombiScheme::containsActive(const LevelVector& level) const{
  return std::find(std::begin(activeNodes), std::end(activeNodes), level) != std::end(activeNodes);
}

void DimAdaptiveCombiScheme::generateActiveNodes(){
	for(const LevelVector& grid: levels_){
		if(!containsAllFwdNeighbours(grid) && grid >= lmin_){
			if(!containsAllBwdNeighbours(grid)){
				std::cout << grid << '\n';
				assert(false);
			}
			activeNodes.push_back(grid);
		}
	}
	std::cout << "active:\n";
	printActive(std::cout);
}


bool DimAdaptiveCombiScheme::addExpansion(const LevelVector& grid){
  std::cout << "isExpansion: " << isExpansion(grid) << "\n";
  if(!isExpansion(grid)){
	  return false;
  } else {
	  for(int i = 0; i < grid.size(); ++i){
		  if(grid.at(i) == lmin_.at(i)){
			  for(int j = 1; j < lmin_.at(i); ++j){
				  auto temp {grid};
				  temp.at(i) = j;
				  levels_.push_back(temp);
			  }
		  }
	  }
	  //TODO add if on border

	  activeNodes.push_back(grid);
	  levels_.push_back(grid);
	  //check for all bwd neighbours if they are still active or not
	  activeNodes.erase(std::remove_if(std::begin(activeNodes), std::end(activeNodes), [this](const LevelVector& vec){
		  return containsAllFwdNeighbours(vec);
	  }), std::end(activeNodes));

	  computeCombiCoeffs();
	  return true;
  }
}

} // end namespace combigrid
