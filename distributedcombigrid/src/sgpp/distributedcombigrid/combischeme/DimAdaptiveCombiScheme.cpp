#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"

namespace combigrid{

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
	return std::find(std::begin(activeGrids_), std::end(activeGrids_), level) != std::end(activeGrids_);
}

void DimAdaptiveCombiScheme::generateActiveGrids(){
	for(const LevelVector& grid: levels_){
		if(!containsFwdNeighbours(grid)){
			assert(grid >= lmin_);
			assert(containsAllBwdNeighbours(grid));
			activeGrids_.push_back(grid);
		}
	}
}

void DimAdaptiveCombiScheme::addExpansion(const LevelVector& grid){
	addExpansionAllDirections(grid);
}

void DimAdaptiveCombiScheme::addBelowMinSubgrids_h(const LevelVector& grid, const int dim, bool toAdd){
	if(dim == -1){
		if(toAdd){
			assert(!contains(grid));
			levels_.push_back(grid);
		}
	} else {
		if(grid.at(dim) == lmin_.at(dim) && !isDummyDim(dim)){
			LevelVector temp {grid};
			temp.at(dim) = 1;
			for(; temp.at(dim) < lmin_.at(dim); ++temp.at(dim)){
				addBelowMinSubgrids_h(temp, dim - 1, true);
			}
			addBelowMinSubgrids_h(temp, dim - 1, toAdd);
		} else {
			addBelowMinSubgrids_h(grid, dim - 1, toAdd);
		}
	}
}

void DimAdaptiveCombiScheme::addBelowMinSubgrids(const LevelVector& grid){
	assert(grid.size() >= 1);
	addBelowMinSubgrids_h(grid, grid.size() - 1, false);
}

bool DimAdaptiveCombiScheme::hasExpansionNeighbour(const LevelVector& grid) const{
	for(int i = 0; i < dim(); ++i){
		LevelVector fwdNeigh {grid};
		++fwdNeigh.at(i);
		if(isExpansion(fwdNeigh)){
			return true;
		}
	}
	return false;
}

LevelVector DimAdaptiveCombiScheme::errorMeasurePartner(const LevelVector& grid) const{
	assert(grid.size() == dim());
	assert(isActive(grid));
	//find a sensible search direction. This is any direction
	//in which at least one backward neighbour exists (that is greater than lmin)
	int searchDim = std::numeric_limits<int>::max();
	for(DimType i = 0; i < dim(); ++i){
		LevelVector bwdNeigh {grid};
		--bwdNeigh.at(i);
		if(contains(bwdNeigh) && bwdNeigh.at(i) >= lmin_.at(i)){
			searchDim = i;
			break;
		}
	}
	assert(searchDim < dim());

	//find the concrete partner
	size_t partnerIndex = std::numeric_limits<size_t>::max();
	LevelType maxLevel = std::numeric_limits<int>::min();
	for(size_t i = 0; i < combiSpaces_.size(); ++i){
		const LevelVector& combiSpace = combiSpaces_.at(i);
		if(coefficients_.at(i) <= -1 && equalsExceptDim(grid, combiSpace, searchDim) && combiSpace.at(searchDim) > maxLevel){
			maxLevel = combiSpace.at(searchDim);
			partnerIndex = i;
		}
	}
	assert(partnerIndex < combiSpaces_.size());

	return combiSpaces_.at(partnerIndex);
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
		addBelowMinSubgrids(exp);
		activeGrids_.push_back(exp);
		levels_.push_back(exp);
	}

	//remove the grid that is now no longer active
	auto gridPos = std::find(std::begin(activeGrids_), std::end(activeGrids_), grid);
	assert(gridPos != std::end(activeGrids_));
	activeGrids_.erase(gridPos);

	computeCombiCoeffs();
}

} // end namespace combigrid
