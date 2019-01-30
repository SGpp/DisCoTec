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

	DimAdaptiveCombiScheme(LevelVector lmin, LevelVector lmax)
	: StaticCombiScheme{StaticCombiScheme::createClassicalScheme(lmin, lmax)} {
		generateActiveGrids();
		printActive(std::cout);
	}

	/**
	 * returns true if and only if l1 and l2 are equal except for the value at position dim
	 */
	static bool equalsExceptDim(const LevelVector& l1, const LevelVector& l2, DimType dim){
		assert(l1.size() == l2.size());
		for(int i = 0; i < l1.size(); ++i){
			if(l1.at(i) != l2.at(i) && i != dim){
				return false;
			}
		}

		return true;
	}

	bool containsAllBwdNeighboursInv(const LevelVector& grid) const{
		for(std::size_t i = 0; i < dim(); ++i){
			LevelVector bwdNeigh{grid};
			--bwdNeigh.at(i);
			if(bwdNeigh.at(i) >= 1 && !contains(bwdNeigh)){
				return false;
			}
		}

		return true;
	}

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

	/**
	 * Returns true if the scheme contains any forward neighbour of this grid
	 */
	bool containsFwdNeighbours(const LevelVector& grid) const;

	/**
	 * Returns true if the scheme contains the given grid
	 */
	bool contains(const LevelVector& level) const;

	bool isActive(const LevelVector& level) const;

	void expand(std::size_t index);

	/*
	 * Expands the scheme around the given grid.
	 * Currently calls addExpansionAllDirections to do the work.
	 * This structure should allow for a easier switch between different
	 * expansion algorithms.
	 */
	void addExpansion(const LevelVector& grid);

	/*
	 * Adds all subgrids of this grid, which are below the minimum
	 * in at least one dimension and which are not currently containd in the grid.
	 * The reason for these grids to be stored explicitly, is that they are
	 * needed in the sparse grid classes and generating them later on
	 * is much more complicated.
	 */
	void addBelowMinSubgrids_h(const LevelVector& grid, int dim, bool toAdd);
	void addBelowMinSubgrids(const LevelVector& grid);

	/*
	 * Returns true if the given grid has the same index as lmin in the
	 * given dimension
	 */
	bool isBorder(const LevelVector& grid, DimType dimension){
		assert(grid.size() == dim());
		assert(dimension < dim());
		return grid.at(dimension) == lmin_.at(dimension);
	}

	/*
	 * Returns true if the given grid has at least one forward neighbour that
	 * is a possible expansion
	 */
	bool hasExpansionNeighbour(const LevelVector& grid) const;

	LevelVector errorMeasurePartner(const LevelVector& grid) const;

	const std::vector<LevelVector>& getActiveGrids() const noexcept{
		return activeGrids_;
	}

	bool isExpansion(const LevelVector& grid) const{
		bool isPotExp = containsOneActiveBwdNeighbour(grid) && containsAllBwdNeighbours(grid);
		assert(!contains(grid));

		if(!isPotExp){
			return false;
		}

		//One musn't expand in a dummy dimension
		for(int i = 0; i < dim(); ++i){
			if(isDummyDim(i) && grid.at(i) > lmax_.at(i)){
				return false;
			}
		}

		return true;
	}

	void printActive(std::ostream& os) const{
		for (uint i = 0; i < activeGrids_.size(); ++i)
			os << "\t" << i << ". "<< activeGrids_[i] << std::endl;

		os << std::endl;
	}

	void checkInvariant(){
		for(const LevelVector grid : levels_){
			assert(containsAllBwdNeighboursInv(grid));
		}

		for(const LevelVector grid : activeGrids_){
			assert(contains(grid));
			assert(!containsFwdNeighbours(grid));
		}
	}

private:
	std::vector<LevelVector> activeGrids_;

	/*
	 * The implementation of addExpansion which adds all neighbours of the
	 * given grid that are valid expansions to the scheme
	 */
	void addExpansionAllDirections(const LevelVector& grid);


	/**
	 * Generates the active nodes of this scheme from all grids
	 * that are contained in this scheme
	 */
	void generateActiveGrids();
};

}
#endif /* SRC_SGPP_COMBIGRID_COMBISCHEME_COMBIMINMAXSCHEME_HPP_ */
