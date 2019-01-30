#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <complex>
#include <cstdarg>
#include <vector>
#include <random>

#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"
#include "test_helper.hpp"

template <class T>
void cmp_vecs1(std::vector<T> vec1, std::vector<T> vec2){
	BOOST_REQUIRE(vec1.size() == vec2.size());
	for(size_t i = 0; i < vec1.size(); ++i){
		BOOST_REQUIRE(vec1.at(i) == vec2.at(i));
	}
}

template <class T>
void print1(std::vector<T> vec1, std::ostream& stream){
	stream << '[';
	for(size_t i = 0; i < vec1.size(); ++i){
		stream << vec1[i] << ", ";
	}
	stream << ']';
}

/*
 * This is a series of expansion which resulted in a bug in
 * one of the previous versions of the adaptive code.
 */
BOOST_AUTO_TEST_CASE(expansion_test1) {
	combigrid::DimAdaptiveCombiScheme scheme {{4,4,4}, {5,5,5}};
	scheme.addExpansion({5,4,4});
	scheme.checkInvariant();
	scheme.addExpansion({6,4,4});
	scheme.checkInvariant();
	scheme.addExpansion({7,4,4});
	scheme.checkInvariant();
	scheme.addExpansion({8,4,4});
	scheme.checkInvariant();
	scheme.addExpansion({4,5,4});
	scheme.checkInvariant();
	scheme.addExpansion({9,4,4});
	scheme.checkInvariant();
	scheme.addExpansion({5,5,4});
	scheme.checkInvariant();
	scheme.addExpansion({4,4,5});
	scheme.checkInvariant();
	scheme.addExpansion({5,4,5});
	scheme.checkInvariant();
	scheme.addExpansion({6,4,5});
	scheme.checkInvariant();
	scheme.addExpansion({7,4,5});
	scheme.checkInvariant();
	scheme.addExpansion({8,4,5});
	scheme.checkInvariant();
	scheme.addExpansion({6,5,4});
	scheme.checkInvariant();
	scheme.addExpansion({7,5,4});
	scheme.checkInvariant();
	scheme.addExpansion({8,5,4});
	scheme.checkInvariant();
	scheme.addExpansion({4,5,5});
	scheme.checkInvariant();
	for(const auto& activeGrid : scheme.getActiveGrids()){
		print1(activeGrid, std::cout);
		std::cout << "\n";
		if(!scheme.hasExpansionNeighbour(activeGrid)){
			continue;
		}
		assert(scheme.isActive(activeGrid));

		//This is the critical part. Here the buggy version didn't find a partner
		//even though the combischeme was valid.
		const combigrid::LevelVector partnerGrid = scheme.errorMeasurePartner(activeGrid);
	}
}
