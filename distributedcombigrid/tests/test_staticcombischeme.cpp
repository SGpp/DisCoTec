#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mpi.h>
#include <iostream>
#include <complex>
#include <cstdarg>
#include <vector>
#include <random>

#include "sgpp/distributedcombigrid/combischeme/DimAdaptiveCombiScheme.hpp"
#include "sgpp/distributedcombigrid/combischeme/StaticCombiScheme.hpp"
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"
#include "test_helper.hpp"

template <class T>
void cmp_vecs(std::vector<T> vec1, std::vector<T> vec2){
	BOOST_REQUIRE(vec1.size() == vec2.size());
	for(size_t i = 0; i < vec1.size(); ++i){
		BOOST_REQUIRE(vec1.at(i) == vec2.at(i));
	}
}

template <class T>
void print(std::vector<T> vec1, std::ostream& stream){
	stream << '[';
	for(size_t i = 0; i < vec1.size(); ++i){
		stream << vec1[i] << ", ";
	}
	stream << ']';
}

void test_helper(combigrid::LevelVector currMin, combigrid::LevelVector currMax, combigrid::LevelType max, combigrid::DimType currDim){
	for(int imin = 1; imin <= max; ++imin){
		for(int imax = imin; imax <= max; ++imax){
			combigrid::LevelVector newMin {currMin};
			combigrid::LevelVector newMax {currMax};
			newMin[currDim] = imin;
			newMax[currDim] = imax;
			if(currDim != 0){
				test_helper(newMin, newMax, max, currDim - 1);
			} else {
				/*
				std::cout << "lmin: ";
				print(newMin, std::cout);
				std::cout << "\n";
				std::cout << "lmax: ";
				print(newMax, std::cout);
				std::cout << "\n";
				 */
				combigrid::StaticCombiScheme staticScheme {combigrid::StaticCombiScheme::createAdaptiveScheme(newMin, newMax)};
				//The comparison with the "classical" generation was not made since it is utterly broken in the CombiMinMaxScheme
				//combigrid::StaticCombiScheme staticScheme {combigrid::StaticCombiScheme::createClassicalScheme(newMin, newMax)};
				//staticScheme.print(std::cout);

				combigrid::CombiMinMaxScheme oldScheme {newMin.size(), newMin, newMax};
				oldScheme.createAdaptiveCombischeme();
				//oldScheme.createClassicalCombischeme();
				//oldScheme.printLevels(std::cout);
				//oldScheme.print(std::cout);

				cmp_vecs(staticScheme.getCombiSpaces(), oldScheme.getCombiSpaces());
				cmp_vecs(staticScheme.getCoeffs(), oldScheme.getCoeffs());
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(same_level_gen) {
	const size_t dim = 4;
	const size_t max = 4;
	test_helper(combigrid::LevelVector(dim, 0), combigrid::LevelVector(dim, 0), max, dim - 1);
}
