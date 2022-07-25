#pragma once
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {
namespace montecarlo {
std::vector<std::vector<real>> getRandomCoordinates(int numCoordinates, size_t dim);
void getNumberSequenceFromSeed(std::vector<real>& randomNumsToBeSet, size_t seed);
real getRandomNumber(real&& a, real&& b);
}
}