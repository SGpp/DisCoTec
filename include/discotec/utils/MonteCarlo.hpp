#pragma once
#include "Types.hpp"

namespace combigrid {

std::vector<real> serializeInterpolationCoords(
    const std::vector<std::vector<real>>& interpolationCoords);

std::vector<std::vector<real>> deserializeInterpolationCoords(
    const std::vector<real>& interpolationCoordsSerial, DimType numDimensions);

namespace montecarlo {
std::vector<std::vector<real>> getRandomCoordinates(int numCoordinates, size_t dim);

void getNumberSequenceFromSeed(std::vector<real>& randomNumsToBeSet, size_t seed);

real getRandomNumber(real&& a, real&& b);
}  // namespace montecarlo
}  // namespace combigrid