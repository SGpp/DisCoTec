#include "utils/MonteCarlo.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <random>

#include "utils/Config.hpp"

namespace combigrid {

std::vector<real> serializeInterpolationCoords(
    const std::vector<std::vector<real>>& interpolationCoords) {
  auto coordsSize = interpolationCoords.size() * interpolationCoords[0].size();
  std::vector<real> interpolationCoordsSerial;
  interpolationCoordsSerial.reserve(coordsSize);
  for (const auto& coord : interpolationCoords) {
    interpolationCoordsSerial.insert(interpolationCoordsSerial.end(), coord.begin(), coord.end());
  }
  return interpolationCoordsSerial;
}

std::vector<std::vector<real>> deserializeInterpolationCoords(
    const std::vector<real>& interpolationCoordsSerial, DimType numDimensions) {
  const int dimInt = static_cast<int>(numDimensions);
  assert(interpolationCoordsSerial.size() % dimInt == 0);
  auto numCoords = interpolationCoordsSerial.size() / dimInt;
  std::vector<std::vector<real>> interpolationCoords;
  interpolationCoords.reserve(numCoords);
  auto it = interpolationCoordsSerial.cbegin();
  for (size_t i = 0; i < numCoords; ++i) {
    interpolationCoords.emplace_back(it, it + dimInt);
    std::advance(it, dimInt);
  }
  assert(interpolationCoords[0].size() == numDimensions);
  assert(interpolationCoords[0].size() * interpolationCoords.size() ==
         interpolationCoordsSerial.size());
  return interpolationCoords;
}

namespace montecarlo {
std::vector<std::vector<real>> getRandomCoordinates(int numCoordinates, size_t dim) {
  std::vector<std::vector<real>> randomCoords (numCoordinates, std::vector<real>(dim));
  // cf. https://stackoverflow.com/a/23143753
  // std::random_device rnd_device;
  static std::mt19937 mersenne_engine {8285545262};  // have 1 seed, for reproducible tests
  static std::uniform_real_distribution<> dist {0., 1.};
  static auto gen = [](){
                  return dist(mersenne_engine);
              };
  for (auto & coord : randomCoords) {
      std::generate(begin(coord), end(coord), gen);
      for ([[maybe_unused]] const auto & c : coord){
        assert( c <= 1. && c >=0. );
      }
  }
  return randomCoords;
}

void getNumberSequenceFromSeed(std::vector<real>& randomNumsToBeSet, size_t seed) {
  std::mt19937 mersenne_engine {seed};
  std::uniform_real_distribution<> dist {0., 1.};
  auto gen = [&](){
                  return dist(mersenne_engine);
              };
  std::generate(std::begin(randomNumsToBeSet), std::end(randomNumsToBeSet), gen);
}

real getRandomNumber(real&& a, real&& b) {
  static std::mt19937 mersenne_engine {8285545262};  // have 1 seed, for reproducible tests
  std::uniform_real_distribution<> dist {std::forward<real>(a), std::forward<real>(b)};
  return dist(mersenne_engine);
}
}  // namespace montecarlo
}  // namespace combigrid