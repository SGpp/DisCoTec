#include <algorithm>
#include <cassert>
#include <functional>
#include <random>

#include "sgpp/distributedcombigrid/utils/Config.hpp"

namespace combigrid {
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
      for (auto & c : coord){
        assert( c <= 1. && c >=0. );
      }
  }
  return randomCoords;
}
}
}