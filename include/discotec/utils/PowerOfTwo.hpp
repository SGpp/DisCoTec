#pragma once

#include <vector>

#include "LevelVector.hpp"
#include "Types.hpp"

namespace combigrid {

// compute powers of two quickly by bit-shifting
inline IndexType powerOfTwoByBitshift(LevelType x) {
  assert(static_cast<long unsigned int>(x) < sizeof(IndexType) * 8);
  assert(x >= 0 && "powerOfTwoByBitshift: negative argument");
  return IndexType{ 1 } << x;
}

/** vector with power two */
const int powerOfTwo[30] = {
    1,       2,       4,       8,       16,       32,       64,       128,       256,      512,
    1024,    2048,    4096,    8192,    16384,    32768,    65536,    131072,    262144,   524288,
    1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456};

/** vector with one over power two */
const double oneOverPowOfTwo[30] = {
    1.0 / 1.0,        1.0 / 2.0,        1.0 / 4.0,         1.0 / 8.0,        1.0 / 16.0,
    1.0 / 32.0,       1.0 / 64.0,       1.0 / 128.0,       1.0 / 256.0,      1.0 / 512.0,
    1.0 / 1024.0,     1.0 / 2048.0,     1.0 / 4096.0,      1.0 / 8192.0,     1.0 / 16384.0,
    1.0 / 32768.0,    1.0 / 65536.0,    1.0 / 131072.0,    1.0 / 262144.0,   1.0 / 524288.0,
    1.0 / 1048576.0,  1.0 / 2097152.0,  1.0 / 4194304.0,   1.0 / 8388608.0,  1.0 / 16777216.0,
    1.0 / 33554432.0, 1.0 / 67108864.0, 1.0 / 134217728.0, 1.0 / 268435456.0};

}  // namespace combigrid
