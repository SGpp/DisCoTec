#include "Stats.hpp"
#include <chrono>

void Stats::startWallclock()
{
  wallStart_ = std::chrono::high_resolution_clock::now();
}

void Stats::stopWallclock()
{
  auto now = std::chrono::high_resolution_clock::now();
  wallTime_ = (now - wallStart_).count();
}

size_t Stats::getWallTime()
{
  return wallTime_;
}

void Stats::addToBytesTransferredInCombination(size_t dataSize)
{
  totalBytesTransferredInCombination_ += dataSize;
}

void Stats::addToBytesTransferredInSizeExchange(size_t dataSize)
{
  totalBytesTransferredInSizeExchange_ += dataSize;
}

void Stats::increaseNumCombinations()
{
  numCombinations_++;
}

size_t Stats::getTotalBytesTransferredInCombination()
{
  return totalBytesTransferredInCombination_;
}

size_t Stats::getTotalBytesTransferredInSizeExchange()
{
  return totalBytesTransferredInSizeExchange_;
}

size_t Stats::getNumCombinations()
{
  return numCombinations_;
}

