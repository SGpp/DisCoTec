#include <chrono>

class Stats
{
  private:
    std::chrono::_V2::high_resolution_clock::time_point wallStart_;
    size_t wallTime_ = 0;
    size_t totalBytesTransferredInCombination_ = 0;
    size_t totalBytesTransferredInSizeExchange_ = 0;
    size_t numCombinations_ = 0;

  public:
    void startWallclock();

    void stopWallclock();

    void addToBytesTransferredInCombination(size_t dataSize);

    void addToBytesTransferredInSizeExchange(size_t dataSize);

    void increaseNumCombinations();

    size_t getWallTime();

    size_t getTotalBytesTransferredInCombination();

    size_t getTotalBytesTransferredInSizeExchange();

    size_t getNumCombinations();
};
