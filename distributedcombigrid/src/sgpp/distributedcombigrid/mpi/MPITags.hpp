#pragma once
#include <limits>

namespace combigrid {
// count down from maximum value, as lower tags indicate
// subspaces in the dfg communication

constexpr int TRANSFER_CLASS_TAG = std::numeric_limits<int>::max() / 2 - 1;
constexpr int TRANSFER_DSGU_DATA_TAG = std::numeric_limits<int>::max() / 2 - 2;
constexpr int TRANSFER_SUBSPACE_DATA_SIZES_TAG = std::numeric_limits<int>::max() / 2 - 3;
constexpr int TRANSFER_TASK_TAG = std::numeric_limits<int>::max() / 2 - 4;
constexpr int TRANSFER_LEVAL_TAG = std::numeric_limits<int>::max() / 2 - 5;
constexpr int TRANSFER_SIGNAL_TAG = std::numeric_limits<int>::max() / 2 - 6;
constexpr int TRANSFER_STATUS_TAG = std::numeric_limits<int>::max() / 2 - 7;
constexpr int TRANSFER_TAG = std::numeric_limits<int>::max() / 2 - 8;
constexpr int TRANSFER__TAG = std::numeric_limits<int>::max() / 2 - 9;

// hangup:
// rescheduling/test_2: 3* `task size0` -> jetzt weiter
// thirdLevel: 7*`Worker with rank 0 processed signal 0` -> jetzt weiter <= most promising
// reduce: run first -> jetzt weiter
// ftolerance

// fine: stats, loadmodel, networkutils,hierarchization,fullgrid,distributedfullgrid, task

}  // namespace combigrid