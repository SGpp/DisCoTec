#pragma once

namespace combigrid {
// count down from maximum value, as lower tags indicate
// subspaces in the dfg communication

//constexpr int MAX_TAG =  MPI_TAG_UB; // avoided, because throws errors with OpenMPI
constexpr int MAX_TAG = 32767; //from MPI Standard
constexpr int TRANSFER_CLASS_TAG = MAX_TAG - 1;
constexpr int TRANSFER_DSGU_DATA_TAG = MAX_TAG - 2;
constexpr int TRANSFER_SUBSPACE_DATA_SIZES_TAG = MAX_TAG - 3;
constexpr int TRANSFER_TASK_TAG = MAX_TAG - 4;
constexpr int TRANSFER_LEVAL_TAG = MAX_TAG - 5;
constexpr int TRANSFER_SIGNAL_TAG = MAX_TAG - 6;
constexpr int TRANSFER_STATUS_TAG = MAX_TAG - 7;
constexpr int TRANSFER_TAG = MAX_TAG - 8;
constexpr int TRANSFER_GHOST_LAYER_TAG = MAX_TAG - 9;
constexpr int TRANSFER__TAG = MAX_TAG - 10;

}  // namespace combigrid