/*
 * ProcessGroupCommands.hpp
 *
 *  Created on: Jun 24, 2014
 *      Author: heenemo
 */

#ifndef PROCESSGROUPSIGNALS_HPP_
#define PROCESSGROUPSIGNALS_HPP_

namespace combigrid {

// attention: changing SignalType might require changing the MPI Type as well
typedef int SignalType;

const SignalType RUN_FIRST = 0;
const SignalType RUN_NEXT = 1;
const SignalType EVAL = 2;
const SignalType GRID_EVAL = 3;
const SignalType COMBINE = 4;
const SignalType EXIT = 5;
const SignalType SYNC_TASKS = 6;
const SignalType TEST_CONVERSION = 7;
const SignalType COMBINE_FG = 8;
const SignalType EV_CALC_FG = 9;
const SignalType EV_CALC_FG_INIT = 10;
const SignalType GRID_GATHER = 11;
const SignalType UPDATE_COMBI_PARAMETERS = 12;
const SignalType ADD_TASK = 13;
const SignalType RECOMPUTE = 14;
const SignalType CHECK_DEAD_PROCS = 15; // check for dead workers
const SignalType RECOVER_COMM = 16;
const SignalType PARALLEL_EVAL = 17;
const SignalType COMBINE_THIRD_LEVEL = 18;
const SignalType COMBINE_TO_FILE_THIRD_LEVEL = 19;

typedef int NormalizationType;
const NormalizationType NO_NORMALIZATION = 0;
const NormalizationType L1_NORMALIZATION = 1;
const NormalizationType L2_NORMALIZATION = 2;
const NormalizationType EV_NORMALIZATION = 3;

typedef int FaultSimulationType;
const FaultSimulationType RANDOM_FAIL = 0;
const FaultSimulationType GROUPS_FAIL = 1;

enum TagType {
  signalTag = 0, statusTag = 1, infoTag = 2
};

// attention: changing StatusType might require changing the MPI Type
typedef int StatusType;

const StatusType PROCESS_GROUP_WAIT = 0;
const StatusType PROCESS_GROUP_BUSY = 1;
const StatusType PROCESS_GROUP_FAIL = 2;


} /* namespace combigrid */

#endif /* PROCESSGROUPCOMMANDS_HPP_ */
