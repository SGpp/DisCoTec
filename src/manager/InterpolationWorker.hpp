#pragma once
#include "io/H5InputOutput.hpp"
#include "manager/TaskWorker.hpp"
#include "mpi/MPISystem.hpp"
#include "utils/Types.hpp"
#include "vtk/DFGPlotFileWriter.hpp"

namespace combigrid {

template <typename CombinableType>
static std::vector<CombinableType> interpolateValues(
    const std::vector<std::unique_ptr<Task>>& tasks,
    const std::vector<std::vector<real>>& interpolationCoords) {
  auto numCoordinates = interpolationCoords.size();

  // call interpolation function on tasks and reduce with combination coefficient
  std::vector<CombinableType> values(numCoordinates, 0.);
  std::vector<CombinableType> kahanTrailingTerm(numCoordinates, 0.);

  for (const auto& t : tasks) {
    const auto coeff = t->getCoefficient();
#pragma omp parallel for simd
    for (size_t i = 0; i < numCoordinates; ++i) {
      auto localValue = t->getDistributedFullGrid().evalLocal(interpolationCoords[i]);
      auto summand = localValue * coeff;
      // cf. https://en.wikipedia.org/wiki/Kahan_summation_algorithm
      volatile auto y = summand - kahanTrailingTerm[i];
      volatile auto t = values[i] + y;
      kahanTrailingTerm[i] = (t - values[i]) - y;
      values[i] = t;
    }
  }
  // reduce interpolated values within process group
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombinableType>()),
                MPI_SUM, theMPISystem()->getLocalComm());
  // TODO is it necessary to correct for the kahan terms across process groups too?
  //  need to reduce across process groups too
  //  these do not strictly need to be allreduce (could be reduce), but it is easier to maintain
  //  that way (all processes end up with valid values)
  MPI_Allreduce(MPI_IN_PLACE, values.data(), static_cast<int>(numCoordinates),
                abstraction::getMPIDatatype(abstraction::getabstractionDataType<CombinableType>()),
                MPI_SUM, theMPISystem()->getGlobalReduceComm());

  // hope for RVO or change
  return values;
}

static void writeInterpolatedValuesPerGrid(
    const std::vector<std::unique_ptr<Task>>& tasks,
    const std::vector<std::vector<real>>& interpolationCoords, const std::string& fileNamePrefix,
    IndexType currentCombinationStep) {
  // call interpolation function on tasks and write out task-wise
  for (size_t i = 0; i < tasks.size(); ++i) {
    auto taskVals = tasks[i]->getDistributedFullGrid().getInterpolatedValues(interpolationCoords);
    // cycle through ranks to write
    if (i % (theMPISystem()->getNumProcs()) == theMPISystem()->getLocalRank()) {
      std::string saveFilePath =
          fileNamePrefix + "_task_" + std::to_string(tasks[i]->getID()) + ".h5";
      std::string groupName = "run_";
      std::string datasetName = "interpolated_" + std::to_string(currentCombinationStep);
      h5io::writeValuesToH5File(taskVals, saveFilePath, groupName, datasetName,
                                tasks[i]->getCurrentTime());
    }
  }
}

template <typename CombinableType>
static void writeInterpolatedValuesSingleFile(
    const std::vector<std::unique_ptr<Task>>& tasks,
    const std::vector<std::vector<real>>& interpolationCoords, const std::string& filenamePrefix,
    IndexType currentCombinationStep) {
  // all processes interpolate
  auto values = interpolateValues<CombinableType>(tasks, interpolationCoords);
  // one process writes
  OTHER_OUTPUT_GROUP_EXCLUSIVE_SECTION {
    MASTER_EXCLUSIVE_SECTION {
      assert(currentCombinationStep >= 0);
      assert(values.size() > 0);
      std::string groupName = "all_grids";
      std::string datasetName = "interpolated_" + std::to_string(currentCombinationStep);
      std::string valuesWriteFilename =
          filenamePrefix + "_values_" + std::to_string(currentCombinationStep) + ".h5";
      h5io::writeValuesToH5File(values, valuesWriteFilename, groupName, datasetName,
                                tasks.front()->getCurrentTime());
    }
  }
}

static void writeVTKPlotFilesOfAllTasks(const std::vector<std::unique_ptr<Task>>& tasks,
                                        int numberOfGrids) {
  for (const auto& task : tasks) {
#ifdef USE_VTK
    for (int g = 0; g < numberOfGrids; ++g) {
      DistributedFullGrid<CombiDataType>& dfg = task->getDistributedFullGrid(g);
      DFGPlotFileWriter<CombiDataType> writer{dfg, g};
      writer.writePlotFile();
    }
#else
    std::cout << "Warning: no vtk output produced as DisCoTec was compiled without VTK."
              << std::endl;
#endif /* USE_VTK */
  }
}
} /* namespace combigrid */
