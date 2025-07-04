#include "discotec/rescheduler/RebalancingTaskRescheduler.hpp"

#include <algorithm>
#include <numeric>
#include <map>

namespace combigrid {

std::vector<std::pair<LevelVector, int>> RebalancingTaskRescheduler::eval(
    const std::map<LevelVector, int>& levelVectorToProcessGroupIndex,
    const std::map<LevelVector, unsigned long>& levelVectorToTaskDuration,
    LoadModel *loadModel) {

  MicrocsecondsLearningLoadModel *sllm = dynamic_cast<MicrocsecondsLearningLoadModel *>(loadModel);

  assert(sllm);
  if (!sllm) {
    // abort scheduling
    return {};
  }
		
	unsigned int iteration = 0;

  // variable containing the result
  std::vector<std::pair<LevelVector, int>> moveTasks{};

  const auto numberOfProcessGroups = levelVectorToProcessGroupIndex.size();

  // we need to mutate this map for multiple iterations
  auto levelVectorToProcessGroupIndexMutCopy = levelVectorToProcessGroupIndex;

  unsigned long averageProcessGroupDuration = 
    std::accumulate(
        begin(levelVectorToTaskDuration), 
        end(levelVectorToTaskDuration), 
        0ul, 
        [](const unsigned long& a, 
           const std::pair<LevelVector, unsigned long>& b) {
          return a + b.second;
        }
      ) / numberOfProcessGroups;

  while (iteration < this->max_iterations_) {

    std::map<int, unsigned long> processGroupIndexToDuration{};
    for (const auto& t : levelVectorToProcessGroupIndexMutCopy) { 
      processGroupIndexToDuration.insert({t.second, 0}); // init process group durations with zero
    }
    for (const auto& t : levelVectorToTaskDuration) { // add up task durations for every process group
      const auto& processGroupIndex = levelVectorToProcessGroupIndexMutCopy[t.first];
      processGroupIndexToDuration[processGroupIndex] += t.second; 
    }

    std::pair<int, unsigned long> fastestGroup;
    std::pair<int, unsigned long> slowestGroup;
    {
      auto minmaxGroup = std::minmax_element(
          std::begin(processGroupIndexToDuration), 
          std::end(processGroupIndexToDuration), 
          [](const std::pair<int, unsigned long>& a, 
             const std::pair<int, unsigned long>& b) {
            return a.second < b.second;
          }
        );
      fastestGroup = *minmaxGroup.first;
      slowestGroup = *minmaxGroup.second;
    }

    if (slowestGroup.second / averageProcessGroupDuration < this->min_inbalance_) {
      return moveTasks;
    }

    std::map<LevelVector, unsigned long> levelVectorToPrognosisOfTaskDurationsOfSlowestProcessGroup;
    for (const auto& t : levelVectorToProcessGroupIndexMutCopy) {
      if (slowestGroup.first == levelVectorToProcessGroupIndexMutCopy[t.first]) {
        levelVectorToPrognosisOfTaskDurationsOfSlowestProcessGroup.insert(
            {t.first, sllm->evalSpecificUOT(t.first).count()}); // use load model for prognosis
      }
    }

    auto moveTask = *std::min_element(
        std::begin(levelVectorToPrognosisOfTaskDurationsOfSlowestProcessGroup), 
        std::end(levelVectorToPrognosisOfTaskDurationsOfSlowestProcessGroup), 
        [&averageProcessGroupDuration, &fastestGroup]
        (const std::pair<LevelVector, unsigned long>& a, 
         const std::pair<LevelVector, unsigned long>& b) { 
          // negative values are possible: need to cast to a signed value for correct results
          return std::abs(static_cast<signed long long>(fastestGroup.second 
                                                        + a.second) 
                          - static_cast<signed long long>(
                                averageProcessGroupDuration)) 
                 < std::abs(static_cast<signed long long>(fastestGroup.second 
                                                          + b.second) 
                            - static_cast<signed long long>(
                                  averageProcessGroupDuration));
        });

    if (fastestGroup.second + moveTask.second < slowestGroup.second) {
      // move task from slowest to fastest group

      bool alreadyMoved = false;
      // check if previously moved: if moved, change group id
      for (size_t i = 0; i < moveTasks.size(); ++i) {
        if (moveTasks[i].first == moveTask.first) {
          alreadyMoved = true;
          moveTasks[i].second = fastestGroup.first;
          break;
        }
      } 
      if (!alreadyMoved) {
        levelVectorToProcessGroupIndexMutCopy[moveTask.first] = fastestGroup.first;
        moveTasks.push_back({moveTask.first, fastestGroup.first});
      }

      ++iteration;
    } else {
      return moveTasks;
    }
  }
  
  return moveTasks;
}

} /* namespace combigrid */
