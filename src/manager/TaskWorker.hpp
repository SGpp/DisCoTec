#pragma once

#include "hierarchization/DistributedHierarchization.hpp"
#include "loadmodel/LearningLoadModel.hpp"
#include "task/Task.hpp"
#include "vtk/DFGPlotFileWriter.hpp"

namespace combigrid {

class TaskWorker {
 public:
  explicit TaskWorker() = default;

  TaskWorker(TaskWorker const&) = delete;
  TaskWorker& operator=(TaskWorker const&) = delete;

  ~TaskWorker() { this->deleteTasks(); };

  inline void dehierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                     const std::vector<bool>& hierarchizationDims,
                                     const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                     const LevelVector& lmin);

  inline void deleteTasks();

  inline std::vector<double> getLpNorms(int p) const;

  inline const TaskContainer& getTasks() const;

  inline void hierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                   const std::vector<bool>& hierarchizationDims,
                                   const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                   const LevelVector& lmin);

  inline void initializeTask(Task* t, LevelVectorList const& taskDecomposition,
                             CommunicatorType taskCommunicator);

  inline void removeTask(size_t index);

  inline void runAllTasks();

 private:
  TaskContainer tasks_{};  /// task storage
  int numGridsPerTask_ = 1;
};

inline void TaskWorker::dehierarchizeFullGrids(
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) {
  bool anyNotBoundary =
      std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });
  for (Task* t : this->getTasks()) {
    for (IndexType g = 0; g < this->numGridsPerTask_; g++) {
      auto& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      if (anyNotBoundary) {
        std::remove_reference_t<decltype(lmin)> zeroLMin(lmin.size(), 0);
        DistributedHierarchization::dehierarchizeDFG(dfg, hierarchizationDims, hierarchicalBases,
                                                     zeroLMin);
      } else {
        DistributedHierarchization::dehierarchizeDFG(dfg, hierarchizationDims, hierarchicalBases,
                                                     lmin);
      }
    }
  }
}

inline void TaskWorker::deleteTasks() {
  for (auto& task : this->tasks_) {
    delete task;
    task = nullptr;
  }
  tasks_.clear();
}

inline std::vector<double> TaskWorker::getLpNorms(int p) const {
  // get Lp norm on every worker; reduce through dfg function
  std::vector<double> lpnorms;
  lpnorms.reserve(this->getTasks().size());
  for (const auto& t : this->getTasks()) {
    auto lpnorm = t->getDistributedFullGrid().getLpNorm(p);
    lpnorms.push_back(lpnorm);
  }
  return lpnorms;
}

inline const TaskContainer& TaskWorker::getTasks() const { return tasks_; }

inline void TaskWorker::hierarchizeFullGrids(
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) {
  bool anyNotBoundary =
      std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });
  for (Task* t : this->getTasks()) {
    for (IndexType g = 0; g < this->numGridsPerTask_; g++) {
      auto& dfg = t->getDistributedFullGrid(static_cast<int>(g));
      // hierarchize dfg
      if (anyNotBoundary) {
        std::remove_reference_t<decltype(lmin)> zeroLMin(lmin.size(), 0);
        DistributedHierarchization::hierarchize(dfg, hierarchizationDims, hierarchicalBases,
                                                zeroLMin);
      } else {
        DistributedHierarchization::hierarchize(dfg, hierarchizationDims, hierarchicalBases, lmin);
      }
    }
  }
}

inline void TaskWorker::initializeTask(Task* t, LevelVectorList const& taskDecomposition,
                                       CommunicatorType taskCommunicator) {
  // add task to task storage
  this->tasks_.push_back(t);
  // initalize task
  tasks_.back()->init(taskCommunicator, taskDecomposition);
}

inline void TaskWorker::removeTask(size_t index) {
  delete (tasks_[index]);
  auto removeIt = tasks_.begin();
  std::advance(removeIt, index);
  tasks_.erase(removeIt);
}

inline void TaskWorker::runAllTasks() {
  if (this->getTasks().empty()) {
    std::cout << "Possible error: No tasks! \n";
  }
  for (auto task : this->getTasks()) {
    task->setFinished(false);  // todo: check if this is necessary or move somewhere else
  }
  for (auto task : this->getTasks()) {
    if (!task->isFinished()) {
      task->run(theMPISystem()->getLocalComm());
    }
  }
}
} /* namespace combigrid */
