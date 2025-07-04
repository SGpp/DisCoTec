#pragma once

#include "hierarchization/DistributedHierarchization.hpp"
#include "../Task.hpp"
#include "../utils/Types.hpp"

namespace combigrid {
/**
 * @class TaskWorker
 *
 * @brief The TaskWorker class is responsible for managing the component grids, and for operations
 * that use them within the process group only
 */
class TaskWorker {
 public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  explicit TaskWorker() = default;

  TaskWorker(TaskWorker const&) = delete;
  TaskWorker& operator=(TaskWorker const&) = delete;

  ~TaskWorker() { this->deleteTasks(); };
#endif  // DOXYGEN_SHOULD_SKIP_THIS

  /**
   * @brief dehierarchize all full grids
   *
   * inverse operation to hierarchizeFullGrids
   *
   * @param boundary the boundary conditions used for hierarchization
   * @param hierarchizationDims the hierarchized dimensions
   * @param hierarchicalBases the hierarchical bases used for hierarchization
   * @param lmin the minimum level vector the grids are hierarchized to
   */
  inline void dehierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                     const std::vector<bool>& hierarchizationDims,
                                     const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                     const LevelVector& lmin);

  /**
   * @brief delete all tasks
   */
  inline void deleteTasks();

  /**
   * @brief get a reference to the last task
   */
  inline std::unique_ptr<Task>& getLastTask();

  /**
   * @brief get the Lp norms of the current component grids
   *
   * @param p the p in Lp norm (0 for maximum norm)
   * @return a vector of results, one for each task
   */
  inline std::vector<double> getLpNorms(int p) const;

  /**
   * @brief get a reference to the tasks container
   */
  inline const std::vector<std::unique_ptr<Task>>& getTasks() const;

  /**
   * @brief hierarchize all full grids
   *
   * inverse operation to dehierarchizeFullGrids
   *
   * @param boundary the boundary conditions used for hierarchization
   * @param hierarchizationDims the dimensions to hierarchize
   * @param hierarchicalBases the hierarchical bases used for hierarchization
   * @param lmin the minimum level vector to hierarchize to
   */
  inline void hierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                   const std::vector<bool>& hierarchizationDims,
                                   const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                   const LevelVector& lmin);

  /**
   * @brief initialize a task, takes over ownership of the task too
   *
   * @param t the task to initialize
   * @param taskDecomposition the decomposition of the task's points to the process grid
   * @param taskCommunicator the communicator for the task
   */
  inline void initializeTask(std::unique_ptr<Task> t, LevelVectorList const& taskDecomposition,
                             CommunicatorType taskCommunicator);

  /**
   * @brief remove a task from the task storage
   */
  inline void removeTask(size_t index);

  /**
   * @brief run all tasks by calling their run method
   */
  inline void runAllTasks();

 private:
  std::vector<std::unique_ptr<Task>> tasks_{};  /// task storage
  int numGridsPerTask_ = 1;
};

inline void TaskWorker::dehierarchizeFullGrids(
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) {
  bool anyNotBoundary =
      std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });
  for (auto& t : this->getTasks()) {
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

inline void TaskWorker::deleteTasks() { tasks_.clear(); }

inline std::unique_ptr<Task>& TaskWorker::getLastTask() {
  assert(!this->tasks_.empty());
  return this->tasks_.back();
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

inline const std::vector<std::unique_ptr<Task>>& TaskWorker::getTasks() const { return tasks_; }

inline void TaskWorker::hierarchizeFullGrids(
    const std::vector<BoundaryType>& boundary, const std::vector<bool>& hierarchizationDims,
    const std::vector<BasisFunctionBasis*>& hierarchicalBases, const LevelVector& lmin) {
  bool anyNotBoundary =
      std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });
  for (const auto& t : this->getTasks()) {
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

inline void TaskWorker::initializeTask(std::unique_ptr<Task> t,
                                       LevelVectorList const& taskDecomposition,
                                       CommunicatorType taskCommunicator) {
  // add task to task storage
  this->tasks_.emplace_back(std::move(t));
  // initalize task
  tasks_.back()->init(taskCommunicator, taskDecomposition);
}

inline void TaskWorker::removeTask(size_t index) {
  auto removeIt = tasks_.begin();
  std::advance(removeIt, index);
  tasks_.erase(removeIt);
}

inline void TaskWorker::runAllTasks() {
  if (this->getTasks().empty()) {
    std::cout << "Possible error: No tasks! \n";
  }
  for (const auto& task : this->getTasks()) {
    task->setFinished(false);  // todo: check if this is necessary or move somewhere else
  }
  for (const auto& task : this->getTasks()) {
    if (!task->isFinished()) {
      task->run(theMPISystem()->getLocalComm());
    }
  }
}
} /* namespace combigrid */
