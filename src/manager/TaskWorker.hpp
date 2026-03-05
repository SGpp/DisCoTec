#pragma once

#include "hierarchization/DistributedHierarchization.hpp"
#include "task/Task.hpp"
#include "utils/Types.hpp"

namespace combigrid {
/**
 * @class TaskWorker
 *
 * @brief The TaskWorker class is responsible for managing the component grids, and for operations
 * that use them within the process group only
 */
template <typename CombiDataType = double>
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
   */
  inline void dehierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                     const std::vector<bool>& hierarchizationDims,
                                     const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                     const LevelVector& lmin) {
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

  /**
   * @brief delete all tasks
   */
  inline void deleteTasks() { tasks_.clear(); }

  /**
   * @brief get a reference to the last task
   */
  inline std::unique_ptr<Task<CombiDataType>>& getLastTask() {
    assert(!this->tasks_.empty());
    return this->tasks_.back();
  }

  /**
   * @brief get the Lp norms of the current component grids
   */
  inline std::vector<double> getLpNorms(int p) const {
    std::vector<double> lpnorms;
    lpnorms.reserve(this->getTasks().size());
    for (const auto& t : this->getTasks()) {
      auto lpnorm = t->getDistributedFullGrid().getLpNorm(p);
      lpnorms.push_back(lpnorm);
    }
    return lpnorms;
  }

  /**
   * @brief get a reference to the tasks container
   */
  inline const std::vector<std::unique_ptr<Task<CombiDataType>>>& getTasks() const {
    return tasks_;
  }

  /**
   * @brief hierarchize all full grids
   */
  inline void hierarchizeFullGrids(const std::vector<BoundaryType>& boundary,
                                   const std::vector<bool>& hierarchizationDims,
                                   const std::vector<BasisFunctionBasis*>& hierarchicalBases,
                                   const LevelVector& lmin) {
    bool anyNotBoundary =
        std::any_of(boundary.cbegin(), boundary.cend(), [](BoundaryType b) { return b == 0; });
    for (const auto& t : this->getTasks()) {
      for (IndexType g = 0; g < this->numGridsPerTask_; g++) {
        auto& dfg = t->getDistributedFullGrid(static_cast<int>(g));
        if (anyNotBoundary) {
          std::remove_reference_t<decltype(lmin)> zeroLMin(lmin.size(), 0);
          DistributedHierarchization::hierarchize(dfg, hierarchizationDims, hierarchicalBases,
                                                  zeroLMin);
        } else {
          DistributedHierarchization::hierarchize(dfg, hierarchizationDims, hierarchicalBases,
                                                  lmin);
        }
      }
    }
  }

  /**
   * @brief initialize a task, takes over ownership of the task too
   */
  inline void initializeTask(std::unique_ptr<Task<CombiDataType>> t,
                             LevelVectorList const& taskDecomposition,
                             CommunicatorType taskCommunicator) {
    this->tasks_.emplace_back(std::move(t));
    tasks_.back()->init(taskCommunicator, taskDecomposition);
  }

  /**
   * @brief remove a task from the task storage
   */
  inline void removeTask(size_t index) {
    auto removeIt = tasks_.begin();
    std::advance(removeIt, index);
    tasks_.erase(removeIt);
  }

  /**
   * @brief run all tasks by calling their run method
   */
  inline void runAllTasks() {
    if (this->getTasks().empty()) {
      std::cout << "Possible error: No tasks! \n";
    }
    for (const auto& task : this->getTasks()) {
      task->setFinished(false);
    }
    for (const auto& task : this->getTasks()) {
      if (!task->isFinished()) {
        task->run(theMPISystem()->getLocalComm());
      }
    }
  }

 private:
  std::vector<std::unique_ptr<Task<CombiDataType>>> tasks_{};  /// task storage
  int numGridsPerTask_ = 1;
};
} /* namespace combigrid */
