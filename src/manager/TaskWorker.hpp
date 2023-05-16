#pragma once

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

  //   void runAllTasks();

  inline void deleteTasks();

  inline const TaskContainer& getTasks() const;

  inline void initializeTask(Task* t, LevelVectorList const& taskDecomposition,
                             CommunicatorType taskCommunicator);

  inline void removeTask(size_t index);

 private:
  TaskContainer tasks_;  /// task storage
};

inline void TaskWorker::deleteTasks() {
  for (auto& task : this->tasks_) {
    delete task;
    task = nullptr;
  }
  tasks_.clear();
}

inline const TaskContainer& TaskWorker::getTasks() const { return tasks_; }

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

} /* namespace combigrid */
