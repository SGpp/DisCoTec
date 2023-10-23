#ifndef TASK_HPP_
#define TASK_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <string>
#include <vector>
#include "fault_tolerance/FaultCriterion.hpp"
#include "fault_tolerance/StaticFaults.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"
#include "mpi/MPISystem.hpp"
#include "utils/LevelVector.hpp"
#include "loadmodel/LoadModel.hpp"


namespace combigrid {

/**
 * Interface for tasks
 */

class Task {
 protected:
  Task();

  Task(DimType dim, const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
       LoadModel* loadModel,
       FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})));

  Task(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
       LoadModel* loadModel,
       FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})));

  // cheapest rule of 5 ever
  Task(const Task& other) = delete;
  Task(Task&& other) = delete;
  Task& operator=(const Task& other) = delete;
  Task& operator=(Task&& other) = delete;

  // fault tolerance info
  FaultCriterion* faultCriterion_;

 public:
  virtual ~Task();
  // pgroup receive task
  static void receive(Task** t, RankType source, CommunicatorType comm);

  // manager send task to pgroup root
  static void send(const Task* const t, RankType dest, CommunicatorType comm);

  // broadcast task
  static void broadcast(Task** t, RankType root, CommunicatorType comm);

  inline DimType getDim() const;

  inline const LevelVector& getLevelVector() const;

  inline size_t getID() const;

  inline real getCoefficient() const;

  /**
   * @brief explicitly set a new ID for the task;
   * may be useful if process groups read their task assignment from file instead of receiving it
   * from the manager
   *
   * make sure ID is continuous and unique!
   */
  inline void setID(size_t ID);

  virtual void run(CommunicatorType lcomm) = 0;

  virtual void changeDir(CommunicatorType lcomm) {
    // do nothing
  }

  virtual void init(CommunicatorType lcomm,
                    std::vector<IndexVector> decomposition = std::vector<IndexVector>()) = 0;

  // inline real estimateRuntime() const;

  inline bool isFinished() const;

  inline void setFinished(bool finished);

  // returns the -th fullgrid gathered from all processors
  virtual void getFullGrid(FullGrid<CombiDataType>& fg, RankType lroot, CommunicatorType lcomm,
                           int n = 0) = 0;

  // This method returns the local part of the n-th distributedFullGrid
  virtual DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) = 0;

  virtual const DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n) const;

  virtual void setZero() = 0;

  // override if there is adaptive timestepping in the solver
  virtual real getCurrentTimestep() const { return 0.; }

  // override to really get the current time in the solver
  virtual real getCurrentTime() const {return 0.;}

  virtual CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const {
    std::cout << "No analytical solution for this task! \n";
    return -0.;
  }

  virtual std::vector<IndexVector> getDecomposition() { return std::vector<IndexVector>(); }

  inline virtual bool isInitialized();

  // do task-specific postprocessing (by default: nothing)
  virtual void doDiagnostics() { std::cout << "doDiagnostics called but not implemented"; }

  // do manager-side task-specific postprocessing, if applicable (by default: nothing)
  virtual void receiveDiagnostics() {}

  inline const std::vector<BoundaryType>& getBoundary() const;

 private:
  friend class boost::serialization::access;

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);

 protected:
  DimType dim_;

  LevelVector l_;  // levelvector of partial solution

  real coeff_;  // coefficient of partial solution

  std::vector<BoundaryType> boundary_;

  size_t id_;  // unique id of task, same on manager and worker

  static size_t count;

  LoadModel* loadModel_;

  bool isFinished_;
};

typedef std::vector<Task*> TaskContainer;

inline const LevelVector& getLevelVectorFromTaskID(const TaskContainer& tasks, size_t task_id){
  auto task = std::find_if(tasks.begin(), tasks.end(), 
    [task_id] (Task* t) {return t->getID() == task_id;}
  );
  assert(task != tasks.end());
  return (*task)->getLevelVector();
}

template <class Archive>
void Task::serialize(Archive& ar, const unsigned int version) {
  ar& faultCriterion_;
  ar& dim_;
  ar& id_;
  ar& l_;
  ar& coeff_;
  ar& boundary_;
  ar& isFinished_;
}

inline DimType Task::getDim() const { return dim_; }

inline const LevelVector& Task::getLevelVector() const { return l_; }

inline const std::vector<BoundaryType>& Task::getBoundary() const { return boundary_; }

inline size_t Task::getID() const { return id_; }

inline real Task::getCoefficient() const { return coeff_; }

inline void Task::setID(size_t id) { id_ = id; }

inline bool Task::isFinished() const { return isFinished_; }

inline void Task::setFinished(bool finished) { isFinished_ = finished; }

inline bool Task::isInitialized() {
  std::cout << "Not implemented!!!";
  return false;
}

// inline real Task::estimateRuntime() const { return loadModel_->eval(l_); }
} /* namespace combigrid */

#endif /* TASK_HPP_ */
