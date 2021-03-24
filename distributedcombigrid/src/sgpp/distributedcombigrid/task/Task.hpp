#ifndef TASK_HPP_
#define TASK_HPP_

#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <string>
#include <vector>
#include "sgpp/distributedcombigrid/fault_tolerance/FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridEnsemble.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include "sgpp/distributedcombigrid/loadmodel/LoadModel.hpp"


namespace combigrid {

/**
 * Interface for tasks
 */

class Task {
 protected:
  Task();

  Task(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
       LoadModel* loadModel, FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})));

  // fault tolerance info
  FaultCriterion* faultCriterion_;

 public:
  virtual ~Task();
  // pgroup receive task
  static void receive(Task** t, RankType source, CommunicatorType comm);

  // manager send task to pgroup root
  static void send(Task** t, RankType dest, CommunicatorType comm);

  // broadcast task
  static void broadcast(Task** t, RankType root, CommunicatorType comm);

  inline DimType getDim() const;

  inline const LevelVector& getLevelVector() const;

  inline const std::vector<bool>& getBoundary() const;

  inline int getID() const;
  virtual size_t getDOFs(){return 1;};

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

  // return the number of grids in a given task; default is one
  virtual size_t getNumGrids() {return 1;}

  virtual void setZero() = 0;

  // override if there is adaptive timestepping in the solver
  virtual real getCurrentTimestep() const { return 0.; }

  // override to really get the current time in the solver
  virtual real getCurrentTime() const {return 0.;}

  virtual void decideToKill() { std::cout << "Kill function not implemented for this task! \n"; }

  virtual CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const {
    std::cout << "No analytical solution for this task! \n";
    return -0.;
  }

  virtual std::vector<IndexVector> getDecomposition() { return std::vector<IndexVector>(); }

  inline virtual bool isInitialized();

  real initFaults(real t_fault, std::chrono::high_resolution_clock::time_point startTimeIteration) {
    return faultCriterion_->init(startTimeIteration, t_fault);
  }

  // override if there is a DFG ensemble in the task
  virtual DFGEnsemble* getDFGEnsemble() const { return nullptr; }

 private:
  friend class boost::serialization::access;

  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version);

 protected:
  DimType dim_;

  LevelVector l_;  // levelvector of partial solution

  std::vector<bool> boundary_;

  int id_;  // unique id of task, same on manager and worker

  static int count;

  LoadModel* loadModel_;

  bool isFinished_;
};

typedef std::vector<Task*> TaskContainer;

inline const LevelVector& getLevelVectorFromTaskID(TaskContainer tasks, int task_id){
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
  ar& boundary_;
  ar& isFinished_;
}

inline DimType Task::getDim() const { return dim_; }

inline const LevelVector& Task::getLevelVector() const { return l_; }

inline const std::vector<bool>& Task::getBoundary() const { return boundary_; }

inline int Task::getID() const { return id_; }

inline bool Task::isFinished() const { return isFinished_; }

inline void Task::setFinished(bool finished) { isFinished_ = finished; }

inline bool Task::isInitialized() {
  std::cout << "Not implemented!!!";
  return false;
}

// inline real Task::estimateRuntime() const { return loadModel_->eval(l_); }
} /* namespace combigrid */

#endif /* TASK_HPP_ */
