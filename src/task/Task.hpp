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
#include "loadmodel/LoadModel.hpp"
#include "mpi/MPISystem.hpp"
#include "utils/LevelVector.hpp"

namespace combigrid {

/**
 * @brief Abstract class for a task that is worked on by a process group (consisting of
 * ProcessGroupWorkers)
 *
 * In the combination technique, this will correspond to the current solution of the PDE for a
 * component grid of a certain level, along with the operators to advance further in time.
 *
 * The task's domain is assumed to be the unit hypercube ((0,1)^d or [0,1)^d [0,1]^d , depending on
 * the boundary conditions)
 */
class Task {
 protected:
  /**
   * @brief Task default constructor, required for serialization
   */
  Task();

  /**
   * @brief Task constructor
   *
   * @param dim dimensionality of the task's domain
   * @param l levelvector of the task
   * @param boundary boundary conditions of the task
   * @param coeff combination coefficient of the task
   * @param loadModel load model for the task, required for dynamic load balancing
   * @param faultCrit fault criterion for the task, required for fault simulation
   */
  Task(DimType dim, const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
       LoadModel* loadModel,
       FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})));

  /**
   * @brief Task constructor
   *
   * @param l levelvector of the task
   * @param boundary boundary conditions of the task
   * @param coeff combination coefficient of the task
   * @param loadModel load model for the task, required for dynamic load balancing
   * @param faultCrit fault criterion for the task, required for fault simulation
   *
   */
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

  /**
   * @brief receive a task from another rank
   */
  static void receive(Task** t, RankType source, CommunicatorType comm);

  /**
   * @brief send a task to another rank
   */
  static void send(const Task* const t, RankType dest, CommunicatorType comm);

  /**
   * @brief broadcast a task to all ranks (typically within the same process group)
   */
  static void broadcast(Task** t, RankType root, CommunicatorType comm);

  /**
   * @brief get the dimensionality of the task's domain
   */
  inline DimType getDim() const;

  /**
   * @brief get the level vector of the task
   */
  inline const LevelVector& getLevelVector() const;

  /**
   * @brief get the unique ID of the task
   */
  inline size_t getID() const;

  /**
   * @brief get the combination coefficient of the task
   */
  inline real getCoefficient() const;

  /**
   * @brief explicitly set a new ID for the task;
   *
   * useful if process groups read their task assignment from file instead of receiving it
   * from the manager
   *
   * make sure ID is continuous and unique!
   */
  inline void setID(size_t ID);

  /**
   * @brief run the task
   *
   * For PDE solvers, this will typically advance the solution to the next timestep
   *
   * @param lcomm local communicator of the process group
   */
  virtual void run(CommunicatorType lcomm) = 0;

  /**
   * @brief change the current working directory
   *
   * useful for tasks that write output files
   *
   * @param lcomm local communicator of the process group
   */
  virtual void changeDir(CommunicatorType lcomm) {
    // do nothing
  }

  /**
   * @brief initialize the task
   *
   * @param lcomm local communicator of the process group
   * @param decomposition optional decomposition of the task's domain
   */
  virtual void init(CommunicatorType lcomm,
                    const std::vector<IndexVector>& decomposition = std::vector<IndexVector>()) = 0;

  /**
   * @brief flag to indicate that the task has finished
   *
   * useful for concurrent setups, needs setFinished() to be called at the end of run()
   */
  inline bool isFinished() const;

  /**
   * @brief set the finished flag
   */
  inline void setFinished(bool finished);

  /**
   * @brief gathers the entire full grid on a single process
   *
   * should only be done for small grids!
   *
   * @param fg the full grid to gather to
   * @param lroot rank of the root process
   * @param lcomm local communicator of the process group
   * @param n index of the distributed full grid to gather
   */
  virtual void getFullGrid(FullGrid<CombiDataType>& fg, RankType lroot, CommunicatorType lcomm,
                           int n = 0);

  /**
   * @brief get the DistributedFullGrid of the task
   *
   * This is where the high-dimensional data is stored. There may be multiple DistributedFullGrid s,
   * to allow for advanced setups with multple fields (e.g. different species in a plasma
   * simulation, concentration and electric field, etc.)
   *
   * @param n index of the distributed full grid to get
   */
  virtual DistributedFullGrid<CombiDataType>& getDistributedFullGrid(size_t n = 0) = 0;

  virtual const DistributedFullGrid<CombiDataType>& getDistributedFullGrid(size_t n) const;

  /**
   * @brief set all DistributedFullGrid values to zero
   */
  virtual void setZero() = 0;

  /**
   * @brief get the current time step of the task
   *
   * useful if there is adaptive timestepping in the solver
   */
  virtual real getCurrentTimestep() const { return 0.; }

  /**
   * @brief get the current time in the task
   *
   * useful if there is adaptive timestepping in the solver
   */
  virtual real getCurrentTime() const { return 0.; }

  /**
   * @brief get the analytical solution of the task
   *
   * only applicable where an analytical solution exists, then it can be used for diagnostics such
   * as error integrals
   *
   * @param coords coordinates in the task's domain
   */
  virtual CombiDataType analyticalSolution(const std::vector<real>& coords, int n = 0) const {
    std::cout << "No analytical solution for this task! \n";
    return -0.;
  }

  /**
   * @brief get the decomposition of task's DistributedFullGrid points to the (Cartesian) process
   * grid
   *
   * @return vector of index vectors, each containing the lower bounds of the points in one
   * dimension
   */
  virtual std::vector<IndexVector> getDecomposition() { return std::vector<IndexVector>(); }

  inline virtual bool isInitialized();

  /**
   * @brief do task-specific postprocessing
   *
   * useful for tasks that need to do additional work after the main task has been completed
   */
  virtual void doDiagnostics() { std::cout << "doDiagnostics called but not implemented"; }

  // do manager-side task-specific postprocessing, if applicable (by default: nothing)

  /**
   * @brief receive diagnostics from the task
   *
   * useful for tasks that need to send additional data to the manager
   */
  virtual void receiveDiagnostics() {}

  /**
   * @brief get the task's boundary flags
   */
  inline const std::vector<BoundaryType>& getBoundary() const;

 private:
  friend class boost::serialization::access;

  /**
   * @brief serialize the task
   *
   * note that the DistributedFullGrids are not serialized here, as this may cause large data
   * volumes. If they need to be serialized, this should be done in the derived classes.
   * If your task needs to serialize additional data, you need to implement this function in your
   * derived class, and put the following into your program:
   * #include "utils/BoostExports.hpp"
   * BOOST_CLASS_EXPORT(MyTask)
   *
   * @param ar archive to serialize to
   * @param version version of the serialization

   */
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

inline const LevelVector& getLevelVectorFromTaskID(const TaskContainer& tasks, size_t task_id) {
  auto task = std::find_if(tasks.begin(), tasks.end(),
                           [task_id](Task* t) { return t->getID() == task_id; });
  assert(task != tasks.end());
  return (*task)->getLevelVector();
}

template <class Archive>
void Task::serialize(Archive& ar, const unsigned int version) {
  ar & faultCriterion_;
  ar & dim_;
  ar & id_;
  ar & l_;
  ar & coeff_;
  ar & boundary_;
  ar & isFinished_;
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
