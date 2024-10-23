#include "Task.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>

namespace combigrid {

Task::Task() : isFinished_(false) {}

Task::Task(DimType dim, const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
           LoadModel* loadModel, FaultCriterion* faultCrit)
    : faultCriterion_(faultCrit),
      dim_(dim),
      l_(l),
      coeff_(coeff),
      boundary_(boundary),
      id_(count++),
      loadModel_(loadModel),
      isFinished_(false) {
  assert(dim_ > 0);
  assert(l_.size() == dim_);
}

Task::Task(const LevelVector& l, const std::vector<BoundaryType>& boundary, real coeff,
           LoadModel* loadModel, FaultCriterion* faultCrit)
    : Task(static_cast<DimType>(l.size()), l, boundary, coeff, loadModel, faultCrit) {}

Task::~Task() { delete faultCriterion_; }

const DistributedFullGrid<CombiDataType>& Task::getDistributedFullGrid(size_t n) const {
  throw std::runtime_error("const getDistributedFullGrid called but not implemented");
}

size_t Task::count = 0;

void Task::send(const Task* const t, RankType dst, CommunicatorType comm) {
  // save data to archive
  std::stringstream ss;
  {
    boost::archive::text_oarchive oa(ss);
    // write class instance to archive
    oa << t;
  }
  // create mpi buffer of archive
  std::string s = ss.str();
  int bsize = static_cast<int>(s.size());
  char* buf = const_cast<char*>(s.c_str());
  MPI_Send(buf, bsize, MPI_CHAR, dst, TRANSFER_TASK_TAG, comm);
}

void Task::receive(Task** t, RankType src, CommunicatorType comm) {
  // receive size of message
  // todo: not really necessary since size known at compile time
  MPI_Status status;
  int bsize;
  MPI_Probe(src, TRANSFER_TASK_TAG, comm, &status);
  MPI_Get_count(&status, MPI_CHAR, &bsize);

  // create buffer of appropriate size and receive
  std::vector<char> buf(bsize);
  MPI_Recv(&buf[0], bsize, MPI_CHAR, src, TRANSFER_TASK_TAG, comm, MPI_STATUS_IGNORE);

  // create and open an archive for input
  std::string s(&buf[0], bsize);
  std::stringstream ss(s);
  {
    boost::archive::text_iarchive ia(ss);
    // read class state from archive
    ia >> *t;
  }
}

void Task::broadcast(Task** t, RankType root, CommunicatorType comm) {
  RankType myID;
  MPI_Comm_rank(comm, &myID);

  char* buf = NULL;
  int bsize;

  // root writes object data into buffer
  std::string s;

  if (myID == root) {
    // save data to archive
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      // write class instance to archive
      oa << *t;
    }
    // create mpi buffer of archive
    s = ss.str();
    bsize = static_cast<int>(s.size());
    buf = const_cast<char*>(s.c_str());
  }

  // root broadcasts object size
  MPI_Bcast(&bsize, 1, MPI_INT, root, comm);

  // non-root procs create buffer which is large enough
  std::vector<char> tmp(bsize);

  if (myID != root) {
    buf = &tmp[0];
  }

  // broadcast of buffer
  MPI_Bcast(buf, bsize, MPI_CHAR, root, comm);

  // non-root procs write buffer to object
  if (myID != root) {
    // create and open an archive for input
    std::string s(buf, bsize);
    std::stringstream ss(s);
    {
      boost::archive::text_iarchive ia(ss);
      // read class state from archive
      ia >> *t;
    }
  }
}

void Task::getFullGrid(FullGrid<CombiDataType>& fg, RankType lroot, CommunicatorType lcomm, int n) {
  this->getDistributedFullGrid(n).gatherFullGrid(fg, lroot);
}

} /* namespace combigrid */
