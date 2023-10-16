#ifndef SRC_SGPP_COMBIGRID_MPI_MPIUTILS_HPP_
#define SRC_SGPP_COMBIGRID_MPI_MPIUTILS_HPP_

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>
#include <string>

#include "mpi/MPITags.hpp"

namespace combigrid {

class MPIUtils {
 public:
  template <typename T>
  static void sendClass(T* t, RankType dst, CommunicatorType comm, int tag = TRANSFER_CLASS_TAG) {
    // save data to archive
    std::stringstream ss;
    {
      boost::archive::text_oarchive oa(ss);
      // write class instance to archive
      oa << *t;
    }
    // create mpi buffer of archive
    std::string s = ss.str();
    int bsize = static_cast<int>(s.size());
    char* buf = const_cast<char*>(s.c_str());
    MPI_Send(buf, bsize, MPI_CHAR, dst, tag, comm);
  }

  template <typename T>
  static void receiveClass(T* t, RankType src, CommunicatorType comm, int tag = TRANSFER_CLASS_TAG) {
    // receive size of message
    // todo: not really necessary since size known at compile time
    MPI_Status status;
    int bsize;
    MPI_Probe(src, tag, comm, &status);
    MPI_Get_count(&status, MPI_CHAR, &bsize);

    // create buffer of appropriate size and receive
    std::vector<char> buf(bsize);
    MPI_Recv(&buf[0], bsize, MPI_CHAR, src, tag, comm, MPI_STATUS_IGNORE);

    // create and open an archive for input
    std::string s(&buf[0], bsize);
    std::stringstream ss(s);
    {
      boost::archive::text_iarchive ia(ss);
      // read class state from archive
      ia >> *t;
    }
  }

  template <typename T>
  static void broadcastContainer(T& v, RankType root, CommunicatorType comm) {
    // root broadcasts object size
    int bsize = static_cast<int>(v.size());
    MPI_Bcast(&bsize, 1, MPI_INT, root, comm);
    // non-root procs resize buffer to large enough
    v.resize(bsize);

    // broadcast of buffer
    MPI_Bcast(
        v.data(), bsize,
        abstraction::getMPIDatatype(abstraction::getabstractionDataType<typename T::value_type>()),
        root, comm);
  }

  template <typename T>
  static void broadcastClass(T* t, RankType root, CommunicatorType comm) {
    RankType myID;
    MPI_Comm_rank(comm, &myID);

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
    }

    // root broadcasts object size
    broadcastContainer(s, root, comm);

    // non-root procs write buffer to object
    if (myID != root) {
      // create and open an archive for input
      std::stringstream ss(s);
      {
        boost::archive::text_iarchive ia(ss);
        // read class state from archive
        ia >> *t;
      }
    }
  }
};
}  // namespace combigrid

#endif /* SRC_SGPP_COMBIGRID_MPI_MPIUTILS_HPP_ */
