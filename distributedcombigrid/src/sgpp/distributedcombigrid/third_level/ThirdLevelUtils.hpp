#ifndef THIRDLEVELUTILSHPP_
#define THIRDLEVELUTILSHPP_

#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/third_level/NetworkUtils.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <stdlib.h>
#include <ctime>
#include <sstream>
#include <mutex>
#include <thread>

// does not have to be group 0
// manager process in globalReduceComm should have the same id as group
#define THIRD_LEVEL_GROUP 0
#define THIRD_LEVEL_MANAGER THIRD_LEVEL_GROUP
#define THIRD_LEVEL_PROCESSGROUP_EXCLUSIVE \
      if ( theMPISystem()->getWorldRank() % theMPISystem()->getNumProcs() == THIRD_LEVEL_GROUP )

namespace combigrid {

class ThirdLevelUtils {
 public:
   ThirdLevelUtils(const std::string& remoteHost, int mesgPort);
   ~ThirdLevelUtils();

   void init();
   void finalize();
   int getCommSize();
   int getRank();

 template<typename FG_ELEMENT>
   void allreduce(std::vector<FG_ELEMENT>& sendrecvbuf, unsigned int numparticipants) const;

 template<typename FG_ELEMENT>
   void isend(std::vector<FG_ELEMENT>& buf) const;

  /*
   * sends data to process on remote system with rank source.
   */
  template<typename FG_ELEMENT>
    void send(const std::vector<FG_ELEMENT>& buff, int dest);

  /*
   * receives data from process on remote system with rank source.
   */
  template<typename FG_ELEMENT>
    void recv(std::vector<FG_ELEMENT>& buff, size_t len, int source);

  /*
   * sends data to intermediary who writes the reduced values to a file.
   */
  template<typename FG_ELEMENT> void
    reduceToFileUniform(const std::vector<FG_ELEMENT>& buff);

  /*
   * Blocks until commSize participants have joined to the communicator on 
   * intermediary.
   */
  void barrier(size_t commSize);

 /*
  * Reduces global first and with remote after.
  */
 template<typename FG_ELEMENT> void
   reduceGlobalRemote(DistributedSparseGridUniform<FG_ELEMENT>& dsg);

 private:
   std::string remoteHost_;
   int remoteMesgPort_, rank_;
   ClientSocket* mesgClient_;

   std::mutex mesgLock_;

   /*
    * Talks to intermediary by sending an operation message and receiving the
    * appropriate response.
    */
   void exchangeMesg(std::string& mesg);

   template<typename FG_ELEMENT> void
   extractSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg, std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const;

   template<typename FG_ELEMENT> void
   setSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg, std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const;
};

}

#endif
