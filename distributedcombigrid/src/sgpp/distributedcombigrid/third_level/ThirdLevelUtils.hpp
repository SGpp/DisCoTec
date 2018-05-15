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
#define RESPONSE_THREAD_TIMEOUT 120
#define AWAIT_RESPONSE_TIMEOUT 120

namespace combigrid {

class ThirdLevelUtils {
 public:
   ThirdLevelUtils(const std::string& remoteHost, int mesgPort, int dataPort);
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
  * Reduces global first and with remote after.
  */
 template<typename FG_ELEMENT> void
   reduceGlobalRemote(DistributedSparseGridUniform<FG_ELEMENT>& dsg);

 private:
   std::string remoteHost_;
   int remoteMesgPort_, remoteDataPort_, rank_;
   ClientSocket *mesgClient_, *dataClient_;

   std::mutex mesgLock_;

   long mesgCount_;

   // TODO better variable names
   std::map<long, std::string> responses_;
   std::mutex responsesLock_;
   std::thread responseThread_;
   int remaining_; // remaining time for responseThread_;
   std::mutex timeoutLock_; // lock to change the remaining timeout

   void fetchResponse();
   bool awaitResponse(long mesgID);
   long sendMessage(const std::string& mesg);

   template<typename FG_ELEMENT> void
   extractSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg, std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const;

   template<typename FG_ELEMENT> void
   setSubspaces(DistributedSparseGridUniform<FG_ELEMENT>& dsg, std::vector<FG_ELEMENT>& buf, std::vector<int>& subspaceSizes) const;

   std::string getResponse() const;
};

}

#endif
