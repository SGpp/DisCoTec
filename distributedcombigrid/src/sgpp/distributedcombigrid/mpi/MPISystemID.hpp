#ifndef MPISYSTEMID_HPP
#define MPISYSTEMID_HPP

#include <tr1/memory>

namespace combigrid {

class MPISystem;

/*!\brief Handle for the MPI communication system.
 // \ingroup mpi
 */
typedef std::shared_ptr<MPISystem> MPISystemID;

/*!\brief Handle for the constant MPI communication system.
 // \ingroup mpi
 */
typedef std::shared_ptr<const MPISystem> ConstMPISystemID;

}  // namespace combigrid

#endif  // MPISYSTEMID_HPP
