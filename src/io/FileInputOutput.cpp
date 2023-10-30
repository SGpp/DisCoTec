#include "io/FileInputOutput.hpp"

#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

namespace combigrid {
void writeConcatenatedFileRootOnly(const char* data, size_t sizeOfData, const std::string& path,
                                   MPI_Comm comm, bool replaceExistingFile) {
  int mpi_size, mpi_rank;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);

  std::vector<int> recvCounts;

  if (mpi_rank == 0) {
    recvCounts.resize(mpi_size);
  }
  int sizeAsInt = static_cast<int>(sizeOfData);
  MPI_Datatype datatype = MPI_INT;
  MPI_Gather(&sizeAsInt, 1, datatype, recvCounts.data(), 1, datatype, 0, comm);

  std::string recvBuffer;
  std::vector<int> displacement;

  if (mpi_rank == 0) {
    displacement.resize(mpi_size);
    std::exclusive_scan(recvCounts.begin(), recvCounts.end(), displacement.begin(), 0);
    recvBuffer.resize(std::accumulate(recvCounts.begin(), recvCounts.end(), 0));
  }

  MPI_Gatherv(data, sizeOfData, MPI_CHAR, &recvBuffer[0], recvCounts.data(), displacement.data(),
              MPI_CHAR, 0, comm);

  if (mpi_rank == 0) {
    if (replaceExistingFile == true) {
      bool fileExist = std::filesystem::exists(path);
      if (fileExist == true) {
        std::remove(path.c_str());
      }
    }
    std::ofstream out(path.c_str());
    out << recvBuffer;
    out.close();
  }
}

bool getFileExistsRootOnly(const std::string& fileName, CommunicatorType comm, RankType myRank,
                           RankType rootRank) {
  bool doesItExist = false;
  if (myRank == rootRank) {
    doesItExist = std::filesystem::exists(fileName);
  }
  MPI_Bcast(&doesItExist, 1, MPI_C_BOOL, rootRank, comm);
  return doesItExist;
}
}  // namespace combigrid