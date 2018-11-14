/*
 * LearningLoadModel.hpp
 *
 *      Author: pollinta
 */

#include <iostream>
#include <string>
#include <chrono>
#include <memory>
#include <thread>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {
  namespace durationsFile{

  //the durationInformation data type needs to be communicated to MPI in every process using it
  MPI_Datatype createMPIDurationType(){
    durationInformation dI;
    MPI_Aint baseaddr, procaddr;
    MPI_Address ( &dI,   &baseaddr);
    MPI_Address ( &dI.nProcesses, &procaddr); 
    MPI_Datatype duration_datatype;
    MPI_Datatype types[] = { MPI_LONG, MPI_UNSIGNED }; 
    int blocklengths[] = { 1, 1 };

    MPI_Aint displacements[] = { 0, procaddr - baseaddr };
    MPI_Type_create_struct(
      2, blocklengths, displacements, types,
      &duration_datatype
    );
    MPI_Type_commit(&duration_datatype);

    //verify that the first field has the expected length
    MPI_Aint extent;
    MPI_Type_extent(MPI_LONG, &extent);
    assert (extent == procaddr - baseaddr);
    assert (extent == sizeof(long int) );
    //and the second field, too
    MPI_Type_extent(MPI_UNSIGNED, &extent);
    assert (extent == sizeof(uint) );
    //and the whole structure, too
    MPI_Type_extent(duration_datatype, &extent);
    assert (extent == sizeof(durationInformation) );

    return duration_datatype;
  }

  std::string getFilename(const LevelVector& levelVector){
    return "./loaddata_" + toString(levelVector) + ".mpi_durations";//TODO which directory
  }

  void DurationsWriteFile::write(const Stats::Event event, uint nProcesses){ 
    std::chrono::milliseconds x = std::chrono::duration_cast<std::chrono::milliseconds>(event.end - event.start);
    long int duration = x.count();
    // long int order = event.end.time_since_epoch().count();
    durationInformation info = {duration, nProcesses};
    write(&info);
  }

  std::vector<durationInformation> DurationsReadFile::readFromBeginning(size_t numberOfItems){
    std::vector<durationInformation> content;
    content.resize(numberOfItems);
    durationInformation buf;
    /* Set the file pointer to 0 */
    MPI_File_seek( fh_, 0, MPI_SEEK_SET ); 
    MPI_Status status;
    int flag;
    MPI_File_read(fh_, &content[0], numberOfItems, durationType_, &status);
    MPI_Test_cancelled(&status, &flag);
    // if flag is set, it may be that the file is not ready for reading (yet)
    for (int i=0; i<100; ++i){
      if (flag){
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        flag = 0;
        MPI_File_read(fh_, &content[0], numberOfItems, durationType_, &status);
        MPI_Test_cancelled(&status, &flag);
      }else{
        break;
      }
    }
    assert (!flag);
    // int count = MPI_Get_count(&status, durationType_, &flag); 
    // assert( count == numberOfItems );
    return content;
  }
  
  } /* namespace durationsFile */
 
} /* namespace combigrid */

