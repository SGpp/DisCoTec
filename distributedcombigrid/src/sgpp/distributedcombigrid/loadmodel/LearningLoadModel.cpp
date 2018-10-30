/*
 * LearningLoadModel.hpp
 *
 *      Author: pollinta
 */

#include <iostream>
#include <string>
#include <chrono>
#include <memory>

#include "sgpp/distributedcombigrid/loadmodel/LearningLoadModel.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {
  namespace durationsFile{

  //the durationInformation data type needs to be communicated to MPI in every process using it
  MPI_Datatype createMPIDurationType(){
    MPI_Datatype duration_datatype;
    MPI_Datatype types[] = { MPI_LONG, MPI_UNSIGNED }; 
    int blocklengths[] = { 1, 1 };
    MPI_Aint extent;
    MPI_Type_extent(MPI_LONG, &extent);
    MPI_Aint displacements[] = {0, extent};
    MPI_Type_create_struct(
      2, blocklengths, displacements, types,
      &duration_datatype
    );
    MPI_Type_commit(&duration_datatype);
    return duration_datatype;
  }

  std::string getFilename(const LevelVector& levelVector){
    return "./loaddata_" + toString(levelVector) + ".mpi_durations";//TODO which directory
  }

  void DurationsWriteFile::write(const Stats::Event event, size_t nProcesses){ //TODO include "metadata" in model: nrg.dat, parameters etc.
    long int duration = (event.end - event.start).count(); 
    // long int order = event.end.time_since_epoch().count();
    durationInformation info = {duration, nProcesses};
    write(&info);
  }

  std::vector<durationInformation> DurationsReadFile::read_all(){
    std::vector<durationInformation> content;
    durationInformation buf;
    // while (){ //TODO
    //   MPI_File_read_all(fh_,buf,1,durationType_,MPI_STATUS_IGNORE);
    // }
    return content;
  }
  
  } /* namespace durationsFile */
 
} /* namespace combigrid */

