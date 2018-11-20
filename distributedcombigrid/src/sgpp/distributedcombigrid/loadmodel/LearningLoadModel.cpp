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

  //the durationInformation data type needs to be communicated to MPI in every process using it
  MPI_Datatype createMPIDurationType(){
    durationInformation dI;
    MPI_Aint baseaddr, duraddr, procaddr;
    MPI_Address ( &dI.task_id,   &baseaddr);
    MPI_Address ( &dI.duration,   &duraddr);
    MPI_Address ( &dI.nProcesses, &procaddr); 
    MPI_Datatype duration_datatype;
    MPI_Datatype types[] = { MPI_INT, MPI_UNSIGNED_LONG, MPI_UNSIGNED }; 
    int blocklengths[] = { 1, 1, 1 };

    MPI_Aint displacements[] = { 0, duraddr - baseaddr, procaddr - baseaddr };
    MPI_Type_create_struct(
      3, blocklengths, displacements, types,
      &duration_datatype
    );
    int err = MPI_Type_commit(&duration_datatype);
    assert(!err); // && std::to_string(err));

    //verify that the fields have the expected length
    MPI_Aint extent;
    MPI_Type_extent(MPI_INT, &extent);
    assert (extent == sizeof(int) );
    // assert (extent == duraddr - baseaddr); // not true, padding
    MPI_Type_extent(MPI_UNSIGNED_LONG, &extent);
    assert (extent == procaddr - duraddr);
    assert (extent == sizeof(long unsigned int) );
    MPI_Type_extent(MPI_UNSIGNED, &extent);
    assert (extent == sizeof(uint) );

    //and the whole structure, too
    MPI_Type_extent(duration_datatype, &extent);
    assert (extent == sizeof(durationInformation) );

    return duration_datatype;
  }

  std::string getFilename(const LevelVector& levelVector){
    return "./loaddata_" + toString(levelVector) + ".durations";//TODO which directory
  }
 
} /* namespace combigrid */

