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
 
} /* namespace combigrid */

