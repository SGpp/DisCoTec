/*
 * FaultCriterion.h
 *
 *  Created on: Feb 24, 2017
 *      Author: oberstei
 */

#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_FAULTCRITERION_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_FAULTCRITERION_HPP_
#include <iostream>
#include "sgpp/distributedcombigrid/utils/Config.hpp"
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <chrono>

//#include "sgpp/distributedcombigrid/fault_tolerance/StaticFaults.hpp"
//class StaticFaults;


namespace combigrid {
/**
 * This class defines the general properties of a fault Criterion which can be used
 * to determine of a processor should be killed
 */
class FaultCriterion {

public:
  FaultCriterion();
  virtual ~FaultCriterion();
  /**
   * This method decides whether a process should fail.
   * The return value indicates whether the process fails.
   * @param ncombi current combination number
   * @param t_iter time spent for the iteration
   * @param globalRank global rank of the process
   */
  virtual bool failNow(int ncombi, real t_iter, int globalRank){
    std::cout << "Use one of the specialized fault criteria for simulating faults!";
    return false;
  }
  /**
   * This method is used to initialize the fault criterion
   * @param startTimeIteration starting time of process -> used to compare to failure time
   * @param t_fault time a process fails (also referred to as failure time)
   */
  virtual real init(std::chrono::high_resolution_clock::time_point  startTimeIteration, real t_fault){
    return -1.0;
  }



private:
  friend class boost::serialization::access;
  /**
   * This method is used for serialization
   */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version){

  }
};

} /* namespace combigrid */

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_FAULTCRITERION_HPP_ */
