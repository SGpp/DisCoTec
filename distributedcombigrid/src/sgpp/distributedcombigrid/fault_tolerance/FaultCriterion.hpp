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

class FaultCriterion {

public:
  FaultCriterion();
  virtual ~FaultCriterion();
  virtual bool failNow(int ncombi, real t_iter, int globalRank){
    std::cout << "Use one of the specialized fault criteria for simulating faults!";
    return false;
  }

  virtual real init(std::chrono::high_resolution_clock::time_point  startTimeIteration, real t_fault){
    return -1.0;
  }



private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version){

  }
};

} /* namespace combigrid */

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_FAULTCRITERION_HPP_ */
