/*
 * WeibullFaults.h
 *
 *  Created on: Feb 27, 2017
 *      Author: oberstei
 */

#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_WEIBULLFAULTS_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_WEIBULLFAULTS_HPP_

#include "FaultCriterion.hpp"

namespace combigrid {

class WeibullFaults: public FaultCriterion {
public:
  WeibullFaults(real k, real lambda, int numberOfCombis, bool faultMaster = false);
  WeibullFaults();

  virtual ~WeibullFaults();
  bool failNow(int ncombi, real t_iter, int globalRank);
  real init(std::chrono::high_resolution_clock::time_point  startTimeIteration, real t_fault_);

private:
  friend class boost::serialization::access;
  //shape parameter
  real k_;
  //scale parameter
  real lambda_;
  //indicates if master of each process group is allowed to fail
  bool faultMaster_;
  //total number of combinations
  int numberOfCombis_;
  //failure time
  real t_fault_;
  //starting time of simulation
  std::chrono::high_resolution_clock::time_point  startTimeIteration_;
  //serialize
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar & boost::serialization::base_object<FaultCriterion>(*this);
    ar & k_;
    ar & lambda_;
    ar & faultMaster_;
    ar & numberOfCombis_;
  }

};

} /* namespace combigrid */

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_WEIBULLFAULTS_HPP_ */
