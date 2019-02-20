/*
 * StaticFaults.h
 *
 *  Created on: Feb 24, 2017
 *      Author: oberstei
 */

#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_STATICFAULTS_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_STATICFAULTS_HPP_

#include "FaultCriterion.hpp"
#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"

namespace combigrid {

class StaticFaults : public FaultCriterion {
 public:
  StaticFaults(FaultsInfo faultsInfo);
  StaticFaults();
  virtual ~StaticFaults();
  bool failNow(int ncombi, real t_iter, int globalRank);

 private:
  friend class boost::serialization::access;

  FaultsInfo faultsInfo_;
  // serialize
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& boost::serialization::base_object<FaultCriterion>(*this);
    ar& faultsInfo_;
  }
};

} /* namespace combigrid */

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_FAULT_TOLERANCE_STATICFAULTS_HPP_ */
