/*
 * StaticFaults.cpp
 *
 *  Created on: Feb 24, 2017
 *      Author: oberstei
 */

#include "StaticFaults.hpp"

namespace combigrid {

StaticFaults::StaticFaults(FaultsInfo faultsInfo) : faultsInfo_(faultsInfo) {}

StaticFaults::StaticFaults() {}

StaticFaults::~StaticFaults() {}
bool StaticFaults::failNow(int ncombi, real t_iter, int globalRank) {
  FaultsInfo faultsInfo = faultsInfo_;
  IndexVector iF = faultsInfo_.iterationFaults_;
  IndexVector rF = faultsInfo_.globalRankFaults_;

  std::vector<IndexType>::iterator it;
  it = std::find(iF.begin(), iF.end(), ncombi);
  IndexType idx = std::distance(iF.begin(), it);
  // std::cout << "faultInfo" << iF[0] << " " << rF[0] << "\n";
  // Check if current iteration is in iterationFaults_
  while (it != iF.end()) {
    // Check if my rank is the one that fails
    if (globalRank == rF[idx]) return true;
    it = std::find(++it, iF.end(), ncombi);
    idx = std::distance(iF.begin(), it);
  }
  return false;
}
} /* namespace combigrid */
