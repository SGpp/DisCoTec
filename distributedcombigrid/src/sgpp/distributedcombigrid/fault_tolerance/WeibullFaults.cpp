/*
 * WeibullFaults.cpp
 *
 *  Created on: Feb 27, 2017
 *      Author: oberstei
 */

#include "WeibullFaults.hpp"
#include <math.h>
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include <random>
namespace combigrid {

WeibullFaults::WeibullFaults(real k, real lambda, int numberOfCombis, bool faultMaster): k_(k), lambda_(lambda), numberOfCombis_(numberOfCombis), faultMaster_(faultMaster) {

}

WeibullFaults::WeibullFaults(): k_(0.7), lambda_(1000), faultMaster_(false) {

}

WeibullFaults::~WeibullFaults() {
  // TODO Auto-generated destructor stub
}

bool WeibullFaults::failNow(int ncombi, real t_iter, int globalRank){
  if(!faultMaster_){
    MASTER_EXCLUSIVE_SECTION {
      return false;
    }
  }
  if(ncombi == numberOfCombis_ - 1){ //do not fault in last iteration
    return false;
  }
  //random probability p
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  //std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  /*std::uniform_real_distribution<> dis(0, 1);
  real p = dis(gen);
  //t= lambda*(-ln(1-p))^(1/k)
   *   real t_fault = lambda_ * pow(-log(1 - p), 1.0/k_);
   *
  */
  std::weibull_distribution<real> distribution(k_,lambda_);
  real t_fault = distribution(rd);
  std::cout << "Sampled time for fault: " << t_fault << "\n";
  if(t_fault < t_iter){
    //std::cout << "Sampled time for fault: " << t_fault  << " with p: " << p << "\n";

    return true;
  }
  return false;
}
} /* namespace combigrid */
