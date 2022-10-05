#include "WeibullFaults.hpp"
#include <math.h>
#include <random>
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
namespace combigrid {

WeibullFaults::WeibullFaults(real k, real lambda, int numberOfCombis, bool faultMaster)
    : k_(k),
      lambda_(lambda),
      faultMaster_(faultMaster),
      numberOfCombis_(numberOfCombis),
      t_fault_(-1.0) {}

WeibullFaults::WeibullFaults() : k_(0.7), lambda_(1000), faultMaster_(false), t_fault_(-1.0) {}

WeibullFaults::~WeibullFaults() {
  // TODO Auto-generated destructor stub
}

bool WeibullFaults::failNow(int ncombi, real t_iter, int globalRank) {
  if (!faultMaster_) {
    MASTER_EXCLUSIVE_SECTION { return false; }
  }
  if (ncombi == numberOfCombis_ - 1) {  // do not fault in last iteration
    return false;
  }
  // random probability p
  // Will be used to obtain a seed for the random number engine
  // std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  /*std::uniform_real_distribution<> dis(0, 1);
  real p = dis(gen);
  //t= lambda*(-ln(1-p))^(1/k)
   *   real t_fault = lambda_ * pow(-log(1 - p), 1.0/k_);
   *
  */

  std::chrono::duration<real> dur = std::chrono::high_resolution_clock::now() - startTimeIteration_;
  real t_running = dur.count();
  std::cout << "Sampled time for fault: " << t_fault_ << " simulation time needed " << t_running
            << "\n";
  if (t_fault_ < t_running) {
    // std::cout << "Sampled time for fault: " << t_fault  << " with p: " << p << "\n";

    return true;
  }
  return false;
}

real WeibullFaults::init(std::chrono::high_resolution_clock::time_point startTimeIteration,
                         real t_fault) {
  if (t_fault < 0) {
    std::random_device rd;
    std::weibull_distribution<real> distribution(k_, lambda_);
    t_fault_ = distribution(rd);
    std::cout << "Sampled time for fault: " << t_fault_ << "\n";

  } else {
    t_fault_ = t_fault;
  }
  startTimeIteration_ = startTimeIteration;
  return t_fault_;
}
} /* namespace combigrid */
