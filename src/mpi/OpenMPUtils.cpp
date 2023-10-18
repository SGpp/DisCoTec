#include "mpi/OpenMPUtils.hpp"

#include <cassert>

#ifdef _OPENMP
// OpenMP header
#include <omp.h>
#endif

namespace combigrid {

bool OpenMPUtils::isOmpEnabled() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

int OpenMPUtils::getNumThreads() {
  // get the number of used OpenMP threads // gratefully borrowed from PLSSVM
  int numOMPthreads = 1;
#ifdef _OPENMP
#pragma omp parallel default(none) shared(numOMPthreads)
  {
#pragma omp master
    numOMPthreads = omp_get_num_threads();
  }
#endif
  return numOMPthreads;
}

int OpenMPUtils::getMaximumActiveLevels() {
  int maxActiveLevels = 0;
#pragma omp parallel default(none) shared(maxActiveLevels)
#ifdef _OPENMP
  maxActiveLevels = omp_get_max_active_levels();
#endif
  return maxActiveLevels;
}

bool OpenMPUtils::setMaximumActiveLevels(int numberOfLevels) {
  bool setSuccessfully = false;
#pragma omp parallel default(none) shared(setSuccessfully, numberOfLevels)
#ifdef _OPENMP
  omp_set_max_active_levels(numberOfLevels);
  if (OpenMPUtils::getMaximumActiveLevels() == numberOfLevels) {
    setSuccessfully = true;
  }
#endif
  return setSuccessfully;
}

int OpenMPUtils::getNumberOfTeams(int numberOfTeamsToTry) {
  int numTeams = 1;
#ifdef _OPENMP
// this will always be one team prior to OpenMP 5.0
#pragma omp teams num_teams(numberOfTeamsToTry) default(none) shared(numTeams)
  numTeams = omp_get_num_teams();
#endif
  assert(numTeams == 1 && "if this works it's time to further work on the ticket lock algorithm");
  return numTeams;
}

// to be called from within a teams block
int OpenMPUtils::getTeamNumber() {
  int myTeamNumber = 0;
#ifdef _OPENMP
  myTeamNumber = omp_get_team_num();
#endif
  return myTeamNumber;
}

}  // namespace combigrid