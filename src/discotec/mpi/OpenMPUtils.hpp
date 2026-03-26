#pragma once

namespace combigrid {

namespace OpenMPUtils {

bool isOmpEnabled();

int getNumThreads();

int getMaximumActiveLevels();

[[nodiscard]] bool setMaximumActiveLevels(int numberOfLevels);

int getNumberOfTeams(int numberOfTeamsToTry);

int getTeamNumber();

}  // namespace OpenMPUtils

}  // namespace combigrid