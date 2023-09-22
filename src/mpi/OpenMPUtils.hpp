#pragma once

namespace combigrid {

namespace OpenMPUtils {

bool isOmpEnabled();

int getNumThreads();

int getMaximumActiveLevels();

[[nodiscard]] bool setMaximumActiveLevels(int numberOfLevels);

}  // namespace OpenMPUtils

}  // namespace combigrid