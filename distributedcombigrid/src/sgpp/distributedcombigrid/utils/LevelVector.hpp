#ifndef LEVELVECTOR_HPP_
#define LEVELVECTOR_HPP_

// #include <sstream>
// #include <string>
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

    typedef IndexVector LevelVector;
    std::string toString(combigrid::LevelVector const& l);
}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
