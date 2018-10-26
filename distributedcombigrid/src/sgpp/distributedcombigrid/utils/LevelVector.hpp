/*
 * LevelVector.hpp
 *
 *  Created on: May 14, 2013
 *      Author: heenemo
 */

#ifndef LEVELVECTOR_HPP_
#define LEVELVECTOR_HPP_

// #include <sstream>
// #include <string>
#include "sgpp/distributedcombigrid/utils/IndexVector.hpp"

namespace combigrid {

    typedef IndexVector LevelVector;
    std::string toString(combigrid::LevelVector const& l);
}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
