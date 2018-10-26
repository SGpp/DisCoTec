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

    // std::string toString(LevelVector const& l){
    //     std::stringstream ss;
    //     for(size_t i = 0; i < l.size(); ++i)
    //     {
    //         if(i != 0)
    //             ss << ",";
    //         ss << l[i];
    //     }
    //     return ss.str();
    // }
}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
