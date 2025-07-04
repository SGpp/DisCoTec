/*
 * GeneGrid.hpp
 *
 *  Created on: Jan 14, 2014
 *      Author: heenemo
 */

#ifndef GENEGRID_HPP_
#define GENEGRID_HPP_

#include <vector>
#include "../../../include/discotec/utils/Types.hpp"
#include "boost/multi_array.hpp"

namespace combigrid {

struct GeneComplex{
  real r;
  real i;
};

typedef boost::multi_array<GeneComplex, 6> GeneGrid;
typedef boost::multi_array_ref<GeneComplex, 6> GeneGridRef;

} /* namespace combigrid */

#endif /* GENEGRID_HPP_ */
