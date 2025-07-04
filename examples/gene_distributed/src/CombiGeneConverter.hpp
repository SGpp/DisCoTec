/*
 * CombiGeneConverter.hpp
 *
 *  Created on: Jan 14, 2014
 *      Author: heenemo
 */

#ifndef COMBIGENECONVERTER_HPP_
#define COMBIGENECONVERTER_HPP_

#include "GeneGrid.hpp"
#include "../../../include/discotec/fullgrid/FullGrid.hpp"
#include "boost/multi_array.hpp"
#include "../../../include/discotec/utils/Types.hpp"

namespace combigrid {

class CombiGeneConverter {

public:
  CombiGeneConverter();
  virtual ~CombiGeneConverter();

  static void GeneGridToFullGrid( GeneGrid& gg, FullGrid<complex> &fg );

  static void FullGridToGeneGrid( FullGrid<complex> &fg, GeneGrid& gg );
};

} /* namespace combigrid */

#endif /* COMBIGENECONVERTER_HPP_ */
