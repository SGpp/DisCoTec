/*
 * GeneCheckpoint.hpp
 *
 *  Created on: Jul 23, 2014
 *      Author: heenemo
 */

#ifndef GENELOCALCHECKPOINT_HPP_
#define GENELOCALCHECKPOINT_HPP_

#include "GeneGrid.hpp"
#include <vector>

namespace combigrid {


class GeneLocalCheckpoint {
public:
  GeneLocalCheckpoint();

  virtual ~GeneLocalCheckpoint();

  inline GeneComplex* getData();

  // todo: probably size is redundant as size can be determint from sizes or
  // bounds
  void writeCheckpoint( GeneComplex* data, size_t size,
                  std::vector<size_t>& sizes,
                  std::vector<size_t>& bounds );

  void initCheckpoint(size_t size,
                    std::vector<size_t>& sizes,
                    std::vector<size_t>& bounds );

  inline size_t getSize() const;

  inline const std::vector<size_t>& getSizes() const;

  inline const std::vector<size_t>& getBounds() const;

private:
  std::vector<GeneComplex> data_;

  size_t size_;

  std::vector<size_t> sizes_;

  std::vector<size_t> bounds_;
};


inline GeneComplex* GeneLocalCheckpoint::getData(){
  return &data_[0];
}


inline size_t GeneLocalCheckpoint::getSize() const{
  return size_;
}


inline const std::vector<size_t>& GeneLocalCheckpoint::getSizes() const{
  return sizes_;
}


inline const std::vector<size_t>& GeneLocalCheckpoint::getBounds() const{
  return bounds_;
}

} /* namespace combigrid */

#endif /* GENECHECKPOINT_HPP_ */
