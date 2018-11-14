/*
 * GeneCheckpoint.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: heenemo
 */

#include "GeneLocalCheckpoint.hpp"
#include <cstring> //memcpy

namespace combigrid {

void GeneLocalCheckpoint::writeCheckpoint( GeneComplex* data, size_t size,
                  std::vector<size_t>& sizes,
                  std::vector<size_t>& bounds )
{
  size_ = size;
  data_.resize(size_);
  sizes_ = sizes;
  bounds_ = bounds;

  // copy data
  memcpy( &data_[0], data, size*sizeof(GeneComplex));
}


void GeneLocalCheckpoint::initCheckpoint(size_t size,
                  std::vector<size_t>& sizes,
                  std::vector<size_t>& bounds ){
  size_ = size;
  data_.resize(size_);
  sizes_ = sizes;
  bounds_ = bounds;
}

GeneLocalCheckpoint::GeneLocalCheckpoint() :
    size_(0)
{
}


GeneLocalCheckpoint::~GeneLocalCheckpoint() {
  // TODO Auto-generated destructor stub
}

} /* namespace combigrid */
