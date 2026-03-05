#ifndef MULTIARRAY_HPP_
#define MULTIARRAY_HPP_

#include "boost/multi_array.hpp"
#include "utils/Types.hpp"

namespace combigrid {
template <typename T, DimType D>
using MultiArray = boost::multi_array<T, D>;

template <typename T, DimType D>
using MultiArrayRef = boost::multi_array_ref<T, D>;

/* create a general multiarray ref.
 *
 * shape is given in the same order x_d, ..., x_2, x_1 as the usual
 * multiarrary access notation, i.e. [index_d][...][index_2][index_1]
 *
 * note that this is the reverse order as e.g. the sizes of FullGrid or
 * DistributedFullGrid are given */
template <typename T, DimType D>
static MultiArrayRef<T, D> createMultiArrayRef(T* data, IndexVector shape) {
  return MultiArrayRef<T, D>(data, shape);
}

/* create a general multiarray.
 *
 * shape is given in the same order x_d, ..., x_2, x_1 as the usual
 * multiarrary access notation, i.e. [index_d][...][index_2][index_1]
 *
 * note that this is the reverse order as e.g. the sizes of FullGrid or
 * DistributedFullGrid are given */
template <typename T, DimType D>
static MultiArray<T, D> createMultiArray(IndexVector shape) {
  return MultiArray<T, D>(shape);
}
}  // namespace combigrid

#endif /* MULTIARRAY_HPP_ */
