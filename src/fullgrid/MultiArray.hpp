#ifndef MULTIARRAY_HPP_
#define MULTIARRAY_HPP_

#include "boost/multi_array.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "fullgrid/FullGrid.hpp"

namespace combigrid {
template <typename T, DimType D>
using MultiArray = boost::multi_array<T, D>;

template <typename T, DimType D>
using MultiArrayRef = boost::multi_array_ref<T, D>;

/* create a multiarray ref to the local data in a distributed fullgrid */
template <typename T, size_t D>
static MultiArrayRef<T, D> createMultiArrayRef(DistributedFullGrid<T>& dfg) {
  /* note that the sizes of dfg are in reverse order compared to shape */
  std::vector<size_t> shape(dfg.getLocalExtents().rbegin(), dfg.getLocalExtents().rend());
  if (!reverseOrderingDFGPartitions) {
    assert(false && "this is not adapted to normal ordering of DFG partitions yet");
  }

  return MultiArrayRef<T, D>(dfg.getData(), shape);
}

/* create a multiarray ref to the data in a fullgrid */
template <typename T, size_t D>
static MultiArrayRef<T, D> createMultiArrayRef(FullGrid<T>& fg) {
  assert(fg.isGridCreated());

  /* note that the sizes of dfg are in reverse order compared to shape */
  std::vector<size_t> shape(fg.getSizes().rbegin(), fg.getSizes().rend());
  if (!reverseOrderingDFGPartitions) {
    assert(false && "this is not adapted to normal ordering of DFG partitions yet");
  }

  return MultiArrayRef<T, D>(fg.getData(), shape);
}

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

typedef boost::multi_array<CombiDataType, 2> MultiArray2;
typedef boost::multi_array<CombiDataType, 3> MultiArray3;
typedef boost::multi_array<CombiDataType, 4> MultiArray4;
typedef boost::multi_array<CombiDataType, 5> MultiArray5;
typedef boost::multi_array<CombiDataType, 6> MultiArray6;
typedef boost::multi_array<CombiDataType, 7> MultiArray7;
typedef boost::multi_array<CombiDataType, 8> MultiArray8;
typedef boost::multi_array<CombiDataType, 9> MultiArray9;
typedef boost::multi_array<CombiDataType, 10> MultiArray10;

typedef boost::multi_array_ref<CombiDataType, 2> MultiArrayRef2;
typedef boost::multi_array_ref<CombiDataType, 3> MultiArrayRef3;
typedef boost::multi_array_ref<CombiDataType, 4> MultiArrayRef4;
typedef boost::multi_array_ref<CombiDataType, 5> MultiArrayRef5;
typedef boost::multi_array_ref<CombiDataType, 6> MultiArrayRef6;
typedef boost::multi_array_ref<CombiDataType, 7> MultiArrayRef7;
typedef boost::multi_array_ref<CombiDataType, 8> MultiArrayRef8;
typedef boost::multi_array_ref<CombiDataType, 9> MultiArrayRef9;
typedef boost::multi_array_ref<CombiDataType, 10> MultiArrayRef10;
}

#endif /* MULTIARRAY_HPP_ */
