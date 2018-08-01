/*
 * Types.hpp
 *
 *  Created on: May 14, 2013
 *      Author: heenemo
 */
#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <mpi.h>
#include <cassert>
#include <stdexcept>
#include "sgpp/distributedcombigrid/utils/Config.hpp"

/* only change these types of if you know what you're doing! */

namespace combigrid {

// IndexType must be signed! because we use -1 in some functions as a return
// value. and large enough so that all grid points can be
// numbered. i.e levelsum (with boundary) < 30 for int32 and < 62 for int64
typedef int64_t IndexType;

typedef IndexType LevelType;

typedef size_t DimType;

typedef MPI_Comm CommunicatorType;

typedef int RankType;
}

namespace abstraction {

typedef enum {
  type_unknown = 0,
  type_float,
  type_double,
  type_double_complex,
  type_float_complex
} DataType;

template <class T>
DataType getabstractionDataType() {
  throw new std::invalid_argument("Datatype is not supported!");
}

template <>
inline DataType getabstractionDataType<float>() {
  return abstraction::type_float;
}

template <>
inline DataType getabstractionDataType<double>() {
  return abstraction::type_double;
}

template <>
inline DataType getabstractionDataType<std::complex<double> >() {
  return abstraction::type_double_complex;
}

template <>
inline DataType getabstractionDataType<std::complex<float> >() {
  return abstraction::type_float_complex;
}

inline MPI_Datatype getMPIDatatype(abstraction::DataType type) {
  switch (type) {
    case abstraction::type_float:
      return MPI_FLOAT;

    case abstraction::type_double:
      return MPI_DOUBLE;

    case abstraction::type_double_complex:
      return MPI_DOUBLE_COMPLEX;

    case abstraction::type_float_complex:
      return MPI_COMPLEX;

    case abstraction::type_unknown:
      throw new std::invalid_argument("Type unknown ConvertType!");

    default:
      assert(false && "you should never get here");
  };

  throw new std::invalid_argument("MPI_Datatype Convert(abstraction::DataType) failed!");
}
}

#endif /* TYPES_HPP_ */
