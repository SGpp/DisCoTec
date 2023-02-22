#ifndef TYPES_HPP_
#define TYPES_HPP_

// to resolve https://github.com/open-mpi/ompi/issues/5157
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#include <cassert>
#include <vector>
#include <stdexcept>
#include "utils/Config.hpp"
#include <stdint.h>
#include <limits.h>
#include <boost/cstdint.hpp>

#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "size_t is some strange datatype"
#endif

/* only change these types of if you know what you're doing! */

namespace combigrid {

  // IndexType must be signed! because we use -1 in some functions as a return
  // value. and large enough so that all grid points can be
  // numbered. i.e levelsum (with boundary) < 30 for int32 and < 62 for int64
  typedef int64_t IndexType;

  typedef IndexType LevelType;

  typedef uint8_t DimType;

  // type used to denote how many boundary points there are in a dimension
  // 0 -> no boundary points;
  // 1 -> boundary point on one side (assume the "left" side with coordinate 0);
  // 2 -> boundary points on each side
  // may be changed to enum class if needed
  // (eg to distinguish one-sided left and right boundaries,
  // or to implement GENE's twisted boundary conditions into the hierarchization operator
  // directly...)
  typedef uint8_t BoundaryType;

  typedef MPI_Comm CommunicatorType;

  typedef int RankType;

  // type used for widely-distributed reduction of subspace sizes
  // it is easily enough to fit the largest subspace (19,1,1,1,1,1) in the current scenario
  // (= 2^19 * 3 * 3 * 3 * 3 * 3 = 2^19 * 3^5 = 127401984)
  typedef uint32_t SubspaceSizeType;
  }  // namespace combigrid

namespace abstraction {

typedef enum {
  type_unknown = 0,
  type_char,
  type_float,
  type_double,
  type_long_double,
  type_double_complex,
  type_float_complex,
  type_long_long,
  type_long,
  type_size_t,
  type_uint32_t
} DataType;

template <class T>
DataType getabstractionDataType() {
  throw new std::invalid_argument("Datatype is not supported!");
}

template <>
inline DataType getabstractionDataType<char>() {
  return abstraction::type_char;
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
inline DataType getabstractionDataType<long double>() {
  return abstraction::type_long_double;
}

template <>
inline DataType getabstractionDataType<std::complex<double> >() {
  return abstraction::type_double_complex;
}

template <>
inline DataType getabstractionDataType<std::complex<float> >() {
  return abstraction::type_float_complex;
}

template <>
inline DataType getabstractionDataType<long long>() {
  return abstraction::type_long_long;
}

template <>
inline DataType getabstractionDataType<long>() {
  return abstraction::type_long;
}

template <>
inline DataType getabstractionDataType<size_t>() {
  return abstraction::type_size_t;
}

template <>
inline DataType getabstractionDataType<uint32_t>() {
  return abstraction::type_uint32_t;
}

inline MPI_Datatype getMPIDatatype(abstraction::DataType type) {
  switch (type) {
    case abstraction::type_char:
      return MPI_CHAR;

    case abstraction::type_float:
      return MPI_FLOAT;

    case abstraction::type_double:
      return MPI_DOUBLE;

    case abstraction::type_long_double:
      return MPI_LONG_DOUBLE;

    case abstraction::type_double_complex:
      return MPI_DOUBLE_COMPLEX;

    case abstraction::type_float_complex:
      return MPI_COMPLEX;

    case abstraction::type_long_long:
      return MPI_LONG_LONG;

    case abstraction::type_long:
      return MPI_LONG;

    case abstraction::type_size_t:
      return MPI_SIZE_T;

    case abstraction::type_uint32_t:
      return MPI_UINT32_T;

    case abstraction::type_unknown:
      throw new std::invalid_argument("Type unknown ConvertType!");

    default:
      assert(false && "you should never get here");
      throw new std::invalid_argument("Type unknown ConvertType!");
  };

  throw new std::invalid_argument("MPI_Datatype Convert(abstraction::DataType) failed!");
}
  template <typename T>
  inline std::string toString(std::vector<T> const& v){
      std::stringstream ss;
      ss << "[[";
      for(size_t i = 0; i < v.size(); ++i)
      {
          if(i != 0)
              ss << ",";
          ss << v[i];
      }
      ss << "]]";
      return ss.str();
  }
}  // namespace abstraction

#endif /* TYPES_HPP_ */
