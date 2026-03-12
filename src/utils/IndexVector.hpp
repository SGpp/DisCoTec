#ifndef INDEXVECTOR_HPP_
#define INDEXVECTOR_HPP_

#include <assert.h>

#include <algorithm>
#include <array>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <vector>

#include "utils/Types.hpp"

namespace combigrid {

typedef std::vector<IndexType> IndexVector;

template <DimType NumDimensions>
using IndexArray = std::array<IndexType, NumDimensions>;

template <DimType NumDimensions>
using BoundaryArray = std::array<BoundaryType, NumDimensions>;

// vector-to-array and array-to-vector conversions
template <DimType Dim, typename T>
inline std::array<T, Dim> toArray(const std::vector<T>& v) {
  assert(v.size() == Dim);
  std::array<T, Dim> a;
  std::copy_n(v.begin(), Dim, a.begin());
  return a;
}

template <DimType Dim, typename T>
inline std::vector<T> toVector(const std::array<T, Dim>& a) {
  return std::vector<T>(a.begin(), a.end());
}

// array arithmetic operators
template <typename T, size_t N>
inline std::array<T, N> operator+(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  std::array<T, N> tmp;
  for (std::size_t i = 0; i < N; ++i) tmp[i] = static_cast<T>(l1[i] + l2[i]);
  return tmp;
}

template <typename T, size_t N>
inline std::array<T, N> operator-(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  std::array<T, N> tmp;
  for (std::size_t i = 0; i < N; ++i) tmp[i] = static_cast<T>(l1[i] - l2[i]);
  return tmp;
}

// Element-wise comparison functions for arrays (not operators, to avoid
// ambiguity with std::array's lexicographic comparison operators)
template <typename T, size_t N>
inline bool elementwiseLessOrEqual(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  for (std::size_t i = 0; i < N; ++i) {
    if (l1[i] > l2[i]) return false;
  }
  return true;
}

template <typename T, size_t N>
inline bool elementwiseLess(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  for (std::size_t i = 0; i < N; ++i) {
    if (!(l1[i] < l2[i])) return false;
  }
  return true;
}

template <typename T, size_t N>
inline bool elementwiseGreater(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  for (std::size_t i = 0; i < N; ++i) {
    if (!(l1[i] > l2[i])) return false;
  }
  return true;
}

template <typename T, size_t N>
inline bool elementwiseGreaterOrEqual(const std::array<T, N>& l1, const std::array<T, N>& l2) {
  for (std::size_t i = 0; i < N; ++i) {
    if (l1[i] < l2[i]) return false;
  }
  return true;
}

// helper function to compare any vector
template <typename T>
inline bool operator==(const std::vector<T>& u, const std::vector<T>& v) {
  using namespace std;
  for (std::size_t i = 0; i < u.size(); ++i) {
    if (u[i] != v[i]) return false;
  }
  return true;
}

template <typename T>
inline bool equals_with_size_check(const std::vector<T>& u, const std::vector<T>& v) {
  if (u.size() != v.size()) {
    return false;
  }
  for (std::size_t i = 0; i < u.size(); ++i) {
    if (u[i] != v[i]) return false;
  }
  return true;
}

// a l1 < l2 if each entry l1,i < l2,i
inline bool operator<=(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] > l2[i]) return false;
  }

  return true;
}

// a l1 < l2 if each entry l1,i < l2,i
inline bool operator<(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (!(l1[i] < l2[i])) return false;
  }

  return true;
}

// a l1 > l2 if each entry l1,i > l2,i
inline bool operator>(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (!(l1[i] > l2[i])) return false;
  }

  return true;
}

// a l1 >= l2 if each entry l1,i >= l2,i
// todo replace some of these operators by using contrast operators
inline bool operator>=(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] < l2[i]) return false;
  }

  return true;
}

template <typename T>
inline std::vector<T> operator+(const std::vector<T>& l1, const std::vector<T>& l2) {
  assert(l1.size() == l2.size());

  std::vector<T> tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i) tmp[i] = static_cast<T>(l1[i] + l2[i]);

  return tmp;
}

inline std::ostream& operator<<(std::ostream& os, const IndexVector& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i) os << l[i] << " ";

  os << "]";

  return os;
}

template <typename T>
inline std::vector<T> operator-(const std::vector<T>& l1, const std::vector<T>& l2) {
  assert(l1.size() == l2.size());
  //  assert(l1 >= l2);

  std::vector<T> tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    tmp[i] = static_cast<T>(l1[i] - l2[i]);
  }

  return tmp;
}

template <typename T>
inline T l1(const std::vector<T>& l) {
  T lsum(0);

  for (std::size_t d = 0; d < l.size(); ++d) lsum += std::abs(l[d]);

  return lsum;
}

/* read in vector from string where ' ' serves as delimiter
 * the vector described by str MUST NOT have a different size than ivec */
template <typename T>
inline std::vector<T>& operator>>(std::string str, std::vector<T>& ivec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(ivec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i) ivec[i] = boost::lexical_cast<int>(strs[i]);

  return ivec;
}

// helper function to output small types (may otherwise not be printed)
inline std::ostream& operator<<(std::ostream& os, const DimType& v) {
  using namespace std;
  os << static_cast<int>(v);
  return os;
}

// helper function to output any vector
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  using namespace std;
  os << "[";
  // copy(v.begin(), v.end(), ostream_iterator<T>(os, ", "));
  for (const auto& any : v) {
    os << any << " ";
  }
  os << "]";
  return os;
}

// helper function to output any map
template <typename T, typename U>
inline std::ostream& operator<<(std::ostream& os, const std::map<T, U>& m) {
  using namespace std;
  for (const auto& any : m) {
    os << "(" << any.first << ") : " << any.second << "; ";
  }
  return os;
}

// helper function to output any set
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {
  using namespace std;
  os << "[";
  for (const auto& any : s) {
    os << any << ", ";
  }
  os << "]";
  return os;
}

// helper function to output any array
template <typename T, size_t N>
inline std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a) {
  using namespace std;
  os << "[";
  for (const auto& any : a) {
    os << any << ", ";
  }
  os << "]";
  return os;
}

}  // namespace combigrid

#endif /* LEVELVECTOR_HPP_ */
