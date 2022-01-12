#ifndef INDEXVECTOR_HPP_
#define INDEXVECTOR_HPP_

#include <assert.h>
#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include <set>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "sgpp/distributedcombigrid/utils/Types.hpp"

namespace combigrid {

typedef std::vector<IndexType> IndexVector;

inline IndexType sum(const IndexVector& l) {
  IndexType lsum(0);

  for (std::size_t d = 0; d < l.size(); ++d) lsum += l[d];

  return lsum;
}

inline bool operator==(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  // std::cout << "== at line number " << __LINE__ << " in file " << __FILE__ << std::endl;

  for (std::size_t i = 0; i < l1.size(); ++i) {
    if (l1[i] != l2[i]) return false;
  }

  return true;
}

// helper function to compare any vector
template<typename T>
inline bool operator==(const std::vector<T> &u, const std::vector<T> &v) {
  using namespace std;
  if(u.size() != v.size()) {
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

inline IndexVector operator+(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());

  IndexVector tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i) tmp[i] = l1[i] + l2[i];

  return tmp;
}

inline std::ostream& operator<<(std::ostream& os, const IndexVector& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i) os << l[i] << " ";

  os << "]";

  return os;
}

inline IndexVector operator-(const IndexVector& l1, const IndexVector& l2) {
  assert(l1.size() == l2.size());
  //  assert(l1 >= l2);

  IndexVector tmp(l1.size());

  for (std::size_t i = 0; i < l1.size(); ++i) {
    tmp[i] = l1[i] - l2[i];
  }

  return tmp;
}

inline IndexType l1(const IndexVector& l) {
  IndexType lsum(0);

  for (std::size_t d = 0; d < l.size(); ++d) lsum += std::abs(l[d]);

  return lsum;
}

/* read in indexvector from string where ' ' serves as delimiter
 * the vector described by str MUST NOT have a different size than ivec */
inline IndexVector& operator>>(std::string str, IndexVector& ivec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(ivec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i) ivec[i] = boost::lexical_cast<int>(strs[i]);

  return ivec;
}

// helper function to output any vector
template<typename T>
inline std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
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
template<typename T, typename U>
inline std::ostream &operator <<(std::ostream &os, const std::map<T, U> &m) {
  using namespace std;
  for (const auto& any : m) {
    os << "(" << any.first << ") : " << any.second << "; ";
  }
  return os;
}

// helper function to output any set
template<typename T>
inline std::ostream &operator <<(std::ostream &os, const std::set<T> &s) {
  using namespace std;
  os << "[";
  for (const auto& any : s) {
    os << any << ", ";
  }
  os << "]";
  return os;
}

// helper function to output any array
template<typename T, size_t N>
inline std::ostream &operator <<(std::ostream &os, const std::array<T, N> &a) {
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
