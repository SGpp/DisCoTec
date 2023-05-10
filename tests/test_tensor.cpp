#define BOOST_TEST_DYN_LINK

#include "fullgrid/Tensor.hpp"
#include "utils/Types.hpp"
#include "test_helper.hpp"

using namespace combigrid;

void checkPrint() {}

BOOST_FIXTURE_TEST_SUITE(tensor, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(30))

BOOST_AUTO_TEST_CASE(test_print_1d) {
  std::vector<double> v(10U);
  std::array<combigrid::IndexType,1> extents = {10};

  Tensor<double, 1U> t(v.data(), std::move(extents));
  for (IndexType i = 0U; i < t.getExtents()[0]; ++i) {
    t(i) = static_cast<double>(i);
  }

  BOOST_CHECK_NO_THROW(print(t));
}

BOOST_AUTO_TEST_CASE(test_print_2d) {
  std::vector<double> v(50U);

  Tensor<double, 2U> t(v.data(), std::array<IndexType,2>({5U, 10U}));
  for (IndexType i = 0U; i < t.getExtents()[0]; ++i) {
    for (IndexType j = 0U; j < t.getExtents()[1]; ++j) {
      t({i, j}) = static_cast<double>(i * 10 + j);
    }
  }

  BOOST_CHECK_NO_THROW(print(t));
}

BOOST_AUTO_TEST_CASE(test_print_3d) {
  std::vector<double> v(24U);

  Tensor<double, 3U> t(v.data(), std::array<IndexType,3>({2U, 3U, 4U}));
  for (IndexType i = 0U; i < t.getExtents()[0]; ++i) {
    for (IndexType j = 0U; j < t.getExtents()[1]; ++j) {
      for (IndexType k = 0U; k < t.getExtents()[2]; ++k) {
        t({i, j, k}) = static_cast<double>(i * 12U + j * 4U + k);
      }
    }
  }

  BOOST_CHECK_NO_THROW(print(t));
}

BOOST_AUTO_TEST_SUITE_END()
