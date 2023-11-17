#define BOOST_TEST_DYN_LINK

#include "fullgrid/Tensor.hpp"
#include "test_helper.hpp"
#include "utils/Types.hpp"

using namespace combigrid;

void checkPrint() {}

BOOST_FIXTURE_TEST_SUITE(tensor, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(30))

BOOST_AUTO_TEST_CASE(test_print_1d) {
  std::vector<double> v(10U);
  IndexVector extents = {10};

  BOOST_TEST_CHECKPOINT("constructing tensor");
  Tensor<double> t(v.data(), std::move(extents));
  BOOST_CHECK_EQUAL(t.size(), v.size());
  BOOST_CHECK_EQUAL(t.getExtentsVector()[0], v.size());
  BOOST_TEST_CHECKPOINT("tensor constructed");

  for (IndexType i = 0U; i < t.getExtentsVector()[0]; ++i) {
    t({i}) = static_cast<double>(i);
  }
  BOOST_TEST_CHECKPOINT("tensor filled");

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    BOOST_CHECK_NO_THROW(print(t));
  }
}

BOOST_AUTO_TEST_CASE(test_print_2d) {
  std::vector<double> v(50U);

  Tensor<double> t(v.data(), IndexVector{5U, 10U});
  for (IndexType i = 0U; i < t.getExtentsVector()[0]; ++i) {
    for (IndexType j = 0U; j < t.getExtentsVector()[1]; ++j) {
      t({i, j}) = static_cast<double>(i * 10 + j);
    }
  }

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    BOOST_CHECK_NO_THROW(print(t));
  }
}

BOOST_AUTO_TEST_CASE(test_print_3d) {
  std::vector<double> v(24U);

  Tensor<double> t(v.data(), IndexVector{2U, 3U, 4U});
  for (IndexType i = 0U; i < t.getExtentsVector()[0]; ++i) {
    for (IndexType j = 0U; j < t.getExtentsVector()[1]; ++j) {
      for (IndexType k = 0U; k < t.getExtentsVector()[2]; ++k) {
        t({i, j, k}) = static_cast<double>(i * 12U + j * 4U + k);
      }
    }
  }
  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    BOOST_CHECK_NO_THROW(print(t));
  }
}

BOOST_AUTO_TEST_CASE(test_6d) {
  std::vector<double> v(5040U, -1.0);
  IndexVector extents = {2U, 3U, 4U, 5U, 6U, 7U};
  auto extentsCopy = extents;
  Tensor<double> t(v.data(), std::move(extentsCopy));
  BOOST_CHECK_EQUAL_COLLECTIONS(t.getExtentsVector().begin(), t.getExtentsVector().end(),
                                extents.begin(), extents.end());
  double increasingValue = 0.0;
  // make sure that Fortran ordering is the default
  for (IndexType i = 0U; i < t.getExtentsVector()[5]; ++i) {
    for (IndexType j = 0U; j < t.getExtentsVector()[4]; ++j) {
      for (IndexType k = 0U; k < t.getExtentsVector()[3]; ++k) {
        for (IndexType l = 0U; l < t.getExtentsVector()[2]; ++l) {
          for (IndexType m = 0U; m < t.getExtentsVector()[1]; ++m) {
            for (IndexType n = 0U; n < t.getExtentsVector()[0]; ++n) {
              t({n, m, l, k, j, i}) = increasingValue;
              increasingValue += 1.0;
            }
          }
        }
      }
    }
  }
  std::vector<double> compare(v.size());
  std::iota(compare.begin(), compare.end(), 0.0);
  BOOST_CHECK_EQUAL_COLLECTIONS(v.begin(), v.end(), compare.begin(), compare.end());
  BOOST_CHECK_EQUAL(t.getExtentsVector()[0], t.dim<0>());
}

BOOST_AUTO_TEST_CASE(test_iterate_lower_boundaries_6d) {
  constexpr DimType dimensionality = 6;
  IndexVector extents = {2U, 3U, 5U, 7U, 11U, 13U};
  auto numElements = std::accumulate(extents.begin(), extents.end(), 1U, std::multiplies<>());
  std::vector<double> v(numElements, -1.0);
  Tensor<double> tensor(v.data(), std::move(extents));

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    for (DimType d = 0; d < dimensionality; ++d) {
      IndexType numberOfPointsVisited = 0;
      const auto strideInThisDimension = tensor.getOffsetsVector()[d];
      const IndexType jump = strideInThisDimension * tensor.getExtentsVector()[d];
      const IndexType numberOfPolesHigherDimensions = static_cast<IndexType>(tensor.size()) / jump;
      IndexType tensorLowestLayerIteratedIndex = 0;

      BOOST_TEST_CHECKPOINT("iterate lowest layer of tensor");
      for (IndexType nHigher = 0; nHigher < numberOfPolesHigherDimensions; ++nHigher) {
        tensorLowestLayerIteratedIndex = nHigher * jump;  // local linear index
        for (IndexType nLower = 0; nLower < tensor.getOffsetsVector()[d];
             ++nLower && ++tensorLowestLayerIteratedIndex && ++numberOfPointsVisited) {
          
          auto arrayIndex = tensor.getVectorIndex(tensorLowestLayerIteratedIndex);
          BOOST_CHECK_EQUAL(arrayIndex[d], 0);
        }
      }
      BOOST_CHECK_EQUAL(numberOfPointsVisited, tensor.size() / tensor.getExtentsVector()[d]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
