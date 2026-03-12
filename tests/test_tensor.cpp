#define BOOST_TEST_DYN_LINK

#include "fullgrid/Tensor.hpp"
#include "test_helper.hpp"
#include "utils/Types.hpp"

using namespace combigrid;

BOOST_FIXTURE_TEST_SUITE(tensor, TestHelper::BarrierAtEnd, *boost::unit_test::timeout(30))

BOOST_AUTO_TEST_CASE(test_print_1d) {
  std::vector<double> v(10U);
  IndexArray<1> extents = {10};

  BOOST_TEST_CHECKPOINT("constructing tensor");
  TensorDim<double, 1> t(v.data(), extents);
  BOOST_CHECK_EQUAL(t.size(), v.size());
  BOOST_CHECK_EQUAL(t.getExtentsArray()[0], static_cast<IndexType>(v.size()));
  BOOST_TEST_CHECKPOINT("tensor constructed");

  for (IndexType i = 0; i < t.template dim<0>(); ++i) {
    t({i}) = static_cast<double>(i);
  }
  BOOST_TEST_CHECKPOINT("tensor filled");

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    BOOST_CHECK_NO_THROW(print(t));
  }
}

BOOST_AUTO_TEST_CASE(test_print_2d) {
  std::vector<double> v(50U);

  TensorDim<double, 2> t(v.data(), IndexArray<2>{5, 10});
  for (IndexType i = 0; i < t.template dim<0>(); ++i) {
    for (IndexType j = 0; j < t.template dim<1>(); ++j) {
      t({i, j}) = static_cast<double>(i * 10 + j);
    }
  }

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    BOOST_CHECK_NO_THROW(print(t));
  }
}

BOOST_AUTO_TEST_CASE(test_print_3d) {
  std::vector<double> v(24U);

  TensorDim<double, 3> t(v.data(), IndexArray<3>{2, 3, 4});
  for (IndexType i = 0; i < t.template dim<0>(); ++i) {
    for (IndexType j = 0; j < t.template dim<1>(); ++j) {
      for (IndexType k = 0; k < t.template dim<2>(); ++k) {
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
  IndexArray<6> extents = {2, 3, 4, 5, 6, 7};
  TensorDim<double, 6> t(v.data(), extents);
  BOOST_CHECK_EQUAL_COLLECTIONS(t.getExtentsArray().begin(), t.getExtentsArray().end(),
                                extents.begin(), extents.end());
  double increasingValue = 0.0;
  // make sure that Fortran ordering is the default
  for (IndexType i = 0U; i < t.template dim<5>(); ++i) {
    for (IndexType j = 0U; j < t.template dim<4>(); ++j) {
      for (IndexType k = 0U; k < t.template dim<3>(); ++k) {
        for (IndexType l = 0U; l < t.template dim<2>(); ++l) {
          for (IndexType m = 0U; m < t.template dim<1>(); ++m) {
            for (IndexType n = 0U; n < t.template dim<0>(); ++n) {
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
  BOOST_CHECK_EQUAL(t.getExtentsArray()[0], t.template dim<0>());
}

BOOST_AUTO_TEST_CASE(test_iterate_lower_boundaries_6d) {
  constexpr DimType dimensionality = 6;
  IndexArray<dimensionality> extents = {2, 3, 5, 7, 11, 13};
  auto numElements =
      std::accumulate(extents.begin(), extents.end(), IndexType{1}, std::multiplies<>());
  std::vector<double> v(numElements, -1.0);
  TensorDim<double, dimensionality> tensor(v.data(), extents);

  if (TestHelper::getRank(MPI_COMM_WORLD) == 0) {
    for (DimType d = 0; d < dimensionality; ++d) {
      IndexType numberOfPointsVisited = 0;
      const auto strideInThisDimension = tensor.getOffsetsArray()[d];
      const IndexType jump = strideInThisDimension * tensor.getExtentsArray()[d];
      const IndexType numberOfPolesHigherDimensions = static_cast<IndexType>(tensor.size()) / jump;
      IndexType tensorLowestLayerIteratedIndex = 0;

      BOOST_TEST_CHECKPOINT("iterate lowest layer of tensor");
      for (IndexType nHigher = 0; nHigher < numberOfPolesHigherDimensions; ++nHigher) {
        tensorLowestLayerIteratedIndex = nHigher * jump;  // local linear index
        for (IndexType nLower = 0; nLower < tensor.getOffsetsArray()[d];
             ++nLower && ++tensorLowestLayerIteratedIndex && ++numberOfPointsVisited) {
          auto arrayIndex = tensor.getVectorIndex(tensorLowestLayerIteratedIndex);
          BOOST_CHECK_EQUAL(arrayIndex[d], 0);
        }
      }
      BOOST_CHECK_EQUAL(numberOfPointsVisited,
                        static_cast<IndexType>(tensor.size()) / tensor.getExtentsArray()[d]);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
