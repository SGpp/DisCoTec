#include "fullgrid/Tensor.hpp"
#include "utils/Config.hpp"

// #include <boost/preprocessor.hpp> //TODO consider using this for automatic instantiations

namespace combigrid {

template class TensorIndexer<1>;
template class TensorIndexer<2>;
template class TensorIndexer<3>;
template class TensorIndexer<4>;
template class TensorIndexer<5>;
template class TensorIndexer<6>;

template class Tensor<float, 1>;
template class Tensor<float, 2>;
template class Tensor<float, 3>;
template class Tensor<float, 4>;
template class Tensor<float, 5>;
template class Tensor<float, 6>;

template class Tensor<double, 1>;
template class Tensor<double, 2>;
template class Tensor<double, 3>;
template class Tensor<double, 4>;
template class Tensor<double, 5>;
template class Tensor<double, 6>;

template class Tensor<combigrid::complex, 1>;
template class Tensor<combigrid::complex, 2>;
template class Tensor<combigrid::complex, 3>;
template class Tensor<combigrid::complex, 4>;
template class Tensor<combigrid::complex, 5>;
template class Tensor<combigrid::complex, 6>;

SomeTensorIndexer makeTensorIndexer(IndexVector extents) {
  SomeTensorIndexer indexer;
  DimType dim = static_cast<DimType>(extents.size());
  switch (dim) {
    case 1: {
      IndexArray<1> extentsArray;
      std::copy_n(extents.begin(), 1, extentsArray.begin());
      indexer = TensorIndexer<1>(std::move(extentsArray));
    } break;
    case 2: {
      IndexArray<2> extentsArray;
      std::copy_n(extents.begin(), 2, extentsArray.begin());
      indexer = TensorIndexer<2>(std::move(extentsArray));
    } break;
    case 3: {
      IndexArray<3> extentsArray;
      std::copy_n(extents.begin(), 3, extentsArray.begin());
      indexer = TensorIndexer<3>(std::move(extentsArray));
    } break;
    case 4: {
      IndexArray<4> extentsArray;
      std::copy_n(extents.begin(), 4, extentsArray.begin());
      indexer = TensorIndexer<4>(std::move(extentsArray));
    } break;
    case 5: {
      IndexArray<5> extentsArray;
      std::copy_n(extents.begin(), 5, extentsArray.begin());
      indexer = TensorIndexer<5>(std::move(extentsArray));
    } break;
    case 6: {
      IndexArray<6> extentsArray;
      std::copy_n(extents.begin(), 6, extentsArray.begin());
      indexer = TensorIndexer<6>(std::move(extentsArray));
    } break;
    default:
      throw std::runtime_error("makeTensorIndexer: unsupported dimensionality");
  }
  return indexer;
}

}  // namespace