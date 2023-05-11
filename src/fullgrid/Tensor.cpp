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

}  // namespace