// Initial draft by Klaus Iglberger -- thank you!
#pragma once

#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <variant>
#include <vector>

#include "utils/IndexVector.hpp"
#include "utils/Types.hpp"

namespace combigrid {

template <DimType NumDimensions>
class TensorIndexer {
 public:
  TensorIndexer() = default;
  explicit TensorIndexer(std::array<IndexType, NumDimensions>&& extents)
      : extents_{extents.begin(), extents.end()}, localOffsets_(extents_.size()) {
    // : extents_{std::move(extents)} {
    IndexType nrElements = 1;
    // cf. https://en.wikipedia.org/wiki/Row-_and_column-major_order#Address_calculation_in_general
    // -> column-major order
    for (DimType j = 0; j < NumDimensions; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
    assert(this->size() == nrElements);
  }

  // have only move constructors for now
  TensorIndexer(TensorIndexer const&) = delete;
  TensorIndexer(TensorIndexer&&) = default;
  TensorIndexer& operator=(TensorIndexer const&) = delete;
  TensorIndexer& operator=(TensorIndexer&&) = default;

  template <DimType Dimension>
  IndexType dim() const {
    return this->getExtentsVector()[Dimension];
  }

  // IndexArray<NumDimensions> getExtents() const {
  //   IndexArray<NumDimensions> extentsArray;
  //       std::copy_n(this->getExtents().begin(), NumDimensions, extentsArray.begin());
  //       return extentsArray;
  // }

  // IndexArray<NumDimensions> getOffsets() const {
  //   IndexArray<NumDimensions> offsetsArray;
  //       std::copy_n(this->getOffsets().begin(), NumDimensions, offsetsArray.begin());
  //       return offsetsArray;
  // }

  size_t size() const {
    auto size = std::accumulate(this->extents_.begin(), this->extents_.end(), 1U,
                                std::multiplies<size_t>());
    assert(size < std::numeric_limits<IndexType>::max());
    return size;
  }

  IndexType sequentialIndex(IndexArray<NumDimensions> indexArray) const {
    return std::inner_product(indexArray.begin(), indexArray.end(),
                              this->getOffsetsVector().begin(), 0);
  }

  IndexArray<NumDimensions> getArrayIndex(IndexType index) const {
    IndexArray<NumDimensions> indexArray;
    for (auto j = NumDimensions; j > 0; --j) {
      auto dim_i = static_cast<DimType>(j - 1);
      const auto quotient = index / localOffsets_[dim_i];
      const auto remainder = index % localOffsets_[dim_i];

      indexArray[dim_i] = quotient;
      index = remainder;
    }
    return indexArray;
  }

  const IndexVector& getExtentsVector() const { return this->extents_; }

  const IndexVector& getOffsetsVector() const { return this->localOffsets_; }

 protected:
  // IndexArray<NumDimensions> extents_{};
  // IndexArray<NumDimensions> localOffsets_{};
  // TODO make these arrays (again)
  IndexVector extents_{};
  IndexVector localOffsets_{};
};

// Implementation of a non-owning tensor
template <typename Type, DimType NumDimensions>
class Tensor : public TensorIndexer<NumDimensions> {
 public:
  Tensor() = default;
  explicit Tensor(Type* data, std::array<IndexType, NumDimensions>&& extents)
      : data_(data), TensorIndexer<NumDimensions>(std::move(extents)) {}

  // have only move constructors for now
  Tensor(Tensor const&) = delete;
  Tensor(Tensor&&) = default;
  Tensor& operator=(Tensor const&) = delete;
  Tensor& operator=(Tensor&&) = default;

  Type* getData() {
    if (this->data_ == nullptr) {
      throw std::runtime_error("Data pointer must not be null!");
    }
    return this->data_;
  }

  const Type* getData() const {
    if (this->data_ == nullptr) {
      throw std::runtime_error("Data pointer must not be null!");
    }
    return this->data_;
  }

  void setData(Type* newData) { this->data_ = newData; }

  Type& operator[](IndexType a) { return this->data_[a]; }
  const Type& operator[](IndexType a) const { return this->data_[a]; }

  Type& operator()(IndexArray<NumDimensions> index) {
    return this->operator[](this->sequentialIndex(index));
  }
  Type const& operator()(IndexArray<NumDimensions> index) const {
    return this->operator[](this->sequentialIndex(index));
  }

 private:
  Type* data_ = nullptr;
};

template <typename Type, DimType NumDimensions>
void print(Tensor<Type, NumDimensions> const& T) {
  if constexpr (NumDimensions == 1) {
    std::cout << '(';
    for (IndexType i = 0U; i < T.getExtentsVector()[0]; ++i) {
      std::cout << ' ' << T[i];
    }
    std::cout << " )\n\n";
  } else if constexpr (NumDimensions == 2) {
    for (IndexType i = 0U; i < T.getExtentsVector()[0]; ++i) {
      std::cout << '(';
      for (IndexType j = 0U; j < T.getExtentsVector()[1]; ++j) {
        std::cout << ' ' << T({i, j});
      }
      std::cout << " )\n";
    }
    std::cout << '\n';
  } else if constexpr (NumDimensions == 3) {
    for (IndexType i = 0U; i < T.getExtentsVector()[0]; ++i) {
      std::cout << "[\n";
      for (IndexType j = 0U; j < T.getExtentsVector()[1]; ++j) {
        std::cout << '(';
        for (IndexType k = 0U; k < T.getExtentsVector()[2]; ++k) {
          std::cout << ' ' << T({i, j, k});
        }
        std::cout << " )\n";
      }
      std::cout << "]\n";
    }
  } else {
    static_assert(NumDimensions == 1 || NumDimensions == 2 || NumDimensions == 3,
                  "refusing to print for dimensions > 3!");
  }
}

using SomeTensorIndexer =
    std::variant<TensorIndexer<0>, TensorIndexer<1>, TensorIndexer<2>, TensorIndexer<3>,
                 TensorIndexer<4>, TensorIndexer<5>, TensorIndexer<6>>;

SomeTensorIndexer makeTensorIndexer(IndexVector extents);

template <typename Type>
using SomeTensor = std::variant<Tensor<Type, 0>, Tensor<Type, 1>, Tensor<Type, 2>, Tensor<Type, 3>,
                                Tensor<Type, 4>, Tensor<Type, 5>, Tensor<Type, 6>>;

template <typename Type>
SomeTensor<Type> makeTensor(Type* data, IndexVector extents) {
  SomeTensor<Type> tensor;
  DimType dim = static_cast<DimType>(extents.size());
  switch (dim) {
    case 1: {
      IndexArray<1> extentsArray;
      std::copy_n(extents.begin(), 1, extentsArray.begin());
      tensor = Tensor<Type, 1>(data, std::move(extentsArray));
    } break;
    case 2: {
      IndexArray<2> extentsArray;
      std::copy_n(extents.begin(), 2, extentsArray.begin());
      tensor = Tensor<Type, 2>(data, std::move(extentsArray));
    } break;
    case 3: {
      IndexArray<3> extentsArray;
      std::copy_n(extents.begin(), 3, extentsArray.begin());
      tensor = Tensor<Type, 3>(data, std::move(extentsArray));
    } break;
    case 4: {
      IndexArray<4> extentsArray;
      std::copy_n(extents.begin(), 4, extentsArray.begin());
      tensor = Tensor<Type, 4>(data, std::move(extentsArray));
    } break;
    case 5: {
      IndexArray<5> extentsArray;
      std::copy_n(extents.begin(), 5, extentsArray.begin());
      tensor = Tensor<Type, 5>(data, std::move(extentsArray));
    } break;
    case 6: {
      IndexArray<6> extentsArray;
      std::copy_n(extents.begin(), 6, extentsArray.begin());
      tensor = Tensor<Type, 6>(data, std::move(extentsArray));
    } break;
    default:
      throw std::runtime_error("makeTensor: unsupported dimensionality");
  }
  return tensor;
}

namespace tensor {

template <typename TensorLike>
inline size_t size(const TensorLike& s) {
  size_t size = 0;
  std::visit([&](auto&& arg) { size = arg.size(); }, s);
  assert(size > 0);
  return size;
}

template <typename TensorLike>
inline const IndexVector& getExtents(const TensorLike& s) {
  return std::visit([&](auto&& arg) -> const IndexVector& { return arg.getExtentsVector(); }, s);
}

template <typename TensorLike>
inline const IndexVector& getOffsets(const TensorLike& s) {
  return std::visit([&](auto&& arg) -> const IndexVector& { return arg.getOffsetsVector(); }, s);
}

template <typename Type>
Type* getData(SomeTensor<Type>& s) {
  return std::visit([&](auto&& arg) -> Type* { return arg.getData(); }, s);
}

template <typename Type>
const Type* getData(const SomeTensor<Type>& s) {
  return std::visit([&](auto&& arg) -> const Type* { return arg.getData(); }, s);
}

template <typename Type>
void setData(SomeTensor<Type>& s, Type* data) {
  return std::visit([&](auto&& arg) -> void { return arg.setData(data); }, s);
}
}  // namespace tensor
}  // namespace combigrid