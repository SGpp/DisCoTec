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
    for (DimType j = 0; j < NumDimensions; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
    assert(this->size() == nrElements);
  }

  // default copy and move constructors for now
  TensorIndexer(TensorIndexer const&) = default;
  TensorIndexer(TensorIndexer&&) = default;
  TensorIndexer& operator=(TensorIndexer const&) = default;
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
    return std::inner_product(indexArray.begin(), indexArray.end(), this->getOffsetsVector().begin(),
                              0);
  }

  const IndexVector& getExtentsVector() const { return this->extents_; }

  const IndexVector& getOffsetsVector() const { return this->localOffsets_; }

 protected:
  // IndexArray<NumDimensions> extents_{}; 
  // IndexArray<NumDimensions> localOffsets_{};
  //TODO make these arrays (again)
  IndexVector extents_{};
  IndexVector localOffsets_{};
};

// Implementation of a non-owning tensor
template <typename Type, DimType NumDimensions>
class Tensor : public TensorIndexer<NumDimensions> {
 public:
  Tensor() = default;
  explicit Tensor(Type* data, std::array<IndexType, NumDimensions>&& extents)
      : data_(data), TensorIndexer<NumDimensions>(std::move(extents)) {
    if (this->data_ == nullptr) {
      throw std::runtime_error("Data pointer must not be null!");
    }
  }

  // delete copy and move constructors for now
  Tensor(Tensor const&) = delete;
  Tensor(Tensor&&) = delete;
  Tensor& operator=(Tensor const&) = delete;
  Tensor& operator=(Tensor&&) = delete;

  Type* getData() { return this->data_; }
  const Type* getData() const { return this->data_; }

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

}  // namespace combigrid