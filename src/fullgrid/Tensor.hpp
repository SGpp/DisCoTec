// Initial draft by Klaus Iglberger -- thank you!
#pragma once

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <variant>
#include <vector>

#include "utils/IndexVector.hpp"
#include "utils/Types.hpp"

namespace combigrid {

class TensorIndexer {
 public:
  TensorIndexer() = default;
  explicit TensorIndexer(IndexVector&& extents)
      : extents_{std::move(extents)}, localOffsets_(extents_.size()) {
    IndexType nrElements = 1;
    // cf. https://en.wikipedia.org/wiki/Row-_and_column-major_order#Address_calculation_in_general
    // -> column-major order
    for (DimType j = 0; j < extents_.size(); j++) {
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

  IndexType sequentialIndex(IndexVector indexVector) const {
    return std::inner_product(indexVector.begin(), indexVector.end(),
                              this->getOffsetsVector().begin(), 0);
  }

  IndexVector getVectorIndex(IndexType index) const {
    IndexVector indexVector(this->extents_.size());
    for (auto j = extents_.size(); j > 0; --j) {
      auto dim_i = static_cast<DimType>(j - 1);
      const auto quotient = index / localOffsets_[dim_i];
      const auto remainder = index % localOffsets_[dim_i];

      indexVector[dim_i] = quotient;
      index = remainder;
    }
    return indexVector;
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
template <typename Type>
class Tensor : public TensorIndexer {
 public:
  Tensor() = default;
  explicit Tensor(Type* data, IndexVector&& extents)
      : data_(data), TensorIndexer(std::move(extents)) {}

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

  Type& operator()(IndexVector index) { return this->operator[](this->sequentialIndex(index)); }
  Type const& operator()(IndexVector index) const {
    return this->operator[](this->sequentialIndex(index));
  }

 private:
  Type* data_ = nullptr;
};

template <typename Type>
void print(Tensor<Type> const& T) {
  // if constexpr (NumDimensions == 1) {
  DimType numDimensions = static_cast<DimType>(T.getExtentsVector().size());
  if (numDimensions == 1) {
    std::cout << '(';
    for (IndexType i = 0U; i < T.getExtentsVector()[0]; ++i) {
      std::cout << ' ' << T[i];
    }
    std::cout << " )\n\n";
  } else if (numDimensions == 2) {
    for (IndexType i = 0U; i < T.getExtentsVector()[0]; ++i) {
      std::cout << '(';
      for (IndexType j = 0U; j < T.getExtentsVector()[1]; ++j) {
        std::cout << ' ' << T({i, j});
      }
      std::cout << " )\n";
    }
    std::cout << '\n';
  } else if (numDimensions == 3) {
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
    // static_
    assert(numDimensions == 1 || numDimensions == 2 || numDimensions == 3);
    //  "refusing to print for dimensions > 3!");
  }
}

}  // namespace combigrid