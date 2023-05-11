// Initial draft by Klaus Iglberger -- thank you!
#pragma once

#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/preprocessor.hpp>

#include "utils/Types.hpp"
#include "utils/IndexVector.hpp"

namespace combigrid {

template <typename Type, DimType NumDimensions>
class BaseTensor {
 public:
  explicit BaseTensor(Type* data, std::array<IndexType, NumDimensions>&& extents)
      : data_{data}, extents_{std::move(extents)} {
    if (this->data_ == nullptr) {
      throw std::runtime_error("Data pointer must not be null!");
    }
    IndexType nrElements = 1;
    for (DimType j = 0; j < NumDimensions; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
  }

  // delete copy and move constructors for now
  BaseTensor(BaseTensor const&) = delete;
  BaseTensor(BaseTensor&&) = delete;
  BaseTensor& operator=(BaseTensor const&) = delete;
  BaseTensor& operator=(BaseTensor&&) = delete;

  const std::array<combigrid::IndexType, NumDimensions>& getExtents() const {
    return this->extents_;
  }
  const std::array<combigrid::IndexType, NumDimensions>& getLocalOffsets() const {
    return this->localOffsets_;
  }

  Type* getData() { return this->data_; }
  const Type* getData() const { return this->data_; }

  Type& operator[](IndexType a) { return this->data_[a]; }
  const Type& operator[](IndexType a) const { return this->data_[a]; }

  size_t size() const {
    return std::accumulate(this->extents_.begin(), this->extents_.end(), 1U,
                           std::multiplies<size_t>());
  }

 protected:
  Type* data_;
  std::array<combigrid::IndexType, NumDimensions> extents_;
  std::array<combigrid::IndexType, NumDimensions> localOffsets_;
};

// Implementation of a non-owning tensor
template <typename Type, DimType NumDimensions>
class Tensor {
 public:
  explicit Tensor(Type* data, std::array<IndexType, NumDimensions>&& extents)
      : base_{BaseTensor<Type, NumDimensions>(data, std::move(extents))} {}

  const std::array<combigrid::IndexType, NumDimensions>& getExtents() const {
    return this->base_.getExtents();
  }

  //   template <DimType Dimension>
  //   IndexType dim() const { return this->getExtents()[Dimension]; }

  Type& operator[](IndexType a) { return this->base_.getData()[a]; }
  Type const& operator[](IndexType a) const { return this->base_.getData()[a]; }

  Type& operator()(std::array<combigrid::IndexType, NumDimensions> index) {
    return this->operator[](
        std::inner_product(index.begin(), index.end(), this->base_.getLocalOffsets().begin(), 0));
  }
  Type const& operator()(std::array<combigrid::IndexType, NumDimensions> index) const {
    return this->operator[](
        std::inner_product(index.begin(), index.end(), this->base_.getLocalOffsets().begin(), 0));
  }

 private:
  BaseTensor<Type, NumDimensions> base_;
};

template <typename Type, DimType NumDimensions>
void print(Tensor<Type, NumDimensions> const& T) {
  if constexpr (NumDimensions == 1) {
    std::cout << '(';
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << ' ' << T[i];
    }
    std::cout << " )\n\n";
  } else if constexpr (NumDimensions == 2) {
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << '(';
      for (IndexType j = 0U; j < T.getExtents()[1]; ++j) {
        std::cout << ' ' << T({i, j});
      }
      std::cout << " )\n";
    }
    std::cout << '\n';
  } else if constexpr (NumDimensions == 3) {
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << "[\n";
      for (IndexType j = 0U; j < T.getExtents()[1]; ++j) {
        std::cout << '(';
        for (IndexType k = 0U; k < T.getExtents()[2]; ++k) {
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