// Initial draft by Klaus Iglberger -- thank you!
#pragma once

#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include "utils/Types.hpp"

namespace combigrid {

template <typename Type, DimType NumDimensions>
class BaseTensor {
 public:
  explicit BaseTensor(Type* data, std::array<IndexType, NumDimensions>&& extents)
      : data_{data}, extents_{std::move(extents)} {}

  // delete copy and move constructors for now
  BaseTensor(BaseTensor const&) = delete;
  BaseTensor(BaseTensor&&) = delete;
  BaseTensor& operator=(BaseTensor const&) = delete;
  BaseTensor& operator=(BaseTensor&&) = delete;

  std::array<combigrid::IndexType, NumDimensions>& getExtents() { return this->extents_; }
  const std::array<combigrid::IndexType, NumDimensions>& getExtents() const { return this->extents_; }

  Type* getData() { return this->data_; }
  const Type* getData() const { return this->data_; }

  Type& operator[](IndexType a) { return this->data_[a]; }
  const Type& operator[](IndexType a) const { return this->data_[a]; }

 protected:
  Type* data_;
  std::array<combigrid::IndexType, NumDimensions> extents_;
};

// Implementation of a non-owning tensor
template <typename Type, DimType NumDimensions>
class Tensor;

template <typename Type>
class Tensor<Type,1> {
 public:
  explicit Tensor(Type* data, std::array<IndexType, 1>&& extents)
      : base_{BaseTensor<Type, 1>(data, std::move(extents))} {}


  std::array<combigrid::IndexType, 1>& getExtents() { return this->base_.getExtents(); }
  const std::array<combigrid::IndexType, 1>& getExtents() const { return this->base_.getExtents(); }

  Type& operator[](IndexType a) { return this->base_.getData()[a]; }
  Type const& operator[](IndexType a) const { return this->base_.getData()[a]; }

  Type& operator()(IndexType a) { return this->operator[](a); }
  Type const& operator()(IndexType a) const { return this->operator[](a); }

  private:
   BaseTensor<Type, 1> base_;
};

template <typename Type>
class Tensor<Type,2> {
 public:
  explicit Tensor(Type* data, std::array<IndexType, 2>&& extents)
      : base_{BaseTensor<Type, 2>(data, std::move(extents))} {}

  std::array<combigrid::IndexType, 2>& getExtents() { return this->base_.getExtents(); }
  const std::array<combigrid::IndexType, 2>& getExtents() const { return this->base_.getExtents(); }

//   template <DimType Dimension>
//   IndexType dim() const { return this->getExtents()[Dimension]; }

  Type& operator[](IndexType a) { return this->base_.getData()[a]; }
  Type const& operator[](IndexType a) const { return this->base_.getData()[a]; }

  Type& operator()(std::array<combigrid::IndexType, 2> index) { return this->operator[](index[1]*getExtents()[0]+index[0]); }
  Type const& operator()(std::array<combigrid::IndexType, 2>index) const { return this->operator[](index[1]*getExtents()[0]+index[0]); }

  private:
   BaseTensor<Type, 2> base_;
};


template <typename Type>
class Tensor<Type,3> {
 public:
  explicit Tensor(Type* data, std::array<IndexType, 3>&& extents)
      : base_{BaseTensor<Type, 3>(data, std::move(extents))} {}

  std::array<combigrid::IndexType, 3>& getExtents() { return this->base_.getExtents(); }
  const std::array<combigrid::IndexType, 3>& getExtents() const { return this->base_.getExtents(); }

  Type& operator[](IndexType a) { return this->base_.getData()[a]; }
  Type const& operator[](IndexType a) const { return this->base_.getData()[a]; }

  Type& operator()(std::array<combigrid::IndexType, 3> index) { return this->operator[]((index[1])*getExtents()[0]+index[0]); }
  Type const& operator()(std::array<combigrid::IndexType, 3>index) const { return this->operator[]((index[1])*getExtents()[0]+index[0]); }

  private:
   BaseTensor<Type, 3> base_;
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
  }
}

}  // namespace combigrid