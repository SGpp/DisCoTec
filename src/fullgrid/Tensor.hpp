// Initial draft by Klaus Iglberger -- thank you!
#pragma once

#include <boost/multi_array.hpp>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

#include "utils/IndexVector.hpp"
#include "utils/Types.hpp"

namespace combigrid {

// Implementation of a non-owning tensor
template <typename Type>
class Tensor : public TensorIndexer {
 public:
  TensorIndexerDim() = default;
  explicit TensorIndexerDim(IndexArray<DIM> extents) : extents_{extents}, localOffsets_{} {
    IndexType nrElements = 1;
    for (DimType j = 0; j < DIM; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
    assert(this->size() == static_cast<size_t>(nrElements));
  }

  // Allow construction from IndexVector for convenience
  explicit TensorIndexerDim(const IndexVector& extents) : localOffsets_{} {
    assert(extents.size() == DIM);
    std::copy_n(extents.begin(), DIM, extents_.begin());
    IndexType nrElements = 1;
    for (DimType j = 0; j < DIM; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
    assert(this->size() == static_cast<size_t>(nrElements));
  }

  TensorIndexerDim(TensorIndexerDim const&) = default;
  TensorIndexerDim(TensorIndexerDim&&) = default;
  TensorIndexerDim& operator=(TensorIndexerDim const&) = default;
  TensorIndexerDim& operator=(TensorIndexerDim&&) = default;

  template <DimType Dimension>
  IndexType dim() const {
    static_assert(Dimension < DIM);
    return extents_[Dimension];
  }

  size_t size() const {
    if constexpr (DIM == 0) {
      return 0;
    }
    auto s =
        std::accumulate(extents_.begin(), extents_.end(), size_t{1}, std::multiplies<size_t>());
    assert(s < static_cast<size_t>(std::numeric_limits<IndexType>::max()));
    return s;
  }

  IndexType sequentialIndex(const IndexArray<DIM>& indexArray) const {
    return std::inner_product(indexArray.begin(), indexArray.end(), localOffsets_.begin(),
                              IndexType{0});
  }

  IndexType sequentialIndex(const IndexVector& indexVector) const {
    assert(indexVector.size() == DIM);
    return std::inner_product(indexVector.begin(), indexVector.end(), localOffsets_.begin(),
                              IndexType{0});
  }

  IndexArray<DIM> getArrayIndex(IndexType index) const {
    IndexArray<DIM> indexArray{};
    for (auto j = DIM; j > 0; --j) {
      auto dim_i = static_cast<DimType>(j - 1);
      indexArray[dim_i] = index / localOffsets_[dim_i];
      index = index % localOffsets_[dim_i];
    }
    return indexArray;
  }

  const IndexVector& getVectorIndex(IndexType index) const {
    static thread_local IndexVector indexVector;
    indexVector.resize(DIM);
    for (auto j = DIM; j > 0; --j) {
      auto dim_i = static_cast<DimType>(j - 1);
      indexVector[dim_i] = index / localOffsets_[dim_i];
      index = index % localOffsets_[dim_i];
    }
    return indexVector;
  }

  const IndexArray<DIM>& getExtentsArray() const { return extents_; }
  const IndexArray<DIM>& getOffsetsArray() const { return localOffsets_; }

  const IndexVector& getExtentsVector() const {
    static thread_local IndexVector v;
    v.assign(extents_.begin(), extents_.end());
    return v;
  }

  const IndexVector& getOffsetsVector() const {
    static thread_local IndexVector v;
    v.assign(localOffsets_.begin(), localOffsets_.end());
    return v;
  }

 protected:
  IndexArray<DIM> extents_{};
  IndexArray<DIM> localOffsets_{};
};

// Compile-time dimensioned non-owning tensor
template <typename Type, DimType DIM>
class TensorDim : public TensorIndexerDim<DIM> {
 public:
  TensorDim() = default;
  explicit TensorDim(Type* data, IndexArray<DIM> extents)
      : TensorIndexerDim<DIM>(extents), data_(data) {}

  explicit TensorDim(Type* data, const IndexVector& extents)
      : TensorIndexerDim<DIM>(extents), data_(data) {}

  TensorDim(TensorDim const&) = delete;
  TensorDim(TensorDim&&) = default;
  TensorDim& operator=(TensorDim const&) = delete;
  TensorDim& operator=(TensorDim&&) = default;

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

  Type& operator()(IndexArray<DIM> index) { return this->operator[](this->sequentialIndex(index)); }
  Type const& operator()(IndexArray<DIM> index) const {
    return this->operator[](this->sequentialIndex(index));
  }

  template <std::size_t NumDims>
  inline boost::multi_array_ref<Type, NumDims> getAsMultiArrayRef() {
    static_assert(NumDims == DIM);
    return boost::multi_array_ref<Type, NumDims>(this->getData(), this->getExtentsVector(),
                                                 boost::fortran_storage_order());
  }

  template <std::size_t NumDims>
  inline boost::const_multi_array_ref<Type, NumDims> getAsConstMultiArrayRef() const {
    static_assert(NumDims == DIM);
    return boost::const_multi_array_ref<Type, NumDims>(this->getData(), this->getExtentsVector(),
                                                       boost::fortran_storage_order());
  }

 private:
  Type* data_ = nullptr;
};

template <typename Type, DimType DIM>
void print(TensorDim<Type, DIM> const& T) {
  if constexpr (DIM == 1) {
    std::cout << '(';
    for (IndexType i = 0; i < T.template dim<0>(); ++i) {
      std::cout << ' ' << T[i];
    }
    std::cout << " )\n\n";
  } else if constexpr (DIM == 2) {
    for (IndexType i = 0; i < T.template dim<0>(); ++i) {
      std::cout << '(';
      for (IndexType j = 0; j < T.template dim<1>(); ++j) {
        std::cout << ' ' << T({i, j});
      }
      std::cout << " )\n";
    }
    std::cout << '\n';
  } else if constexpr (DIM == 3) {
    for (IndexType i = 0; i < T.template dim<0>(); ++i) {
      std::cout << "[\n";
      for (IndexType j = 0; j < T.template dim<1>(); ++j) {
        std::cout << '(';
        for (IndexType k = 0; k < T.template dim<2>(); ++k) {
          std::cout << ' ' << T({i, j, k});
        }
        std::cout << " )\n";
      }
      std::cout << "]\n";
    }
  } else {
    static_assert(DIM <= 3, "print not supported for dimensions > 3");
  }
}

// Dim-templated TensorIndexer: stores extents and offsets as std::array
template <DimType DIM>
class TensorIndexerDim {
 public:
  TensorIndexerDim() = default;
  explicit TensorIndexerDim(IndexArray<DIM> extents) : extents_{extents} {
    IndexType nrElements = 1;
    for (DimType j = 0; j < DIM; j++) {
      localOffsets_[j] = nrElements;
      nrElements = nrElements * extents_[j];
    }
    assert(this->size() == static_cast<size_t>(nrElements));
  }

  TensorIndexerDim(TensorIndexerDim const&) = delete;
  TensorIndexerDim(TensorIndexerDim&&) = default;
  TensorIndexerDim& operator=(TensorIndexerDim const&) = delete;
  TensorIndexerDim& operator=(TensorIndexerDim&&) = default;

  size_t size() const {
    if constexpr (DIM == 0) {
      return 0;
    }
    auto s = std::accumulate(extents_.begin(), extents_.end(), static_cast<size_t>(1),
                             std::multiplies<size_t>());
    assert(s < static_cast<size_t>(std::numeric_limits<IndexType>::max()));
    return s;
  }

  IndexType sequentialIndex(const IndexArray<DIM>& indexVector) const {
    return std::inner_product(indexVector.begin(), indexVector.end(), localOffsets_.begin(),
                              static_cast<IndexType>(0));
  }

  IndexType sequentialIndex(const IndexVector& indexVector) const {
    assert(indexVector.size() == DIM);
    return std::inner_product(indexVector.begin(), indexVector.end(), localOffsets_.begin(),
                              static_cast<IndexType>(0));
  }

  IndexArray<DIM> getVectorIndex(IndexType index) const {
    IndexArray<DIM> indexVector{};
    for (auto j = static_cast<int>(DIM); j > 0; --j) {
      auto dim_i = static_cast<DimType>(j - 1);
      indexVector[dim_i] = index / localOffsets_[dim_i];
      index = index % localOffsets_[dim_i];
    }
    return indexVector;
  }

  const IndexArray<DIM>& getExtents() const { return extents_; }
  const IndexArray<DIM>& getOffsets() const { return localOffsets_; }

  // backward-compatible accessors returning vectors
  IndexVector getExtentsVector() const { return toVector<DIM>(extents_); }
  IndexVector getOffsetsVector() const { return toVector<DIM>(localOffsets_); }

 protected:
  IndexArray<DIM> extents_{};
  IndexArray<DIM> localOffsets_{};
};

// Dim-templated non-owning tensor
template <typename Type, DimType DIM>
class TensorDim : public TensorIndexerDim<DIM> {
 public:
  TensorDim() = default;
  explicit TensorDim(Type* data, IndexArray<DIM> extents)
      : TensorIndexerDim<DIM>(extents), data_(data) {}

  TensorDim(TensorDim const&) = delete;
  TensorDim(TensorDim&&) = default;
  TensorDim& operator=(TensorDim const&) = delete;
  TensorDim& operator=(TensorDim&&) = default;

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

  Type& operator()(const IndexArray<DIM>& index) {
    return this->operator[](this->sequentialIndex(index));
  }
  Type const& operator()(const IndexArray<DIM>& index) const {
    return this->operator[](this->sequentialIndex(index));
  }

  template <std::size_t NumDims>
  inline boost::multi_array_ref<Type, NumDims> getAsMultiArrayRef() {
    static_assert(NumDims == DIM, "getAsMultiArrayRef: NumDims != Dim");
    auto extentsVec = this->getExtentsVector();
    return boost::multi_array_ref<Type, NumDims>(this->getData(), extentsVec,
                                                 boost::fortran_storage_order());
  }

  template <std::size_t NumDims>
  inline boost::const_multi_array_ref<Type, NumDims> getAsConstMultiArrayRef() const {
    static_assert(NumDims == DIM, "getAsConstMultiArrayRef: NumDims != Dim");
    auto extentsVec = this->getExtentsVector();
    return boost::const_multi_array_ref<Type, NumDims>(this->getData(), extentsVec,
                                                       boost::fortran_storage_order());
  }

 private:
  Type* data_ = nullptr;
};

template <typename Type, DimType DIM>
void print(TensorDim<Type, DIM> const& T) {
  if constexpr (DIM == 1) {
    std::cout << '(';
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << ' ' << T[i];
    }
    std::cout << " )\n\n";
  } else if constexpr (DIM == 2) {
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << '(';
      for (IndexType j = 0U; j < T.getExtents()[1]; ++j) {
        std::cout << ' ' << T(IndexArray<2>{i, j});
      }
      std::cout << " )\n";
    }
    std::cout << '\n';
  } else if constexpr (DIM == 3) {
    for (IndexType i = 0U; i < T.getExtents()[0]; ++i) {
      std::cout << "[\n";
      for (IndexType j = 0U; j < T.getExtents()[1]; ++j) {
        std::cout << '(';
        for (IndexType k = 0U; k < T.getExtents()[2]; ++k) {
          std::cout << ' ' << T(IndexArray<3>{i, j, k});
        }
        std::cout << " )\n";
      }
      std::cout << "]\n";
    }
  } else {
    // print nothing for dimensions > 3
  }
}

}  // namespace combigrid
