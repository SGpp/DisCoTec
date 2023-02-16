#pragma once
#include <cassert>

#include "utils/IndexVector.hpp"
#include "utils/Types.hpp"

// template <typename FG_ELEMENT>
// using SliceIterator = typename boost::numeric::ublas::tensor<FG_ELEMENT>::iterator;
namespace combigrid {
template <typename FG_ELEMENT>
class SliceIterator {
 public:
  // cf. https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = FG_ELEMENT;
  using pointer = FG_ELEMENT*;
  using reference = FG_ELEMENT&;

  SliceIterator(const std::vector<int>& subsizes, const std::vector<int>& starts,
                const IndexVector& offsets, std::vector<FG_ELEMENT>& data)
      : currentLocalIndex_(starts),
        subsizes_(subsizes),
        starts_(starts),
        offsets_(offsets),
        dataPointer_(data.data(),
        dataSize_(data.size())) {
    this->setEndIndex();
    this->validateSizes();
  }

  // cheap rule of 5
  explicit SliceIterator() = default;
  SliceIterator(const SliceIterator& other) = delete;
  SliceIterator& operator=(const SliceIterator&) = delete;
  SliceIterator(SliceIterator&& other) = delete;
  SliceIterator& operator=(SliceIterator&& other) = delete;

  reference operator*() const {
#ifndef NDEBUG
    // make sure to only dereference when we actually have the data mapped
    if (linearize(currentLocalIndex_) > dataSize_) {
      std::cout << " subsizes " << subsizes_ << " starts " << starts_ << " end " << endIndex()
                << std::endl;
      std::cout << " ref waah currentLocalIndex_" << currentLocalIndex_ << " linearized "
                << linearize(currentLocalIndex_) << " numElements " << dataSize_
                << std::endl;
    }
    assert(linearize(currentLocalIndex_) <= dataSize_);
    assert(linearize(currentLocalIndex_) < endIndex());
#endif  // ndef NDEBUG
    return dataPointer_[linearize(currentLocalIndex_)];
  }

  pointer operator->() const { return &(dataPointer_[linearize(currentLocalIndex_)]); }
  pointer getPointer() const { return &(dataPointer_[linearize(currentLocalIndex_)]); }

  const std::vector<int>& getVecIndex() const { return currentLocalIndex_; }
  IndexType getIndex() const { return linearize(currentLocalIndex_); }

  const SliceIterator& operator++() {
#ifndef NDEBUG
    assert(linearize(currentLocalIndex_) <= dataSize_);
    auto cpLocalIdx = currentLocalIndex_;
    auto idxbefore = linearize(currentLocalIndex_);
#endif  // ndef NDEBUG
    // increment
    // Fortran ordering
    currentLocalIndex_[0] += 1;
    for (DimType d = 0; d < this->getDimension() - 1; ++d) {
      // wrap around
#ifndef NDEBUG
      if (currentLocalIndex_[d] > starts_[d] + subsizes_[d]) {
        std::cout << " subsizes " << subsizes_ << " starts " << starts_ << " end " << endIndex()
                  << " idxbefore " << idxbefore << std::endl;
        std::cout << "waah currentLocalIndex_" << currentLocalIndex_ << " before " << cpLocalIdx
                  << std::endl;
      }
      assert(currentLocalIndex_[d] >= starts_[d]);
      assert(!(currentLocalIndex_[d] > starts_[d] + subsizes_[d]));
#endif  // ndef NDEBUG
      // expect that it is more likely that we do not have to wrap around
      // if (__builtin_expect((currentLocalIndex_[d] == starts_[d] + subsizes_[d]), false)) {
      // (not sure if this is always the case)
      if (currentLocalIndex_[d] == starts_[d] + subsizes_[d]) {
        currentLocalIndex_[d] = starts_[d];
        ++currentLocalIndex_[d + 1];
      }
    }
    assert(idxbefore < linearize(currentLocalIndex_));
    assert(linearize(currentLocalIndex_) <= endIndex());
    return *this;
  }
  friend bool operator==(const SliceIterator& a, const SliceIterator& b) {
    return a.currentLocalIndex_ == b.currentLocalIndex_;
  };
  friend bool operator!=(const SliceIterator& a, const SliceIterator& b) {
    return a.currentLocalIndex_ != b.currentLocalIndex_;
  };
  pointer begin() { return &(dataPointer_[firstIndex()]); }
  pointer end() { return &(dataPointer_[endIndex()]); }
  bool isAtEnd() const { return (linearize(currentLocalIndex_) == endIndex()); }
  int size() const {
    return std::accumulate(subsizes_.begin(), subsizes_.end(), 1, std::multiplies<int>());
  }

  IndexType firstIndex() const { return linearize(starts_); }
  void setEndIndex() {
    std::vector<int> endVectorIndex = this->starts_;
    endVectorIndex[this->getDimension() - 1] += this->subsizes_[this->getDimension() - 1];
    linearEndIndex_ = linearize(endVectorIndex);
  }
  IndexType endIndex() const {
    assert(linearEndIndex_ != -1);
    return linearEndIndex_;
  }

 protected:
  inline DimType getDimension() const { return static_cast<DimType>(subsizes_.size()); }

  inline IndexType linearize(const std::vector<int>& indexVector) const {
    IndexType index = 0;
    // Fortran ordering
    for (DimType d = 0; d < subsizes_.size(); ++d) {
      index += offsets_[d] * indexVector[d];
    }
    return index;
  }

  void validateSizes() {
    assert(offsets_[0] == 1);
    assert(std::accumulate(this->subsizes_.begin(), this->subsizes_.end(), 1,
                           std::multiplies<int>()) <= this->dataSize_);
    assert(linearize(this->currentLocalIndex_) <= this->dataSize_);
    assert(currentLocalIndex_.size() == this->getDimension());
    assert(subsizes_.size() == this->getDimension());
    assert(starts_.size() == this->getDimension());
  }

  std::vector<int> currentLocalIndex_;
  std::vector<int> subsizes_;
  std::vector<int> starts_;
  IndexVector offsets_;
  FG_ELEMENT* dataPointer_;
  size_t dataSize_;
  IndexType linearEndIndex_ = -1;
};
}  // namespace combigrid