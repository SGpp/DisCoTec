// @author Theresa Pollinger
#ifndef DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_
#define DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_

#include <bitset>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

// https://simsgs.informatik.uni-stuttgart.de:8444/heenemo/combi/issues/39

namespace combigrid {

enum class MultipleGrids {onlyOne, isHigherNeighbor, isLowerNeighbor};

class GridEnumeration{
  // encode the grid number as bits: if bit is set at position d, it means it is the
  // "second" point at this position in dimension d.

  static bool isHigherInDimension(char dimension, unsigned long testnumber){
      std::bitset<sizeof(unsigned long) * 8> b(testnumber);
      return b.test(dimension);
  }

  static unsigned long getLowerNeighborInDimension(char dimension, unsigned long testnumber){
      assert(isHigherInDimension(dimension, testnumber));
      std::bitset<sizeof(unsigned long) * 8> b(testnumber);
      b.set(dimension, false);
      return b.to_ulong();
  }

  static unsigned long getUpperNeighborInDimension(char dimension, unsigned long testnumber){
      assert(!isHigherInDimension(dimension, testnumber));
      std::bitset<sizeof(unsigned long) * 8> b(testnumber);
      b.set(dimension, true);
      return b.to_ulong();
  }
};


class DFGEnsemble{
public:
  using DFG = DistributedFullGrid<CombiDataType>;
  using DFGPointer = std::unique_ptr<DFG>;

  template<class... U>
  DFGEnsemble(int dimensions, U&&... restDFGArgs){
    size_t numGridsInEnsemble = powerOfTwo[dimensions];

    ensemble_ = std::vector<DFGPointer>();
    ensemble_.reserve(numGridsInEnsemble);
    for (size_t i=0; i < numGridsInEnsemble; ++i){
      ensemble_.push_back(DFGPointer(new DFG(dimensions, std::forward<U>(restDFGArgs)...)));
    }
  }

  inline DFG& getDFG(size_t number) { return *(ensemble_[number]); }

  inline const DFG& getDFG(size_t number) const { return *(ensemble_[number]); }

  size_t getNumFullGrids() const { return ensemble_.size(); }

private:
  std::vector<DFGPointer> ensemble_;
}; // class DFGEnsemble

} // namespace combigrid


#endif /* DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_ */
