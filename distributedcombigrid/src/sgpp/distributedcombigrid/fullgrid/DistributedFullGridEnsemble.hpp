// @author Theresa Pollinger
#ifndef DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_
#define DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

// https://simsgs.informatik.uni-stuttgart.de:8444/heenemo/combi/issues/39

class DFGEnsemble{
public:
  using DFG = DistributedFullGrid<CombiDataType>;
  using DFGPointer = std::unique_ptr<DFG>;

  template<class... U>
  DFGEnsemble(int dimensions, U&&... restDFGArgs){
    int numGridsInEnsemble = powerOfTwo[dimensions];

    ensemble_ = std::vector<DFGPointer>();
    ensemble_.reserve(numGridsInEnsemble);
    for (int i=0; i < numGridsInEnsemble; ++i){
      ensemble_.push_back(DFGPointer(new DFG(dimensions, std::forward<U>(restDFGArgs)...)));
    }
  }

  inline DFG& getDFG(int number) { return *(ensemble_[number]); }

  inline const DFG& getDFG(int number) const { return *(ensemble_[number]); }

  const int getNumFullGrids() const { return ensemble_.size(); }

private:
  std::vector<DFGPointer> ensemble_;
}; // class DFGEnsemble


#endif /* DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_ */
