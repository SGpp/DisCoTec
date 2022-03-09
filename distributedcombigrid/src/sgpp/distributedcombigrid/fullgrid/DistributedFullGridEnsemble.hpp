// @author Theresa Pollinger
#ifndef DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_
#define DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_

#include <bitset>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/sparsegrid/DistributedSparseGridUniform.hpp"

// https://simsgs.informatik.uni-stuttgart.de:8444/heenemo/combi/issues/39

namespace combigrid {

class GridEnumeration{
public:
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

  static size_t getNumberOfHigherDimensions(unsigned long testnumber){
      std::bitset<sizeof(unsigned long) * 8> b(testnumber);
      return b.count();
  }
};


class DFGEnsemble{
public:
  using DFG = DistributedFullGrid<CombiDataType>;
  using DFGPointer = std::unique_ptr<DFG>;

  template<class... U>
  DFGEnsemble(DimType dimensions, U&&... restDFGArgs){
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

  std::vector<size_t> getIndicesOfLowerGridsInDimension(DimType dimension){
    std::vector<size_t> lowerGrids;
    for (size_t i = 0; i < this->getNumFullGrids(); ++i) {
      if (!GridEnumeration::isHigherInDimension(static_cast<char>(dimension), i)) {
        lowerGrids.push_back(i);
      }
    }
    assert(lowerGrids.size() == this->getNumFullGrids()/2);
    return lowerGrids;
  }

  DimType getDimension() const {
    return this->getDFG(0).getDimension();
  }

  std::vector<DFGPointer>::iterator begin() {
    return ensemble_.begin();
  }

  std::vector<DFGPointer>::iterator end() {
    return ensemble_.end();
  }

  /**
   * @brief recursive call to evaluate all neighbor points' contributions to the coordinate (on this
   * part of the grid)
   *
   * @param localIndex the (in-all-dimensions lower) neighbor of coords
   * @param dim the current dimension to split on (start with 0)
   * @param coords the coordinate to interpolate on
   * @param onEnsembleIndices the set of ensemble's dfgs to interpolate on; is bisected on each recursive iteration
   * @return CombiDataType the interpolated value at coords
   */
  CombiDataType evalMultiindexRecursively(const IndexVector& localIndex, DimType dim,
                                          const std::vector<real>& coords,
                                          const std::set<unsigned long>& onEnsembleIndices) const {
    assert(!(dim > this->getDimension()));
    if (dim == this->getDimension()) {
      assert(onEnsembleIndices.size() == 1);
      // std::cout << "eval on " << onEnsembleIndices << " at " << localIndex << std::endl;
      return this->getDFG(*(onEnsembleIndices.cbegin())).evalLocalIndexOn(localIndex, coords);
    } else {
      CombiDataType sum = 0.;
      IndexVector localIndexDimPlusOne = localIndex;
      localIndexDimPlusOne[dim] += 1;
      // std::cout << localIndex << localIndexDimPlusOne << std::endl;
      std::set<unsigned long> onEnsembleIndicesLower, onEnsembleIndicesUpper;
      for (const auto& ensembleIndex : onEnsembleIndices) {
        if (GridEnumeration::isHigherInDimension(static_cast<char>(dim), ensembleIndex)) {
          onEnsembleIndicesLower.insert(ensembleIndex);
        } else {
          onEnsembleIndicesUpper.insert(ensembleIndex);
        }
      }

      sum += evalMultiindexRecursively(localIndex, dim+1, coords, onEnsembleIndicesLower);
      sum += evalMultiindexRecursively(localIndexDimPlusOne, dim+1, coords, onEnsembleIndicesUpper);
      return sum;
    }
  }

  // basically stolen from DistributedFullGrid!
  void eval(const std::vector<real>& coords, CombiDataType& value, MPI_Request* request = nullptr) const {
    assert(coords.size() == this->getDimension());

    // get the lowest-index point of the points
    // whose basis functions contribute to the interpolated value
    auto lowerCoords = this->getDFG(0).getLowerBoundsCoords();
    auto h = this->getDFG(0).getGridSpacing();
    IndexVector localIndexLowerNonzeroNeighborPoint (this->getDimension());
    for (DimType d = 0 ; d < this->getDimension() ; ++d){
      assert(coords[d] >= 0. && coords[d] <= 1.);
      localIndexLowerNonzeroNeighborPoint[d] = static_cast<IndexType>(std::floor((coords[d] - lowerCoords[d]) / h[d]));
    }

    // evaluate at those points and sum up according to the basis function
    // needs to be recursive in order to be dimensionally adaptive
    std::vector<unsigned long> rangeOfEnsembleIndices(this->getNumFullGrids());
    std::iota(rangeOfEnsembleIndices.begin(), rangeOfEnsembleIndices.end(), 0);
    std::set<unsigned long> onEnsembleIndices(rangeOfEnsembleIndices.begin(), rangeOfEnsembleIndices.end());

    value = this->evalMultiindexRecursively(localIndexLowerNonzeroNeighborPoint, 0, coords, onEnsembleIndices);

    if (request == nullptr) {
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, this->getDFG(0).getMPIDatatype(), MPI_SUM, this->getDFG(0).getCommunicator());
    } else {
      MPI_Iallreduce(MPI_IN_PLACE, &value, 1, this->getDFG(0).getMPIDatatype(), MPI_SUM, this->getDFG(0).getCommunicator(), request);
    }
  }

  /** evaluates the ensemble at the specified coordinate
   * @param coords ND coordinates on the unit hypercube [0,1]^D*/
  CombiDataType eval(const std::vector<real>& coords) const {
    CombiDataType value;
    eval(coords, value);
    return value;
  }

  /** evaluates the ensemble on the specified coordinates
   * @param interpolationCoords vector of ND coordinates on the unit hypercube [0,1]^D*/
  std::vector<CombiDataType> getInterpolatedValues(const std::vector<std::vector<real>>& interpolationCoords) const {
    auto numValues = interpolationCoords.size();
    std::vector<CombiDataType> values;
    values.resize(numValues);
    for (size_t i = 0; i < numValues; ++i) {
      this->eval(interpolationCoords[i], values[i]);
    }
    return values;
  }

private:
  std::vector<DFGPointer> ensemble_;
}; // class DFGEnsemble

} // namespace combigrid


#endif /* DISTRIBUTEDCOMBIFULLGRIDENSEMBLE_HPP_ */
