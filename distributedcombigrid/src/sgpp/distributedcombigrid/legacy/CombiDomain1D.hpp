#ifndef COMBIDOMAIN1D_HPP_
#define COMBIDOMAIN1D_HPP_

#include <sgpp/distributedcombigrid/legacy/AbstractStretchingMaker.hpp>
#include <sgpp/distributedcombigrid/legacy/CombiEquidistantStretching.hpp>
#include <sgpp/distributedcombigrid/legacy/combigrid_utils.hpp>

#include <vector>

namespace combigrid {

/** */
class Domain1D {
 public:
  /** */
  Domain1D(int level, double min, double max, AbstractStretchingMaker& stretching);

  /**
   *  copy constructor!
   */
  Domain1D(const Domain1D& domain);

  virtual ~Domain1D() { ; }

  /** return if the axis is scaled */
  inline bool isAxisScaled() const { return _isStretched; }

  /** return the minimum of the domain */
  inline double getMinDomain() const { return _min; }

  /** return the maximum of the domain */
  inline double getMaxDomain() const { return _max; }

  /** return if the axis is scaled */
  const std::vector<double>& axisScaling() const { return _stretching; }
  /**
   *  @return returns a vector of doubles with the derivative of the coordinate
   *transform evaluated at the corresponding grid points...
   *
   */
  const std::vector<double>& axisJacobian() const { return _jacobian; }

  Stretching getStretchingType() const { return _stretching_type; }

  /** return the level of the domain , then the number of points are (2^L) + 1*/
  inline int getLevel() const { return _level; }

  /** transform the real coordinates to unit coordinates
   * @param coordReal real non-unit coordinates
   * @param coordUnit coordinates in the unit cube
   * @param level_in level of the resolution required
   * @param noBoundary make extrapolation for the boundary cells*/
  void transformRealToUnit(double coordReal, double& coordUnit, int level_in = 0,
                           bool noBoundary = false) const;

  /** transform from unit index to real coordinates
   * @param level input level
   * @param index the index 0..2^level
   * @param realCoord the real coordinate */
  void transformUnitToReal(int level, int index, double& realCoord) const;

  /** flocated the point on one axis
   * @param coordReal the coord on real domain
   * @param level_in input level
   * @param startIndex the left index of the cell
   * @param intersect the intersection of the cell
   */
  void findEntry(double coordReal, int level_in, int& startIndex, double& intersect) const;

  /** returns the mesh width /scaling
   * @param index point index
   * @param level_in the level resolution
   * @param h0 the first mesh width
   * @param h1 the second mesh width */
  void getMeshWidth(int index, int level_in, double& h0, double& h1) const;

 private:
  /** the level in the case of stretched */
  int _level;

  /** if stretching is needed*/
  bool _isStretched;

  /** minimum value of the axis */
  double _min;

  /** minimum value of the axis */
  double _max;

  /** if the axis is scaled then here we store the scaling */
  std::vector<double> _stretching;
  std::vector<double> _jacobian;
  Stretching _stretching_type;
};
}  // namespace combigrid

#endif /* COMBIDOMAIN1D_HPP_ */
