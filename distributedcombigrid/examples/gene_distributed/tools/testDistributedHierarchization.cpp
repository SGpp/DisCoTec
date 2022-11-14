/* ****************************************************************************
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author mh
#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/hierarchization/Hierarchization.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include <cassert>
#include <complex.h>

using namespace combigrid;

void test2D();
void test5D();
void test6D();

static complex testFunction2D(double x1, double x2) {
  //return complex((-4 * x1 * (x1 - 1)) * (-1 * x2 * (x2 - 1)), 1.0);
  return complex(1.0,1.0);
}

static complex testFunction5D(double x1, double x2, double x3, double x4,
    double x5) {
  //return complex(x1 * x2 * x3 * x4, x2 * x3 * x4 * x5);
  return complex(1.0, 1.0);
}

template<typename T>
bool isequal(T val1, T val2, real eps) {
  bool ret = ( std::abs(val1 - val2) < eps );

  return ret;
}


// helper function to output double vector
inline std::ostream& operator<<(std::ostream& os, const std::vector<double>& l) {
  os << "[";

  for (size_t i = 0; i < l.size(); ++i)
    os << l[i] << " ";

  os << "]";

  return os;
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  //make sure number of mpi procs fits to desired test case
  //test2D();

  test6D();

  //testGENE();

  MPI_Finalize();

  return 0;
}



void test2D() {
  DimType dim = 2;
  LevelVector l = { 2, 3 };
  CommunicatorType comm = MPI_COMM_WORLD;
  std::vector<BoundaryType> boundary(dim, true);
  IndexVector gpd = { 1,2 };

  std::vector<IndexVector> decomposition(dim);
  decomposition[0] = { 0 };
  decomposition[1] = { 0, 4 };

  // create distributed fg and fill with test function
  DistributedFullGrid<complex> dfg(dim, l, comm, boundary, gpd, true, decomposition);
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] = testFunction2D(coords[0], coords[1]);
  }

  // print coordinates
  int myrank = dfg.getMpiRank();
  for (int r = 0; r < dfg.getMpiSize(); ++r) {
    if (r == myrank) {
      IndexVector pcoords(2);
      dfg.getPartitionCoords(pcoords);
      std::cout << "coordinates of rank " << myrank << " " << pcoords
          << std::endl;
      for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
        std::vector<double> coords(dim);
        dfg.getCoordsLocal(li, coords);

        std::cout << li << " (" << coords[0] << "," << coords[1] << ")"
            << std::endl;
      }
    }
    MPI_Barrier(dfg.getCommunicator());
  }

  // create fg and fill with test function
  FullGrid<complex> fg(dim, l, boundary);
  fg.createFullGrid();
  for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);

    fg.getData()[i] = testFunction2D(coords[0], coords[1]);
  }

  // hierarchize fg and distributed fg
  Hierarchization::hierarchize(fg);
  DistributedHierarchization::hierarchize(dfg);

  // compare fg and distributed fg
  /* alt
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);
    assert(dfg.getData()[li] == fg.getData()[gi] && "Hierarchization");
  }*/
    // compare fg and distributed fg
  bool kaputt = false;
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    assert(coords_dfg == coords_fg);

    complex val1 = dfg.getData()[li];
    complex val2 = fg.getData()[gi];


    if( !(isequal<complex>( val1, val2, 1e-3 ) ) ) {
      std::cout << "coords = " << coords_dfg
                << "fg val = " << val2
                << "dfg val = " << val1
                << std::endl;
      kaputt = true;
    }
  }

  assert( !kaputt );

  // dehiarchize fg and distributed fg
  Hierarchization::dehierarchize(fg);
  DistributedHierarchization::dehierarchize(dfg);

  // compare fg and distributed fg to test function
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    assert(coords_dfg == coords_fg);

    //assert( dfg.getData()[li] == fg.getData()[gi] );
    assert(
        isequal(fg.getData()[gi], testFunction2D(coords_fg[0], coords_fg[1]),
            1e-12));
    assert(
        isequal(dfg.getData()[li], testFunction2D(coords_dfg[0], coords_dfg[1]),
            1e-10));

  }
}


void test5D() {
  CommunicatorType comm = MPI_COMM_WORLD;
  DimType dim = 5;
  LevelVector l(dim);
  l[4] = 4;
  l[3] = 4;
  l[2] = 4;
  l[1] = 1;
  l[0] = 2;
  std::vector<BoundaryType> boundary(dim, false);
  boundary[0] = true; //x
  boundary[2] = true; //z
  boundary[3] = true; //v
  boundary[4] = true; //w
  IndexVector gpd(dim, 1);
  gpd[0] = 1;
  gpd[1] = 1;
  gpd[2] = 1;
  gpd[3] = 1;
  gpd[4] = 2;

  // create distributed fg and fill with test function
  DistributedFullGrid<complex> dfg(dim, l, comm, boundary, gpd);
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] = testFunction5D(coords[0], coords[1], coords[2],
        coords[3], coords[4]);
  }

  // create fg and fill with test function
  FullGrid<complex> fg(dim, l, boundary);
  fg.createFullGrid();
  for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);

    fg.getData()[i] = testFunction5D(coords[0], coords[1], coords[2], coords[3],
        coords[4]);
  }

  // hierarchize fg and distributed fg
  Hierarchization::hierarchize(fg);
  DistributedHierarchization::hierarchize(dfg);

  // compare fg and distributed fg
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);
    assert(dfg.getData()[li] == fg.getData()[gi]);
  }

  // dehiarchize fg and distributed fg
  Hierarchization::dehierarchize(fg);
  DistributedHierarchization::dehierarchize(dfg);

  // compare fg and distributed fg to test function
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    assert(coords_dfg == coords_fg);

    //assert( dfg.getData()[li] == fg.getData()[gi] );
    assert(
        isequal(fg.getData()[gi],
            testFunction5D(coords_fg[0], coords_fg[1], coords_fg[2],
                coords_fg[3], coords_fg[4]), 1e-12));
    assert(
        isequal(dfg.getData()[li],
            testFunction5D(coords_dfg[0], coords_dfg[1], coords_dfg[2],
                coords_dfg[3], coords_dfg[4]), 1e-10));
  }
}


void test6D() {
  CommunicatorType comm = MPI_COMM_WORLD;
  DimType dim = 6;

  LevelVector l = { 2, 1, 4, 4, 3, 1 };
  std::vector<BoundaryType> boundary = { true, false, true, true, true, false };
  IndexVector gpd = { 1, 1, 1, 1, 2, 1 };

  std::vector<IndexVector> decomposition(6);
  decomposition[0] = { 0 };
  decomposition[1] = { 0 };
  decomposition[2] = { 0 };
  decomposition[3] = { 0 };
  decomposition[4] = { 0, 4 };
  decomposition[5] = { 0 };

  // create distributed fg and fill with test function
  DistributedFullGrid<complex> dfg(dim, l, comm, boundary, gpd, true, decomposition );
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    std::vector<double> coords(dim);
    dfg.getCoordsLocal(li, coords);

    dfg.getData()[li] = testFunction5D(coords[0], coords[1], coords[2],
        coords[3], coords[4]);
  }

  // create fg and fill with test function
  FullGrid<complex> fg(dim, l, boundary);
  fg.createFullGrid();
  for (size_t i = 0; i < static_cast<size_t>(fg.getNrElements()); ++i) {
    std::vector<double> coords(dim);
    fg.getCoords(i, coords);

    fg.getData()[i] = testFunction5D(coords[0], coords[1], coords[2], coords[3],
        coords[4]);
  }

  // hierarchize fg and distributed fg
  Hierarchization::hierarchize(fg);
  DistributedHierarchization::hierarchize(dfg);

  // compare fg and distributed fg
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    assert(coords_dfg == coords_fg);

    complex val1 = dfg.getData()[li];
    complex val2 = fg.getData()[gi];

    if( !(isequal<complex>( val1, val2, 1e-3 ) ) ) {
      std::cout << "coords = " << coords_dfg
                << "fg val = " << val2
                << "dfg val = " << val1
                << std::endl;
    }
  }

  // dehiarchize fg and distributed fg
  Hierarchization::dehierarchize(fg);
  DistributedHierarchization::dehierarchize(dfg);

  // compare fg and distributed fg to test function
  for (IndexType li = 0; li < dfg.getNrLocalElements(); ++li) {
    IndexType gi = dfg.getGlobalLinearIndex(li);

    std::vector<double> coords_fg(dim);
    fg.getCoords(gi, coords_fg);

    std::vector<double> coords_dfg(dim);
    dfg.getCoordsLocal(li, coords_dfg);

    assert(coords_dfg == coords_fg);

    //assert( dfg.getData()[li] == fg.getData()[gi] );
    assert(
        isequal(fg.getData()[gi],
            testFunction5D(coords_fg[0], coords_fg[1], coords_fg[2],
                coords_fg[3], coords_fg[4]), 1e-12));
    assert(
        isequal(dfg.getData()[li],
            testFunction5D(coords_dfg[0], coords_dfg[1], coords_dfg[2],
                coords_dfg[3], coords_dfg[4]), 1e-10));
  }
}
