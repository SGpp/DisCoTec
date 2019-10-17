/*
 * hierarchiSeriell.hpp
 *
 *  Created on: 02.07.2014
 *      Author: P. Butz
 */
#ifndef TESTDISTRIBUTEDHIERARCHIZATION_HPP_
#define TESTDISTRIBUTEDHIERARCHIZATION_HPP_

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/hierarchization/DistributedHierarchization.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include <sys/time.h>

namespace combigrid {

/** Class to test distributed hierarchization functionalities
 * run this test with 8 mpi processes to test dim = 1,2,3
 * run with 32 mpi processes for dim = 5 */
class TestDistributedHierarchization {
public:

  static void test_all_cases() {
    int size, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 1D without boundary: ";
      {
        const DimType dim = 1;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        boundary[0] = false;
        gpd[0] = 8;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 1D with boundary: ";
      {
        const DimType dim = 1;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        boundary[0] = true;
        gpd[0] = 8;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base_boundary(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 2D without boundary: ";
      {
        const DimType dim = 2;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        levels[1] = 5;

        boundary[0] = false;
        boundary[1] = false;

        gpd[0] = 2;
        gpd[1] = 4;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 2D with boundary: ";
      {
        const DimType dim = 2;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        levels[1] = 5;

        boundary[0] = true;
        boundary[1] = true;

        gpd[0] = 2;
        gpd[1] = 4;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base_boundary(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 3D without boundary: ";
      {
        const DimType dim = 3;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        levels[1] = 5;
        levels[2] = 3;

        boundary[0] = false;
        boundary[1] = false;
        boundary[2] = false;

        gpd[0] = 2;
        gpd[1] = 2;
        gpd[2] = 2;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 8) {
      if (rank == 0)
        std::cout << "testing 3D with mixed boundaries: ";
      {
        const DimType dim = 3;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        levels[1] = 5;
        levels[2] = 3;

        boundary[0] = true;
        boundary[1] = false;
        boundary[2] = true;

        gpd[0] = 2;
        gpd[1] = 2;
        gpd[2] = 2;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base_boundary(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;
    }

    if (size == 32) {
      if (rank == 0)
        std::cout << "testing 5D with mixed boundaries: ";
      {
        const DimType dim = 5;
        LevelVector levels(dim);
        std::vector<bool> boundary(dim, false);
        std::vector<size_t> gpd(dim);

        levels[0] = 4;
        levels[1] = 4;
        levels[2] = 4;
        levels[3] = 4;
        levels[4] = 4;

        boundary[0] = true;
        boundary[1] = false;
        boundary[2] = true;
        boundary[3] = false;
        boundary[4] = true;

        gpd[0] = 2;
        gpd[1] = 2;
        gpd[2] = 2;
        gpd[3] = 2;
        gpd[4] = 2;

        CommunicatorType comm = MPI_COMM_WORLD;

        test_base_boundary(dim, levels, boundary, gpd, comm);
      }
      if (rank == 0)
        std::cout << "OK" << std::endl;

      theStatsContainer()->save(
          "stats" + boost::lexical_cast<std::string>(rank) + ".dat");
    }
  }
  ;

  // d-dim parabola
  static double testFunction(std::vector<double>& coords) {
    size_t dim(coords.size());

    double res = 1.0;
    for (size_t i = 0; i < dim; ++i) {
      res *= -4.0 * coords[i] * (coords[i] - 1);
    }

    return res;
    //return 1.0;
  }

  // d-dim parabola with offset
  static double testFunctionBoundary(std::vector<double>& coords) {
    size_t dim(coords.size());

    //return res + 1.0;
    return 1.0;
  }

  static double referenceSolutionParabola(LevelVector& levels) {
    const DimType dim = levels.size();
    LevelType lsum = 0;
    for (size_t i = 0; i < levels.size(); ++i)
      lsum += levels[i];

    return std::pow(0.25, lsum - dim);
  }

  // only works properly if ALT_LEVEL_VECTOR in FullGrid.hpp not defined
  static void checkReferenceSolution(DistributedFullGrid<double>& dfg) {
    for (IndexType i = 0; i < dfg.getNrLocalElements(); ++i) {
      LevelVector level(dfg.getDimension());
      LevelVector index(dfg.getDimension());
      IndexType globalIdx = dfg.getGlobalLinearIndex(i);
      dfg.getGlobalLI(globalIdx, level, index);

      double coeff = 1.0;
      for (int k = 0; k < dfg.getDimension(); ++k)
        if (level[k] > 1)
          coeff *= 0.5;

      assert(
          dfg.getElementVector()[i] == 0.0
              || dfg.getElementVector()[i] == coeff);
    }
  }

  static void checkReferenceSolutionParabola(DistributedFullGrid<double>& dfg) {
    for (IndexType i = 0; i < dfg.getNrLocalElements(); ++i) {
      LevelVector level(dfg.getDimension());
      LevelVector index(dfg.getDimension());
      IndexType globalIdx = dfg.getGlobalLinearIndex(i);
      dfg.getGlobalLI(globalIdx, level, index);

      COMBIGRID_ERROR_TEST_EQUAL(dfg.getElementVector()[i],
          referenceSolutionParabola(level), 1e-12, "check parabola");
    }
  }

  static void test_base(DimType dim, LevelVector& levels,
      std::vector<bool>& boundary, std::vector<size_t>& gpd,
      CommunicatorType comm) {
    assert(levels.size() == dim);
    assert(boundary.size() == dim);
    assert(gpd.size() == dim);
    size_t size = 1;
    for (size_t i = 0; i < gpd.size(); ++i)
      size *= gpd[i];

    DistributedFullGrid<double> dfg(dim, levels, comm, boundary, gpd);

    assert(size == dfg.getCommunicatorSize());

    // initialize with testfunction
    std::vector<double> coords(dim, 0.0);
    for (int i = 0; i < dfg.getNrLocalElements(); i++) {
      dfg.getCoordsLocal(i, coords);
      dfg.getElementVector()[i] = testFunction(coords);
    }

    DistributedHierarchization::hierarchize<double>(dfg);

    // first check: compare to reference solution
    checkReferenceSolutionParabola(dfg);

    /*
     Hierarchization::dehierarchize<double>(fg1);

     // second check: compare to test function
     for (int i = 0 ; i < dfg.getNrElements() ; i++) {
     dfg.getCoords( i , coords );
     COMBIGRID_ERROR_TEST_EQUAL( fg1.getElementVector()[i],
     testFunction(coords) , 1e-12 , "fg1" );
     }
     */

    // print data
    /*
     {
     int rank = dfg.getRank();
     int size = dfg.getCommunicatorSize();
     for( int r=0; r < size; ++r ){
     if( r == rank ){
     std::cout << "rank " << r << ":" << std::endl;
     std::cout << dfg;
     }

     MPI_Barrier(comm);
     }
     }
     */
  }

  static void test_base_boundary(DimType dim, LevelVector& levels,
      std::vector<bool>& boundary, std::vector<size_t>& gpd,
      CommunicatorType comm) {
    assert(levels.size() == dim);
    assert(boundary.size() == dim);
    assert(gpd.size() == dim);
    size_t size = 1;
    for (size_t i = 0; i < gpd.size(); ++i)
      size *= gpd[i];

    DistributedFullGrid<double> dfg(dim, levels, comm, boundary, gpd);

    assert(size == dfg.getCommunicatorSize());

    // initialize with testfunction
    std::vector<double> coords(dim, 0.0);
    for (int i = 0; i < dfg.getNrLocalElements(); i++) {
      dfg.getCoordsLocal(i, coords);
      dfg.getElementVector()[i] = testFunctionBoundary(coords);
    }

    DistributedHierarchization::hierarchize<double>(dfg);

    // first check: compare to reference solution
    checkReferenceSolution(dfg);

    /*
     Hierarchization::dehierarchize<double>(fg1);

     // second check: compare to test function
     for (int i = 0 ; i < fg1.getNrElements() ; i++) {
     fg1.getCoords( i , coords );
     COMBIGRID_ERROR_TEST_EQUAL( fg1.getElementVector()[i],
     testFunctionBoundary(coords) , 1e-12 , "fg1" );
     COMBIGRID_ERROR_TEST_EQUAL( fg2.getElementVector()[i],
     testFunctionBoundary(coords) , 1e-12 , "fg2" );
     }
     */

    // print data
    /*
     {
     int rank = dfg.getRank();
     int size = dfg.getCommunicatorSize();
     for( int r=0; r < size; ++r ){
     if( r == rank ){
     std::cout << "rank " << r << ":" << std::endl;
     std::cout << dfg;
     }

     MPI_Barrier(comm);
     }
     }
     */
  }
};
}

#endif /* TESTFULLGRID_HPP_ */
