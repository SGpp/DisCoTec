/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKENSEMPLE_HPP_
#define TASKENSEMPLE_HPP_

#include <deal.II/base/timer.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <hyper.deal.combi/include/functionalities/dynamic_convergence_table.h>
#include <hyper.deal.combi/include/functionalities/vector_dummy.h>
#include <math.h>

#include <cassert>
#include <sstream>

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridEnsemble.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"
//#include <hyper.deal.combi/applications/advection_reference_dealii/include/application.h>

namespace combigrid {

class TaskEnsemble : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskEnsemble(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
               LoadModel* loadModel, std::string filename, real dt, bool makeExactSolution,
               IndexVector p = IndexVector(0),
               FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit) {
    dt_ = (dt), p_ = (p), make_exact_ = (makeExactSolution), initialized_ = (false),
    stepsTotal_ = (0), dfgEnsemble_ = (nullptr), _filename = (filename);
  }

  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);
    assert(dfgEnsemble_ == nullptr);

    int lrank;
    MPI_Comm_rank(lcomm, &lrank);

    /* create distributed full grid. we try to find a balanced ratio between
     * the number of grid points and the number of processes per dimension
     * by this very simple algorithm. to keep things simple we require powers
     * of two for the number of processes here. */
    int np;
    MPI_Comm_size(lcomm, &np);

    // check if power of two
    if (!((np > 0) && ((np & (~np + 1)) == np)))
      assert(false && "number of processes not power of two");

    DimType dim = this->getDim();
    IndexVector p(dim, 1);
    const LevelVector& l = this->getLevelVector();

    if (p_.size() == 0) {
      // compute domain decomposition
      IndexType prod_p(1);

      while (prod_p != static_cast<IndexType>(np)) {
        DimType dimMaxRatio = 0;
        real maxRatio = 0.0;

        for (DimType k = 0; k < dim; ++k) {
          real ratio = std::pow(2.0, l[k]) / p[k];

          if (ratio > maxRatio) {
            maxRatio = ratio;
            dimMaxRatio = k;
          }
        }

        p[dimMaxRatio] *= 2;
        prod_p = 1;

        for (DimType k = 0; k < dim; ++k) prod_p *= p[k];
      }
    } else {
      p = p_;
    }

    // create local subgrid on each process
    dfgEnsemble_ = new DFGEnsemble(dim, l, lcomm, this->getBoundary(), p);

    // for visualization of grids in Paraview
    Point<Problem::dim_> offset;
    for (unsigned int i = 0; i < dim; i++) offset[i] = l[i] * 1.1;

    // if(!make_exact_){
    this->problem = std::make_shared<Problem>(lcomm, table, offset);
    this->problem->reinit(_filename);

    std::vector<std::array<Number, Problem::dim_ + 2>> coords_dealii = problem->get_result();
    dealii_size = coords_dealii.size();

    // array containing the coordinates of all grids.
    std::vector<std::vector<std::array<Number, Problem::dim_>>> element_coords(
        dfgEnsemble_->getNumFullGrids());

    // loop over each fullgrid to set coordinates of each point in element_coords.
    // furthermore the DisCoTec grids are initialized with NaN
    for (unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); ++i) {
      auto& dfg = dfgEnsemble_->getDFG(i);
      std::vector<CombiDataType>& elements = dfg.getElementVector();

      element_coords[i].resize(elements.size());
      /* loop over local subgrids and set initial values */
      for (size_t j = 0; j < elements.size(); ++j) {
        IndexType globalLinearIndex = dfg.getGlobalLinearIndex(j);
        std::vector<real> globalCoords(dim);
        dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
        for (size_t l = 0; l < dim; ++l) {
          element_coords[i][j][l] = globalCoords[l];
        }

        elements[j] = NAN;  // TaskEnsemble::myfunction(globalCoords, 0.0);
      }
    }

    // corners is an array that contains exactly the one corner point at index i that is set by
    // deal.ii at Degree of freedom i.
    // So you know which direction this degree of freedom i is pointing
    corners.resize(dfgEnsemble_->getNumFullGrids());

    int partitions = p[0] * p[1];

    // set grid sizes
    int Nx = std::pow(2, l[0]) + 1;
    int Ny = std::pow(2, l[1]) + 1;
    int Nz = 1;
    if (dim == 3) {
      Nz = std::pow(2, l[2]) + 1;
      partitions *= p[2];
    }

    assert(partitions == 1);  // partitioning isnt correct when in parallel

#pragma region Index_mapping_init
    // indexmapping[i], initialized with -1 contains at index j (=index_mapping[i][j]) the index of
    // the deal.ii point that corresponds to the DisCoTec point at index j
    index_mapping.resize(dfgEnsemble_->getNumFullGrids());

    // loop over each degree of freedom and init index_mapping with -1 such that we later can
    // control which index is not set yet.
    for (unsigned int dof = 0; dof < dfgEnsemble_->getNumFullGrids(); ++dof) {
      index_mapping[dof].resize(element_coords[dof].size());
      for (unsigned int j = 0; j < element_coords[dof].size(); ++j) {
        index_mapping[dof][j] = -1;
      }
    }

    // this may cause memory problems when working with high discretizations,
    index_sub.resize(Nx * Ny * Nz);
    int x = 0, y = 0, z = 0, linearized_index = 0;

    // index mapping is done via hashing:
    // each point is mapped to its linearized index in the grid -> O(n)
    for (unsigned int el_index = 0; el_index < (element_coords[0].size()); el_index++) {
      x = element_coords[0][el_index][0] * (Nx - 1);
      y = element_coords[0][el_index][1] * (Ny - 1);

      linearized_index = y * (Nx) + x;
      if (dim == 3) {
        z = element_coords[0][el_index][2] * (Nz - 1);
        linearized_index = x + (Nx)*y + z * (Nx) * (Ny);
      }
      // remember the index of the point and its linearized index
      index_sub[linearized_index] = el_index;
    }
    // same for the deal.II points ->O(2n)
    for (unsigned int x2 = 0; x2 < (coords_dealii.size()); x2++) {
      // x,y and z are the indices of the point with respect to its corresponding axis, so it is not
      // the deal.ii index
      x = std::round(coords_dealii[x2][0] * (Nx - 1));
      y = std::round(coords_dealii[x2][1] * (Ny - 1));

      // this bool indicates wether this point is on the boundary in z-direction
      bool z_corner = true;
      if (dim == 3) {
        z = coords_dealii[x2][2] * (Nz - 1);
        if (!(z == (Nz - 1) || z == 0)) z_corner = false;
      }

      // check here whether cornerpoint, if so, save this point in corners
      if ((x == (Nx - 1) || x == 0) && (y == (Ny - 1) || y == 0) && z_corner) {
        corners[coords_dealii[x2][dim + 1]].resize(dim);
        for (unsigned int d = 0; d < dim; ++d)
          corners[coords_dealii[x2][dim + 1]][d] = coords_dealii[x2][d];
      }
      // claculate the linearized index
      linearized_index = (Nx)*y + x + z * (Nx) * (Ny);

      // If we already had an DisCoTec point with this linearized index we know these points need to
      // be mapped
      if (index_sub[linearized_index] != -1) {
        // coords_dealii[x2][dim + 1] decodes the degree of freedom.
        // then in this DOF save at the index of the DisCoTec point the index of the deal.ii point.
        // can be understood as elements[index_sub[linearized_index]]=coords_dealii[x2]
        index_mapping[coords_dealii[x2][dim + 1]][index_sub[linearized_index]] = x2;
      }
    }
#pragma endregion Index_mapping_init

#pragma region Index_mapping_boundary
    // check the index_mapping for Completeness, beacuse the boundary values, namely corners, edges
    // and (in 3d) boundary areas need special treatment
    for (unsigned int point_index = 0; point_index < element_coords[0].size(); ++point_index) {
      for (unsigned int dof = 0; dof < dfgEnsemble_->getNumFullGrids(); ++dof) {
        if (index_mapping[dof][point_index] == -1) {
          // x,y and z are the indices of the point with respect to its corresponding axis, so it is
          // not the discotec index
          //those values are used to detect the type of boundary point: corner, edge or boundary area
          x = std::round(element_coords[dof][point_index][0] * (Nx - 1));
          y = std::round(element_coords[dof][point_index][1] * (Ny - 1));
          bool z_corner = true;
          if (dim == 3) {
            z = element_coords[dof][point_index][2] * (Nz - 1);
            if (!(z == (Nz - 1) || z == 0)) z_corner = false;
          }

#pragma region Corner
          if ((x == (Nx - 1) || x == 0) && (y == (Ny - 1) || y == 0) && z_corner) {
            // this point is a corner
            // on one DisCoTec grid all corner points have the same value as it displays 1 Degree of
            // freedom and we have periodic boundary values

            //calculate the corner that is stored in this degree of freedom
            x = std::round(corners[dof][0] * (Nx - 1));
            y = std::round(corners[dof][1] * (Ny - 1));
            if (dim == 3) z = corners[dof][2] * (Nz - 1);

            //this is the linearized index of the corner point.
            linearized_index = (Nx)*y + x + z * (Nx) * (Ny);
            
            //this point needs the same value as the corner point and the corner points DisCoTec index ist stored in index_sub[linearized_index]
            index_mapping[dof][point_index] = index_mapping[dof][index_sub[linearized_index]];
          }
#pragma endregion Corner
#pragma region Edge
          else {
            //this is a edge
            //adding a part of the negated corner gives the point coordinates like in Deal.ii such that we know which dof of a point it is.
            double x_moved_point = element_coords[dof][point_index][0] + 1.0 / 100.0 * ((int)(!((bool)corners[dof][0])) - 0.5);
            double y_moved_point = element_coords[dof][point_index][1] + 1.0 / 100.0 * ((int)(!((bool)corners[dof][1])) - 0.5);
            double z_moved_point=0;
            //check whether the x_ coordinate is on the x-axis boundary
            if (element_coords[dof][point_index][0] == 0.0 || element_coords[dof][point_index][0] == 1.0) {
              if (x_moved_point< 0.0 || x_moved_point > 1.0) {
                //the moved point is outside the domain in x direction so copy the value from the negated x coordinate
                x_moved_point = ((int)!(element_coords[dof][point_index][0]));
              } else {
                //x_ is inside the domain so keep this value
                x_moved_point = element_coords[dof][point_index][0];
              }
            } else {
              x_moved_point = element_coords[dof][point_index][0];
            }

            //same for the y coordinate: If y is at y-axis boundary 
            if (element_coords[dof][point_index][1] == 0.0 || element_coords[dof][point_index][1] == 1.0) {
              if (y_moved_point < 0.0 || y_moved_point > 1.0) {
                //the moved point is outside the domain in y direction so copy the value from the negated y coordinate
                y_moved_point = ((int)!(element_coords[dof][point_index][1]));
              } else {
                y_moved_point = element_coords[dof][point_index][1];
              }
            } else {
              y_moved_point = element_coords[dof][point_index][1];
            }

            //the same for the z coordinate if we are three dimensional
            if (dim == 3) {

              //move the point in z-direction with a little offset according to the "direction" of the degree of freedom
              z_moved_point = element_coords[dof][point_index][2] + 1.0 / 100.0 * ((int)(!((bool)corners[dof][2])) - 0.5);
              if (element_coords[dof][point_index][2] == 0 || element_coords[dof][point_index][2] == 1) {
                if (z_moved_point < 0 || z_moved_point > 1) {
                  //the moved point is outside the domain in z direction so copy the value from the negated y coordinate
                  z_moved_point = (int)!(element_coords[dof][point_index][2]);
                } else {
                  z_moved_point = element_coords[dof][point_index][2];
                }
              } else {
                z_moved_point = element_coords[dof][point_index][2];
              }
            }
            linearized_index = (Nx)*y_moved_point * (Ny - 1) + x_moved_point * (Nx - 1) + z_moved_point * (Nz - 1) * (Nx) * (Ny);
            
            index_mapping[dof][point_index] = index_mapping[dof][index_sub[linearized_index]];
          }
#pragma endregion Edge
        }
      }
    }
#pragma endregion Index_mapping_boundary

    initialized_ = true;
  }

  /* this is were the application code kicks in and all the magic happens.
   * do whatever you have to do, but make sure that your application uses
   * only lcomm or a subset of it as communicator.
   * important: don't forget to set the isFinished flag at the end of the computation.
   */
  void run(CommunicatorType lcomm) {
    assert(initialized_);
    table.start("time->fullgrid");
    if (!make_exact_) {
      std::vector<std::array<Number, Problem::dim_ + 2>> old_result(dealii_size);
      if (stepsTotal_ == 1) {
        auto& dfg = dfgEnsemble_->getDFG(0);
        std::vector<CombiDataType>& elements = dfg.getElementVector();
        int counter = 0;
        for (size_t j = 0; j < elements.size(); ++j) {
          IndexType globalLinearIndex = dfg.getGlobalLinearIndex(j);
          std::vector<real> globalCoords(dim);
          dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
          std::cout << elements[j] << " ";
          counter++;
          if (counter == 5) {
            std::cout << std::endl;
            counter = 0;
          }
        }
      }
      if (stepsTotal_ > 0 && do_combine) {
        // iterates over all dfgs of the ensemble
        for (unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++) {
          auto& dfg = dfgEnsemble_->getDFG(i);
          std::vector<CombiDataType>& elements = dfg.getElementVector();

          // iterate here over all points of this grid
          for (unsigned int l = 0; l < index_mapping[i].size(); l++) {
            old_result[index_mapping[i][l]][Problem::dim_] = elements[l];
          }
        }
        problem->set_result(old_result);
      }
      // Stats::startEvent("Task "+std::to_string(this->getID()));
      problem->reinit_time_integration(stepsTotal_ * dt_, (stepsTotal_ + 1) * dt_);

      // process problem
      problem->solve();
      // Stats::stopEvent("Task "+std::to_string(this->getID()));

      // process problem
      problem->solve();
      // Stats::stopEvent("Task "+std::to_string(this->getID()));

      std::vector<std::array<Number, Problem::dim_ + 2>> result = problem->get_result();
      // std::cout << "result " << result.size() << " " << result[0].size() << std::endl;
      // for (auto& r : result) {
      //   for (auto& s : r) {
      //       std::cout << s << " ";
      //   }
      //   std::cout << "," << std::endl;
      // }
      // std::cout << std::endl;

      // iterate over all dfgs
      for (unsigned int i_dash = 0; i_dash < dfgEnsemble_->getNumFullGrids(); i_dash++) {
        auto& dfg = dfgEnsemble_->getDFG(i_dash);
        auto i = powerOfTwo[Problem::dim_] - i_dash - 1;
        std::vector<CombiDataType>& elements = dfg.getElementVector();

        // iterate over the elements of this dfg
        for (unsigned int l = 0; l < index_mapping[i].size(); l++) {
          elements[l] = result[index_mapping[i][l]][Problem::dim_];
        }
      }

    } else {
      std::vector<std::array<Number, Problem::dim_ + 2>> coords_dealii = problem->get_result();
      // iterate over all points here and compute the value here.
      for (unsigned int i = 0; i < coords_dealii.size(); ++i) {
        std::array<Number, Problem::dim_> coord;
        for (int dim = 0; dim < Problem::dim_; ++dim) {
          coord[dim] = coords_dealii[i][dim];
        }
        coords_dealii[i][Problem::dim_] = exactSolution(coord, (stepsTotal_ + 1) * dt_);
      }

      for (unsigned int i_dash = 0; i_dash < dfgEnsemble_->getNumFullGrids(); i_dash++) {
        auto& dfg = dfgEnsemble_->getDFG(i_dash);
        auto i = powerOfTwo[Problem::dim_] - i_dash - 1;
        std::vector<CombiDataType>& elements = dfg.getElementVector();

        // iterate over the elements of this dfg
        for (unsigned int l = 0; l < index_mapping[i].size(); l++) {
          elements[l] = coords_dealii[index_mapping[i][l]][Problem::dim_];
        }
      }
    }
    stepsTotal_++;
    table.stop("time->fullgrid");
    this->setFinished(true);
    // for(unsigned int i = 0; i< dfgEnsemble_->getNumFullGrids(); ++i) {
    //   std::cout << GridEnumeration::isHigherInDimension(static_cast<char>(0), i) << " " <<
    //   GridEnumeration::isHigherInDimension(static_cast<char>(1), i) << std::endl; std::cout <<
    //   dfgEnsemble_->getDFG(i) << std::endl;
    // }
  }

  double exactSolution(std::array<Number, Problem::dim_> coordinates, double time) {
    double result = 1.0;
    Tensor<1, 3> advection;
    advection[0] = 1;
    advection[1] = 0.15;
    advection[2] = -0.05;
    for (unsigned int d = 0; d < Problem::dim_; ++d)
      result *= std::pow(std::sin((coordinates[d] - time * advection[d]) * dealii::numbers::PI), 2);

    return result;
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) {
    return dfgEnsemble_->getDFG(n);
  }

  DFGEnsemble* getDFGEnsemble() const override { return dfgEnsemble_; }

  /* this function evaluates the combination solution on a given full grid.
   * here, a full grid representation of your task's solution has to be created
   * on the process of lcomm with the rank r.
   * typically this would require gathering your (in whatever way) distributed
   * solution on one process and then converting it to the full grid representation.
   * the DistributedFullGrid class offers a convenient function to do this.
   */
  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    assert(fg.getLevels() == getDistributedFullGrid(n).getLevels());

    getDistributedFullGrid(n).gatherFullGrid(fg, r);
  }

  // return the number of grids; here it is the number of grids in ensemble
  size_t getNumGrids() override { return powerOfTwo[dim_]; }

  size_t getDOFs() override {
    // std::vector<std::array<Number, Problem::dim_ + 2>> coordalii = problem->get_result();
    // dealii_size=coordalii.size();
    return dealii_size;
  }

  void setZero() {}

  ~TaskEnsemble() {
    if (dfgEnsemble_ != nullptr) delete dfgEnsemble_;
  }

 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskEnsemble() : make_exact_(false), initialized_(false), stepsTotal_(0), dfgEnsemble_(nullptr) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t dealii_size = 0;
  std::vector<std::vector<int>> index_mapping;
  std::vector<int> index_sub;
  std::vector<std::vector<Number>> corners;
  std::vector<std::vector<Number>> index_mapping2;
  bool make_exact_ = false;
  IndexVector p_;
  std::shared_ptr<Problem> problem;
  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  // DistributedFullGrid<CombiDataType>* dfg_;
  DFGEnsemble* dfgEnsemble_;
  std::string _filename;
  /**
   * The serialize function has to be extended by the new member variables.
   * However this concerns only member variables that need to be exchanged
   * between manager and workers. We do not need to add "local" member variables
   * that are only needed on either manager or worker processes.
   * For serialization of the parent class members, the class must be
   * registered with the BOOST_CLASS_EXPORT macro.
   */
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    // handles serialization of base class
    ar& boost::serialization::base_object<Task>(*this);

    // add our new variables
    ar& dt_;

    ar& p_;
    ar& _filename;
    ar& make_exact_;
  }
};

}  // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
