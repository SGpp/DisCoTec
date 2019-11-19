/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKEXAMPLE_HPP_
#define TASKEXAMPLE_HPP_

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGrid.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

#include <deal.II/base/timer.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <hyper.deal.combi/include/functionalities/dynamic_convergence_table.h>
#include <hyper.deal.combi/include/functionalities/vector_dummy.h>
#include <hyper.deal.combi/applications/advection_reference_dealii/include/application.h>

namespace combigrid {

class TaskExample : public Task {
 public:
  typedef Application<dim, degree, degree + 1, Number, VectorizedArrayType, VectorType> Problem;

  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskExample(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
              LoadModel* loadModel,std::string filename,  real dt, size_t nsteps, IndexVector p = IndexVector(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        nsteps_(nsteps),
        p_(p),
        initialized_(false),
        stepsTotal_(0),dfg_(NULL),
        _filename(filename){}
        

  void init(CommunicatorType lcomm,
            std::vector<IndexVector> decomposition = std::vector<IndexVector>()) {
    assert(!initialized_);
    assert(dfg_ == NULL);

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

    if (lrank == 0) {
      std::cout << "init task " << this->getID() << " with l = " << this->getLevelVector()
                << " and p = " << p << std::endl;
    }

    // create local subgrid on each process
    dfg_ = new DistributedFullGrid<CombiDataType>(dim, l, lcomm, this->getBoundary(), p);

    /* loop over local subgrid and set initial values */
    this->problem = std::make_shared<Problem>(lcomm,  table);
    this->problem->reinit(_filename);

    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    std::vector<std::array<Number, Problem::dim_ >> element_coords(elements.size());
    std::vector<std::array<Number, Problem::dim_ + 1>> coords_dealii = problem->get_result();
    size_result=coords_dealii.size();
    for (size_t i = 0; i < elements.size(); ++i) {
      IndexType globalLinearIndex = dfg_->getGlobalLinearIndex(i);
      std::vector<real> globalCoords(dim);
      dfg_->getCoordsGlobal(globalLinearIndex, globalCoords);
      for(size_t l=0; l<dim;++l)
        element_coords[i][l]=globalCoords[l];
      
      elements[i] = 0;

    }
   
    int i=0;
    double eps=5*10e-16;
    index_mapping.resize(element_coords.size());

    for(unsigned int el_index=0;el_index<(element_coords.size());el_index++){
      //Create the point for the point of combi
      i=0;
      Point<Problem::dim_> combi;
      for(auto y:element_coords[el_index])
        combi[i++]=y;
      //loop over dealii points
      for(unsigned int x2=(coords_dealii.size())-1;x2>=numbers::invalid_unsigned_int;x2--){
        //create Point
        Point<Problem::dim_> dealii;
        i=0;
        for(auto y:coords_dealii[x2])
          dealii[i++]=y;

        //check if the points are equal
        if(combi.distance(dealii)<=eps)
        {
          index_mapping[el_index]=x2;
          //std::cout << x2 << " ";
          //std::cout <<"Combi Point: "<<combi << " Dealii: "<<dealii<<std::endl;
          break;//break is bad
        }
      }
    }
      
    initialized_ = true;
  }

  /* this is were the application code kicks in and all the magic happens.
   * do whatever you have to do, but make sure that your application uses
   * only lcomm or a subset of it as communicator.
   * important: don't forget to set the isFinished flag at the end of the computation.
   */

  
  void run(CommunicatorType lcomm) {
    assert(initialized_);
    std::cout << "Run of Task"<< this->getID()<<std::endl;
    std::vector<CombiDataType>& elements = dfg_->getElementVector();
    
    
    std::vector<std::array<Number, Problem::dim_ + 1>> old_result(size_result);
   // 
   if(stepsTotal_>0)
    for(unsigned int i = 0; i < index_mapping.size(); i++)
      old_result[index_mapping[i]][Problem::dim_]=elements[i].real();
    
    problem->set_result(old_result);
    problem->reinit_time_integration(stepsTotal_*dt_, (stepsTotal_ + nsteps_)*dt_);

    //process problem
    problem->solve();

    std::vector<std::array<Number, Problem::dim_ + 1>> result = problem->get_result();
    
    for(unsigned int i = 0; i < index_mapping.size(); i++)
      elements[i]=result[index_mapping[i]][Problem::dim_];
    
    stepsTotal_ += nsteps_;

    
    this->setFinished(true);
  }

  /* this function evaluates the combination solution on a given full grid.
   * here, a full grid representation of your task's solution has to be created
   * on the process of lcomm with the rank r.
   * typically this would require gathering your (in whatever way) distributed
   * solution on one process and then converting it to the full grid representation.
   * the DistributedFullGrid class offers a convenient function to do this.
   */
  void getFullGrid(FullGrid<CombiDataType>& fg, RankType r, CommunicatorType lcomm, int n = 0) {
    assert(fg.getLevels() == dfg_->getLevels());

    dfg_->gatherFullGrid(fg, r);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return *dfg_; }

  void setZero() {}

  ~TaskExample() {
    if (dfg_ != NULL) delete dfg_;
  }


 protected:
  /* if there are local variables that have to be initialized at construction
   * you have to do it here. the worker processes will create the task using
   * this constructor before overwriting the variables that are set by the
   * manager. here we need to set the initialized variable to make sure it is
   * set to false. */
  TaskExample() : initialized_(false), stepsTotal_(0), dfg_(NULL) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;       // TODO
  size_t nsteps_; // TODO
  size_t size_result;
  std::vector<Number> index_mapping;
  IndexVector p_;
  std::shared_ptr<Problem> problem;

  // pure local variables that exist only on the worker processes
  bool initialized_;
  size_t stepsTotal_;
  DistributedFullGrid<CombiDataType>* dfg_;
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
    ar& nsteps_;
    ar& p_;
    ar& _filename;
  }
};

}  // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
