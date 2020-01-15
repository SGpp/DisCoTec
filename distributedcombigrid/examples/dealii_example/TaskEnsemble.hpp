/*
 * TaskExample.hpp
 *
 *  Created on: Sep 25, 2015
 *      Author: heenemo
 */

#ifndef TASKENSEMPLE_HPP_
#define TASKENSEMPLE_HPP_

#include "sgpp/distributedcombigrid/fullgrid/DistributedFullGridEnsemble.hpp"
#include "sgpp/distributedcombigrid/task/Task.hpp"

#include <deal.II/base/timer.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <hyper.deal.combi/include/functionalities/dynamic_convergence_table.h>
#include <hyper.deal.combi/include/functionalities/vector_dummy.h>
//#include <hyper.deal.combi/applications/advection_reference_dealii/include/application.h>


namespace combigrid {

class TaskEnsemble : public Task {
 public:
  /* if the constructor of the base task class is not sufficient we can provide an
   * own implementation. here, we add dt, nsteps, and p as a new parameters.
   */
  TaskEnsemble(DimType dim, LevelVector& l, std::vector<bool>& boundary, real coeff,
              LoadModel* loadModel, std::string filename,  real dt,  IndexVector p = IndexVector(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit),
        dt_(dt),
        p_(p),
        initialized_(false),
        stepsTotal_(0),
        dfgEnsemble_(nullptr),
        _filename(filename) {}

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

    if (lrank == 0) {
      std::cout << "init task " << this->getID() << " with l = " << this->getLevelVector()
                << " and p = " << p << std::endl;
    }

    // create local subgrid on each process
    dfgEnsemble_ = new DFGEnsemble(dim, l, lcomm, this->getBoundary(), p);


    //for visualization of grids in Paraview
    Point<Problem::dim_> offset;
    for(unsigned int i=0;i<dim;i++)
      offset[i]=l[i]*1.1;

    
    this->problem = std::make_shared<Problem>(lcomm, table,offset);
    this->problem->reinit(_filename);


      
    
    if(do_combine){
    std::vector<std::array<Number, Problem::dim_ + 2>> coords_dealii = problem->get_result();
    size_result=coords_dealii.size();

    std::vector<std::vector<std::array<Number, Problem::dim_ >>> element_coords(dfgEnsemble_->getNumFullGrids()); 

    /* loop over local subgrids and set initial values */
    for (unsigned int i = 0; i< dfgEnsemble_->getNumFullGrids(); ++i){
      auto& dfg = dfgEnsemble_->getDFG(i);
      std::vector<CombiDataType>& elements = dfg.getElementVector();

      element_coords[i].resize(elements.size());

      for (size_t j = 0; j < elements.size(); ++j) {
        IndexType globalLinearIndex = dfg.getGlobalLinearIndex(j);
        std::vector<real> globalCoords(dim);
        dfg.getCoordsGlobal(globalLinearIndex, globalCoords);
        for(size_t l=0; l<dim;++l)
          element_coords[i][j][l]=globalCoords[l];
      
        elements[j] =0;// TaskEnsemble::myfunction(globalCoords, 0.0);
      }
    }

    int vec=0;
    double eps=5*10e-8;
    
    index_mapping.resize(dfgEnsemble_->getNumFullGrids());
    
    for (unsigned int i = 0; i< dfgEnsemble_->getNumFullGrids(); ++i){
      auto& dfg = dfgEnsemble_->getDFG(i);
      std::vector<CombiDataType>& elements = dfg.getElementVector();

      index_mapping[i].resize(elements.size());

      for(unsigned int el_index=0;el_index<elements.size();el_index++){
      //Create the point for the point of combi
        vec=0;
        Point<Problem::dim_> combi;
        for(auto y:element_coords[i][el_index])
          combi[vec++]=y;
      //loop over dealii points
        for(unsigned int x2=(coords_dealii.size())-1;x2<numbers::invalid_unsigned_int;x2--){
        //for(unsigned int x2=0;x2<(coords_dealii.size());x2++){
          //create Point
          Point<Problem::dim_> dealii;
          vec=0;
          for(unsigned int vec=0;vec<dim;vec++)
            dealii[vec]=coords_dealii[x2][vec];

          //check if the points are equal
          if(combi.distance(dealii)<=eps && coords_dealii[x2][dim+1]==i )
          {
            index_mapping[i][el_index]=x2;
            //std::cout << x2 << " ";
            //std::cout <<"Combi Point: "<<combi << " Dealii: "<<dealii<<std::endl;
            break;//break is bad
          }
        }
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

    int lrank;
    MPI_Comm_rank(lcomm, &lrank);
    // TODO if your Example uses another data structure, you need to copy
    // the data from elements to that data structure

    /* pseudo timestepping to demonstrate the behaviour of your typical
     * time-dependent simulation problem. */
    // TODO replace by your time time stepping algorithm

    
   // 
    if(stepsTotal_>0 && do_combine)
    {
      std::vector<std::array<Number, Problem::dim_ + 2>> old_result(size_result);

      //iterates over all dfgs of the ensemble
      for(unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++){
        auto& dfg = dfgEnsemble_->getDFG(i);
        std::vector<CombiDataType>& elements = dfg.getElementVector();

        //iterate here over all points of this grid
        for(unsigned int l = 0; l < elements.size(); l++) {       

          old_result[index_mapping[i][l]][Problem::dim_]=elements[l];
          //old_result[index_mapping[i][l]][Problem::dim_+1]=i;
        }
      }
      problem->set_result(old_result);
    }
    problem->reinit_time_integration(stepsTotal_*dt_, (stepsTotal_ + 1)*dt_);

    //process problem
    problem->solve();
    
    if(do_combine){
      std::vector<std::array<Number, Problem::dim_ + 2>> result = problem->get_result();

      //iterate over all dfgs
      for(unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++){         
        auto& dfg = dfgEnsemble_->getDFG(i);
        std::vector<CombiDataType>& elements = dfg.getElementVector();
  
        //iterate over the elements of this dfg
        for(unsigned int l=0; l<elements.size(); l++){
          //TODO: missing values at boundary
          elements[l]=result[index_mapping[i][l]][Problem::dim_];

        }
      }          
    }
    stepsTotal_ ++;

    this->setFinished(true);
  }

  DistributedFullGrid<CombiDataType>& getDistributedFullGrid(int n = 0) { return dfgEnsemble_->getDFG(n); }

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

  static real myfunction(std::vector<real>& coords, real t) {
    real u = std::cos(M_PI * t);

    for (size_t d = 0; d < coords.size(); ++d) u *= std::cos(2.0 * M_PI * coords[d]);

    return u;

    /*
    double res = 1.0;
    for (size_t i = 0; i < coords.size(); ++i) {
      res *= -4.0 * coords[i] * (coords[i] - 1);
    }


    return res;
    */
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
  TaskEnsemble() : initialized_(false), stepsTotal_(0), dfgEnsemble_(nullptr) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t size_result;
  std::vector<std::vector<Number>> index_mapping;
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
  }
};

}  // namespace combigrid

#endif /* TASKEXAMPLE_HPP_ */
