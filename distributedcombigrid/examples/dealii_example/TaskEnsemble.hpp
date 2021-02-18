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
#include <sstream>
#include <cassert>
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
              LoadModel* loadModel, std::string filename,  real dt, bool makeExactSolution, IndexVector p = IndexVector(0),
              FaultCriterion* faultCrit = (new StaticFaults({0, IndexVector(0), IndexVector(0)})))
      : Task(dim, l, boundary, coeff, loadModel, faultCrit){
        dt_=(dt),
        p_=(p),
        make_exact_=(makeExactSolution),
        initialized_=(false),
        stepsTotal_=(0),
        dfgEnsemble_=(nullptr),
        _filename=(filename); }

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


    //for visualization of grids in Paraview
    Point<Problem::dim_> offset;
    for(unsigned int i=0;i<dim;i++)
      offset[i]=l[i]*1.1;

    
    if(!make_exact_){
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
        for(size_t l=0; l<dim;++l){
          element_coords[i][j][l]=globalCoords[l];
        }
        
        elements[j] =0;// TaskEnsemble::myfunction(globalCoords, 0.0);
      }
    }
    std::stringstream ss;
    int Nx=std::pow(2,l[0])+1;
    int Ny=std::pow(2,l[1])+1;
    int Nz=1;
    int partitions=p[0]*p[1];
    if(dim==3){
      Nz=std::pow(2,l[2])+1;
      partitions*=p[2];
    }
      
    assert(partitions==1);//partitioning isnt correct when in parallel
    index_mapping.resize(dfgEnsemble_->getNumFullGrids());
    for (unsigned int i = 0; i< dfgEnsemble_->getNumFullGrids(); ++i)
      index_mapping[i].resize(element_coords[0].size());
    //this may cause memory problems when working with high discretizations, 
    std::vector<int> index_sub(Nx*Ny*Nz,-1);
    
    int z=0;
    int x=0,y=0,linearized_index=0;
    //index mapping is done via hashing:
    //each point is mapped to its linearized index in the grid -> O(n)
    
    for(unsigned int el_index=0;el_index<(element_coords[0].size());el_index++){
      
      x=element_coords[0][el_index][0]*(Nx-1);
      y=element_coords[0][el_index][1]*(Ny-1);      
      if(dim==2){               
        linearized_index=y*(Nx)+x;
      }
      else if(dim==3){
        z=element_coords[0][el_index][2]*(Nz-1); 
        linearized_index=x+(Nx)*y+z*(Nx)*(Ny);
      }
      //ss << "["<<linearized_index<<"]";
      index_sub[linearized_index]=el_index;
    }
    
    //same for the deal.II points ->O(2n)
    for(unsigned int x2=0; x2<(coords_dealii.size());x2++){
      x=std::round(coords_dealii[x2][0]*(Nx-1));
      y= std::round(coords_dealii[x2][1]*(Ny-1));    
      if(dim==3)
        z=coords_dealii[x2][2]*(Nz-1);  
      
      
      linearized_index=(Nx)*y+x+z*(Nx)*(Ny);
      
      if(index_sub[linearized_index]!=-1)        
        index_mapping[coords_dealii[x2][dim+1]][index_sub[linearized_index]]=x2;
      
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

  table.start("time->fullgrid");
   if(!make_exact_){
     std::vector<std::array<Number, Problem::dim_ + 2>> old_result(size_result);
    
    if(stepsTotal_>0)
    {
      

      //iterates over all dfgs of the ensemble
      for(unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++){
        auto& dfg = dfgEnsemble_->getDFG(i);
        std::vector<CombiDataType>& elements = dfg.getElementVector();

        //iterate here over all points of this grid
        for(unsigned int l = 0; l < index_mapping[i].size(); l++) {       

          old_result[index_mapping[i][l]][Problem::dim_]=elements[l];
          //old_result[index_mapping[i][l]][Problem::dim_+1]=i;
        }
      }
      problem->set_result(old_result);
    }
    Stats::startEvent("Task "+std::to_string(this->getID()));
    problem->reinit_time_integration(stepsTotal_*dt_, (stepsTotal_ + 1)*dt_);

    //process problem
    problem->solve();
    Stats::stopEvent("Task "+std::to_string(this->getID()));
    
    std::vector<std::array<Number, Problem::dim_ + 2>> result = problem->get_result();

      //iterate over all dfgs
      for(unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++){         
        auto& dfg = dfgEnsemble_->getDFG(i);
        std::vector<CombiDataType>& elements = dfg.getElementVector();
  
        //iterate over the elements of this dfg
        for(unsigned int l=0; l<index_mapping[i].size(); l++){
          //TODO: missing values at boundary
          elements[l]=result[index_mapping[i][l]][Problem::dim_];

        }
      }          
    
  }
  else{
    for(unsigned int i = 0; i < dfgEnsemble_->getNumFullGrids(); i++){
      auto& dfg_ = dfgEnsemble_->getDFG(i);
      std::vector<CombiDataType>& elements = dfg_.getElementVector();
      std::vector<std::array<Number, Problem::dim_ >> element_coords(elements.size());
      for (size_t j = 0; j < elements.size(); ++j) {
        IndexType globalLinearIndex = dfg_.getGlobalLinearIndex(j);
        std::vector<real> globalCoords(dim);
        dfg_.getCoordsGlobal(globalLinearIndex, globalCoords);
        for(size_t l=0; l<dim;++l)
          element_coords[j][l]=globalCoords[l];
        //the coordinates are now stored in elemen_coords[i]
        //compute the exact solution on this point and store it in the grid;
        elements[j] = exactSolution(element_coords[j],(stepsTotal_+1)*dt_);
      }
      }
  }
    stepsTotal_ ++;
  table.stop("time->fullgrid");
    this->setFinished(true);
  }

   double exactSolution(std::array<Number, Problem::dim_ > coordinates, double time){
    double result=1.0;
    Tensor<1, 3> advection;
    advection[0]=1;
    advection[1]=0.15;
    advection[2]=-0.05;
    for(unsigned int d = 0; d < Problem::dim_; ++d)
      result *= std::pow(std::sin((coordinates[d]-time*advection[d])*dealii::numbers::PI),2);
    
    return result;    
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
  TaskEnsemble() : make_exact_(false),initialized_(false), stepsTotal_(0), dfgEnsemble_(nullptr) {}

 private:
  friend class boost::serialization::access;

  // new variables that are set by manager. need to be added to serialize
  real dt_;
  size_t size_result;
  std::vector<std::vector<Number>> index_mapping;
  std::vector<std::vector<Number>> index_mapping2;
  bool make_exact_=false;
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