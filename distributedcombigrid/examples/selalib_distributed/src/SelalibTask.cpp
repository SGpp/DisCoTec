/*
 * SelalibTask.cpp
 *
 *  Created on: Nov 19, 2020
 *      Author: obersteiner
 */

#include "SelalibTask.hpp"
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <unistd.h>
#include <fstream>
//#include "CombiGeneConverter.hpp"
#include "sgpp/distributedcombigrid/mpi/MPISystem.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupSignals.hpp"
#include <math.h>

//#include "sgpp/distributedcombigrid/utils/StatsContainer.hpp"

namespace combigrid
{

SelalibTask::SelalibTask( DimType dim, LevelVector& l,
                    std::vector<bool>& boundary, real coeff, LoadModel* loadModel,
                    std::string& path, real dt, real combitime, size_t nsteps,
                    IndexVector p , FaultCriterion *faultCrit,
                    IndexType numSpecies, size_t checkpointFrequency, size_t offsetForDiagnostics)
    : Task( dim, l, boundary, coeff, loadModel, faultCrit),
      path_( path ),
      dt_( dt ),
      combitime_(combitime),
      nsteps_( nsteps ),
      stepsTotal_(0),
      combiStep_(0),
      p_(p),
      checkpoint_(), initialized_(false),
      checkpointInitialized_(false),
      nspecies_(numSpecies),
      currentTime_(0.0),
      checkpointFrequency_(checkpointFrequency),
      offsetForDiagnostics_(offsetForDiagnostics)
{

// currently all boundaries are set
assert( boundary[0] == true );//x
assert( boundary[1] == true );//y
assert( boundary[2] == true );//z
assert( boundary[3] == true );//v
assert( boundary[4] == true );//w
assert( boundary[5] == true );//z
}

SelalibTask::SelalibTask() :
    checkpoint_(),
    dfgVector_(0),
    initialized_(false),
    checkpointInitialized_(false),
    currentTime_(0.0)
{
  ;
}

SelalibTask::~SelalibTask()
{
  for(unsigned int g = 0; g < dfgVector_.size(); g++){
    delete dfgVector_[g];
  }
  dfgVector_.clear();
}

/**
 * This routine is used to initialize everything needed for running GENE.
 * This routine does NOT actually run GENE. It only changes directories
 */
void
SelalibTask::run( CommunicatorType lcomm )
{
  using namespace std::chrono;

  // change dir to wdir
  if( chdir( path_.c_str() ) ){


    printf( "could not change to directory %s \n", path_.c_str() );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }


  // it is more save to wait here until all procs are in the right directory
  MPI_Barrier( lcomm );

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "run task " << this->getID() << std::endl;
  }

}
/**
 * This routine is used to change the directory to the directory of the current task
 */
void SelalibTask::changeDir(CommunicatorType lcomm){
  // change dir to wdir
  if( chdir( path_.c_str() ) ){
    printf( "could not change to directory %s \n", path_.c_str() );
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  // it is more save to wait here until all procs are in the right directory
  MPI_Barrier( lcomm );

  MASTER_EXCLUSIVE_SECTION{
    std::cout << "changed to task " << this->getID() << std::endl;
  }
}
/**
 * This method is used to kill a process according to the faultCriterion.
 * Failing processes call the Sim_FT_kill_me() function
 */
void SelalibTask::decideToKill(){ //toDo check if combiStep should be included in task and sent to process groups in case of reassignment
  using namespace std::chrono;

  int globalRank;
  // MPI_Comm_rank(lcomm, &lrank);
  MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

  //check if killing necessary
  if (combiStep_ != 0 && faultCriterion_->failNow(combiStep_, -1.0, globalRank)){
        std::cout<<"Rank "<< globalRank <<" failed at iteration "<<combiStep_<<std::endl;
        StatusType status=PROCESS_GROUP_FAIL;
        MASTER_EXCLUSIVE_SECTION{
          simft::Sim_FT_MPI_Send( &status, 1, MPI_INT,  theMPISystem()->getManagerRank(), statusTag,
                            theMPISystem()->getGlobalCommFT() );
        }
        theMPISystem()->sendFailedSignal();
        simft::Sim_FT_kill_me();
  }
  combiStep_++;
}

/**
 * This routine initializes SelalibTask; currently it only sets a bool value.
 */
void SelalibTask::init(CommunicatorType lcomm, std::vector<IndexVector> decomposition){
//  if( dfg_ == NULL ){
//      dfg_ = new DistributedFullGrid<CombiDataType>( dim_, l_, lcomm,
//          this->getBoundary(), p_, false);
//  }
  initialized_ = true;
}

/**
 * This routine writes the grid values from gene (data) to a checkpoint file
 * @param data GENE grid
 * @param size size of the data grid (total)
 * @param sizes size of the data grid in each dimension
 * @param bounds bounds of the grid in distributed implementation
 */
void
SelalibTask::writeLocalCheckpoint( GeneComplex* data, size_t size,
                                std::vector<size_t>& sizes,
                                std::vector<size_t>& bounds )
{
  std::cout << "Number of species in checkpoint: " << sizes[0] << "\n";
  // todo: doing it like this will require two times copying
  for(unsigned int i= 0; i < sizes.size(); i++){
    int index_l = sizes.size()- 1 - i ; // sizes is reversed order of l; i.e. l is x y z v w spec and sizes spec, w, v, z, y, x
    assert(sizes[i] == pow(2,l_[index_l]));
      
  }
  checkpoint_.writeCheckpoint( data, size, sizes, bounds );
  checkpointInitialized_= true;

}
/**
 * Initializes local checkpoint.
 * This method is needed in case tasks are redistributed or recomputed
 * and the checkpoint is not initialized during read memory
 * @param size size of data (total)
 * @param sizes size of data in each dimension
 * @param bounds bounds of data in distributed array
 */
void SelalibTask::InitLocalCheckpoint(size_t size,
    std::vector<size_t>& sizes,
    std::vector<size_t>& bounds ){
  for(unsigned int i= 0; i < sizes.size(); i++){
    int index_l = sizes.size()- 1 - i ; // sizes is reversed order of l; i.e. l is x y z v w spec and sizes spec, w, v, z, y, x
    assert(sizes[i] == pow(2,l_[index_l])); //check if parameters match
  }
  checkpoint_.initCheckpoint(size,sizes,bounds);
  checkpointInitialized_= true;
}

/**
 * This routine returns the complete fullgrid on rank lroot
 * Therefore the fullgrid is collected from the other ranks of the local communicator
 */
void SelalibTask::getFullGrid( FullGrid<CombiDataType>& fg, RankType lroot,
                            CommunicatorType lcomm, int species)
{
  dfgVector_[species]->gatherFullGrid(fg, lroot);
}

/**
 * This routine returns the local part of the fullgrid
 */
DistributedFullGrid<CombiDataType>& SelalibTask::getDistributedFullGrid(int specie){
  return *dfgVector_[specie];
}

void
SelalibTask::saveCheckpoint( FullGrid<CombiDataType>& fg, const char* filename )
{
// first make sure we have a valid fullgrid
assert( fg.isGridCreated() );
assert( fg.getDimension() == 6 );
const std::vector<bool>& boundary = fg.returnBoundaryFlags();
assert( boundary[0] == true ); //x
assert( boundary[1] == false ); //y
assert( boundary[2] == true ); // z
assert( boundary[3] == true ); // v
assert( boundary[4] == true ); // w
assert( boundary[5] == false ); // nspec
assert( fg.getLevels()[1] == 1 ); // y
assert( fg.getLevels()[5] == 1 ); // nspec

size_t nspec = 1;
size_t nw = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[4] ) );
size_t nv = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[3] ) );
size_t nz = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[2] ) );
size_t ny = 1;
size_t nx = static_cast<size_t>( std::pow( 2.0, fg.getLevels()[0] ) );

// create MultiArray to store Gene grid temporarily
// note reverse notation of shape vector
IndexVector selalibShape = { (IndexType) nspec, (IndexType) nw, (IndexType) nv, (IndexType) nz, (IndexType) ny, (IndexType) nx };
MultiArray<CombiDataType,6> selalibGrid = createMultiArray<CombiDataType,6>( selalibShape );

// create MultiArrayRef to fg
MultiArrayRef6 fgData = createMultiArrayRef<CombiDataType,6>( fg );

// copy data to gene grid
// note that the last grid points in x,z,v,w dimension are ignored
for( size_t n=0; n < selalibShape[0]; ++n ){ //n_spec
  for( size_t m=0; m < selalibShape[1]; ++m ){ //w
    for( size_t l=0; l < selalibShape[2]; ++l ){ //v
      for( size_t k=0; k < selalibShape[3]; ++k ){ //z
        for( size_t j=0; j < selalibShape[4]; ++j ){ //y
          for( size_t i=0; i < selalibShape[5]; ++i ){ //x
            selalibGrid[n][m][l][k][j][i] = fgData[n][m][l][k][j][i];
          }
        }
      }
    }
  }
}

// save to file
std::ofstream cpFile( filename, std::ofstream::binary );

std::cout << "saving full grid as GENE checkpoint" << std::endl;

// write prec flag
const char* prec = "DOUBLE";
cpFile.write( prec, 6 );
std::cout << "prec:(6)" << prec << std::endl;

// write time and dt
// todo: use real values here
double mytime( 0 ), dt( 0 );
cpFile.write( (char*) &mytime, sizeof(double) );
cpFile.write( (char*) &dt, sizeof(double) );
std::cout << "mytime(" << sizeof(double) << "): " << mytime << std::endl;
std::cout << "dt(" << sizeof(double) << "): " << dt << std::endl;

// write resolution
int res[6];
res[0] = static_cast<int>(nx);
res[1] = static_cast<int>(ny);
res[2] = static_cast<int>(nz);
res[3] = static_cast<int>(nv);
res[4] = static_cast<int>(nw);
res[5] = static_cast<int>(nspec);
cpFile.write( (char*) res, sizeof(int) * 6 );
std::cout << "resolution:" << "\n" << "\t nx(" << sizeof(int) << ") " << nx
    << "\n" << "\t ny(" << sizeof(int) << ") " << ny << "\n" << "\t nz("
    << sizeof(int) << ") " << nz << "\n" << "\t nv(" << sizeof(int) << ") "
    << nv << "\n" << "\t nw(" << sizeof(int) << ") " << nw << "\n"
    << "\t n_spec(" << sizeof(int) << ") " << nspec << std::endl;

// write data
size_t dsize( selalibGrid.num_elements() );
cpFile.write( (char*) selalibGrid.origin(), sizeof(GeneComplex) * dsize );

cpFile.close();
}

/**
 * This routine sets the dfg to zero
 */
void SelalibTask::setZero(){
  if(dfgVector_.size() != 0){
    for(int i=0; i< dfgVector_.size(); i++){
      std::vector<CombiDataType>& data = dfgVector_[i]->getElementVector();

      for( size_t i=0; i<data.size(); ++i ){
        data[i] = CombiDataType(0.0);
      }
    }
  }
}

/**
 * This routine initializes the dfg only if it doesn't exist so far
 * The content of the dfg is set to zero in case it is initialized.
 */
void SelalibTask::initDFG( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){
  // todo: keep in mind
  // in this version the dfg is only created once. this only works if always exactly
  // the same set of processes is used by gene
  // this will probably not work, when the task is moved to another group.
  if( dfgVector_.size() == 0 ){
    dfgVector_.resize(nspecies_,NULL);
    for(int i=0; i<nspecies_; i++){
      dfgVector_[i] = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
          this->getBoundary(), p_, false, decomposition );
    }
    setZero();

  }
  //std::cout << "initDFG \n";

}

/**
 * This routine creates a new dfg (and initializes it). This is needed when a task is redistributed or recomputed
 * The content of the dfg is set to zero in case it is initialized.
 */
void SelalibTask::initDFG2( CommunicatorType comm,
                        std::vector<IndexVector>& decomposition ){

  // this is the clean version. however requires creation of dfg before each
  // combination step

  if(dfgVector_.size() != nspecies_){
    dfgVector_.resize(nspecies_,NULL);
  }
  for(int i=0; i<nspecies_; i++){
    if(dfgVector_[i] != NULL ){
      delete dfgVector_[i];
    }
    dfgVector_[i] = new DistributedFullGrid<CombiDataType>( dim_, l_, comm,
        this->getBoundary(), p_, false, decomposition );
  }
  setZero();
}


/** copy and convert data from local checkpoint to dfg
 * This is used for importing data from GENE to our framework
 * */
void SelalibTask::setDFG(){
  // step one: copy data of localcheckpoint into dfg

  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = {(IndexType) b[0],(IndexType) b[2],(IndexType) b[4],(IndexType) b[6],(IndexType) b[8],(IndexType) b[10] };
  IndexVector upperBounds = {(IndexType) b[1], (IndexType)b[3],(IndexType) b[5],(IndexType) b[7],(IndexType) b[9],(IndexType) b[11] };
  IndexVector lcpSizes = upperBounds - lowerBounds;

  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<CombiDataType,6>( checkpoint_.getData(), lcpSizes );
  for(int d=0; d < dim_; d++){
    for(int species=0; species<nspecies_; species++){

      MultiArrayRef6 dfgData =
          createMultiArrayRef<CombiDataType,6>( *dfgVector_[species] );

      // some checks
      const IndexVector p( dfgVector_[species]->getParallelization().rbegin(),
                            dfgVector_[species]->getParallelization().rend() );
      IndexVector tmp( dfgVector_[species]->getDimension() );
      dfgVector_[species]->getPartitionCoords( tmp );

      // get partition coords of the current process in the ordering as
      // the shape of the multi arrays, i.e. spec, w, v, z, y, x
      IndexVector coords( tmp.rbegin(), tmp.rend() );

      const size_t* dfgShape = dfgData.shape();
      const size_t* lcpShape = lcpData.shape();

      for( DimType d=0; d<dfgVector_[species]->getDimension(); ++d ){
        // for the last rank in the dimension w(d=1), v(2), z(3), y(4) (only in non-linear), x(5) the number of elements in
        // dfg and lcp differ
        //last process in line has boundary points not included in gene
        if( coords[d] == p[d] - 1){
          assert( dfgShape[d] == lcpShape[d] + 1 );
        } else{

              assert( dfgShape[d] == lcpShape[d] );
        }
      }
      // copy data from local checkpoint to dfg
      // note that on the last process in some dimensions dfg is larger than the
      // local checkpoint
      for( size_t n=0; n < lcpShape[0]; ++n){ //w
        for( size_t m=0; m < lcpShape[1]; ++m ){ //v
          for( size_t l=0; l < lcpShape[2]; ++l ){ //u
            for( size_t k=0; k < lcpShape[3]; ++k ){ //z
              for( size_t j=0; j < lcpShape[4]; ++j ){ //y
                for( size_t i=0; i < lcpShape[5]; ++i ){ //x
                  dfgData[0][m][l][k][j][i]= lcpData[species][m][l][k][j][i];
                }
              }
            }
          }
        }
      }



    /*
    // check if last grid point of x is zero
    if( coords[1] == p[1] - 1 ){
      for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
        for( size_t m=0; m < dfgShape[1]; ++m ) //w
          for( size_t l=0; l < dfgShape[2]; ++l ) //v
            for( size_t k=0; k < dfgShape[3]; ++k ) //z
              for( size_t j=0; j < dfgShape[4]; ++j ){ //y
                  assert( dfgData[n][ m ][l][k][j][ dfgShape[5]-1 ]
                          == complex(0.0) );
                }
    }

    // check if last grid point of w is zero
    if( coords[1] == p[1] - 1 ){
      for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
          for( size_t l=0; l < dfgShape[2]; ++l ) //v
            for( size_t k=0; k < dfgShape[3]; ++k ) //z
              for( size_t j=0; j < dfgShape[4]; ++j ) //y
                for( size_t i=0; i < dfgShape[5]; ++i ){ //x
                  assert( dfgData[n][ dfgShape[1]-1 ][l][k][j][i]
                          == complex(0.0) );
                }
    }

    // check if last grid point of v is zero
    if( coords[2] == p[2] - 1 ){
      for( size_t n=0; n < dfgShape[0]; ++n ) //n_spec
        for( size_t m=0; m < dfgShape[1]; ++m ) //w
            for( size_t k=0; k < dfgShape[3]; ++k ) //z
              for( size_t j=0; j < dfgShape[4]; ++j ) //y
                for( size_t i=0; i < dfgShape[5]; ++i ){ //x
                  assert( dfgData[n][ m ][ dfgShape[2]-1 ][k][j][i]
                          == complex(0.0) );
                }
    }
    */

    // adapt z-boundary
    adaptBoundaries(species, d);

    }
  } 
}

/*
 * This routine applies the z boundary condition
 */
void SelalibTask::adaptBoundaries(int species, int d){
  if( dfgVector_[species]->getParallelization()[d] > 1 )
    adaptBoundariesGlobal(species, d);
  else
    adaptBoundariesLocal(species, d);
}


/** the term local may be misleading in the context of GENE
 *  this simply means that the copy operations are done locally without any
 *  MPI communication. this is only possible in the case that there is only
 *  one process in z direction (and in x, what we require anyway)
 */
void SelalibTask::adaptBoundariesLocal(int species, int d){
  // make sure no parallelization in z
  assert( dfgVector_[species]->getParallelization()[2] == 1 );

  MultiArrayRef6 dfgData = createMultiArrayRef<CombiDataType,6>( *dfgVector_[species] );

  adaptBoundariesKernel(dfgData, dfgData, species, d);
}
/**
 * This method updates the last z component of target data according to the z boundary conditions.
 * The value at z[0] is taken from source data which might be a separate array in case of domain decomposition in z direction
 * @param sourceData contains the z[0] elements at z[0]
 * @param targetData needs to be updated at last z coordinate using boundary condition
 *
 */
void SelalibTask::adaptBoundariesKernel(MultiArrayRef6& sourceData, MultiArrayRef6& targetData, int species, int d){
  // calculate offset and factor
    const IndexVector p( dfgVector_[species]->getParallelization().rbegin(),
                          dfgVector_[species]->getParallelization().rend() );
    IndexVector tmp( dfgVector_[species]->getDimension() );
    dfgVector_[species]->getPartitionCoords( tmp );

    // get partition coords of the current process in the ordering as
    // the shape of the multi arrays, i.e. spec, w, v, z, y, x
    IndexVector coords( tmp.rbegin(), tmp.rend() );

    const size_t* targetShape = targetData.shape();

    std::vector<int> sizes(dim_);
    std::vector<int> offset(dim_);

    for(int i; i < dim_; i++){
      sizes[i] = i == d ? 1 : targetShape[i];
      offset[i] = i == d ? targetShape[i] - 1 : 0;
    }
    

    for( size_t n=0; n < sizes[0]; ++n ){ //w
      for( size_t m=0; m < sizes[1]; ++m ){ //v
        for( size_t l=0; l < sizes[2]; ++l ){ //u
          for( size_t k=0; k < sizes[3]; ++k ){ //z
            for( size_t j=0; j < sizes[4]; ++j ){ //y
              for( size_t i=0; i < sizes[5]; ++i ){ //x
                
                targetData[n + offset[0]][m + offset[1]][l + offset[2]][k + offset[3]][j + offset[4]][i + offset[5]] =
                    sourceData[n][m][l][k][j][i];
              } 
            }
          }
        }
      }
    }
}


/* name may be misleading, this is for local iv compuations
 * global refers to a global adaption of the z-boundary conditions, which is
 * the case when there is more than one process in this direction
 */
void SelalibTask::adaptBoundariesGlobal(int species, int d){
  // make sure at least two processes in z-direction
  assert( dfgVector_[species]->getParallelization()[d] > 1 );

  // everything in normal dfg notation here except MPI subarrays!

  // get parallelization and coordinates of process
  const IndexVector& p = dfgVector_[species]->getParallelization();
  IndexVector coords( dfgVector_[species]->getDimension() );
  dfgVector_[species]->getPartitionCoords( coords );

  // x range of current process in global coordinates
  IndexVector lBounds( dfgVector_[species]->getLowerBounds() );
  IndexVector uBounds( dfgVector_[species]->getUpperBounds() );

  IndexType xoffset;
  CombiDataType factor;

  bool sendingProc = false;
  if( coords[d] == 0 ){
    sendingProc = true;
  }
  else if( coords[d] == p[d]-1 ){
    sendingProc = false;
  } else{
    // other processes don't have to do anything here
    return;
  }

  if(species==0){
    requestArray_ = new MPI_Request[nspecies_];
    receiveBufferArray_.resize(nspecies_);
  }
  // set lower bounds of subarray
  IndexVector subarrayLowerBounds = dfgVector_[species]->getLowerBounds();

  // set upper bounds of subarray
  IndexVector subarrayUpperBounds = dfgVector_[species]->getUpperBounds();

  if( sendingProc ){
    subarrayLowerBounds[d] = 0;
    subarrayUpperBounds[d] = 1;
  } else{
    subarrayLowerBounds[d] = dfgVector_[species]->getGlobalSizes()[2] - 1;
    subarrayUpperBounds[d] = dfgVector_[species]->getGlobalSizes()[2];
  }

  // create MPI datatype
  IndexVector sizes( dfgVector_[species]->getLocalSizes() );
  IndexVector subsizes = subarrayUpperBounds - subarrayLowerBounds;
  // the starts are local indices
  IndexVector starts = subarrayLowerBounds - dfgVector_[species]->getLowerBounds();

  // convert to mpi notation
  // also we have to use int as type for the indizes
  std::vector<int> csizes(sizes.rbegin(), sizes.rend());
  std::vector<int> csubsizes(subsizes.rbegin(), subsizes.rend());
  std::vector<int> cstarts(starts.rbegin(), starts.rend());

  // create subarray view on data
  MPI_Datatype mysubarray;
  MPI_Type_create_subarray(static_cast<int>(dfgVector_[species]->getDimension()),
                           &csizes[0], &csubsizes[0], &cstarts[0],
                           MPI_ORDER_C, dfgVector_[species]->getMPIDatatype(), &mysubarray);
  MPI_Type_commit(&mysubarray);
  if(sendingProc){
    coords[d] = p[d]-1 ;
  }
  else{
    coords[d] = 0;
  }
  RankType r = dfgVector_[species]->getRank( coords );
  if( sendingProc ){
    MPI_Isend( dfgVector_[species]->getData(), 1, mysubarray, r, 1000, dfgVector_[species]->getCommunicator(),
               &requestArray_[species]);

  } else{
    IndexType numElements = 1;
    for(int i=0; i<subsizes.size(); i++){
      numElements *= subsizes[i];
    }
    receiveBufferArray_[species] = new CombiDataType[numElements];
    //receiver does not need mysubarray
    MPI_Irecv( receiveBufferArray_[species], numElements, dfgVector_[species]->getMPIDatatype(), r, 1000, dfgVector_[species]->getCommunicator(),
               &requestArray_[species]);
  }



  if(species == nspecies_ - 1){
    MPI_Waitall(nspecies_, requestArray_, MPI_STATUSES_IGNORE );
    //toDo check if correct
    //nur letze processe in z-richtung
    if(!sendingProc){
      for(int i=0; i<nspecies_; i++){
        MultiArrayRef6 dfgData = createMultiArrayRef<CombiDataType,6>( *dfgVector_[i] );
        //dfg data is stored in reverse ordering of shape
        std::vector<size_t> shape( subsizes.rbegin(),
                                       subsizes.rend() );
        MultiArrayRef6 receivedData =  MultiArrayRef<CombiDataType,6>(receiveBufferArray_[i],shape);
        adaptBoundariesKernel(receivedData,dfgData, i, d);
        delete[] receiveBufferArray_[i];
      }
    }
    delete[] requestArray_;
    receiveBufferArray_.clear();
  }
}


/* copy dfg data to local checkpoint
 * This is used to get the data from our framework to GENE
 * */
void SelalibTask::getDFG(){
  // get multiarrayref to lcp
  // todo: put this into a function as it is used multiple times
  const std::vector<size_t>& b = checkpoint_.getBounds();
  IndexVector lowerBounds = {(IndexType) b[0],(IndexType) b[2],(IndexType) b[4],(IndexType) b[6],(IndexType) b[8],(IndexType) b[10] };
  IndexVector upperBounds = {(IndexType) b[1], (IndexType)b[3],(IndexType) b[5],(IndexType) b[7],(IndexType) b[9],(IndexType) b[11] };

  IndexVector lcpSizes = upperBounds - lowerBounds;
  MultiArrayRef<GeneComplex,6> lcpData =
    createMultiArrayRef<GeneComplex,6>( checkpoint_.getData(), lcpSizes );



  // copy data back to lcp
  // note that the last grid points in x,z,v,w dimension are ignored
      // get multiarrayref to dfg
  int species = 0;
  MultiArrayRef6 dfgDataN =
  createMultiArrayRef<CombiDataType,6>( *dfgVector_[species] );
  const size_t* lcpShape = lcpData.shape();
  for( size_t n=0; n < lcpShape[0]; ++n ){ //w

    for( size_t m=0; m < lcpShape[1]; ++m ){ //v
      for( size_t l=0; l < lcpShape[2]; ++l ){ //u
        for( size_t k=0; k < lcpShape[3]; ++k ){ //z
          for( size_t j=0; j < lcpShape[4]; ++j ){ //y
            for( size_t i=0; i < lcpShape[5]; ++i ){ //x
              lcpData[n][m][l][k][j][i] = dfgDataN[n][m][l][k][j][i];
            }
          }
        }
      }
    }
  }
}


} /* namespace combigrid */
