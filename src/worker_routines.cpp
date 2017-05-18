#include "worker_routines.hpp"
#include "GeneGrid.hpp"
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include "sgpp/distributedcombigrid/fullgrid/FullGrid.hpp"
#include "CombiGeneConverter.hpp"
#include "boost/multi_array.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "GeneLocalCheckpoint.hpp"
#include "GeneTask.hpp"

double calc_time_start;

using namespace combigrid;

void checkpoint_write_memory_(GeneComplex* g_1, double *timep, double *dtp,
                              int *li1p, int *li2p, int *lj1p, int *lj2p,
                              int *lk1p, int *lk2p, int *ll1p, int *ll2p,
                              int *lm1p, int *lm2p, int *ln1p, int *ln2p,
                              int *ni0p, int *nj0p, int *nz0p,
                              int *nv0p, int *nw0p, int *n_specp,
                              MPI_Fint* comm_gene_f ) {
  double tstart = MPI_Wtime();

  Task* tt = pgroup->getCurrentTask();
  GeneTask* t = static_cast< GeneTask* >(tt);

  // this information is also contained in gene checkpoint files
  // we don't use it at the moment
	double mytime = *timep;
	double dt = *dtp;

	// calculate local size of checkpoint
  int d6 = *ln2p - *ln1p +1;
  int d5 = *lm2p - *lm1p +1;
  int d4 = *ll2p - *ll1p +1;
  int d3 = *lk2p - *lk1p +1;
  int d2 = *lj2p - *lj1p +1;
  int d1 = *li2p - *li1p +1;
  size_t size = static_cast<size_t>(d1*d2*d3*d4*d5*d6);

  /* store bounds of local arrays in following order
   * ln1/2, lm1/2, ll1/2, lk1/2, lj1/2, li1/2
   * upper bounds are +1 compared to gene
   * which corresponds to the dimensions
   * spec,w,v,z,y,x
   */
  std::vector<size_t> bounds(12);
  bounds[0] = static_cast<size_t>(*ln1p);
  bounds[1] = static_cast<size_t>(*ln2p + 1);
  bounds[2] = static_cast<size_t>(*lm1p);
  bounds[3] = static_cast<size_t>(*lm2p + 1);
  bounds[4] = static_cast<size_t>(*ll1p);
  bounds[5] = static_cast<size_t>(*ll2p + 1);
  bounds[6] = static_cast<size_t>(*lk1p);
  bounds[7] = static_cast<size_t>(*lk2p + 1);
  bounds[8] = static_cast<size_t>(*lj1p);
  bounds[9] = static_cast<size_t>(*lj2p + 1);
  bounds[10] = static_cast<size_t>(*li1p);
  bounds[11] = static_cast<size_t>(*li2p + 1);

  // get global sizes
  // order: spec,w,v,z,y,x
  std::vector<size_t> sizes(6);
  sizes[0] = static_cast<size_t>(*n_specp);
  sizes[1] = static_cast<size_t>(*nw0p);
  sizes[2] = static_cast<size_t>(*nv0p);
  sizes[3] = static_cast<size_t>(*nz0p);
  sizes[4] = static_cast<size_t>(*nj0p);
  sizes[5] = static_cast<size_t>(*ni0p);

  // todo: create decomposition vectors
  // gather lowerbounds of all procs in dfg typical ordering
  int lsize;
  MPI_Comm_size( *comm_gene_f, &lsize );
  std::vector<int> myLowerBounds = { *li1p, *lj1p, *lk1p, *ll1p, *lm1p, *ln1p };
  std::vector<int> allLowerBounds( lsize * myLowerBounds.size() );
  MPI_Allgather( &myLowerBounds[0], static_cast<int>( myLowerBounds.size() ), MPI_INT,
              &allLowerBounds[0], static_cast<int>( myLowerBounds.size() ), MPI_INT,
              *comm_gene_f );

  // extract decomposition from alllowerbounds
  std::vector<IndexVector> decomposition(6);
  for( int r=0; r<lsize; ++r ){
    for( int i=0; i<6; ++i )
      decomposition[i].push_back( allLowerBounds[ r*6 + i ] );
  }

  // decomposition contains only unique elements
  for( int d=0; d<6; ++d ){
    std::sort( decomposition[d].begin(),
               decomposition[d].end() );
    IndexVector::iterator last = std::unique( decomposition[d].begin(),
                                              decomposition[d].end());
    decomposition[d].erase( last, decomposition[d].end() );
  }

  // check if dfg created and create if necessary
  t->initDFG( *comm_gene_f, decomposition );

  // set decomposition in combiparameters
  CombiParameters& param = pgroup->getCombiParameters();
  param.setApplicationComm( *comm_gene_f );

  t->writeLocalCheckpoint( g_1, size, sizes, bounds );
  //t->setTimeCPMem( MPI_Wtime() - tstart );

  t->setDFG();
}

void checkpoint_read_memory_(GeneComplex* g_1, int *li1p, int *li2p,
                             int *lj1p, int *lj2p, int *lk1p, int *lk2p,
                             int *ll1p, int *ll2p, int *lm1p, int *lm2p,
                             int *ln1p, int *ln2p) {
  Task* tt = pgroup->getCurrentTask();
  GeneTask* t = static_cast< GeneTask* >(tt);

  // write dfg data to local checkpoint
  t->getDFG();

  GeneLocalCheckpoint& cp = t->getLocalCheckpoint();

  assert( cp.getSize() > 0 );

  // check local size of checkpoint
  int d6 = *ln2p - *ln1p +1;
  int d5 = *lm2p - *lm1p +1;
  int d4 = *ll2p - *ll1p +1;
  int d3 = *lk2p - *lk1p +1;
  int d2 = *lj2p - *lj1p +1;
  int d1 = *li2p - *li1p +1;
  size_t size = static_cast<size_t>(d1*d2*d3*d4*d5*d6);

  /*
  std::cout << "expected cp size " << size
            << " actual cp size " << cp.getSize() << std::endl;
            */

  assert( cp.getSize() == size );

  // check bounds
  std::vector<size_t> bounds(12);
  bounds[0] = static_cast<size_t>(*ln1p);
  bounds[1] = static_cast<size_t>(*ln2p + 1);
  bounds[2] = static_cast<size_t>(*lm1p);
  bounds[3] = static_cast<size_t>(*lm2p + 1);
  bounds[4] = static_cast<size_t>(*ll1p);
  bounds[5] = static_cast<size_t>(*ll2p + 1);
  bounds[6] = static_cast<size_t>(*lk1p);
  bounds[7] = static_cast<size_t>(*lk2p + 1);
  bounds[8] = static_cast<size_t>(*lj1p);
  bounds[9] = static_cast<size_t>(*lj2p + 1);
  bounds[10] = static_cast<size_t>(*li1p);
  bounds[11] = static_cast<size_t>(*li2p + 1);

  const std::vector<size_t>& cp_bounds = cp.getBounds();
  for( size_t i=0; i<bounds.size(); ++i )
    assert( bounds[i] == cp_bounds[i] );

  memcpy ( g_1, cp.getData(), cp.getSize() * sizeof(GeneComplex) );
}


volatile int endStallForDebugger=1;

void stallForDebugger()
{
        while (!endStallForDebugger) ;
}


void worker_wait( MPI_Comm lcomm, int* worker_stat, int nprocs, int ngroup ){
  // first time wait function is called: mpi system and pgroup are initialized
  if( !theMPISystem()->isInitialized() )
    theMPISystem()->init( ngroup, nprocs, lcomm );

  if( pgroup == NULL )
    pgroup = new ProcessGroupWorker();

  //when unfinished tasks -> run next task
  // sonst pgroup wait

  SignalType signal = pgroup->wait();

  *worker_stat = signal;
}


//fortran interface
void worker_wait_(MPI_Fint* lcomm_f, int* worker_stat, int* nprocs, int* ngroup ){
  worker_wait( (MPI_Comm) *lcomm_f, worker_stat, *nprocs, *ngroup );
}


void worker_ready(double wtime, double time_perf,
                  double time_iv, double time_cp){
  MASTER_EXCLUSIVE_SECTION {
    // copy parameters.dat to parameters
    rename( "parameters.dat", "parameters" );
  }

  // get current task
  Task* t = pgroup->getCurrentTask();

  // set finished status
  t->setFinished(true);

  pgroup->ready();
}


// fortran interface
void worker_ready_( double* wtime, double* time_perf,
                    double* time_iv, double* time_cp ){
  worker_ready(*wtime, *time_perf, *time_iv, *time_cp );
}


void write_omega_(int* itime, double* gamma, double* omega){
	//open omega_out with append
	//std::ofstream myfile( "omega_out2.dat", std::ios::app );
	//myfile << *itime << " " << *gamma << " " << *omega << std::endl;
	//myfile.close();

  //Task* t = pgroup->getCurrentTask();
  //t->setLambda( *itime, combigrid::complex( *gamma, *omega ) );
}


/* get access to the nrg 0 value of the first species when it is written
 * by gene. to get the correct value at the combination step,
 * make sure 'istep_nrg' in the parameters divides through the number
 * of time steps */
void set_nrg_(double* time, double* nrg0){
  Task* tt = pgroup->getCurrentTask();
  GeneTask* t = static_cast< GeneTask* >(tt);

  t->setNrg(*nrg0);
}


void init_stats_(){
  Stats::initialize();
}

void finalize_stats_(){
  Stats::finalize();

  /* write stats to json file for postprocessing */
  Stats::write( "timers.json" );
}

void set_group_id_(int* color){
  Stats::setAttribute("group", std::to_string(*color));
}


void gene_time_start_(){
  Stats::startEvent("worker run");
}

void gene_time_stop_(){
  Stats::stopEvent("worker run");
}
