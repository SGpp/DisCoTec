#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include <unistd.h>
#include "GeneTask.hpp"
#include "sgpp/distributedcombigrid/manager/ProcessGroupWorker.hpp"
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "GeneGrid.hpp"

#define LOCAL_ROOT 0
#define DEBUG 0

using namespace combigrid;

static ProcessGroupWorker* pgroup = NULL;

extern "C" void checkpoint_write_memory_(
                                  GeneComplex* g_1, double *timep, double *dtp,
                                  int *li1p, int *li2p, int *lj1p, int *lj2p,
                                  int *lk1p, int *lk2p, int *ll1p, int *ll2p,
                                  int *lm1p, int *lm2p, int *ln1p, int *ln2p,
                                  int *ni0p, int *nj0p, int *nz0p,
                                  int *nv0p, int *nw0p, int *n_spec,
                                  MPI_Fint* comm_gene_f );

extern "C" void checkpoint_read_memory_(GeneComplex* g_1, int *li1, int *li2,
                                        int *lj1, int *lj2, int *lk1, int *lk2,
                                        int *ll1, int *ll2, int *lm1, int *lm2,
                                        int *ln1, int *ln2, int *ni0p, int *nj0p, int *nz0p,
                                        int *nv0p, int *nw0p, int *n_specp);

void write_to_file(MPI_Comm lcomm);

//combigrid::complex eval( MPI_Comm lcomm, std::vector<combigrid::real>& x );

// fortran interfaces
extern "C" void worker_wait_(MPI_Fint* comm_gene_f, int* worker_stat,
                             int* nprocs, int* ngroup );
extern "C" void worker_ready_(double* wtime, double* time_perf,
                              double* time_iv, double* time_cp );

// append gamma and omega to omega_out.dat
extern "C" void write_omega_(int* itime, double* gamma, double* omega);

extern "C" void set_nrg_(double* time, double* nrg0 );

extern "C" void mpi_ft_init_();

extern "C" void mpi_ft_finalize_();


extern "C" void decide_to_kill_();

// c interfaces
void worker_wait(MPI_Comm comm, int* worker_stat, int nprocs, int ngroup );
void worker_ready(double wtime, double time_perf,
                  double time_iv, double time_cp );

