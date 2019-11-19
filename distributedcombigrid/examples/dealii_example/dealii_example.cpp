/*#include <deal.II/base/mpi.h>
#include <hyper.deal/include/mpi.h>

using namespace dealii;
using namespace hyperdealii;
int main(int argc, char **argv){
  Utilities::MPI::MPI_InitFinalize mpi(argc,argv,1);
}
*/
// deal.II
#include <deal.II/base/mpi.h>
#include <deal.II/base/revision.h>

#include <hyper.deal.combi/include/functionalities/dynamic_convergence_table.h>

#ifdef LIKWID_PERFMON
#  include <likwid.h>
#else
#  define LIKWID_MARKER_INIT
#  define LIKWID_MARKER_THREADINIT
#  define LIKWID_MARKER_SWITCH
#  define LIKWID_MARKER_REGISTER(regionTag)
#  define LIKWID_MARKER_START(regionTag)
#  define LIKWID_MARKER_STOP(regionTag)
#  define LIKWID_MARKER_CLOSE
#  define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

using namespace dealii;

const unsigned int dim    = 2;
const unsigned int degree = 1; /*dummy value*/
#define degree_min 1
#define degree_max 15
const MPI_Comm comm = MPI_COMM_WORLD;

typedef double                     Number;
typedef VectorizedArray<Number, 1> VectorizedArrayType;

// include is only possible here since constants are also used in those files (WIP)
//#include "../../include/functionalities/vector.h"
//#include "../../include/functionalities/vector_mf.h"
#include <hyper.deal.combi/include/functionalities/vector_dummy.h>

// application
#include <hyper.deal.combi/applications/advection_reference_dealii/include/application.h>
#include <hyper.deal.combi/applications/advection_reference_dealii/include/parameters_driver.h>

typedef LinearAlgebra::distributed::Vector<Number> VectorType;
// typedef LinearAlgebra::sharedmpi::Vector<Number> VectorType;

template<int dim, int degree>
void
run_degree(const std::string & file_name)
{
  typedef Application<dim, degree, degree + 1, Number, VectorizedArrayType, VectorType> Problem;

  // create table for measurements
  hyperdeal::DynamicConvergenceTable table;

  // process problem
  Problem problem(comm, table);
  problem.reinit(file_name);
  problem.solve();
  // print results
  table.print(false);
}

int
main(int argc, char ** argv)
{
  try
  {
    Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

    if(Utilities::MPI::this_mpi_process(comm) == 0)
    {
      printf("deal.II git version %s on branch %s\n", DEAL_II_GIT_SHORTREV, DEAL_II_GIT_BRANCH);
      printf("lanes = %d\n\n", VectorizedArrayType::n_array_elements);
    }

    deallog.depth_console(0);

    // process command line arguments
    AssertThrow(argc >= 2, ExcMessage("You have not provided enough command-line arguments."));

    ParametersDriver param;

    std::string file_name(argv[1]);
    param.deserialize(file_name);
  std::cout << param.dim << file_name << std::endl;
    // check input parameters
    AssertThrow(dim == param.dim, ExcMessage("dim is currently fixed at compile time!"));


    LIKWID_MARKER_INIT;
    LIKWID_MARKER_THREADINIT;

    switch(param.degree)
    {
#if degree_min <= 1 && 1 <= degree_max
      case 1:
        run_degree<dim, 1>(file_name);
        break;
#endif
#if degree_min <= 2 && 2 <= degree_max
      case 2:
        run_degree<dim, 2>(file_name);
        break;
#endif
#if degree_min <= 3 && 3 <= degree_max
      case 3:
        run_degree<dim, 3>(file_name);
        break;
#endif
#if degree_min <= 4 && 4 <= degree_max
      case 4:
        run_degree<dim, 4>(file_name);
        break;
#endif
#if degree_min <= 5 && 5 <= degree_max
      case 5:
        run_degree<dim, 5>(file_name);
        break;
#endif
#if degree_min <= 6 && 6 <= degree_max
      case 6:
        run_degree<dim, 6>(file_name);
        break;
#endif
#if degree_min <= 7 && 7 <= degree_max
      case 7:
        run_degree<dim, 7>(file_name);
        break;
#endif
#if degree_min <= 8 && 8 <= degree_max
      case 8:
        run_degree<dim, 8>(file_name);
        break;
#endif
#if degree_min <= 9 && 9 <= degree_max
      case 9:
        run_degree<dim, 9>(file_name);
        break;
#endif
#if degree_min <= 10 && 10 <= degree_max
      case 10:
        run_degree<dim, 10>(file_name);
        break;
#endif
#if degree_min <= 11 && 11 <= degree_max
      case 11:
        run_degree<dim, 11>(file_name);
        break;
#endif
#if degree_min <= 12 && 12 <= degree_max
      case 12:
        run_degree<dim, 12>(file_name);
        break;
#endif
#if degree_min <= 13 && 13 <= degree_max
      case 13:
        run_degree<dim, 13>(file_name);
        break;
#endif
#if degree_min <= 14 && 14 <= degree_max
      case 14:
        run_degree<dim, 14>(file_name);
        break;
#endif
#if degree_min <= 15 && 15 <= degree_max
      case 15:
        run_degree<dim, 15>(file_name);
        break;
#endif
      default:
        AssertThrow(false, ExcMessage("This degree is not implemented!"));
    }

    LIKWID_MARKER_CLOSE;
  }
  catch(std::exception & exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }
  return 0;
}
