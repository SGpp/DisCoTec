/*
 * Config.hpp
 *
 *  Created on: Feb 26, 2016
 *      Author: heenemo
 */

#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_

#include <complex>

#define ASSERT(cond, msg) {\
    if(!(cond))\
    {\
        std::stringstream str;\
        str << msg;\
        std::cerr << msg;\
        BOOST_ASSERT_MSG(cond, str.str().c_str());\
    }\
}

namespace combigrid{

/* With this config class the distributed combigrid module can be configured
 * for a specific application.
 */

/* set the datatype for floating point numbers. usually this would be float or
 * double.
 */
typedef double real;

// the datatype for complex numbers will change accordingly. do not modify this.
typedef std::complex<real> complex;

/* set the datatype for the values stored in any type of grid. essentially you
 * have two options: real values or complex numbers. other datatypes like int
 * have not been tested and operations on the grids like evaluation or
 * hierarchization might produce unexpected results.
 */
//typedef real CombiDataType;
typedef complex CombiDataType;

/* nonblocking mpi collective calls (MPI_Iallreduce and the likes) usually yield
 * better performance in some of the operations in CombiCom. if you observe
 * problems with these functions uncomment to fall back to blocking counterpart
 * of the function.
 */
constexpr bool USE_NONBLOCKING_MPI_COLLECTIVE = true;

/* for some applications it is necessary to send the ready signal while the
 * process is in the application code. in this case this flag can be set to
 * true to avoid that the ready signal is sent automatically.
 */
constexpr bool omitReadySignal = true;

/* using a uniform domain decomposition for all component grids (the same
 * number of processes in each dimension) yields a significantly better performance
 * for the combination and eval operation.
 * to enable the uniform operations set this to true.
 * so far, only the uniform operations are properly implemented
 */
constexpr bool uniformDecomposition = true;

/* switch on fault tolerance functionality */
const bool ENABLE_FT = false;
const bool isGENE = true;
//const bool GENE_Global = true;
//const bool GENE_Linear = true;
}

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_ */
