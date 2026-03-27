#ifndef DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_
#define DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_

#include <complex>

#ifndef NDEBUG
#include <boost/assert.hpp>
#define ASSERT(cond, msg) {\
    if(!(cond))\
    {\
        std::stringstream str;\
        str << msg;\
        std::cout << msg;\
        std::cerr << msg;\
        BOOST_ASSERT_MSG(cond, str.str().c_str());\
    }\
}
#else // NDEBUG
#define ASSERT(cond, msg)
#endif

namespace combigrid {

/* With this config class the distributed combigrid module can be configured
 * for a specific application.
 */

/* set the datatype for floating point numbers. usually this would be float or
 * double.
 */
typedef double real;

// the datatype for complex numbers will change accordingly. do not modify this.
typedef std::complex<real> complex;

/* using a uniform domain decomposition for all component grids (the same
 * number of processes in each dimension) yields a significantly better performance
 * for the combination and eval operation.
 * to enable the uniform operations set this to true.
 * so far, only the uniform operations are properly implemented
 */
/* switch on fault tolerance functionality */
#ifdef ENABLEFT
constexpr bool ENABLE_FT = true;  // TODO move this switch to a more sensible place
#else
constexpr bool ENABLE_FT = false;
#endif
typedef real CombiDataType;

}  // namespace combigrid

#endif /* DISTRIBUTEDCOMBIGRID_SRC_SGPP_DISTRIBUTEDCOMBIGRID_UTILS_CONFIG_HPP_ */
