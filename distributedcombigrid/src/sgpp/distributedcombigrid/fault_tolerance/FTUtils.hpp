#ifndef HELPER_HPP_
#define HELPER_HPP_

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric> 
#include <vector>
#include <map>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <valarray>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"

namespace combigrid {

typedef std::vector<std::vector<real>> matrix;
typedef std::map<int, LevelVector> IdToLevelDict;
typedef std::map<LevelVector, real> CombigridDict;
typedef std::vector<LevelVector> LevelVectorList;

/* used to convert a string to a number in any format */
template<typename T>
T str_to_number(const std::string& no);

/* used to remove a vector element of vec from position pos */
template<typename T>
void remove(std::vector<T>& vec, size_t pos);

/* used for calling the python code as python script_name level_min level_max */
std::string python_code_caller(const std::string& script_name,
    const LevelVectorList& levels, const int& dim);

/* used to get data for the GCP when minimizing the interpolation error */
CombigridDict get_python_data(const std::string& script_run, const int& dim);

/* used to create the M matrix for the interpolation based problem */
matrix create_M_matrix(const CombigridDict& aux_downset, const int& dim);

/* used to create the inverse of M */
matrix get_inv_M(const CombigridDict& aux_downset, const int& dim);

/* used to create the inverse of M matrix for the interpolation based problem */
CombigridDict set_entire_downset_dict(const LevelVectorList levels,
    const CombigridDict& received_dict, const int& dim);

/* used to get a vector of the entire downset indices */
LevelVectorList get_downset_indices(const CombigridDict& entire_downset,
    const int& dim);

/* used to filter the input LevelVectorList faults such that only faults from the partial downset
 (input from python) are considered */
LevelVectorList filter_faults(const LevelVectorList& faults_input, const IndexType& l_max,
    const CombigridDict& received_dict);

/* used to create an entire downset dictionary used for setting up the M matrix */
CombigridDict create_aux_entire_dict(const CombigridDict& entire_downset,
    const int& dim);

/* used to print the new dictionary after the optimization is performed */
CombigridDict create_out_dict(const CombigridDict& given_downset,
    const std::vector<real>& new_c, const int& dim);

/* used to set the row and column variables for the optimization problem */
std::string set_aux_var_name(const std::string& var_name, const int& index);

/* generates no_of_levels random faults */
int generate_random_fault(const int& no_of_levels);

/* used to generate random variables for the W matrix in the optimization problem */
std::vector<real> gen_rand(const int& size);

/* used to compute the size of the downset */
int get_size_downset(const LevelVector& level_max, const int& dim);

/* used to compute the L1 norm of a vector */
int l1_norm(const LevelVector& u);

/* used to compute factorial; needed to compute size of the downset */
int factorial(const int& dim);

/* test whether b >= a */
bool test_greater(const LevelVector& b, const LevelVector& a);

/* used to create a multi-index based on maximum level */
LevelVectorList mindex(const int& dimension, const LevelVector& level_max);

/* used to check input levels dimensionality */
LevelVectorList check_dimensionality(const LevelVectorList& input_levels,
    LevelVector& ignored_dimensions);

/* used to ignore certain dimensions of the input faults based on input levels & ignored dimension */
LevelVectorList check_faults(const LevelVectorList& input_faults,
    const LevelVector& ignored_dimensions);

/* used to create a new dictionary, based on the given dictionary and the ignored dimensions */
CombigridDict set_new_given_dict(const CombigridDict& given_dict,
    const LevelVector& ignored_dimensions, const int& dim);

/* used to check whether the input levels are correct */
/* i.e. they satisfy: l_max - l_min = c*ones(dim) */
void check_input_levels(const LevelVectorList& levels);

/* used to extract only the coefficients output by the optimization problem */
/* which correspond to levels from the given downset */
std::vector<real> select_coeff_downset(const std::vector<real>& all_c,
    const CombigridDict& given_downset, const CombigridDict& aux_downset);


/* contains information about fault simulation:
** numFaults_: how many faults occur
** iterationFaults_: vector of timesteps at which processes fail
** globalRankFaults_: global rank of process that fails
**/
struct FaultsInfo {
public:
  int numFaults_;
  IndexVector iterationFaults_;
  IndexVector globalRankFaults_;
  IndexVector taskFaults_;
private:
  friend class boost::serialization::access;
  // serialize
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
};


template<class Archive>
void FaultsInfo::serialize(Archive& ar, const unsigned int version) {
  ar& numFaults_;
  ar& iterationFaults_;
  ar& globalRankFaults_;
  ar& taskFaults_;
}

} // namespace combigrid
#endif /* HELPER_HPP_ */
