#ifndef LPOPTINTERP_HPP_
#define LPOPTINTERP_HPP_

#include "sgpp/distributedcombigrid/fault_tolerance/LPOptimization.hpp"

namespace combigrid {

class LP_OPT_INTERP : public LP_OPT {
 private:
  /* levels of grid indices */
  LevelVectorList i_levels;
  /* top level of grid indices*/
  LevelVector level_max;
  /* dimension of the problem */
  int i_dim;
  /* total size of the optimization problem */
  int total_size;

  /* new dimensionality after the input levels are checked */
  int new_dim;
  /* new levels based on new_dim */
  LevelVectorList new_levels;
  /* if it is the case, the dimension(s) that is/are ignored */
  LevelVector ignored_dimensions;
  /* if it is the case, new faults, based on ignored dimensions */
  LevelVectorList new_faults;

  /* no. of constraints */
  int no_faults;
  /* level max sum */
  IndexType l_max;

  //  total size of the optimization problem
  /* down set size */
  int size_downset;

  /* inverse of M s.t. w = Mc*/
  matrix inv_M;

  /* given dictionary */
  CombigridDict given_downset;
  /* modified given downset based on ignored dimensions */
  CombigridDict new_given_downset;
  /* entire donwset with corresponding indices */
  CombigridDict entire_downset;
  /* auxiliary dictionary, used to create M and inv(M) */
  CombigridDict aux_entire_dict;
  /* downset indices as a 2d vector*/
  LevelVectorList downset_indices;
  /* faults input by user */
  LevelVectorList input_faults;
  /* faults that are in the given downset from the set of input faults */
  LevelVectorList valid_input_faults;
  /* faults that have to be recomputed */
  mutable LevelVectorList recompute_faults;

 public:
  LP_OPT_INTERP();

  LP_OPT_INTERP(const LevelVectorList& _levels, const int& _dim, const int& _opt_type,
                const CombigridDict& _given_downset, const LevelVectorList& _input_faults);

  LP_OPT_INTERP(const LP_OPT_INTERP& obj);

  LP_OPT_INTERP& operator=(const LP_OPT_INTERP& rhs);

  virtual void init_opti_prob(const std::string& prob_name);

  virtual void set_constr_matrix();

  virtual void set_constr_matrix(const std::vector<real>& W);

  virtual int getNumFaults();

  virtual void solve_opti_problem() const;

  virtual CombigridDict get_results(LevelVectorList& recomp_faults) const;

  virtual ~LP_OPT_INTERP();
};
}
#endif /* LPOPTINTERP_HPP_ */
