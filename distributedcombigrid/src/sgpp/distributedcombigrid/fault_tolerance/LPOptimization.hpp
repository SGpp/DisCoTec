#ifndef LPOPT_HPP_
#define LPOPT_HPP_

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "sgpp/distributedcombigrid/fault_tolerance/FTUtils.hpp"

#include "glpk.h"

namespace combigrid {

class LP_OPT {
 protected:
  /* optimization type: GLP_MIN or GLP_MAX */
  int opt_type;

  /* glp problem; used in every glpk function */
  glp_prob* i_lp_prob;

  /* constraint matrix */
  std::vector<double> constr_mat;
  /* row index vector for the constraint matrix */
  std::vector<int> row_index;
  /* colums index vector for the constraint matrix */
  std::vector<int> col_index;

 public:
  /* used for LP optimization problem initialization */
  /* number of auxiliary and structural variables is set */
  /* as well as objective function and the constraints */
  virtual void init_opti_prob(const std::string& prob_name) = 0;

  /* used to set up the constraint matrix and its row and column indices*/
  virtual void set_constr_matrix() = 0;

  /* used to set up the constraint matrix and its row and column indices when it depends on matrix
   * W*/
  virtual void set_constr_matrix(const std::vector<real>& W) = 0;

  /* used to solve the linear programming problem, using the simplex algorithm */
  virtual void solve_opti_problem() const = 0;

  /* used to get the output of the LP problem, i.e. c and d vectors */
  virtual CombigridDict get_results(LevelVectorList& faults) const = 0;

  /* destructor; used to deallocate all the allocated memory */
  virtual ~LP_OPT() {}
};
}
#endif /* LPOPT_HPP_ */
