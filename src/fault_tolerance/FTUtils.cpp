#include "fault_tolerance/FTUtils.hpp"
#include <random>
#include <numeric>
#include <valarray>

namespace combigrid {
template <typename T>
T str_to_number(const std::string& no) {
  T value;
  std::stringstream stream(no);
  stream >> value;

  if (stream.fail()) {
    std::runtime_error e(no);
    std::cout << "Error in the conversion of " << no << "!" << std::endl;
    throw e;
  }

  return value;
}

std::string python_code_caller(const std::string& script_name, const LevelVectorList& levels,
                               const int& dim) {
  LevelType levels_no = static_cast<LevelType>(levels.size());
  LevelType level_size = static_cast<LevelType>(levels[0].size());
  LevelType one_level = 0;
  std::stringstream caller;

  caller << "python " << script_name << " " << dim << " ";

  for (int i = 0; i < levels_no; ++i) {
    for (int j = 0; j < level_size; ++j) {
      one_level = levels[i][j];
      caller << one_level << " ";
    }
  }

  return caller.str();
}

CombigridDict get_python_data(const std::string& script_run, const int& dim) {
  FILE* stream;
  char buffer[256];
  std::string level_x_str;
  std::string level_y_str;
  std::string coeff_str;

  double coeff = 0.0;

  CombigridDict dict;

  stream = popen(script_run.c_str(), "r");

  if (stream) {
    while (!feof(stream)) {
      if (fgets(buffer, sizeof(buffer), stream) != NULL) {
        std::string one_level_str;
        int one_level = 0;
        LevelVector levels;
        std::stringstream temp(buffer);

        for (int i = 0; i < dim; ++i) {
          temp >> one_level_str;
          one_level = str_to_number<int>(one_level_str);
          levels.push_back(one_level);
        }

        temp >> coeff_str;
        coeff = str_to_number<double>(coeff_str);

        dict.insert(std::make_pair(levels, coeff));
      }
    }

    pclose(stream);
  } else {
    throw "Error reading script output!";
  }

  return dict;
}

matrix get_inv_M(const CombigridDict& aux_downset, const int& dim) {
  int size_downset = static_cast<int>(aux_downset.size());
  int i = 0;
  int j = 0;

  std::valarray<int> c(dim);
  std::valarray<int> w(dim);
  std::valarray<int> diff(dim);
  std::valarray<int> zeros(0, dim);

  matrix M_inv(size_downset, std::vector<real>(size_downset, 0.0));

  for (auto ii = aux_downset.begin(); ii != aux_downset.end(); ++ii) {
    i = static_cast<int>(ii->second);
    for (int it = 0; it < dim; ++it) {
      w[it] = static_cast<int>(ii->first[it]);
    }

    for (auto jj = ii; jj != aux_downset.end(); ++jj) {
      j = static_cast<int>(jj->second);
      for (int it = 0; it < dim; ++it) {
        c[it] = static_cast<int>(jj->first[it]);
      }

      diff = c - w;

      if (((diff.sum() > 0) || (diff.sum() <= dim)) &&
          ((diff.max() <= 1) && (diff >= zeros).min())) {
        M_inv[i][j] = pow(-1, diff.sum());
      }
    }
  }

  return M_inv;
}

CombigridDict set_entire_downset_dict(const LevelVectorList levels,
                                      const CombigridDict& received_dict, const int& dim) {
  LevelVector level_min = levels.front();
  LevelVector level_max = levels.back();
  CombigridDict active_downset;
  CombigridDict output;

  LevelVector level_active_downset;

  LevelVector level;
  LevelVectorList all_levels;
  LevelVectorList feasible_levels;

  double key = 0.0;

  all_levels = mindex(dim, level_max);

  for (auto ii = received_dict.begin(); ii != received_dict.end(); ++ii) {
    if (ii->second > 0.0) {
      active_downset.insert(std::make_pair(ii->first, ii->second));
    }
  }

  for (auto ii = active_downset.begin(); ii != active_downset.end(); ++ii) {
    level_active_downset = ii->first;

    for (unsigned int i = 0; i < all_levels.size(); ++i) {
      if (test_greater(level_active_downset, all_levels[i]) &&
          test_greater(all_levels[i], level_min)) {
        feasible_levels.push_back(all_levels[i]);
      }
    }
  }

  for (unsigned int i = 0; i < feasible_levels.size(); ++i) {
    level = feasible_levels[i];
    auto ii = received_dict.find(level);

    if (ii != received_dict.end()) {
      key = ii->second;
      output.insert(std::make_pair(level, key));
    } else {
      key = 0.0;
      output.insert(std::make_pair(level, key));
    }
  }

  return output;
}

CombigridDict create_aux_entire_dict(const CombigridDict& entire_downset, const int& dim) {
  real key = 0;
  int i = 0;

  CombigridDict aux_dict;

  for (auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii) {
    LevelVector levels;
    key = static_cast<real>(i);

    for (int i = 0; i < dim; ++i) {
      levels.push_back(ii->first[i]);
    }

    ++i;
    aux_dict.insert(std::make_pair(levels, key));
  }

  return aux_dict;
}

LevelVectorList get_downset_indices(const CombigridDict& entire_downset, const int& dim) {
  LevelVectorList indices;

  for (auto ii = entire_downset.begin(); ii != entire_downset.end(); ++ii) {
    LevelVector index;

    for (int i = 0; i < dim; ++i) {
      index.push_back(ii->first[i]);
    }

    indices.push_back(index);
  }

  return indices;
}

LevelVectorList filter_faults(const LevelVectorList& faults_input, const IndexType& l_max,
                              const CombigridDict& received_dict) {
  int no_faults = 0;
  int level_fault = 0;
  LevelVectorList faults_output;

  no_faults = static_cast<int>(faults_input.size());

  for (int i = 0; i < no_faults; ++i) {
    auto it = received_dict.find(faults_input[i]);

    if (it != received_dict.end()) {
      level_fault = std::accumulate(faults_input[i].begin(), faults_input[i].end(), 0);

      if ((level_fault == l_max) || (level_fault == (l_max - 1))) {
        faults_output.push_back(faults_input[i]);
      }
    }
  }

  return faults_output;
}

CombigridDict create_out_dict(const CombigridDict& given_downset, const std::vector<real>& new_c,
                              const int& dim) {
  real key = 0;
  int i = 0;

  CombigridDict out_dict;

  for (auto ii = given_downset.begin(); ii != given_downset.end(); ++ii) {
    LevelVector levels;
    key = new_c[i];

    for (int i = 0; i < dim; ++i) {
      levels.push_back(ii->first[i]);
    }

    ++i;
    out_dict.insert(std::make_pair(levels, key));
  }

  return out_dict;
}

std::string set_aux_var_name(const std::string& var_name, const int& index) {
  std::stringstream aux_var;
  aux_var << var_name << index;

  return aux_var.str();
}
// generate integer random numbers between 0 and #levels-1
int generate_random_fault(const int& no_of_levels) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<int> rand_num(0, no_of_levels - 1);

  return rand_num(rng);
}

int l1_norm(const LevelVector& u) {
  int norm = 0;

  for (auto elem : u) {
    norm += abs(static_cast<int>(elem));
  }

  return norm;
}

bool test_greater(const LevelVector& b, const LevelVector& a) {
  int dim = static_cast<int>(a.size());
  bool test = true;

  for (int i = 0; i < dim; ++i) {
    test *= (b[i] >= a[i]) ? true : false;
  }

  return test;
}

LevelVectorList mindex(const int& dimension, const LevelVector& level_max) {
  int j = 0;
  int norm = 0;
  IndexType sum_level_max = 0;

  LevelVector temp(dimension, 1);
  LevelVectorList mindex_result;

  auto upper_limit = std::max_element(level_max.begin(), level_max.end());

  sum_level_max = std::accumulate(level_max.begin(), level_max.end(), 0);

  while (true) {
    norm = l1_norm(temp);

    if (norm <= sum_level_max) {
      mindex_result.push_back(temp);
    }

    for (j = dimension - 1; j >= 0; --j) {
      if (++temp[j] <= *upper_limit)
        break;
      else
        temp[j] = 1;
    }

    if (j < 0) break;
  }

  return mindex_result;
}

LevelVectorList check_dimensionality(const LevelVectorList& input_levels,
                                     LevelVector& ignored_dimensions) {
  LevelVector l_min = input_levels[0];
  LevelVector l_max = input_levels[1];

  LevelVector new_l_min;
  LevelVector new_l_max;
  LevelVectorList new_levels;

  for (size_t i = 0; i < l_min.size(); ++i) {
    if (l_max[i] == l_min[i]) {
      ignored_dimensions.push_back(i);
    } else {
      new_l_min.push_back(l_min[i]);
      new_l_max.push_back(l_max[i]);
    }
  }

  new_levels.push_back(new_l_min);
  new_levels.push_back(new_l_max);

  return new_levels;
}

LevelVectorList check_faults(const LevelVectorList& input_faults,
                             const LevelVector& ignored_dimensions) {
  LevelVectorList new_faults;

  for (unsigned int i = 0; i < input_faults.size(); ++i) {
    LevelVector new_fault;

    for (unsigned int j = 0; j < input_faults[0].size(); ++j) {
      if (std::find(ignored_dimensions.begin(), ignored_dimensions.end(), j) ==
          ignored_dimensions.end()) {
        new_fault.push_back(input_faults[i][j]);
      }
    }

    new_faults.push_back(new_fault);
  }

  return new_faults;
}

CombigridDict set_new_given_dict(const CombigridDict& given_dict,
                                 const LevelVector& ignored_dimensions, const int& dim) {
  real key = 0.0;
  CombigridDict new_given_dict;

  for (const auto & ii : given_dict) {
    LevelVector new_level;

    for (int i = 0; i < dim; ++i) {
      if (std::find(ignored_dimensions.begin(), ignored_dimensions.end(), i) ==
          ignored_dimensions.end()) {
        new_level.push_back(ii.first[i]);
      }
    }

    key = ii.second;
    new_given_dict.insert(std::make_pair(new_level, key));
  }

  return new_given_dict;
}

void check_input_levels(const LevelVectorList& levels) {
  LevelVector l_min = levels[0];
  LevelVector l_max = levels[1];
  LevelVector c;

  // cf. https://stackoverflow.com/questions/3642700/vector-addition-operation
  c.reserve( l_min.size() );
  // for (unsigned int i = 0; i < l_min.size(); ++i) {
  //   c.push_back(l_max[i] - l_min[i]);
  // }
  std::transform(l_max.begin(), l_max.end(), l_min.begin(), std::back_inserter(c), std::minus<IndexType>());

  assert (std::adjacent_find(c.begin(), c.end(), std::not_equal_to<int>()) == c.end() &&
    "Input levels are incorrect!" &&
    "Please input them of the form: l_max = l_min + c*ones(dim), c>=1, integer"
  );
}

std::vector<double> select_coeff_downset(const std::vector<double>& all_c,
                                         const CombigridDict& given_downset,
                                         const CombigridDict& aux_downset) {
  int given_downset_index = 0;
  std::vector<double> donwset_c;

  for (const auto & ii : aux_downset) {
    if (given_downset.find(ii.first) != given_downset.end()) {
      given_downset_index = static_cast<int>(ii.second);
      donwset_c.push_back(all_c.at(given_downset_index));
    }
  }

  return donwset_c;
}

}  // namespace combigrid
