#include <assert.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <vector>
#include "boost/lexical_cast.hpp"
#include "fullgrid/FullGrid.hpp"
#include "fullgrid/MultiArray.hpp"
#include "utils/Config.hpp"
#include "utils/LevelVector.hpp"
#include "utils/Types.hpp"
// #include "timing.h"

using namespace combigrid;

void readPlotFile(const char* pltFileName, std::vector<CombiDataType>& data,
                  IndexVector& resolution);

void calcNorms(std::vector<CombiDataType>& dleft, std::vector<CombiDataType>& dright,
               std::vector<real>& Norms);

real l2Norm(std::vector<CombiDataType>& data);

int main(int argc, char** argv) {
  std::cout << argc << "\n";
  assert(argc == 6);

  // mode normalize abs : either normalize w.r.t. l2 Norm of grids, or don't
  char* mode = argv[1];

  // files
  std::string filenameLeft(argv[2]);
  std::string filenameRight(argv[3]);
  std::string filenameError(argv[4]);
  std::string prefix(argv[5]);

  std::vector<CombiDataType> data1;
  std::vector<CombiDataType> data2;

  IndexVector res1;
  IndexVector res2;

  // get data and res of first file
  readPlotFile(filenameLeft.c_str(), data1, res1);

  // get data and res of second file
  readPlotFile(filenameRight.c_str(), data2, res2);

  // check sizes
  assert(data1.size() == data2.size());
  assert(res1 == res2);

  // normalize with l2 norm
  real l2norm1 = l2Norm(data1);
  real l2norm2 = l2Norm(data2);
  std::cout << "l2 norm grid 1: " << l2norm1 << " l2 norm grid 2: " << l2norm2 << "\n";

  real tmp1 = 1.0 / l2norm1;
    real tmp2 = 1.0 / l2norm2;
  if (mode[0] == 'n') {
    for (auto i = 0; i < data1.size(); ++i) data1[i] *= tmp1;
    for (auto i = 0; i < data2.size(); ++i) data2[i] *= tmp2;
  } else if (mode[0] == 'a') {
    ;
  } else {
    assert(!"wrong parameter");
  }
//   for( auto i=0; i<data1.size(); ++i ){
//     std::cout << data1[i] << " ";
//   }
//   std::cout << "\n data2";
//   for( auto i=0; i<data2.size(); ++i ){
//     std::cout << data2[i] << " ";
//   }
//   std::cout << "\n";

  // calc l2 norm of values
  real err = 0.0;
  for (auto i = 0; i < data1.size(); ++i) {
    real tmp = std::abs(data1[i] - data2[i]);

    if (std::abs(tmp) / std::abs(data1[i]) > 1e-12)
      std::cout << " i: " << i << " values: " << std::abs(data2[i]) << " " << std::abs(data1[i]) << " ";
    err += tmp * tmp;
  }
  std::cout << "\n";
  err = std::sqrt(err);
  // open file in append mode
  std::ofstream ofs(filenameError.c_str(), std::ofstream::app);

  // write prefix and err
  ofs << prefix << " " << err << std::endl;

  std::cout << "error " << err << std::endl;

  return 0;
}

real l2Norm(std::vector<CombiDataType>& data) {
  real nrm = 0.0;
  for (auto d : data) nrm += std::abs(d) * std::abs(d);

  return sqrt(nrm);
}

void readPlotFile(const char* pltFileName, std::vector<CombiDataType>& data,
                  IndexVector& resolution) {
  // load gene grid from cp file
  std::cout << "reading plot file " << pltFileName << std::endl;

  //   double tstart = timing();

  // check if file exists
  struct stat buffer;
  assert(stat(pltFileName, &buffer) == 0);

  std::ifstream pltFile(pltFileName);

  // read dim and resolution
  int dim;
  pltFile.read((char*)&dim, sizeof(int));

  std::vector<int> res(dim);
  for (size_t i = 0; i < dim; ++i) pltFile.read((char*)&res[i], sizeof(int));

  resolution.assign(res.begin(), res.end());
  std::cout << "resolution " << resolution << std::endl;

  // calc data size
  size_t dsize = 1;
  for (auto r : resolution) dsize *= r;

  // read data
  std::vector<CombiDataType> tmp(dsize);
  pltFile.read((char*)&tmp[0], sizeof(CombiDataType) * dsize);

  pltFile.close();

  //   std::cout << "time to load plot file: " << timing() - tstart << "s" << std::endl;

  // create multiarray view on tmp
  IndexVector shape(resolution.begin(), resolution.end());
  if (dim == 6) {
    auto grid = createMultiArrayRef<CombiDataType, 6>(&tmp[0], shape);

    // copy tmp to data
    // copy data from local checkpoint to dfg
    // note that on the last process in some dimensions dfg is larger than the
    // local checkpoint
    for (size_t n = 0; n < shape[0]; ++n) {
      for (size_t m = 0; m < shape[1]; ++m) {
        for (size_t l = 0; l < shape[2]; ++l) {
          for (size_t k = 0; k < shape[3]; ++k) {
            for (size_t j = 0; j < shape[4]; ++j) {
              for (size_t i = 0; i < shape[5]; ++i) {
                data.push_back(grid[n][m][l][k][j][i]);
              }
            }
          }
        }
      }
    }
  } else if (dim == 5) {
    auto grid = createMultiArrayRef<CombiDataType, 5>(&tmp[0], shape);
    for (size_t n = 0; n < shape[0]; ++n) {
      for (size_t m = 0; m < shape[1]; ++m) {
        for (size_t l = 0; l < shape[2]; ++l) {
          for (size_t k = 0; k < shape[3]; ++k) {
            for (size_t j = 0; j < shape[4]; ++j) {
              data.push_back(grid[n][m][l][k][j]);
            }
          }
        }
      }
    }
  } else if (dim == 4) {
    auto grid = createMultiArrayRef<CombiDataType, 4>(&tmp[0], shape);
    for (size_t n = 0; n < shape[0]; ++n) {
      for (size_t m = 0; m < shape[1]; ++m) {
        for (size_t l = 0; l < shape[2]; ++l) {
          for (size_t k = 0; k < shape[3]; ++k) {
            data.push_back(grid[n][m][l][k]);
          }
        }
      }
    }
  } else if (dim == 3) {
    auto grid = createMultiArrayRef<CombiDataType, 3>(&tmp[0], shape);
    for (size_t n = 0; n < shape[0]; ++n) {
      for (size_t m = 0; m < shape[1]; ++m) {
        for (size_t l = 0; l < shape[2]; ++l) {
          data.push_back(grid[n][m][l]);
        }
      }
    }
  } else if (dim == 2) {
    auto grid = createMultiArrayRef<CombiDataType, 2>(&tmp[0], shape);
    for (size_t n = 0; n < shape[0]; ++n) {
      for (size_t m = 0; m < shape[1]; ++m) {
        data.push_back(grid[n][m]);
      }
    }
  } else {
    assert(!"wrong number of dimensions");
  }
}

void calcNorms(std::vector<CombiDataType>& dleft, std::vector<CombiDataType>& dright,
               std::vector<real>& Norms) {
  Norms.resize(9);

  real el1(0.0), el2(0.0), emax(0.0);
  real ll1(0.0), ll2(0.0), lmax(0.0);
  real rl1(0.0), rl2(0.0), rmax(0.0);

  assert(dleft.size() == dright.size());

  for (size_t i = 0; i < dleft.size(); ++i) {
    const CombiDataType ei = dleft[i] - dright[i];
    const real eiabs = std::abs(ei);
    const real liabs = std::abs(dleft[i]);
    const real riabs = std::abs(dright[i]);

    el1 += eiabs;
    el2 += eiabs * eiabs;
    if (eiabs > emax) emax = eiabs;

    ll1 += liabs;
    ll2 += liabs * liabs;
    if (liabs > lmax) lmax = liabs;

    rl1 += riabs;
    rl2 += riabs * riabs;
    if (riabs > rmax) rmax = riabs;
  }

  el2 = std::sqrt(el2);
  ll2 = std::sqrt(ll2);
  rl2 = std::sqrt(rl2);

  Norms[0] = el1;
  Norms[1] = el2;
  Norms[2] = emax;
  Norms[3] = ll1;
  Norms[4] = ll2;
  Norms[5] = lmax;
  Norms[6] = rl1;
  Norms[7] = rl2;
  Norms[8] = rmax;
}
