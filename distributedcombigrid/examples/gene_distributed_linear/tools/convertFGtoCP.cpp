/*
 * errorCalc.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: heenemo
 */
#include <assert.h>
#include <fstream>
#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "sgpp/distributedcombigrid/utils/LevelVector.hpp"
#include <iostream>
#include <vector>
#include "timing.h"
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"
#include "../src/GeneTask.hpp"

using namespace combigrid;

int main( int argc, char** argv ){
  assert( argc == 3 );

  const char* in = argv[1];
  const char* out = argv[2];

  std::cout << "loading FG " << in << std::endl;
  real tstart = timing();
  FullGrid<complex> fg( in );
  std::cout << "time to load FG " << timing() - tstart << "s" << std::endl;
  const LevelVector& fg_l = fg.getLevels();
  std::cout << "fg levels = "
	    << LevelVector( fg_l.begin(), fg_l.end() )
	    << std::endl;

  std::cout << "writing GENE checkpoint " << out << std::endl;

  GeneTask::saveCheckpoint( fg, out );

  return 0;
}
