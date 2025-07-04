#include <iostream>
#include <string>
#include <chrono>
#include <memory>
#include <thread>
#include <sys/stat.h>

#include "discotec/loadmodel/AveragingLoadModel.hpp"
#include "discotec/utils/LevelVector.hpp"

namespace combigrid {

  unsigned long int getAverageOfFirstColumn(std::istream& fs){
    std::string line, val;
    unsigned long int sum=0, count=0;
    while (std::getline (fs, line)) { 
      std::stringstream s (line);
      /* get first value (',' delimited) */
      getline (s, val, ',');       
      sum += std::stoul (val);
      ++count;
    }
    return sum/count;
  }
 
  void AveragingLoadModel::addLevelVectorToLoadModel(const LevelVector& lv) {
    durationsOfLevels_->insert(make_pair(lv, std::deque<DurationInformation>()));

    struct stat buf;
    // if file exists, read contents
    if (stat(getFilename(lv).c_str(), &buf) != -1){
      std::fstream fs;
      fs.open(getFilename(lv), std::fstream::in | std::fstream::app);
      fs.seekg(0, std::ios::beg);
      fs.clear();
      if (fs.peek() != EOF) {
        last_durations_avg_.insert(make_pair(lv, getAverageOfFirstColumn(fs)));
      }
      fs.clear();
      fs.close();
    }
    // else, make new empty file
    else{
      std::fstream fs;
      fs.open(getFilename(lv), std::fstream::out);    
      fs.close();
    }

  }
} /* namespace combigrid */

