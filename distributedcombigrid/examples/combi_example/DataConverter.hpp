/*
this function reads in the data from a json file formated like that one from hyperdeal and returns the data 
we need for a Combi

*/

#ifndef DATACONVERTER_H
#define DATACONVERTER_H

#include <string>
#include <vector>
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"

#include "sgpp/distributedcombigrid/utils/Types.hpp"

using namespace combigrid;

class Converter{
    public:
    Converter(std::string filename):_fileName(filename){};
    Converter();
    void print();
    //reads the data from the file and stores it in the given variables.
    void readFile(DimType dim,LevelVector lmin, LevelVector lmax, LevelVector leval, IndexVector p,size_t ncombi, std::vector<bool> boundary);

    private:
    const std::string _fileName;
};
#endif