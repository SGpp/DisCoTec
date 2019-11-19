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
    
    Converter();
    //reads the data from a json file and stores it in param file
    
    void toParam(std::string, std::string);

    void toJSONForDealII(std::string, std::string, LevelVector);

    //
    void toJSON(std::string, std::string, LevelVector );

};
#endif