#include <string>
#include <vector>
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "DataConverter.hpp"
#include "../../../tools/json.hpp"

using namespace combigrid;



    void Converter::readFile(DimType dim,LevelVector lmin, LevelVector lmax, LevelVector leval, IndexVector p,size_t ncombi, std::vector<bool> boundary){
        namespace pt = boost::property_tree;

        // Create a root
        pt::ptree root;

        // Load the json file in this ptree
        try{
            pt::read_json(_fileName, root);
            int dimx=root.get<int>("General.DimX");
            int dimv=root.get<int>("General.DimV");
            dim*=dimx+dimv;

            //lmax:
            int refineX=root.get<int>("Case.NRefinementsX");
            int refineY=root.get<int>("Case.NRefinementsY");
            int subX=root.get<int>("Case.NSubdivisionsX");
            int subV=root.get<int>("Case.NSubdivisionsV");

            //lmin:


            //ncombi?:

            //p for partitions
            int partX=root.get<int>("General.PartitionX");
            int partY=root.get<int>("General.PartitionY");

            //boundary
            bool periodicX=root.get<bool>("Case.PeriodicX",true);
            bool periodicV=root.get<bool>("Case.PeriodicV",true);

            leval =lmax;
        }
        catch(pt::json_parser::json_parser_error error){

        }
        catch(...){

        }
    }
    void Converter::print(){
        std::cout << "Name:" <<_fileName;
    }

