#include <string>
#include <vector>
#include <math.h>
#include <sstream>
#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"

#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"
#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "sgpp/distributedcombigrid/utils/Types.hpp"
#include "DataConverter.hpp"

using namespace combigrid;


    Converter::Converter(){};
    void Converter::toParam(std::string _filename, std::string newname){
        namespace pt = boost::property_tree;

        // Create a root
        pt::ptree root,ini_root,ini_ct, ini_application, ini_manager;

        // Load the json file in this ptree
        try{
            pt::read_json(_filename, root);
            int dimx=root.get<int>("General.DimX");
            int dimv=root.get<int>("General.DimV");
            int dim=dimx+dimv;
            //lmax:
            int refineX=root.get<int>("Case.NRefinementsX");
            int refineV=root.get<int>("Case.NRefinementsV");
            int subX=root.get<int>("Case.NSubdivisionsX.X");
            int subV=root.get<int>("Case.NSubdivisionsV.X");
            LevelVector lma;

            //zur Basis 2 hoch beziehen!
            for(int i=0; i<dimx;i++){
                lma.push_back(pow(subX*(refineX+1),1.0/3.0)); 
           }
            for(int i=0; i<dimx;i++){
                lma.push_back(pow(1.0/3.0,subV*(refineV+1)));
            }
            
            //lmin=lmax


            //ncombi=10

            //p for partitions
            int partX=root.get<int>("General.PartitionX");
            int partY=root.get<int>("General.PartitionV");
            LevelVector p_{partX,partY};
            
            //boundary
            
            LevelVector bound;
            for(int i=0;i<dimx;i++){
                bound.push_back(true);
            }
            for(int i=0;i<dimv;i++){
                bound.push_back(true);
            }
            std::string levelString;
            std::ostringstream vec,p;
            vec << lma;
            levelString=vec.str();
            p<<p_[0]<<" "<<p_[1]<<std::endl;
            //leval=lmax;
            ini_ct.put("dim",dim);
            ini_ct.put("lmin",levelString);
            ini_ct.put("lmax",levelString);
            ini_ct.put("leval",levelString);
            ini_ct.put("p",p.str());
            ini_ct.put("ncombi",10);

            ini_root.add_child("ct",ini_ct);
            ini_application.put("dt",0.001);
            ini_application.put("steps",100);
            ini_root.add_child("application",ini_application);
            ini_manager.put("ngroup",2);
            ini_manager.put("nprocs",1);
            ini_root.add_child("manager",ini_manager);
            pt::write_ini(newname,ini_root);

        }
        catch(const pt::json_parser::json_parser_error& error){
            std::cout<< "Not a valid json File";
        }
        catch(...){
            std::cout << "Something went wrong";
        }
    }
    
    //converts the data from the ini file to a json format and stores the data in newname.json

    void Converter::toJSON(std::string _fileName, std::string newname,LevelVector l){
        namespace pt = boost::property_tree;

        // Create a root
        pt::ptree root;
        pt::ptree json_root, json_general, json_case, json_SpatialDiscretization, json_TemporalDiscretization, json_postprocessing;
        pt::ptree jsubX, jsubV,jsout, jVTK;
        
        // Load the ini file in this ptree
        try{
            pt::read_ini(_fileName, root);
            std::cout << "Name "<<_fileName;
            int dim=root.get<int>("ct.dim",2);
            int dimx=root.get<int>("ct.dimx",dim/2);
            int dimv=root.get<int>("ct.dimv",dim/2);

            //degree is 1 per default
            int degx=root.get<int>("ct.degreex",1);
            int degv=root.get<int>("ct.degreev",1);

            // case is per default hyperrectangle and periodic is per default 
            //nrefinments is set to 0 and the distribution is controlled via subdivisions

            //lmax:
            LevelVector lmax(dimx+dimv),p(dimx+dimv);
            lmax=l;
            root.get<std::string>("ct.p") >> p;
            
            json_general.put("DimX", dimx);
            json_general.put("DimV", dimv);
            json_general.put("DegreeX", degx);
            json_general.put("DegreeV", degv);
            json_general.put("PartitionX", p[0]);
            json_general.put("PartitionV", p[1]);
            json_general.put("Case","hyperrectangle");
            json_root.add_child("General", json_general);

            json_case.put("NRefinementsX",0);
            json_case.put("NRefinementsV",0);
            json_case.put("PeriodicX",true);
            json_case.put("PeriodicV",true);
            jsubX.put("X",pow(2,lmax[0]));
            jsubX.put("Y",pow(2,lmax[0]));
            jsubX.put("Z",pow(2,lmax[0]));
            jsubV.put("X",pow(2,lmax[1]));
            jsubV.put("Y",pow(2,lmax[1]));
            jsubV.put("Z",pow(2,lmax[1]));
            json_case.add_child("NSubdivisionsX",jsubX);
            json_case.add_child("NSubdivisionsV",jsubV);
            json_root.add_child("Case", json_case);

            json_SpatialDiscretization.put("TriangulationType","FullyDistributed");
            json_SpatialDiscretization.put("MappingX",1);
            json_SpatialDiscretization.put("MappingV",1);
            json_root.add_child("SpatialDiscretization",json_SpatialDiscretization);

            json_TemporalDiscretization.put("FinalTime",2);
            json_TemporalDiscretization.put("CFLNumber",0.1);
            json_root.add_child("TemporalDiscretization",json_TemporalDiscretization);
            
            jsout.put("Tick",0.1);
            json_postprocessing.add_child("StandardOutput",jsout);
            jVTK.put("Enabled",false);
            json_postprocessing.add_child("VTK",jVTK);
            
            json_root.add_child("Postprocessing", json_postprocessing);
            pt::write_json(newname,json_root);
           
            
        }
        catch(const pt::ini_parser::ini_parser_error& error){
            std::cout<< "Not a valid File";
        }
        catch(...){
            std::cout << "Something went wrong";
        }
    }

    void Converter::toJSONForDealII(std::string _fileName, std::string newname,LevelVector l, int output_prefix){
        namespace pt = boost::property_tree;

        // Create a root
        pt::ptree root;
        pt::ptree json_root, json_general, json_case, json_SpatialDiscretization, json_TemporalDiscretization, json_postprocessing;
        pt::ptree jsubX, jsubV,jsout, jVTK;
        
        // Load the ini file in this ptree
        try{
            pt::read_ini(_fileName, root);
            std::cout << "Name "<<_fileName;
            int dim=root.get<int>("ct.dim",2);

            //degree is 1 per default
            int deg=root.get<int>("ct.degree",1);


            double dt=root.get<double>("application.dt",1);
            std::cout << "ct.dt="<<dt<<std::endl;
            int nsteps=root.get<int>("ct.nsteps",1);
            std::string isDG=root.get<std::string>("ct.FE","FE_Q");

            // case is per default hyperrectangle and periodic is per default 
            //nrefinments is set to 0 and the distribution is controlled via subdivisions

            //lmax:
            LevelVector lmax(dim),p(dim);
            lmax=l;
            root.get<std::string>("ct.p") >> p;
            
            json_general.put("Dim", dim);
            json_general.put("Degree", deg);
            json_general.put("PartitionX", p[std::min(std::max(dim-1, 0),0)]);
            json_general.put("PartitionY", p[std::min(std::max(dim-1, 1),1)]);
            json_general.put("PartitionZ", p[std::min(dim-1, 2)]);
            json_general.put("Case","hyperrectangle");
            json_root.add_child("General", json_general);

            json_case.put("NRefinements",0);
            json_case.put("Periodic",true);
            jsubX.put("X",pow(2,lmax[std::min(std::max(dim-1, 0),0)]));
            jsubX.put("Y",pow(2,lmax[std::min(std::max(dim-1, 1),1)]));
            //watch out here
            
            jsubX.put("Z",pow(2,lmax[std::min(dim-1, 2)]));
            json_case.add_child("NSubdivisions",jsubX);
            json_root.add_child("Case", json_case);

            json_SpatialDiscretization.put("TriangulationType","FullyDistributed");
            json_SpatialDiscretization.put("Mapping",1);
            json_SpatialDiscretization.put("FE",isDG);
            json_root.add_child("SpatialDiscretization",json_SpatialDiscretization);

            json_TemporalDiscretization.put("FinalTime",1);
            json_TemporalDiscretization.put("CFLNumber",0.15);
            json_root.add_child("TemporalDiscretization",json_TemporalDiscretization);
            
            jsout.put("Tick",dt);
            json_postprocessing.add_child("StandardOutput",jsout);
            jVTK.put("Enabled",false);
            //jVTK.put("Prefix","solution/solution_"+std::to_string(output_prefix));
            json_postprocessing.add_child("VTK",jVTK);
            
            json_root.add_child("Postprocessing", json_postprocessing);
            pt::write_json(newname,json_root);
           
            
        }
        catch(const pt::ini_parser::ini_parser_error& error){
            std::cout<< "Not a valid File";
        }
        catch(...){
            std::cout << "Something went wrong";
        }
    }

