#include <boost/filesystem.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/multi_array.hpp>
#include <filesystem>
#include <iomanip>  // std::setprecision()
#include <iostream>
#include <string>
#include <vector>

#include <highfive/H5File.hpp>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "utils/Types.hpp"

using namespace combigrid;
namespace fs = boost::filesystem;  // for file operations

int main(int argc, char** argv) {
  [[maybe_unused]] auto mpiOnOff = MpiOnOff(&argc, &argv);

  std::chrono::high_resolution_clock::time_point init_time =
      std::chrono::high_resolution_clock::now();

  // read in parameter file
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

  DimType dim = cfg.get<DimType>("ct.dim");
  DimType dim_x = dim/2;
  assert(dim_x == 3);
  LevelVector lmin(dim), lmax(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;  // minimal level vector for each grid
  cfg.get<std::string>("ct.lmax") >> lmax;  // maximum level vector -> level vector of target grid
  // for some reason, the coordinates are reversed in selalib's h5 files
  LevelVector lmax_x = {lmax[2], lmax[1], lmax[0]};
  std::string basename = cfg.get<std::string>("preproc.basename");
  std::string nameDiagnostics = cfg.get<std::string>("application.name_diagnostics_phi", "diagnostics3d_add-");

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createClassicalCombischeme();
  // combischeme.makeFaultTolerant();
  auto levels = combischeme.getCombiSpaces();
  auto coeffs = combischeme.getCoeffs();
  std::cout << " combi levels " << levels << std::endl;

  std::string firstTaskFolder = basename + std::to_string(0);
  std::vector<std::string> diagnosticsFiles;
  // iterate files in folder starting with nameDiagnostics
  for (const auto& entry : fs::directory_iterator(firstTaskFolder)) {
    if (entry.path().filename().string().find(nameDiagnostics) == 0) {
      diagnosticsFiles.push_back(entry.path().filename().string());
    }
  }

  MPI_Comm commSelfPeriodic;
  MPI_Cart_create(MPI_COMM_SELF, 3, new int[3]{1,1,1}, new int[3]{1,1,1}, 1, &commSelfPeriodic);

  for (const auto& df : diagnosticsFiles) {
    auto combinedValues = combigrid::OwningDistributedFullGrid<double>(dim_x, lmax_x, commSelfPeriodic, {1,1,1}, {1,1,1}, false );
    auto kahanValues = combigrid::OwningDistributedFullGrid<double>(dim_x, lmax_x, commSelfPeriodic, {1,1,1}, {1,1,1}, false );
    for (size_t i = 0; i < levels.size(); i++) {
      IndexVector level_reversed {levels[i][2], levels[i][1], levels[i][0]};
      // path to task folder
      std::string taskFolder = basename + std::to_string(i);
      // assert that diagnostics are there
      std::string fullPathString = taskFolder + "/" + df;
      if (!fs::exists(fullPathString)) {
        throw std::runtime_error("Not found: " + fullPathString);
      }
      std::cout << " processing " << fullPathString << std::endl;
      std::ifstream inputFileStream(fullPathString, std::ifstream::in);
      // read values as h5 with highfive
      IndexVector extents {1<<level_reversed[0], 1<<level_reversed[1], 1<<level_reversed[2]};
      std::cout << "level " << level_reversed[0] << " " << level_reversed[1] << " " << level_reversed[2] << std::endl;
      boost::multi_array<double, 3> values{ extents, boost::fortran_storage_order()};
      HighFive::File h5_file(fullPathString, HighFive::File::ReadOnly);
      // assume one group in file
      HighFive::Group group = h5_file.getGroup("/");
      auto dataset = group.getDataSet("phi");
      dataset.read(values);
      if (!(values.storage_order() == boost::fortran_storage_order())) {
        throw std::runtime_error("storage order not fortran");
      }

      // pass values to DistributedFullGrid
      auto componentGrid = combigrid::DistributedFullGrid<double>(dim_x, level_reversed, commSelfPeriodic, {1,1,1}, values.data(), {1,1,1}, false);
      if (i == 0) {
        std::cout << "corner values " << fullPathString << std::endl;
        std::cout << values.shape()[0] << " " << values.shape()[1] << " " << values.shape()[2] << " " << std::endl;
        std::cout << values[0][0][0] << " " << values[0][0][extents[2]-1] << " " << values[0][extents[1]-1][0] << " " << values[0][extents[1]-1][extents[2]-1] << std::endl;
        std::cout << values[extents[0]-1][0][0] << " " << values[extents[0]-1][0][extents[2]-1] << " " << values[extents[0]-1][extents[1]-1][0] << " " << values[extents[0]-1][extents[1]-1][extents[2]-1] << std::endl;
      }
      {
        assert(values.shape()[0] == extents[0]);
        assert(values.shape()[1] == extents[1]);
        assert(values.shape()[2] == extents[2]);
        assert(values[0][0][0] == componentGrid.getData()[0]);
        assert(values[1][0][0] == componentGrid.getData()[1]);
        assert(values[2][0][0] == componentGrid.getData()[2]);
        assert(values[3][0][0] == componentGrid.getData()[3]);
        assert(values[extents[0]-1][0][0] == componentGrid.getData()[extents[0]-1]);
        assert(level_reversed[0] <= lmax_x[0]);
        assert(level_reversed[1] <= lmax_x[1]);
        assert(level_reversed[2] <= lmax_x[2]);
      }

      // interpolate on all points of combinedValues
      std::vector<real> coordinates_j = {0., 0., 0.};
      size_t j = 0;
      //three-fold loop, increase coordinate by grid spacing
      for (coordinates_j[2] = 0.; coordinates_j[2] < 1.; coordinates_j[2] += 1. / (1 << lmax_x[2])) {
        for (coordinates_j[1] = 0.; coordinates_j[1] < 1.; coordinates_j[1] += 1. / (1 << lmax_x[1])) {
          for (coordinates_j[0] = 0.; coordinates_j[0] < 1.; coordinates_j[0] += 1. / (1 << lmax_x[0])) {
            double summand = 0.;
            componentGrid.evalLocal(coordinates_j, summand);
            summand *= coeffs[i];
            auto y = summand - kahanValues.getData()[j];
            auto t = combinedValues.getData()[j] + y;
            kahanValues.getData()[j] = (t - combinedValues.getData()[j]) - y;
            combinedValues.getData()[j] = t;
            ++j;
          }
        }
      }
    }
    // output interpolated values to h5 file
    std::string outputFileName = "combined_phi_" + df;
    HighFive::File outputFile(outputFileName, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    HighFive::Group group = outputFile.createGroup("/");
    std::vector<size_t> extentsMax {static_cast<size_t>(1<<lmax_x[0]), static_cast<size_t>(1<<lmax_x[1]), static_cast<size_t>(1<<lmax_x[2])};
    HighFive::DataSet dataset = group.createDataSet<double>("phi", HighFive::DataSpace(extentsMax));
    dataset.write(combinedValues.getData());
    outputFile.flush();
  }

  std::chrono::high_resolution_clock::time_point end_time =
      std::chrono::high_resolution_clock::now();
  std::cout << "total "
            << std::chrono::duration_cast<std::chrono::seconds>(end_time - init_time).count()
            << " seconds" << std::endl;

  return 0;
}
