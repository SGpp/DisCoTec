#include <boost/multi_array.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <iostream>
#include <string>
#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "fullgrid/DistributedFullGrid.hpp"
#include "hierarchization/DistributedHierarchization.hpp"
#include "sparsegrid/DistributedSparseGridUniform.hpp"
#include "utils/Stats.hpp"
#include "utils/Types.hpp"

using namespace combigrid;
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
  DimType dim_x = dim / 2;
  assert(dim_x == 3);
  LevelVector lmin(dim), lmax(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;  // minimal level vector for each grid
  cfg.get<std::string>("ct.lmax") >> lmax;  // maximum level vector -> level vector of target grid
  // for some reason, the coordinates are reversed in selalib's h5 files
  LevelVector lmax_x = {lmax[2], lmax[1], lmax[0]};
  LevelVector lmin_x = {lmin[2], lmin[1], lmin[0]};
  std::string basename = cfg.get<std::string>("preproc.basename");
  std::string nameDiagnostics =
      cfg.get<std::string>("application.name_diagnostics_phi", "diagnostics3d_add-");

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createClassicalCombischeme();
  // combischeme.makeFaultTolerant();
  auto levels = combischeme.getCombiSpaces();
  auto coeffs = combischeme.getCoeffs();
  std::cout << " combi levels " << levels << std::endl;

  std::string firstTaskFolder = basename + std::to_string(0);
  std::vector<std::string> diagnosticsFiles;
  // iterate files in folder starting with nameDiagnostics
  for (const auto& entry : std::filesystem::directory_iterator(firstTaskFolder)) {
    if (entry.path().filename().string().find(nameDiagnostics) == 0) {
      diagnosticsFiles.push_back(entry.path().filename().string());
    }
  }

  MPI_Comm commSelfPeriodic;
  MPI_Cart_create(MPI_COMM_SELF, 3, new int[3]{1, 1, 1}, new int[3]{1, 1, 1}, 1, &commSelfPeriodic);

  for (const auto& df : diagnosticsFiles) {
    combigrid::OwningDistributedFullGrid<double> combinedFullGrid(dim_x, lmax_x, commSelfPeriodic,
                                                                  {1, 1, 1}, {1, 1, 1}, false);
    combigrid::DistributedSparseGridUniform<double> combinedSparseGrid(dim_x, lmax_x, lmin_x,
                                                                       commSelfPeriodic);
    combinedSparseGrid.registerDistributedFullGrid(combinedFullGrid);
    combinedSparseGrid.createSubspaceData();

    for (size_t i = 0; i < levels.size(); i++) {
      IndexVector level_reversed{levels[i][2], levels[i][1], levels[i][0]};
      // path to task folder
      std::string taskFolder = basename + std::to_string(i);
      // assert that diagnostics are there
      std::string fullPathString = taskFolder + "/" + df;
      if (!std::filesystem::exists(fullPathString)) {
        throw std::runtime_error("Not found: " + fullPathString);
      }
      std::cout << " processing " << fullPathString << std::endl;

      Stats::startEvent("read file");
      std::ifstream inputFileStream(fullPathString, std::ifstream::in);
      // read values as h5 with highfive
      IndexVector extents{1 << level_reversed[0], 1 << level_reversed[1], 1 << level_reversed[2]};
      std::cout << "level " << level_reversed[0] << " " << level_reversed[1] << " "
                << level_reversed[2] << std::endl;
      boost::multi_array<double, 3> values{extents, boost::fortran_storage_order()};
      HighFive::File h5_file(fullPathString, HighFive::File::ReadOnly);
      // assume one group in file
      HighFive::Group group = h5_file.getGroup("/");
      auto dataset = group.getDataSet("phi");
      dataset.read(values);
      if (!(values.storage_order() == boost::fortran_storage_order())) {
        throw std::runtime_error("storage order not fortran");
      }
      Stats::stopEvent("read file");
      auto durationRead = Stats::getDuration("read file");
      std::cout << "read file in " << durationRead << " milliseconds" << std::endl;

      // pass values to DistributedFullGrid
      auto componentGrid = combigrid::DistributedFullGrid<double>(
          dim_x, level_reversed, commSelfPeriodic, {1, 1, 1}, values.data(), {1, 1, 1}, false);
      if (i == 0) {
        std::cout << "corner values " << fullPathString << std::endl;
        std::cout << values.shape()[0] << " " << values.shape()[1] << " " << values.shape()[2]
                  << " " << std::endl;
        std::cout << values[0][0][0] << " " << values[0][0][extents[2] - 1] << " "
                  << values[0][extents[1] - 1][0] << " "
                  << values[0][extents[1] - 1][extents[2] - 1] << std::endl;
        std::cout << values[extents[0] - 1][0][0] << " "
                  << values[extents[0] - 1][0][extents[2] - 1] << " "
                  << values[extents[0] - 1][extents[1] - 1][0] << " "
                  << values[extents[0] - 1][extents[1] - 1][extents[2] - 1] << std::endl;
      }
      {
        assert(values.shape()[0] == extents[0]);
        assert(values.shape()[1] == extents[1]);
        assert(values.shape()[2] == extents[2]);
        assert(values[0][0][0] == componentGrid.getData()[0]);
        assert(values[1][0][0] == componentGrid.getData()[1]);
        assert(values[2][0][0] == componentGrid.getData()[2]);
        assert(values[3][0][0] == componentGrid.getData()[3]);
        assert(values[extents[0] - 1][0][0] == componentGrid.getData()[extents[0] - 1]);
        assert(level_reversed[0] <= lmax_x[0]);
        assert(level_reversed[1] <= lmax_x[1]);
        assert(level_reversed[2] <= lmax_x[2]);
      }

      Stats::startEvent("hierarchize");
      DistributedHierarchization::hierarchizeHierarchicalHat<double>(componentGrid,
                                                                     {true, true, true}, lmin_x);
      Stats::stopEvent("hierarchize");
      auto durationHierarchize = Stats::getDuration("hierarchize");
      std::cout << "hierarchized in " << durationHierarchize << " milliseconds" << std::endl;

      Stats::startEvent("kahan sum");
      combinedSparseGrid.addDistributedFullGrid(componentGrid, coeffs[i]);
      Stats::stopEvent("kahan sum");
      auto durationSum = Stats::getDuration("kahan sum");
      std::cout << "summed in " << durationSum << " milliseconds" << std::endl;
    }
    // scatter values to target grid
    Stats::startEvent("scatter");
    combinedFullGrid.extractFromUniformSG(combinedSparseGrid);
    Stats::stopEvent("scatter");
    Stats::startEvent("dehierarchize");
    DistributedHierarchization::dehierarchizeHierarchicalHat<double>(combinedFullGrid,
                                                                     {true, true, true}, lmin_x);
    Stats::stopEvent("dehierarchize");
    std::cout << "scattered in " << Stats::getDuration("scatter") << " and dehierarchized in "
              << Stats::getDuration("dehierarchize") << " milliseconds" << std::endl;

    boost::const_multi_array_ref<double, 3> combinedValues =
        boost::const_multi_array_ref<double, 3>(combinedFullGrid.getData(),
                                                combinedFullGrid.getGlobalSizes());
    boost::multi_array<double, 3> combinedValuesCopy = combinedValues;

    // output interpolated values to h5 file
    std::string outputFileName = "combined_phi_" + df;
    HighFive::File outputFile(outputFileName, HighFive::File::ReadWrite |
                                                  HighFive::File::OpenOrCreate |
                                                  HighFive::File::Truncate);
    HighFive::Group group;
    if (outputFile.exist("/")) {
      group = outputFile.getGroup("/");
    } else {
      group = outputFile.createGroup("/");
    }
    HighFive::DataSet dataset =
        group.createDataSet<double>("phi", HighFive::DataSpace::From(combinedValuesCopy));
    dataset.write(combinedValuesCopy);
    outputFile.flush();
  }

  std::chrono::high_resolution_clock::time_point end_time =
      std::chrono::high_resolution_clock::now();
  std::cout << "total "
            << std::chrono::duration_cast<std::chrono::seconds>(end_time - init_time).count()
            << " seconds" << std::endl;

  return 0;
}
