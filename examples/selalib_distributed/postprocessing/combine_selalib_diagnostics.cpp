#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/serialization/export.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/tokenizer.hpp>
#include <chrono>
#include <filesystem>
#include <iomanip>  // std::setprecision()
#include <iostream>
#include <string>
#include <vector>

#include "combischeme/CombiMinMaxScheme.hpp"
#include "utils/Types.hpp"

using namespace combigrid;
namespace qi = boost::spirit::qi;  // for parsing

int main(int argc, char** argv) {
  std::chrono::high_resolution_clock::time_point init_time =
      std::chrono::high_resolution_clock::now();

  // read in parameter file
  std::string paramfile = "ctparam";
  if (argc > 1) paramfile = argv[1];
  boost::property_tree::ptree cfg;
  boost::property_tree::ini_parser::read_ini(paramfile, cfg);

  DimType dim = cfg.get<DimType>("ct.dim");
  LevelVector lmin(dim), lmax(dim);
  cfg.get<std::string>("ct.lmin") >> lmin;  // minimal level vector for each grid
  cfg.get<std::string>("ct.lmax") >> lmax;  // maximum level vector -> level vector of target grid
  std::string basename = cfg.get<std::string>("preproc.basename");
  std::string nameDiagnostics = cfg.get<std::string>("application.name_diagnostics", "vp_B2_3d3v");

  CombiMinMaxScheme combischeme(dim, lmin, lmax);
  combischeme.createClassicalCombischeme();
  auto levels = combischeme.getCombiSpaces();
  auto coeffs = combischeme.getCoeffs();

  std::vector<std::string> diagnosticsFiles{nameDiagnostics +
                                            ".dat"};  // or nameDiagnostics + "add.dat"};

  std::string suffix = "";
  std::string baseFolder = "./" + basename;
  // read in all diagnostics files
  for (const auto& df : diagnosticsFiles) {
    auto combinedValues = std::vector<std::vector<double>>();
    auto kahanValues = std::vector<std::vector<double>>();
    for (size_t i = 0; i < levels.size(); i++) {
      // path to task folder
      std::string taskFolder = baseFolder + suffix + std::to_string(i);
      // assert that diagnostics are there
      std::string fullPathString = taskFolder + "/" + df;
      if (!std::filesystem::exists(fullPathString)) {
        throw std::runtime_error("Not found: " + fullPathString);
      }
      std::ifstream inputFileStream(fullPathString, std::ifstream::in);
      std::string line;
      size_t numLines = 0;
      while (std::getline(inputFileStream, line)) {
        ++numLines;
        std::istringstream iss(line);
        std::string token;
        // std::cout << "line" << line << std::endl;
        int column = 0;
        while (std::getline(iss, token, ' ')) {
          // std::cout << token << std::endl;
          double parsedDouble = 0.;
          try {
            parsedDouble = std::stod(token);
          } catch (std::invalid_argument& e) {
            continue;
          }
          if (combinedValues.size() < column + 1) {
            combinedValues.emplace_back();
            kahanValues.emplace_back();
          }
          if (combinedValues[column].size() < numLines) {
            combinedValues[column].push_back(parsedDouble * coeffs[i]);
            kahanValues[column].push_back(0.);
          } else {
            auto summand = parsedDouble * coeffs[i];
            auto y = summand - kahanValues[column][numLines - 1];
            auto t = combinedValues[column][numLines - 1] + y;
            kahanValues[column][numLines - 1] = (t - combinedValues[column][numLines - 1]) - y;
            combinedValues[column][numLines - 1] = t;
          }
          ++column;
        }
      }
    }
    std::ofstream outputFileStream("combined_" + df);
    outputFileStream << std::fixed << std::setprecision(12) << std::scientific;
    // write out combined data
    for (size_t j = 0; j < combinedValues[0].size(); ++j) {
      for (size_t k = 0; k < combinedValues.size(); ++k) {
        outputFileStream << " " << combinedValues[k][j];
      }
      outputFileStream << std::endl;
    }
  }

  std::chrono::high_resolution_clock::time_point end_time =
      std::chrono::high_resolution_clock::now();
  std::cout << "total "
            << std::chrono::duration_cast<std::chrono::microseconds>(end_time - init_time).count()
            << " " << std::endl;

  return 0;
}
