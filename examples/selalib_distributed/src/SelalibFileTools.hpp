#pragma once

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "../../../include/discotec/utils/Types.hpp"

// helper funtion to read a bool vector from string
inline std::vector<bool>& operator>>(const std::string& str, std::vector<bool>& vec) {
  std::vector<std::string> strs;
  boost::split(strs, str, boost::is_any_of(" "));

  assert(vec.size() == strs.size());

  for (size_t i = 0; i < strs.size(); ++i) vec[i] = boost::lexical_cast<bool>(strs[i]);

  return vec;
}

// cf https://stackoverflow.com/questions/37931691/replace-a-word-in-text-file-using-c
std::string getFile(std::ifstream& is) {
  std::stringstream buffer;
  buffer << is.rdbuf();
  return buffer.str();
}

// cf https://stackoverflow.com/questions/5878775/how-to-find-and-replace-string/5878802
std::string replaceFirstOccurrence(std::string& s, const std::string& toReplace,
                                   const std::string& replaceWith) {
  std::size_t pos = s.find(toReplace);
  if (pos == std::string::npos) return s;
  return s.replace(pos, toReplace.length(), replaceWith);
}

void setCheckpointRestart(const std::string& basename,
                          const std::vector<combigrid::LevelVector>& levels,
                          const std::string& suffix = "") {
  std::string baseFolder = "./" + basename;
  for (size_t i = 0; i < levels.size(); i++) {
    // path to task folder
    std::string taskFolder = baseFolder + suffix + std::to_string(i);
    // assert that checkpoint is there
    std::string checkpointString = taskFolder + "/distribution_end-0000.h5";
    if (!std::filesystem::exists(checkpointString)) {
      throw std::runtime_error("No checkpoint to re-start from " + checkpointString);
    }
    {
      // adapt each parameter file
      std::ifstream inputFileStream(taskFolder + "/param.nml", std::ifstream::in);
      auto contents = getFile(inputFileStream);
      std::string newRestartInfo =
          "restart = .true.\n        restart_filename = \"distribution_end-0000.h5\"";
      contents = replaceFirstOccurrence(contents, "restart = .false.", newRestartInfo);
      std::ofstream outputFileStream(taskFolder + "/param_new.nml");
      outputFileStream << contents;
    }
    std::rename((taskFolder + "/param_new.nml").c_str(), (taskFolder + "/param.nml").c_str());
  }
}

bool adaptParameterFile(const std::string& infilename, const std::string& outfilename,
                        const std::vector<int>& resolution, const std::vector<int>& p,
                        size_t nsteps, double dt, size_t n_diagnostics,
                        const std::string& name_diagnostics) {
  assert(resolution.size() == p.size());
  std::ifstream inputFileStream(infilename, std::ifstream::in);
  auto contents = getFile(inputFileStream);
  contents = replaceFirstOccurrence(contents, "$nsteps", std::to_string(nsteps));
  contents = replaceFirstOccurrence(contents, "$dt", std::to_string(dt));
  contents = replaceFirstOccurrence(contents, "$n_diagnostics", std::to_string(n_diagnostics));
  contents = replaceFirstOccurrence(contents, "$name_diagnostics", name_diagnostics);
  for (combigrid::DimType d = 0; d < resolution.size(); ++d) {
    contents = replaceFirstOccurrence(contents, "$nx" + std::to_string(d + 1),
                                      std::to_string(resolution[d]));
    contents = replaceFirstOccurrence(contents, "$p" + std::to_string(d + 1), std::to_string(p[d]));
  }
  std::ofstream outputFileStream(outfilename);
  outputFileStream << contents;
  return true;
}

bool adaptParameterFileFirstFolder(const std::string& basename, const std::vector<int>& resolution,
                                   const std::vector<int>& p, size_t nsteps, double dt,
                                   size_t n_diagnostics, const std::string& name_diagnostics,
                                   const std::string& suffix = "") {
  std::string baseFolder = "./" + basename;
  std::string taskFolder = baseFolder + suffix + std::to_string(0);
  std::string templateFolder = "./template";

  bool yes = adaptParameterFile(templateFolder + "/param.nml", taskFolder + "/param.nml",
                                resolution, p, nsteps, dt, n_diagnostics, name_diagnostics);
  assert(yes);
  return yes;
}

bool createTaskFolders(const std::string& basename,
                       const std::vector<combigrid::LevelVector>& levels,
                       const std::vector<size_t>& taskNumbers, const std::vector<int>& p,
                       size_t nsteps, double dt, size_t n_diagnostics,
                       const std::string& name_diagnostics, const std::string& suffix = "") {
  assert(levels.size() == taskNumbers.size());
  std::string baseFolder = "./" + basename;
  std::string templateFolder = "./template";
  bool yes = std::filesystem::exists(templateFolder);
  assert(yes);
  for (size_t i = 0; i < levels.size(); i++) {
    // path to task folder
    std::string taskFolder = baseFolder + suffix + std::to_string(taskNumbers[i]);
    // copy template directory
    if (!std::filesystem::exists(taskFolder) && !std::filesystem::create_directory(taskFolder)) {
      throw std::runtime_error("Cannot create destination directory " + taskFolder);
    }
    // adapt each parameter file
    std::vector<int> resolutions;
    for (combigrid::DimType d = 0; d < levels[0].size(); ++d) {
      resolutions.push_back(static_cast<combigrid::IndexType>(std::pow(2, levels[i][d])));
    }
    yes = adaptParameterFile(templateFolder + "/param.nml", taskFolder + "/param.nml", resolutions,
                             p, nsteps, dt, n_diagnostics, name_diagnostics);
    assert(yes);
    // copy all other files
    for (const auto& dirEnt : std::filesystem::recursive_directory_iterator{templateFolder}) {
      const auto& path = dirEnt.path();
      const auto relativePathStr = std::filesystem::relative(path, templateFolder);

      // auto relativePathStr = path.string();
      // boost::replace_first(relativePathStr, templateFolder, "");
      std::filesystem::copy(path, taskFolder / relativePathStr,
                            std::filesystem::copy_options::skip_existing);
    }
  }
  return yes;
}
