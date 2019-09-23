[33mcommit 88d37a7e6f5e302d8ef562cca0b3abb58a7beaa2[m
Author: Marius Goehring <goehrims@ipvs.uni-stuttgart.de>
Date:   Sun Sep 15 21:25:34 2019 +0200

    new Class DataConverter

[1mdiff --git a/compile.sh b/compile.sh[m
[1mindex 4ec902b..43146e3 100755[m
[1m--- a/compile.sh[m
[1m+++ b/compile.sh[m
[36m@@ -1,2 +1,2 @@[m
[31m-scons -j 8 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 VERBOSE=1 RUN_BOOST_TESTS=1 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=mpicc.mpich FC=mpifort.mpich CXX=mpicxx.mpich OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 #DEBUG_OUTPUT=1[m
[32m+[m[32mscons -j 8 SG_ALL=1 SG_DISTRIBUTEDCOMBIGRID=1 VERBOSE=1 RUN_BOOST_TESTS=0 RUN_CPPLINT=0 BUILD_STATICLIB=0 CC=mpicc.mpich FC=mpifort.mpich CXX=mpicxx.mpich OPT=1 TIMING=0 UNIFORMDECOMPOSITION=1 #DEBUG_OUTPUT=1[m
 [m
[1mdiff --git a/distributedcombigrid/examples/combi_example/DataConverter.cpp b/distributedcombigrid/examples/combi_example/DataConverter.cpp[m
[1mnew file mode 100644[m
[1mindex 0000000..dc7b7f9[m
[1m--- /dev/null[m
[1m+++ b/distributedcombigrid/examples/combi_example/DataConverter.cpp[m
[36m@@ -0,0 +1,62 @@[m
[32m+[m[32m#include <string>[m
[32m+[m[32m#include <vector>[m
[32m+[m[32m#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"[m
[32m+[m
[32m+[m[32m#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"[m
[32m+[m[32m#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"[m
[32m+[m[32m#include <boost/property_tree/ptree.hpp>[m
[32m+[m[32m#include <boost/property_tree/json_parser.hpp>[m
[32m+[m
[32m+[m[32m#include "sgpp/distributedcombigrid/utils/Types.hpp"[m
[32m+[m[32m#include "DataConverter.hpp"[m
[32m+[m[32m#include "../../../tools/json.hpp"[m
[32m+[m
[32m+[m[32musing namespace combigrid;[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m[32m    void Converter::readFile(DimType dim,LevelVector lmin, LevelVector lmax, LevelVector leval, IndexVector p,size_t ncombi, std::vector<bool> boundary){[m
[32m+[m[32m        namespace pt = boost::property_tree;[m
[32m+[m
[32m+[m[32m        // Create a root[m
[32m+[m[32m        pt::ptree root;[m
[32m+[m
[32m+[m[32m        // Load the json file in this ptree[m
[32m+[m[32m        try{[m
[32m+[m[32m            pt::read_json(_fileName, root);[m
[32m+[m[32m            int dimx=root.get<int>("General.DimX");[m
[32m+[m[32m            int dimv=root.get<int>("General.DimV");[m
[32m+[m[32m            dim*=dimx+dimv;[m
[32m+[m
[32m+[m[32m            //lmax:[m
[32m+[m[32m            int refineX=root.get<int>("Case.NRefinementsX");[m
[32m+[m[32m            int refineY=root.get<int>("Case.NRefinementsY");[m
[32m+[m[32m            int subX=root.get<int>("Case.NSubdivisionsX");[m
[32m+[m[32m            int subV=root.get<int>("Case.NSubdivisionsV");[m
[32m+[m
[32m+[m[32m            //lmin:[m
[32m+[m
[32m+[m
[32m+[m[32m            //ncombi?:[m
[32m+[m
[32m+[m[32m            //p for partitions[m
[32m+[m[32m            int partX=root.get<int>("General.PartitionX");[m
[32m+[m[32m            int partY=root.get<int>("General.PartitionY");[m
[32m+[m
[32m+[m[32m            //boundary[m
[32m+[m[32m            bool periodicX=root.get<bool>("Case.PeriodicX",true);[m
[32m+[m[32m            bool periodicV=root.get<bool>("Case.PeriodicV",true);[m
[32m+[m
[32m+[m[32m            leval =lmax;[m
[32m+[m[32m        }[m
[32m+[m[32m        catch(pt::json_parser::json_parser_error error){[m
[32m+[m
[32m+[m[32m        }[m
[32m+[m[32m        catch(...){[m
[32m+[m
[32m+[m[32m        }[m
[32m+[m[32m    }[m
[32m+[m[32m    void Converter::print(){[m
[32m+[m[32m        std::cout << "Name:" <<_fileName;[m
[32m+[m[32m    }[m
[32m+[m
[1mdiff --git a/distributedcombigrid/examples/combi_example/DataConverter.hpp b/distributedcombigrid/examples/combi_example/DataConverter.hpp[m
[1mnew file mode 100644[m
[1mindex 0000000..5918400[m
[1m--- /dev/null[m
[1m+++ b/distributedcombigrid/examples/combi_example/DataConverter.hpp[m
[36m@@ -0,0 +1,32 @@[m
[32m+[m[32m/*[m
[32m+[m[32mthis function reads in the data from a json file formated like that one from hyperdeal and returns the data[m[41m [m
[32m+[m[32mwe need for a Combi[m
[32m+[m
[32m+[m[32m*/[m
[32m+[m
[32m+[m[32m#ifndef DATACONVERTER_H[m
[32m+[m[32m#define DATACONVERTER_H[m
[32m+[m
[32m+[m[32m#include <string>[m
[32m+[m[32m#include <vector>[m
[32m+[m[32m#include "sgpp/distributedcombigrid/combischeme/CombiMinMaxScheme.hpp"[m
[32m+[m
[32m+[m[32m#include "sgpp/distributedcombigrid/loadmodel/LinearLoadModel.hpp"[m
[32m+[m[32m#include "sgpp/distributedcombigrid/manager/CombiParameters.hpp"[m
[32m+[m
[32m+[m[32m#include "sgpp/distributedcombigrid/utils/Types.hpp"[m
[32m+[m
[32m+[m[32musing namespace combigrid;[m
[32m+[m
[32m+[m[32mclass Converter{[m
[32m+[m[32m    public:[m
[32m+[m[32m    Converter(std::string filename):_fileName(filename){};[m
[32m+[m[32m    Converter();[m
[32m+[m[32m    void print();[m
[32m+[m[32m    //reads the data from the file and stores it in the given variables.[m
[32m+[m[32m    void readFile(DimType dim,LevelVector lmin, LevelVector lmax, LevelVector leval, IndexVector p,size_t ncombi, std::vector<bool> boundary);[m
[32m+[m
[32m+[m[32m    private:[m
[32m+[m[32m    const std::string _fileName;[m
[32m+[m[32m};[m
[32m+[m[32m#endif[m
\ No newline at end of file[m
[1mdiff --git a/distributedcombigrid/examples/combi_example/combi_example b/distributedcombigrid/examples/combi_example/combi_example[m
[1mnew file mode 100755[m
[1mindex 0000000..fcfd569[m
Binary files /dev/null and b/distributedcombigrid/examples/combi_example/combi_example differ
[1mdiff --git a/distributedcombigrid/examples/combi_example/combi_example.cpp b/distributedcombigrid/examples/combi_example/combi_example.cpp[m
[1mindex 06e7d3a..6f637b1 100644[m
[1m--- a/distributedcombigrid/examples/combi_example/combi_example.cpp[m
[1m+++ b/distributedcombigrid/examples/combi_example/combi_example.cpp[m
[36m@@ -26,6 +26,8 @@[m
 #include "sgpp/distributedcombigrid/utils/Types.hpp"[m
 // include user specific task. this is the interface to your application[m
 #include "TaskExample.hpp"[m
[32m+[m[32m#include "DataConverter.cpp"[m
[32m+[m[32m#include "DataConverter.hpp"[m
 [m
 using namespace combigrid;[m
 [m
[36m@@ -83,6 +85,8 @@[m [mint main(int argc, char** argv) {[m
     dt = cfg.get<combigrid::real>("application.dt");[m
     nsteps = cfg.get<size_t>("application.nsteps");[m
 [m
[32m+[m[32m    Converter converter("ct.param");[m
[32m+[m[32m    converter.print();[m
     // TODO: read in boundary vector from ctparam[m
     std::vector<bool> boundary(dim, true);[m
 [m
[1mdiff --git a/distributedcombigrid/examples/combi_example/ctparam b/distributedcombigrid/examples/combi_example/ctparam[m
[1mindex 29e5624..e3927c0 100644[m
[1m--- a/distributedcombigrid/examples/combi_example/ctparam[m
[1m+++ b/distributedcombigrid/examples/combi_example/ctparam[m
[36m@@ -4,7 +4,7 @@[m [mlmin = 3 3[m
 lmax = 10 10  [m
 leval = 5 5 [m
 p = 1 1[m
[31m-ncombi = 10 t[m
[32m+[m[32mncombi = 10[m
 [m
 [application][m
 dt = 1e-3[m
[1mdiff --git a/lib/__init__.py b/lib/__init__.py[m
[1mnew file mode 100644[m
[1mindex 0000000..e69de29[m
[1mdiff --git a/lib/pysgpp/__init__.py b/lib/pysgpp/__init__.py[m
[1mnew file mode 100644[m
[1mindex 0000000..65d4f85[m
[1m--- /dev/null[m
[1m+++ b/lib/pysgpp/__init__.py[m
[36m@@ -0,0 +1,2 @@[m
[32m+[m[32mfrom pysgpp_swig import *[m
[32m+[m[32mimport extensions[m
[1mdiff --git a/lib/pysgpp/extensions/__init__.py b/lib/pysgpp/extensions/__init__.py[m
[1mnew file mode 100644[m
[1mindex 0000000..e69de29[m
[1mdiff --git a/tools/json.hpp b/tools/json.hpp[m
[1mnew file mode 100644[m
[1mindex 0000000..2a32a82[m
[1m--- /dev/null[m
[1m+++ b/tools/json.hpp[m
[36m@@ -0,0 +1,22684 @@[m
[32m+[m[32m/*[m
[32m+[m[32m    __ _____ _____ _____[m
[32m+[m[32m __|  |   __|     |   | |  JSON for Modern C++[m
[32m+[m[32m|  |  |__   |  |  | | | |  version 3.7.0[m
[32m+[m[32m|_____|_____|_____|_|___|  https://github.com/nlohmann/json[m
[32m+[m
[32m+[m[32mLicensed under the MIT License <http://opensource.org/licenses/MIT>.[m
[32m+[m[32mSPDX-License-Identifier: MIT[m
[32m+[m[32mCopyright (c) 2013-2019 Niels Lohmann <http://nlohmann.me>.[m
[32m+[m
[32m+[m[32mPermission is hereby  granted, free of charge, to any  person obtaining a copy[m
[32m+[m[32mof this software and associated  documentation files (the "Software"), to deal[m
[32m+[m[32min the Software  without restriction, including without  limitation the rights[m
[32m+[m[32mto  use, copy,  modify, merge,  publish, distribute,  sublicense, and/or  sell[m
[32m+[m[32mcopies  of  the Software,  and  to  permit persons  to  whom  the Software  is[m
[32m+[m[32mfurnished to do so, subject to the following conditions:[m
[32m+[m
[32m+[m[32mThe above copyright notice and this permission notice shall be included in all[m
[32m+[m[32mcopies or substantial portions of the Software.[m
[32m+[m
[32m+[m[32mTHE SOFTWARE  IS PROVIDED "AS  IS", WITHOUT WARRANTY  OF ANY KIND,  EXPRESS OR[m
[32m+[m[32mIMPLIED,  INCLUDING BUT  NOT  LIMITED TO  THE  WARRANTIES OF  MERCHANTABILITY,[m
[32m+[m[32mFITNESS FOR  A PARTICULAR PURPOSE AND  NONINFRINGEMENT. IN NO EVENT  SHALL THE[m
[32m+[m[32mAUTHORS  OR COPYRIGHT  HOLDERS  BE  LIABLE FOR  ANY  CLAIM,  DAMAGES OR  OTHER[m
[32m+[m[32mLIABILITY, WHETHER IN AN ACTION OF  CONTRACT, TORT OR OTHERWISE, ARISING FROM,[m
[32m+[m[32mOUT OF OR IN CONNECTION WITH THE SOFTWARE  OR THE USE OR OTHER DEALINGS IN THE[m
[32m+[m[32mSOFTWARE.[m
[32m+[m[32m*/[m
[32m+[m
[32m+[m[32m#ifndef INCLUDE_NLOHMANN_JSON_HPP_[m
[32m+[m[32m#define INCLUDE_NLOHMANN_JSON_HPP_[m
[32m+[m
[32m+[m[32m#define NLOHMANN_JSON_VERSION_MAJOR 3[m
[32m+[m[32m#define NLOHMANN_JSON_VERSION_MINOR 7[m
[32m+[m[32m#define NLOHMANN_JSON_VERSION_PATCH 0[m
[32m+[m
[32m+[m[32m#include <algorithm> // all_of, find, for_each[m
[32m+[m[32m#include <cassert> // assert[m
[32m+[m[32m#include <ciso646> // and, not, or[m
[32m+[m[32m#include <cstddef> // nullptr_t, ptrdiff_t, size_t[m
[32m+[m[32m#include <functional> // hash, less[m
[32m+[m[32m#include <initializer_list> // initializer_list[m
[32m+[m[32m#include <iosfwd> // istream, ostream[m
[32m+[m[32m#include <iterator> // random_access_iterator_tag[m
[32m+[m[32m#include <memory> // unique_ptr[m
[32m+[m[32m#include <numeric> // accumulate[m
[32m+[m[32m#include <string> // string, stoi, to_string[m
[32m+[m[32m#include <utility> // declval, forward, move, pair, swap[m
[32m+[m[32m#include <vector> // vector[m
[32m+[m
[32m+[m[32m// #include <nlohmann/adl_serializer.hpp>[m
[32m+[m
[32m+[m
[32m+[m[32m#include <utility>[m
[32m+[m
[32m+[m[32m// #include <nlohmann/detail/conversions/from_json.hpp>[m
[32m+[m
[32m+[m
[32m+[m[32m#include <algorithm> // transform[m
[32m+[m[32m#include <array> // array[m
[32m+[m[32m#include <ciso646> // and, not[m
[32m+[m[32m#include <forward_list> // forward_list[m
[32m+[m[32m#include <iterator> // inserter, front_inserter, end[m
[32m+[m[32m#include <map> // map[m
[32m+[m[32m#include <string> // string[m
[32m+[m[32m#include <tuple> // tuple, make_tuple[m
[32m+[m[32m#include <type_traits> // is_arithmetic, is_same, is_enum, underlying_type, is_convertible[m
[32m+[m[32m#include <unordered_map> // unordered_map[m
[32m+[m[32m#include <utility> // pair, declval[m
[32m+[m[32m#include <valarray> // valarray[m
[32m+[m
[32m+[m[32m// #include <nlohmann/detail/exceptions.hpp>[m
[32m+[m
[32m+[m
[32m+[m[32m#include <exception> // exception[m
[32m+[m[32m#include <stdexcept> // runtime_error[m
[32m+[m[32m#include <string> // to_string[m
[32m+[m
[32m+[m[32m// #include <nlohmann/detail/input/position_t.hpp>[m
[32m+[m
[32m+[m
[32m+[m[32m#include <cstddef> // size_t[m
[32m+[m
[32m+[m[32mnamespace nlohmann[m
[32m+[m[32m{[m
[32m+[m[32mnamespace detail[m
[32m+[m[32m{[m
[32m+[m[32m/// struct to capture the start position of the current token[m
[32m+[m[32mstruct position_t[m
[32m+[m[32m{[m
[32m+[m[32m    /// the total number of characters read[m
[32m+[m[32m    std::size_t chars_read_total = 0;[m
[32m+[m[32m    /// the number of characters read in the current line[m
[32m+[m[32m    std::size_t chars_read_current_line = 0;[m
[32m+[m[32m    /// the number of lines read[m
[32m+[m[32m    std::size_t lines_read = 0;[m
[32m+[m
[32m+[m[32m    /// conversion to size_t to preserve SAX interface[m
[32m+[m[32m    constexpr operator size_t() const[m
[32m+[m[32m    {[m
[32m+[m[32m        return chars_read_total;[m
[32m+[m[32m    }[m
[32m+[m[32m};[m
[32m+[m
[32m+[m[32m} // namespace detail[m
[32m+[m[32m} // namespace nlohmann[m
[32m+[m
[32m+[m[32m// #include <nlohmann/detail/macro_scope.hpp>[m
[32m+[m
[32m+[m
[32m+[m[32m#include <utility> // pair[m
[32m+[m[32m// #include <nlohmann/thirdparty/hedley/hedley.hpp>[m
[32m+[m[32m/* Hedley - https://nemequ.github.io/hedley[m
[32m+[m[32m * Created by Evan Nemerson <evan@nemerson.com>[m
[32m+[m[32m *[m
[32m+[m[32m * To the extent possible under law, the author(s) have dedicated all[m
[32m+[m[32m * copyright and related and neighboring rights to this software to[m
[32m+[m[32m * the public domain worldwide. This software is distributed without[m
[32m+[m[32m * any warranty.[m
[32m+[m[32m *[m
[32m+[m[32m * For details, see <http://creativecommons.org/publicdomain/zero/1.0/>.[m
[32m+[m[32m * SPDX-License-Identifier: CC0-1.0[m
[32m+[m[32m */[m
[32m+[m
[32m+[m[32m#if !defined(JSON_HEDLEY_VERSION) || (JSON_HEDLEY_VERSION < 9)[m
[32m+[m[32m#if defined(JSON_HEDLEY_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_VERSION 9[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_STRINGIFY_EX)[m
[32m+[m[32m    #undef JSON_HEDLEY_STRINGIFY_EX[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_STRINGIFY_EX(x) #x[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_STRINGIFY)[m
[32m+[m[32m    #undef JSON_HEDLEY_STRINGIFY[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_STRINGIFY(x) JSON_HEDLEY_STRINGIFY_EX(x)[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CONCAT_EX)[m
[32m+[m[32m    #undef JSON_HEDLEY_CONCAT_EX[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_CONCAT_EX(a,b) a##b[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CONCAT)[m
[32m+[m[32m    #undef JSON_HEDLEY_CONCAT[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_CONCAT(a,b) JSON_HEDLEY_CONCAT_EX(a,b)[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_VERSION_ENCODE)[m
[32m+[m[32m    #undef JSON_HEDLEY_VERSION_ENCODE[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_VERSION_ENCODE(major,minor,revision) (((major) * 1000000) + ((minor) * 1000) + (revision))[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_VERSION_DECODE_MAJOR)[m
[32m+[m[32m    #undef JSON_HEDLEY_VERSION_DECODE_MAJOR[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_VERSION_DECODE_MAJOR(version) ((version) / 1000000)[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_VERSION_DECODE_MINOR)[m
[32m+[m[32m    #undef JSON_HEDLEY_VERSION_DECODE_MINOR[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_VERSION_DECODE_MINOR(version) (((version) % 1000000) / 1000)[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_VERSION_DECODE_REVISION)[m
[32m+[m[32m    #undef JSON_HEDLEY_VERSION_DECODE_REVISION[m
[32m+[m[32m#endif[m
[32m+[m[32m#define JSON_HEDLEY_VERSION_DECODE_REVISION(version) ((version) % 1000)[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__GNUC__) && defined(__GNUC_PATCHLEVEL__)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_VERSION JSON_HEDLEY_VERSION_ENCODE(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)[m
[32m+[m[32m#elif defined(__GNUC__)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_VERSION JSON_HEDLEY_VERSION_ENCODE(__GNUC__, __GNUC_MINOR__, 0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_GNUC_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_MSVC_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_MSVC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(_MSC_FULL_VER) && (_MSC_FULL_VER >= 140000000)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION JSON_HEDLEY_VERSION_ENCODE(_MSC_FULL_VER / 10000000, (_MSC_FULL_VER % 10000000) / 100000, (_MSC_FULL_VER % 100000) / 100)[m
[32m+[m[32m#elif defined(_MSC_FULL_VER)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION JSON_HEDLEY_VERSION_ENCODE(_MSC_FULL_VER / 1000000, (_MSC_FULL_VER % 1000000) / 10000, (_MSC_FULL_VER % 10000) / 10)[m
[32m+[m[32m#elif defined(_MSC_VER)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION JSON_HEDLEY_VERSION_ENCODE(_MSC_VER / 100, _MSC_VER % 100, 0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_MSVC_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_MSVC_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(_MSC_VER)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#elif defined(_MSC_VER) && (_MSC_VER >= 1400)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION_CHECK(major,minor,patch) (_MSC_FULL_VER >= ((major * 10000000) + (minor * 100000) + (patch)))[m
[32m+[m[32m#elif defined(_MSC_VER) && (_MSC_VER >= 1200)[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION_CHECK(major,minor,patch) (_MSC_FULL_VER >= ((major * 1000000) + (minor * 10000) + (patch)))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_MSVC_VERSION_CHECK(major,minor,patch) (_MSC_VER >= ((major * 100) + (minor)))[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_INTEL_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_INTEL_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__INTEL_COMPILER) && defined(__INTEL_COMPILER_UPDATE)[m
[32m+[m[32m    #define JSON_HEDLEY_INTEL_VERSION JSON_HEDLEY_VERSION_ENCODE(__INTEL_COMPILER / 100, __INTEL_COMPILER % 100, __INTEL_COMPILER_UPDATE)[m
[32m+[m[32m#elif defined(__INTEL_COMPILER)[m
[32m+[m[32m    #define JSON_HEDLEY_INTEL_VERSION JSON_HEDLEY_VERSION_ENCODE(__INTEL_COMPILER / 100, __INTEL_COMPILER % 100, 0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_INTEL_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_INTEL_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_INTEL_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_INTEL_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_INTEL_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_INTEL_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PGI_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_PGI_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__PGI) && defined(__PGIC__) && defined(__PGIC_MINOR__) && defined(__PGIC_PATCHLEVEL__)[m
[32m+[m[32m    #define JSON_HEDLEY_PGI_VERSION JSON_HEDLEY_VERSION_ENCODE(__PGIC__, __PGIC_MINOR__, __PGIC_PATCHLEVEL__)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PGI_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_PGI_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_PGI_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_PGI_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_PGI_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_PGI_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_SUNPRO_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_SUNPRO_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__SUNPRO_C) && (__SUNPRO_C > 0x1000)[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION JSON_HEDLEY_VERSION_ENCODE((((__SUNPRO_C >> 16) & 0xf) * 10) + ((__SUNPRO_C >> 12) & 0xf), (((__SUNPRO_C >> 8) & 0xf) * 10) + ((__SUNPRO_C >> 4) & 0xf), (__SUNPRO_C & 0xf) * 10)[m
[32m+[m[32m#elif defined(__SUNPRO_C)[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION JSON_HEDLEY_VERSION_ENCODE((__SUNPRO_C >> 8) & 0xf, (__SUNPRO_C >> 4) & 0xf, (__SUNPRO_C) & 0xf)[m
[32m+[m[32m#elif defined(__SUNPRO_CC) && (__SUNPRO_CC > 0x1000)[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION JSON_HEDLEY_VERSION_ENCODE((((__SUNPRO_CC >> 16) & 0xf) * 10) + ((__SUNPRO_CC >> 12) & 0xf), (((__SUNPRO_CC >> 8) & 0xf) * 10) + ((__SUNPRO_CC >> 4) & 0xf), (__SUNPRO_CC & 0xf) * 10)[m
[32m+[m[32m#elif defined(__SUNPRO_CC)[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION JSON_HEDLEY_VERSION_ENCODE((__SUNPRO_CC >> 8) & 0xf, (__SUNPRO_CC >> 4) & 0xf, (__SUNPRO_CC) & 0xf)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_SUNPRO_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_SUNPRO_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_SUNPRO_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_SUNPRO_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_SUNPRO_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_EMSCRIPTEN_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_EMSCRIPTEN_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__EMSCRIPTEN__)[m
[32m+[m[32m    #define JSON_HEDLEY_EMSCRIPTEN_VERSION JSON_HEDLEY_VERSION_ENCODE(__EMSCRIPTEN_major__, __EMSCRIPTEN_minor__, __EMSCRIPTEN_tiny__)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_EMSCRIPTEN_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_EMSCRIPTEN_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_EMSCRIPTEN_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_EMSCRIPTEN_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_EMSCRIPTEN_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_EMSCRIPTEN_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_ARM_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_ARM_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__CC_ARM) && defined(__ARMCOMPILER_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_ARM_VERSION JSON_HEDLEY_VERSION_ENCODE(__ARMCOMPILER_VERSION / 1000000, (__ARMCOMPILER_VERSION % 1000000) / 10000, (__ARMCOMPILER_VERSION % 10000) / 100)[m
[32m+[m[32m#elif defined(__CC_ARM) && defined(__ARMCC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_ARM_VERSION JSON_HEDLEY_VERSION_ENCODE(__ARMCC_VERSION / 1000000, (__ARMCC_VERSION % 1000000) / 10000, (__ARMCC_VERSION % 10000) / 100)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_ARM_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_ARM_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_ARM_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_ARM_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_ARM_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_ARM_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_IBM_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_IBM_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__ibmxl__)[m
[32m+[m[32m    #define JSON_HEDLEY_IBM_VERSION JSON_HEDLEY_VERSION_ENCODE(__ibmxl_version__, __ibmxl_release__, __ibmxl_modification__)[m
[32m+[m[32m#elif defined(__xlC__) && defined(__xlC_ver__)[m
[32m+[m[32m    #define JSON_HEDLEY_IBM_VERSION JSON_HEDLEY_VERSION_ENCODE(__xlC__ >> 8, __xlC__ & 0xff, (__xlC_ver__ >> 8) & 0xff)[m
[32m+[m[32m#elif defined(__xlC__)[m
[32m+[m[32m    #define JSON_HEDLEY_IBM_VERSION JSON_HEDLEY_VERSION_ENCODE(__xlC__ >> 8, __xlC__ & 0xff, 0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_IBM_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_IBM_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_IBM_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_IBM_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_IBM_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_IBM_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_TI_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_TI_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__TI_COMPILER_VERSION__)[m
[32m+[m[32m    #define JSON_HEDLEY_TI_VERSION JSON_HEDLEY_VERSION_ENCODE(__TI_COMPILER_VERSION__ / 1000000, (__TI_COMPILER_VERSION__ % 1000000) / 1000, (__TI_COMPILER_VERSION__ % 1000))[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_TI_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_TI_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_TI_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_TI_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_TI_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_TI_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CRAY_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_CRAY_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(_CRAYC)[m
[32m+[m[32m    #if defined(_RELEASE_PATCHLEVEL)[m
[32m+[m[32m        #define JSON_HEDLEY_CRAY_VERSION JSON_HEDLEY_VERSION_ENCODE(_RELEASE_MAJOR, _RELEASE_MINOR, _RELEASE_PATCHLEVEL)[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_CRAY_VERSION JSON_HEDLEY_VERSION_ENCODE(_RELEASE_MAJOR, _RELEASE_MINOR, 0)[m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CRAY_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_CRAY_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_CRAY_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_CRAY_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_CRAY_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_CRAY_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_IAR_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_IAR_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__IAR_SYSTEMS_ICC__)[m
[32m+[m[32m    #if __VER__ > 1000[m
[32m+[m[32m        #define JSON_HEDLEY_IAR_VERSION JSON_HEDLEY_VERSION_ENCODE((__VER__ / 1000000), ((__VER__ / 1000) % 1000), (__VER__ % 1000))[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_IAR_VERSION JSON_HEDLEY_VERSION_ENCODE(VER / 100, __VER__ % 100, 0)[m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_IAR_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_IAR_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_IAR_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_IAR_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_IAR_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_IAR_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_TINYC_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_TINYC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__TINYC__)[m
[32m+[m[32m    #define JSON_HEDLEY_TINYC_VERSION JSON_HEDLEY_VERSION_ENCODE(__TINYC__ / 1000, (__TINYC__ / 100) % 10, __TINYC__ % 100)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_TINYC_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_TINYC_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_TINYC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_TINYC_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_TINYC_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_TINYC_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DMC_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_DMC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__DMC__)[m
[32m+[m[32m    #define JSON_HEDLEY_DMC_VERSION JSON_HEDLEY_VERSION_ENCODE(__DMC__ >> 8, (__DMC__ >> 4) & 0xf, __DMC__ & 0xf)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DMC_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_DMC_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_DMC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_DMC_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_DMC_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DMC_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_COMPCERT_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_COMPCERT_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__COMPCERT_VERSION__)[m
[32m+[m[32m    #define JSON_HEDLEY_COMPCERT_VERSION JSON_HEDLEY_VERSION_ENCODE(__COMPCERT_VERSION__ / 10000, (__COMPCERT_VERSION__ / 100) % 100, __COMPCERT_VERSION__ % 100)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_COMPCERT_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_COMPCERT_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_COMPCERT_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_COMPCERT_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_COMPCERT_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_COMPCERT_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PELLES_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_PELLES_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__POCC__)[m
[32m+[m[32m    #define JSON_HEDLEY_PELLES_VERSION JSON_HEDLEY_VERSION_ENCODE(__POCC__ / 100, __POCC__ % 100, 0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PELLES_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_PELLES_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_PELLES_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_PELLES_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_PELLES_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_PELLES_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_VERSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    defined(JSON_HEDLEY_GNUC_VERSION) && \[m
[32m+[m[32m    !defined(__clang__) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_INTEL_VERSION) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_PGI_VERSION) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_ARM_VERSION) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_TI_VERSION) && \[m
[32m+[m[32m    !defined(__COMPCERT__)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_VERSION JSON_HEDLEY_GNUC_VERSION[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_VERSION_CHECK)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_VERSION_CHECK[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch) (JSON_HEDLEY_GCC_VERSION >= JSON_HEDLEY_VERSION_ENCODE(major, minor, patch))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_ATTRIBUTE(attribute) __has_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_ATTRIBUTE(attribute) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_ATTRIBUTE(attribute,major,minor,patch) __has_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_ATTRIBUTE(attribute,major,minor,patch) __has_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_CPP_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_CPP_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_cpp_attribute) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_CPP_ATTRIBUTE(attribute) __has_cpp_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_CPP_ATTRIBUTE(attribute) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_CPP_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_CPP_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_cpp_attribute) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_CPP_ATTRIBUTE(attribute,major,minor,patch) __has_cpp_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_CPP_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_CPP_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_CPP_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_cpp_attribute) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_CPP_ATTRIBUTE(attribute,major,minor,patch) __has_cpp_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_CPP_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_BUILTIN)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_BUILTIN[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_builtin)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_BUILTIN(builtin) __has_builtin(builtin)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_BUILTIN(builtin) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_BUILTIN)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_BUILTIN[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_builtin)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_BUILTIN(builtin,major,minor,patch) __has_builtin(builtin)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_BUILTIN(builtin,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_BUILTIN)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_BUILTIN[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_builtin)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_BUILTIN(builtin,major,minor,patch) __has_builtin(builtin)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_BUILTIN(builtin,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_FEATURE)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_FEATURE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_feature)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_FEATURE(feature) __has_feature(feature)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_FEATURE(feature) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_FEATURE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_FEATURE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_feature)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_FEATURE(feature,major,minor,patch) __has_feature(feature)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_FEATURE(feature,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_FEATURE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_FEATURE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_feature)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_FEATURE(feature,major,minor,patch) __has_feature(feature)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_FEATURE(feature,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_EXTENSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_EXTENSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_extension)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_EXTENSION(extension) __has_extension(extension)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_EXTENSION(extension) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_EXTENSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_EXTENSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_extension)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_EXTENSION(extension,major,minor,patch) __has_extension(extension)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_EXTENSION(extension,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_EXTENSION)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_EXTENSION[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_extension)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_EXTENSION(extension,major,minor,patch) __has_extension(extension)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_EXTENSION(extension,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_DECLSPEC_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_DECLSPEC_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_declspec_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_DECLSPEC_ATTRIBUTE(attribute) __has_declspec_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_DECLSPEC_ATTRIBUTE(attribute) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_DECLSPEC_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_DECLSPEC_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_declspec_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_DECLSPEC_ATTRIBUTE(attribute,major,minor,patch) __has_declspec_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_DECLSPEC_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_DECLSPEC_ATTRIBUTE)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_DECLSPEC_ATTRIBUTE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_declspec_attribute)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_DECLSPEC_ATTRIBUTE(attribute,major,minor,patch) __has_declspec_attribute(attribute)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_DECLSPEC_ATTRIBUTE(attribute,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_HAS_WARNING)[m
[32m+[m[32m    #undef JSON_HEDLEY_HAS_WARNING[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_warning)[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_WARNING(warning) __has_warning(warning)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_HAS_WARNING(warning) (0)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GNUC_HAS_WARNING)[m
[32m+[m[32m    #undef JSON_HEDLEY_GNUC_HAS_WARNING[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_warning)[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_WARNING(warning,major,minor,patch) __has_warning(warning)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GNUC_HAS_WARNING(warning,major,minor,patch) JSON_HEDLEY_GNUC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_GCC_HAS_WARNING)[m
[32m+[m[32m    #undef JSON_HEDLEY_GCC_HAS_WARNING[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__has_warning)[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_WARNING(warning,major,minor,patch) __has_warning(warning)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_GCC_HAS_WARNING(warning,major,minor,patch) JSON_HEDLEY_GCC_VERSION_CHECK(major,minor,patch)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if \[m
[32m+[m[32m    (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || \[m
[32m+[m[32m    defined(__clang__) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(18,4,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(6,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_CRAY_VERSION_CHECK(5,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TINYC_VERSION_CHECK(0,9,17) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) && defined(__C99_PRAGMA_OPERATOR))[m
[32m+[m[32m    #define JSON_HEDLEY_PRAGMA(value) _Pragma(#value)[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(15,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_PRAGMA(value) __pragma(value)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_PRAGMA(value)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DIAGNOSTIC_PUSH)[m
[32m+[m[32m    #undef JSON_HEDLEY_DIAGNOSTIC_PUSH[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_DIAGNOSTIC_POP)[m
[32m+[m[32m    #undef JSON_HEDLEY_DIAGNOSTIC_POP[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__clang__)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("clang diagnostic push")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("clang diagnostic pop")[m
[32m+[m[32m#elif JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("warning(push)")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("warning(pop)")[m
[32m+[m[32m#elif JSON_HEDLEY_GCC_VERSION_CHECK(4,6,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("GCC diagnostic push")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("GCC diagnostic pop")[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(15,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH __pragma(warning(push))[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP __pragma(warning(pop))[m
[32m+[m[32m#elif JSON_HEDLEY_ARM_VERSION_CHECK(5,6,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("push")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("pop")[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(8,1,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("diag_push")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("diag_pop")[m
[32m+[m[32m#elif JSON_HEDLEY_PELLES_VERSION_CHECK(2,90,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH _Pragma("warning(push)")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP _Pragma("warning(pop)")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_PUSH[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_POP[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED)[m
[32m+[m[32m    #undef JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED[m
[32m+[m[32m#endif[m
[32m+[m[32m#if JSON_HEDLEY_HAS_WARNING("-Wdeprecated-declarations")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")[m
[32m+[m[32m#elif JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("warning(disable:1478 1786)")[m
[32m+[m[32m#elif JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("diag_suppress 1215,1444")[m
[32m+[m[32m#elif JSON_HEDLEY_GCC_VERSION_CHECK(4,3,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(15,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED __pragma(warning(disable:4996))[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("diag_suppress 1291,1718")[m
[32m+[m[32m#elif JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,13,0) && !defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("error_messages(off,E_DEPRECATED_ATT,E_DEPRECATED_ATT_MESS)")[m
[32m+[m[32m#elif JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,13,0) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("error_messages(off,symdeprecated,symdeprecated2)")[m
[32m+[m[32m#elif JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("diag_suppress=Pe1444,Pe1215")[m
[32m+[m[32m#elif JSON_HEDLEY_PELLES_VERSION_CHECK(2,90,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED _Pragma("warn(disable:2241)")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_DEPRECATED[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS)[m
[32m+[m[32m    #undef JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS[m
[32m+[m[32m#endif[m
[32m+[m[32m#if JSON_HEDLEY_HAS_WARNING("-Wunknown-pragmas")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("clang diagnostic ignored \"-Wunknown-pragmas\"")[m
[32m+[m[32m#elif JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("warning(disable:161)")[m
[32m+[m[32m#elif JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("diag_suppress 1675")[m
[32m+[m[32m#elif JSON_HEDLEY_GCC_VERSION_CHECK(4,3,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(15,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS __pragma(warning(disable:4068))[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("diag_suppress 163")[m
[32m+[m[32m#elif JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS _Pragma("diag_suppress=Pe161")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_UNKNOWN_PRAGMAS[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL)[m
[32m+[m[32m    #undef JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL[m
[32m+[m[32m#endif[m
[32m+[m[32m#if JSON_HEDLEY_HAS_WARNING("-Wcast-qual")[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL _Pragma("clang diagnostic ignored \"-Wcast-qual\"")[m
[32m+[m[32m#elif JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL _Pragma("warning(disable:2203 2331)")[m
[32m+[m[32m#elif JSON_HEDLEY_GCC_VERSION_CHECK(3,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL _Pragma("GCC diagnostic ignored \"-Wcast-qual\"")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_DEPRECATED)[m
[32m+[m[32m    #undef JSON_HEDLEY_DEPRECATED[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_DEPRECATED_FOR)[m
[32m+[m[32m    #undef JSON_HEDLEY_DEPRECATED_FOR[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus) && (__cplusplus >= 201402L)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) [[deprecated("Since " #since)]][m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) [[deprecated("Since " #since "; use " #replacement)]][m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_HAS_EXTENSION(attribute_deprecated_with_message) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,5,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(5,6,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,13,0) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,3,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) __attribute__((__deprecated__("Since " #since)))[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__("Since " #since "; use " #replacement)))[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(deprecated) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) __attribute__((__deprecated__))[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) __attribute__((__deprecated__))[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(14,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) __declspec(deprecated("Since " # since))[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) __declspec(deprecated("Since " #since "; use " #replacement))[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_MSVC_VERSION_CHECK(13,10,0) || \[m
[32m+[m[32m    JSON_HEDLEY_PELLES_VERSION_CHECK(6,50,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) _declspec(deprecated)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) __declspec(deprecated)[m
[32m+[m[32m#elif JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since) _Pragma("deprecated")[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement) _Pragma("deprecated")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED(since)[m
[32m+[m[32m    #define JSON_HEDLEY_DEPRECATED_FOR(since, replacement)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_UNAVAILABLE)[m
[32m+[m[32m    #undef JSON_HEDLEY_UNAVAILABLE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(warning) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,3,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_UNAVAILABLE(available_since) __attribute__((__warning__("Not available until " #available_since)))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_UNAVAILABLE(available_since)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_WARN_UNUSED_RESULT)[m
[32m+[m[32m    #undef JSON_HEDLEY_WARN_UNUSED_RESULT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus) && (__cplusplus >= 201703L)[m
[32m+[m[32m    #define JSON_HEDLEY_WARN_UNUSED_RESULT [[nodiscard]][m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(warn_unused_result) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,4,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__)) || \[m
[32m+[m[32m    (JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,15,0) && defined(__cplusplus)) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))[m
[32m+[m[32m#elif defined(_Check_return_) /* SAL */[m
[32m+[m[32m    #define JSON_HEDLEY_WARN_UNUSED_RESULT _Check_return_[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_WARN_UNUSED_RESULT[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_SENTINEL)[m
[32m+[m[32m    #undef JSON_HEDLEY_SENTINEL[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(sentinel) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(5,4,0)[m
[32m+[m[32m    #define JSON_HEDLEY_SENTINEL(position) __attribute__((__sentinel__(position)))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_SENTINEL(position)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_NO_RETURN)[m
[32m+[m[32m    #undef JSON_HEDLEY_NO_RETURN[m
[32m+[m[32m#endif[m
[32m+[m[32m#if JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __noreturn[m
[32m+[m[32m#elif JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __attribute__((__noreturn__))[m
[32m+[m[32m#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN _Noreturn[m
[32m+[m[32m#elif defined(__cplusplus) && (__cplusplus >= 201103L)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN [[noreturn]][m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(noreturn) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,2,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(18,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(17,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __attribute__((__noreturn__))[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(13,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __declspec(noreturn)[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(6,0,0) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN _Pragma("FUNC_NEVER_RETURNS;")[m
[32m+[m[32m#elif JSON_HEDLEY_COMPCERT_VERSION_CHECK(3,2,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __attribute((noreturn))[m
[32m+[m[32m#elif JSON_HEDLEY_PELLES_VERSION_CHECK(9,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN __declspec(noreturn)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_NO_RETURN[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_UNREACHABLE)[m
[32m+[m[32m    #undef JSON_HEDLEY_UNREACHABLE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_UNREACHABLE_RETURN)[m
[32m+[m[32m    #undef JSON_HEDLEY_UNREACHABLE_RETURN[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    (JSON_HEDLEY_HAS_BUILTIN(__builtin_unreachable) && (!defined(JSON_HEDLEY_ARM_VERSION))) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,5,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(13,1,5)[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE() __builtin_unreachable()[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(13,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE() __assume(0)[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(6,0,0)[m
[32m+[m[32m    #if defined(__cplusplus)[m
[32m+[m[32m        #define JSON_HEDLEY_UNREACHABLE() std::_nassert(0)[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_UNREACHABLE() _nassert(0)[m
[32m+[m[32m    #endif[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE_RETURN(value) return value[m
[32m+[m[32m#elif defined(EXIT_FAILURE)[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE() abort()[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE()[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE_RETURN(value) return value[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(JSON_HEDLEY_UNREACHABLE_RETURN)[m
[32m+[m[32m    #define JSON_HEDLEY_UNREACHABLE_RETURN(value) JSON_HEDLEY_UNREACHABLE()[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_ASSUME)[m
[32m+[m[32m    #undef JSON_HEDLEY_ASSUME[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_MSVC_VERSION_CHECK(13,10,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_ASSUME(expr) __assume(expr)[m
[32m+[m[32m#elif JSON_HEDLEY_HAS_BUILTIN(__builtin_assume)[m
[32m+[m[32m    #define JSON_HEDLEY_ASSUME(expr) __builtin_assume(expr)[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(6,0,0)[m
[32m+[m[32m    #if defined(__cplusplus)[m
[32m+[m[32m        #define JSON_HEDLEY_ASSUME(expr) std::_nassert(expr)[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_ASSUME(expr) _nassert(expr)[m
[32m+[m[32m    #endif[m
[32m+[m[32m#elif \[m
[32m+[m[32m    (JSON_HEDLEY_HAS_BUILTIN(__builtin_unreachable) && !defined(JSON_HEDLEY_ARM_VERSION)) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,5,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(13,1,5)[m
[32m+[m[32m    #define JSON_HEDLEY_ASSUME(expr) ((void) ((expr) ? 1 : (__builtin_unreachable(), 1)))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_ASSUME(expr) ((void) (expr))[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m
[32m+[m[32mJSON_HEDLEY_DIAGNOSTIC_PUSH[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_WARNING("-Wvariadic-macros") || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,0,0)[m
[32m+[m[32m    #if defined(__clang__)[m
[32m+[m[32m        #pragma clang diagnostic ignored "-Wvariadic-macros"[m
[32m+[m[32m    #elif defined(JSON_HEDLEY_GCC_VERSION)[m
[32m+[m[32m        #pragma GCC diagnostic ignored "-Wvariadic-macros"[m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_NON_NULL)[m
[32m+[m[32m    #undef JSON_HEDLEY_NON_NULL[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(nonnull) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,3,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NON_NULL(...) __attribute__((__nonnull__(__VA_ARGS__)))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_NON_NULL(...)[m
[32m+[m[32m#endif[m
[32m+[m[32mJSON_HEDLEY_DIAGNOSTIC_POP[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PRINTF_FORMAT)[m
[32m+[m[32m    #undef JSON_HEDLEY_PRINTF_FORMAT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__MINGW32__) && JSON_HEDLEY_GCC_HAS_ATTRIBUTE(format,4,4,0) && !defined(__USE_MINGW_ANSI_STDIO)[m
[32m+[m[32m    #define JSON_HEDLEY_PRINTF_FORMAT(string_idx,first_to_check) __attribute__((__format__(ms_printf, string_idx, first_to_check)))[m
[32m+[m[32m#elif defined(__MINGW32__) && JSON_HEDLEY_GCC_HAS_ATTRIBUTE(format,4,4,0) && defined(__USE_MINGW_ANSI_STDIO)[m
[32m+[m[32m    #define JSON_HEDLEY_PRINTF_FORMAT(string_idx,first_to_check) __attribute__((__format__(gnu_printf, string_idx, first_to_check)))[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(format) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(5,6,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_PRINTF_FORMAT(string_idx,first_to_check) __attribute__((__format__(__printf__, string_idx, first_to_check)))[m
[32m+[m[32m#elif JSON_HEDLEY_PELLES_VERSION_CHECK(6,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_PRINTF_FORMAT(string_idx,first_to_check) __declspec(vaformat(printf,string_idx,first_to_check))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_PRINTF_FORMAT(string_idx,first_to_check)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CONSTEXPR)[m
[32m+[m[32m    #undef JSON_HEDLEY_CONSTEXPR[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus)[m
[32m+[m[32m    #if __cplusplus >= 201103L[m
[32m+[m[32m        #define JSON_HEDLEY_CONSTEXPR constexpr[m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(JSON_HEDLEY_CONSTEXPR)[m
[32m+[m[32m    #define JSON_HEDLEY_CONSTEXPR[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PREDICT)[m
[32m+[m[32m    #undef JSON_HEDLEY_PREDICT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_LIKELY)[m
[32m+[m[32m    #undef JSON_HEDLEY_LIKELY[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_UNLIKELY)[m
[32m+[m[32m    #undef JSON_HEDLEY_UNLIKELY[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_UNPREDICTABLE)[m
[32m+[m[32m    #undef JSON_HEDLEY_UNPREDICTABLE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if JSON_HEDLEY_HAS_BUILTIN(__builtin_unpredictable)[m
[32m+[m[32m    #define JSON_HEDLEY_UNPREDICTABLE(expr) __builtin_unpredictable(!!(expr))[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m  JSON_HEDLEY_HAS_BUILTIN(__builtin_expect_with_probability) || \[m
[32m+[m[32m  JSON_HEDLEY_GCC_VERSION_CHECK(9,0,0)[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT(expr, value, probability) __builtin_expect_with_probability(expr, value, probability)[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_TRUE(expr, probability) __builtin_expect_with_probability(!!(expr), 1, probability)[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_FALSE(expr, probability) __builtin_expect_with_probability(!!(expr), 0, probability)[m
[32m+[m[32m#  define JSON_HEDLEY_LIKELY(expr) __builtin_expect(!!(expr), 1)[m
[32m+[m[32m#  define JSON_HEDLEY_UNLIKELY(expr) __builtin_expect(!!(expr), 0)[m
[32m+[m[32m#if !defined(JSON_HEDLEY_BUILTIN_UNPREDICTABLE)[m
[32m+[m[32m    #define JSON_HEDLEY_BUILTIN_UNPREDICTABLE(expr) __builtin_expect_with_probability(!!(expr), 1, 0.5)[m
[32m+[m[32m#endif[m
[32m+[m[32m#elif \[m
[32m+[m[32m  JSON_HEDLEY_HAS_BUILTIN(__builtin_expect) || \[m
[32m+[m[32m  JSON_HEDLEY_GCC_VERSION_CHECK(3,0,0) || \[m
[32m+[m[32m  JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m  (JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,15,0) && defined(__cplusplus)) || \[m
[32m+[m[32m  JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m  JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m  JSON_HEDLEY_TI_VERSION_CHECK(6,1,0) || \[m
[32m+[m[32m  JSON_HEDLEY_TINYC_VERSION_CHECK(0,9,27)[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT(expr, expected, probability) \[m
[32m+[m[32m    (((probability) >= 0.9) ? __builtin_expect(!!(expr), (expected)) : (((void) (expected)), !!(expr)))[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_TRUE(expr, probability) \[m
[32m+[m[32m    (__extension__ ({ \[m
[32m+[m[32m        JSON_HEDLEY_CONSTEXPR double hedley_probability_ = (probability); \[m
[32m+[m[32m        ((hedley_probability_ >= 0.9) ? __builtin_expect(!!(expr), 1) : ((hedley_probability_ <= 0.1) ? __builtin_expect(!!(expr), 0) : !!(expr))); \[m
[32m+[m[32m    }))[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_FALSE(expr, probability) \[m
[32m+[m[32m    (__extension__ ({ \[m
[32m+[m[32m        JSON_HEDLEY_CONSTEXPR double hedley_probability_ = (probability); \[m
[32m+[m[32m        ((hedley_probability_ >= 0.9) ? __builtin_expect(!!(expr), 0) : ((hedley_probability_ <= 0.1) ? __builtin_expect(!!(expr), 1) : !!(expr))); \[m
[32m+[m[32m    }))[m
[32m+[m[32m#  define JSON_HEDLEY_LIKELY(expr)   __builtin_expect(!!(expr), 1)[m
[32m+[m[32m#  define JSON_HEDLEY_UNLIKELY(expr) __builtin_expect(!!(expr), 0)[m
[32m+[m[32m#else[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT(expr, expected, probability) (((void) (expected)), !!(expr))[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_TRUE(expr, probability) (!!(expr))[m
[32m+[m[32m#  define JSON_HEDLEY_PREDICT_FALSE(expr, probability) (!!(expr))[m
[32m+[m[32m#  define JSON_HEDLEY_LIKELY(expr) (!!(expr))[m
[32m+[m[32m#  define JSON_HEDLEY_UNLIKELY(expr) (!!(expr))[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(JSON_HEDLEY_UNPREDICTABLE)[m
[32m+[m[32m    #define JSON_HEDLEY_UNPREDICTABLE(expr) JSON_HEDLEY_PREDICT(expr, 1, 0.5)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_MALLOC)[m
[32m+[m[32m    #undef JSON_HEDLEY_MALLOC[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(malloc) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(12,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_MALLOC __attribute__((__malloc__))[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(14, 0, 0)[m
[32m+[m[32m    #define JSON_HEDLEY_MALLOC __declspec(restrict)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_MALLOC[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PURE)[m
[32m+[m[32m    #undef JSON_HEDLEY_PURE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(pure) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(2,96,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__)) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_PURE __attribute__((__pure__))[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(6,0,0) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_PURE _Pragma("FUNC_IS_PURE;")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_PURE[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CONST)[m
[32m+[m[32m    #undef JSON_HEDLEY_CONST[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(const) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(2,5,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__)) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_CONST __attribute__((__const__))[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_CONST JSON_HEDLEY_PURE[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_RESTRICT)[m
[32m+[m[32m    #undef JSON_HEDLEY_RESTRICT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) && !defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_RESTRICT restrict[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_MSVC_VERSION_CHECK(14,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_PGI_VERSION_CHECK(17,10,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,14,0) && defined(__cplusplus)) || \[m
[32m+[m[32m    JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    defined(__clang__)[m
[32m+[m[32m    #define JSON_HEDLEY_RESTRICT __restrict[m
[32m+[m[32m#elif JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,3,0) && !defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_RESTRICT _Restrict[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_RESTRICT[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_INLINE)[m
[32m+[m[32m    #undef JSON_HEDLEY_INLINE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)) || \[m
[32m+[m[32m    (defined(__cplusplus) && (__cplusplus >= 199711L))[m
[32m+[m[32m    #define JSON_HEDLEY_INLINE inline[m
[32m+[m[32m#elif \[m
[32m+[m[32m    defined(JSON_HEDLEY_GCC_VERSION) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(6,2,0)[m
[32m+[m[32m    #define JSON_HEDLEY_INLINE __inline__[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_MSVC_VERSION_CHECK(12,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_INLINE __inline[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_INLINE[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_ALWAYS_INLINE)[m
[32m+[m[32m    #undef JSON_HEDLEY_ALWAYS_INLINE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(always_inline) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_ALWAYS_INLINE __attribute__((__always_inline__)) JSON_HEDLEY_INLINE[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(12,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_ALWAYS_INLINE __forceinline[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(7,0,0) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_ALWAYS_INLINE _Pragma("FUNC_ALWAYS_INLINE;")[m
[32m+[m[32m#elif JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_ALWAYS_INLINE _Pragma("inline=forced")[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_ALWAYS_INLINE JSON_HEDLEY_INLINE[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_NEVER_INLINE)[m
[32m+[m[32m    #undef JSON_HEDLEY_NEVER_INLINE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(noinline) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(10,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m    (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE __attribute__((__noinline__))[m
[32m+[m[32m#elif JSON_HEDLEY_MSVC_VERSION_CHECK(13,10,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE __declspec(noinline)[m
[32m+[m[32m#elif JSON_HEDLEY_PGI_VERSION_CHECK(10,2,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE _Pragma("noinline")[m
[32m+[m[32m#elif JSON_HEDLEY_TI_VERSION_CHECK(6,0,0) && defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE _Pragma("FUNC_CANNOT_INLINE;")[m
[32m+[m[32m#elif JSON_HEDLEY_IAR_VERSION_CHECK(8,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE _Pragma("inline=never")[m
[32m+[m[32m#elif JSON_HEDLEY_COMPCERT_VERSION_CHECK(3,2,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE __attribute((noinline))[m
[32m+[m[32m#elif JSON_HEDLEY_PELLES_VERSION_CHECK(9,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE __declspec(noinline)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_NEVER_INLINE[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_PRIVATE)[m
[32m+[m[32m    #undef JSON_HEDLEY_PRIVATE[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_PUBLIC)[m
[32m+[m[32m    #undef JSON_HEDLEY_PUBLIC[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_IMPORT)[m
[32m+[m[32m    #undef JSON_HEDLEY_IMPORT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(_WIN32) || defined(__CYGWIN__)[m
[32m+[m[32m    #define JSON_HEDLEY_PRIVATE[m
[32m+[m[32m    #define JSON_HEDLEY_PUBLIC   __declspec(dllexport)[m
[32m+[m[32m    #define JSON_HEDLEY_IMPORT   __declspec(dllimport)[m
[32m+[m[32m#else[m
[32m+[m[32m    #if \[m
[32m+[m[32m        JSON_HEDLEY_HAS_ATTRIBUTE(visibility) || \[m
[32m+[m[32m        JSON_HEDLEY_GCC_VERSION_CHECK(3,3,0) || \[m
[32m+[m[32m        JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,11,0) || \[m
[32m+[m[32m        JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m        JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m        JSON_HEDLEY_IBM_VERSION_CHECK(13,1,0) || \[m
[32m+[m[32m        JSON_HEDLEY_TI_VERSION_CHECK(8,0,0) || \[m
[32m+[m[32m        (JSON_HEDLEY_TI_VERSION_CHECK(7,3,0) && defined(__TI_EABI__) && defined(__TI_GNU_ATTRIBUTE_SUPPORT__))[m
[32m+[m[32m        #define JSON_HEDLEY_PRIVATE __attribute__((__visibility__("hidden")))[m
[32m+[m[32m        #define JSON_HEDLEY_PUBLIC  __attribute__((__visibility__("default")))[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_PRIVATE[m
[32m+[m[32m        #define JSON_HEDLEY_PUBLIC[m
[32m+[m[32m    #endif[m
[32m+[m[32m    #define JSON_HEDLEY_IMPORT    extern[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_NO_THROW)[m
[32m+[m[32m    #undef JSON_HEDLEY_NO_THROW[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(nothrow) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,3,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_THROW __attribute__((__nothrow__))[m
[32m+[m[32m#elif \[m
[32m+[m[32m    JSON_HEDLEY_MSVC_VERSION_CHECK(13,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0)[m
[32m+[m[32m    #define JSON_HEDLEY_NO_THROW __declspec(nothrow)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_NO_THROW[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_FALL_THROUGH)[m
[32m+[m[32m    #undef JSON_HEDLEY_FALL_THROUGH[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    defined(__cplusplus) && \[m
[32m+[m[32m    (!defined(JSON_HEDLEY_SUNPRO_VERSION) || JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,15,0)) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_PGI_VERSION)[m
[32m+[m[32m    #if \[m
[32m+[m[32m        (__cplusplus >= 201703L) || \[m
[32m+[m[32m        ((__cplusplus >= 201103L) && JSON_HEDLEY_HAS_CPP_ATTRIBUTE(fallthrough))[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH [[fallthrough]][m
[32m+[m[32m    #elif (__cplusplus >= 201103L) && JSON_HEDLEY_HAS_CPP_ATTRIBUTE(clang::fallthrough)[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH [[clang::fallthrough]][m
[32m+[m[32m    #elif (__cplusplus >= 201103L) && JSON_HEDLEY_GCC_VERSION_CHECK(7,0,0)[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH [[gnu::fallthrough]][m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(JSON_HEDLEY_FALL_THROUGH)[m
[32m+[m[32m    #if JSON_HEDLEY_GNUC_HAS_ATTRIBUTE(fallthrough,7,0,0) && !defined(JSON_HEDLEY_PGI_VERSION)[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH __attribute__((__fallthrough__))[m
[32m+[m[32m    #elif defined(__fallthrough) /* SAL */[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH __fallthrough[m
[32m+[m[32m    #else[m
[32m+[m[32m        #define JSON_HEDLEY_FALL_THROUGH[m
[32m+[m[32m    #endif[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_RETURNS_NON_NULL)[m
[32m+[m[32m    #undef JSON_HEDLEY_RETURNS_NON_NULL[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_ATTRIBUTE(returns_nonnull) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(4,9,0)[m
[32m+[m[32m    #define JSON_HEDLEY_RETURNS_NON_NULL __attribute__((__returns_nonnull__))[m
[32m+[m[32m#elif defined(_Ret_notnull_) /* SAL */[m
[32m+[m[32m    #define JSON_HEDLEY_RETURNS_NON_NULL _Ret_notnull_[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_RETURNS_NON_NULL[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_ARRAY_PARAM)[m
[32m+[m[32m    #undef JSON_HEDLEY_ARRAY_PARAM[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m    defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) && \[m
[32m+[m[32m    !defined(__STDC_NO_VLA__) && \[m
[32m+[m[32m    !defined(__cplusplus) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_PGI_VERSION) && \[m
[32m+[m[32m    !defined(JSON_HEDLEY_TINYC_VERSION)[m
[32m+[m[32m    #define JSON_HEDLEY_ARRAY_PARAM(name) (name)[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_ARRAY_PARAM(name)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_IS_CONSTANT)[m
[32m+[m[32m    #undef JSON_HEDLEY_IS_CONSTANT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_REQUIRE_CONSTEXPR)[m
[32m+[m[32m    #undef JSON_HEDLEY_REQUIRE_CONSTEXPR[m
[32m+[m[32m#endif[m
[32m+[m[32m/* Note the double-underscore. For internal use only; no API[m
[32m+[m[32m * guarantees! */[m
[32m+[m[32m#if defined(JSON_HEDLEY__IS_CONSTEXPR)[m
[32m+[m[32m    #undef JSON_HEDLEY__IS_CONSTEXPR[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if \[m
[32m+[m[32m    JSON_HEDLEY_HAS_BUILTIN(__builtin_constant_p) || \[m
[32m+[m[32m    JSON_HEDLEY_GCC_VERSION_CHECK(3,4,0) || \[m
[32m+[m[32m    JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TINYC_VERSION_CHECK(0,9,19) || \[m
[32m+[m[32m    JSON_HEDLEY_ARM_VERSION_CHECK(4,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_IBM_VERSION_CHECK(13,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_TI_VERSION_CHECK(6,1,0) || \[m
[32m+[m[32m    JSON_HEDLEY_SUNPRO_VERSION_CHECK(5,10,0) || \[m
[32m+[m[32m    JSON_HEDLEY_CRAY_VERSION_CHECK(8,1,0)[m
[32m+[m[32m    #define JSON_HEDLEY_IS_CONSTANT(expr) __builtin_constant_p(expr)[m
[32m+[m[32m#endif[m
[32m+[m[32m#if !defined(__cplusplus)[m
[32m+[m[32m#  if \[m
[32m+[m[32m       JSON_HEDLEY_HAS_BUILTIN(__builtin_types_compatible_p) || \[m
[32m+[m[32m       JSON_HEDLEY_GCC_VERSION_CHECK(3,4,0) || \[m
[32m+[m[32m       JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m       JSON_HEDLEY_IBM_VERSION_CHECK(13,1,0) || \[m
[32m+[m[32m       JSON_HEDLEY_CRAY_VERSION_CHECK(8,1,0) || \[m
[32m+[m[32m       JSON_HEDLEY_ARM_VERSION_CHECK(5,4,0) || \[m
[32m+[m[32m       JSON_HEDLEY_TINYC_VERSION_CHECK(0,9,24)[m
[32m+[m[32m#if defined(__INTPTR_TYPE__)[m
[32m+[m[32m    #define JSON_HEDLEY__IS_CONSTEXPR(expr) __builtin_types_compatible_p(__typeof__((1 ? (void*) ((__INTPTR_TYPE__) ((expr) * 0)) : (int*) 0)), int*)[m
[32m+[m[32m#else[m
[32m+[m[32m    #include <stdint.h>[m
[32m+[m[32m    #define JSON_HEDLEY__IS_CONSTEXPR(expr) __builtin_types_compatible_p(__typeof__((1 ? (void*) ((intptr_t) ((expr) * 0)) : (int*) 0)), int*)[m
[32m+[m[32m#endif[m
[32m+[m[32m#  elif \[m
[32m+[m[32m       (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) && !defined(JSON_HEDLEY_SUNPRO_VERSION) && !defined(JSON_HEDLEY_PGI_VERSION)) || \[m
[32m+[m[32m       JSON_HEDLEY_HAS_EXTENSION(c_generic_selections) || \[m
[32m+[m[32m       JSON_HEDLEY_GCC_VERSION_CHECK(4,9,0) || \[m
[32m+[m[32m       JSON_HEDLEY_INTEL_VERSION_CHECK(17,0,0) || \[m
[32m+[m[32m       JSON_HEDLEY_IBM_VERSION_CHECK(12,1,0) || \[m
[32m+[m[32m       JSON_HEDLEY_ARM_VERSION_CHECK(5,3,0)[m
[32m+[m[32m#if defined(__INTPTR_TYPE__)[m
[32m+[m[32m    #define JSON_HEDLEY__IS_CONSTEXPR(expr) _Generic((1 ? (void*) ((__INTPTR_TYPE__) ((expr) * 0)) : (int*) 0), int*: 1, void*: 0)[m
[32m+[m[32m#else[m
[32m+[m[32m    #include <stdint.h>[m
[32m+[m[32m    #define JSON_HEDLEY__IS_CONSTEXPR(expr) _Generic((1 ? (void*) ((intptr_t) * 0) : (int*) 0), int*: 1, void*: 0)[m
[32m+[m[32m#endif[m
[32m+[m[32m#  elif \[m
[32m+[m[32m       defined(JSON_HEDLEY_GCC_VERSION) || \[m
[32m+[m[32m       defined(JSON_HEDLEY_INTEL_VERSION) || \[m
[32m+[m[32m       defined(JSON_HEDLEY_TINYC_VERSION) || \[m
[32m+[m[32m       defined(JSON_HEDLEY_TI_VERSION) || \[m
[32m+[m[32m       defined(__clang__)[m
[32m+[m[32m#    define JSON_HEDLEY__IS_CONSTEXPR(expr) ( \[m
[32m+[m[32m        sizeof(void) != \[m
[32m+[m[32m        sizeof(*( \[m
[32m+[m[32m                  1 ? \[m
[32m+[m[32m                  ((void*) ((expr) * 0L) ) : \[m
[32m+[m[32m((struct { char v[sizeof(void) * 2]; } *) 1) \[m
[32m+[m[32m                ) \[m
[32m+[m[32m              ) \[m
[32m+[m[32m                                            )[m
[32m+[m[32m#  endif[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY__IS_CONSTEXPR)[m
[32m+[m[32m    #if !defined(JSON_HEDLEY_IS_CONSTANT)[m
[32m+[m[32m        #define JSON_HEDLEY_IS_CONSTANT(expr) JSON_HEDLEY__IS_CONSTEXPR(expr)[m
[32m+[m[32m    #endif[m
[32m+[m[32m    #define JSON_HEDLEY_REQUIRE_CONSTEXPR(expr) (JSON_HEDLEY__IS_CONSTEXPR(expr) ? (expr) : (-1))[m
[32m+[m[32m#else[m
[32m+[m[32m    #if !defined(JSON_HEDLEY_IS_CONSTANT)[m
[32m+[m[32m        #define JSON_HEDLEY_IS_CONSTANT(expr) (0)[m
[32m+[m[32m    #endif[m
[32m+[m[32m    #define JSON_HEDLEY_REQUIRE_CONSTEXPR(expr) (expr)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_BEGIN_C_DECLS)[m
[32m+[m[32m    #undef JSON_HEDLEY_BEGIN_C_DECLS[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_END_C_DECLS)[m
[32m+[m[32m    #undef JSON_HEDLEY_END_C_DECLS[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(JSON_HEDLEY_C_DECL)[m
[32m+[m[32m    #undef JSON_HEDLEY_C_DECL[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_BEGIN_C_DECLS extern "C" {[m
[32m+[m[32m    #define JSON_HEDLEY_END_C_DECLS }[m
[32m+[m[32m    #define JSON_HEDLEY_C_DECL extern "C"[m
[32m+[m[32m#else[m
[32m+[m[32m    #define JSON_HEDLEY_BEGIN_C_DECLS[m
[32m+[m[32m    #define JSON_HEDLEY_END_C_DECLS[m
[32m+[m[32m    #define JSON_HEDLEY_C_DECL[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_STATIC_ASSERT)[m
[32m+[m[32m    #undef JSON_HEDLEY_STATIC_ASSERT[m
[32m+[m[32m#endif[m
[32m+[m[32m#if \[m
[32m+[m[32m  !defined(__cplusplus) && ( \[m
[32m+[m[32m      (defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)) || \[m
[32m+[m[32m      JSON_HEDLEY_HAS_FEATURE(c_static_assert) || \[m
[32m+[m[32m      JSON_HEDLEY_GCC_VERSION_CHECK(6,0,0) || \[m
[32m+[m[32m      JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0) || \[m
[32m+[m[32m      defined(_Static_assert) \[m
[32m+[m[32m    )[m
[32m+[m[32m#  define JSON_HEDLEY_STATIC_ASSERT(expr, message) _Static_assert(expr, message)[m
[32m+[m[32m#elif \[m
[32m+[m[32m  (defined(__cplusplus) && (__cplusplus >= 201703L)) || \[m
[32m+[m[32m  JSON_HEDLEY_MSVC_VERSION_CHECK(16,0,0) || \[m
[32m+[m[32m  (defined(__cplusplus) && JSON_HEDLEY_TI_VERSION_CHECK(8,3,0))[m
[32m+[m[32m#  define JSON_HEDLEY_STATIC_ASSERT(expr, message) static_assert(expr, message)[m
[32m+[m[32m#elif defined(__cplusplus) && (__cplusplus >= 201103L)[m
[32m+[m[32m#  define JSON_HEDLEY_STATIC_ASSERT(expr, message) static_assert(expr)[m
[32m+[m[32m#else[m
[32m+[m[32m#  define JSON_HEDLEY_STATIC_ASSERT(expr, message)[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_CONST_CAST)[m
[32m+[m[32m    #undef JSON_HEDLEY_CONST_CAST[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus)[m
[32m+[m[32m#  define JSON_HEDLEY_CONST_CAST(T, expr) (const_cast<T>(expr))[m
[32m+[m[32m#elif \[m
[32m+[m[32m  JSON_HEDLEY_HAS_WARNING("-Wcast-qual") || \[m
[32m+[m[32m  JSON_HEDLEY_GCC_VERSION_CHECK(4,6,0) || \[m
[32m+[m[32m  JSON_HEDLEY_INTEL_VERSION_CHECK(13,0,0)[m
[32m+[m[32m#  define JSON_HEDLEY_CONST_CAST(T, expr) (__extension__ ({ \[m
[32m+[m[32m        JSON_HEDLEY_DIAGNOSTIC_PUSH \[m
[32m+[m[32m        JSON_HEDLEY_DIAGNOSTIC_DISABLE_CAST_QUAL \[m
[32m+[m[32m        ((T) (expr)); \[m
[32m+[m[32m        JSON_HEDLEY_DIAGNOSTIC_POP \[m
[32m+[m[32m    }))[m
[32m+[m[32m#else[m
[32m+[m[32m#  define JSON_HEDLEY_CONST_CAST(T, expr) ((T) (expr))[m
[32m+[m[32m#endif[m
[32m+[m
[32m+[m[32m#if defined(JSON_HEDLEY_REINTERPRET_CAST)[m
[32m+[m[32m    #undef JSON_HEDLEY_REINTERPRET_CAST[m
[32m+[m[32m#endif[m
[32m+[m[32m#if defined(__cplusplus)[m
[32m+[m[32m    #define JSON_HEDLEY_REINTERP