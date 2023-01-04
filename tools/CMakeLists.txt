cmake_minimum_required(VERSION 3.24.2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../distributedcombigrid/tools)

add_executable(errorCalc errorCalc.cpp)

find_package(Boost REQUIRED)
target_link_libraries(errorCalc discotec Boost::boost)
