cmake_minimum_required(VERSION 3.24.2)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(errorCalc errorCalc.cpp)
target_compile_features(errorCalc PRIVATE cxx_std_17)

find_package(Boost REQUIRED)
target_link_libraries(errorCalc discotec Boost::boost)

add_subdirectory(subspace_writer)
