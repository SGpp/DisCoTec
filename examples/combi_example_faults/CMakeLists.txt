# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed combi faults example"
        LANGUAGES CXX
        DESCRIPTION "TODO")

if (NOT TARGET discotec)
    add_subdirectory(../../src discotec)
endif ()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
find_package(Boost REQUIRED)

add_executable(combi_example_faults combi_example_faults.cpp)
target_compile_features(combi_example_faults PRIVATE cxx_std_17)
target_link_libraries(combi_example_faults discotec Boost::boost)

install(TARGETS combi_example_faults DESTINATION examples/combi_example_faults)
install(FILES ctparam run.sh DESTINATION examples/combi_example_faults)