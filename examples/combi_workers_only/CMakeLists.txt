# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed combi example without manager rank"
        LANGUAGES CXX
        DESCRIPTION "TODO")

if (NOT TARGET discotec)
    add_subdirectory(../../src discotec)
endif ()
find_package(Boost REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(combi_example_worker_only combi_example_worker_only.cpp)
target_compile_features(combi_example_worker_only PRIVATE cxx_std_17)
if(DISCOTEC_USE_HIGHFIVE)
    enable_language(C)
    find_package(HDF5 REQUIRED)
    find_package(HighFive REQUIRED)
    target_include_directories(discotec PRIVATE ${HighFive_INCLUDE_DIRS})
    target_link_libraries(combi_example_worker_only PUBLIC HighFive discotec Boost::boost)
    target_compile_definitions(combi_example_worker_only PUBLIC DISCOTEC_USE_HIGHFIVE)
else()
target_link_libraries(combi_example_worker_only discotec Boost::boost)

endif()

install(TARGETS combi_example_worker_only DESTINATION examples/combi_example_worker_only)
install(FILES ctparam run.sh DESTINATION examples/combi_example_worker_only)