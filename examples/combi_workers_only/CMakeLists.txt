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
target_link_libraries(combi_example_worker_only discotec Boost::boost)

install(TARGETS combi_example_worker_only DESTINATION examples/combi_example_worker_only)
install(FILES ctparam run.sh DESTINATION examples/combi_example_worker_only)