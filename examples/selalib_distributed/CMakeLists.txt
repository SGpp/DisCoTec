# TODO# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.18)

project("DisCoTec distributed with selalib"
        LANGUAGES C CXX Fortran
        DESCRIPTION "Run SeLaLib simulations with the combination technique using DisCoTec")

if (NOT TARGET discotec)
    add_subdirectory(../../src discotec)
endif ()
find_package(Boost REQUIRED COMPONENTS filesystem program_options system serialization)

find_package(SELALIB REQUIRED)
set( SELALIB_INTERFACE_LIBRARY_NAME "libsll_m_sim_bsl_vp_3d3v_cart_dd_slim_interface" )

if(NOT EXISTS ${SELALIB_DIR}/lib/${SELALIB_INTERFACE_LIBRARY_NAME}.a)
    message(FATAL_ERROR "SELALIB: ${SELALIB_DIR}/lib/${SELALIB_INTERFACE_LIBRARY_NAME}.a not found")
endif()

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(combine_selalib_diagnostics src/combine_selalib_diagnostics.cpp)
add_executable(selalib_distributed src/selalib_distributed.cpp)
add_executable(selalib_distributed_workers_only src/selalib_distributed_workers_only.cpp)

set( HDF5_USE_STATIC_LIBRARIES ON )
find_package( HDF5 COMPONENTS C Fortran REQUIRED )

target_compile_features(combine_selalib_diagnostics PRIVATE cxx_std_17)
target_link_libraries(combine_selalib_diagnostics discotec ${Boost_LIBRARIES}  ${SELALIB_DIR}/lib/${SELALIB_INTERFACE_LIBRARY_NAME}.a)

target_compile_features(selalib_distributed PRIVATE cxx_std_17)
set_target_properties(selalib_distributed PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(selalib_distributed discotec ${Boost_LIBRARIES}  ${SELALIB_DIR}/lib/${SELALIB_INTERFACE_LIBRARY_NAME}.a  ${SELALIB_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES} ${OpenMP_CXX_LIBRARIES} ${OpenMP_C_LIBRARIES})

target_compile_features(selalib_distributed_workers_only PRIVATE cxx_std_17)
set_target_properties(selalib_distributed_workers_only PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(selalib_distributed_workers_only discotec ${Boost_LIBRARIES}  ${SELALIB_DIR}/lib/${SELALIB_INTERFACE_LIBRARY_NAME}.a  ${SELALIB_LIBRARIES} ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES} ${OpenMP_CXX_LIBRARIES} ${OpenMP_C_LIBRARIES})
