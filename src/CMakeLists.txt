# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "DisCoTec is a code for running the distributed sparse grid combination technique with MPI parallelization.")


set(DISCOTEC_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/combischeme/CombiMinMaxScheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/combischeme/CombiThirdLevelScheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/FaultCriterion.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/FTUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/LPOptimizationInterpolation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/StaticFaults.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/WeibullFaults.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/io/H5InputOutput.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/loadmodel/AverageOfLastNLoadModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/loadmodel/AveragingLoadModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/loadmodel/LinearLoadModel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/manager/ProcessGroupManager.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/manager/ProcessGroupWorker.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/manager/ProcessManager.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi/MPIMemory.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi/MPISystem.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/MPI-FT_bitflips.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/MPI-FT_ULFM.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/MPI-FT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_abort.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_allreduce.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_barrier.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_bcast.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_cart_coords.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_cart_create.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_cart_rank.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_create.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_free.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_group.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_rank.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_set_errhandler.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_size.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_comm_split.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_finalize.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_get_count.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_iallreduce.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_ibarrier.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_ibcast.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_init.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_irecv.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_isend.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_probe.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_recv.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_reduce.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_request_free.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_send.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_wait.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mpi_fault_simulator/Sim_FT_waitall.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rescheduler/RebalancingTaskRescheduler.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/task/Task.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_level/NetworkUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/third_level/ThirdLevelUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/LevelSetUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/LevelVector.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/MonteCarlo.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/utils/Stats.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/vtk/PlotFileWriter.cpp
        )

add_library(discotec ${DISCOTEC_SOURCES})
target_compile_features(discotec PRIVATE cxx_std_17)
target_include_directories(discotec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Defining options
include(CMakeDependentOption)

option(DISCOTEC_BUILD_MISSING_DEPS "Automatic fetch and build first layer of missing dependencies" ON)
option(DISCOTEC_TIMING "Build with high-resolution timing" ON)
if (DISCOTEC_TIMING)
    target_compile_definitions(discotec PUBLIC TIMING)
endif ()

option(DISCOTEC_USE_HIGHFIVE "Interpolation output with HighFive/HDF5" ON)
option(DISCOTEC_UNIFORMDECOMPOSITION "Use uniform decomposition" ON) # TODO: @polinta: does not compile if off
if (DISCOTEC_UNIFORMDECOMPOSITION)
    target_compile_definitions(discotec PUBLIC UNIFORMDECOMPOSITION) # has to be PUBLIC for tests, rename to DISCOTEC_UNIFORMDECOMPOSITION?
endif ()
option(DISCOTEC_GENE "Build for GENE (as task library)" OFF) #TODO: handle
if (DISCOTEC_GENE)
    message(FATAL_ERROR "Building with GENE is not supported yet")
endif ()
option(DISCOTEC_OPENMP "Parallelize with OpenMP" OFF)
if (DISCOTEC_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(discotec PUBLIC OpenMP::OpenMP_CXX)
endif ()
option(DISCOTEC_ENABLEFT "DisCoTec with algorithm-based fault tolerance" OFF)
if (DISCOTEC_ENABLEFT)
    add_compile_definitions(ENABLEFT)
endif ()
option(DISCOTEC_USE_VTK "Use VTK for visualization" OFF)
option(DISCOTEC_USE_LTO "Use link-time optimization" ON)
if (DISCOTEC_USE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output)
    if (result)
        set_property(TARGET discotec PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else ()
        message(FATAL "IPO is not supported: ${output}")
    endif ()
endif ()
#TODO: find or fetch GENE if enabled

option(DISCOTEC_OMITREADYSIGNAL "Omit ready signal" OFF)
if (DISCOTEC_OMITREADYSIGNAL)
    target_compile_definitions(discotec PRIVATE OMITREADYSIGNAL)
endif ()

option(DISCOTEC_USENONBLOCKINGMPICOLLECTIVE "Use non-blocking MPI collective operations" OFF)
if (DISCOTEC_USENONBLOCKINGMPICOLLECTIVE)
    target_compile_definitions(discotec PRIVATE USENONBLOCKINGMPICOLLECTIVE)
endif ()

#ISGENE #TODO: handle if access to GENE

# Handle dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)
find_package(MPI REQUIRED)
target_link_libraries(discotec PUBLIC MPI::MPI_CXX)

include(FetchContent)

find_package(Boost REQUIRED COMPONENTS serialization headers program_options unit_test_framework system filesystem)

target_include_directories(discotec PRIVATE Boost::headers)
target_link_libraries(discotec PRIVATE Boost::serialization Boost::program_options Boost::system )

find_package(GLPK REQUIRED)
target_link_libraries(discotec PUBLIC ${GLPK_LIBRARIES})


if (DISCOTEC_USE_HIGHFIVE)
    enable_language(C)
    find_package(HighFive REQUIRED)  # TODO: add minimum version and components
    target_include_directories(discotec PRIVATE ${HighFive_INCLUDE_DIRS})
    target_link_libraries(discotec PRIVATE HighFive)
    target_compile_definitions(discotec PRIVATE HAVE_HIGHFIVE)

endif ()


if (DISCOTEC_USE_VTK)
    enable_language(C)
    find_package(VTK REQUIRED) # TODO: add minimum version and components

    target_link_libraries(discotec PRIVATE  VTK::CommonCore VTK::CommonDataModel VTK::IOParallelXML)

    target_include_directories(discotec PRIVATE ${VTK_INCLUDE_DIRS})
    target_compile_definitions(discotec PRIVATE USE_VTK)
endif ()
install(TARGETS discotec
        EXPORT discotecTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
        )