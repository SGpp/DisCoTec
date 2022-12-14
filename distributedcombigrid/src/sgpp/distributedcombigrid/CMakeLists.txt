# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.24)

project("DisCoTec"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "DisCoTec is a code for running the distributed sparse grid combination technique with MPI parallelization.")


#TODO: find or fetch GENE if enabled

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(DISCOTEC_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/combicom/CombiCom.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/combischeme/CombiMinMaxScheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/combischeme/CombiThirdLevelScheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/FaultCriterion.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/FTUtils.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/LPOptimizationInterpolation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/StaticFaults.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fault_tolerance/WeibullFaults.cpp
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
target_include_directories(discotec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../) #TODO: fix paths
find_package(MPI REQUIRED)
target_link_libraries(discotec PUBLIC MPI::MPI_CXX)

include(FetchContent)

find_package(Boost COMPONENTS serialization headers) #TODO: add all required components
if (Boost_FOUND)
    message(STATUS "Boost found")
else ()
    message(WARNING "Boost not found. Building own Boost!")
    FetchContent_Declare(
            Boost
            GIT_REPOSITORY https://github.com/boostorg/boost.git
            GIT_TAG boost-1.80.0
            GIT_SHALLOW TRUE
            QUIET
    )
    FetchContent_MakeAvailable(Boost)
endif ()

target_include_directories(discotec PUBLIC Boost::headers)
target_link_libraries(discotec PRIVATE Boost::serialization) #TODO: add all required components


FetchContent_Declare(
        glpk
        URL http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz # TODO: use newer and not hardcoded Version
)
FetchContent_MakeAvailable(glpk)
execute_process(COMMAND ${glpk_SOURCE_DIR}/configure --prefix=${glpk_BINARY_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND make -j ${nproc} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND make install OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${glpk_BINARY_DIR})
target_include_directories(discotec PUBLIC ${glpk_BINARY_DIR}/include)
target_link_libraries(discotec PUBLIC ${glpk_BINARY_DIR}/lib/libglpk.a)


option(DISCOTEC_TIMING "Build with high-resolution timing" ON)
if (DISCOTEC_TIMING)
    target_compile_definitions(discotec PUBLIC TIMING)
endif ()

option(DISCOTEC_HDF5 "Interpolation output with HDF5" ON) # ON
include(CMakeDependentOption)
cmake_dependent_option(DISCOTEC_USE_HIGHFIVE "Use HighFive for HDF5" ON DISCOTEC_HDF5 OFF)
if (DISCOTEC_HDF5)
    #TODO: find or fetch highfive if hdf5 is enabled ??
    find_package(HDF5 COMPONENTS CXX HINTS HDF5_PREFER_PARALLEL=True)
    if (HDF5_FOUND AND HDF5_IS_PARALLEL)
        message(STATUS "Found HDF5: ${HDF5_CXX_LIBRARIES}")
    else ()
        message(WARNING "Parallel HDF5 not found, building own HDF5!")
        FetchContent_Declare(
                HDF5
                GIT_REPOSITORY git@github.com:HDFGroup/hdf5.git
                GIT_TAG hdf5-1_12_0
                GIT_SHALLOW TRUE
                QUIET
        )
        set(HDF5_ENABLE_PARALLEL ON CACHE BOOL "" FORCE)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_FRAMEWORKS OFF CACHE BOOL "" FORCE)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
        set(TEST_SHELL_SCRIPTS OFF CACHE BOOL "" FORCE)
        set(HDF5_TEST_TOOLS OFF CACHE BOOL "" FORCE)
        set(HDF5_TEST_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(HDF5_TEST_SWMR OFF CACHE BOOL "" FORCE)
        set(HDF5_TESTS_PARALLEL OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(HDF5)
    endif ()
    target_link_libraries(discotec PUBLIC ${HDF5_CXX_LIBRARIES})

    if (DISCOTEC_USE_HIGHFIVE)
        find_package(HighFive QUIET)
        if (HighFive_FOUND)
            message(STATUS "Found HighFive: ${HighFive_INCLUDE_DIRS}")
        else ()
            message(WARNING "HighFive not found, building own HighFive!")
            FetchContent_Declare(
                    HighFive
                    GIT_REPOSITORY git@github.com:BlueBrain/HighFive.git
                    GIT_TAG v2.6.2
                    GIT_SHALLOW TRUE
                    QUIET
            )
            set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" FORCE)
            set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" FORCE)
            set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" FORCE)
            set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" FORCE)
            FetchContent_MakeAvailable(HighFive)
        endif ()
        target_include_directories(discotec PUBLIC ${HighFive_INCLUDE_DIRS})
        target_link_libraries(discotec PUBLIC HighFive)
    endif ()
endif ()


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


option(DISCOTEC_ENABLEFT "DisCoTec with algorithm-based fault tolerance" ON)
if (DISCOTEC_ENABLEFT)
    add_compile_definitions(ENABLEFT)
endif ()

option(DISCOTEC_DEBUG "Build with assertions and debug symbols" OFF) #TODO: Move to build type
if (DISCOTEC_DEBUG OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_definitions(discotec PUBLIC DEBUG_OUTPUT)
endif ()

option(DISCOTEC_USE_VTK "Use VTK for visualization" ON)
if (DISCOTEC_USE_VTK)
    enable_language(C)
    find_package(VTK) # TODO: add minimum version and components
    if (VTK_FOUND)
        message(STATUS "Found VTK: ${VTK_VERSION}")
    else ()
        message(WARNING "VTK not found, building own VTK!")
        FetchContent_Declare(
                VTK
                GIT_REPOSITORY git@github.com:Kitware/VTK.git
                GIT_TAG v9.0.1
                GIT_SHALLOW TRUE
        )
        add_library(hdf5::hdf5 INTERFACE IMPORTED)
        set(VTK_USE_MPI ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(HDF5 VTK)
    endif ()


    target_link_libraries(discotec PRIVATE ${VTK_LIBRARIES})
    target_include_directories(discotec PRIVATE ${VTK_INCLUDE_DIRS})
    target_compile_definitions(discotec PRIVATE USE_VTK)
endif ()

option(DISCOTEC_OMITREADYSIGNAL "Omit ready signal" OFF)
if (DISCOTEC_OMITREADYSIGNAL)
    target_compile_definitions(discotec PRIVATE OMITREADYSIGNAL)
endif ()

option(DISCOTEC_USENONBLOCKINGMPICOLLECTIVE "Use non-blocking MPI collective operations" OFF)
if (DISCOTEC_USENONBLOCKINGMPICOLLECTIVE)
    target_compile_definitions(discotec PRIVATE USENONBLOCKINGMPICOLLECTIVE)
endif ()

#USE_NON_BLOCKING_COLLECTIVES #wird nirgends in SCONS definiert also ka wann definieren
#ALT_LEVEL_VECTOR #wird nirgends in SCONS definiert also ka wann definieren
#HAVE_HIGHFIVE #wird nirgends in SCONS definiert also ka wann definieren
#SIMFTMPI #TODO: @polinta: what is this? in confg.hpp the same is set no matter if defined or not ...
#ISGENE #TODO: handle if access to GENE