# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.24)

project("DisCoTec"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "TODO")

#TODO: find or fetch highfive if hdf5 is enabled ??
#TODO: find or fetch GENE if enabled
#TODO: find or fetch hdf5 if is enabled


set (DISCOTEC_SOURCES
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
target_include_directories(discotec PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../) #TODO: fix paths
find_package(MPI REQUIRED)
target_link_libraries(discotec PUBLIC MPI::MPI_CXX)

find_package(Boost REQUIRED COMPONENTS serialization headers) #TODO: add all required components
target_include_directories(discotec PUBLIC Boost::headers)
target_link_libraries(discotec PRIVATE Boost::serialization)
include(FetchContent)


if(DISCOTEC_TIMING)
    target_compile_definitions(discotec PUBLIC TIMING)
endif()

if(DISCOTEC_HDF5)
    find_package(HDF5 COMPONENTS CXX HINTS HDF5_PREFER_PARALLEL=True)
    if(HDF5_FOUND AND HDF5_IS_PARALLEL)
        message(STATUS "Found HDF5: ${HDF5_CXX_LIBRARIES}")
        #TODO: handle
    else()
        message(WARNING "HDF5 not found, building own HDF5!")
        FetchContent_Declare(
                HDF5
                HDF5_ENABLE_PARALLEL=ON
                CMAKE_BUILD_TYPE:STRING=Release
                BUILD_SHARED_LIBS:BOOL=OFF
                BUILD_TESTING:BOOL=OFF
                DHDF5_BUILD_TOOLS:BOOL=OFF
                TEST_SHELL_SCRIPTS:BOOL=OFF
                GIT_REPOSITORY git@github.com:HDFGroup/hdf5.git
                GIT_TAG        hdf5-1_12_0)
        FetchContent_MakeAvailable(HDF5)
    endif()
    target_link_libraries(discotec PUBLIC ${HDF5_CXX_LIBRARIES})
endif()

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
