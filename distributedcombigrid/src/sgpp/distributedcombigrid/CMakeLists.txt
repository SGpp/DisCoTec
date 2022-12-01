# Copyright (C) 2008-today The SG++ Project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

cmake_minimum_required(VERSION 3.24)

project("DisCoTec"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "TODO")


option(DISCOTEC_ENABLEFT "DisCoTec with algorithm-based fault tolerance" OFF)
option(DISCOTEC_GENE "Build for GENE (as task library)" OFF)
option(DISCOTEC_TIMING "Build with high-resolution timing" ON)
option(DISCOTEC_TEST "Build with Bost tests" ON) #TODO: sinvoller Bescreiben oder in Boost auslagern
option(DISCOTEC_OPENMP "Parallelize with OpenMP" OFF)
option(DISCOTEC_HDF5 "Interpolation output with HDF5" ON)
option(DISCOTEC_DEBUG "Build with assertions and debug symbols" OFF) #TODO: Move to build type


#TODO: find or fetch highfive if hdf5 is enabled ??
#TODO: find or fetch GENE if enabled
#TODO: find or fetch hdf5 if is enabled

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

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

include(ExternalProject)

ExternalProject_Add( #TODO: first search for glpk
    glpk
    URL http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz # TODO: use newer and not hardcoded Version
    DOWNLOAD_EXTRACT_TIMESTAMP false
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
    BUILD_COMMAND make -j htopinstall # TODO: get thread count and build parallel
    )
add_dependencies(discotec glpk)


ExternalProject_Get_property(glpk INSTALL_DIR) # TODO: better
target_include_directories(discotec PUBLIC ${INSTALL_DIR}/include)
target_link_libraries(discotec PUBLIC ${INSTALL_DIR}/lib/libglpk.a)
find_package(Boost REQUIRED COMPONENTS serialization)
target_link_libraries(discotec PUBLIC ${Boost_LIBRARIES})
