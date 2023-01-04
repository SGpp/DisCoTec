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
        combischeme/CombiMinMaxScheme.cpp
        combischeme/CombiThirdLevelScheme.cpp
        fault_tolerance/FaultCriterion.cpp
        fault_tolerance/FTUtils.cpp
        fault_tolerance/LPOptimizationInterpolation.cpp
        fault_tolerance/StaticFaults.cpp
        fault_tolerance/WeibullFaults.cpp
        loadmodel/AverageOfLastNLoadModel.cpp
        loadmodel/AveragingLoadModel.cpp
        loadmodel/LinearLoadModel.cpp
        manager/ProcessGroupManager.cpp
        manager/ProcessGroupWorker.cpp
        manager/ProcessManager.cpp
        mpi/MPIMemory.cpp
        mpi/MPISystem.cpp
        mpi_fault_simulator/MPI-FT_bitflips.cpp
        mpi_fault_simulator/MPI-FT_ULFM.cpp
        mpi_fault_simulator/MPI-FT.cpp
        mpi_fault_simulator/Sim_FT_abort.cpp
        mpi_fault_simulator/Sim_FT_allreduce.cpp
        mpi_fault_simulator/Sim_FT_barrier.cpp
        mpi_fault_simulator/Sim_FT_bcast.cpp
        mpi_fault_simulator/Sim_FT_cart_coords.cpp
        mpi_fault_simulator/Sim_FT_cart_create.cpp
        mpi_fault_simulator/Sim_FT_cart_rank.cpp
        mpi_fault_simulator/Sim_FT_comm_create.cpp
        mpi_fault_simulator/Sim_FT_comm_free.cpp
        mpi_fault_simulator/Sim_FT_comm_group.cpp
        mpi_fault_simulator/Sim_FT_comm_rank.cpp
        mpi_fault_simulator/Sim_FT_comm_set_errhandler.cpp
        mpi_fault_simulator/Sim_FT_comm_size.cpp
        mpi_fault_simulator/Sim_FT_comm_split.cpp
        mpi_fault_simulator/Sim_FT_finalize.cpp
        mpi_fault_simulator/Sim_FT_get_count.cpp
        mpi_fault_simulator/Sim_FT_iallreduce.cpp
        mpi_fault_simulator/Sim_FT_ibarrier.cpp
        mpi_fault_simulator/Sim_FT_ibcast.cpp
        mpi_fault_simulator/Sim_FT_init.cpp
        mpi_fault_simulator/Sim_FT_irecv.cpp
        mpi_fault_simulator/Sim_FT_isend.cpp
        mpi_fault_simulator/Sim_FT_probe.cpp
        mpi_fault_simulator/Sim_FT_recv.cpp
        mpi_fault_simulator/Sim_FT_reduce.cpp
        mpi_fault_simulator/Sim_FT_request_free.cpp
        mpi_fault_simulator/Sim_FT_send.cpp
        mpi_fault_simulator/Sim_FT_test.cpp
        mpi_fault_simulator/Sim_FT_wait.cpp
        mpi_fault_simulator/Sim_FT_waitall.cpp
        rescheduler/RebalancingTaskRescheduler.cpp
        task/Task.cpp
        third_level/NetworkUtils.cpp
        third_level/ThirdLevelUtils.cpp
        utils/LevelSetUtils.cpp
        utils/LevelVector.cpp
        utils/MonteCarlo.cpp
        utils/Stats.cpp
        vtk/PlotFileWriter.cpp
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
option(DISCOTEC_HDF5 "Interpolation output with HDF5" ON) # ON
cmake_dependent_option(DISCOTEC_USE_HIGHFIVE "Use HighFive for HDF5" ON DISCOTEC_HDF5 OFF)
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
option(DISCOTEC_DEBUG "Build with assertions and debug symbols" OFF) #TODO: Move to build type
if (DISCOTEC_DEBUG OR CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_definitions(discotec PUBLIC DEBUG_OUTPUT)
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
find_package(MPI REQUIRED)
target_link_libraries(discotec PUBLIC MPI::MPI_CXX)

include(FetchContent)

find_package(Boost QUIET COMPONENTS serialization headers program_options unit_test_framework system filesystem)
if (Boost_FOUND)
    message(STATUS "Boost found")
else ()
    if (NOT DISCOTEC_BUILD_MISSING_DEPS)
        message(FATAL_ERROR "Boost not found")
    endif ()
    message(WARNING "Boost not found. Building own Boost!")
    list(APPEND CMAKE_MESSAGE_INDENT "Boost: ")
    set(Boost_VERSION 1.80.0 CACHE STRING "Boost version to build")
    set(Boost_URL https://boostorg.jfrog.io/artifactory/main/release/${Boost_VERSION}/source/boost_${Boost_VERSION//./_}.tar.gz)
    FetchContent_Declare(
            Boost_binaries
            URL ${Boost_URL}
            OVERRIDE_FIND_PACKAGE
            )
    FetchContent_MakeAvailable(Boost_binaries)
    message(STATUS "Cloning Boost")
    FetchContent_GetProperties(Boost_binaries)
    message(STATUS "Bootstraping Boost")
    execute_process(COMMAND ./bootstrap.sh --prefix=${boost_binaries_BINARY_DIR} --with-libraries=serialization,program_options,system,math,test,headers WORKING_DIRECTORY  ${boost_binaries_SOURCE_DIR}  OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    message(STATUS "Building Boost")
    execute_process(COMMAND ./b2 WORKING_DIRECTORY  ${boost_binaries_SOURCE_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND ./b2 headers WORKING_DIRECTORY  ${boost_binaries_SOURCE_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
    execute_process(COMMAND ./b2 install WORKING_DIRECTORY  ${boost_binaries_SOURCE_DIR} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)


    set(BOOST_ROOT ${boost_binaries_BINARY_DIR})
    set(Boost_NO_SYSTEM_PATHS ON)
    find_package(Boost ${Boost_VERSION} REQUIRED COMPONENTS serialization program_options system unit_test_framework)

    list(POP_BACK CMAKE_MESSAGE_INDENT)
endif ()

target_include_directories(discotec PRIVATE Boost::headers)
target_link_libraries(discotec PRIVATE Boost::serialization Boost::program_options Boost::system )

message(STATUS "Fetching glpk")
FetchContent_Declare(
        glpk
        URL http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz # TODO: use newer and not hardcoded Version
        QUIET
)
FetchContent_MakeAvailable(glpk)
message(STATUS "GLPK: configure ...")
execute_process(COMMAND ${glpk_SOURCE_DIR}/configure --prefix=${glpk_BINARY_DIR} OUTPUT_QUIET ERROR_QUIET COMMAND_ERROR_IS_FATAL ANY)
message(STATUS "GLPK: build ...")
execute_process(COMMAND make -j ${nproc} OUTPUT_QUIET ERROR_QUIET COMMAND_ERROR_IS_FATAL ANY)
message(STATUS "GLPK: install ...")
execute_process(COMMAND make install OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${glpk_BINARY_DIR})
target_include_directories(discotec PUBLIC ${glpk_BINARY_DIR}/include)
target_link_libraries(discotec PUBLIC ${glpk_BINARY_DIR}/lib/libglpk.a)


if (DISCOTEC_HDF5 OR DISCOTEC_USE_VTK)
    enable_language(C)

    set(HDF5_PREFER_PARALLEL ON)
    find_package(HDF5) # TODO: add minimum version and components


    if (HDF5_FOUND AND HDF5_IS_PARALLEL)
        message(STATUS "Found HDF5")
    else ()
        if (NOT DISCOTEC_BUILD_MISSING_DEPS)
            message(FATAL_ERROR "HDF5 not found")
        endif ()
        list(APPEND CMAKE_MESSAGE_INDENT "HDF5: ")

        message(WARNING "Parallel HDF5 not found, building own HDF5!")
        FetchContent_Declare(
                HDF5
                GIT_REPOSITORY git@github.com:HDFGroup/hdf5.git
                GIT_TAG hdf5-1_13_0
                GIT_SHALLOW TRUE
                QUIET
                OVERRIDE_FIND_PACKAGE
        )
        set(HDF5_ENABLE_PARALLEL ON CACHE BOOL "" FORCE)
        set(HDF5_IS_PARALLEL ON CACHE BOOL "" FORCE)
        set(HDF5_EXTERNALLY_CONFIGURED ON CACHE BOOL "" FORCE)
        set(HDF5_EXTERNALLY_CONFIGURED ON CACHE BOOL "" FORCE)
        set(HDF5_BUILD_PARALLEL_TOOLS OFF CACHE BOOL "" FORCE)
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
        set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "" FORCE) # TODO: enable fortran when using selalib
        set(HDF5_ENABLE_THREADSAFE OFF CACHE BOOL "" FORCE)
        set(HDF5_BUILD_JAVA OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(HDF5)

        if (IS_DIRECTORY "${HDF5_SOURCE_DIR}")
            set_property(DIRECTORY ${HDF5_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
        endif ()
        list(POP_BACK CMAKE_MESSAGE_INDENT)

        target_include_directories(discotec PRIVATE ${HDF5_SOURCE_DIR}/src) # TODO
        target_include_directories(discotec PRIVATE ${HDF5_BINARY_DIR}/src) # TODO:

        target_link_libraries(discotec PUBLIC hdf5-static)

        if("${HDF5_LIBRARIES}" MATCHES ".*-NOTFOUND")
         set(HDF5_LIBRARIES "")
        endif()
        if ("${HDF5_INCLUDE_DIRS}" MATCHES ".*-NOTFOUND")
         set(HDF5_INCLUDE_DIRS "")
        endif()
    endif ()


    if (DISCOTEC_USE_HIGHFIVE)
        find_package(HighFive QUIET)  # TODO: add minimum version and components
        if (HighFive_FOUND)
            message(STATUS "Found HighFive: ${HighFive_INCLUDE_DIRS}")
        else ()
            if (NOT DISCOTEC_BUILD_MISSING_DEPS)
                message(FATAL_ERROR "HighFive not found")
            endif ()
            message(WARNING "HighFive not found, building own HighFive!")
            list(APPEND CMAKE_MESSAGE_INDENT "HighFive: ")
            FetchContent_Declare(
                    HighFive
                    GIT_REPOSITORY git@github.com:BlueBrain/HighFive.git
                    GIT_TAG v2.6.2
                    GIT_SHALLOW TRUE
                    QUIET
                    OVERRIDE_FIND_PACKAGE
            )
            set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "" FORCE)
            set(HIGHFIVE_PARALLEL_HDF5 ON CACHE BOOL "" FORCE)
            set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "" FORCE)
            set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "" FORCE)
            set(HIGHFIVE_USE_BOOST OFF CACHE BOOL "" FORCE)

            FetchContent_MakeAvailable(HighFive)
            if (IS_DIRECTORY "${HighFive_SOURCE_DIR}")
                set_property(DIRECTORY ${HighFive_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
            endif ()
        endif ()
        target_include_directories(discotec PRIVATE ${HighFive_INCLUDE_DIRS})
        target_link_libraries(discotec PRIVATE HighFive)
        target_compile_definitions(discotec PRIVATE HAVE_HIGHFIVE)
        list(POP_BACK CMAKE_MESSAGE_INDENT)
    endif ()
endif ()


if (DISCOTEC_USE_VTK)
    enable_language(C)
    find_package(VTK) # TODO: add minimum version and components
    if (VTK_FOUND)
        message(STATUS "Found VTK: ${VTK_VERSION}")
    else ()
        if (NOT DISCOTEC_BUILD_MISSING_DEPS)
            message(FATAL_ERROR "VTK not found")
        endif ()
        message(WARNING "VTK not found, building own VTK!")
        list(APPEND CMAKE_MESSAGE_INDENT "VTK: ")
        FetchContent_Declare(
                VTK
                GIT_REPOSITORY git@github.com:Kitware/VTK.git
                GIT_TAG v9.1.0
                GIT_SHALLOW TRUE
        )
        set(VTK_USE_MPI ON CACHE BOOL "" FORCE)
        set(VTK_MODULE_USE_EXTERNAL_VTK_hdf5 ON CACHE BOOL "" FORCE)

        FetchContent_MakeAvailable(VTK)

        if (IS_DIRECTORY "${VTK_SOURCE_DIR}")
            set_property(DIRECTORY ${VTK_SOURCE_DIR} PROPERTY EXCLUDE_FROM_ALL YES)
        endif ()
    endif ()




    target_link_libraries(discotec PRIVATE  VTK::CommonCore VTK::CommonDataModel VTK::IOParallelXML)

    target_include_directories(discotec PRIVATE ${VTK_INCLUDE_DIRS})
    target_compile_definitions(discotec PRIVATE USE_VTK)
    list(POP_BACK CMAKE_MESSAGE_INDENT)
endif ()