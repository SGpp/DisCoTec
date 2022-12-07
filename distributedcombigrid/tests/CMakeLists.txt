cmake_minimum_required(VERSION 3.5)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")


enable_testing()
find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})


file(GLOB UNIT_TESTS_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(test_distributedcombigrid_boost ${UNIT_TESTS_SRC_FILES})
add_dependencies(test_distributedcombigrid_boost glpk)

target_link_libraries(test_distributedcombigrid_boost ${Boost_LIBRARIES} discotec)
target_include_directories(test_distributedcombigrid_boost PUBLIC ${Boost_INCLUDE_DIRS})
set_property(TARGET test_distributedcombigrid_boost PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9')

# TODO: automate test split up and split up further
add_test(NAME distributedfullgrid COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=distributedfullgrid)
add_test(NAME distributedsparsegrid COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=distributedsparsegrid)
add_test(NAME ftolerance COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=ftolerance)
add_test(NAME integration COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=integration)
add_test(NAME fullgrid COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=fullgrid)
add_test(NAME hierarchization COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=hierarchization)
add_test(NAME loadmodel COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=loadmodel)
add_test(NAME reduce COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=reduce)
add_test(NAME rescheduling COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=rescheduling)
add_test(NAME stats COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=stats)
add_test(NAME task COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=task)
add_test(NAME mpisystem COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=mpisystem)
add_test(NAME networkutils COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=networkutils)
add_test(NAME thirdlevel COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 9 ${CMAKE_CURRENT_SOURCE_DIR}/test_distributedcombigrid_boost --run_test=thirdlevel)
