cmake_minimum_required(VERSION 3.5)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")


enable_testing()
find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB UNIT_TESTS_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(test_distributedcombigrid_boost ${UNIT_TESTS_SRC_FILES})
target_link_libraries(test_distributedcombigrid_boost ${Boost_LIBRARIES} discotec)
target_include_directories(test_distributedcombigrid_boost PUBLIC ${Boost_INCLUDE_DIRS})
enable_testing()
add_test(DemoTest test_distributedcombigrid_boost)