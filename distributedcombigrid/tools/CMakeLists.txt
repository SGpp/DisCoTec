cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(errorCalc errorCalc.cpp)
target_link_libraries(errorCalc discotec)
