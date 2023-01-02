cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec main project"
        LANGUAGES CXX
        DESCRIPTION "Main CMake for DiscoTec. Use this to configure the complete project.")

add_subdirectory(src/sgpp/distributedcombigrid/)
add_subdirectory(examples)
add_subdirectory(tools)

option(DISCOTEC_TEST "Build with Boost tests" ON) #TODO: more sensible description
if(DISCOTEC_TEST)
    enable_testing()
    add_subdirectory(tests)
endif()