cmake_minimum_required(VERSION 3.18)

project("DisCoTec main project"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "TODO")

add_subdirectory(src/sgpp/distributedcombigrid/)
add_subdirectory(examples)
add_subdirectory(third_level_manager)
add_subdirectory(tools)

option(DISCOTEC_TEST "Build with Bost tests" ON) #TODO: sinvoller Bescreiben oder in Boost auslagern
if(DISCOTEC_TEST)
    enable_testing()
    add_subdirectory(tests)
endif()









