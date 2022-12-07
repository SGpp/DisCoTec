cmake_minimum_required(VERSION 3.18)

project("DisCoTec"
        VERSION 0.1.0
        LANGUAGES CXX
        DESCRIPTION "TODO")

add_subdirectory(examples)
add_subdirectory(src/sgpp/distributedcombigrid/)
enable_testing()
add_subdirectory(tests)
add_subdirectory(third_level_manager)
add_subdirectory(tools)