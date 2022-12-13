cmake_minimum_required(VERSION 3.18)

project("DisCoTec"
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




option(DISCOTEC_GENE "Build for GENE (as task library)" OFF)
option(DISCOTEC_TIMING "Build with high-resolution timing" ON) # Done
option(DISCOTEC_OPENMP "Parallelize with OpenMP" OFF)

option(DISCOTEC_HDF5 "Interpolation output with HDF5" OFF) # ON


option(DISCOTEC_ENABLEFT "DisCoTec with algorithm-based fault tolerance" OFF)
if(DISCOTEC_ENABLEFT)
    # TODO: package handeling
    add_compile_definitions(ENABLEFT)
endif()

option(DISCOTEC_DEBUG "Build with assertions and debug symbols" OFF) #TODO: Move to build type
if(DISCOTEC_DEBUG)
    add_compile_definitions(DEBUG_OUTPUT) #TODO: Add only to specific targets if they are sorted
endif()



#USE_NON_BLOCKING_COLLECTIVES #wird nirgends in SCONS definiert also ka wann definieren
#ALT_LEVEL_VECTOR #wird nirgends in SCONS definiert also ka wann definieren
#HAVE_HIGHFIVE #wird nirgends in SCONS definiert also ka wann definieren
#USE_VTK
#SIMFTMPI
#USENONBLOCKINGMPICOLLECTIVE
#OMITREADYSIGNAL
#UNIFORMDECOMPOSITION
#ISGENE
