cmake_minimum_required(VERSION 3.24.2)

project("DisCoTec Examples"
        LANGUAGES CXX
        DESCRIPTION "Use this to build all DisCoTec examples")


add_subdirectory(combi_example)
add_subdirectory(combi_example_faults)
add_subdirectory(combi_workers_only)
add_subdirectory(distributed_advection)
add_subdirectory(distributed_third_level)
#add_subdirectory(gene_distributed) # TODO: add if access to gene
#add_subdirectory(gene_distributed_linear) # TODO: add if access to gene
#add_subdirectory(gene_distributed_third_level) # TODO: add if access to gene
option(DISCOTEC_WITH_SELALIB "Enable SeLaLib example" OFF)
if(DISCOTEC_WITH_SELALIB)
    add_subdirectory(selalib_distributed)
endif()