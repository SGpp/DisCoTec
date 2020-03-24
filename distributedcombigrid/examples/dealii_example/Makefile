CC=mpicxx
CFLAGS=-std=c++17 -march=native -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall -DUNIFORMDECOMPOSITION

SGPP_DIR=/home/sgsscratch/goehrims/combi
GLPK_DIR=/home/sgsscratch/goehrims/combi/glpk
HYPER_DIR=/home/sgsscratch/goehrims/
DEBUG_DIR=/home/sgsscratch/goehrims/DEBUG

DEAL_II_DIR=/home/sgsscratch/goehrims/matrixfree-install
DEAL_II_INC_DIR=/home/sgsscratch/goehrims/matrixfree


LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
INC_GLPK=-I$(GLPK_DIR)/include
LD_GLPK=-L$(GLPK_DIR)/lib
INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 

LD_DEALII=-L$(DEAL_II_DIR)/lib
INC_DEALII=-I$(DEAL_II_DIR)/include 
INC_BOOST=-I$(DEAL_II_INC_DIR)/bundled/boost-1.70.0/include
INC_P4EST= -I$(DEBUG_DIR)/include

INC_HYPERDEAL=-I$(HYPER_DIR)

LDIR=$(LD_SGPP) $(LD_GLPK) $(LD_DEALII)
INC=$(INC_SGPP) $(INC_GLPK) $(INC_DEALII) $(INC_BOOST) $(INC_HYPERDEAL) $(INC_P4EST)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -ldeal_II

all: dealii_example 

dealii_example: dealii_example.cpp TaskExample.hpp DataConverter.cpp parameter_handler.cc TaskEnsemble.hpp
	mpicxx $(CFLAGS) $(LDIR) $(INC) -o dealii_example dealii_example.cpp parameter_handler.cc $(LIBS)

clean:
	rm -f *.o out/* dealii_example
