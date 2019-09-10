CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall 
#CFLAGS+=-DSIMFTMPI

SGPP_DIR=/home/sgsscratch/goehrims/combi
GLPK_DIR=/home/sgsscratch/goehrims/combi/glpk

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
INC_GLPK=-I$(GLPK_DIR)/include
LD_GLPK=-L$(GLPK_DIR)/lib
INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)



LIBS=-lsgppdistributedcombigrid -lboost_serialization

all: combi_example_faults

combi_example_faults: combi_example_faults.cpp TaskExample.hpp 
	mpicxx.mpich $(CFLAGS) $(LDIR) $(INC) -o combi_example_faults combi_example_faults.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_example_faults
