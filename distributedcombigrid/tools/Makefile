CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall 
CFLAGS+=-rdynamic
#you may need to add some more compile flags

SGPP_DIR=$(realpath $(shell pwd)/../../)

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp -L$(SGPP_DIR)/glpk/lib -lglpk

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ -I../src

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

LIBS=-Wl,-Bdynamic -lsgppdistributedcombigrid -lboost_serialization

errorCalc: errorCalc.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o errorCalc errorCalc.cpp $(LIBS)