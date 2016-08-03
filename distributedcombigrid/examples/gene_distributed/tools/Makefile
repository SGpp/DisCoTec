CC=mpic++.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall 
CFLAGS+=-DSIMFTMPI -rdynamic

SGPP_DIR=$(HOME)/path/to/SGpp

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ -I$(SGPP_DIR)/combigrid/src/

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

LIBS=-Wl,-Bdynamic -lsgppdistributedcombigrid -lsgppcombigrid -lboost_serialization 

all: testBoundaryZ

testBoundaryZ: testBoundaryZ.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o testBoundaryZ testBoundaryZ.cpp $(LIBS)


