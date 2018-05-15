CC=mpic++.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(HOME)/path/to/SGpp

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

LIBS=-lsgppdistributedcombigrid -lboost_serialization 

all: combi_third_level_example

combi_example: combi_third_level_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_third_level_example combi_third_level_example.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_third_level_example
