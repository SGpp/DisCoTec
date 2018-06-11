CC=mpic++
CFLAGS=-std=c++14 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(HOME)/Thesis/combi

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -lglpk

all: combi_example

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_example
