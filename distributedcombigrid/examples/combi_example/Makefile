CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wall -Wextra -Wno-unused-local-typedefs 

SGPP_DIR=../../..
GLPK_DIR=/home_local/oberstei/Documents/ExaHD/hlrs-tools/glpk

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 
LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
INC_GLPK=-I$(GLPK_DIR)/include
LD_GLPK=-L$(GLPK_DIR)/lib

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)

# LDFLAGS=-Wl,-V

LIBS=-lsgppdistributedcombigrid -lboost_serialization 

all: combi_example 

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o out/* combi_example
