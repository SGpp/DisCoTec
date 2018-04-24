CC=mpic++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=/home_local/oberstei/Documents/ExaHD/combi-ft-main/combi-ft/

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
INC_GLPK=-I/home_local/oberstei/Documents/ExaHD/hlrs-tools/glpk/include
LD_GLPK=-L/home_local/oberstei/Documents/ExaHD/hlrs-tools/glpk/lib
INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)

LIBS=-lsgppdistributedcombigrid -lboost_serialization 

all: combi_example 

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_example
