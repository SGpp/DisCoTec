CC=/opt/mpich/bin/mpic++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=../../../
GLPK_DIR=$(SGPP_DIR)/glpk/

LD_SGPP=$(SGPP_DIR)/lib/sgpp
INC_SGPP=$(SGPP_DIR)/distributedcombigrid/src/
LD_GLPK=$(GLPK_DIR)/lib
INC_GLPK=$(GLPK_DIR)/include
LD_BOOST=
INC_BOOST=

LDIR=-L$(LD_SGPP) -L$(LD_GLPK) -L$(LD_BOOST)
INC=-I$(INC_SGPP) -I$(INC_GLPK) -I$(INC_BOOST)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -lglpk

all: combi_example

combi_example: combi_example.cpp TaskConstParaboloid.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LIBS) -DUNIFORMDECOMPOSITION -DTIME

clean:
	rm -f *.o out/* combi_example
