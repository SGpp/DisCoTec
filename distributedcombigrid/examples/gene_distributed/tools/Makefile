CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall 
CFLAGS+=-rdynamic
#you may need to add some more compile flags

SGPP_DIR=$(realpath $(shell pwd)/../../../../)

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp -L../lib/ -L$(SGPP_DIR)/glpk/lib -lglpk

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ -I../src/ 

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

LIBS=-Wl,-Bdynamic -lsgppdistributedcombigrid -lboost_serialization timing.o

all: testDistributedHierarchization

testBoundaryZ: testBoundaryZ.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o testBoundaryZ testBoundaryZ.cpp -lgenetask $(LIBS)

errorCalc: errorCalc.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o errorCalc errorCalc.cpp $(LIBS)

errorCalcNew: errorCalcNew.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o errorCalcNew errorCalcNew.cpp $(LIBS)

errorCalcGene: errorCalcGene.cpp timing.o
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o errorCalcGene errorCalcGene.cpp -lgenetask $(LIBS)

convertFGtoCP: convertFGtoCP.cpp timing.o
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o convertFGtoCP convertFGtoCP.cpp $(LIBS)

testHierarchization: testHierarchization.cpp timing.o
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o testHierarchization testHierarchization.cpp $(LIBS)

testDistributedHierarchization: testDistributedHierarchization.cpp timing.o
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o testDistributedHierarchization testDistributedHierarchization.cpp $(LIBS)

timing.o: timing.c timing.h
	$(CC) $(CFLAGS) -c timing.c


