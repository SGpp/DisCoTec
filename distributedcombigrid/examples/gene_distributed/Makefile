CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall -fpic -ggdb

SGPP_DIR=../../..

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
LD_GLPK=-L$(SGPP_DIR)/../hlrs-tools/glpk/lib
INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 
INC_GLPK=-I$(SGPP_DIR)/../hlrs-tools/glpk/include

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -lglpk #-DDEBUG_OUTPUT

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: manager lib errorCalc

manager: $(OBJ_FILES)
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o $@ $^ $(LIBS)

obj/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -c -o $@ $< $(LIBS)

lib: $(OBJ_FILES)
	$(CC) -shared -o lib/libgenetask.so $^

errorCalc: tools/errorCalcNew.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o errorCalc tools/errorCalcNew.cpp $(LIBS) tools/timing.o

clean:
	rm obj/*.o out/* manager lib/*
