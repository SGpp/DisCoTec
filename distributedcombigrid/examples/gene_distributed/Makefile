CC=mpicxx.mpich
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall -fpic -ggdb #-DTIMING

SGPP_DIR=/home/sgsscratch/goehrims/combi
GLPK_DIR=/home/sgsscratch/goehrims/combi/glpk

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
INC_GLPK=-I$(GLPK_DIR)/include
LD_GLPK=-L$(GLPK_DIR)/lib
INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 

LDIR=$(LD_SGPP) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_GLPK)

LIBS=-lsgppdistributedcombigrid -lboost_serialization #-lglpk#-DDEBUG_OUTPUT

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: manager lib timing errorCalc

manager: $(OBJ_FILES)
	mpicxx.mpich $(CFLAGS) $(LDIR) $(INC) -o $@ $^ $(LIBS)

obj/%.o: src/%.cpp
	mpicxx.mpich $(CFLAGS) $(LDIR) $(INC) -c -o $@ $< $(LIBS)

lib: $(OBJ_FILES)
	mpicxx.mpich -shared -o lib/libgenetask.so $^

timing: tools/timing.c tools/timing.h
	mpicxx.mpich $(CFLAGS) -c tools/timing.c -o tools/timing.o

errorCalc: tools/errorCalcNew.cpp
	mpicxx.mpich $(CFLAGS) $(LDIR) $(INC) -o errorCalc tools/errorCalcNew.cpp $(LIBS) tools/timing.o

clean:
	rm obj/*.o out/* manager lib/*
