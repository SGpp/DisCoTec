CC=mpic++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall -fpic -ggdb -rdynamic

SGPP_DIR=$(HOME)/combi-gene-faults

USRNAME=ipvalf
LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
LD_BOOST=-L/zhome/academic/HLRS/ipv/$(USRNAME)/hlrs-tools/boost_1_58_0/stage/lib
LD_GLPK=-L/zhome/academic/HLRS/ipv/$(USRNAME)/hlrs-tools/glpk/lib

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/ 
INC_BOOST=-I/zhome/academic/HLRS/ipv/$(USRNAME)/hlrs-tools/boost_1_58_0
INC_GLPK=-I/zhome/academic/HLRS/ipv/$(USRNAME)/hlrs-tools/glpk/include

LDIR=$(LD_SGPP) $(LD_BOOST) $(LD_GLPK)
INC=$(INC_SGPP) $(INC_BOOST) $(INC_GLPK)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -lglpk -lboost_exception -lboost_mpi

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

all: manager lib

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
