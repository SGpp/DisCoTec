CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=/home/marci/UNI/HIWI/combi/

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp
LD_MPI=-L/opt/mpich/lib
LD_SIMPLE_AMQP=-L$(SGPP_DIR)/distributedcombigrid/SimpleAmqpClient/build

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/
INC_MPI=-I/opt/mpich/include
INC_SIMPLE_AMQP=-I$(SGPP_DIR)/distributedcombigrid/SimpleAmqpClient/src/

LDIR=$(LD_SGPP) $(LD_MPI) $(LD_SIMPLE_AMQP)
INC=$(INC_SGPP) $(INC_MPI) $(INC_SIMPLE_AMQP)

LIBS=-lsgppdistributedcombigrid -lboost_serialization -lmpi -lSimpleAmqpClient

all: combi_example

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LIBS)

clean:
	rm -f *.o out/* combi_example
