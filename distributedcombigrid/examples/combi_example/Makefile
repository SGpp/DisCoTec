CC=mpic++
CFLAGS=-std=c++11 -g -fopenmp -Wall -Wextra

SGPP_DIR=$(HOME)/path/to/SGpp

LD_SGPP=-L$(SGPP_DIR)/lib/sgpp

INC_SGPP=-I$(SGPP_DIR)/distributedcombigrid/src/

LDIR=$(LD_SGPP)
INC=$(INC_SGPP)

# LDFLAGS=-Wl,-V

LIBS=-lsgppdistributedcombigrid -lboost_serialization 

all: combi_example 

combi_example: combi_example.cpp TaskExample.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o combi_example combi_example.cpp $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o out/* combi_example
