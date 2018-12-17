CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(SGPP)
GLPK_DIR=$(GLPK)

RABBIT_MQ_CLIENT_LIB=SimpleAmqpClient/build/
RABBIT_MQ_CLIENT_INC=SimpleAmqpClient/src/

SGPP_SRC=../../../src/
NETWORK_UTILS=$(SGPP_SRC)/sgpp/distributedcombigrid/third_level/

MPI_LIB=/opt/mpich/lib/
MPI_INC=/opt/mpich/include/

LDIR=-L$(RABBIT_MQ_CLIENT_LIB) -L$(MPI_LIB)
INC=-I$(RABBIT_MQ_CLIENT_INC) -I$(SGPP_SRC) -I$(MPI_INC)

LIBS=-lSimpleAmqpClient #-lboost

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp Params.cpp Params.hpp System.cpp System.hpp $(NETWORK_UTILS)/NetworkUtils.hpp $(NETWORK_UTILS)/NetworkUtils.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o thirdLevelManager thirdLevelManager.cpp Params.cpp System.cpp $(NETWORK_UTILS)/NetworkUtils.cpp $(LIBS) #> build.log 2>&1

clean:
	rm -f *.o out/* thirdLevelManager
