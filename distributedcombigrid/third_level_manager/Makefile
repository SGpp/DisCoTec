CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(SGPP)

RABBIT_MQ_CLIENT_LIB=../../../SimpleAmqpClient/build/
RABBIT_MQ_CLIENT_INC=../../../SimpleAmqpClient/src/

SGPP_SRC=../../../src/
THIRD_LEVEL=$(SGPP_SRC)/sgpp/distributedcombigrid/third_level/

MPI_LIB=/opt/mpich/lib/
MPI_INC=/opt/mpich/include/

LDIR=-L$(RABBIT_MQ_CLIENT_LIB) -L$(MPI_LIB)
INC=-I$(RABBIT_MQ_CLIENT_INC) -I$(SGPP_SRC) -I$(MPI_INC)

LIBS=-lSimpleAmqpClient

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp Params.cpp Params.hpp System.cpp System.hpp $(THIRD_LEVEL)/MessageUtils.cpp $(THIRD_LEVEL)/MessageUtils.hpp $(THIRD_LEVEL)/NetworkUtils.hpp $(THIRD_LEVEL)/NetworkUtils.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o thirdLevelManager thirdLevelManager.cpp Params.cpp System.cpp $(THIRD_LEVEL)/MessageUtils.cpp $(THIRD_LEVEL)/NetworkUtils.cpp $(LIBS) #> build.log 2>&1

clean:
	rm -f *.o out/* thirdLevelManager
