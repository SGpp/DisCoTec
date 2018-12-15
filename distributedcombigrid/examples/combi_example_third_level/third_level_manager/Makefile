CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_DIR=$(SGPP)
GLPK_DIR=$(GLPK)

RABBIT_MQ_CLIENT_LIB=../SimpleAmqpClient/build/
RABBIT_MQ_CLIENT_INC=../SimpleAmqpClient/src/

LDIR=-L$(RABBIT_MQ_CLIENT_LIB)
INC=-I$(RABBIT_MQ_CLIENT_INC)

LIBS=-lSimpleAmqpClient

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o thirdLevelManager thirdLevelManager.cpp $(LIBS)

clean:
	rm -f *.o out/* thirdLevelManager
