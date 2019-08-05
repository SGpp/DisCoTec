CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

RABBIT_MQ_CLIENT_LIB=../SimpleAmqpClient/build/
RABBIT_MQ_CLIENT_INC=../SimpleAmqpClient/src/

SGPP_SRC=../src/
THIRD_LEVEL_INC=$(SGPP_SRC)/sgpp/distributedcombigrid/third_level/

LDIR=-L$(RABBIT_MQ_CLIENT_LIB)
INC=-I$(RABBIT_MQ_CLIENT_INC) -I$(SGPP_SRC)

LIBS=-lSimpleAmqpClient

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp Params.cpp\
                   Params.hpp System.cpp System.hpp\
                   $(THIRD_LEVEL_INC)/MessageUtils.cpp\
                   $(THIRD_LEVEL_INC)/MessageUtils.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.cpp
	$(CC) $(CFLAGS) $(LDIR) $(INC) -o thirdLevelManager thirdLevelManager.cpp\
                                    Params.cpp System.cpp\
                                    $(THIRD_LEVEL_INC)/MessageUtils.cpp\
                                    $(THIRD_LEVEL_INC)/NetworkUtils.cpp $(LIBS) #> build.log 2>&1

clean:
	rm -f *.o out/* thirdLevelManager
