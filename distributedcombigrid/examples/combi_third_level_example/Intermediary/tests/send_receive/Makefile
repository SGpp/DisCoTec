CC=g++
THIRD_LEVEL=../../../../../src/sgpp/distributedcombigrid/third_level
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall
INC = -I. -I$(THIRD_LEVEL)

all: srTest

srTest: srTest.cpp $(THIRD_LEVEL)/NetworkUtils.cpp
	$(CC) $(INC) $(CFLAGS) -o srTest srTest.cpp $(THIRD_LEVEL)/NetworkUtils.cpp -pthread

clean:
	rm -f *.o srTest
