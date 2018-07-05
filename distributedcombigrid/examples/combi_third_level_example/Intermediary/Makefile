CC=g++
THIRD_LEVEL=../../../src/sgpp/distributedcombigrid/third_level
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

INC = -I. -I$(THIRD_LEVEL)

all: intermediary

intermediary: main.cpp Intermediary.cpp Communicator.cpp Participant.cpp $(THIRD_LEVEL)/NetworkUtils.cpp
	$(CC) $(INC) $(CFLAGS) -o Intermediary main.cpp Intermediary.cpp Communicator.cpp Participant.cpp $(THIRD_LEVEL)/NetworkUtils.cpp -pthread

clean:
	rm -f *.o Intermediary
