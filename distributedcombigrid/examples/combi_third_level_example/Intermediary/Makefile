CC=g++
CFLAGS=-std=c++11 -g -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

INC = -I.

all: intermediary

intermediary: intermediary.cpp NetworkUtils.cpp
	$(CC) $(INC) $(CFLAGS) -o intermediary intermediary.cpp NetworkUtils.cpp

clean:
	rm -f *.o out/* intermediary
