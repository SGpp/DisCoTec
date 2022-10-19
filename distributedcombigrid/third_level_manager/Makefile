CC=mpicxx
CFLAGS=-std=c++17 -O3 -DNDEBUG -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall -DUNIFORMDECOMPOSITION -DTIMING

SGPP_SRC=../src/
THIRD_LEVEL_INC=$(SGPP_SRC)/sgpp/distributedcombigrid/third_level/

BOOST_INC=
INC=-I$(SGPP_SRC) $(BOOST_INC)

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp\
                   Params.cpp Params.hpp\
                   System.cpp System.hpp Stats.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.cpp -lboost_program_options
	$(CC) $(CFLAGS) $(INC) -o thirdLevelManager thirdLevelManager.cpp\
			    -DBROKER_ON_SYSTEM=0\
                            Params.cpp System.cpp Stats.cpp\
                            $(THIRD_LEVEL_INC)/NetworkUtils.cpp -lboost_program_options -DDEBUG_ON # > build.log 2>&1

testUFTPsystem: test_uftp.cpp
	$(CC) $(CFLAGS) $(INC) -o testUFTPsystem test_uftp.cpp \
                            -L$(SGPP_SRC)/../../lib/sgpp -lsgppdistributedcombigrid

clean:
	rm -f *.o out/* thirdLevelManager testUFTPsystem
