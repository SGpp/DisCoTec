CC=mpicxx
CFLAGS=-std=c++11 -O3 -fopenmp -Wno-deprecated-declarations -Wno-unused-local-typedefs -Wno-deprecated -Wno-uninitialized -Wall

SGPP_SRC=../src/
THIRD_LEVEL_INC=$(SGPP_SRC)/sgpp/distributedcombigrid/third_level/

BOOST_INC=
INC=-I$(SGPP_SRC) $(BOOST_INC)

all: thirdLevelManager

thirdLevelManager: thirdLevelManager.cpp thirdLevelManager.hpp Params.cpp\
                   Params.hpp System.cpp System.hpp Stats.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.hpp\
                   $(THIRD_LEVEL_INC)/NetworkUtils.cpp
	$(CC) $(CFLAGS) $(INC) -o thirdLevelManager thirdLevelManager.cpp\
                            Params.cpp System.cpp Stats.cpp\
                            $(THIRD_LEVEL_INC)/NetworkUtils.cpp # -DDEBUG_ON > build.log 2>&1

clean:
	rm -f *.o out/* thirdLevelManager
