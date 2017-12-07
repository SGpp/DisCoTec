#!/bin/bash
. module.sh
USRNAME="ipvober"
scons -j 16 OPT=1 VERBOSE=1 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 CXX=CC COMPILE_BOOST_TESTS=0 RUN_BOOST_TESTS=0 BOOST_LIBRARY_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/boost_1_58_0/stage/lib BOOST_INCLUDE_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/boost_1_58_0 GLPK_LIBRARY_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/glpk/lib GLPK_INCLUDE_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/glpk/include BUILD_STATICLIB=1

scons -j 16 OPT=1 VERBOSE=1 SG_ALL=0 SG_DISTRIBUTEDCOMBIGRID=1 CXX=CC COMPILE_BOOST_TESTS=0 RUN_BOOST_TESTS=0 BOOST_LIBRARY_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/boost_1_58_0/stage/lib BOOST_INCLUDE_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/boost_1_58_0 GLPK_LIBRARY_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/glpk/lib GLPK_INCLUDE_PATH=/zhome/academic/HLRS/ipv/$USRNAME/hlrs-tools/glpk/include BUILD_STATICLIB=0 
