#!/bin/bash

SGPP_DIR=/home/marci/UNI/HIWI/combi/
export LD_LIBRARY_PATH=$SGPP_DIR/lib/sgpp:$LD_LIBRARY_PATH

THIRD_LEVEL_MANAGER=$SGPP_DIR/distributedcombigrid/third_level_manager

#while true; do
#  $THIRD_LEVEL_MANAGER/thirdLevelManager $THIRD_LEVEL_MANAGER/example.ini
#  sleep 3
#done

#gdb -ex run --args $THIRD_LEVEL_MANAGER/thirdLevelManager $THIRD_LEVEL_MANAGER/example.ini
$THIRD_LEVEL_MANAGER/thirdLevelManager $THIRD_LEVEL_MANAGER/example.ini
