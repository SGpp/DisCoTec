#!/bin/bash

LIB_SGPP_DIR=/home/marci/UNI/HIWI/combi/
LIB_SIMPLE_AMQP_DIR=$LIB_SGPP_DIR/distributedcombigrid/SimpleAmqpClient/build
export LD_LIBRARY_PATH=$LIB_SGPP_DIR/lib/sgpp:$LIB_SIMPLE_AMQP_DIR:$LD_LIBRARY_PATH

MANAGER_DIR=$LIB_SGPP_DIR/distributedcombigrid/examples/gene_distributed_third_level/third_level_manager

$MANAGER_DIR/thirdLevelManager $MANAGER_DIR/example.ini
