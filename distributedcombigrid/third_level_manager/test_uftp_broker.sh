#!/bin/bash

FILELRZ=/hppfs/work/pn34mi/di39qun2/DisCoTec-third-level/distributedcombigrid/third_level_manager/dsg_part_0_0
FILEHAWK=/lustre/hpe/ws10/ws10.3/ws/TODO/dsg_part_1_0
#TODO test
HAWKURL=gridftp-fr1.hww.hlrs.de:9000

# start with one iteration of copying file back and forth
while true
do
if test -f "$FILELRZ" then
    # ./copy_sng_to_hawk.sh
    # try to copy the file
    uftp cp $FILELRZ $FILEHAWK
fi
done

while true
if test -f "$FILEHAWK" then
    # ./copy_hawk_to_sng.sh
    # try to copy the file
    uftp cp $HAWKURL/$FILEHAWK $FILELRZ
fi
sleep 20s
done

