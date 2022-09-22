#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client

PATHLRZ=/hppfs/work/pn34mi/di39qun2/DisCoTec-third-level/distributedcombigrid/third_level_manager/
PATHHAWK=/lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely/DisCoTec-third-level/distributedcombigrid/third_level_manager/
FILELRZ=$PATHLRZ/dsg_part_0_0
FILEHAWK=$PATHHAWK/dsg_part_1_0
#TODO test
HAWKURL=gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS

# start with one iteration of copying file back and forth
while true
do
if test -f "$FILELRZ" then
    # ./copy_sng_to_hawk.sh
    # try to copy the file
    uftp cp  -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $FILELRZ $HAWKURL:$PATHHAWK
    break
fi
done

while true
do
if test -f "$FILEHAWK" then
    # ./copy_hawk_to_sng.sh
    # try to copy the file
    uftp cp  -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$FILEHAWK $PATHLRZ/
fi
sleep 20s
done
