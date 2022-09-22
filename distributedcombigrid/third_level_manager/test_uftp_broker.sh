#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client

PATHLRZ=//hppfs/scratch/02/di39qun2/third_level_manager/
PATHHAWK=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely/DisCoTec-third-level/distributedcombigrid/third_level_manager/
FILELRZ=dsg_part_0_*
FILEHAWK=dsg_part_1_*
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=ipvpolli

TOKEN_TRANSFER_FORWARD=uftp_transfer_0.txt
TOKEN_TRANSFER_BACKWARD=uftp_transfer_1.txt
TOKEN_STOP=uftp_transfer_stop.txt

echo "$FILELRZ" "$FILEHAWK"

# loop that is infinite until $TOKEN_STOP is created
while true
do
    # one iteration of copying file back and forth
    while true
    do
    if test -f "$PATHLRZ/$TOKEN_TRANSFER_FORWARD"; then
        # ./copy_sng_to_hawk.sh
        # try to copy the file
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $PATHLRZ/$FILELRZ $HAWKURL:$PATHHAWK
        rm "$PATHLRZ/$TOKEN_TRANSFER_FORWARD"
        echo copied "$FILELRZ"
        break
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break 2
    fi
    done
    
    while true
    do
    # test if the file is there by trying to copy it (maybe there is a better way?)
    # the part after `>` is to ignore the potential (error) output
    uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK "$HAWKURL:$PATHHAWK/$TOKEN_TRANSFER_BACKWARD" /dev/null > /dev/null 2>&1
    cp_status=$?
    if (exit $cp_status); then
        # ./copy_hawk_to_sng.sh
        # try to copy the file
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$PATHHAWK/$FILEHAWK $PATHLRZ
        uftp rm --quiet -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK "$HAWKURL:$PATHHAWK/$TOKEN_TRANSFER_BACKWARD"
        echo copied "$FILEHAWK"
        break
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break 2
    fi
    #sleep 20s
    done
done
