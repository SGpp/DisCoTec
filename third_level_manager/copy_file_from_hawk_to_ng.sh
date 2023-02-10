#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client

PATHLRZ=//hppfs/scratch/0F/di93yuw/discotec_widely_distributed/numerical_target_64_both/
PATHHAWK=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely-distributed-ct/scenarios/numerical_target_64_both/
FILEHAWK=${PATHHAWK}/dsg_0_step*_0
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=ipvpolli

TOKEN_TRANSFER_BACKWARD=${FILEHAWK:0:-2}_token.txt
TOKEN_STOP=uftp_transfer_stop.txt

echo "$FILEHAWK"
echo "$TOKEN_TRANSFER_BACKWARD"

step=0

# loop that is infinite until $TOKEN_STOP is created
while true
do
    FILEHAWK=${PATHHAWK}/dsg_0_step${step}_0
    TOKEN_TRANSFER_BACKWARD=${FILEHAWK:0:-2}_token.txt
    # test if the file is there by trying to copy it (maybe there is a better way?)
    # the part after `>` is to ignore the potential (error) output
    uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$TOKEN_TRANSFER_BACKWARD /dev/null > /dev/null 2>&1
    cp_status=$?
    if (exit $cp_status); then
        # ./copy_hawk_to_sng.sh
        # try to copy the file
        FILEHAWK_INSTANCE=$(echo $TOKEN_TRANSFER_BACKWARD)
        FILEHAWK_INSTANCE=${FILEHAWK_INSTANCE:0:-10}_0
        echo copy $FILEHAWK_INSTANCE from hawk to NG
        starttime=`date +%s`
        uftp cp -t 20 -n 4 -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$FILEHAWK_INSTANCE $PATHLRZ
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$TOKEN_TRANSFER_BACKWARD $PATHLRZ
        uftp rm --quiet -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK "$HAWKURL:$TOKEN_TRANSFER_BACKWARD"
        endtime=`date +%s`
        echo copied "$FILEHAWK_INSTANCE" in `expr $endtime - $starttime` seconds.
        step=$(($step+1))
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break 
    fi
    sleep 1
done
