#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client
#set -x
set -e
PATHLRZ=${PATHLRZ:=//hppfs/scratch/0F/di93yuw/discotec_widely_distributed/numerical_target_64_both/}
PATHHAWK=${PATHHAWK:=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely-distributed-ct/scenarios/numerical_target_64_both/}

HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=${USERHAWK:=ipvpolli}


#TRANSFER_BACKWARD=${FILETOHAWK}
#TOKEN_TRANSFER_BACKWARD=${FILETOHAWK:0:-2}_token.txt
#TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt

TOKEN_STOP=uftp_transfer_stop.txt

step=0


# loop that is infinite until $TOKEN_STOP is created
while true
do
    FILETOHAWK=${PATHLRZ}/dsg_0_step${step}_0
    FILELRZ=${FILELRZ:=${PATHLRZ}/dsg_1_step${step}_0}
    TRANSFER_BACKWARD=${FILETOHAWK}
    TOKEN_TRANSFER_BACKWARD=${FILETOHAWK:0:-2}_token.txt
    TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt

    if test -f $TRANSFER_BACKWARD; then
      echo "Transfer backward"
      touch $TOKEN_TRANSFER_FORWARD
      sleep 10
    fi
    if test -f $TOKEN_TRANSFER_BACKWARD; then
      echo "Transfer forward finished"
      rm $TRANSFER_BACKWARD
      rm $TOKEN_TRANSFER_BACKWARD
      step=$(($step+1))
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break
    fi
    sleep 1
done
