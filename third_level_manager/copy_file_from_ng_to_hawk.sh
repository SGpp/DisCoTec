#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client

PATHLRZ=${PATHLRZ:=//hppfs/scratch/0F/di93yuw/discotec_widely_distributed/numerical_target_64_both/}
PATHHAWK=${PATHHAWK:=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely-distributed-ct/scenarios/numerical_target_64_both/}
FILELRZ=${FILELRZ:=${PATHLRZ}/dsg_1_step*_0}
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=${USERHAWK:=ipvpolli}

PROCS=${PROCS:=8}
THREADS_PER_PROC=${THREADS_PER_PROC:=8}
STREAMS=${STREAMS:=2}


TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt
TOKEN_STOP=uftp_transfer_stop.txt

echo "$FILELRZ"

# loop that is infinite until $TOKEN_STOP is created
while true
do
    #echo "$TOKEN_TRANSFER_FORWARD"
    if test -f $TOKEN_TRANSFER_FORWARD; then
        # ./copy_sng_to_hawk.sh
        # try to copy the file
        FILELRZ_INSTANCE=$(echo $TOKEN_TRANSFER_FORWARD)
        FILELRZ_INSTANCE=${FILELRZ_INSTANCE:0:-10}_0
        echo instance "$FILELRZ_INSTANCE"
        size=$(stat --printf="%s" $FILELRZ_INSTANCE)
        startblock=0
        endblock=$((size/PROCS))
        starttime=`date +%s`
        for ((i=1; i<=PROCS; i++)); do
          echo "Block: $startblock-$endblock of $size"
          uftp cp -t $THREADS_PER_PROC -n $STREAMS -C -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $FILELRZ_INSTANCE $HAWKURL:$PATHHAWK &
          startblock=$((endblock+1))
          if [ $i -eq $((PROCS)) ]; then
              endblock=$((size-1))
          else
              endblock=$((endblock+size/PROCS))
          fi
        done
        wait
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $TOKEN_TRANSFER_FORWARD $HAWKURL:$PATHHAWK
        rm $TOKEN_TRANSFER_FORWARD
        endtime=`date +%s`
        echo copied "$FILELRZ_INSTANCE" in `expr $endtime - $starttime` seconds.
        throughput=$( echo "scale=4;($size/1024/1024)/(($endtime-$starttime))" | bc )
        throughput_bits=$( echo "scale=4;($throughput*8)" | bc )
        echo "Average throughput: $throughput MB/s; $throughput_bits Mbit/s"
        rm $FILELRZ_INSTANCE
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break
    fi
    sleep 1
done
