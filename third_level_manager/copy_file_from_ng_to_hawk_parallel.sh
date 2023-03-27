#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)

module load uftp-client

set -e

# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

ID=${ID:=0}
PATHLRZ=${PATHLRZ:=//hppfs/work/pn36xu/di93yuw/uftp}
PATHHAWK=${PATHHAWK:=/lustre/hpe/ws10/ws10.3/ws/ipvsavcr-discotec/wd_alex/scaled_down_target_to_try_64/}
FILELRZ=${FILELRZ:=${PATHLRZ}/dsg_${ID}_step*_0}
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=${USERHAWK:=ipvsavcr}
num_hosts=3


TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt
TOKEN_STOP=uftp_transfer_stop.txt

PROCS=${PROCS:=16}
THREADS_PER_PROC=${THREADS_PER_PROC:=8}
STREAMS=${STREAMS:=3}

echo "$FILELRZ"
echo "$TOKEN_TRANSFER_FORWARD"

step=0

# loop that is infinite until $TOKEN_STOP is created
while true
do
    FILELRZ=${PATHLRZ}/dsg_${ID}_step${step}_0
    TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt

    if test -f $TOKEN_TRANSFER_FORWARD; then
        # ./copy_sng_to_hawk.sh
        # try to copy the file
        FILELRZ_INSTANCE=$(echo $TOKEN_TRANSFER_FORWARD)
        FILELRZ_INSTANCE=${FILELRZ_INSTANCE:0:-10}_0
        echo instance "$FILELRZ_INSTANCE"
        size=$(stat --printf="%s" $FILELRZ_INSTANCE)
        local_size=$((size/num_hosts))
        host=$(hostname)
        if [ $host == "login05" ]; then
            local_start=0
        elif [ $host == "login06" ]; then
            local_start=$((local_size))
        elif [ $host == "login07" ]; then
            local_start=$((local_size*2))
        fi
        startblock=$local_start
        endblock=$((local_start+(local_size/PROCS)))
        starttime=`date +%s`
        pids=( )
        for ((i=1; i<=PROCS; i++)); do
          echo "Block: $startblock-$endblock of $size"
          uftp cp -B "${startblock}-${endblock}-p" -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $FILELRZ_INSTANCE $HAWKURL:$PATHHAWK/ &
          pids+=($!) # store background pids
          startblock=$((endblock))
          if [ $i -eq $((PROCS-1)) ]; then
              if [ $host ==  "login07" ]; then
                        endblock=$((size))
              else
                endblock=$((local_start+local_size))
              fi
          else
              endblock=$((endblock+local_size/PROCS))
          fi
        done
        # wait for all pids
        for pid in "${pids[@]}"; do
          wait "$pid" || exit
        done
        touch "${TOKEN_TRANSFER_FORWARD}_${host}"
        if [ $host == "login05" ]; then
            while [ ! -f ${TOKEN_TRANSFER_FORWARD}_login05 ] || [ ! -f ${TOKEN_TRANSFER_FORWARD}_login06 ] || [ ! -f ${TOKEN_TRANSFER_FORWARD}_login07 ]; do
                sleep 1
            done
            rm ${TOKEN_TRANSFER_FORWARD}_login05
            rm ${TOKEN_TRANSFER_FORWARD}_login06
            rm ${TOKEN_TRANSFER_FORWARD}_login07
            uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $TOKEN_TRANSFER_FORWARD $HAWKURL:$PATHHAWK/
            rm -f $TOKEN_TRANSFER_FORWARD
        fi

        endtime=`date +%s`
        echo "copied $local_size bites from file  "$FILELRZ_INSTANCE" in `expr $endtime - $starttime` seconds."

        throughput=$( echo "scale=4;($local_size/1024/1024)/(($endtime-$starttime))" | bc )
        throughput_bits=$( echo "scale=4;($throughput*8)" | bc )
        echo "Local average  throughput: $throughput MB/s; $throughput_bits Mbit/s"
        if [ $host == "login05" ]; then
            throughput=$( echo "scale=4;($size/1024/1024)/(($endtime-$starttime))" | bc )
            throughput_bits=$( echo "scale=4;($throughput*8)" | bc )
            echo "Approx total average throughput: $throughput MB/s; $throughput_bits Mbit/s"
            rm -f $FILELRZ_INSTANCE
        fi
        date
        step=$((step+1))
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break
    fi
    sleep 1
done
