#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)
# example call: PATHLRZ=/hppfs/work/pn36xu/di39qun2/scenarios/weak_widely/weak_asymmetric_2  USERHAWK=ipvsavcr PATHHAWK=/lustre/hpe/ws10/ws10.3/ws/ipvsavcr-discotec/scenarios/weak_widely/weak_asymmetric_2 /hppfs/scratch/0F/di93yuw/DisCoTec/third_level_manager/copy_file_from_hawk_to_ng.sh

module load uftp-client


# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# Paths and filenames with default values
# (can be overwritten by environment variables or command line arguments)
PATHLRZ=${PATHLRZ:=//hppfs/scratch/0F/di93yuw/discotec_widely_distributed/numerical_target_64_both/}
PATHHAWK=${PATHHAWK:=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely-distributed-ct/scenarios/numerical_target_64_both/}
FILEHAWK=${FILEHAWK:=${PATHHAWK}/dsg_0_step*_0}
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=${USERHAWK:=ipvpolli}

# Tokens for the file transfer. If the file is not there, the script will wait until it is there.
TOKEN_TRANSFER_BACKWARD=${FILEHAWK:0:-2}_token.txt
# Or stops if this file is created.
TOKEN_STOP=uftp_transfer_stop.txt

# Default values for the parallel file transfer
# (can be overwritten by environment variables or command line arguments)
PROCS=${PROCS:=4} # number of independent parallel transfers (each transfer uses $THREADS_PER_PROC threads)
THREADS_PER_PROC=${THREADS_PER_PROC:=8} # number of uftp threads per proc
STREAMS=${STREAMS:=3} # number of parallel uftp streams per thread

echo "$FILEHAWK"
echo "$TOKEN_TRANSFER_BACKWARD"

step=0

# loop that is infinite until $TOKEN_STOP is found
while true
do
    # check for the current ready token
    FILEHAWK=${PATHHAWK}/dsg_0_step${step}_0
    TOKEN_TRANSFER_BACKWARD=${FILEHAWK:0:-2}_token.txt
    # test if the file is there by trying to copy it (maybe there is a better way?)
    # the part after `>` is to ignore the potential (error) output
    uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$TOKEN_TRANSFER_BACKWARD /dev/null > /dev/null 2>&1
    cp_status=$?
    pids=( )
    if (exit $cp_status); then
        # ./copy_hawk_to_sng.sh
        # try to copy the file
        FILEHAWK_INSTANCE=$(echo $TOKEN_TRANSFER_BACKWARD)
        FILEHAWK_INSTANCE=${FILEHAWK_INSTANCE:0:-10}_0
        echo copy $FILEHAWK_INSTANCE from hawk to NG
        size=$(uftp ls -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK  $HAWKURL:$FILEHAWK_INSTANCE | awk '{print $2 }')
        startblock=0
        endblock=$((size/PROCS))
        starttime=`date +%s`
        for ((i=1; i<=PROCS; i++)); do
          echo "Block: $startblock-$endblock of $size"
          uftp cp -n $STREAMS -B "${startblock}-${endblock}-p" -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$FILEHAWK_INSTANCE $PATHLRZ/ &
          pids+=($!) # store background pids
          startblock=$((endblock))
          if [ $i -eq $((PROCS)) ]; then
              endblock=$((size))
          else
              endblock=$((endblock+size/PROCS))
          fi
        done
        # wait for all processes to finish
        for pid in "${pids[@]}"; do
            wait $pid || exit
        done
        # copy the finish token file witch signals that the copie is finished
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$TOKEN_TRANSFER_BACKWARD $PATHLRZ/
        # remove the start token file on hawk
        uftp rm --quiet -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK "$HAWKURL:$TOKEN_TRANSFER_BACKWARD"

        # print the throughput
        endtime=`date +%s`
        echo copied "$FILEHAWK_INSTANCE" in `expr $endtime - $starttime` seconds.
        throughput=$( echo "scale=4;($size/1024/1024)/(($endtime-$starttime))" | bc )
        throughput_bits=$( echo "scale=4;($throughput*8)" | bc )
        echo "Average throughput: $throughput MB/s; $throughput_bits Mbit/s"
        step=$(($step+1))

        # remove the copied file on hawk
        uftp rm -q -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $HAWKURL:$FILEHAWK_INSTANCE
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break
    fi
    # wait for 1 second to reduce the load on the system
    sleep 1
done
