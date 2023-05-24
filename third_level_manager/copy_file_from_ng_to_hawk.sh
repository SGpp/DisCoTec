#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de
# as described on https://doku.lrz.de/display/PUBLIC/Data+Transfer+Options+on+SuperMUC-NG#DataTransferOptionsonSuperMUCNG-UNICOREFileTransfer(UFTP)
#example call: PATHLRZ=/hppfs/work/pn36xu/di39qun2/scenarios/weak_widely/weak_asymmetric_2  USERHAWK=ipvsavcr PATHHAWK=/lustre/hpe/ws10/ws10.3/ws/ipvsavcr-discotec/scenarios/weak_widely/weak_asymmetric_2 /hppfs/scratch/0F/di93yuw/DisCoTec/third_level_manager/copy_file_from_ng_to_hawk.sh

module load uftp-client

# exit on error
set -e

# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# Paths and filenames with default values
# (can be overwritten by environment variables or command line arguments)
PATHLRZ=${PATHLRZ:=//hppfs/scratch/0F/di93yuw/discotec_widely_distributed/numerical_target_64_both/}
PATHHAWK=${PATHHAWK:=//lustre/hpe/ws10/ws10.1/ws/ipvpolli-widely-distributed-ct/scenarios/numerical_target_64_both/}
FILELRZ=${FILELRZ:=${PATHLRZ}/dsg_1_step*_0}
HAWKURL=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
USERHAWK=${USERHAWK:=ipvpolli}

# Tokens for the file transfer. If the file is not there, the script will wait until it is there.
TOKEN_TRANSFER_FORWARD=${FILELRZ:0:-2}_token.txt
# Or stops if this file is created.
TOKEN_STOP=uftp_transfer_stop.txt

# Default values for the parallel file transfer
# (can be overwritten by environment variables or command line arguments)
PROCS=${PROCS:=4} # number of independent parallel transfers (each transfer uses $THREADS_PER_PROC threads)
THREADS_PER_PROC=${THREADS_PER_PROC:=8} # number of uftp threads per proc
STREAMS=${STREAMS:=3} # number of parallel uftp streams per thread




echo "$FILELRZ"

# loop that is infinite until $TOKEN_STOP is found
while true
do
    # check for the current ready token
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
        pids=( )
        for ((i=1; i<=PROCS; i++)); do
          echo "Block: $startblock-$endblock of $size"
          uftp cp -n $STREAMS -B "${startblock}-${endblock}-p" -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $FILELRZ_INSTANCE $HAWKURL:$PATHHAWK/ &
          pids+=($!) # store background pids
          startblock=$((endblock))
          if [ $i -eq $((PROCS)) ]; then
              endblock=$((size))
          else
              endblock=$((endblock+size/PROCS))
          fi
        done
        # wait for all pids
        for pid in "${pids[@]}"; do
          wait "$pid" || exit
        done

        # copy the finish token file witch signals that the copie is finished
        uftp cp -i ~/.uftp/id_uftp_to_hlrs -u $USERHAWK $TOKEN_TRANSFER_FORWARD $HAWKURL:$PATHHAWK/
        # remove the start token
        rm -f $TOKEN_TRANSFER_FORWARD

        # print the throughput
        endtime=`date +%s`
        echo copied "$FILELRZ_INSTANCE" in `expr $endtime - $starttime` seconds.
        throughput=$( echo "scale=4;($size/1024/1024)/(($endtime-$starttime))" | bc )
        throughput_bits=$( echo "scale=4;($throughput*8)" | bc )
        echo "Average throughput: $throughput MB/s; $throughput_bits Mbit/s"

        # remove the copied file
        rm $FILELRZ_INSTANCE
    fi
    if test -f "$PATHLRZ/$TOKEN_STOP"; then
        break
    fi
    # wait for 1 second to reduce the load on the system
    sleep 1
done
