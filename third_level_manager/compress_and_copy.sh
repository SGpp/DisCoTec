#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de

module load uftp-client

# exit on error
set -e -x

# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# Paths and filenames with default values
# (can be overwritten by environment variables or command line arguments)
SOURCE_PATH=${SOURCE_PATH:=//p/scratch/widediscotecjsc/widely/widely_measurements/widely_14x512/}
TARGET_PATH=${TARGET_PATH:=//hppfs/work/pn36xu/di39qun2/widely/widely_measurements/widely_150x512/}
SOURCE_FILE=${SOURCE_FILE:=${SOURCE_PATH}/dsg_1_step*_0} #TODO: modify the number accordingly
TARGET_URL=https://datagw03.supermuc.lrz.de:9000/rest/auth/DATAGW
TARGET_USER=${TARGET_USER:=di39qun2}
IDENTITYFILE=${IDENTITYFILE:=~/.uftp/id_uftp_to_lrz}
NUM_FILES=${NUM_FILES:=8}

LZ4_PATH=${LZ4_PATH:=/p/scratch/widediscotecjsc/lz4yume/lz4}

# host_nr=$(( $(hostname | tr -dc '0-9') - 5 ))

# Tokens for the file transfer. If the file is not there, the script will wait until it is there.
TOKEN_TRANSFER_FORWARD_WILDCARD=${SOURCE_FILE:0:-2}_token.txt
# Or stops if this file is created.
TOKEN_STOP=uftp_transfer_stop.txt

# Default values for the parallel file transfer
# (can be overwritten by environment variables or command line arguments)
PROCS=${PROCS:=4} # number of independent parallel transfers (each transfer uses $THREADS_PER_PROC threads)
THREADS_PER_PROC=${THREADS_PER_PROC:=8} # number of uftp threads per proc
STREAMS=${STREAMS:=3} # number of parallel uftp streams per thread
STEP=0
TOKEN_TRANSFER_FORWARD=${TOKEN_TRANSFER_FORWARD_WILDCARD/\*/$STEP}

echo "$SOURCE_FILE"


# loop that is infinite until $TOKEN_STOP is found
while true
do
    # check for the current ready token
    TOKEN_TRANSFER_FORWARD=${TOKEN_TRANSFER_FORWARD_WILDCARD/\*/$STEP}
    if test -f $TOKEN_TRANSFER_FORWARD; then
        # ./copy_sng_to_hawk.sh
        # try to copy the file

        starttime=`date +%s`
        for ((i=0; i<NUM_FILES; i++)); do
          SOURCE_FILE_INSTANCE=$(echo $TOKEN_TRANSFER_FORWARD)
          SOURCE_FILE_INSTANCE=${SOURCE_FILE_INSTANCE:0:-10}_0.part${i}
          echo instance "$SOURCE_FILE_INSTANCE"

          {
            # compress the file
            ${LZ4_PATH} -1 -f $SOURCE_FILE_INSTANCE $SOURCE_FILE_INSTANCE.lz4 # --rm TODO: ADD if you want to remove the original file
            # copy the compressed file
            uftp cp -n ${STREAMS} -t ${THREADS_PER_PROC} -i ${IDENTITYFILE} -u ${TARGET_USER} ${SOURCE_FILE_INSTANCE}.lz4 ${TARGET_URL}:${TARGET_PATH}
            # generate and copy the finish token file witch signals that the copie of this part is finished
            touch ${TOKEN_TRANSFER_FORWARD}.part${i}.lz4
            uftp cp -n ${STREAMS}  -i ${IDENTITYFILE} -u ${TARGET_USER} ${TOKEN_TRANSFER_FORWARD}.part${i}.lz4 ${TARGET_URL}:${TARGET_PATH}
            rm -f ${TOKEN_TRANSFER_FORWARD}.part${i}.lz4
          } &

        done
        wait

        # remove the start token
        rm -f ${TOKEN_TRANSFER_FORWARD}

        # print the needed time
        endtime=`date +%s`
        echo copied in `expr $endtime - $starttime` seconds.


        STEP=$((STEP+1))
    fi
    if test -f "$SOURCE_PATH/$TOKEN_STOP"; then
        break
    fi
    # wait for 1 second to reduce the load on the system
    sleep 1
done

