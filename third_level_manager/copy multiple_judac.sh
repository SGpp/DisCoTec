#!/bin/bash
# execute this script on skx-arch.supermuc.lrz.de

module load uftp-client

# exit on error
#set -e #-x

# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# Paths and filenames with default values
# (can be overwritten by environment variables or command line arguments)
SOURCE_PATH=${SOURCE_PATH:=/p/scratch/widediscotecjsc/3-systems/measurements/3-systems_3x8}
DSG_NUM=${DSG_NUM:=1}
SOURCE_FILE=${SOURCE_FILE:=${SOURCE_PATH}/dsg_${DSG_NUM}_step*_0}

TARGET_PATH_NG=${TARGET_PATH_NG:=//hppfs/work/pn36xu/di39qun2/3-systems/measurements/3-systems_8x8}
TARGET_URL_NG=https://datagw03.supermuc.lrz.de:9000/rest/auth/DATAGW
TARGET_USER_NG=${TARGET_USER_NG:=di39qun2}
IDENTITYFILE_NG=${IDENTITYFILE_NG:=~/.uftp/id_uftp_to_lrz}

TARGET_PATH_HLRS=${TARGET_PATH_HLRS:=//lustre/hpe/ws11/ws11.1/ws/ipvpolli-discotec_weak_scaling/3-systems/measurements/3-system_3x8}
TARGET_URL_HLRS=https://gridftp-fr1.hww.hlrs.de:9000/rest/auth/HLRS
TARGET_USER_HLRS=${TARGET_USER_HLRS:=ipvpolli}
IDENTITYFILE_HLRS=${IDENTITYFILE_HLRS:=~/.uftp/id_uftp_to_hlrs}
DSG_NUM_HAWK=${DSG_NUM_HAWK:=0}
RECIEVE_FILE=${RECIEVE_FILE:=${TARGET_PATH_HLRS}/dsg_${DSG_NUM_HAWK}_step*_0}

NUM_FILES=${NUM_FILES:=8}
HOST_NR_NG=$(( $(hostname | tr -dc '0-9') % 5 - 1 ))
NUM_HOSTS=4

# Tokens for the file transfer. If the file is not there, the script will wait until it is there.
TOKEN_TRANSFER_FORWARD_WILDCARD=${SOURCE_FILE:0:-1}token.txt
TOKEN_TRANSFER_BACKWARD_WILDCARD=${RECIEVE_FILE:0:-1}token.txt
# Or stops if this file is created.
TOKEN_STOP=uftp_transfer_stop.txt

# Default values for the parallel file transfer
# (can be overwritten by environment variables or command line arguments)

THREADS_PER_PROC=${THREADS_PER_PROC:=8} # number of uftp threads per proc
STREAMS=${STREAMS:=3} # number of parallel uftp streams per thread
STEP_SEND=0
STEP_RECEIVE=0

LOCAL_SEQUENCE=$(seq $HOST_NR_NG $NUM_HOSTS $((NUM_FILES-1)))

echo File Pattern: "$SOURCE_FILE"
echo Target parts: $LOCAL_SEQUENCE
NUM_LOCAL_FILES=$(echo $LOCAL_SEQUENCE | wc -w)
copie=0
receive=0




# loop that is infinite until $TOKEN_STOP is found
while true
do
    for part in $LOCAL_SEQUENCE; do
      # send
      # check for the current ready token
      TOKEN_TRANSFER_FORWARD=${TOKEN_TRANSFER_FORWARD_WILDCARD/\*/$STEP_SEND}
      TOKEN_TRANSFER_FORWARD=${TOKEN_TRANSFER_FORWARD}.part${part}
      if test -f $TOKEN_TRANSFER_FORWARD; then
          rm $TOKEN_TRANSFER_FORWARD
          SOURCE_FILE_INSTANCE=${TOKEN_TRANSFER_FORWARD/token.txt/0}
          copie=$((copie+1))

          {
            starttime=`date +%s`
            echo "$(date +%T) start copying $(basename $SOURCE_FILE_INSTANCE) to NG as copie $copie / $NUM_LOCAL_FILES"
            # copy the compressed file
            uftp cp -n ${STREAMS} -t ${THREADS_PER_PROC} -i ${IDENTITYFILE_NG} -u ${TARGET_USER_NG} ${SOURCE_FILE_INSTANCE} ${TARGET_URL_NG}:${TARGET_PATH_NG}
            # generate TOKEN_TRANSFER_FORWARD on remote system
            uftp cp -i ${IDENTITYFILE_NG} -u ${TARGET_USER_NG} -B0 /dev/zero ${TARGET_URL_NG}:${TARGET_PATH_NG}/$(basename ${TOKEN_TRANSFER_FORWARD})
            endtime=`date +%s`
            echo $(date +%T) copied $(basename $SOURCE_FILE_INSTANCE) in $((endtime-starttime)) seconds to NG.
          } &
          {
            starttime=`date +%s`
            echo "$(date +%T) start copying $(basename $SOURCE_FILE_INSTANCE) to HAWK as copie $copie / $NUM_LOCAL_FILES"
            # copy the compressed file
            uftp cp -n ${STREAMS} -t ${THREADS_PER_PROC} -i ${IDENTITYFILE_HLRS} -u ${TARGET_USER_HLRS} ${SOURCE_FILE_INSTANCE} ${TARGET_URL_HLRS}:${TARGET_PATH_HLRS}
            # generate TOKEN_TRANSFER_FORWARD on remote system
            uftp cp  -i ${IDENTITYFILE_HLRS} -u ${TARGET_USER_HLRS} -B0 /dev/zero ${TARGET_URL_HLRS}:${TARGET_PATH_HLRS}/$(basename ${TOKEN_TRANSFER_FORWARD})
            endtime=`date +%s`
            echo $(date +%T) copied $(basename $SOURCE_FILE_INSTANCE) in $((endtime-starttime)) seconds to HAWK.
          } &
      fi
      # receive
      TOKEN_TRANSFER_BACKWARD=${TOKEN_TRANSFER_BACKWARD_WILDCARD/\*/$STEP_RECEIVE}
      TOKEN_TRANSFER_BACKWARD=${TOKEN_TRANSFER_BACKWARD}.part${part}
      # if [ ! -f ${SOURCE_PATH}/$(basename $TOKEN_TRANSFER_BACKWARD) ]; then
        uftp cp -i ${IDENTITYFILE_HLRS} -u ${TARGET_USER_HLRS} $TARGET_URL_HLRS:$TOKEN_TRANSFER_BACKWARD /dev/null > /dev/null 2>&1
        cp_status=$?
        if (exit $cp_status); then
          uftp rm -q -i ${IDENTITYFILE_HLRS} -u ${TARGET_USER_HLRS} $TARGET_URL_HLRS:$TOKEN_TRANSFER_BACKWARD
          SOURCE_FILE_INSTANCE=${TOKEN_TRANSFER_BACKWARD/token.txt/0}
          receive=$((receive+1))
          {
              starttime=`date +%s`
              echo "$(date +%T) start copying $(basename $SOURCE_FILE_INSTANCE) from HAWK as copie $copie / $NUM_LOCAL_FILES"
              # copy the compressed file
              uftp cp -n ${STREAMS} -t ${THREADS_PER_PROC} -i ${IDENTITYFILE_HLRS} -u ${TARGET_USER_HLRS} ${TARGET_URL_HLRS}:${SOURCE_FILE_INSTANCE} ${SOURCE_PATH}
              # generate TOKEN_TRANSFER_BACKWARD
              touch ${SOURCE_PATH}/$(basename $TOKEN_TRANSFER_BACKWARD)
              endtime=`date +%s`
              echo $(date +%T) copied $(basename $SOURCE_FILE_INSTANCE) in $((endtime-starttime)) seconds from HAWK.
          } &

        fi
      # fi
    done


    if [ $copie -eq $NUM_LOCAL_FILES ]; then
      copie=0
      STEP_SEND=$((STEP_SEND+1))
    fi

    if [ $receive -eq $NUM_LOCAL_FILES ]; then
      receive=0
      STEP_RECEIVE=$((STEP_RECEIVE+1))
    fi

    if test -f "$SOURCE_PATH/$TOKEN_STOP"; then
        break
    fi
    # wait for 1 second to reduce the load on the system
    sleep 1
done

