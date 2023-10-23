#!/bin/bash

# exit on error
set -e

# kill all processes from this script if the script is terminated
trap "trap - SIGTERM && kill -- -$$" SIGINT SIGTERM EXIT

# Paths and filenames with default values
TARGET_PATH=${TARGET_PATH:=//p/scratch/widediscotecjsc/widely/widely_measurements/widely_14x512/}
number_of_files=${number_of_files:=8}
DSG_NUM=${DSG_NUM:=0}
LZ4_PATH=${LZ4_PATH:=/p/scratch/widediscotecjsc/lz4yume/lz4}

cd $TARGET_PATH


while true; do
  # check for current ready tokens
  for i in dsg_${DSG_NUM}*_token.txt.part*.lz4; do
    [ -f "$i" ] || break
    echo $i
    # remove the token and start asynchronus decompression
    rm $i
    {
    ${LZ4_PATH} -f ${i/token.txt/0}
    # generate ready token for the decompressed part
    touch ${i/.lz4/}
    } &
    sleep 1
  done


  for STEP in {0..10}; do
    count=0
    # check for the current ready token
    if [ -f "dsg_${DSG_NUM}_step${STEP}_token.txt" ]; then
      # this step was already done
      continue
    fi
    # count the number of decrompressed files
    for i in dsg_${DSG_NUM}*${STEP}_token.txt.part*; do
      [ -f "$i" ] || break
      if ! [[ $i =~ ".lz4" ]]; then
        count=$((count+1))
      fi
    done
    # if all files are decompressed, generate the global ready token and cleanup
    if [ $count -eq $number_of_files ]; then
      touch  dsg_${DSG_NUM}_step${STEP}_token.txt
      for file in ${STEP}_token.txt.part*.lz4; do
        [ -f "$file" ] || break
        rm $file
      done
      for file in dsg_${DSG_NUM}*step${STEP}_0.part*.lz4; do
        [ -f "$file" ] || break
        rm $file
      done
    fi
  done
  sleep 1
done

