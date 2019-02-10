#! /bin/bash
if [ $# -eq 0 ]
then
	echo "Usage: forward.sh fromPort toPort"
	exit 0
fi

mkfifo backpipe
nc -l $1 0<backpipe | nc localhost $2 1>backpipe
