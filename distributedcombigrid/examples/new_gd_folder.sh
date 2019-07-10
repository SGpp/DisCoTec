#!/bin/bash

if [ -z "$1" ]; then
    echo "Argument: name of the example directory to be created"
    exit 0
fi

mkdir $1

cd $1

export CODE_DIR=../gene_distributed/

ln -s $CODE_DIR/src 
ln -s $CODE_DIR/lib
ln -s $CODE_DIR/manager
ln -s $CODE_DIR/setenv.sh
ln -s $CODE_DIR/preproc.py
ln -s $CODE_DIR/updateOffset.py

# these things can and should be adapted, especially ctparam and the batch file
cp -r $CODE_DIR/template .
cp $CODE_DIR/run.sh .
cp $CODE_DIR/ctparam .
cp $CODE_DIR/batch.pbs .
#etc

cd -
