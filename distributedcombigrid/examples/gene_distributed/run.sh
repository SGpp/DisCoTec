#!/bin/bash
#make
cp ginstance/loaddata_* template/
rm -r ginstance*
rm out/*
python3 preproc.py
cd ginstance
echo "offset 0\n" > offset.txt
source start.bat
cd ..
