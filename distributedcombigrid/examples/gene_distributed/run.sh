#!/bin/bash
#make
cp ginstance/l_* template/
rm -r ginstance*
rm out/*
python3 preproc.py
cd ginstance
echo "offset 0" > offset.txt
source start.bat
cd ..
