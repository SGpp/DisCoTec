#!/bin/bash
#make
cp ginstance/loaddata_* template/
rm -r ginstance*
rm out/*
python3 preproc.py
cd ginstance
source start.bat
cd ..
