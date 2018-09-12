#!/bin/bash
#make
rm -r ginstance*
rm out/*
python3 preproc.py
cd ginstance
source start.bat
cd ..
