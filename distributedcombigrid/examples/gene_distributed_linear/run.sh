#!/bin/bash
#make
rm -r ginstance*
rm out/*
python preproc.py
cd ginstance
source start.bat
cd ..
