import numpy as np
import matplotlib.pyplot as plt
from ParIO import *
from grids import *
import os
from energylib import *


suffix=raw_input("Enter suffix for parameters file (e.g., .dat or _0001:\n")

#Select directory and make link to parameters.dat file
if os.path.isfile('./parameters'+suffix):
    par0=Parameters()
    par0.Read_Pars('./parameters'+suffix)
    print 'Using diagdir:',par0.pardict['diagdir']
    cd=raw_input("Change directory?(y/n):\n")
    if cd=='y':
        valid_dir=False
        while not valid_dir:
            diagdir=raw_input("Enter new directory:\n")
            if os.path.isfile(diagdir+'/parameters'+suffix):
                os.system('rm parameters'+suffix)
                os.system('ln -s '+diagdir+'/parameters'+suffix)
                valid_dir=True
            else:
                print "No parameters file in directory:", diagdir
                diagdir=raw_input("Enter new directory:\n")
    else:
        pass
else:
    valid_dir=False
    while not valid_dir:
        diagdir=raw_input("Enter diag directory:\n")
        if os.path.isfile(diagdir+'/parameters'+suffix):
            os.system('ln -s '+diagdir+'/parameters'+suffix)
            valid_dir=True
        else:
            print "No parameters file in directory:", diagdir+parfile
            diagdir=raw_input("Enter diag directory:\n")


par0.Read_Pars('./parameters'+suffix)
par=par0.pardict

#cm=common_data(par)



