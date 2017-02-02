#!/usr/bin/env python
from os import listdir,chdir,curdir
import os.path as osp
#from ParIO import *

####### This script sums up all the CPU-hours spent on the present
####### directory and its subdirectories. 
####### USAGE: Start script from the directory to be analyzed. 
####### (the script itself can be placed at an arbitrary location)

#create a list of lines from a given parameters file
def prep_parfile(file):
    f=open(file)
    lines=f.readlines()
    return lines

#return the line containing a given parameter from the parfile list
def get_par(lines,par):
    for line in lines:
        if par in line:
            return line

#Recursively examine all parts of the present directory
def traverse(dir):
    misscount=0
    cpuh=0
    numjobs=0
    ls=listdir(dir)
    for entry in ls:
        #if entry is directory, traverse it
        if osp.isdir(entry):
            chdir(entry)
            print 'Entering ',entry
            dcpuh,dmisscount,dnumjobs=traverse(osp.curdir)
            cpuh+=dcpuh;misscount+=dmisscount;numjobs+=dnumjobs
            chdir('..')
        #else, if it is a parameters file, extract the CPU-hours (only if nonlinear)
        elif 'parameters' in entry:
            lines=prep_parfile(entry)
            try: 
                if get_par(lines,'nonlinear').split()[2]=='T':
                   try:
                       np=int(get_par(lines,'n_procs_sim').split()[2])
                       tim=float(get_par(lines,'time for initial value solver').split()[6])
                       #the following line raises an exception if tim==0, 
                       #so that the run is counted as 'missing the necessary data'
                       1./tim
                       cpuh+=np*tim/3600
                       numjobs+=1
                       print '%s: n_procs: %d, Runtime: %7.1f s, CPUh (cumulative per directory): %10.1f' %(entry,np,tim,cpuh)
                   except:
                       print 'missing data in %s' %entry
                       misscount+=1
            except:
                pass
    return cpuh,misscount,numjobs 

#main part
sum=0;misscount=0;numjobs=0
dsum,dmisscount,dnumjobs=traverse('.')
sum+=dsum;misscount+=dmisscount;numjobs+=dnumjobs
print '\nTotal CPU-hours spent on present subdir: %14.6g' %(sum)
print 'Data was extracted from %d successful jobs.' %(numjobs)
if misscount!=0:
    print 'Info: %d runs did not contain the necessary data to estimate their CPUh consumption!' %(misscount)
