#!/usr/bin/env python
from sys import exit,path
path.append('..')
from ParIO import Parameters
from fieldlib import *
from momlib import *
from os.path import join
from common import common_data
from averages import *
from geom import *
from fluxes import *

##### INPUT #####
time1=0 #beginning and end of time window
time2=50
rundir='' #directory containing simulation data
diags=['hfspec'] #which diags to use (currently there's only a half one)
need_moms= (set(['hfspec']) & set(diags))!=None  #this is supposed to check if one of the diags that will be done requires mom data ('diags' is compared to a fixed set of such diagnostics, which consists only of hfspec atm.))
runs=['dat',1,2] #runs to be examined, can also be given as, e.g. range(1,5)
spec='Electrons'

##### FIND RUNS OF INTEREST (ROI) #####
roi=[]
mins=[0]
times=[]
oldmax=0.
for run in runs:
    #read parfile for each run
    par=Parameters()
    fileext='_%s' %(str(run))
    if run in ['dat','act']: fileext='.dat'
    par.Read_Pars(join(rundir,'./parameters'+fileext))
    pars=par.pardict
    #istep_field and istep_mom synchronous?
    isf=int(pars['istep_field'])
    ism=int(pars['istep_mom'])
    if isf!=ism: exit('Error: Differing istep_field and istep_mom currently not supported!')
    
    #is the current run in the selected time window?
    curr_fld=fieldfile(join(rundir,'field'+fileext),pars)
    min,max=curr_fld.get_minmaxtime()  
    if min>time2:
        break
    #if time windows overlap, use the data from the continuation run
    if oldmax>min:
        for i in range(len(times)):
            if times[i]>=min: 
                del times[i:]
                break
    oldmax=max
    mins.append(min)
    #generate array of all field/mom times within the time window
    if not (time2<min or time1>max): 
        roi.append(run)
        times.extend([time for time in curr_fld.tfld if time>=time1 and time<=time2])

##### INIT DERIVED VALUES #####
#reset to first run
run=roi[0]
common=common_data(pars)
common.rundir=rundir
common.fileext='_%s' %(str(run))
if run in ['dat','act']: fileext='.dat'
par.Read_Pars(join(rundir,'./parameters'+fileext))
count=0

#set file handles for field and mom files and initialize their classes
fieldfilename=join(rundir,'field'+common.fileext)
momfilename=join(rundir,'mom_%s'%(spec)+common.fileext)

curr_fld=fieldfile(fieldfilename,pars)
if need_moms: curr_mom=momfile(momfilename,pars)

#make time array available to common module
common.tlen=len(times)
common.times=np.array(times,dtype=np.float64)
common.set_profiles(spec)

##### INIT DIAGS #####
#instantiate geometry and averaging classes
geom=geometry(common,pars)
avg=averages(common,geom)

#rudimentary heat flux module
fl=heat_fluxes(common,curr_fld,curr_mom,avg,pars)


##### DIAGNOSTIC LOOP #####
for time in times:
    #time stepping
    print 'Time: %1.10g' %(time)
    common.tind=times.index(time)

    #set file pointer to next file at the appropriate time
    #'mins' contains the starting time for each run
    if time in mins[1:] and len(roi)>1:
        count+=1
        run=roi[count]

        common.fileext='_%s' %(str(run))
        if run in ['dat','act']: common.fileext='.dat'

        common.set_profiles(spec)

        fieldfilename=join(rundir,'field'+common.fileext)
        momfilename=join(rundir,'mom_%s'%(spec)+common.fileext)

        #redirect makes the field and mom classes read a new file
        #it's done that way since when just reinitializing the classes,
        #all diags that were initialized before still try to work on
        #the old files
        curr_fld.redirect(fieldfilename)
        if need_moms:
            curr_mom.redirect(momfilename)

    #tell the field/mom classes which timestep to look for
    curr_fld.set_time(time)
    if need_moms: curr_mom.set_time(time)

    #run a diag
    fl.high_ky_fraction()

#run plotting
fl.plot()    
    
    
