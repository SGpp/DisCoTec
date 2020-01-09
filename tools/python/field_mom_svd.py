import numpy as np
import matplotlib.pyplot as plt
from sys import path
import sys
#GENE_path="/afs/ipp/home/d/davidhat/gene0/tools/python"
#path.append(GENE_path)
from ParIO import Parameters
import os

#The SCALAPACK SVD routine is embedded in GENE.
#This script creates a GENE parameters file that has all the correct parameters
#set to perform an SVD analysis of the field and moment dat.  The routine in GENE 
#(in the file diag_fmsvd.F90)
#reads in phi,  and the species(1) Tpar and Tperp from the fields and mom_* files and
#produces an SVD decomposition of the data.  Each SVD mode is defined by a 3D mode
#structure for phi, Tpar, and Tperp.  As such, one can associate a heat flux with
#each SVD mode.  The mode structures are output in the directory 'diagdir/field_mom_svd',
# which is created by this script.  The files are in the standard GENE format (field_svd and mom_*_svd) 
#so that the mode structures can be visualized and analyzed using the GENE 
#diagnostic tool.  The mode number takes the place of the time coordinate in the
#output files.  Additionally, the SVD time modes are output in 'c_tot_svd.dat'
#Instructions:
#Modify the parameters below in the 'User Input' section and execute the script
#Output: a GENE parameters file with the appropriate parameters set for the SVD
#analysis


##########User Input########
##########User Input########
##########User Input########

#use every sparse_factor'th time point
sparse_factor = 1 
start_time = 100.0
data_path='/gpfs/ipp/davidhat/yglobal_svd/'
#example: file_suffixes=['_0001','_0002','_0003','_0004']
file_suffixes=['.dat']
#Total number of runs
nruns=len(file_suffixes)
nprocs=64
#Use swap_endian=True if the data you want to analyze has opposite endianness
#Warning: swap_endian may not work on all compilers
swap_endian=False

##########User Input########
##########User Input########
##########User Input########

ntimes=np.zeros(nruns,dtype=np.int32)
fsuf_out=''

for i in range(nruns):
    fsuf_out=fsuf_out+' \''+file_suffixes[i]+'\''

for i in range(nruns):
    par=Parameters()
    par.Read_Pars(data_path+'parameters'+file_suffixes[i])
    #print par['number of computed time steps']
    #print par.pardict
    #print par.pardict['step_time']
    print par.pardict['number of computed time steps']
    ntstep=float(par.pardict['number of computed time steps'])
    if ntstep==-1:
        sys.exit("Error!! Modify this script to account for parameters.dat files from unfinished runs!")
    else:
        ntimes[i]=int(ntstep/float(par.pardict['istep_field']))

print 'ntimes',ntimes 

mat_size=np.sum(ntimes)*float(par.pardict['nx0'])*float(par.pardict['nky0'])*float(par.pardict['nz0'])*3
mem_requirement=mat_size*16/1000000.0
print "Matrix is ", mem_requirement, " MB"
print mem_requirement/float(nprocs), " MB per processor (if start_time=0.0)."
#print mem_requirement/float(nprocs)/, " MB per processor."

#print "par.nmldict",par.nmldict

if 'extended_diags' in par.namelists:
    print "extended_diags is a namelist"
else:
    par.namelists.append('extended_diags')

if 'perf_vec' in par.pardict:
    par.pardict['perf_vec']='1 1 1 1 1 1 1 1 1'
else:
    par.pardict['perf_vec']='1 1 1 1 1 1 1 1 1'
    par.nmldict['perf_vec']='general'

if 'comp_type' in par.pardict:
    par.pardict['comp_type']="\'SV\'"
else:
    par.pardict['comp_type']="\'SV\'"
    par.nmldict['comp_type']='general'

par.pardict['SVD_field_mom']='T'
par.nmldict['SVD_field_mom']='extended_diags'
par.pardict['SVD_n_restarts']=str(nruns)
par.nmldict['SVD_n_restarts']='extended_diags'
par.pardict['SVD_df_n_time']=str(ntimes)[1:-1]
par.nmldict['SVD_df_n_time']='extended_diags'
par.pardict['SVD_start_time']=str(start_time)
par.nmldict['SVD_start_time']='extended_diags'
par.pardict['SVD_sparse_factor']=str(sparse_factor)
par.nmldict['SVD_sparse_factor']='extended_diags'
par.pardict['SVD_df_file_path']='\''+data_path+'\''
par.nmldict['SVD_df_file_path']='extended_diags'
par.pardict['SVD_df_file_suffix']=fsuf_out
par.nmldict['SVD_df_file_suffix']='extended_diags'

os.system('mkdir '+data_path+'/field_mom_svd/')
par.pardict['diagdir']='\''+data_path+'/field_mom_svd/'+'\''
print "Double check processor assignment in parameters file!"

par.pardict['n_procs_sim']=str(nprocs)
par.pardict['n_procs_s']=str(1)
par.pardict['n_procs_x']=str(1)
par.pardict['n_procs_y']=str(1)
par.pardict['n_procs_z']=str(1)
par.pardict['n_procs_v']=str(1)
par.pardict['n_procs_w']=str(min(int(float(par.pardict['nw0'])),nprocs))
procs_left=nprocs/int(float(par.pardict['n_procs_w']))
if procs_left > 1:
    par.pardict['n_procs_v']=str(min(procs_left,int(float(par.pardict['nv0']))/2))
    procs_left=procs_left/int(float(par.pardict['n_procs_v']))
if procs_left > 1:
    par.pardict['n_procs_z']=str(min(procs_left,int(float(par.pardict['nz0']))/2))
    procs_left=procs_left/int(float(par.pardict['n_procs_z']))
if procs_left > 1:
    print "Distribute your own processors!"


print par.pardict['n_procs_s']
print par.pardict['n_procs_x']
print par.pardict['n_procs_y']
print par.pardict['n_procs_z']
print par.pardict['n_procs_v']
print par.pardict['n_procs_w']

par.Write_Pars('parameters')







