#!/usr/bin/env python
import os, sys, math, pickle
from optparse import OptionParser
from math import ceil

from parameters_scan import *
from rungene import *
#from runtest import *

from dense_grids import *
from random_grids import *
try:
    from sparse_grids import *
except:
    print 'SG++ not found'

#read parameters file
par=parameters_scan()
par.read_scanpars('./parameters')

if par.scandim==0:
    sys.exit('\nERROR: no scan defined in parameters file!\n')

#compute total number of processors
procnum=int(par.pardict['n_procs_sim'])*int(par.pardict['n_parallel_sims'])

#input arguments
parser = OptionParser(usage='SG++ interface for GENE')
parser.add_option("-s", "--syscall", dest="syscall", type='string',
                  default='gmake -f ../makefile run N_PES='+str(procnum),
                  help="System call to start scangene")
parser.add_option("-m", "--mode", dest="mode", default='sparse', type='string',
                  help="Sets mode of operation. Possible values: dense, sparse, testsample")
parser.add_option("-p", "--points", dest="npoint_in", default=1e40, type='int', metavar='NP',
                  help="For --mode=sparse, the grid is refined while until the number of gridpoints reaches NP. For --mode=testsample, NP is the exact number of sample points")
parser.add_option("-e", "--errorlimit", dest="errlim", default=1e-40, type='float',
                  help="Grid is refined until the error is less than ERRLIM")
parser.add_option("-n", "--nrefpoints", dest="nref", type='int',
                  default=int(ceil(float(par.pardict['n_parallel_sims'].strip("'"))/par.scandim)),
                  help="Number of grid points that are refined per stage")
parser.add_option("-c", "--continue", dest="cont", action="store_true",default=False,
                  help="Continue an old scan")
parser.add_option("-l", "--limit", dest="limit", default=0., type='float',
                  help="Only points with estimated value of less than limit are used for testsample")
parser.add_option("-d", "--directory", dest="limdir", type='string',
                  default='',help="Directory containing the scan that is used for estimating the value of the sample points")
(opts, args) = parser.parse_args()

if not opts.cont:
#delete old status/log files 
    for entry in os.listdir(par.ddir):
        if os.path.isfile(par.ddir+'/'+entry):
            os.remove(par.ddir+'/'+entry)

#initialize gene
gene=rungene(par,opts.syscall,opts.cont)
#gene=runtest(par,opts.syscall,opts.cont)

if opts.mode=='rectest':
    #only to test speedup due to old checkpoints
    print "testing subspace recycling" 
    parlist=[]
    center=[]
    gene.chpt.rectest=True
    for ind in range(par.scandim):
        center.append(0.5*(par.scanrange[ind][0]+par.scanrange[ind][1]))
    for isamp in range(opts.npoint_in):
        parlist.append(center)
    reslist=gene.comp(parlist)
else:
    #do a hierarchical parameters scan
    if opts.mode=='dense':
        grid=dense_grids(par)
    elif opts.mode=='sparse':
        grid=sparse_grids(par,opts.nref,opts.errlim)
    elif opts.mode=='testsample':
        grid=random_grids(par,opts.nref,opts.limdir,opts.limit)

    if opts.cont:
        #continuation of old scan
        grid.nadapt=-1
        for logfile in os.listdir(grid.logdir):
            if 'obj' in logfile:
                 grid.nadapt=max(grid.nadapt,int(logfile.split('obj_')[1]))

        sfile = open('%s/%s.obj_%d' % (grid.logdir,grid.type,grid.nadapt), 'rb')
        oldscanpars=pickle.load(sfile)
        oldscanrange=pickle.load(sfile)
        grid.serialgrid=pickle.load(sfile)
        grid.serialvec=pickle.load(sfile)
        npointslist=pickle.load(sfile)
        sfile.close()
        if oldscanpars != par.scanpars:
            print 'scan parameters of old and new parameter file do not agree, continuation not possible'
            sys.exit()
        if oldscanrange != par.scanrange:
            print 'scan range of old and new parameter file do not agree, continuation not possible'
            sys.exit()
        grid.reconstruct_grid()
        if par.use_checkpoints:
            gene.chpt.update(grid.get_parlist(0,grid.npoints))
            for points in npointslist:
                gene.chpt.init_metric(points)
            while os.path.exists('%s/%dk' % (gene.chptfiles,gene.kmax+1)):
                gene.kmax+=1
        noldpoints=npointslist[-1]
        gene.noldpars=noldpoints
        gene.refnum=grid.nadapt+1
        if not par.external_file: 
            while os.path.exists('%s/%dk' % (gene.genefiles,gene.kmax+1)):
                gene.kmax+=1
        grid.refine()

    else:
    #old logfiles
        if os.path.exists(par.ddir+'/DGlogs'):
            shutil.rmtree(par.ddir+'/DGlogs')
        if os.path.exists(par.ddir+'/SGlogs'):
            shutil.rmtree(par.ddir+'/SGlogs')
        if os.path.exists(par.ddir+'/RDlogs'):
            shutil.rmtree(par.ddir+'/RDlogs')
        os.mkdir(grid.logdir)
        open('%s/start_tag' % grid.logdir, 'w').close()
        if os.path.exists(gene.done):
            shutil.rmtree(gene.done)
        grid.create_initial_grid()
        noldpoints=0
        npointslist=[]

    while (grid.npoints!=noldpoints and grid.npoints<=opts.npoint_in):
        npointslist.append(grid.npoints)
        print "refinement level %d, number of new grid points %d" % (grid.nadapt, grid.npoints-noldpoints)
        parlist=grid.get_parlist(noldpoints,grid.npoints)
       
        #compute results
        reslist= gene.comp(parlist)

        #update number of computed points 
        noldpoints=grid.npoints

        grid.add_res(reslist)

        #save status
        sfile = open('%s/%s.obj_%d' % (grid.logdir,grid.type,grid.nadapt), 'wb')

        pickle.dump(par.scanpars,sfile)
        pickle.dump(par.scanrange,sfile)
        pickle.dump(grid.serialgrid, sfile)
        pickle.dump(grid.serialvec,sfile)
        pickle.dump(npointslist,sfile)
        sfile.close()
        
        #initialize metric for checkpoints
        if par.use_checkpoints:
            gene.chpt.init_metric(grid.npoints)

        grid.refine()
        
