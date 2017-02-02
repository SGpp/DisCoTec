#!/usr/bin/env python
import os, sys, pickle, shutil
from optparse import OptionParser

from analyze_dir import *
from parameters_scan import *

from dense_grids import *
from random_grids import *
try:
    from sparse_grids import *
except:
    print 'SG++ not found'

#input arguments
parser = OptionParser(usage='SG++ interface for GENE')
parser.add_option("-t", "--target", dest="target", default=1, type='int',
                  help="Target function for the reanalyzed data (0: QL model, 1: numerator, 2: weight, 3: denominator, 4: difference of eigenvalues)")
(opts, args) = parser.parse_args()

oldpath=args[0]
newpath=args[1]
#read parameters file
par=parameters_scan()
par.read_scanpars('%s/parameters'%oldpath)

if os.path.exists(oldpath+'/DGlogs'):
    print 'reanalyzing dense grid'
    mode='DGlogs'
elif os.path.exists(oldpath+'/SGlogs'):
    print 'reanalyzing sparse grid'
    mode='SGlogs'
else:
    print 'reanalyzing test sample'
    mode='RDlogs'

if mode=='RDlogsa':
    print 'da'

else:
    if mode=='DGlogs':
        grid=dense_grids(par)
    elif mode=='SGlogs':
        grid=sparse_grids(par,1,0.)
    elif mode=='RDlogs':
        grid=random_grids(par,1)

    grid.nadapt=-1
    for logfile in os.listdir('%s/%s' % (oldpath,mode)):
        if 'obj' in logfile:
            grid.nadapt=max(grid.nadapt,int(logfile.split('obj_')[1]))

    if grid.nadapt>=0:
        sfile = open('%s/%s/%s.obj_%d' % (oldpath,mode,grid.type,grid.nadapt), 'rb')
        scanpars=pickle.load(sfile)
        scanrange=pickle.load(sfile)
        grid.serialgrid=pickle.load(sfile)
        grid.serialvec=pickle.load(sfile)
        npointslist=pickle.load(sfile)
        sfile.close()
    elif mode=='RDlogs':
        qlfile=open(oldpath+'/quasilinear.log','r')
        head=qlfile.readline()
        parlist=[]
        for line in qlfile:
            spline=line.split('|')
            parlist.append([])
            for ind in range(par.scandim):
                parlist[-1].append(float(spline[1+ind]))
        grid.serialgrid=parlist
        npointslist=[len(parlist)]
        grid.nadapt=0
    else:
        print 'no valid data in %s' % oldpath
        sys.exit()
    
    grid.reconstruct_grid()
    noldpoints=npointslist[-1]

    par.ddir=newpath

    newres= analyze_dir(1000,False,par=par)
   

    bounds=npointslist[:]
    bounds.insert(0,0)
    newreslist=[]
    for ref in range(grid.nadapt+1):
        newreslist.extend(newres.ql( '%s/genefiles/ref_%d'%(oldpath,ref),grid.get_parlist(bounds[ref],bounds[ref+1]),bounds[ref],opts.target))

    grid.delete_old_res()
    grid.add_res(newreslist)

    newlogdir=newpath+'/'+mode

    if os.path.exists(newlogdir):
        shutil.rmtree(newlogdir)
    os.mkdir(newlogdir)

    #save status
    sfile = open('%s/%s.obj_%d' % (newlogdir,grid.type,grid.nadapt), 'wb')

    pickle.dump(par.scanpars,sfile)
    pickle.dump(par.scanrange,sfile)
    pickle.dump(grid.serialgrid, sfile)
    pickle.dump(grid.serialvec,sfile)
    pickle.dump(npointslist,sfile)
    sfile.close()
        
