import os, sys, pickle, shutil, numpy

#sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/sg++_p2.6/release_090/lib/pysgpp')
#sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/release_090/lib/pysgpp')
sys.path.append('/lustre/jhome6/fsngeul/fsgeul01/soft/release_090/lib/pysgpp')
from pysgpp import *

class analyze_SG:
    def __init__(self,path,sgo_num=-1):
        self.path=path+'/SGlogs'
        self.logname='sparse_grid.obj'
        #if no number is given, use the latest (=most precise) sparse_grid.obj

        refnum=-1
        for file in os.listdir(self.path):
            fs=file.split('_')
            if fs[1]=='grid.obj':
                refnum= max(refnum, int(fs[2]))
        refnum=refnum+1

        if refnum==0:
            sys.exit('\nERROR: no %s in path!\n'%self.logname)

        self.nstages=refnum

        if sgo_num==-1:
            sgo_num=self.nstages-1

        sfile = open('%s/%s_%d' % (self.path,self.logname,sgo_num), 'rb')
        self.scanpars=pickle.load(sfile)
        self.scanrange=pickle.load(sfile)
        self.serialgrid=pickle.load(sfile)
        self.serialvec=pickle.load(sfile)

        self.points_per_stage=pickle.load(sfile)

        sfile.close()
        #convert back to SG++
        self.alpha=DataVector(len(self.serialvec))
        for i in range(len(self.serialvec)):
            self.alpha[i]=self.serialvec[i]
        self.grid = Grid.unserialize(self.serialgrid)
        self.ndim=self.grid.getStorage().dim()
      
        self.interpol=self.grid.createOperationEval()

        self.reflev=-1
       
#        print "Dimensionality of the scan:   %d" % self.ndim
#        print "Number of refinement levels:  %d" % refnum


    def ninterp(self,pvals,reflev_in=-1):
 
        if reflev_in==-1:
            reflev_in=self.nstages-1
        coordshape=pvals.shape

        dm = DataMatrix(coordshape[1],coordshape[0])
        coordvec=numpy.zeros(coordshape[1])       
        for ind2 in range(self.ndim):
            coordvec[:]=(pvals[ind2,:]-self.scanrange[ind2][0])/self.scanrange[ind2][2]
            for ind1 in range(coordshape[1]):
                dm.set(ind1,ind2,coordvec[ind1])

        if self.reflev != reflev_in:
            #re-initialize reduced alpha
            self.redalpha = DataVector(0)
            for i in range(self.points_per_stage[reflev_in]):
                self.redalpha.append(self.alpha[i])
            for i in range(self.points_per_stage[reflev_in],self.points_per_stage[-1]):
                self.redalpha.append(0.)
            self.reflev=reflev_in
     
        dv = DataVector(coordshape[1])
        self.grid.createOperationB().multTranspose(self.redalpha, dm, dv)

        return dv


