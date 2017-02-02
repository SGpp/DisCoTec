import os, sys, pickle, shutil, copy, numpy
try:
    from scipy import ndimage, reshape, mgrid, ogrid, zeros
    with_scipy= True
except ImportError:
    with_scipy = False

class analyze_DG:
    def __init__(self,path,sgo_num=-1):
        self.path=path+'/DGlogs'
        self.logname='dense_grid.obj'
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

        self.ndim=len(self.scanpars)

        self.reflev=-1
       
#        print "Dimensionality of the scan:   %d" % self.ndim
#        print "Number of refinement levels:  %d" % refnum


    def ninterp(self,pvals,reflev_in=-1):
        if not with_scipy:
            sys.exit('Scipy is not available on this system, dense grid interpolation not possible')
        if reflev_in==-1:
            reflev_in=self.nstages-1
        
        #transform pval
        points_per_dim=2**reflev_in+1     
        coordarr=numpy.zeros(pvals.shape)       
        for ind in range(self.ndim):
            #convert to [0,points_per_dim-1]
            coordarr[ind,:]=(pvals[ind,:]-self.scanrange[ind][0])/self.scanrange[ind][2]*(points_per_dim-1)
            
        #transform data
        datshape=[]
        for int in range(self.ndim):
            datshape.append(points_per_dim)
        if self.reflev != reflev_in:
            #re-initialize data array
            self.datarr=zeros(datshape)
            for ind in range(self.points_per_stage[reflev_in]):
                #the coordinates have to be rescaled for reflev_in < maximum reflev
                coord=[]
                for dim in range(self.ndim):
                    coord.append([self.serialgrid[ind][dim]/2**(self.nstages-reflev_in-1)])
                if (self.serialvec[ind]==[]):
                    print ind
                    print self.serialvec[ind]
                    self.datarr[coord]=0.
                    print 'WARNING: missing data point has been set to zero'
                    sys.exit()
                else:
                    self.datarr[coord]=self.serialvec[ind]
            self.reflev=reflev_in

        points=ndimage.map_coordinates(self.datarr, coordarr, order=1) 
        
        return points
