import os, sys, pickle, shutil, copy, numpy
try:
    from scipy import ndimage, reshape, mgrid, ogrid, zeros
    with_scipy= True
except ImportError:
    with_scipy = False

class analyze_RD:
    def __init__(self,path,sgo_num=-1):
        self.path=path+'/RDlogs'
        self.logname='random_grid.obj'
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
       


