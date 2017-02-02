import sys, random, numpy
from parameters_scan import *

class random_grids:
    
    def __init__(self,par_in,nref_in,limdir,limit):
        global par
        par=par_in
        print "random distribution in %s" % par.scanpars         
        self.logdir=par.ddir+'/RDlogs'
        self.type='random_grid'
        self.nadapt=0
        self.serialvec=[]
        self.serialgrid=[]

        if nref_in == -1:
            self.nref=100
        else:
            self.nref=nref_in

        self.limit=limit
        if limit>0.:
            if os.path.exists(limdir+'/SGlogs'):
                self.limvals=analyze_SG(limdir)
            elif os.path.exists(limdir+'/DGlogs'):
                self.limvals=analyze_DG(limdir)
            else:
                print 'Directory for the interpolation of test point values is not valid'
                sys.exit()

    def create_initial_grid(self):
        #initial grid: first set of random parameter points

        self.serialgrid.extend(self.get_newpoints())
        self.npoints=len(self.serialgrid)

    def get_newpoints(self):
        newpoints=[]
        arr=numpy.zeros((3,1))
        while len(newpoints)<self.nref:
            newpoints.append([])
            for ind in range(par.scandim):
                rnum=random.uniform(par.scanrange[ind][0], par.scanrange[ind][1])
                newpoints[-1].append(rnum)
                arr[ind,0]=rnum
            if self.limit>0.:
                if self.limvals.ninterp(arr)[0]>self.limit:
                    newpoints.pop()

        return newpoints

    def reconstruct_grid(self):
        #initial grid: corners of scan interval
        self.npoints=len(self.serialgrid)
        return

    def get_parlist(self,lower,upper):
        #create parlist
        parlist=[]
        for i in range(lower,upper):
            parlist.append(self.serialgrid[i])

        return parlist

    def delete_old_res(self):
        self.serialvec=[]

    def add_res(self,reslist):
        self.serialvec.extend(reslist)

    def refine(self):
        self.serialgrid.extend(self.get_newpoints())
        self.npoints=len(self.serialgrid)
        self.nadapt+=1
