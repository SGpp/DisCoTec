import sys, math

from parameters_scan import *
sys.path.append('/lustre/jhome6/fsngeul/fsgeul01/soft/release_090/lib/pysgpp')
sys.path.append('/lustre/jhome6/fsngeul/fsgeul01/soft/release_090/bin')
#sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/sg++_p2.6/release_090/lib/pysgpp')
#sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/sg++_p2.6/release_090/bin')

from pysgpp import *
import tools

class sparse_grids():
    
    def __init__(self,par_in,nref_in,errlim_in):
        global par
        par=par_in

        if nref_in == -1:
            self.nref=int(ceil(float(par.pardict['n_parallel_sims'].strip("'"))/par.scandim))
        else:
            self.nref=nref_in

        self.errlim=errlim_in
        print "sparse scan in %s" % par.scanpars 
        
        self.resvec = DataVector(0)

        #bbox= grid.getBoundingBox()
        #for i in range(par.scandim):
        #    bbobj=bbox.getBoundary(i)
        #    bbobj.leftBoundary=par.scanrange[i][0]
        #    bbobj.rightBoundary=par.scanrange[i][1]
        #    bbox.setBoundary(i,bbobj)
        self.logdir=par.ddir+'/SGlogs'
        self.type='sparse_grid'
        self.nadapt=0
        

    def create_initial_grid(self):
        #initialize SGpp
        self.grid = Grid.createModLinearGrid(par.scandim)
#        self.grid = Grid.createLinearGrid(par.scandim)
#        self.grid = Grid.createLinearBoundaryGrid(par.scandim)
#        self.grid = Grid.createLinearTrapezoidBoundaryGrid(par.scandim)

        self.gridStorage = self.grid.getStorage()

        self.gridGen = self.grid.createGridGenerator()
        self.gridGen.regular(1)
        self.npoints=self.gridStorage.size()


    def reconstruct_grid(self):
        #use grid from previous computation as initial grid
        self.grid = Grid.unserialize(self.serialgrid)
        self.gridGen = self.grid.createGridGenerator()  
        self.gridStorage = self.grid.getStorage()
        self.npoints=self.gridStorage.size()

        self.alpha=DataVector(self.npoints)
        for i in xrange(self.npoints):
            self.alpha[i]=self.serialvec[i]
        #de-hierarchize to initialize resvec
        self.grid.createOperationHierarchisation().doDehierarchisation(self.alpha)
        for i in range(self.npoints):
            self.resvec.append(self.alpha[i])
        #hierarchize again for refinement etc
        self.grid.createOperationHierarchisation().doHierarchisation(self.alpha)
        
    def get_parlist(self,lower,upper):
        #create parlist
        parlist=[]
        for i in range(lower,upper):
            gp = self.gridStorage.get(i)
            parlist.append([])
            for j in range(par.scandim):
                #map from [0,1]^d to scanrange
                pval=gp.abs(j)*par.scanrange[j][2]+par.scanrange[j][0]
                parlist[-1].append(pval)

        return parlist

    def delete_old_res(self):
        self.resvec= DataVector(0)

    def add_res(self,reslist):
       # convert to coefficient vector
        for res in reslist:
            self.resvec.append(res)

        self.alpha=DataVector(self.npoints)
        for i in xrange(self.npoints):
            self.alpha[i]=self.resvec[i]

        #hierarchize
        self.grid.createOperationHierarchisation().doHierarchisation(self.alpha)

        self.serialgrid=self.grid.serialize()
        self.serialvec=[self.alpha[i] for i in range(len(self.alpha))]


    def refine(self):
        p = DataVector(par.scandim)
        refcrit=3

        for k in range(len(self.alpha)):
            ind = self.gridStorage.seq(self.gridStorage.get(k))
            gp=self.gridStorage.get(k)
            if refcrit==1:
                #weighting by value
                self.alpha[ind] = self.alpha[ind] * math.exp(-10*self.resvec[ind]**2)
                #weighting by volume
                self.alpha[ind] = self.alpha[ind] * (2**(par.scandim-gp.getLevelSum()))
                #points next to the boundary get extra factor
                for d in range(par.scandim):
                    (l,i) = gp.get(d)
                    if (i == 1 or i == 2**l-1):
                        self.alpha[ind] = self.alpha[ind]*2.
            elif refcrit==2:
                #weighting by volume
                self.alpha[ind] = self.alpha[ind] * 2**(par.scandim-gp.getLevelSum())
            elif refcrit==3:
                #weighting by volume
                self.alpha[ind] = self.alpha[ind] * (2**(par.scandim-gp.getLevelSum()))
                #points next to the boundary get extra factor
                for d in range(par.scandim):
                    (l,i) = gp.get(d)    
                    if (i == 1 or i == 2**l-1):
                        self.alpha[ind] = self.alpha[ind]*2.

        self.gridGen.refine(SurplusRefinementFunctor(self.alpha,self.nref,self.errlim))
        self.npoints=self.gridStorage.size()
        self.nadapt+=1
