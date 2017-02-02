import sys
from parameters_scan import *

class dense_grids:
    
    def __init__(self,par_in):
        global par
        par=par_in
        print "dense scan in %s" % par.scanpars         
        self.logdir=par.ddir+'/DGlogs'
        self.type='dense_grid'
        self.nadapt=0
        self.cvalgrid=[]
        self.serialvec=[]
        self.serialgrid=[]
        self.nintervals=0

    def create_initial_grid(self):
        #initial grid: corners of scan interval
        self.nintervals=1

        linset=[]
        for ind in range(par.scandim):
            linset.append(range(self.nintervals+1))

        self.serialgrid=self.create_parset(linset)
        self.cvalgrid=self.ind2range(self.serialgrid)
        self.npoints=len(self.cvalgrid)
        return 

    def reconstruct_grid(self):
        #initial grid: corners of scan interval
        self.nintervals=2**self.nadapt

        linset=[]
        for ind in range(par.scandim):
            linset.append(range(self.nintervals+1))

        self.cvalgrid=self.ind2range(self.serialgrid)
        self.npoints=len(self.cvalgrid)
        return

    def ind2range(self,indgrid):
        rangegrid=[]
        for cind in indgrid:
            cval=[]
            for dim in range(par.scandim):
                cval.append(par.scanrange[dim][0]+float(cind[dim])/self.nintervals*par.scanrange[dim][2])
            rangegrid.append(cval)

        return rangegrid

    def get_parlist(self,lower,upper):
        #create parlist
        parlist=[]
        for i in range(lower,upper):
            parlist.append(self.cvalgrid[i])

        return parlist

    def delete_old_res(self):
        self.serialvec=[]

    def add_res(self,reslist):
        self.serialvec.extend(reslist)

    def refine(self):
        #double the number of parameter intervals in each dimension
        self.nintervals=self.nintervals*2
        
        linset=[]
        for ind in range(par.scandim):
            linset.append(range(self.nintervals+1))

        newinds=self.create_parset(linset)

        #transform the old grid indices to the new range of points
        for ind1 in range(len(self.serialgrid)):
            for ind2 in range(len(self.serialgrid[ind1])):
                self.serialgrid[ind1][ind2]=self.serialgrid[ind1][ind2]*2

        #remove the gridpoints that have been computed before
        for element in self.serialgrid:
            newinds.remove(element)


        self.serialgrid.extend(newinds)  
        self.cvalgrid.extend(self.ind2range(newinds))
        self.npoints+=len(newinds)
        self.nadapt+=1
         
    def create_parset(self,linset):
        nparsets=[1]
        for ind in range(len(linset)):
            nparsets.append(nparsets[ind]*len(linset[ind]))

        parset=[]
        for pset in range(nparsets[-1]):
            parset.append([])
            rest=pset
            for ind in range(len(linset)):
                div=rest/nparsets[-ind-2]
                parset[pset].insert(0,linset[-ind-1][div])
                rest=rest-div*nparsets[-ind-2]

        return parset


    def createdenselist(self):
#        print par.stepwidth
        linset=[]
        for ind in range(par.scandim):
            linset.append([])
            i=0
            while (par.scanrange[ind][0]+i*par.stepwidth[ind]) <= par.scanrange[ind][1]:
                linset[ind].append(par.scanrange[ind][0]+i*par.stepwidth[ind])
                i=i+1
        parset=self.create_parset(linset)
        return(parset)
