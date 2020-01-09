import numpy as np
from os.path import join

class common_data:
    def __init__(self,p):
        global pars
        pars=p
        self.nx=int(pars['nx0'])
        self.ny=int(pars['nky0'])
        self.nz=int(pars['nz0'])
        self.kymin=float(pars['kymin'])
        self.kygrid=np.array([self.kymin*j for j in range(self.ny)],dtype=np.float64)
        try: 
            lx=float(pars['lx'])
        except:	
            lx=1
        self.zgrid=np.array([-np.pi+k*2*np.pi/self.nz for k in range(self.nz)],dtype=np.float64)
        self.x_local=self.isTrue('x_local')
        if self.x_local:
            self.kxmin=2*np.pi/lx
            self.kxgridp=[self.kxmin*i for i in range(int(self.nx/2+1))]
            self.kxgridn=[self.kxmin*i for i in range(int(-self.nx/2+1),0)]
            self.kxgrid=np.array(self.kxgridp+self.kxgridn,dtype=np.float64)
            #ordered kx array
            self.kxgrid_ord=np.array(self.kxgridn+self.kxgridp[0:-1],dtype=np.float64)
        else:
            self.xgrid=np.array([-lx/2.+lx/(self.nx-1)*i for i in range(0,self.nx)])
            try:
                self.minor_r=float(pars['minor_r'])
            except:
                pass
            self.rhostar=float(pars['rhostar'])
            self.x0=float(pars['x0'])
            self.x_agrid=self.xgrid*self.rhostar+self.x0

    def isTrue(self,var):
        try: 
            return pars[var] in ['T','t','.t.','.T.']
        except:
            return True

    def set_profiles(self,spec):
        try:
            dat=np.genfromtxt(join(self.rundir,'profiles_'+spec+self.fileext))
            xgrid=dat[:,0]
            self.temp_prof=dat[:,2]
            self.dens_prof=dat[:,3]
        except:
            self.temp_prof=1.
            self.dens_prof=1
