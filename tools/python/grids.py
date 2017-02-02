import numpy as np
from ParIO import *

def get_grids(pars='get'):
    #if par['x_local']==False or par['y_local']==False:
    #    print "Need to set up get_grids for x,y global!"
    #    return
    #else:
        if pars=='get':
            par0=Parameters()
            parfile=raw_input('Enter parameter file name:')
            par0.Read_Pars(parfile)
            pars=par0.pardict
            self.pars=pars
            pars['suffix']=suffix

        par0=Parameters()
        parfile=raw_input('Eneter parameter file name:')
        par0.Read_Pars(parfile)
        par=par0.pardict
        print par
        par['nx0']=int(float(par['nx0']))
        par['nky0']=int(float(par['nky0']))
        par['nz0']=int(float(par['nz0']))
        par['kymin']=float(par['kymin'])
        kxmin=2.0*np.pi/float(par['lx'])
        kxmax=par['nx0']/2*kxmin
        kx=np.arange(par['nx0'],dtype='float')/(par['nx0']-1)*kxmin*(par['nx0']-1)-(kxmax-kxmin)
        kx=np.roll(kx,par['nx0']/2+1)
        kx=kx-kx[0]
        #print kx
        ky=np.arange(par['nky0'],dtype='float')*par['kymin']
        ky=ky-ky[0]
        #print ky
        kzmin=2.0*np.pi/par['nz0']
        zgrid=np.arange(par['nz0'],dtype='float')/(par['nz0']-1)*(2.0*np.pi-kzmin)-np.pi
        #print zgrid
        return kx,ky,zgrid

        par['nx0']=int(float(par['nx0']))
        par['nky0']=int(float(par['nky0']))
        par['nz0']=int(float(par['nz0']))
        par['kymin']=float(par['kymin'])
        kxmin=2.0*np.pi/float(par['lx'])
        kxmax=par['nx0']/2*kxmin
        kx=np.arange(par['nx0'],dtype='float')/(par['nx0']-1)*kxmin*(par['nx0']-1)-(kxmax-kxmin)
        kx=np.roll(kx,par['nx0']/2+1)
        kx=kx-kx[0]
        #print kx
        ky=np.arange(par['nky0'],dtype='float')*par['kymin']
        ky=ky-ky[0]
        #print ky
        kzmin=2.0*np.pi/par['nz0']
        zgrid=np.arange(par['nz0'],dtype='float')/(par['nz0']-1)*(2.0*np.pi-kzmin)-np.pi
        #print zgrid
        return kx,ky,zgrid



