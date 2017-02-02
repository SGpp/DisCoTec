#!/usr/bin/env python
#Extracts profile data from HDF5 files and writes GENE-readable
#profile files for Ions and Electrons, assuming equal densities.

try: 
    import h5py as h5
except:
    print 'This script requires the Python module h5py in order to open HDF5 files.'
    exit(1)

import numpy as np
from sys import argv,exit
#from matplotlib.pyplot import *


if len(argv)!=3:
    print 'Need two arguments: filename and desired coordinate.'
    exit(1)
else:
    filename=argv[1]
    coord=argv[2]
coords=('rho_v','rho_t','rho_p')
if coord not in coords:
    print 'Second argument has to be %s, %s or %s.' %coords
    exit(1)

#read relevant contents of h5 file
h5dat=h5.File(filename)
psigrid=h5dat['data']['grid']['PSI'][:]
Te=h5dat['data']['var1d']['Te'][:]
dTedpsi=h5dat['data']['var1d']['dTedpsi'][:]
Ti=h5dat['data']['var1d']['Ti'][:]
dTidpsi=h5dat['data']['var1d']['dTidpsi'][:]
dqdpsi=h5dat['data']['var1d']['dqdpsi'][:]
q=h5dat['data']['var1d']['q'][:]
ne=h5dat['data']['var1d']['ne'][:]
dnedpsi=h5dat['data']['var1d']['dnedpsi'][:]
ni=h5dat['data']['var1d']['ni'][:]
dnidpsi=h5dat['data']['var1d']['dnidpsi'][:]
dpsidrhotor=h5dat['data']['var1d']['dpsidrhotor'][:]
rho_tor=h5dat['data']['var1d']['rho_tor'][:]
volume=h5dat['data']['var1d']['Volume'][:]
dVdpsi=h5dat['data']['var1d']['dVdpsi'][:]
rho_psi=np.sqrt(psigrid/psigrid[-1])
rho_vol=np.sqrt(volume/volume[-1])

#select radial variable for output
coords={'rho_v': rho_vol, 'rho_p': rho_psi, 'rho_t': rho_tor}
radius=coords[coord]

#write ion profiles
proffile=open('profile_Ions','w')
proffile.write('#%11s %12s %12s %12s\n' %(coord,'rho_psi','Ti/keV','ne/1e19m^-3'))
proffile.write('#\n')
for i in range(len(psigrid)):
    proffile.write('%10.6e %10.6e %10.6e %10.6e\n' % (radius[i],rho_psi[i],Ti[i]/1e3,ne[i]/1e19))
proffile.close()

#write electron profiles
proffile=open('profile_Electrons','w')
proffile.write('#%11s %12s %12s %12s\n' %(coord,'rho_psi','Te/keV','ne/1e19m^-3'))
proffile.write('#\n')
for i in range(len(psigrid)):
    proffile.write('%10.6e %10.6e %10.6e %10.6e\n' % (radius[i],rho_psi[i],Te[i]/1e3,ne[i]/1e19))
proffile.close()

#fig1=Figure()
#plot(psigrid[0:-1],np.diff(psigrid)/np.diff(rho_tor))
#plot(psigrid,dpsidrhotor)
#psigrid=rho_psi
#plot(abs(psigrid),rho_tor,label='rho_tor')
#plot(abs(psigrid),rho_psi,label='rho_psi')
#plot(abs(psigrid),rho_vol,label='rho_vol')
#legend(loc=4)
#plot(psigrid,ni)
#show()

