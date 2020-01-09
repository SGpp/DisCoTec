from averages import *
from pylab import *

class heat_fluxes():
    def __init__(self,cm,f,m,a,p):
        global par, common, fld, mom, avg, nx,ny,nz
        common=cm;fld=f;mom=m;avg=a;par=p
        nx=common.nx;ny=common.ny;nz=common.nz
        tlen=common.tlen
        npct=fld.npct
        nprt=fld.nprt

	hf3d=np.empty((nz,ny,nx),dtype=npct)
	#self.hf4d=np.empty((tlen,nz,ny,nx),dtype=npct)
	#hf_kx_z=np.empty((tlen,nz,nx-1),dtype=npct)
	#hf_ky_z=np.empty((tlen,nz,ny),dtype=npct)
	#hf_kx=np.empty((tlen,nx-1),dtype=npct)
	#hf_ky=np.empty((tlen,ny),dtype=npct)
	#hf_kxky=np.empty((tlen,ny,nx-1),dtype=npct)
	#hf_zprof=np.empty((tlen,nz),dtype=npct)
        self.hf_tyx=np.empty((tlen,ny,nx),dtype=npct)
	#hf_avg=np.empty((tlen),dtype=npct)
        self.hf_xprof=np.empty((tlen,nx),dtype=nprt)
        self.hf_tind=-1

    def hf_spec(which='y',zarg=-1):
        var=self.hf3d
        if which=='x':
            xspec=2*np.sum(var,axis=1)
            xspec[:,0]=var[:,0,0]+2*np.sum(var[:,:,0],axis=1)
            spec=orderzx(xspec)
        elif which=='y':	
            yspec=2*np.sum(var,axis=2)
            yspec[:,0]=var[:,0,0]+2*np.sum(var[:,0,0:nx/2-1],axis=1)
            spec=yspec
        elif which=='xy':
            spec=2*orderzyx(var)
        elif which=='avg':
            spec=2*np.sum(np.sum(hf3d,axis=1),axis=1)
        if zarg==-1:
            var_out=np.average(spec,weights=geom.jacobian,axis=0)
        elif zarg>=0:
            var_out=spec[zarg]
        elif zarg==-2:
            var_out=spec
        return var_out

    def calc_hf(self):
        dp3d=common.dens_prof.reshape(1,1,nx)
        tp3d=common.temp_prof.reshape(1,1,nx)
        ky3d=common.kygrid.reshape(1,ny,1)
        if common.tind!=self.hf_tind:
            self.hf3d=(-1j*ky3d*fld.phi()*np.conj(
                    1.5*mom.dens()*dp3d+
                    0.5*mom.tpar()*dp3d+mom.tperp()*tp3d)) 
        #remember that hf has already been calculated for this time
        self.hf_tind=common.tind
            
    def store_4d(self):
        self.hf4d[common.tind]=self.hf3d

    def store_zavg(self):
        self.hf_tyx[common.tind]=avg.z(self.hf3d)

    def high_ky_fraction(self):
        self.calc_hf()
        self.store_zavg()
        
    def plot(self):
        hf_yx=avg.t(self.hf_tyx)
        contourf(hf_yx)
        show()

