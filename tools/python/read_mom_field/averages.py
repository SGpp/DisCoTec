import numpy as np

class averages:
    def __init__(self,cm,gm):
        global common,geom,nx,ny,nz
        common=cm; geom=gm
        nx=common.nx;ny=common.ny;nz=common.nz
        
        #init weights for time averages
        common.tweights=np.empty(common.tlen,dtype=np.float64)
        if common.tlen > 1:
            for i in range(1,common.tlen-1):
                common.tweights[i]=0.5*(common.times[i+1]-common.times[i-1])
            common.tweights[0]=common.times[1]-common.times[0]
            common.tweights[-1]=common.times[-1]-common.times[-2]
        else:
            common.tweights[0] = 1
        if not common.x_local:
            self.jaco3d=np.repeat(geom.jacobian[:,:],ny,axis=1).reshape(nz,ny,nx)

    def t(self,var):
        return np.average(var,weights=common.tweights,axis=0)

    def z(self,var):
        if common.x_local:
            return np.average(var,weights=geom.jacobian,axis=0)
        else:
            return np.average(var,weights=self.jaco3d,axis=0)
