import numpy as np
from os.path import join

class geometry:
    def __init__(self,cm,p):
        global common,pars
        common=cm; pars=p

        self.geomtype=pars['magn_geometry']
        geom=self.getgeom()
        self.gxx=geom[0]
        self.gxy=geom[1]
        self.gxz=geom[2]
        self.gyy=geom[3]
        self.gyz=geom[4]
        self.gzz=geom[5]
        self.Bfield=geom[6]
        self.dBdx=geom[7]
        self.dBdy=geom[8]
        self.dBdz=geom[9]
        self.jacobian=geom[10]
        self.R=geom[11]
        self.Z=geom[12]
        self.drhodr=geom[13]
        self.drhodz=geom[14]
        self.Cy=geom[15]
        if not common.x_local:
            self.Cxy=geom[16]
            self.q=geom[17]
        
        
        

    #returns the geometry from a non-hdf5 file
    def getgeom(self):
        try: 
            if pars['x_local'] in ['.t.','.T.','t']: 
                local=1
            else:
                local=0
        except:
            local=1
        try:
            if pars['y_local'] in ['.t.','.T.','t']: 
                ylocal=1
            else:
                ylocal=0
        except:
            ylocal=1
        if not local:
            var=common.nx
        if not ylocal:
            var=common.ny
        if local and ylocal:
            geom=np.empty((16,common.nz),dtype=np.float64)
        else:
            geom=np.empty((18,common.nz,var),dtype=np.float64)

        g=open(join(common.rundir,self.geomtype.strip("'"))+common.fileext)
        i=0
        k=0
        if not (local and ylocal):
            dict={'gxx': 0, 'gxy': 1, 'gxz': 2, 'gyy': 3, 'gyz':4, 'gzz': 5,
                  'Bfield': 6, 'dBdx': 7, 'dBdy': 8, 'dBdz': 9, 'jacobian': 10,
                  'geo_R':11, 'geo_Z': 12, 'geo_c1': 13, 'geo_c2': 14, 'C_y': 15,
                  'C_xy': 16, 'q': 17}
            for line in g.readlines():
                if len(line.split())!=16:
                    try: 
                            num=dict[line.strip()]
                            i=0
                            k=0
                            ik=0
                    except: pass
                else:
                    if num in [15,16,17]:
                            geom[num,k,i:i+16]=line.split()[:]
                            i=(i+16)%var
                    else:
                            arr=line.split()[:]
                            lb=ik%var
                            ub=min(lb+16,var)
                            k=ik/var
                            dim=ub-lb
                            geom[num,k,lb:ub]=arr[0:dim]
                            if dim<16:
                                geom[num,k+1,0:16-dim]=arr[dim:]
                            ik+=16
        else:
            for line in g.readlines():
                if len(line.split())==16:
                    geom[:,k]=line.split()[:]
                    k+=1
        return geom
