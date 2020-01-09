#!/usr/bin/env python
import struct
from os.path import getsize, join
import numpy as np

class fieldfile():
    #class constructor
    def __init__(self,file,pars):
        self.pars=pars
        self.file=file
        self.set_gridcounts()
        self.set_sizes()
        self.set_datatypes()
        self.te,self.tesize=self.TimeEntry()
        self.define_arrays()
        self.redirect(self.file)

    #call this routine to read from a new field file
    def redirect(self,file):
        self.file=file
        try: 
            close(self.f)
        except:
            pass
        self.f=open(file,'rb')
        self.tfld=[]
        self.get_timearray()
        self.reset_tinds()

    #set resolution
    def set_gridcounts(self):
        self.nx=int(self.pars['nx0'])
        self.ny=int(self.pars['nky0'])
        self.nz=int(self.pars['nz0'])

    #entry sizes
    def set_sizes(self):
        self.nfields=int(self.pars['n_fields'])
        self.intsize=4
        realsize=4+(self.pars['PRECISION']=='DOUBLE')*4
        complexsize=2*realsize
        self.entrysize=self.nx*self.ny*self.nz*complexsize
        #jumps in bytes in field/mom files
        self.leapfld=self.nfields*(self.entrysize+2*self.intsize)

    #real and complex datatypes according to endianness
    def set_datatypes(self):
        try: self.bigendian=pars['ENDIANNESS']=='BIG'
        except:
            self.bigendian=False
        if self.bigendian:
            self.nprt=(np.dtype(np.float64)).newbyteorder()
            self.npct=(np.dtype(np.complex128)).newbyteorder()
        else:
            self.nprt=np.dtype(np.float64)
            self.npct=np.dtype(np.complex128)

    def define_arrays(self):
        self.phi3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
        self.apar3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
            

    def get_timearray(self):
        #get time arrays for field file
        for i in range(int(getsize(self.file)/(self.leapfld+self.tesize))):
            self.tfld.append(float(self.te.unpack(self.f.read(self.tesize))[1]))
            self.f.seek(self.leapfld,1)

    def get_minmaxtime(self):
        if not self.tfld:
            self.get_timearray()
        return self.tfld[0],self.tfld[-1]

    #defines the struct for a time entry in field
    def TimeEntry(self):
            if self.bigendian:
                    timeentry=struct.Struct('>idi')
            else:
                    timeentry=struct.Struct('=idi')
            return timeentry,timeentry.size

    #calculate offset in field file for a given timestep and variable
    def offset(self,var):
            if var in [i for i in range(self.nfields)]:
                    return self.tesize+self.tind*(self.tesize+self.leapfld)+var*(self.entrysize+2*self.intsize)+self.intsize

    #returns field for given timestep
    def readvar(self,var):
        self.f.seek(self.offset(var))	
        var3d=np.fromfile(self.f,count=self.nx*self.ny*self.nz,dtype=self.npct).reshape(self.nz,self.ny,self.nx)
        return var3d

    #return time index for given time, if it's present
    #otherwise, and exception is raised
    def get_tind(self):
        return self.tfld.index(self.time)

    #set current timestep
    def set_time(self,time):
        self.time=time
        self.tind=self.get_tind()

    #reset time indices
    def reset_tinds(self):
        self.tind=0
        self.ftind=[-1]*self.nfields

    #return phi. this will only read from file if necessary
    def phi(self):
        if self.ftind[0]==self.tind:
            pass
        else:
            self.ftind[0]=self.tind
            self.phi3d=self.readvar(0)
        return self.phi3d

    #same for apar
    def apar(self):
        if self.nfields>1:
            if self.ftind[1]==self.tind:
                pass
            else:
                self.ftind[1]=self.tind
                self.apar3d=self.readvar(1)
            return self.apar3d
