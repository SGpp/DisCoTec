#!/usr/bin/env python
import struct
from os.path import getsize, join
import numpy as np

class momfile():
	def __init__(self,file,pars):
            self.pars=pars
            self.file=file
            self.set_gridcounts()
            self.set_sizes()
            self.set_datatypes()
            self.define_arrays()
            self.redirect(self.file)

        def redirect(self,file):
	    self.file=file	 
            try: 
                close(self.m)
            except:
                pass
            self.m=open(file,'rb')
            self.tmom=[]
            self.te,self.tesize=self.TimeEntry()
            self.get_timearray()
            self.reset_tinds()
            
        
        def set_gridcounts(self):
            self.nx=int(self.pars['nx0'])
            self.ny=int(self.pars['nky0'])
            self.nz=int(self.pars['nz0'])

        def set_sizes(self):
            self.nmoms=int(self.pars['n_moms'])
            self.intsize=4
            realsize=4+(self.pars['PRECISION']=='DOUBLE')*4
            complexsize=2*realsize
            self.entrysize=self.nx*self.ny*self.nz*complexsize
            #jumps in bytes in field/mom files
            self.leapmom=self.nmoms*(self.entrysize+2*self.intsize)

           
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
		self.dens3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
		self.tpar3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
		self.tperp3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
		self.qpar3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
		self.qperp3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
		self.upar3d=np.empty((self.nz,self.ny,self.nx),dtype=self.npct)
            
           
        def get_timearray(self):
#get time arrays for field file
            for i in range(getsize(self.file)/(self.leapmom+self.tesize)):
                self.tmom.append(float(self.te.unpack(self.m.read(self.tesize))[1]))
                self.m.seek(self.leapmom,1)

        def get_minmaxtime(self):
            if not self.tmom:
                self.get_timearray()
            return self.tmom[0],self.tmom[-1]

#defines the struct for a time entry in field
	def TimeEntry(self):
		if self.bigendian:
			timeentry=struct.Struct('>idi')
		else:
			timeentry=struct.Struct('=idi')
		return timeentry,timeentry.size


	def offset(self,var):
		if var in [i for i in range(self.nmoms)]:
			leap=self.leapmom
			return self.tesize+self.tind*(self.tesize+leap)+var*(self.entrysize+2*self.intsize)+self.intsize

#returns field for given timestep
	def readvar(self,var):
            self.m.seek(self.offset(var))	
            var3d=np.fromfile(self.m,count=self.nx*self.ny*self.nz,dtype=self.npct).reshape(self.nz,self.ny,self.nx)
            return var3d

        def reset_tinds(self):
            self.tind=0
            self.mtind=[-1]*self.nmoms

        def get_tind(self):
            if not self.tmom:
		    self.get_timearray()
            return self.tmom.index(self.time)

        def set_time(self,time):
            self.time=time
            self.tind=self.get_tind()

	def dens(self):
		if self.mtind[0]==self.tind:
			pass
		else:
			self.mtind[0]=self.tind
			self.dens3d=self.readvar(0)
		return self.dens3d


	def tpar(self):
		if self.mtind[1]==self.tind:
			pass
		else:
			self.mtind[1]=self.tind
			self.tpar3d=self.readvar(1)
		return self.tpar3d

	def tperp(self):
		if self.mtind[2]==self.tind:
			pass
		else:
			self.mtind[2]=self.tind
			self.tperp3d=self.readvar(2)
		return self.tperp3d

	def qpar(self):
		if self.mtind[3]==self.tind:
			pass
		else:
			self.mtind[3]=self.tind
			self.qpar3d=self.readvar(3)
		return self.qpar3d

	def qperp(self):
		if self.mtind[4]==self.tind:
			pass
		else:
			self.mtind[4]=self.tind
			self.qperp3d=self.readvar(4)
		return self.qperp3d

	def upar(self):
		if self.mtind[5]==self.tind:
			pass
		else:
			self.mtind[5]=self.tind
			self.upar3d=self.readvar(5)
		return self.upar3d

