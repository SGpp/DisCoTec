#!/usr/bin/env python
import struct
import numpy as np
from os.path import getsize, join
from ParIO import *
from grids import *
import matplotlib.pyplot as plt
import os

class energyfile():
    def __init__(self,pars='get',suffix='.dat',verbose=False):
        if pars=='get':
            par0=Parameters()
            par0.Read_Pars('./parameters'+suffix)
            pars=par0.pardict
            self.pars=pars
            pars['suffix']=suffix
      
        self.verbose=verbose
        kx,ky,zgrid=get_grids(pars=self.pars)
        self.kx=kx
        self.ky=ky
        self.zgrid=zgrid
        self.suffix=suffix
        if os.path.isfile(pars['diagdir'][1:-1]+'/energy'+suffix):
            self.file=pars['diagdir'][1:-1]+'/energy'+suffix
            self.file3d=pars['diagdir'][1:-1]+'/energy3d'+suffix
            self.diagdir=self.pars['diagdir'][1:-1]
        else:
            diagdir=raw_input("Enter directory:\n")
            self.diagdir=diagdir
            self.file=diagdir+'/energy'+suffix
            self.file3d=diagdir+'/energy3d'+suffix
        self.set_gridcounts()
        self.set_sizes()
        self.set_datatypes()
        self.te,self.tesize=self.TimeEntry()
        self.define_arrays()
        self.redirect(self.file3d)
        temp=np.genfromtxt(self.file)
        self.time=temp[:,0]
        self.energy1d=temp[:,1]
        self.dedt1d=temp[:,2]
        self.drive1d=temp[:,3]
        self.sources1d=temp[:,4]
        self.coll1d=temp[:,5]
        self.hypv1d=temp[:,6]
        self.hypz1d=temp[:,7]
        self.hypxy1d=temp[:,8]
        self.nl1d=temp[:,9]
        self.zv1d=temp[:,10]
        self.curv1d=temp[:,11]
        self.ctest1d=temp[:,12]
        self.dedttest1d=temp[:,13]

    #call this routine to read from a new energy file
    def redirect(self,file):
        self.file3d=file
        try: 
            close(self.f)
        except:
            pass
        self.f=open(file,'rb')
        self.tfld=np.array([])
        self.get_timearray()
        self.reset_tinds()


    def set_sizes(self):
        self.nentries=6
        self.intsize=4
        realsize=4+(self.pars['PRECISION']=='DOUBLE')*4
        complexsize=2*realsize
        self.entrysize=self.nx*self.ny*self.nz*realsize
        #jumps in bytes in field/mom files
        self.leapfld=self.nentries*(self.entrysize+2*self.intsize)

    #set resolution
    def set_gridcounts(self):
        self.nx=int(self.pars['nx0'])
        self.ny=int(self.pars['nky0'])
        self.nz=int(self.pars['nz0'])


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
        self.energy3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)
        #All nonconservative terms
        self.dedtnc3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)
        #Drive term
        self.drive3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)
        #Collisions
        self.coll3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)
        #All artificial dissipation (hyp_x,y,z,v)
        self.diss3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)
        #Nonlinearity
        self.nl3d=np.empty((self.nz,self.ny,self.nx),dtype=self.nprt)

    def get_timearray(self):
        #get time arrays for field file
        #print getsize(self.file3d)/(self.leapfld+self.tesize)
        for i in range(getsize(self.file3d)/(self.leapfld+self.tesize)):
            if self.verbose:
                print "i",i
            self.tfld=np.append(self.tfld,(float(self.te.unpack(self.f.read(self.tesize))[1])))
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
            if var in [i for i in range(self.nentries)]:
                    return self.tesize+self.tind*(self.tesize+self.leapfld)+var*(self.entrysize+2*self.intsize)+self.intsize

    #returns field for given timestep
    def readvar(self,var):
        self.f.seek(self.offset(var))	
        var3d=np.fromfile(self.f,count=self.nx*self.ny*self.nz,dtype=self.nprt).reshape(self.nz,self.ny,self.nx)
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
        self.ftind=[-1]*self.nentries

    #return energy3d. this will only read from file if necessary
    def energy(self):
        if self.ftind[0]==self.tind:
            pass
        else:
            self.ftind[0]=self.tind
            self.energy3d=self.readvar(0)
        return self.energy3d

    #return dedtnc3d. this will only read from file if necessary
    def dedtnc(self):
        if self.ftind[1]==self.tind:
            pass
        else:
            self.ftind[1]=self.tind
            self.dedtnc3d=self.readvar(1)
        return self.dedtnc3d

    #return drive3d. this will only read from file if necessary
    def drive(self):
        if self.ftind[2]==self.tind:
            pass
        else:
            self.ftind[2]=self.tind
            self.drive3d=self.readvar(2)
        return self.drive3d

    #return coll3d. this will only read from file if necessary
    def coll(self):
        if self.ftind[3]==self.tind:
            pass
        else:
            self.ftind[3]=self.tind
            self.coll3d=self.readvar(3)
        return self.coll3d

    #return diss3d. this will only read from file if necessary
    def diss(self):
        if self.ftind[4]==self.tind:
            pass
        else:
            self.ftind[4]=self.tind
            self.diss3d=self.readvar(4)
        return self.diss3d

    #return nl3d. this will only read from file if necessary
    def nl(self):
        if self.ftind[5]==self.tind:
            pass
        else:
            self.ftind[5]=self.tind
            self.nl3d=self.readvar(5)
        return self.nl3d

    def plot1d(self,which_term=-1):
        if which_term==-1:
            plt.plot(self.time,self.drive1d,label='drive')
            plt.plot(self.time,self.coll1d,label='coll')
            plt.plot(self.time,self.hypv1d,label='hypv')
            plt.plot(self.time,self.hypz1d,label='hypz')
            plt.plot(self.time,self.hypxy1d,label='hypxy')
            plt.xlabel('time')
            plt.legend()
            plt.show()
        else:
            temp=np.genfromtxt(self.file)
            plt.plot(temp[:,0],temp[:,which_term],label=str(which_term))
            plt.xlabel('time')
            plt.legend()
            plt.show()

    def test1d3d(self):
        nt3=len(self.tfld)
        entest=np.empty(nt3,dtype='float')
        detest=np.empty(nt3,dtype='float')
        drtest=np.empty(nt3,dtype='float')
        cotest=np.empty(nt3,dtype='float')
        dstest=np.empty(nt3,dtype='float')
        nltest=np.empty(nt3,dtype='float')
        print "!!!Warning: need jacobian!!!"
        for i in range(nt3):
            print i
            self.tind=i
            self.energy3d=self.energy()
            self.dedtnc3d=self.dedtnc()
            self.drive3d=self.drive()
            self.coll3d=self.coll()
            self.diss3d=self.diss()
            self.nl3d=self.nl()
            entest[i]=(2.0*np.sum(self.energy3d[:,1:,:])+np.sum(self.energy3d[:,0,:]))/float(self.pars['nz0'])
            detest[i]=(2.0*np.sum(self.dedtnc3d[:,1:,:])+np.sum(self.dedtnc3d[:,0,:]))/float(self.pars['nz0'])
            drtest[i]=(2.0*np.sum(self.drive3d[:,1:,:])+np.sum(self.drive3d[:,0,:]))/float(self.pars['nz0'])
            cotest[i]=(2.0*np.sum(self.coll3d[:,1:,:])+np.sum(self.coll3d[:,0,:]))/float(self.pars['nz0'])
            dstest[i]=(2.0*np.sum(self.diss3d[:,1:,:])+np.sum(self.diss3d[:,0,:]))/float(self.pars['nz0'])
            nltest[i]=(2.0*np.sum(self.nl3d[:,1:,:])+np.sum(self.nl3d[:,0,:]))/float(self.pars['nz0'])
        plt.plot(self.time,self.energy1d,label='energy1d')
        plt.plot(self.tfld,entest,label='energy3d')
        plt.xlabel('time')
        plt.legend()
        plt.show()
        plt.plot(self.time,self.coll1d,label='coll1d')
        plt.plot(self.tfld,cotest,label='coll3d')
        plt.xlabel('time')
        plt.legend()
        plt.show()
        plt.plot(self.time,self.drive1d,label='drive1d')
        plt.plot(self.tfld,drtest,label='drive3d')
        plt.xlabel('time')
        plt.legend()
        plt.show()
        plt.plot(self.time,self.coll1d+self.hypv1d+self.hypz1d+self.hypxy1d,label='diss1d')
        plt.plot(self.tfld,dstest,label='diss3d')
        plt.xlabel('time')
        plt.legend()
        plt.show()

    def dissipation_kg1_vs_kl1(self,start_time=-1.0,end_time=-1.0):

        ctot=0.0
        cgt=0.0
        cls=0.0
        dtot=0.0
        dgt=0.0
        dls=0.0
        if start_time==-1.0:
            start_time=self.tfld[0]
        if end_time==-1.0:
            end_time=self.tfld[-1]
        if start_time >= end_time:
            stop

        print "!!!!Warning using s-alpha (alpha=0) geometry . . . need to generalize!!!"
        print "Using shat = ",float(self.pars['shat'])
        print "Using epsilon_t = ",float(self.pars['trpeps'])
        gxx=np.zeros(self.nz,dtype='float')
        gxy=np.zeros(self.nz,dtype='float')
        gyy=np.zeros(self.nz,dtype='float')
        jacobian=np.zeros(self.nz,dtype='float')
        gxx[:]=1.0
        for i in range(self.nz):
            gyy[i]=1.0+(self.zgrid[i]*float(self.pars['shat']))**2
            gxy[i]=float(self.pars['shat'])*self.zgrid[i]
            jacobian[i]=1.0+float(self.pars['trpeps'])*np.cos(self.zgrid[i])
        avg_jac=np.sum(jacobian)/float(self.nz)
        print "avg_jac",avg_jac
        
        plt.plot(self.zgrid,gxx,label='gxx')
        plt.plot(self.zgrid,gyy,label='gyy')
        plt.plot(self.zgrid,gxy,label='gxy')
        plt.plot(self.zgrid,jacobian,label='jacobian')
        plt.legend()
        plt.show()

        istart=np.argmin(abs(self.tfld-start_time))
        iend=np.argmin(abs(self.tfld-end_time))
        ntime=iend-istart+1
        print "istart", istart
        print "iend", iend
        print "ntime", ntime

        cavg=np.zeros((self.nz,self.ny,self.nx),dtype='float')
        davg=np.zeros((self.nz,self.ny,self.nx),dtype='float')
        for i in range(istart,iend+1):
            print i, " of ", iend
            self.tind=i
            cavg+=self.coll()
            davg+=self.diss()

        cavg=cavg/float(ntime)
        davg=davg/float(ntime)

        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    kperp=np.sqrt(gxx[k]*self.kx[i]**2+gyy[k]*self.ky[j]**2+2.0*gxy[k]*self.kx[i]*self.ky[j])
                    if j==0:
                        factor=1.0
                    else:
                        factor=2.0
                    c0=factor*(jacobian[k]*cavg[k,j,i])/(avg_jac*float(self.pars['nz0']))
                    d0=factor*(jacobian[k]*davg[k,j,i])/(avg_jac*float(self.pars['nz0']))
                    ctot+=c0
                    dtot+=d0
                    if kperp > 1:
                    #    print i,j,self.kx[i],self.ky[j],kperp,self.zgrid[k], "gt"
                        cgt+=c0
                        dgt+=d0
                    else:
                    #    print i,j,self.kx[i],self.ky[j],kperp,self.zgrid[k], "lt"
                        cls+=c0
                        dls+=d0


        print "total dissipation:", ctot+dtot
        print "dissipation at kperp<1:", cls+dls
        print "dissipation at kperp>1:", cgt+dgt

        print "total collisional diss:", ctot
        print "collisional dissipation at kperp<1:", cls
        print "collisional dissipation at kperp>1:", cgt

        print "total hyp diss:", dtot
        print "hyp dissipation at kperp<1:", dls
        print "hyp dissipation at kperp>1:", dgt


    def energy_slice(self,ikx,iky,start_time=-1.0,end_time=-1.0):

        if start_time==-1.0:
            start_time=self.tfld[0]
        if end_time==-1.0:
            end_time=self.tfld[-1]
        if start_time >= end_time:
            stop

        print "!!!!Warning using s-alpha (alpha=0) geometry . . . need to generalize!!!"
        print "Using shat = ",float(self.pars['shat'])
        print "Using epsilon_t = ",float(self.pars['trpeps'])
        gxx=np.zeros(self.nz,dtype='float')
        gxy=np.zeros(self.nz,dtype='float')
        gyy=np.zeros(self.nz,dtype='float')
        jacobian=np.zeros(self.nz,dtype='float')
        gxx[:]=1.0
        for i in range(self.nz):
            gyy[i]=1.0+(self.zgrid[i]*float(self.pars['shat']))**2
            gxy[i]=float(self.pars['shat'])*self.zgrid[i]
            jacobian[i]=1.0+float(self.pars['trpeps'])*np.cos(self.zgrid[i])
        avg_jac=np.sum(jacobian)/float(self.nz)
        print "avg_jac",avg_jac

        istart=np.argmin(abs(self.tfld-start_time))
        iend=np.argmin(abs(self.tfld-end_time))
        ntime=iend-istart+1
        print "istart", istart
        print "iend", iend
        print "ntime", ntime

        #cavg=np.zeros((self.nz,self.ny,self.nx),dtype='float')
        #davg=np.zeros((self.nz,self.ny,self.nx),dtype='float')
        eslice=np.zeros((ntime,7),dtype='float')
        for i in range(istart,iend+1):
            print i, " of ", iend
            self.tind=i
            temp=self.energy()
            eslice[i-istart,1]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))
            temp=self.dedtnc()
            eslice[i-istart,2]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))
            temp=self.drive()
            eslice[i-istart,3]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))
            temp=self.coll()
            eslice[i-istart,4]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))
            temp=self.diss()
            eslice[i-istart,5]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))
            temp=self.nl()
            eslice[i-istart,6]=np.sum(temp[:,iky,ikx]*jacobian[:])/(avg_jac*float(self.pars['nz0']))

        eslice[:,0]=self.tfld[istart:iend+1]
        eslice_average=np.sum(eslice,axis=0)/float(ntime)
        file_name=self.diagdir+'/eslice_ikx'+str(ikx)+'iky'+str(iky)+self.suffix
        f=open(file_name,'w')
        header_string = '#ikx='+str(ikx)+'\n'
        header_string += '#iky='+str(iky)+'\n'
        header_string += '#1:Time 2:Energy 3:dedtnc 4:drive 5:coll 6:hyp 7:NL \n'
        header_string += '#Averages (including time): \n'
        header_string += '#'+str(eslice_average)+'\n'
        f.write(header_string)
        np.savetxt(f,eslice)
        f.close()



