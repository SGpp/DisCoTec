import os, sys, struct, scipy, numpy
from math import sqrt,ceil
from logfile import *
from scipy.optimize import leastsq

class checkpoints:
    def __init__(self,path_in,par_in):
        global par,path,distlog,corrlog
        path=path_in
        par=par_in
        self.g_len=int(par.pardict['n_spec'])*int(par.pardict['nx0'])*int(par.pardict['nky0'])*int(par.pardict['nz0'])
        self.g_len=self.g_len*int(par.pardict['nv0'])*int(par.pardict['nw0'])

        #initialize oldpars[0] with 0
        self.oldpars=numpy.zeros((1,par.scandim))
        self.metric_initialized=False
        distlog=logfile(par,'distances.log',False)
        corrlog=logfile(par,'correlation.log',False)
        self.noldpoints=0
        self.rectest=False
        print 'using numpy version ', numpy.version.version

    def get_closest_checkpoint(self,parpoint):
        if self.rectest:
            #this is only for testing purposes
            self.noldpoints+=1
            return self.noldpoints
        else:
            neighbour=0
            closestdist=1000.
            if len(self.oldpars)>1:
                diffvecs=numpy.array(parpoint)-self.oldpars[1:]
                dist=self.distance(self.metric,diffvecs)
                neighbour=numpy.argmin(dist)+1
            return neighbour

    def distance(self,metric,diffveclist):
        if int(numpy.version.version.split('.')[1])>5:
            distlist=numpy.einsum('ij,jk,ik->i',diffveclist,metric,diffveclist)
            return distlist
        else:
            prodlist=numpy.dot(diffveclist,metric)
            distlist=[]
            for ind in range(len(diffveclist)):
                dist=numpy.dot(prodlist[ind],diffveclist[ind])
                distlist.append(dist)        
            return numpy.array(distlist)

    def update(self,parlist):        
        self.oldpars=numpy.append(self.oldpars,numpy.array(parlist),axis=0)

    def residuals(self,metlist, y, x):
        metric=self.get_metric(metlist)
        err = y-self.distance(metric,x)
        return err

    def get_metric(self,coefflist):
        metric=numpy.zeros((par.scandim,par.scandim))   
        count=0
        for ind1 in range(par.scandim):
            #lower triangle
            for ind2 in range(ind1+1):
                metric[ind1,ind2]=coefflist[count]
                metric[ind2,ind1]=coefflist[count]
                count+=1
        return metric

    def init_metric(self,points_for_metric):
        if self.metric_initialized:
            return

        #initialize unit matrix
        metcoeffs=[]
        for ind1 in range(par.scandim):
            for ind2 in range(ind1+1):                        
                if ind1==ind2:
                    metcoeffs.append(1.)
                else:
                    metcoeffs.append(0.)
        metcoeffs=numpy.array(metcoeffs)
        self.metric=self.get_metric(metcoeffs)

        if par.scandim==1:
            #unit matrix is perfect
            self.metric_initialized=True
            return

        if points_for_metric > par.scandim:
            co=[]
            diffvecs=[]
            for ind1 in range(1,points_for_metric+1):
                ch1='%s/0k/checkpoint_%d' % (path,ind1)
                for ind2 in range(ind1,points_for_metric+1):
                    ch2='%s/0k/checkpoint_%d'%(path,ind2)
                    co.append(1.-self.correlate(ch1,ch2))
                    diffvecs.append(self.oldpars[ind1]-self.oldpars[ind2])
            lsq = leastsq(self.residuals, metcoeffs, args=(numpy.array(co),numpy.array(diffvecs)))
            metcoeffs=lsq[0]

            self.metric=self.get_metric(metcoeffs)
            self.metric_initialized=True
            print self.metric
        
        return

    def correlate(self,ch1,ch2):

        corrlist=[]
        ch1f=open(ch1,'rb')
        s1 = ch1f.read(6)            #remove precision header
        while ch1f.read(40)!='':     #checks for another eigenvector by looking for prec/time header
            filepos1=ch1f.tell()
            ch2f=open(ch2,'rb')      #remove precision header
            s2 = ch2f.read(6)
            while ch2f.read(40)!='': #checks for another eigenvector by looking for prec/time header
                ch1f.seek(filepos1)
                g1=scipy.fromfile(file=ch1f, dtype=complex, count=self.g_len)
                g2=scipy.fromfile(file=ch2f, dtype=complex, count=self.g_len)
                corr=abs(numpy.vdot(g1,g2))
                norm=sqrt(numpy.vdot(g1,g1).real*numpy.vdot(g2,g2).real)
                corrlist.append(corr/norm)
            ch2f.close()
        ch1f.close()
        corrlist.sort()
        avcorr=0.
        for ind in range(-int(par.pardict['n_ev']),0):
            avcorr+=corrlist[ind]
        avcorr/=int(par.pardict['n_ev'])
        return avcorr
#        return corrlist[-1]
