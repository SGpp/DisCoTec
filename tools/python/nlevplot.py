#!/usr/bin/env python
import sys
import numpy as np
from scipy.optimize import fsolve
import pylab as pl
from matplotlib.pyplot import *
from matplotlib.font_manager import fontManager, FontProperties
from ParIO import Parameters
import os
import re


#input arguments
evplot=1
dtevplot=0
collevplot=0
vlasevplot=0
infiles=[]
fextout=[]
fext=[]
nfiles=0
count=0
annotate=0
verbose=0
my_xylim = 0
tmin=0
tmax=999999999
itmin=1
itmax=0
auto_tmin=1
input_error = 0
iarg=0
legpos=-1111
legyoff=0
plotmodel = []
ymax = -1
ymin = -1
xmin = 1111
xmax = 1111
plotRKC1=0
plotRKC2=0
plotRKC3=0
plotRKC4=0
plotRK2=0
plotRK3=0
plotRK4=0
plotRK4M=0
dt_max = -1

def showhelp():
   print "usage:"
   print "nlevplot.py <eigenvalues file>:"
   print "options: "
   print "  -tmin <number>  : set tmin (time window)"
   print "  -tmax <number>  : set tmax (time window)"
   print "  -h              : show help"
   print "  -leg            : position of legend (0-10)"
   print "  -v              : verbose"
   print "  -collev         : plot coll ev's"
   print "  -vlasev         : plot vlasov ev's"
   print "  -dt             : show timetrace of dt_nlev"
   print "  -ev             : show timetrace of nlev eigenvalues"
   print "  -pm  <number>   : number of ev model to plot: 1,2 (0=none)"
   print "  -RK4            : show stability contour of the RK4 scheme in ev plot"
   print "                    (also works with -RK2 -RK4M -RKC1 -RKC2 -RKC3 and/or -RKC4)"
   print "  ->see source for more details. many figuer options are hardcoded"
   exit(0)

while iarg < (len(sys.argv)-1):
   iarg+=1 #skip first arg==skript name
   arg=sys.argv[iarg]
   valid=re.match(r'.*eigenvalues(_nl)?(.*)',arg)
   if valid:
      infiles.append(str(arg))
      ext=str(valid.group(2))
      fext.append(ext)
      if re.match(r'[.]dat$',ext):
         ext='_dat'
      fextout.append(ext)
      nfiles+=1
   if str(arg)=='-np':
      evplot=0
   if str(arg)=='-leg':
      iarg+=1
      legpos=int(float(sys.argv[iarg]))
      legyoff = float(sys.argv[iarg])-legpos
   if str(arg)=='-pm':
      iarg+=1
      plotmodel.append(int(sys.argv[iarg]))
   if str(arg)=='-h' or str(arg)=='--h' or str(arg)=='-help' or str(arg)=='--help':
      showhelp()
   if str(arg)=='-dt':
      dtevplot=1
      evplot=0
      dt2ev=1
   if str(arg)=='-ev':
      dtevplot=1
      evplot=0
      dt2ev=-1
   if str(arg)=='-collev':
      evplot=0
      collevplot=1
   if str(arg)=='-vlasev':
      evplot=0
      vlasevplot=1
   if str(arg)=='-ymax':
      iarg+=1
      ymax=float(sys.argv[iarg])
   if str(arg)=='-annotate':
      iarg+=1
      annotate=int(sys.argv[iarg])
   if str(arg)=='-ymin':
      iarg+=1
      ymin=float(sys.argv[iarg])
   if str(arg)=='-xmin':
      iarg+=1
      xmin=float(sys.argv[iarg])
   if str(arg)=='-xmax':
      iarg+=1
      xmax=float(sys.argv[iarg])
   if str(arg)=='-v':
      verbose=1
   if str(arg)=='-tmin':
      iarg+=1  #read following arg (should be number)
      tmin = float(sys.argv[iarg])
      auto_tmin=0
      if (tmin < 0):
         print "tmin ill defined"
         showhelp()
   if str(arg)=='-tmax':
      iarg+=1  #read following arg (should be number)
      tmax = float(sys.argv[iarg])
      if (tmax < 0):
         print "tmax ill defined"
         showhelp()
   if str(arg)=='-RKC1':
      plotRKC1=1
   if str(arg)=='-RKC2':
      plotRKC2=1
   if str(arg)=='-RKC3':
      plotRKC3=1
   if str(arg)=='-RKC4':
      plotRKC4=1
   if str(arg)=='-RK2':
      plotRK2=1
   if str(arg)=='-RK3':
      plotRK3=1
   if str(arg)=='-RK4':
      plotRK4=1
   if str(arg)=='-RK4M':
      plotRK4M=1
      
if (plotmodel==[]): plotmodel.append(2)
if (plotmodel[0] < 1): plotmodel=[]
print "plot models: ", plotmodel

if nfiles==0:
   print "no iniput files"
   showhelp()

for filename in infiles:
   if (not os.path.isfile(filename)):
      print "%s is not a file: exit\n" %(filename)
      showhelp()


def file_len(fname):
   with open(fname) as f:
      for i, l in enumerate(f):
         pass
   return i + 1

def parse_nlevfile(fname):
   with open(fname) as f:
      nblocks = 0
      nnlev = 0
      nmodel=np.zeros((2))
      max_ev_ind = 0
      ev_ind = 0
      i = 0
      for i, l in enumerate(f):
         isnlev=re.search(r'#.*nlev',l)
         if isnlev:
            nnlev+=1
         ismodel1=re.search(r'#.*model vemax',l)
         ismodel2=re.search(r'#.*model lin\+ve\s*$',l)
         if ismodel1: nmodel[0]+=1
         if ismodel2: nmodel[1]+=1
         ishead=re.search(r'#\S*linear',l)
         if ishead:
            max_ev_ind = max(max_ev_ind,ev_ind)
            nblocks+=1
            ev_ind = 0
         iseigenvalue=re.match(r'(?![#])\s*(\S+)\s*(\S+)\s*$',l)
         if iseigenvalue:
            ev_ind+=1
   return [i + 1,nblocks,nnlev,nmodel,max_ev_ind]


#taylor series expansion coefficients for various schemes
#first order RKC
RKC1_taylor=[1.0,0.0,0.0,0.0]   
RKC2_taylor=[1,0.128067367779301943819803000224]
RKC3_taylor=[1,0.1519274412957031,0.005791682236923549]
RKC4_taylor=[1,0.1602892682704187,0.00825224619254291,0.0001329321183124034]
#second order
RK2_taylor=[1,1./2]
RKC23_taylor=[1,1./2,0.06272956665954434]
RKC24_taylor=[1,1./2,0.0802977344285056,0.004034842618250562]
RKC25_taylor=[1,1./2,0.0878276180032986,0.006304771569814507,0.0001584977157570301]
#third order
turbo_taylor=[1,1./2,1./6,1./32]
RK3ssp_taylor=[1,1./2,1./6, 1.938478371506826E-002]
#4th order
RK4_taylor=[1,1./2,1./6,1./24]
RK4M_taylor=[1,1./2,1./6,1./24,0.005562334,0.0009340186]
#even higher order
RK5_taylor=[1,1./2,1./6,1./24,1./120]
RK6_taylor=[1,1./2,1./6,1./24,1./120,1./720]
def RKpoly(x,y,kin):
   #RK stability polynomial for up to six stages
   ## stability boundary: poly==0
   nstages=len(kin)
   k=np.zeros(6)
   for i in range(nstages):
      k[i]=kin[i]
   k1=k[0]
   k2=k[1]
   k3=k[2]
   k4=k[3]
   k5=k[4]
   k6=k[5]
   #poly = 1+ 2*x +x**2+y**2
   #poly = (x-1)**2+(y)**2
   #poly= 1. + 2*k1*x + k1**2*x**2 + 2*k2*x**2 + 2*k1*k2*x**3 + k2**2*x**4 + \
   #k1**2*y**2 - 2*k2*y**2 + 2*k1*k2*x*y**2 + 2*k2**2*x**2*y**2 + k2**2*y**4

   poly = \
   1 + 2*k1*x + k1**2*x**2 + 2*k2*x**2 + 2*k1*k2*x**3 + 2*k3*x**3 + \
   k2**2*x**4 + 2*k1*k3*x**4 + 2*k4*x**4 + 2*k2*k3*x**5 + \
   2*k1*k4*x**5 + 2*k5*x**5 + k3**2*x**6 + 2*k2*k4*x**6 + 2*k1*k5*x**6 + \
   2*k6*x**6 + 2*k3*k4*x**7 + 2*k2*k5*x**7 + \
   2*k1*k6*x**7 + k4**2*x**8 + 2*k3*k5*x**8 + 2*k2*k6*x**8 + 2*k4*k5*x**9 + \
   2*k3*k6*x**9 + k5**2*x**10 + 2*k4*k6*x**10 + \
   2*k5*k6*x**11 + k6**2*x**12 + k1**2*y**2 - 2*k2*y**2 + 2*k1*k2*x*y**2 - \
   6*k3*x*y**2 + 2*k2**2*x**2*y**2 - \
   12*k4*x**2*y**2 + 4*k2*k3*x**3*y**2 - 4*k1*k4*x**3*y**2 - \
   20*k5*x**3*y**2 + 3*k3**2*x**4*y**2 + 2*k2*k4*x**4*y**2 - \
   10*k1*k5*x**4*y**2 - 30*k6*x**4*y**2 + 6*k3*k4*x**5*y**2 - \
   2*k2*k5*x**5*y**2 - 18*k1*k6*x**5*y**2 + 4*k4**2*x**6*y**2 + \
   4*k3*k5*x**6*y**2 - 8*k2*k6*x**6*y**2 + 8*k4*k5*x**7*y**2 + \
   5*k5**2*x**8*y**2 + 6*k4*k6*x**8*y**2 + 10*k5*k6*x**9*y**2 + \
   6*k6**2*x**10*y**2 + k2**2*y**4 - 2*k1*k3*y**4 + 2*k4*y**4 + \
   2*k2*k3*x*y**4 - 6*k1*k4*x*y**4 + 10*k5*x*y**4 + \
   3*k3**2*x**2*y**4 - 2*k2*k4*x**2*y**4 - 10*k1*k5*x**2*y**4 + \
   30*k6*x**2*y**4 + 6*k3*k4*x**3*y**4 - 10*k2*k5*x**3*y**4 - \
   10*k1*k6*x**3*y**4 + 6*k4**2*x**4*y**4 - 20*k2*k6*x**4*y**4 + \
   12*k4*k5*x**5*y**4 - 12*k3*k6*x**5*y**4 + \
   10*k5**2*x**6*y**4 + 4*k4*k6*x**6*y**4 + 20*k5*k6*x**7*y**4 + \
   15*k6**2*x**8*y**4 + k3**2*y**6 - 2*k2*k4*y**6 + \
   2*k1*k5*y**6 - 2*k6*y**6 + 2*k3*k4*x*y**6 - 6*k2*k5*x*y**6 + \
   10*k1*k6*x*y**6 + 4*k4**2*x**2*y**6 - 4*k3*k5*x**2*y**6 - \
   8*k2*k6*x**2*y**6 + 8*k4*k5*x**3*y**6 - 16*k3*k6*x**3*y**6 + \
   10*k5**2*x**4*y**6 - 4*k4*k6*x**4*y**6 + 20*k5*k6*x**5*y**6 + \
   20*k6**2*x**6*y**6 + k4**2*y**8 - 2*k3*k5*y**8 + 2*k2*k6*y**8 + \
   2*k4*k5*x*y**8 - 6*k3*k6*x*y**8 + 5*k5**2*x**2*y**8 - \
   6*k4*k6*x**2*y**8 + 10*k5*k6*x**3*y**8 + 15*k6**2*x**4*y**8 + k5**2*y**10 - \
   2*k4*k6*y**10 + 2*k5*k6*x*y**10 + 6*k6**2*x**2*y**10 + k6**2*y**12
   return poly


#output figures:
if evplot==1:
   if (legpos==-1111): legpos=2
   #fig1 = pl.figure(1,figsize=(3.5, 6))
   #plot11 = fig1.add_subplot(1,1,1)
   #ax1 = fig1.gca()
   #subplots_adjust(left=0.2 , bottom=0.1, right=0.95, top=0.95, wspace=0.2, hspace=0.1)
   fig1 = pl.figure(1,figsize=(3.7, 4.5))
   plot11 = fig1.add_subplot(1,1,1)
   ax1 = fig1.gca()
   subplots_adjust(left=0.23 , bottom=0.12, right=0.97, top=0.98, wspace=0.0, hspace=0.0)
   
   label1=[]
   leg1=[]
if (dtevplot==1):
   if (legpos==-1111): legpos=0
   label2=[]
   leg2=[]
   if (nfiles==1):
      fig2 = pl.figure(1,figsize=(6,4.5))
      plot21 = fig2.add_subplot(1,1,1)
      ax2 = fig2.gca()
      subplots_adjust(left=0.13 , bottom=0.13, right=0.95, top=0.96, wspace=0.2, hspace=0.2)
   else:
      fig2 = pl.figure(1,figsize=(5.5,3.8))
      plot21 = fig2.add_subplot(1,1,1)
      ax2 = fig2.gca()
      subplots_adjust(left=0.16, bottom=0.16, right=0.98, top=0.96, wspace=0.2, hspace=0.2)
if (vlasevplot==1 or collevplot==1):
   fig3 = pl.figure(1,figsize=(4.5, 5))
   label3=[]
   leg3=[]
   plot31 = fig3.add_subplot(1,1,1)
   ax3 = fig3.gca()
   subplots_adjust(left=0.18 , bottom=0.15, right=0.95, top=0.95, wspace=0.2, hspace=0.2)

#allocate contour arrays
res=100
if (plotRKC1==1) : ZRKC1 = np.zeros((res,res))
if (plotRKC2==1) : ZRKC2 = np.zeros((res,res))
if (plotRKC3==1) : ZRKC3 = np.zeros((res,res))
if (plotRKC4==1) : ZRKC4 = np.zeros((res,res))
if (plotRK2==1)  : ZRK2  = np.zeros((res,res))
if (plotRK3==1)  : ZRK3  = np.zeros((res,res))
if (plotRK4==1)  : ZRK4  = np.zeros((res,res))
if (plotRK4M==1) : ZRK4M = np.zeros((res,res))
   

dtsim_L=[]

par=Parameters()

if (evplot==1 or dtevplot==1):  #nlev settings
   print evplot, dtevplot
   for n in range(nfiles):
      filename=infiles[n]
      par.Read_Pars('parameters%s' %(fext[n]))
      n_nlev=float(par.pardict['n_nlev'])
      [filelen,nblocks,nnlev,nmodel,max_ev_ind]=parse_nlevfile(filename)
      if (verbose==1) :
         print filelen,nblocks,nnlev,max_ev_ind

      print nmodel,len(nmodel)
      nm=len(nmodel)

      #ntimes=(nblocks-1)/2+1 #divide nonlin. entries by 2 (linearized nl and nl_model)
      ntimes=int(np.max(nnlev)+1) #divide nonlin. entries by 2 (linearized nl and nl_model) or 3 with nlonly

      if nblocks==0:
        print "skipping file %s (no (%i) data blocks)" %(filename,nblocks)
        break
      else:
        print "reading %i data blocks" %(nblocks)
        print "containing %i time entries" %(ntimes)
        print "and maximal  %i ev entries..." %(max_ev_ind)
       
      evs_nlev=np.zeros((ntimes,max_ev_ind,2))+1000
      evs_model=np.zeros((ntimes,max_ev_ind,2,nm))+1000
      evs_nl_only=np.zeros((ntimes,max_ev_ind,2))+1000
      dt_sim=np.empty((ntimes))
      dt_nlev=np.empty((ntimes))
      dt_ExB=np.empty((ntimes))
      dt_model=np.empty((ntimes,nm))
      dt_combo=np.empty((ntimes))
      dt_nl_only=np.empty((ntimes))
      time=np.empty((ntimes))
      minreal_nlev =0.0
      maxreal_nlev =0.0
      maximag_nlev = 0.0
      maximag_lin = 0.0
      minreal_model = np.zeros((nm))
      maxreal_model = np.zeros((nm))
      maximag_model = np.zeros((nm))

      file=open(filename)
      lines=file.readlines()
      iblock = -1
      ev_ind = 0
      read_ev_nl = 0
      read_ev_model = 0
      dt_vlasov = -1
      
      for i in range(filelen):
         islinear=re.search(r'#linear(.*)dt_max\s*=\s*(\S*)\s*(\s*,\s*dt_vlasov\s*=\s*(\S*))?\s*\)',lines[i])
         if islinear:
            iblock+=1
            dt_max = float(islinear.group(2))
            if (islinear.group(4)): dt_vlasov = float(islinear.group(4))
            time[iblock]=0.0
            dt_nlev[iblock]= dt_max
            dt_ExB[iblock]= 3*dt_max
            dt_model[iblock,:]= dt_max
            dt_nl_only[iblock]= dt_max
            dt_sim[iblock] = dt_max
            ev_ind = 0
            read_ev_lin = 1
            read_ev_nl = 0
            read_ev_model = 0

         isnonlinear=re.search(r'#nonlinear.*time\s*(\S+).*=\s*\(\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*\)',lines[i])
         if (isnonlinear):
            ev_ind = 0
            is_nl_model=re.search(r'model',lines[i]) #lin+ev
            is_nl_model1=re.search(r'model.*vemax',lines[i]) #lin+ev
            is_nl_model2=re.search(r'model\s*lin',lines[i]) #lin+ev
            is_nlev=re.search(r'nlev',lines[i])
            read_ev_lin = 0
            if (is_nl_model):
               #settings for ev block read
               if(is_nl_model1):read_ev_model=1
               if(is_nl_model2):read_ev_model=2
               read_ev_nl = 0
               dt_model[iblock,read_ev_model-1]=float(isnonlinear.group(2))
               if (verbose==1) :
                  #in the second block we should have collected all dt's:
                  print "dt_nlev = %s\n" %(str(dt_nlev[iblock]))
                  print "dt_ExB = %s\n" %(str(dt_ExB[iblock]))
                  print "dt_model = %s\n" %(str(dt_model[iblock,read_ev_model-1]))

            if (is_nlev):
               #first block:  nlev
               iblock+=1
               time[iblock]=float(isnonlinear.group(1))
               dt_nlev[iblock]=float(isnonlinear.group(2))
               dt_ExB[iblock]=float(isnonlinear.group(3))
               dt_sim[iblock]=float(isnonlinear.group(4))
               #dt_nl is dt_model in this case
               
               if (auto_tmin==0): 
                  if itmin==1 and time[iblock]>=tmin:
                     itmin=iblock #substract linear
               #else leave tmin= 0 and i11tmin= 0
               if time[iblock]<=tmax:
                  itmax=iblock+1
               if (verbose==1): print "reading time[%i] = %f\n" %(iblock,time[iblock])

               #settings for ev block read
               read_ev_model = 0
               read_ev_nl = 1

            ev_ind = 0
         
         iseigenvalue=re.match(r'(?![#])\s*(\S+)\s+(\S+)\s*$',lines[i])
         if (iseigenvalue):
            real = float(iseigenvalue.group(1))
            imag = float(iseigenvalue.group(2))
            if (read_ev_lin==1):
               maximag_lin = np.max([maximag_lin,np.abs(imag)])
               if (ev_ind<max_ev_ind):
                  evs_nlev[0,ev_ind,0]=real
                  evs_nlev[0,ev_ind,1]=imag
               if (verbose==1): print 'lin',ev_ind,evs_nlev[0,ev_ind,:]
            
            if (read_ev_nl==1):
               minreal_nlev = np.min([minreal_nlev,real])
               maxreal_nlev = np.max([maxreal_nlev,real])
               maximag_nlev = np.max([maximag_nlev,np.abs(imag)])
               if (ev_ind<max_ev_ind):
                  evs_nlev[iblock,ev_ind,0]=real
                  evs_nlev[iblock,ev_ind,1]=imag
                  if (verbose==1): print 'nlev',iblock,ev_ind,evs_nlev[iblock,ev_ind,:]
               else:
                  if (verbose==1): print 'skipping ev',real,imag
            
            if (read_ev_model>0):
               minreal_model[read_ev_model-1] = np.min([minreal_model[read_ev_model-1],real])
               maxreal_model[read_ev_model-1] = np.max([maxreal_model[read_ev_model-1],real])
               maximag_model[read_ev_model-1] = np.max([maximag_model[read_ev_model-1],np.abs(imag)])
               if (ev_ind<max_ev_ind):
                  evs_model[iblock,ev_ind,0,read_ev_model-1]=real
                  evs_model[iblock,ev_ind,1,read_ev_model-1]=imag
                  if (verbose==1): print 'mod',iblock,ev_ind,read_ev_model,evs_model[iblock,ev_ind,:,read_ev_model-1]
            ev_ind +=1

      print itmin,itmax,ntimes
      print time
      mindt_nlev=np.min(dt_nlev)
      mindt=np.min(dt_sim)
      mindt_model = 10000000
      for pm in plotmodel:
         tmp=np.min(dt_model[:,pm-1])
         if (tmp>0):
            mindt_model=np.min([tmp,mindt_model])


      mindt_ExB=np.min(dt_ExB)
      
      mindt_tot = mindt_nlev

      for pm in plotmodel:
         print "max imag model %d = %f" %(pm,maximag_model[pm-1])
      print "max imag nlev = %f" %(maximag_nlev)
      print "max imag lin = %f" %(maximag_lin)
      print "dt_max = %f" %(dt_max)
      if (dt_vlasov>0): print 'dt_vlasov = ',dt_vlasov
      if (dt_vlasov>0): dt_max=dt_vlasov
      print "min dt_nlev = %f" %(mindt_nlev)
      print "min dt_model = %f" %(mindt_model)
      print "min dt_sim = %f" %(mindt)
      print "min dt_ExB = %f" %(mindt_ExB)

      if (mindt_model>0):
          dtscale=mindt_model
      else:
          dtscale=mindt
      #for pm in plotmodel:
      #   print dt_model[:,pm-1]

      minreal=-3.5 #looks nice for RK4
      maxreal=0.25 #looks nice for RK4
      tmp = -1
      for pm in plotmodel:
         print 'plot model',plotmodel
         tmp=np.max([maximag_model[pm-1],tmp])
      maximag=1.06*np.max([maximag_nlev,tmp])*dtscale
     
      print minreal, maxreal, maximag

      x=np.linspace(minreal,maxreal,num=res)
      y=np.linspace(-maximag,maximag,num=res)
      for i in range(res):
         for j in range(res):
            #WARNING: we need the transposed Z
            if (plotRKC1==1) : ZRKC1[j,i]=float(RKpoly(x[i],y[j],RKC1_taylor))
            if (plotRKC2==1) : ZRKC2[j,i]=float(RKpoly(x[i],y[j],RKC2_taylor))
            if (plotRKC3==1) : ZRKC3[j,i]=float(RKpoly(x[i],y[j],RKC3_taylor))
            if (plotRKC4==1) : ZRKC4[j,i]=float(RKpoly(x[i],y[j],RKC4_taylor))
            if (plotRK2==1)  : ZRK2[j,i]=float(RKpoly(x[i],y[j],RK2_taylor))
            if (plotRK3==1)  : ZRK3[j,i]=float(RKpoly(x[i],y[j],RK3_taylor))
            if (plotRK4==1)  : ZRK4[j,i]=float(RKpoly(x[i],y[j],RK4_taylor))
            if (plotRK4M==1) : ZRK4M[j,i]=float(RKpoly(x[i],y[j],RK4M_taylor))

      for iblock in range(1,ntimes):
         if (time[iblock] >= tmin and time[iblock] <=tmax and evplot!=0):
            if (verbose==1): print "plotting time = %f\n" %(time[iblock])
            lw2 = 3
            lw = 2
            sz = 10
            fsz = 15
            lsz = 14
            msz = 9
            font= FontProperties(size=lsz)
            style=['m^','bo','rd','bo','ks','bo']
            plot11.set_xlabel(r'Re($\lambda\, \Delta t)$ ,$\Delta t=$'+str(0.0001*np.floor(10000*dtscale)),fontsize=fsz)
            plot11.set_ylabel(r'Im($\lambda\, \Delta t)$',fontsize=fsz)
     
            plot11.set_ylim([-maximag,maximag])
            plot11.set_xlim([minreal,maxreal])
     
            plot11.xaxis.set_tick_params ( which='major', labelsize=lsz ) 
            plot11.yaxis.set_tick_params ( which='major', labelsize=lsz ) 
            plot11.xaxis.set_tick_params(width=1)
            plot11.set_xticks([-3.5,-2.5,-1.5,-0.5])
           
            
            plot11.plot(evs_nlev[0,:,0]*dtscale,evs_nlev[0,:,1]*dtscale,'g^',markeredgewidth=lw,markersize=msz)
            plot11.plot(evs_nlev[iblock,:,0]*dtscale,evs_nlev[iblock,:,1]*dtscale,'bo',markeredgewidth=lw,markersize=msz)
            for pm in plotmodel:
               plot11.plot(evs_model[iblock,:,0,pm-1]*dtscale,evs_model[iblock,:,1,pm-1]*dtscale,'rs',markeredgewidth=lw,markersize=msz)
            plot11.plot(x,0*x,'k--')
            plot11.plot(0*y,y,'k--')
            if (plotRKC1==1): plot11.contour(x,y,ZRKC1,levels=[1],cmap=get_cmap('jet'),linewidths=lw2)
            if (plotRKC2==1): plot11.contour(x,y,ZRKC2,levels=[1],cmap=cm.gray,linewidths=lw2)
            if (plotRKC3==1): plot11.contour(x,y,ZRKC3,levels=[1],cmap=get_cmap('jet'),linewidths=lw2)
            if (plotRKC4==1): plot11.contour(x,y,ZRKC4,levels=[1],cmap=get_cmap('jet'),linewidths=lw2)
            if (plotRK2==1) : plot11.contour(x,y,ZRK2,levels=[1],cmap=get_cmap('PRgn'),linewidths=lw2)
            if (plotRK3==1) : plot11.contour(x,y,ZRK3,levels=[1],cmap=get_cmap('summer'),linewidths=lw2)
            if (plotRK4==1) : plot11.contour(x,y,ZRK4,levels=[1],cmap=get_cmap('PRgn'),linewidths=lw2)
            if (plotRK4M==1): plot11.contour(x,y,ZRK4M,levels=[1],cmap=get_cmap('rainbow'),linewidths=lw2)
            leg1=[r'$\lambda_L$','$\lambda_{\mathrm{nlev}}$',r'$\lambda_L+\lambda_{\chi,\mathrm{max}}$']
            if (legpos>=0): legend(leg1,loc=legpos,prop=font)
     
            max_nlev_i=np.max(abs(evs_nlev[iblock,:,1]))*dtscale
            mid=-(max_nlev_i+maximag_lin*dtscale)/2.
            if (annotate==1):
               plot11.annotate('', xy=(-0.3,-max_nlev_i), xytext=(-0.3,-maximag_lin*dtscale),
                     arrowprops=dict(edgecolor='black',facecolor='white', frac=0.4,linewidth=2))
               plot11.annotate('nonlinear\n shift', xy=(-1.8, mid), xytext=(-1.8, mid),fontsize=14)
               plot11.annotate('RK4', xy=(-2.6, 0.2), xytext=(-2.6, 0.2),fontsize=fsz)
               #show()
            if (annotate==3):
               plot11.annotate('', xy=(-0.3,-max_nlev_i), xytext=(-0.3,-maximag_lin*dtscale),
                     arrowprops=dict(edgecolor='black',facecolor='white', frac=0.1,linewidth=2))
               plot11.annotate('nonlinear\n shift', xy=(-1.8, mid), xytext=(-1.8, mid),fontsize=fsz)
               plot11.annotate('RK4M', xy=(-2.99, 0.2), xytext=(-2.99, 0.2),fontsize=fsz)
               #show()
            if (annotate==2):
               plot11.annotate('', xy=(-0.4,-1.66), xytext=(-0.4,-0.7),
                     arrowprops=dict(edgecolor='black',facecolor='white', shrink=.0,linewidth=2))
               plot11.annotate('nonlinear\n shift', xy=(-1.8, -1.2), xytext=(-1.8, -1.2),fontsize=fsz)
            savefig('nlev%s_%.4i.png'%(fextout[n],iblock))   
            savefig('nlev%s_%.4i.pdf'%(fextout[n],iblock))   
            plot11.clear()

      if (dtevplot==1):
         font= FontProperties(size=14);
         fsz = 16
         lw=2
         msz = 6
         if (ntimes>1):
            #set [tmin,tmax] (if not user-specified)
            if tmin==0:
               if time[1]-time[0] > 0.1*(time[-1]-time[1]):
                  if (verbose==1):
                     print "set tmin to time(1)"
                  tmin=time[1]
               else:
                  tmin=min(time)
            if tmax==999999999:
               tmax=max(time)

         style=['k-','k--','k--','4-','x-','<-']
         dt_L = dt_max
         legpm=['','','']
         if (dt2ev==-1):
             name = 'ev'
             legnlev= r'nlev: $\lambda_{\mathrm{nlev}}$'
             leglin = r'linear: $\lambda_{L}$'
             legsim = r'simulation: $\Delta t_L/\Delta t$'
             legpm[1]=r'estimate: $\lambda_{\chi,\mathrm{max}}$'
             legpm[2]=r'estimate: $\lambda_L+\lambda_{\chi,\mathrm{max}}$'
             plot21.set_ylabel(r'max. eigenvalue $\lambda/\lambda_L$',fontsize=fsz)
         if (dt2ev==1):
             name = 'dt'
             legnlev= r'nlev: $\Delta t_{\mathrm{nlev}}$'
             leglin = r'nlev: $\Delta t_{L}$'
             legsim = r'simulation: $\Delta t/\Delta t_L$'
             legpm[1]=r'estimate: $\Delta t_{\mathrm{comb}}$'
             legpm[2]=r'estimate: $C_{RK}/(\lambda_L+\lambda_{\chi,\mathrm{max}})$'
             plot21.set_ylabel(r'$(\Delta t/\Delta t_L)$',fontsize=fsz)
         
         if (nfiles==1):
            if (ymax==-1): ymax=1.5*dt_max/dt_L
            if (ymin==-1): ymin=0
            plot21.plot(\
            time[itmin:itmax],np.ones(itmax-itmin)*dt_max/dt_L,'g^-',\
            time[itmin:itmax],(dt_nlev[itmin:itmax]/dt_L)**(dt2ev),'bo-',\
            #make less dense points:
            time[itmin:itmax],((1./(1./dt_ExB[itmin:itmax]+1./dt_L))/dt_L)**(dt2ev),'r--s',\
            time[itmin:itmax],(dt_sim[itmin:itmax]/dt_L)**(dt2ev),'k-',\
            linewidth=lw,markersize=msz)
            leg1=[leglin,legnlev]
            leg1.append(legpm[2])
            leg1.append(legsim)
            for pm in plotmodel:
               if (pm>0 and pm<=2):
                  print "here",(dt_model[itmin:itmax,pm-1]/dt_L)**dt2ev
                  plot21.plot(time[itmin:itmax],(dt_model[itmin:itmax,pm-1]/dt_L)**(dt2ev),style[pm+2],\
                  linewidth=lw)
                  leg1.append(legpm[pm])
            plot21.set_ylim([ymin,ymax])
         else:
            if (ymax==-1): ymax=0.55
            if (ymin==-1): ymin=0.25
            #dt_L=0.276E-2   #set this manually
            dt_L=5.18E-03    #dt_vlasov!!
            name=name+'stability'
            dtsim_L.append(dt_sim[-1]/dt_L)
            if n==0: 
               # this is for the paper: timetrace with fixed dt and instabilities
               plot21.plot(time[itmin:itmax],dt_nlev[itmin:itmax]/dt_L,'o-',\
               linewidth=lw)
               leg1=[r'$\Delta t_{\rm nlev}$',\
                     r'$\Delta t_{\mathrm{sim}}=$%.2f$\Delta t _L$'%(dt_max/dt_L)]
               plot21.plot(time[itmin:itmax],dt_sim[itmin:itmax]/dt_L,style[n],linewidth=lw+1)
               plot21.set_ylim([ymin,ymax])
            else:
               ##from nrg files that are added to the argument list...
               if n==2: time[itmax-1]=378.090111 ##insert this manually
               if n==1: time[itmax-1]=379.872111 ##insert this manually
               plot21.plot(time[itmin:itmax],dt_sim[itmin:itmax]/dt_L,style[n],linewidth=lw)
               plot21.plot(time[itmax-1],dt_sim[itmax-1]/dt_L,'rx',markersize=10,markeredgewidth=4)
         plot21.set_xlabel(r'time $[c_s/L_{\mathrm{ref}}]$',fontsize=fsz)
         plot21.set_xlim([tmin,tmax])

   if (dtevplot==1):
      if (legpos>=0):
         legend(leg1,bbox_to_anchor=(0, legyoff, 1, 1-legyoff),loc=legpos,prop=font)
      if (annotate==1):
         print dtsim_L
         plot21.annotate(r'$\Delta t_{\mathrm{sim}}=%.2f\Delta t_L$'%(dtsim_L[-1]), xy=(373, dtsim_L[-1]), xytext=(374, 0.47),
               arrowprops=dict(edgecolor='black',facecolor='white', shrink=.001),fontsize=fsz)
         plot21.annotate(r'$\Delta t_{\mathrm{sim}}=%.2f\Delta t_L$'%(dtsim_L[-2]), xy=(374, dtsim_L[-2]), xytext=(377, 0.27),
               arrowprops=dict(edgecolor='black',facecolor='white', shrink=.001),fontsize=fsz)

      savefig('timetrace_%s%s.pdf'%(name,fextout[n]))   
      show()

if (collevplot==1 or vlasevplot==1):
   for n in range(nfiles):
      filename=infiles[n]
      filelen=file_len(filename)
      evs=np.zeros((filelen,2))
      file=open(filename)
      lines=file.readlines()
      for i in range(filelen):
         iseigenvalue=re.match(r'(?![#])\s*(\S+)\s*(\S+)\s*$',lines[i])
         if iseigenvalue:
            evs[i,0]=iseigenvalue.group(1)
            evs[i,1]=iseigenvalue.group(2)
      file.close()

      font= FontProperties(size=14);
      lw = 3
      lw2 = 2
      sz = 12
      fsz = 16
      lsz = 14
 

      minreal_ev=np.min(evs[:,0])
      maxreal_ev=np.max(evs[:,0])
      maximag_ev=np.max(np.abs(evs[:,1]))

      betaRKC4=31.039324044870003
      betaRK4M=4.8469
      betaRK4=2.78
      betaRK=betaRK4
      if (plotRK4M): betaRK=betaRK4M

      dt_real=0.0001*np.floor(10000*abs(betaRKC4)/abs(minreal_ev))
      dt_imag = abs(betaRK/abs(maximag_ev))
      

      print 'minreal= ', minreal_ev, 'dt_real= ', dt_real
      print 'maximag= ', maximag_ev, 'dt_imag= ', dt_imag
      print 'dt_max from ev = ', dt_max
      print 'stability boundary RKC: ', -betaRKC4
      print 'stability boundary RK: ', betaRK
      
      if (dt_max ==-1): dt_max = min(dt_real,dt_imag)

      if (collevplot==1): 
          name='coll'; floorfk=1;nticks=5
          minreal=np.min([-3.5,minreal_ev*1.05*dt_max])
          maxreal = -minreal/32
      else:
          name='vlas'; floorfk=10000;nticks=4
          minreal = minreal_ev*1.3*dt_max
          maxreal = maxreal_ev*2*dt_max
      if (xmin != 1111 ): minreal=xmin
      if (xmax != 1111 ): maxreal=xmax
      
      #this holds for coll and vlas. evs
      if (ymax==-1): 
         maximag=max(1.1*maximag_ev*dt_max,4)
      else: 
         maximag=ymax

      if (minreal !=1111 or maxreal!=1111):
         plot31.set_xticks(np.linspace(-np.floor(-floorfk*minreal)/floorfk,0,num=nticks))
      else:
         plot31.set_xticks([-3.5,-2.5,-1.5,-0.5])

      res=100
      x=np.linspace(minreal,maxreal,num=res)
      y=np.linspace(-maximag,maximag,num=res)

      for i in range(res):
         for j in range(res):
            #WARNING: we need the transposed Z
            if (plotRKC1==1) : ZRKC1[j,i]=float(RKpoly(x[i],y[j],RKC1_taylor))
            if (plotRKC2==1) : ZRKC2[j,i]=float(RKpoly(x[i],y[j],RKC2_taylor))
            if (plotRKC3==1) : ZRKC3[j,i]=float(RKpoly(x[i],y[j],RKC3_taylor))
            if (plotRKC4==1) : ZRKC4[j,i]=float(RKpoly(x[i],y[j],RKC4_taylor))
            if (plotRK2==1)  : ZRK2[j,i]=float(RKpoly(x[i],y[j],RK2_taylor))
            if (plotRK3==1)  : ZRK3[j,i]=float(RKpoly(x[i],y[j],RK3_taylor))
            if (plotRK4==1)  : ZRK4[j,i]=float(RKpoly(x[i],y[j],RK4_taylor))
            if (plotRK4M==1) : ZRK4M[j,i]=float(RKpoly(x[i],y[j],RK4M_taylor))

      style=['m^','bo','rd','bo','ks','bo']
      plot31.plot(evs[:,0]*dt_max,evs[:,1]*dt_max,'rx',linewidth=lw2)
      plot31.plot(x,0*x,'k--')
      plot31.plot(0*y,y,'k--')
      if (plotRKC1==1): plot31.contour(x,y,ZRKC1,levels=[1],cmap=get_cmap('autumn'),linewidths=lw2)
      if (plotRKC2==1): plot31.contour(x,y,ZRKC2,levels=[1],cmap=cm.gray,linewidths=lw2)
      if (plotRKC3==1): plot31.contour(x,y,ZRKC3,levels=[1],cmap=get_cmap('jet'),linewidths=lw2)
      if (plotRKC4==1): plot31.contour(x,y,ZRKC4,levels=[1],cmap=get_cmap('summer'),linewidths=lw2)
      if (plotRK2==1) : plot31.contour(x,y,ZRK2,levels=[1],cmap=get_cmap('summer'),linewidths=lw2)
      if (plotRK3==1) : plot31.contour(x,y,ZRK3,levels=[1],cmap=get_cmap('summer'),linewidths=lw2)
      if (plotRK4==1) : plot31.contour(x,y,ZRK4,levels=[1],cmap=get_cmap('PRgn'),linewidths=lw2)
      if (plotRK4M==1): plot31.contour(x,y,ZRK4M,levels=[1],cmap=get_cmap('autumn'),linewidths=lw2)

      plot31.set_xlabel(r'Re($\lambda\, \Delta t)$',fontsize=fsz)
      plot31.set_ylabel(r'Im($\lambda\, \Delta t)$',fontsize=fsz)
      
      plot31.set_ylim([-maximag,maximag])
      plot31.set_xlim([minreal,maxreal])
      
      plot31.xaxis.set_tick_params ( which='major', labelsize=lsz ) 
      plot31.yaxis.set_tick_params ( which='major', labelsize=lsz ) 
      leg1=[r'scaled ev: $\Delta t$='+str(0.0001*np.floor(10000*dt_max))]
      if (legpos>=0): legend(leg1,loc=legpos,prop=font)
      if (annotate==1):
         plot31.annotate('', xy=(-0.0021, 2), xytext=(0, 2.008),
               arrowprops=dict(edgecolor='black',facecolor='white', shrink=.0))
         plot31.annotate('num. damping', xy=(-0.0022, 2), xytext=(-0.0022, 2.4))
            
      savefig('evplot_%s_%d%s.eps'%(name,annotate,fextout[n]))   
      savefig('evplot_%s_%d%s.png'%(name,annotate,fextout[n]))   
      show()

