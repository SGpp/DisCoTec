#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 09:32:25 2013

@author: dtold
"""

from pylab import *
from sys import argv,exit,stdout
import optparse as op

nr=100;ntheta=150
parser=op.OptionParser(description='Extract Miller shaping parameters from EQDSK files.')
parser.add_option('--rovera','-r',action='store_const',const=1)
parser.add_option('--conv','-c',action='store_const',const=1)
options,args=parser.parse_args()
use_r_a=options.rovera
write_rhotor_rhopol_conversion=options.conv
if not write_rhotor_rhopol_conversion:
    if len(args)!=2:
        exit("""
Please give two arguments: <EQDSK filename> <Position in rho_tor>
optional: -r <Position in r/a>
-c: Write rhotor/rhopol conversion to file\n""")
    
filename=args[0]
file=open(filename,'r')

def find(val,arr):
    return argmin(abs(arr-val))


    
eqdsk=file.readlines()
print 'Header: %s' %eqdsk[0]
#set resolutions
nw=int(eqdsk[0].split()[-2]);nh=int(eqdsk[0].split()[-1])
pw=(nw/8/2)*2 #psi-width, number of flux surfaces around position of interest
print 'Resolution: %d x %d' %(nw,nh)

entrylength=16
try:
    rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
except:
    entrylength=15
    try:
        rdim,zdim,rctr,rmin,zmid=[float(eqdsk[1][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[1])/entrylength)]
    except:
        exit('Error reading EQDSK file, please check format!')
rmag,zmag,psiax,psisep,Bctr=[float(eqdsk[2][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[2])/entrylength)]
dum,psiax2,dum,rmag2,dum=[float(eqdsk[3][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[3])/entrylength)]
zmag2,dum,psisep2,dum,dum=[float(eqdsk[4][j*entrylength:(j+1)*entrylength]) for j in range(len(eqdsk[4])/entrylength)]
if rmag!=rmag2: sys.exit('Inconsistent rmag: %7.4g, %7.4g' %(rmag,rmag2))
if psiax2!=psiax: sys.exit('Inconsistent psiax: %7.4g, %7.4g' %(psiax,psiax2))
if zmag!=zmag2: sys.exit('Inconsistent zmag: %7.4g, %7.4g' %(zmag,zmag2) )
if psisep2!=psisep: sys.exit('Inconsistent psisep: %7.4g, %7.4g' %(psisep,psisep2))
F=empty(nw,dtype=float)
p=empty(nw,dtype=float)
ffprime=empty(nw,dtype=float)
pprime=empty(nw,dtype=float)
qpsi=empty(nw,dtype=float)
psirz_1d=empty(nw*nh,dtype=float)
start_line=5
lines=range(nw/5)
if nw%5!=0: lines=range(nw/5+1)
for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    F[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    p[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    ffprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    pprime[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

lines_twod=range(nw*nh/5)
if nw*nh%5!=0: lines_twod=range(nw*nh/5+1)
for i in lines_twod:
    n_entries=len(eqdsk[i+start_line])/entrylength
    psirz_1d[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1
psirz=psirz_1d.reshape(nh,nw)

for i in lines:
    n_entries=len(eqdsk[i+start_line])/entrylength
    qpsi[i*5:i*5+n_entries]=[float(eqdsk[i+start_line][j*entrylength:(j+1)*entrylength]) for j in range(n_entries)]
start_line=i+start_line+1

#invert sign of psi if necessary to guarantee increasing values for interpolation
if psisep<psiax:
    psirz=-psirz
    ffprime=-ffprime
    pprime=-pprime
    psiax*=-1
    psisep*=-1

#ignore limiter data etc. for the moment
dw=rdim/(nw-1)
dh=zdim/(nh-1)
rgrid=array([rmin+i*dw for i in range(nw)])
zgrid=array([zmid-zdim/2.+i*dh for i in range(nh)])
#contourf(rgrid,zgrid,psirz,70);gca().set_aspect('equal')
#show()

#create 5th order 2D spline representation of Psi(R,Z)
from scipy.interpolate import RectBivariateSpline as RBS
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline as US
interpol_order=3
psi_spl=RBS(zgrid,rgrid,psirz,kx=interpol_order,ky=interpol_order)

#linear grid of psi, on which all 1D fields are defined
linpsi=linspace(psiax,psisep,nw)
#create rho_tor grid
x_fine=linspace(psiax,psisep,nw*10)
phi_fine=empty((nw*10),dtype=float)
phi_fine[0]=0.
#spline of q for rhotor grid
q_spl_psi=US(linpsi,qpsi,k=interpol_order,s=1e-5)

for i in range(1,nw*10):
    x=x_fine[:i+1]
    y=q_spl_psi(x)
    phi_fine[i]=trapz(y,x)
rho_tor_fine=sqrt(phi_fine/phi_fine[-1])
rho_tor_spl=US(x_fine,rho_tor_fine,k=interpol_order,s=1e-5)
rho_tor=empty(nw,dtype=float)
for i in range(nw):
    rho_tor[i]=rho_tor_spl(linpsi[i])

if write_rhotor_rhopol_conversion:
    rt_rp_filename='rt_rp_%s' %filename
    rt_rp_file=open(rt_rp_filename,'w')
    rt_rp_file.write('# rho_tor          rho_pol\n')
    for i in range(len(x_fine)):
        rho_pol=sqrt((x_fine[i]-psiax)/(psisep-psiax))
        rt_rp_file.write('%16.8e %16.8e\n' %(rho_tor_fine[i],rho_pol))
    rt_rp_file.close()
    exit('\nWrote rhotor/rhopol conversion to %s.' %rt_rp_filename)

t1=arctan2(zmid-zdim/2.-zmag,rmin-rmag)
t2=arctan2(zmid-zdim/2.-zmag,rmin+rdim-rmag)
t3=arctan2(zmid+zdim/2.-zmag,rmin+rdim-rmag)
t4=arctan2(zmid+zdim/2.-zmag,rmin-rmag)

theta_arr=linspace(-pi,pi,ntheta)
#for i in range(nw):
#    curr_psi=linpsi[i]

print 'Finding flux surface shapes...'
R=empty((nw,ntheta),dtype=float)
Z=empty((nw,ntheta),dtype=float)
dr=rdim*cos(theta_arr)
dz=rdim*sin(theta_arr)
for j in range(len(theta_arr)):
    stdout.write('\r Finished %4.1f%%.' %(j*100./(ntheta-1)))
    stdout.flush()
    theta=theta_arr[j]
    r_pol=linspace(rmag,rmag+dr[j],nr) #array([rmag+i*dr for i in range(nr)])
    z_pol=linspace(zmag,zmag+dz[j],nr) #array([zmag+i*dz for i in range(nr)])
    psi_rad=psi_spl.ev(z_pol,r_pol)
    psi_rad[0]=psiax
    #must restrict interpolation range because of non-monotonic psi around coils
    cutoff=0
    for i in range(1,len(psi_rad)):
        if psi_rad[i]<psi_rad[i-1]:
            cutoff=i
            break
    psi_rad=psi_rad[:i]
    end_ind=argmin(abs(psi_rad-psisep))
    end_ind+=(1 if (psi_rad[end_ind]<psisep) else 0)
    indsep=end_ind+1
    R_int=interp1d(psi_rad[:indsep],r_pol[:indsep],kind=interpol_order)
    R[:,j]=R_int(linpsi)
    Z_int=interp1d(psi_rad[:indsep],z_pol[:indsep],kind=interpol_order)
    Z[:,j]=Z_int(linpsi)

print '\nFinding flux surface centers...'
#find average elevation for all flux surfaces
Z_avg=empty(nw,dtype=float)
ds=empty(ntheta,dtype=float)
for i in range(nw):
    ds[1:ntheta-1]=0.5*sqrt((R[i,2:ntheta]-R[i,0:ntheta-2])**2+(Z[i,2:ntheta]-Z[i,0:ntheta-2])**2)
    ds[0]=0.5*sqrt((R[i,1]-R[i,-1])**2+(Z[i,1]-Z[i,-1])**2)
    ds[-1]=0.5*sqrt((R[i,0]-R[i,-2])**2+(Z[i,0]-Z[i,-2])**2)
    Z_avg[i]=average(Z[i,:],weights=ds)

#find R0 for all flux surfaces
R0=empty(nw,dtype=float)
R0[0]=rmag
r_avg=empty(nw,dtype=float)
r_avg[0]=0.
r_maxmin=empty(nw,dtype=float)
r_maxmin[0]=0.
for i in range(1,nw):
    stdout.write('\r Finished %4.1f%%.' %(i*100./(nw-1)))
    stdout.flush()
    R_array=R[i,ntheta/4:3*ntheta/4]
    Z_array=Z[i,ntheta/4:3*ntheta/4]
    #low field side
    Z_int=interp1d(Z_array,range(ntheta/2),kind=interpol_order)
    ind_Zavg=Z_int(Z_avg[i])
    R_int=interp1d(range(ntheta/2),R_array,kind=interpol_order)
    R_out=R_int(ind_Zavg)
    R_max=amax(R_array)
    #high field side
    R_array=roll(R[i,:-1],ntheta/2)[ntheta/4:3*ntheta/4]
    Z_array=roll(Z[i,:-1],ntheta/2)[ntheta/4:3*ntheta/4]

    #have to use negative Z_array here to have increasing order
    Z_int=interp1d(-Z_array,range(ntheta/2),kind=interpol_order)
    #again negative
    ind_Zavg=Z_int(-Z_avg[i])

    R_int=interp1d(range(ntheta/2),R_array,kind=interpol_order)
    R_in=R_int(ind_Zavg)
    R_min=amin(R_array)
    R0[i]=0.5*(R_out+R_in)
    r_avg[i]=0.5*(R_out-R_in)
    r_maxmin[i]=0.5*(R_max-R_min)


radpos=float(args[1])
if use_r_a:
    r_a=radpos
    #find psi index of interest (for the specified r/a position)
    poi_ind=find(radpos,r_avg/r_avg[-1])
else:
    #find psi index of interest (for the specified rho_tor position)
    poi_ind=find(radpos,rho_tor)


#modified theta grid for each flux surface
#arrays equidistant on modified theta grid are marked by 'tm' index!!!
linpsi_spl=US(r_avg,linpsi,k=interpol_order,s=1e-5)
ravg_spl=US(linpsi,r_avg,k=interpol_order,s=1e-5)
if not use_r_a:
    ravg_rho_spl=US(rho_tor,r_avg/r_avg[-1],k=interpol_order,s=1e-5)
    r_a=ravg_rho_spl(radpos)
    print '\nExamine %d flux surfaces around position r/a=%7.4g...' %(pw,r_a)
else:
    print '\nExamine %d flux surfaces around position rho_tor=%7.4g...' %(pw,radpos)
rmaxmin_spl=US(linpsi,r_maxmin,k=interpol_order,s=1e-5)
q_spl=US(r_avg,qpsi,k=interpol_order,s=1e-5)
R0_spl=US(r_avg,R0,k=interpol_order,s=1e-5)
F_spl=US(r_avg,F,k=interpol_order,s=1e-5)
r=r_a*r_avg[-1]
psi=linpsi_spl(r)
psi_N=(psi-psiax)/(psisep-psiax)
R0_pos=R0_spl(r)
F_pos=F_spl(r)
Bref_miller=F_pos/R0_pos



print 'Coordinates: r=%8.5g, psi=%8.5g, psi_N=%8.5g, r/R0=%8.5g, rho_tor=%8.5g, r_maxmin=%8.5g' %(r,psi,psi_N,r/R0_pos,rho_tor_spl(psi),rmaxmin_spl(psi))
psi_stencil=range(poi_ind-pw/2,poi_ind+pw/2)
if psi_stencil[0]<1: psi_stencil=[psi_stencil[i]+1-psi_stencil[0] for i in range(len(psi_stencil))]
if psi_stencil[-1]>nw-1: psi_stencil=[psi_stencil[i]-(psi_stencil[-1]-nw+1) for i in range(len(psi_stencil))]
R_tm=empty((pw,ntheta),dtype=float)
Z_tm=empty((pw,ntheta),dtype=float)
R_extended=empty(2*ntheta-1,dtype=float)
Z_extended=empty(2*ntheta-1,dtype=float)
#theta_mod[0]=theta_arr
#R_tm[0]=R[0]
#Z_tm[0]=Z[0]
theta_tmp=linspace(-2.*pi,2*pi,2*ntheta-1)

print 'Interpolating to flux-surface dependent (proper) theta grid...'
for i in psi_stencil:
    stdout.write('\r Finished %4.1f%%.' %(psi_stencil.index(i)*100./(len(psi_stencil)-1)))
    stdout.flush()
    imod=i-psi_stencil[0]
    #print 'Finished %4.1f%%.' %(float(i)/(pw-1)*100)   
    R_extended[0:(ntheta-1)/2]=R[i,(ntheta+1)/2:-1]
    R_extended[(ntheta-1)/2:(3*ntheta-3)/2]=R[i,:-1]
    R_extended[(3*ntheta-3)/2:]=R[i,0:(ntheta+3)/2]
    Z_extended[0:(ntheta-1)/2]=Z[i,(ntheta+1)/2:-1]
    Z_extended[(ntheta-1)/2:(3*ntheta-3)/2]=Z[i,:-1]
    Z_extended[(3*ntheta-3)/2:]=Z[i,0:(ntheta+3)/2]
    #for j in range(ntheta):
    theta_mod_ext=arctan2(Z_extended-Z_avg[i],R_extended-R0[i])
    #introduce 2pi shifts to theta_mod_ext
    for ind in range(ntheta):
        if theta_mod_ext[ind+1]<0. and theta_mod_ext[ind]>0. and abs(theta_mod_ext[ind+1]-theta_mod_ext[ind])>pi:
            lshift_ind=ind
        if theta_mod_ext[-ind-1]>0. and theta_mod_ext[-ind]<0. and abs(theta_mod_ext[-ind-1]-theta_mod_ext[-ind])>pi:
            rshift_ind=ind
    theta_mod_ext[-rshift_ind:]+=2.*pi            
    theta_mod_ext[:lshift_ind+1]-=2.*pi
    #print theta_mod, theta_arr
#    plot(theta_mod_ext)
#    plot(theta_tmp)
#    show()
    theta_int=interp1d(theta_mod_ext,theta_tmp,kind=interpol_order)
    theta_orig_tm=theta_int(theta_arr)
    R_int=interp1d(theta_mod_ext,R_extended,kind=interpol_order)
    Z_int=interp1d(theta_mod_ext,Z_extended,kind=interpol_order)
    R_tm[imod]=R_int(theta_arr)
    Z_tm[imod]=Z_int(theta_arr)
#    plot(R_tm[imod],Z_tm[imod])
#gca().set_aspect('equal')

#now we have the flux surfaces on a symmetric grid in theta (with reference to R0(r), Z0(r))
#symmetrize flux surfaces
#figure()
R_sym=empty((pw,ntheta),dtype=float)
Z_sym=empty((pw,ntheta),dtype=float)
for i in psi_stencil:
    imod=i-psi_stencil[0]
    Z_sym[imod,:]=0.5*(Z_tm[imod,:]-Z_tm[imod,::-1])+Z_avg[i]
    R_sym[imod,:]=0.5*(R_tm[imod,:]+R_tm[imod,::-1])
#    plot(R_sym[imod],Z_sym[imod])
#gca().set_aspect('equal')
#show()
figure()
dq_dr_avg=empty(pw,dtype=float)
dq_dpsi=empty(pw,dtype=float)
drR=empty(pw,dtype=float)
kappa=empty(pw,dtype=float)
delta=empty(pw,dtype=float)
s_kappa=empty(pw,dtype=float)
s_delta=empty(pw,dtype=float)
delta_upper=empty(pw,dtype=float)
delta_lower=empty(pw,dtype=float)
zeta_arr=empty((pw,4),dtype=float)
zeta=empty(pw,dtype=float)
s_zeta=empty(pw,dtype=float)
for i in psi_stencil:
    imod=i-psi_stencil[0]
    #calculate delta
    stencil_width=ntheta/10
    for o in range(2):
        if o: 
            ind=argmax(Z_sym[imod])
            section=range(ind+stencil_width/2,ind-stencil_width/2,-1)
        else:
            ind=argmin(Z_sym[imod])
            section=range(ind-stencil_width/2,ind+stencil_width/2)
        x=R_sym[imod,section]
        y=Z_sym[imod,section]
        y_int=interp1d(x,y,kind=interpol_order)
        x_fine=linspace(amin(x),amax(x),stencil_width*100)
        y_fine=y_int(x_fine)
        if o:
            x_at_extremum=x_fine[argmax(y_fine)]
            delta_upper[imod]=(R0[i]-x_at_extremum)/r_avg[i]
            Z_max=amax(y_fine)
        else:
            x_at_extremum=x_fine[argmin(y_fine)]
            delta_lower[imod]=(R0[i]-x_at_extremum)/r_avg[i]
            Z_min=amin(y_fine)
    #calculate kappa
    kappa[imod]=(Z_max-Z_min)/2./r_avg[i]
                
#linear extrapolation (in psi) for axis values
#delta_upper[0]=2*delta_upper[1]-delta_upper[2]
#delta_lower[0]=2*delta_lower[1]-delta_lower[2]
#kappa[0]=2*kappa[1]-kappa[2]
#zeta[0]=2*zeta[1]-zeta[2]
delta=0.5*(delta_upper+delta_lower)

#calculate zeta
for i in psi_stencil:
    imod=i-psi_stencil[0]
    x=arcsin(delta[imod])
    #find the points that correspond to Miller-theta=+-pi/4,+-3/4*pi and extract zeta from those
    for o in range(4):
        if o==0:
            val=pi/4.
            searchval=cos(val+x/sqrt(2))
            searcharr=(R_sym[imod]-R0[i])/r_avg[i]
        elif o==1:
            val=3.*pi/4
            searchval=cos(val+x/sqrt(2))
            searcharr=(R_sym[imod]-R0[i])/r_avg[i]
        elif o==2:
            val=-pi/4.
            searchval=cos(val-x/sqrt(2))
            searcharr=(R_sym[imod]-R0[i])/r_avg[i]
        elif o==3:
            val=-3.*pi/4
            searchval=cos(val-x/sqrt(2))
            searcharr=(R_sym[imod]-R0[i])/r_avg[i]
        if o in [0,1]:
            searcharr2=searcharr[ntheta/2:]
            ind=find(searchval,searcharr2)+ntheta/2
        else:
            searcharr2=searcharr[0:ntheta/2]
            ind=find(searchval,searcharr2)
#        print o,ind
        section=range(ind-stencil_width/2,ind+stencil_width/2)
        theta_sec=theta_arr[section]
        if o in [0,1]:
            theta_int=interp1d(-searcharr[section],theta_sec,kind=interpol_order)
            theta_of_interest=theta_int(-searchval)
        else:
            theta_int=interp1d(searcharr[section],theta_sec,kind=interpol_order)            
            theta_of_interest=theta_int(searchval)
        Z_sec=Z_sym[imod,section]
        Z_sec_int=interp1d(theta_sec,Z_sec,kind=interpol_order)
#        print searchval,val, theta_sec
        Z_val=Z_sec_int(theta_of_interest)
        zeta_arr[imod,o]=arcsin((Z_val-Z_avg[i])/kappa[imod]/r_avg[i]) 
    zeta_arr[imod,1]=pi-zeta_arr[imod,1]
    zeta_arr[imod,3]=-pi-zeta_arr[imod,3]
#    print zeta_arr[i]
    zeta[imod]=0.25*(pi+zeta_arr[imod,0]-zeta_arr[imod,1]-zeta_arr[imod,2]+zeta_arr[imod,3])



Bref_efit=abs(F[0]/R0[0])
Lref_efit=sqrt(2*abs(phi_fine[-1])/Bref_efit)

kappa_spl=US(r_avg[psi_stencil],kappa,k=interpol_order,s=1e-5)
delta_spl=US(r_avg[psi_stencil],delta,k=interpol_order,s=1e-5)
zeta_spl=US(r_avg[psi_stencil],zeta,k=interpol_order,s=1e-5)
amhd=empty(pw,dtype=float)
amhd_Miller=empty(pw,dtype=float)
Vprime=empty(pw,dtype=float)
dV_dr=empty(pw,dtype=float)
V=empty(pw,dtype=float)
V_manual=empty(pw,dtype=float)
r_FS=empty(pw,dtype=float)
for i in psi_stencil:
    imod=i-psi_stencil[0]
    Vprime[imod]=abs(sum(qpsi[i]*R_sym[imod]**2/F[i])*4*pi**2/ntheta)
    dV_dr[imod]=abs(sum(qpsi[i]*R_sym[imod]**2/F[i])*4*pi**2/ntheta)/ravg_spl.derivatives(linpsi[i])[1]
#    V[imod]=trapz(Vprime[:imod+1],linpsi[psi_stencil]) 
    r_FS[imod]=average(sqrt((R_sym[imod]-R0[i])**2+(Z_sym[imod]-Z_avg[i])**2),weights=qpsi[i]*R_sym[imod]**2/F[i])
    amhd[imod]=-qpsi[i]**2*R0[i]*pprime[i]*8*pi*1e-7/Bref_miller**2/ravg_spl.derivatives(linpsi[i])[1]
#    amhd_Miller[imod]=-2*Vprime[imod]/(2*pi)**2*(V[imod]/2/pi**2/R0[i])**0.5*4e-7*pi*pprime[i]
    dq_dr_avg[imod]=q_spl.derivatives(r_avg[i])[1]
    dq_dpsi[imod]=q_spl_psi.derivatives(linpsi[i])[1]
    drR[imod]=R0_spl.derivatives(r_avg[i])[1]
    s_kappa[imod]=kappa_spl.derivatives(r_avg[i])[1]*r_avg[i]/kappa[imod]
    s_delta[imod]=delta_spl.derivatives(r_avg[i])[1]*r_avg[i]/sqrt(1-delta[imod]**2)
    s_zeta[imod]=zeta_spl.derivatives(r_avg[i])[1]*r_avg[i]
amhd_spl=US(r_avg[psi_stencil],amhd,k=interpol_order,s=1e-5)
rFS_spl=US(r_avg[psi_stencil],r_FS,k=interpol_order,s=1e-5)
drR_spl=US(r_avg[psi_stencil],drR,k=interpol_order,s=1e-5)
Zavg_spl=US(r_avg,Z_avg,k=interpol_order,s=1e-5)

#plot(r_avg[psi_stencil],V)
#figure()
plot(r_avg[psi_stencil],kappa)
title('Elongation')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$\kappa$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)
figure()
plot(r_avg[psi_stencil],s_kappa)
title('Elongation (Derivative)')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$s_\kappa$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)
figure()
plot(r_avg[psi_stencil],delta)
title('Triangularity')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$\delta$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)
figure()
plot(r_avg[psi_stencil],s_delta)
title('Triangularity (Derivative)')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$s_\delta$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)
#figure()
#plot(r_avg[psi_stencil],delta_diff)
#xlabel(r'$r_{avg}$',fontsize=14)
#ylabel(r'$\delta_u-\delta_l$',fontsize=14)
figure()
plot(r_avg[psi_stencil],zeta)
title('Squareness')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$\zeta$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)
figure()
plot(r_avg[psi_stencil],s_zeta)
title('Squareness (Derivative)')
xlabel(r'$r_{avg}$',fontsize=14)
ylabel(r'$s_\zeta$',fontsize=14)
axvline(r,0,1,ls='--',color='k',lw=2)

#select a given flux surface
ind=poi_ind-psi_stencil[0]#find(r_a,r_avg/r_avg[-1])
Lref=R0_pos
print '\n\nShaping parameters for flux surface r=%9.5g, r/a=%9.5g:' %(r,r_a)
#print 'r_FS= %9.5g (flux-surface averaged radius)\n' %rFS_spl(r)
print 'Copy the following block into a GENE parameters file:\n'
print 'trpeps  = %9.5g' %(r/R0_pos)
print 'q0      = %9.5g' %q_spl(r)
print 'shat    = %9.5g !(defined as r/q*dq_dr)' %(r/q_spl(r)*q_spl.derivatives(r)[1])
#print 'shat=%9.5g (defined as (psi-psiax)/q*dq_dpsi)' %((psi-linpsi[0])/q_spl(r)*q_spl_psi.derivatives(psi)[1])
print 'amhd    = %9.5g' %amhd_spl(r)
#print 'amhd_Miller=%9.5g' %amhd_Miller[ind]
#print test[ind]
print 'drR     = %9.5g' %drR_spl(r)
print 'kappa   = %9.5g' %kappa_spl(r)
print 's_kappa = %9.5g' % (kappa_spl.derivatives(r)[1]*r/kappa_spl(r))
print 'delta   = %9.5g' %delta_spl(r) 
print 's_delta = %9.5g' %(delta_spl.derivatives(r)[1]*r/sqrt(1-delta_spl(r)**2)) 
print 'zeta    = %9.5g' %zeta_spl(r) 
print 's_zeta  = %9.5g' %(zeta_spl.derivatives(r)[1]*r)
print 'minor_r = %9.5g' %(1.0)
print 'major_R = %9.5g' %(R0_pos/r*r_a)
print '\nFor normalization to major radius, set instead:'
print 'minor_r = %9.5g' %(r/R0_pos/r_a)
print 'major_R = %9.5g' %(1.0)
print 'The same conversion factor must be considered for input frequencies and gradients.'
print '\nAdditional information:'
print 'Lref        = %9.5g !for Lref=a convention' %r_avg[-1]
print 'Lref        = %9.5g !for Lref=R0 convention' %R0_pos
print 'Bref        = %9.5g' %Bref_miller
#minor radius at average flux surface elevation (GYRO definition)
print '\na (avg elev)= %9.5g' %r_avg[-1]
print 'R0          = %9.5g' %R0_pos 
print 'Z0          = %9.5g' %Zavg_spl(r)
print 'Lref_efit   = %9.5g' %Lref_efit
print 'Bref_efit   = %9.5g' %Bref_efit
print 'B_unit(GYRO)= %9.5g' %(q_spl(r)/r/ravg_spl.derivatives(psi)[1])
#print 'Vprime=%9.5g; dV_dr=%9.5g' %(Vprime[ind],dV_dr[ind])
#print 'V=%9.5g' %V[ind]
#minor radius defined by 0.5(R_max-R_min), where R_max and R_min can have any elevation
#print 'a_maxmin    = %9.5g' %r_maxmin[-1]
print 'dpsi/dr     = %9.5g' %(1./ravg_spl.derivatives(psi)[1])
#print 'dr_maxmin/dr     = %9.5g' %(rmaxmin_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])
print 'drho_tor/dr = %9.5g' %(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])
print 'Gradient conversion omt(rho_tor) -> a/LT; factor = %9.5g' %(r_avg[-1]*(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])) 
print 'Gradient conversion omt(rho_tor) -> R/LT; factor = %9.5g' %(R0_pos*(rho_tor_spl.derivatives(psi)[1]/ravg_spl.derivatives(psi)[1])) 
figure()
plot(R_tm[ind],Z_tm[ind],'k-',lw=2,label='original')
plot(R_sym[ind],Z_sym[ind],'r-',lw=2,label='symmetrized')
#for v in range(-10,10,4):
#zeta_test=-0.011086#zeta[ind]
plot(R0_pos+r*cos(theta_arr+arcsin(delta_spl(r))*sin(theta_arr)),Zavg_spl(r)+kappa_spl(r)*r*sin(theta_arr+zeta_spl(r)*sin(2*theta_arr)),label='Miller')
title('Flux surface shapes')
xlabel('$R$/m',fontsize=14)
ylabel('$Z$/m',fontsize=14)
gca().set_aspect('equal')
legend(loc=10,prop={'size':10})

provide_conversions=0
if provide_conversions:
    f=open('rho_tor-r_a_conv','w')
    f.write('#r/a    rho_tor\n')
    for i in range(0,nw):
        f.write('%16.8e %16.8e\n' %(ravg_spl(linpsi[i])/r_avg[-1],rho_tor_spl(linpsi[i])))


show()
    



##Fourier series representation of flux surfaces    
#RZ_c=R+1j*Z
#for i in range(0,nw,4):
#    RZ_c_fft=fft(RZ_c[i])
#    #plot(RZ_c_fft)
#    trunc_order=2
#    RZ_c_fft[trunc_order:nw/2]=0.;RZ_c_fft[nw/2+1:-trunc_order+1]=0.
#    RZ_trunc=ifft(RZ_c_fft)
#    R_trunc=real(RZ_trunc)
#    Z_trunc=imag(RZ_trunc)
#    plot(R_trunc,Z_trunc,lw=2,ls='--')
#contour(rgrid,zgrid,psirz,20)
#gca().set_aspect('equal')
#show()
