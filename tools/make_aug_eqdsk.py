#!/usr/bin/env python
#author: D. Told (dtold@ipp.mpg.de)
from ctypes import *
from pylab import *
from sys import exit
from scipy.interpolate import UnivariateSpline
import sys

shotnr =   int(raw_input("Enter shot number: "))
timept = float(raw_input("Enter time in sec: "))
print ''
print 'Hit Enter for defaults given in brackets'
diag = raw_input("Select diag [eqi]: ")
if not diag:
	diag='eqi'
exp =  raw_input("Select exp [augd]: ")
if not exp:
	exp='augd'
edstr=raw_input("Enter edition [0]: ")
if not edstr:
	ed = 0
else:
	ed = int(edstr)

#maximum dimension of AUG magnetic equilibria
rmin = 0.75
rmax = 2.67
zmin = -1.504
zmax = 1.504

#link to libkk
try:
	if sizeof(c_long) == 8:
	    libkkso = '/afs/ipp/aug/ads/lib64/@sys/libkk.so'
	else:
	    libkkso = '/afs/ipp/aug/ads/lib/@sys/libkk.so'

	libkk = cdll.LoadLibrary(libkkso)
except:
	libkk=0
	pass	
if not libkk:
	print 'Could not find libkk for access to ASDEX equilibria.'
	exit()

#cfarr returns a C-type float array with n elements
def cfarr(n):
	return c_float*n 

#convert discharge data to C-types
shot=c_int(shotnr)
edition=c_int(ed)
time=c_float(timept)

#standard AUG resolution
rdim=129
zdim=257

#default values for limiter and boundary points to fulfil the official efit format
nbbs=6
limtr=6

#definitions for libkk calls
m=c_int(rdim-1)
n=c_int(zdim-1)
mdim=c_int(rdim)
error=c_int(0)
lpfx=c_int(10)
pfxx=cfarr(10)(0.)
rpfx=cfarr(10)(0.)
zpfx=cfarr(10)(0.)
Ri=cfarr(rdim)(0.)
zj=cfarr(zdim)(0.)

pfm=empty((zdim,rdim),dtype=float)
fpf=(c_float*1)()
rhopol=(c_float*1)()
lin=c_int(1)

#cylindrical coordinate grid
Ri[:]=array([rmin+i*(rmax-rmin)/(rdim-1) for i in range(rdim)])
zj[:]=array([zmin+i*(zmax-zmin)/(zdim-1) for i in range(zdim)])

#provide arrays with 2000 entries for 1d profiles
#(which should hopefully be sufficient, but can also easily be changed)
lpf_p=2000
lpf=c_int(lpf_p)
pfl=(c_float*lpf_p)()
qpl=(c_float*lpf_p)()
pres=(c_float*lpf_p)()
presp=(c_float*lpf_p)()
f=(c_float*lpf_p)()
fp=(c_float*lpf_p)()
ffp=(c_float*lpf_p)()
#typ=11: calculate 1d profiles if not present in the shotfile
typ=c_int(11)

#get position of magnetic axis, psi_sep and psi_ax
libkk.kkeqpfx(byref(error), exp, diag, shot, byref(edition), byref(time),
	      byref(lpfx), byref(pfxx), byref(rpfx), byref(zpfx))
rmag=rpfx[0]
zmag=zpfx[0]
psi_sep=pfxx[1]
psi_ax=pfxx[0]

#read 2d poloidal flux matrix
for i in range(rdim):
	for j in range(zdim):
		rin=(c_float*1)(Ri[i])
		zin=(c_float*1)(zj[j])
		libkk.kkrzpfn(byref(error), exp, diag, shot, byref(edition), byref(time),
			      rin, zin, lin,
			      byref(fpf), byref(rhopol))
		pfm[j,i]=fpf[0]
#get 1d profiles
libkk.kkeqqpl(byref(error), exp, diag, shot, byref(edition), byref(time),
	      byref(lpf), byref(pfl), byref(qpl))
libkk.kkeqpres(byref(error), exp, diag, shot, byref(edition), byref(time),
	       byref(lpf), byref(pfl), byref(pres), byref(presp))
libkk.kkeqffs(byref(error), exp, diag, shot, byref(edition), byref(time),
	      typ, byref(lpf), byref(pfl), byref(f), byref(fp))
libkk.kkeqffp(byref(error), exp, diag, shot, byref(edition), byref(time),
	      typ, byref(lpf), byref(pfl), byref(ffp))
#get magnetic field at axis
rin=(c_float*1)(rmag)
zin=(c_float*1)(zmag)
lin=c_int(1)
br=(c_float*1)()
bz=(c_float*1)()
bt=(c_float*1)()
fpf=(c_float*1)()
jpol=(c_float*1)()
libkk.kkrzbrzt(byref(error), exp, diag, shot, byref(edition), byref(time),
	       rin, zin, lin,
	       byref(br), byref(bz), byref(bt), byref(fpf), byref(jpol))
Btor=bt[0]

#generate equidistant psi-grid with rdim points
psi_lin=linspace(psi_ax,psi_sep,rdim)
pfl_np=array(pfl[:])
#generate equidistant rho_pol/rho_tor-grid with rdim points
rho_pol=sqrt((psi_lin-psi_ax)/(psi_sep-psi_ax))
rho_tor=empty(rdim)
for i in range(rdim):
	rhopf = (c_float * 1)(rho_pol[i])
	lrho = c_int(1)
	rhotf = (c_float * 1)()
	ftf = (c_float * 1)()
	libkk.kkrhopto(byref(error), exp, diag, shot, byref(edition), byref(time),
		       rhopf, lrho,
		       byref(rhotf), byref(fpf), byref(ftf))
	rho_tor[i]=rhotf[0]
#lpf.value contains the number of actual entries in the 1d profiles
lpf_p=lpf.value+1
#cut arrays to this size
pfl_np=array(pfl_np[:lpf_p])
f=array(f[:lpf_p])
ffp=array(ffp[:lpf_p])
pres=array(pres[:lpf_p])
presp=array(presp[:lpf_p])
qpl=array(qpl[:lpf_p])

#if psi decreases towards the edge, reverse the arrays (requirement
#for standard numpy linear interpolation)
if pfl_np[0]>pfl_np[-1]:
	pfl_np=pfl_np[::-1]
	f=f[::-1]
	ffp=ffp[::-1]
	pres=pres[::-1]
	presp=presp[::-1]
	qpl=qpl[::-1]

#Spline interpolation from the lpf-sized grid to an rdim-sized grid
spl_f=UnivariateSpline(pfl_np,f,s=0)
f_out=spl_f(psi_lin)
spl_ffp=UnivariateSpline(pfl_np,ffp,s=0)
ffp_out=spl_ffp(psi_lin)*2.*pi
spl_p=UnivariateSpline(pfl_np,pres,s=0)
p_out=spl_p(psi_lin)
spl_pp=UnivariateSpline(pfl_np,presp,s=0)
pp_out=spl_pp(psi_lin)*2.*pi
#neglect safety factor at separatrix as this seems to be broken
#in the kk-routines
spl_q = UnivariateSpline(pfl_np[1:lpf_p],qpl[1:lpf_p], s=0)
q_out = spl_q(psi_lin)

#some plotting commands; "show()" will show the plot in a popup window,
#which allows also to save it.
#plot(pfl_np[1:lpf_p],qpl[1:lpf_p])
#plot(psi_lin[:-1],q_out[:-1])
subplot(121,aspect='equal')
contour(Ri,zj,pfm,50)
xlabel(r'R/m')
ylabel(r'Z/m')

subplot(122)
plot(rho_pol,abs(q_out),'b.')
grid(True)
xlabel('rho_pol')
ylabel('q')
ylim(0)

ttitle = "Shot "+str(shotnr)+", diag '"+diag+"', exp '"+exp+"', ed '"+str(ed)+"' at t="+ str(timept)+"s"
suptitle(ttitle)

print 'Check and close plot window to proceed ...'
show()

#create G-EQDSK file
filename='AUG_%d_%.4gs.eqd'%(shot.value,time.value)
f=open(filename,'w')
print 'Writing G-EQDSK file "%s".' %(filename) 
f.write('%8s%8d%7.4gs%8s%8s%8s%4d%4d%4d\n' %(str.upper(exp),shot.value,time.value,str.upper(diag),"","",ed,rdim,zdim))
f.write('%16.9e%16.9e%16.9e%16.9e%16.9e\n' %(rmax-rmin,zmax-zmin,rmag,rmin,(zmax+zmin)/2.))
f.write('%16.9e%16.9e%16.9e%16.9e%16.9e\n' %(rmag,zmag,psi_ax/2./pi,psi_sep/2./pi,Btor))
f.write('%16.9e%16.9e%16.9e%16.9e%16.9e\n' %(0.,psi_ax/2./pi,0.,rmag,0.))
f.write('%16.9e%16.9e%16.9e%16.9e%16.9e\n' %(zmag,0.,psi_sep/2./pi,0.,0.))
for i in range(rdim):
	f.write('%16.9e' %(f_out[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
for i in range(rdim):
	f.write('%16.9e' %(p_out[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0:f.write('\n')
for i in range(rdim):
	f.write('%16.9e' %(ffp_out[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
for i in range(rdim):
	f.write('%16.9e' %(pp_out[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
for j in range(zdim):
	for i in range(rdim):
		f.write('%16.9e' %(pfm[j,i]/2./pi))
		if (i+j*rdim+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
#make q positive
for i in range(rdim):
	f.write('%16.9e' %(abs(q_out[i])))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
f.write('%5d%5d\n' %(nbbs,limtr))
#no boundary/limiter data at the moment -> filled with zeros
for i in range(2*nbbs):
	f.write('%16.9e' %(0.))
	if (i+1)%5==0: f.write('\n')
if ((i+1)%5!=0): f.write('\n')
for i in range(2*limtr):
	f.write('%16.9e' %(0.))
	if (i+1)%5==0: f.write('\n')
if ((i+1)%5!=0): f.write('\n')
#additional: rho_pol and rho_tor
for i in range(rdim):
	f.write('%16.9e' %(rho_pol[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
for i in range(rdim):
	f.write('%16.9e' %(rho_tor[i]))
	if (i+1)%5==0: f.write('\n')
if (i+1)%5!=0: f.write('\n')
f.close()



