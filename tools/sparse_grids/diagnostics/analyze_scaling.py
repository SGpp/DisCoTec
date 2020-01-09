#!/usr/bin/env python
import sys, numpy
from optparse import OptionParser

from gplot import *

try:
    from analyze_SG import *
except:
    print 'SG++ not found'
from analyze_DG import *
from analyze_RD import *

parser = OptionParser(usage='%prog [options] scanpath testpath') #,epilog='hierdertext')
parser.add_option("-a", "--allerrors", action="store_const", const=1, dest="mode", default=0,
                  help="Show errors for all sample points")
parser.add_option("-e", "--extremeerrors", action="store_const", const=2, dest="mode", default=0,
                  help="Show extreme errors (min/max)")
parser.add_option("-l", "--limit", dest="limit", default=0., type='float',
                  help="Limit test sample to points with abs(val) < limit")
(opts, args) = parser.parse_args()

testpath = args[-1] 

#read test points
if os.path.exists(testpath+'/RDlogs'):
    testpoints=analyze_RD(testpath)
    tpars=numpy.transpose(numpy.array(testpoints.serialgrid))
    tres=numpy.array(testpoints.serialvec)
else:
    tfile = open(testpath+'/quasilinear.log', 'r')
    head=tfile.readline()  #remove first line
    ndim=len(head.split('|'))-2
    tpars=[]
    tres=[]
    for line in tfile:
        tpars.append([])
        spline=line.split('|')
        for ind in range(ndim):
            tpars[int(spline[0])-1].append(float(spline[ind+1]))
        tres.append(float(spline[ndim+1]))
    tfile.close()
    tpars=numpy.transpose(numpy.array(tpars))
    tres=numpy.array(tres)

#errors are given in % relative to the range of the test values
scale=100/(max(tres)-min(tres))

if opts.limit>0:
    #limit the test points
    ltres=[]
    ltpars=[]
    for ind in range(tres.shape[0]):
        if abs(tres[ind])<opts.limit:
            ltres.append(tres[ind])
            ltpars.append(tpars[:,ind])
    tres=numpy.array(ltres)
    tpars=numpy.transpose(numpy.array(ltpars))



graphs=[]

for pind in range(len(args)-1):
    path=args[pind]
    if len(args)>2:
        ext='%d' % (pind+1)
    else:
        ext=''

    if os.path.exists(path+'/DGlogs'):
        print 'analyzing dense grid %s' % ext
        grid=analyze_DG(path)
    else:
        print 'analyzing sparse grid %s' % ext
        grid=analyze_SG(path)

    pnum=grid.points_per_stage
    looplim=grid.nstages
    ntest=len(tres)

    inderr=numpy.zeros((ntest,looplim))

    for reflev in range(looplim):
        #compute estimates for a given level
        estres=grid.ninterp(tpars,reflev)
        inderr[:,reflev]=abs(estres-tres)*scale

    averr=numpy.sum(inderr,axis=0)/ntest

    graphs.append(graph('w lp ti "average error %s"'%ext,pnum,averr))
    if opts.mode==1:
        for ind in range(len(tres)):
            graphs.append(graph('w p notitle',pnum,inderr[ind]))
    elif opts.mode==2:
        minv=[]
        maxv=[]
        for reflev in range(grid.nstages):
            maxval=0
            minval=1e40
            for tpoint in range(len(tres)):
                maxval=max(maxval,inderr[tpoint][reflev])
                minval=min(minval,inderr[tpoint][reflev])
            minv.append(minval)
            maxv.append(maxval)
        graphs.append(graph('w lp ti "maximum error %s"'%ext,pnum,maxv))
        graphs.append(graph('w lp ti "minimum error %s"'%ext,pnum,minv))

    plots=[plot(graphs,'number of points','relative error in %',loglog=True)]
    mplot=[multiplot(plots,pause=-1)]
show(mplot,'scaling')
