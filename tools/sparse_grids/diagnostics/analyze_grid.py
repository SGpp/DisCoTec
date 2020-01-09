#!/usr/bin/env python
import sys, numpy
from optparse import OptionParser

sys.path.append('../')
from analyze_scanlog import *
from gplot import *
try:
    from analyze_SG import *
except:
    print 'SG++ not found'
from analyze_DG import *
from analyze_RD import *

parser = OptionParser(usage='%prog [options] scanpath [testpath]') #,epilog='hierdertext')
parser.add_option("-m", "--movie", action="store_true", dest="movie", default=False,
                  help="Show movie of grid refinement")
parser.add_option("-l", "--limit", dest="tlim", default=0,
                  help="Limit for number of test points to be shown (points with the smallest errors are discarded)")
parser.add_option("-s", "--select", dest="plotnr", default=-1,
                  help="Show only 2D projection with the selected number")
(opts, args) = parser.parse_args()

path = args[0] 

testdata=False
if len(args)>1:
    if os.path.exists(args[1]):
        testpath = args[1]
        testdata=True
        if os.path.exists(testpath+'/RDlogs'):
            testpoints=analyze_RD(testpath)
            tpars=numpy.transpose(numpy.array(testpoints.serialgrid))
            tres=numpy.array(testpoints.serialvec)
        else:
            test=analyze_scanlog(testpath)
            tpars=numpy.array(test.coords)

if os.path.exists(path+'/DGlogs'):
    print 'analyzing dense grid'
    grid=analyze_DG(path)
else:
    print 'analyzing sparse grid'
    grid=analyze_SG(path)

refpoints=grid.points_per_stage
scan=analyze_scanlog(path)

if opts.movie:
    startlev=0
else:
    startlev=len(refpoints)-1

mplots=[]
for reflev in range(startlev,len(refpoints)):
    if testdata:
        form=''
        for ind in range(scan.scandim+1):
            form+='f,'
        errlist=numpy.zeros(len(tres),dtype=[('coords','f4',(scan.scandim)),('error','f4')])
        estres=grid.ninterp(tpars,reflev)
        errlist['coords']=numpy.transpose(tpars)
        errlist['error']=abs(estres-tres)

        slist=numpy.sort(errlist,order='error')
        errlist=slist[-int(opts.tlim):]
    plots=[]
    pnum=0
    for ind1 in range(scan.scandim):
        for ind2 in range(ind1+1,scan.scandim):
            pnum+=1
            if pnum==int(opts.plotnr) or int(opts.plotnr)<0:
                graphs=[]
                graphs.append(graph('ls 1 ti "grid points"',scan.coords[ind1][0:refpoints[reflev]],scan.coords[ind2][0:refpoints[reflev]]))
                if testdata:
                    graphs.append(graph('ti "test points" w p pt 2 palette z',errlist['coords'][:,ind1],errlist['coords'][:,ind2],errlist['error']))
                plots.append(plot(graphs,scan.scanpars[ind1],scan.scanpars[ind2],xrange=grid.scanrange[ind1],yrange=grid.scanrange[ind2]))
    if reflev==len(refpoints)-1:
        mplots.append(multiplot(plots,pause=-1))
    else:
        mplots.append(multiplot(plots,pause=0.5))

#print plotdata
show(mplots,'2D-projections')
