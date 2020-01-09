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

parser = OptionParser(usage='%prog [options] scanpath') #,epilog='hierdertext')
parser.add_option("-c", "--center", action="store", dest="center", default=[],
                  help="Set center for slices/lines in all dimensions")
parser.add_option("-m", "--manyplots", action="store_false", dest="multiplot", default=True,
                  help="Switch off multiplot of all slices/lines")
parser.add_option("-l", "--lines", action="store_true", dest="lines", default=False,
                  help="Plot lines instead of 2D slices")
parser.add_option("-r", "--resolution", dest="res", default=0, type='int',
                  help="Resolution for cuts")
(opts, args) = parser.parse_args()

path = args[0] 

#read grid and result vector
if os.path.exists(path+'/DGlogs'):
    print 'analyzing dense grid'
    grid=analyze_DG(path)
else:
    print 'analyzing sparse grid'
    grid=analyze_SG(path)

scan=analyze_scanlog(path)

#choose point for the slices
center=numpy.zeros((grid.ndim,1))

if opts.center==[]:
    for ind in range(grid.ndim):
        center[ind,0]=(grid.scanrange[ind][0]+grid.scanrange[ind][1])/2.
else:
    numbers=opts.center.strip('[]').split(',')
    if len(numbers) != grid.ndim:
        sys.exit('\nERROR: The center option format must be [c1,c2,..,cn],\nit must have %d entries to match the chosen scan directory\n'% grid.ndim)
    for ind in range(len(numbers)):
        center[ind,0]=float(numbers[ind])

print 'showing all cuts through',center[:,0]

def getcut(resolution,x0):
    plotdim=[]
    for ind in range(len(x0)):
        if x0[ind]<0:
            plotdim.append(ind)

    if len(plotdim) ==0:
        #point
        pval1=grid.ninterp(x0,grid.nstages-1)
        pval2=[]
    elif len(plotdim) == 1:
        pval1=numpy.zeros(resolution)
        pval2=numpy.zeros(resolution)
        #line
        sys.exit()
        x=copy.deepcopy(x0)
        for i in xrange(resolution):
            d=plotdim[0]
            x[d] = grid.scanrange[d][0] + float(i)/(resolution-1)*grid.scanrange[d][2]
            pval1[i]=x[d]               
            pval2[i]=grid.ninterp(x,grid.nstages-1)[0]
        pval[3]=[]
    elif len(plotdim) == 2:
        #surface
        pval1=numpy.zeros((resolution,resolution))
        pval2=numpy.zeros((resolution,resolution))
        pval3=numpy.zeros((resolution,resolution))
        x=copy.deepcopy(x0)
        for i in xrange(resolution):
            for j in xrange(resolution):
                d=plotdim[0]
                x[d] = grid.scanrange[d][0] + float(i)/(resolution-1)*grid.scanrange[d][2]
                pval1[i,j]=x[d]
                    
                d=plotdim[1]
                x[d] = grid.scanrange[d][0] + float(j)/(resolution-1)*grid.scanrange[d][2]
                pval2[i,j]=x[d]
                    
                z=grid.ninterp(x,grid.nstages-1)

                pval3[i,j]=grid.ninterp(x,grid.nstages-1)[0]

    else:
        sys.stderr.write("Error! Can't plot grid with dimensionality %d..." % len(plotdim))

    if pval2==[]:
        return pval1
    elif pval3==[]:
        return pval1, pval2
    else:
        return pval1, pval2, pval3

if scan.scandim==1:
    opts.lines=True

#default for resolution
if opts.res==0:
    if opts.lines:
        opts.res=50
    else:
        opts.res=25

cval=getcut(opts.res,center)
plots=[]
for ind1 in range(scan.scandim):
    if opts.lines:
        tempcenter=copy.deepcopy(center)
        tempcenter[ind1]=-1
        x,y=getcut(opts.res,tempcenter)
        if opts.multiplot:
            graphs=[graph("w l ti 'interpolant'",x,y)]
            graphs.append(graph("ls 3 ti 'center'",center[ind1],cval))
            plots.append(plot(graphs,grid.scanpars[ind1],'interpolated value'))
        else:
            graphs=[graph("w l ti 'interpolant'",x,y)]
            graphs.append(graph("ls 3 ti 'center'",center[ind1],cval))
            plot=[plot(graphs,grid.scanpars[ind1],'interpolated value')]
            mplot=[multiplot(plot,pause=-1)]
            show(mplot,'cuts')
    else:
        for ind2 in range(ind1+1,scan.scandim):
            tempcenter=copy.deepcopy(center)
            tempcenter[ind1]=-1
            tempcenter[ind2]=-1
            x,y,z=getcut(opts.res,tempcenter)
            if opts.multiplot:
                graphs=[graph("notitle",x,y,z)]
                plots.append(plot(graphs,grid.scanpars[ind1],grid.scanpars[ind2],'interpolated value', contonly=True))
            else:
                graphs=[graph("notitle",x,y,z)]
                plots=[plot(graphs,grid.scanpars[ind1],grid.scanpars[ind2],'interpolated value')]
                mplot=[multiplot(plots,pause=-1)]
                show(mplot,'cuts')

if opts.multiplot:
    mplots=[multiplot(plots,pause=-1)]
    show(mplots,'cuts')


