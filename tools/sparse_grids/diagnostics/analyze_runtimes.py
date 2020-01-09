#!/usr/bin/env python
import os
from optparse import OptionParser
from numpy import *

from gplot import *
try:
    from analyze_SG import *
except:
    print 'SG++ not found'
from analyze_DG import *
from analyze_RD import *

parser = OptionParser(usage='%prog [options] path')#,epilog='hierdertext')
parser.add_option("-i", "--iterations", action="store_const", const=1, dest="mode", default=0,
                  help="show iterations instead of runtimes")
parser.add_option("-c", "--cutoff", dest="ubound",type='float', default=0.,
                  help="cutoff for runtimes")
parser.add_option("-r", "--reuse",  dest="reuse", action="store_true",default=False,
                  help="analyze subspace recycling")
(opts, args) = parser.parse_args()
path = args[0]

if os.path.exists(path+'/DGlogs'):
    grid=analyze_DG(path)
elif os.path.exists(path+'/SGlogs'):
    grid=analyze_SG(path)
elif os.path.exists(path+'/RDlogs'):
    grid=analyze_RD(path)

pnum=grid.points_per_stage
ndim=grid.ndim

try:
    tfile=open('%s/parameters'% path,'r')
except:
    tfile=open('%s/genefiles/0k/parameters_1'% path,'r')

n_parallel_sims=1
for line in tfile:
    if line.find('n_parallel_sims')>=0:
        n_parallel_sims=int(line.split('=')[1].strip())
        break

command="grep 'Percentage' %s/genefiles/logs/gene_* > temp1" % path
out=os.system(command)
tfile=open('temp1','r')
ifrac=[]
for line in tfile:
    ifrac.append(float(line.split(':')[2].strip().strip('%')))
tfile.close()
command="grep 'Total' %s/genefiles/logs/gene_* > temp1" % path
out=os.system(command)
tfile=open('temp1','r')
rtime=[]
refnum=[]
for line in tfile:
    rtime.append(float(line.split(':')[2].strip().strip('sec')))
    sline=line.split(':')[0].split('_')
    if 'k' in sline[-1]:
        refnum.append(int(sline[-2]))
    else:
        refnum.append(int(sline[-1]))
tfile.close()

#combine runtimes and idletimes from the same refinement stage
runtime = zeros(grid.nstages)
idletime = zeros(grid.nstages)
while len(refnum)>0:
    ind=refnum.pop()
    rt=rtime.pop()
    runtime[ind]+=rt
    idletime[ind]+=ifrac.pop()/100*rt

totalidle=0.
totalgenetime=0.
for ind in range(len(runtime)):
    totalidle=totalidle+idletime[ind]
    totalgenetime=totalgenetime+runtime[ind]

logfile=open('%s/runtimes.log'%path)
logfile.readline()  #discard header
iterations=[]
times=[]
for line in logfile:
    iterations.append(int(line.split('|')[ndim+1]))
    times.append(float(line.split('|')[ndim+2]))
logfile.close()

if opts.reuse:
    distances=[]
    distfile=open('%s/distances.log'%path)
    distfile.readline()  #discard header
    for line in distfile:
        distances.append(sqrt(float(line.split('|')[ndim+1])))
    distfile.close()

maxit=max(iterations)
sumtime=0.
nctime=0.
for ind in range(len(times)):
    sumtime=sumtime+times[ind]
    if iterations[ind]==maxit:
        nctime=nctime+times[ind]
runnrs=range(1,len(times)+1)

try:
    stime=os.stat('%s/start_tag' % grid.path).st_ctime
except:
    stime=os.stat('%s/genefiles/0k/parameters_1' % path).st_ctime

#compute script times
oldtime=stime
reftimes=[]
for stage in range(grid.nstages):
    newtime=os.stat('%s/%s_%d'% (grid.path,grid.logname,stage)).st_ctime
    reftimes.append(newtime-oldtime)
    oldtime=newtime

scripttime=0
for entry in reftimes:
    scripttime+=entry

print "Total time for scan: %f s" % scripttime
print "GENE fraction: %f %%" % (totalgenetime/scripttime*100)
print "-----------------------"
print "Details for gene:"
print "Idle time due to load imbalancing at the end of each scan: %f %%" % (totalidle/totalgenetime*100)
tgruntime=(totalgenetime-totalidle)*n_parallel_sims
print "Average time for initialization: %f %%" % ((tgruntime-sumtime)/tgruntime*100)
print "Time for runs that stopped due to ev_max_it: %f %%" % (nctime/tgruntime*100)

if opts.mode==1:
    data=iterations
    what='iterations'
    whatlong=what
else:
    data=times
    what='time'
    whatlong='runtime [s]'
#compute averages
pnum.insert(0,0)  #insert zero as first index of the times list
avs=[]
evtime=[]
avdist=[]
for stage in range(grid.nstages):
    avs.append(0.)
    evtime.append(0.)
    for ind in range(pnum[stage],pnum[stage+1]):
        avs[stage]+=data[ind]
        evtime[stage]+=times[ind]
    avs[stage]=avs[stage]/(pnum[stage+1]-pnum[stage])
    if opts.reuse:
        avdist.append(0.)
        for ind in range(pnum[stage],pnum[stage+1]):
            avdist[stage]+=distances[ind]
        avdist[stage]=avdist[stage]/(pnum[stage+1]-pnum[stage])

#add 1 to have the axes consistent with numbering of the files
pnum=[entry+1 for entry in pnum]

nbins=min(len(times)/10,50)

if nbins<2:
    sys.exit('\nERROR: not enough runs for histogram!\n')


#histogram
hist, bins=histogram(data, bins=nbins)
hist=concatenate((array([0]),hist,array([0])))
firstbin=2*bins[0]-bins[1]
bins=concatenate((array([firstbin]),bins))
whist=hist*bins

if opts.ubound > 0.:
    ubound=opts.ubound
else:
    ubound=ceil(bins[-1])

genefrac=[]
ifrac=[]
gfrac=[]
evfrac=[]
inifrac=[]
sfrac=[]
for ind in range(grid.nstages):
    ev=evtime[ind]/n_parallel_sims/reftimes[ind]*100
    gene=runtime[ind]/reftimes[ind]*100
    idle=idletime[ind]/reftimes[ind]*100
    ini=gene-ev-idle
    inifrac.append(ini)
    evfrac.append(ini+ev)
    ifrac.append(ini+ev+idle)
    gfrac.append(gene)
    sfrac.append(100.0)
#have last value twice for steps plot
avs.append(avs[-1])
ifrac.append(ifrac[-1])
gfrac.append(gfrac[-1])
evfrac.append(evfrac[-1])
inifrac.append(inifrac[-1])
sfrac.append(sfrac[-1])
ref=range(len(gfrac))

plots=[]
#plot1
graphs=[]
graphs.append(graph("w p title '%s for each gene run'" % what,runnrs,data))
graphs.append(graph("w steps title 'average %s in each stage'" % what,pnum,avs))
plots.append(plot(graphs,'run number',whatlong,yrange=[0,ubound]))

if not opts.reuse:
#plot2
    graphs=[graph("w steps title 'histogram of %s'" % what,bins,hist)]
    plots.append(plot(graphs,whatlong,'number of computations',xrange=[0,ubound]))

#plot3
    graphs=[graph("w steps title 'weighted histogram of %s'" % what,bins,whist)]
    plots.append(plot(graphs,whatlong,'accumulated %s in each bin' % what,xrange=[0,ubound]))
else:
#plot2
    graphs=[graph("w lp title 'average runtime per stage'",ref[:-1],avs[:-1])]
    plots.append(plot(graphs,'refinement stage','average runtime [s]'))
#plot3
    if dense:
        startind=2
    else:
        startind=1
    graphs=[graph("w lp title 'average distance per stage'",ref[startind:-1],avdist[startind:])]
    plots.append(plot(graphs,'refinement stage','average distance'))


#plot4
graphs=[]
graphs.append(graph("w steps title 'gene initialization'",ref,inifrac))
graphs.append(graph("w steps title 'gene EV solver'",ref,evfrac))
graphs.append(graph("w steps title 'gene idle time'",ref,ifrac))
graphs.append(graph("w steps title 'rest: python'",ref,sfrac))

plots.append(plot(graphs,'refinement stage','accumulated fraction of runtime [%]'))

mplots=[multiplot(plots,pause=-1)]
show(mplots,'runtimes')
