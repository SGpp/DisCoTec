#!/usr/bin/env python
import sys, os
from optparse import OptionParser
sys.path.append('./sparse_grids/diagnostics')
from gplot import *

parser = OptionParser(usage='%prog [options] filename')
#parser.add_option("-a", "--allfiles", action="store_const", const=1, dest="all", default=0,
#                  help="Use all files which start with the filename given as argument")
(opts, args) = parser.parse_args()

filelist=[]
for arg in args:
    if os.path.isfile(arg):
        filelist.append(arg)
    else:
        dir=os.path.dirname(arg)
        nlist=os.listdir(dir)
        for filename in nlist:
            if 'autopar' in filename:
                filelist.append(dir+'/'+filename)

plist=[]
tlist=[]

#read test points
for file in filelist:
    print 'analyzing %s'%file
    pfile = open(file)
    for line in pfile:
        if ('best' in line) or ('Choice' in line):
            active=0
        elif 'parallelization' in line:
            active=1
            besttime=0.
        elif active:
            sline=line.split()
            time=float(sline[14])
            if besttime==0.:
                besttime=time
                bestperf=sline[0:9]
            else:
                perfchange=(besttime-time)/besttime*100.
                for i in range(9):
                    if sline[i] != bestperf[i]:
                        changedindex=i
                perfentry=[changedindex+1,int(sline[changedindex])]
                try: 
                    ind=plist.index(perfentry)
                    tlist[ind].append(perfchange)
                except:
                    plist.append(perfentry)
                    tlist.append([])
                    tlist[-1].append(perfchange)
                if time<besttime:
                    besttime=time
                    bestperf=sline[0:9]

    pfile.close()

plots=[]
graphs=[]
for ind in range(len(plist)):
    graphs=[graph("w p title 'perfvec(%i)=%i'" % (plist[ind][0],plist[ind][1]),range(len(tlist[ind])),tlist[ind])]
    plots.append(plot(graphs,'number','relative speedup in percent'))
mplots=[multiplot(plots,pause=-1)]
show(mplots,'perf')

