#!/usr/bin/env python
import sys, os, math
sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/release_090/lib/pysgpp')
sys.path.append('/afs/ipp-garching.mpg.de/home/f/flm/duennegitter/release_090/bin')
from pysgpp import *
import tools

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createModLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (gridStorage.dim())

# create regular grid, level 3
level = 2
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  %d" % (gridStorage.size())

# set function values in alpha
f = lambda x, y: math.sin( (x-0.5)*30 + (y-0.5)*5 )

maxadapt=15
for numadapt in range(maxadapt+1):

    # create coefficient vector
    alpha = DataVector(gridStorage.size())
    print "length of alpha-vector: %d" % (len(alpha))

    for i in xrange(gridStorage.size()):
        gp = gridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))

    # hierarchize
    grid.createOperationHierarchisation().doHierarchisation(alpha)

    tools.writeGnuplot("out%02d.gplt"%(numadapt), grid, alpha, 50)

    # refine
    if numadapt < maxadapt:
        gridGen.refine(SurplusRefinementFunctor(alpha, 1, 0.1))


#print grid.serialize()
#print alpha
