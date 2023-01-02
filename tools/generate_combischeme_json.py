#!/usr/bin/env python3

from tokenize import group
import numpy as np
import itertools as it
from icecream import ic
from scipy.special import binom
import json


class ClassicDiagonalActiveSet():
    def __init__(self, lmax, lmin=None, diagonalIndex=0):
        self.lmax = lmax
        if lmin == None:
            self.lmin = [1 for i in range(len(self.lmax))]
        else:
            self.lmin = lmin
        self.diagonalIndex = diagonalIndex

    def getMinLevelSum(self):
        return self.getLevelMinima().sum()

    def getLevelMinima(self):
        maxInd = np.argmax(np.array(self.lmax)-np.array(self.lmin))
        lm = np.array(self.lmin)
        lm[maxInd] = self.lmax[maxInd]
        return lm

    def getActiveSet(self):
        listOfRanges = [list(range(0, self.lmax[i]+1))
                        for i in range(len(self.lmax))]
        listOfAllGrids = list(it.product(*listOfRanges))
        s = set()
        levelSum = self.getMinLevelSum()+self.diagonalIndex
        for grid in listOfAllGrids:
            if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(self.lmin)).all()):
                s.add(grid)

        return s


class combinationSchemeArbitrary():
    def __init__(self, activeSet):
        self.activeSet = activeSet
        self.updateDict()

    def getDownSet(self, l):
        subs = [list(range(0, x + 1)) for x in l]
        downSet = it.product(*subs)
        return downSet

    def getSubspaces(self):
        subspacesSet = set()
        for l in self.dictOfScheme:
            for subspace in self.getDownSet(l):
                subspacesSet.add(subspace)
        return subspacesSet

    def updateDict(self):
        lmin = self.activeSet.lmin
        lmax = self.activeSet.lmax
        firstLevelDifference = lmax[0] - lmin[0]
        dim = len(lmin)
        uniformLevelDifference = [(lmax[i] - lmin[i]) == firstLevelDifference for i in range(dim)]
        if uniformLevelDifference and firstLevelDifference >= dim:
            ic("binomial")
            # can calculate it by binomial formula
            listOfRanges = [list(range(lmin[i], lmax[i]+1))
                            for i in range(len(lmax))]
            self.dictOfScheme = {}
            for q in range(dim):
                coeff = (-1)**q * binom(dim-1, q)
                levelSum = sum(lmin) + firstLevelDifference - q
                # ic(coeff, levelSum)
                for grid in it.product(*listOfRanges):
                    if (np.sum(grid) == levelSum and (np.array(grid) >= np.array(lmin)).all()):
                        self.dictOfScheme[grid] = coeff

        else:
            self.dictOfScheme = {}
            for l in self.activeSet.getActiveSet():
                self.dictOfScheme[l] = 1

            dictOfSubspaces = {}

            for l in self.dictOfScheme:
                for subspace in self.getDownSet(l):
                    if subspace in dictOfSubspaces:
                        dictOfSubspaces[subspace] += 1
                    else:
                        dictOfSubspaces[subspace] = 1

            # # remove subspaces which are too much
            while(set(dictOfSubspaces.values()) != set([1])):

                for subspace in dictOfSubspaces:
                    currentCount = dictOfSubspaces[subspace]
                    if currentCount != 1:
                        diff = currentCount - 1

                        if subspace in self.dictOfScheme:

                            self.dictOfScheme[subspace] -= diff
                            if self.dictOfScheme[subspace] == 0:
                                del self.dictOfScheme[subspace]

                        else:
                            self.dictOfScheme[subspace] = -diff

                        for l in self.getDownSet(subspace):
                            dictOfSubspaces[l] -= diff

    def getCombinationDictionary(self):
        return self.dictOfScheme


lmin = [2]*6
lmax = [3]*6
activeSet = ClassicDiagonalActiveSet(lmax, lmin)
# scheme = active set und Koeffizienten dazu
scheme = combinationSchemeArbitrary(activeSet)

numGridsOfSize = {}

totalNumPointsCombi = 0
for key, value in scheme.getCombinationDictionary().items():
    totalNumPointsCombi += np.prod([2**l + 1 for l in key])
    levelSum=np.sum([l for l in key])
    if levelSum in numGridsOfSize:
        numGridsOfSize[levelSum] += 1
    else:
        numGridsOfSize[levelSum] = 1

numGrids = len(scheme.getCombinationDictionary())
ic(numGrids)
ic(numGridsOfSize)
ic(totalNumPointsCombi, totalNumPointsCombi/1e13)
print(totalNumPointsCombi)
mem = (totalNumPointsCombi*8) # minimum memory requirement of full grids in scheme in bytes
mem = mem/(2**20) # memory requirement in Mebibytes

totalNumPointsSparse = 0
# use optimized sparse grid, with lmax 1 lower in each dim
reduced_lmax = [max(lmax[d] - 1, lmin[d]) for d in range(len(lmax))]
# reduced_lmax = [9, 5, 8,8,7]
activeSetSparse = ClassicDiagonalActiveSet(reduced_lmax, lmin)
schemeSparse = combinationSchemeArbitrary(activeSetSparse)
for key in schemeSparse.getSubspaces():
    # ic(key)
    totalNumPointsSparse += np.prod([2**(l-1) if l > 0 else 2 for l in key])
print(totalNumPointsSparse)
mem_sg = (totalNumPointsSparse*8)
mem_sg = mem_sg/(2**20) # memory requirement of sparse grid in Mebibytes

mem_total = mem_sg+mem
ic(mem_total, mem_sg/mem_total)

numProcessGroups = 3
assignment = { i : {} for i in range(numProcessGroups)}
roundRobinIndex = 0
for size in numGridsOfSize:
    # ic(size)
    for key, value in scheme.getCombinationDictionary().items():
        levelSum=np.sum([l for l in key])
        if levelSum == size:
            # ic(key, value)
            assignment[roundRobinIndex][key] = value
            roundRobinIndex = (roundRobinIndex+1)%numProcessGroups

schemeList = []
for group_no in assignment.keys():
    # ic(assignment[group_no])
    schemeList += [{"coeff": coeff, "level": list(level), "group_no": group_no}
              for level, coeff in assignment[group_no].items()]

# ic(schemeList)
jsonString = json.dumps(schemeList)

with open('scheme.json', 'w') as f:
    f.write(jsonString)
