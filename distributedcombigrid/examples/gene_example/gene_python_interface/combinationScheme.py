'''
Created on Jan 3, 2014

@author: ckow
'''
import abc,logging
#from reportlab.graphics.shapes import NotImplementedError
import numpy as np
import itertools as it
from collections import OrderedDict

class combinationSchemeBase(object):
    '''
    base class for combination scheme
    '''
    
 
    
    def __init__(self):
        self.log=logging.getLogger(__name__)
        
    def getMaxTupel(self):
        d = self.getCombinationDictionary()
        maxLevel = None
        for k in list(d.keys()):
            if maxLevel==None:
                maxLevel=list(k)
            else:
                for i in range(self.getDim()):
                    if k[i]>maxLevel[i]:
                        maxLevel[i]=k[i]
        return tuple(maxLevel)

    @abc.abstractmethod
    def getCombinationDictionary(self):
        return

    @abc.abstractmethod
    def getDim(self):
        return
    

class combinationSchemeArbitrary(combinationSchemeBase):
    def __init__(self,activeSet):
        combinationSchemeBase.__init__(self)
        self.activeSet = activeSet
        self.updateDict()
        self.log.info('created combischeme with \n'+str(self.getCombinationDictionary()).replace(', (','\n('))
        
    def getDownSet(self,l):
        subs = [list(range(0, x + 1)) for x in l]
        downSet = it.product(*subs)
        return downSet
    
    def getUnifiedDownSet(self):
        s= set([])
        for a in self.activeSet:
            # print a,(self.getBoundedDownSet(a))
            s=s.union((self.getBoundedDownSet(a)))
            self.log.info(s)
        self.log.info(s)
        return s
    
    def getBoundedDownSet(self,l):
        lmin=self.getLmin()
        subs = list(map(lambda x,y:list(range(y, x + 1)),l,lmin))
        downSet = it.product(*subs)
#         self.log.info(list(downSet))
        return list(downSet)
        
    def getLmin(self):
        erg=None
        for s in self.activeSet:
            if erg==None:
                erg=list(s)
            else:
                for i in range(len(s)):
                    if s[i]<erg[i]:
                        erg[i]=s[i]
        return tuple(erg)
    
    def updateDict(self):
        self.dictOfScheme={}
        for l in self.activeSet:
            self.dictOfScheme[l] = 1
            
        dictOfSubspaces={}
        
        for l in self.dictOfScheme:
            for subspace in self.getDownSet(l):
                if subspace in dictOfSubspaces:
                    dictOfSubspaces[subspace] += 1
                else:
                    dictOfSubspaces[subspace] = 1
                    
        # self.log.debug(str(dictOfSubspaces))

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

        self.log.debug('created scheme '+str(self.dictOfScheme))
        
    def getCombinationDictionary(self):
        return self.dictOfScheme
    
    def getDim(self):
        return len(list(self.activeSet)[0])
    
    def removeSubspace(self,subspace):
        if subspace in self.dictOfScheme:
            del self.dictOfScheme[subspace]
            self.log.debug('subspace '+str(subspace)+' removed')
        else:
            raise AssertionError('The subspace does not exist in the scheme.')
    
    def removeZeroSubspaces(self):
        combiDict = self.getCombinationDictionary()
        for key, val in list(combiDict.items()):
            if val == 0:
                del self.dictOfScheme[key]
    
    def getNearestNeighbors(self):
        dictOfNeighbors = {}
        keys = list(self.dictOfScheme.keys())
        keys_search = list(self.dictOfScheme.keys())
        while keys:
            key = keys.pop()
            d = np.inf
            for key_n in keys_search:
                d_new = sum(abs(np.array(key) - np.array(key_n)))
                if d_new < d and d_new > 0:
                    d = d_new
                    key_nearest = key_n
                    dictOfNeighbors[key] = key_n
            if key_nearest in keys:
                keys.remove(key_nearest)
        return dictOfNeighbors

    def getKNearestNeighbors(self,K):
        if K == 1:
            return self.getNearestNeighbors()
        dictOfNeighbors = {}
        keys = list(self.dictOfScheme.keys())
        while keys:
            key = keys.pop()
            keys_search = list(self.dictOfScheme.keys())
            for k in range(K):
                d = np.inf
                for key_n in keys_search:
                    d_new = sum(abs(np.array(key) - np.array(key_n)))
                    if d_new < d and d_new > 0:
                        d = d_new
                        key_nearest = key_n
#                if key_nearest in keys:
#                    keys.remove(key_nearest)
                if k == 0:
                    dictOfNeighbors[key] = [key_nearest]
                else:
                    dictOfNeighbors[key].append(key_nearest)
                keys_search.remove(key_nearest)

        return dictOfNeighbors

    def hierSubs(self,lmin,lmax):
        list_ranges = [list(range(lmin[0],lmax[0]+1)),list(range(lmin[1],lmax[1]+1))]
        return it.product(*list_ranges)

    def keysPerLevel(self,lmin,lmax):
        dictOfKeysPerLevel = {}
        for key in self.hierSubs(lmin,lmax):
            l1 = sum(key)
            if l1 not in dictOfKeysPerLevel:
                dictOfKeysPerLevel[l1] = []
            dictOfKeysPerLevel[l1].append(key)
        return dictOfKeysPerLevel

    def getLargerSubspaces(self,grid,unified=True):
        if unified:
            downSet = self.getUnifiedDownSet()
        else:
            downSet = list(self.dictOfScheme.keys())
        largerSubspaces = []
        dim = len(grid)
        for level in downSet:
            isLarger = True
            for i in range(dim):
                if level[i] < grid[i]:
                    isLarger = False
                    break
            if isLarger:
                largerSubspaces.append(level)
        return largerSubspaces

class combinationSchemeFaultTolerant(combinationSchemeArbitrary):
    def __init__(self,activeSet):
        combinationSchemeArbitrary.__init__(self,activeSet.getActiveSet())
        self.activeSet = activeSet.getActiveSet()
        self.updateDict()
        self.numExtraLayers = 2
        self.extendSchemeExtraLayers(activeSet, self.numExtraLayers)
        self.tmpVec = []
        self.log.info('created combischeme with \n'+str(self.getCombinationDictionary()).replace(', (','\n('))

    def extendSchemeExtraLayers(self,activeSet,numExtraLayers):
        extraGrids = activeSet.getExtraGrids(numExtraLayers)
        for grid in extraGrids:
            self.dictOfScheme[grid] = 0
        
    def getHierarchicalCoefficients(self):
        downSet = self.getUnifiedDownSet()
        w = {}
        for grid in downSet:
            w[grid] = 0
            largerSubspaces = self.getLargerSubspaces(grid)
            for largerSpace in largerSubspaces:
                if largerSpace in self.dictOfScheme:
                    w[grid] += self.dictOfScheme[largerSpace] 
        return w
    
    def evalGCP(self,w):
        Q = 0
        for wi in w:
            Q+= 4**(-sum(wi))*w[wi]
        return Q
    
    def updateFaults(self,faults):
        maxLevelSum = sum(list(self.activeSet)[0]) #all grids in the active set have the same levelSum
        noRecomputationLayers = 2
        removeInds = []
        for i in range(len(faults)):
            if sum(faults[i]) <= maxLevelSum - noRecomputationLayers:
                removeInds.append(i)
        newFaults = list(faults)
        for i in range(len(removeInds)):
            del newFaults[removeInds[i]-i]
        return newFaults

    def divideFaults(self,faults):
        maxLevelSum = sum(list(self.activeSet)[0])
        qZeroFaults = []
        qOneFaults = []
        for fault in faults:
            l1 = sum(fault)
            if l1 == maxLevelSum:
                qZeroFaults.append(fault)
            else:
                qOneFaults.append(fault)
        return qZeroFaults, qOneFaults
    
    def partitionFaults(self,faults, w):
        isPartitioned = {}
        partitionedFaults = []
        if len(faults) == 1:
            isPartitioned[faults[0]] = True
            return isPartitioned
        
        neighborsDict = {}
        for fault in faults:
            isPartitioned[fault] = False
            nbrList = self.getNeighbors(fault, w)
            s = set()
            if len(nbrList) == 0:
                s.add(tuple(fault))
                neighborsDict[fault] = s
            else:
                nbrList.append(fault)
                nbrList = set(nbrList)
                neighborsDict[fault] = nbrList

        for fault in faults:
            if not isPartitioned[fault]:
                isPartitioned[fault] = True
                currentNbrs = neighborsDict[fault]
                self.findCurrentPartition(fault, currentNbrs,faults,
                                     neighborsDict,isPartitioned,partitionedFaults)
        return partitionedFaults
                
    def findCurrentPartition(self,fault,currentNbrs,faults,neighborsDict,isPartitioned,partitionedFaults):
        self.tmpVec.append(fault)
        for flt in faults:
            if not isPartitioned[flt]:
                interSet = currentNbrs.intersection(neighborsDict[flt])
                if len(interSet) != 0:
                    isPartitioned[flt] = True
                    unionSet = currentNbrs.union(neighborsDict[flt])
                    currentNbr = unionSet
                    self.findCurrentPartition(flt,currentNbr,faults,
                                            neighborsDict,isPartitioned,partitionedFaults)
        if len(self.tmpVec) != 0:
            partitionedFaults.append(self.tmpVec)
        self.tmpVec = []
        return partitionedFaults
    
    def generateCasesGCP(self,faults):
        w = OrderedDict()
        for grid in self.getUnifiedDownSet():
            w[grid] = 1
            
        qZeroFaults, qOneFaults = self.divideFaults(faults)
        
        for fault in qZeroFaults:
            del w[fault]
        
        partitionedFaults = self.partitionFaults(faults, w)
        
        allW = []
        allWDicts = []
        if len(faults) == 1:
            oneW = self.generateOneCaseGCP(faults[0], w)
            for case,numCase in zip(oneW,list(range(len(oneW)))):
                wTmp = w.copy()
                for grid in case:
                    wTmp[grid] = oneW[numCase][grid]
                allWDicts.append(wTmp)
            return allWDicts
        else:
            for part in partitionedFaults:
                currentWs = []
                result = []
                onePartW = self.generateOnePartitionCasesGCP(part,w)
                for case in onePartW:
                    currentWs.append(list(case.values()))
                if len(allW) == 0:
                    allW = currentWs
                else:
                    for pairs in it.product(allW,currentWs):
                        result.append(list(np.prod(np.array(pairs),axis=0)))
                    allW = result

        for case,numCase in zip(allW,list(range(len(allW)))):
            wTmp = w.copy()
            for grid,j in zip(w,list(range(len(w)))):
                wTmp[grid] = allW[numCase][j]
            allWDicts.append(wTmp)
        return allWDicts
    
    def generateOnePartitionCasesGCP(self,faults,w):
        qZeroFaults, qOneFaults = self.divideFaults(faults)

        allW = []
        if len(qOneFaults) == 0:
            allW.append(w.copy())
            return allW
        
        numFaults = len(qOneFaults)
        casesGCP = {}
        indexList = []
        for fault in qOneFaults:
            casesGCP[fault] = self.generateOneCaseGCP(fault, w)
            lenCase = len(casesGCP[fault])
            indexList.append(list(range(lenCase)))

        allW = []
        canBeChanged = {}
        for wi in w:
            canBeChanged[wi] = True
        for index in it.product(*indexList):
            tmpW = w.copy()
            changedGrids = []
            for i,fault in zip(list(range(numFaults)),qOneFaults):
                success = True
                currentW = casesGCP[fault][index[i]]
                for cw in currentW:
                    if currentW[cw] != tmpW[cw]:
                        if canBeChanged[cw]:
                            tmpW[cw] = currentW[cw]
                            canBeChanged[cw] = False
                            changedGrids.append(cw)
                        else:
                            success = False
                            break
                    else:
                        canBeChanged[cw] = False
                        changedGrids.append(cw)
                if not success:
                    break
            if success:
                allW.append(tmpW)
            for grid in changedGrids:
                canBeChanged[grid] = True
        
        return allW
    
    def getNeighbors(self,grid,allGrids):
        neighbors = []
        for k in range(self.getDim()):
            nbr = list(grid)
            nbr[k] +=1
            nbr = tuple(nbr)
            if nbr in allGrids:
                neighbors.append(nbr)
        return neighbors
    
    def generateOneCaseGCP(self,fault,w):

        cases = []
        neighbors = self.getNeighbors(fault, w)
        
        # case w_i = 0
        caseGCP = {}
        caseGCP[fault] = 0
        for nbr in neighbors:
            caseGCP[nbr] = 0
        cases.append(caseGCP)
        # cases k=1,...,d
        caseGCP = {}
        caseGCP[fault] = 1
        for nbr in neighbors:
            caseGCP[nbr] = 0
        for nbr in neighbors:
            caseGCP[nbr] = 1
            cases.append(caseGCP.copy())
            caseGCP[nbr] = 0      
        return cases
    
    def chooseBestGCP(self,allW):
        if len(allW) == 1:
            return allW[0]
        Q_max = -np.inf
        for w in allW:
            Q = self.evalGCP(w)
            if Q > Q_max:
                Q_max = Q
                W_max = w
        return W_max
        
    def buildGCPmatrix(self,w,faults):
        maxLevelSum = sum(list(self.activeSet)[0])
        qZeroFaults = []
        for fault in faults:
            l1 = sum(fault)
            if l1 == maxLevelSum:
                qZeroFaults.append(fault)
        indsDict = {}
        i = 0
        for w_i in w:
            indsDict[w_i] = i
            i+=1
        M = np.matlib.identity(len(w))
        dictOfSchemeFaults = set(self.getUnifiedDownSet()).difference(set(qZeroFaults))
        for w_i in w:
            largerGrids = set(self.getLargerSubspaces(w_i)).intersection(dictOfSchemeFaults)
            for grid in largerGrids :
                M[indsDict[w_i],indsDict[grid]] = 1
        return M
    
    def solveGCP(self,M,w):
        return np.linalg.solve(M, w)
    
    def recoverSchemeGCP(self,faults):
        
        # faults below layer two will e recalculated
        newFaults = self.updateFaults(faults)
        
        # if all faults occur below the second diagonal, return these
        if len(newFaults) == 0:
            return faults
        # generate all the w coeff. vectors that lead to 
        # a scheme ommitting the faults 
        allW = self.generateCasesGCP(newFaults)
        
        # out of those, choose the best one
        wBest = self.chooseBestGCP(allW)
        
        # calculate the c coeffs: first generate
        # the coefficient matrix
        M = self.buildGCPmatrix(wBest,faults)
        
        # solve for c (combination coeffs)
        c = self.solveGCP(M,list(wBest.values()))
        
        # calculate new maxLevelSum
        maxLevelSum = 0
        self.activeSet = set()
        for grid,i in zip(wBest,list(range(len(c)))):
            l1 = sum(grid)
            if l1 > maxLevelSum and c[i] != 0:
                maxLevelSum = l1
   
        # update the combischeme
        for grid,i in zip(wBest,list(range(len(c)))):
            l1 = sum(grid)
            if l1==maxLevelSum and c[i] == 1:
                self.activeSet.add(grid)
        
        self.updateDict()
        
        # check that no q=0,1 faults appear on the scheme (this should always hold!)
        # return grids to be recomputed
        if len(set(self.dictOfScheme.keys()).intersection(set(newFaults))) == 0:
            return set(self.dictOfScheme.keys()).intersection(set(faults))
        else:
            return "Error"
class combinationSchemeSGPP(combinationSchemeBase):
    
    def __init__(self, pysgppScheme):
        combinationSchemeBase.__init__(self)
        self.pysgppScheme = pysgppScheme
        self.updateDict()
        
    def updateDict(self):
        levels = self.pysgppScheme.getLevels()
        coeff = self.pysgppScheme.getCoef()
        self.combiDict = {}
        for l,c in zip(levels,coeff):
            self.combiDict[l] = c
        
    def changeCombischeme(self):
        raise NotImplementedError('Please implement this')
        
    def getCombinationDictionary(self):
        return self.combiDict
    
    def getDim(self):
        return self.pysgppScheme.getDim()
    
def showScheme3D(schemeobject,coords=[0,1,2]):
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    scheme = schemeobject.getCombinationDictionary()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    col = ['blue','red','black','green']
    colors = it.cycle(col)
    
    for s in scheme:
        if scheme[s]>0:
            c = col[0]
        else:
            c = col[1]  
        ax.scatter(s[coords[0]],s[coords[1]],s[coords[2]],s=abs(scheme[s])*300.0,c=c)
    
    plt.show()
        
