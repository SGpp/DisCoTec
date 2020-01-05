'''
Created on Jan 5, 2014

@author: ckow
'''
import numpy as np
import itertools as it
import abc
import logging
import random
from argparse import ArgumentError

class ActiveSetBase(object, metaclass=abc.ABCMeta):
    def __init__(self):
        self.log = logging.getLogger(__name__)
        self.log.debug('Created ActiveSet Object')
       
    @abc.abstractmethod
    def getActiveSet(self, *args):
        pass

       
class ActiveSetTiltedPlane(ActiveSetBase):
    def __init__(self, planeNormalVector, anchorVector):
        super(ActiveSetTiltedPlane, self).__init__()
        self.n = np.array(planeNormalVector)
        self.anchor = np.array(anchorVector)
        self.setRoundingMethod(np.floor)
        self.lmax=np.ones(len(self.n))*20
    
    def setRoundingMethod(self, method):
        self.roundingMethod = method
    
    def getIntersectionPointsWithCoordinateAxes(self, currentDim):
        indexSet = list(range(currentDim)) + list(range(currentDim + 1, len(self.n)))
        ret= np.sum(self.n[indexSet] * (self.anchor[indexSet] - 1.)) / self.n[currentDim] + self.anchor[currentDim]
        return ret
    
    def getLmax(self):
        erg = []
        for i in range(len(self.n)):
            erg.append(int(np.ceil(self.getIntersectionPointsWithCoordinateAxes(i))))
        return tuple(erg)
        
    def getIntersectionPoints(self, currentDim):
        a = np.array(self.anchor)
        n = np.array(self.n)
        dim = currentDim
        lmax = self.getLmax()
        
        dimIndices = list(range(dim)) + list(range(dim + 1, len(lmax)))
        maxInd = np.array(lmax)[dimIndices]
        coords = [list(range(1, l + 1)) for l in maxInd]
        tuples = list(it.product(*coords))
        
        erg = {}
        for i in range(len(tuples)):
            t = tuples[i]
            s = np.array(t)
            ind = tuple(s - np.ones(len(s), dtype=int))
            erg[t] = np.sum((a[dimIndices] - s) * n[dimIndices]) / (n[dim] * 1.0) + a[dim]
#             erg[t] = int(np.floor(erg[t]))
            
        ergSet = set()
        for e in erg:
            eList = list(e)
            listAll = eList[:dim] + [erg[e]] + eList[dim:]
            if (np.array(listAll) >= 0).all():
                ergSet.add(tuple(listAll))
        return ergSet

    def getActiveSet(self):
        erg = set()
        for i in range(len(self.n)):
            s = self.getIntersectionPoints(i)
            for grid in s:
                g = (self.roundingMethod(grid)).astype(int)
                if (g >= 1).all():
                    erg.add(tuple(g))
        return erg
            
class RootActiveSet(ActiveSetBase):
    def __init__(self,lmax,lmin):
        super(ActiveSetBase,self).__init__()
        self.lmin = lmin
        self.lmax = lmax

    def getActiveSet(self, *args):
        d = len(self.lmax)
        erg = set()
        for i in range(d):
            buf = list(self.lmin)
            buf[i]=self.lmax[i]
            erg.add(tuple(buf))
        erg.add(self.lmin)
        return erg

class ClassicDiagonalActiveSet(ActiveSetBase):
    def __init__(self,lmax,lmin=None,diagonalIndex=0):
        super(ClassicDiagonalActiveSet, self).__init__()
        self.lmax = lmax
        if lmin==None:
            self.lmin = [1 for i in range(len(self.lmax))]
        else:
            self.lmin = lmin
        self.setDiagonalIndex(diagonalIndex)
        
    def setDiagonalIndex(self,i):
        self.diagonalIndex = i

    def getMinLevelSum(self):
        return self.getLevelMinima().sum()
#         return self.getMaxLevel()-1+len(self.lmax)
    
    def getMaxLevel(self):
        return np.max(self.lmax)
    
    def getLevelMinima(self):
        maxInd = np.argmax(np.array(self.lmax)-np.array(self.lmin))
        lm = np.array(self.lmin)
        lm[maxInd]=self.lmax[maxInd]
        return lm
        
    
    def getActiveSet(self):
        listOfRanges=[list(range(0,self.lmax[i]+1)) for i in range(len(self.lmax))]
        listOfAllGrids = list(it.product(*listOfRanges))
        s = set()
        levelSum = self.getMinLevelSum()+self.diagonalIndex
        for grid in listOfAllGrids:
            if (np.sum(grid)==levelSum and (np.array(grid)>=np.array(self.lmin)).all()):
                s.add(grid)
                
        return s
    
    def getExtraGrids(self,numExtraDiags):
        listOfRanges=[list(range(1,self.lmax[i]+1)) for i in range(len(self.lmax))]
        listOfAllGrids = list(it.product(*listOfRanges))
        diff = np.array(self.lmax)-np.array(self.lmin)
        dim = len(self.lmax)-len(np.where(diff == 0)[0])
        lmin_eff = np.zeros(dim)
        lmax_eff = np.zeros(dim)
        j = 0
        for i in range(len(self.lmin)):
            if diff[i] != 0:
                lmin_eff[j] = self.lmin[i]
                lmax_eff[j] = self.lmax[i]
                j+=1

        s = set()
        for q in range(dim,dim+numExtraDiags):
            levelSum = np.array(lmax_eff).max()-np.array(lmin_eff).max()+np.array(self.lmin).sum()-q
#             print self.getMinLevelSum(),self.getMaxLevel().max(), np.array(self.lmin).max(), np.array(self.lmin).sum(),q,levelSum
            for grid in listOfAllGrids:
                if (np.sum(grid)==levelSum and (np.array(grid)>=np.array(self.lmin)).all()):
                    s.add(grid)
        return s
    
    def getEffLmin(self):
        diff = np.array(self.lmax)-np.array(self.lmin)
        activeInds = np.where(diff != 0)[0]
        lmin = np.array(tuple(self.lmin[i] for i in activeInds))
        lmax = np.array(tuple(self.lmax[i] for i in activeInds))

#         For lmin with extra 1's:
        erg = np.ones(len(self.lmin))
        ergTmp = lmax - np.min(lmax - lmin)*np.ones(len(lmax))
        for i in range(len(activeInds)):
            erg[activeInds[i]] = ergTmp[i]
        return  tuple(map(int,erg))

        # For lmin without extra 1's:
#         erg = lmax - np.min(lmax - lmin)*np.ones(len(lmax))
#         return  tuple(map(int,erg))
    
class ThinnedDiagonalActiveSet(ClassicDiagonalActiveSet):
    def __init__(self,lmax,thinningFactor=None,thinningNumber=None,lmin=None, diagonalIndex=0):
        super(ThinnedDiagonalActiveSet, self).__init__(lmax,lmin,diagonalIndex)
        if thinningFactor==None and thinningNumber==None:
            raise ArgumentError('Please fix a thinningNumber or thinningFactor')
        if thinningFactor!=None and thinningNumber!=None:
            raise ArgumentError('Please fix either thinningNumber or thinningFactor, not both')
        
        if thinningFactor!=None:
            self.thinningFactor = thinningFactor
            self.getActiveSet = self.getActiveSetThinningFactor
        else:
            self.thinningNumber = thinningNumber
            self.getActiveSet = self.getActiveSetFixedNumber
        
    
    def isSetValid(self,activeSet):
        ##get maxima
        it=iter(activeSet)
        maxInd=list(next(it))
        for s in activeSet:
            for i in range(len(s)):
                if s[i]>maxInd[i]:
                    maxInd[i]=s[i]
        for el,l in zip(maxInd,self.lmax):
            if el!=l:
                return False
        return True
    
    def removeNGrids(self,s,nGrids):
        if nGrids>len(s)-len(self.lmax):
            raise ArgumentError('too high thinning, cannot remove more grids than are there')
        counter=0
        while(1):
            counter+=1
            randomSample =set( random.sample(s,len(s)-nGrids))
            if self.isSetValid(randomSample):
                self.log.debug(str(counter)+' draws where required to find the active Set')
                return randomSample
    
    def getActiveSetFixedNumber(self):
        classicSet = ClassicDiagonalActiveSet.getActiveSet(self)
        return self.removeNGrids(classicSet, self.thinningNumber)
                
            
        
    def getActiveSetThinningFactor(self):
        classicSet = ClassicDiagonalActiveSet.getActiveSet(self)
        nGrids = len(classicSet)
        rounds = nGrids-int(round(nGrids*self.thinningFactor))
        return self.removeNGrids(classicSet, rounds)
                
        
        

def showActiveSet3D(listOfsets):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = it.cycle(['black','red','blue','green'])
    
    for sets,color in zip(listOfsets,colors):
        s = np.array(list(sets))
    
        ax.scatter(s[:,0],s[:,1],s[:,2],c=color)
    
    plt.show()


    
    
if __name__ == '__main__':
#     print ActiveSet.getIntersectionPoints((0.5,1), (1,5), 0, (10,10))
#     print ActiveSet.getIntersectionPoints((0.5,1), (1,5), 1, (10,10))

#     print ActiveSet.getIntersectionPoints((0.5, 1, 1), (3, 1, 1) , 0, (10, 10, 10))
    
#     actSet = ActiveSetTiltedPlane((1, 1.,0.5), (2,2,2))
#     actSet.setRoundingMethod(np.ceil)
#     
#     print actSet.getActiveSet()

#     s = ClassicDiagonalActiveSet((10,10,10),lmin=(1,1,1),diagonalIndex=0)
#     s2 = ClassicDiagonalActiveSet((10,10,10),lmin=(3,4,5),diagonalIndex=0)
#     print s.getMinLevelSum()
#     print s2.getMinLevelSum()
#     activeSet = s.getActiveSet()
    log=logging.getLogger(__name__)
    log.addHandler(logging.StreamHandler())
    log.setLevel(logging.DEBUG)
    
#     s= ThinnedDiagonalActiveSet((10,10,10),thinningFactor=0.9,lmin=(3,1,1),diagonalIndex=0)
#    s= ThinnedDiagonalActiveSet((10,10,10),thinningNumber=20,lmin=(3,1,1),diagonalIndex=0)
    s= ClassicDiagonalActiveSet((5,6,5),lmin=(3,5,3))
    activeSet=s.getActiveSet()
    
    
    listOfSets=[activeSet]
#     for i in range(5):
#         s.setDiagonalIndex(i+1)
#         listOfSets.append(s.getActiveSet())
    
    
    
    
#     print activeSet
    showActiveSet3D(listOfSets)
    
