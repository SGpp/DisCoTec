'''
Created on Jan 13, 2014

@author: kowitz
'''
import numpy as np
import logging


class DataManagementBase(object):
    def getSliceTuple6D(self, coordDict):
        printCoord = set(self.coordNames).difference(list(coordDict.keys()))

        self.log.debug('Create Slice for ' + str(printCoord))
        t = [None for i in range(6)]
        for c in coordDict:
            t[self.coordNames.index(c)] = coordDict[c]

        for c in printCoord:
            t[self.coordNames.index(c)] = slice(0, self.data.shape[self.coordNames.index(c)])
#         t = t[::-1]
        self.log.debug(t)
        return t

    def getSliceTupleBallooning(self, coordDict,data):
        printCoord = set(self.coordNamesBallooning).difference(list(coordDict.keys()))
        self.log.debug('Create Slice for ' + str(printCoord))
        t = [None for i in range(self.dim)]
        for c in coordDict:
            t[self.coordNamesBallooning.index(c)] = coordDict[c]

        for c in printCoord:
            t[self.coordNamesBallooning.index(c)] = slice(0, data.shape[self.coordNamesBallooning.index(c)])
#         t = t[::-1]
        self.log.debug('test'+str(t))
        return t

    def getSlice(self, coordDict):
        t = self.getSliceTuple6D(coordDict)
        return self.data[t]


    def getBallooningSlice(self, coordDict,inversion=True):
        data = self.getBallooning(inversion=inversion)
        self.log.debug(str(data.shape))
        t = self.getSliceTupleBallooning(coordDict,data)
        self.log.debug(t)
        return data[t]


class GeneField(DataManagementBase):
    coordNames = ['x', 'y', 'z'][::-1]
    coordNamesBallooning = ['xz', 'y'][::-1]
    def __init__(self,data):
        self.data = data
        self.log = logging.getLogger(__name__)
        self.dim = 3
        if len(self.data.shape) != self.dim:
            raise TypeError('not the right amount of dimensions')

        self.log.debug('Create GeneFields')
        self.shapeDict = {}
        for c,s in zip(self.coordNames,self.data.shape):
            self.shapeDict[c]=s

    def getBallooning(self,inversion=True):
        s = self.data.shape
        d= self.data.copy()
        if inversion:
            d[:,:,::2]*=-1.0
        d= d[  :-1, :, :].transpose(( 1, 2, 0)).reshape(s[1], s[2] * (s[0] - 1))
        return d



class GenePdf(DataManagementBase):
    '''
    a convenience function to access gene results
    '''
    coordNames = ['x', 'y', 'z', 'v', 'w', 's'][::-1]
    coordNamesBallooning = ['xz', 'y', 'v', 'w', 's'][::-1]

    def __init__(self, data):
        self.data = data
        self.dim = 5
        self.log = logging.getLogger(__name__)
        if len(self.data.shape) != 6:
            raise TypeError('array does not have the right amount of dimensions')
        self.log.debug('Created GenePdf')
        self.shapeDict={}
        for c,s in zip(self.coordNames,self.data.shape):
            self.shapeDict[c]=s




    def getBallooning(self,inversion=True):
        s = self.data.shape
        d= self.data.copy()
        if inversion:
            d[:,:,:,:,:,::2]*=-1.0
        d= d[:, :, :, :, :, :].transpose((0, 1, 2, 4, 5, 3)).reshape(s[0], s[1], s[2], s[4], s[5] * (s[3]))
        return d







if __name__ == '__main__':
    
    log = logging.getLogger(__name__)
    log.addHandler(logging.StreamHandler())
    log.setLevel(logging.DEBUG)
    
    a = np.random.rand(3, 4, 5, 6, 7, 8)
    pdf = GenePdf(a)
    
    coordDict = {'x':1, 'y':1, 'w':1, 's':0}
    coordDict2 = {'xz':1}
    print(pdf.getSlice(coordDict) - pdf.data[0, 1, :, :, 1, 1])
    print(pdf.getBallooningSlice(coordDict2).shape)
    
