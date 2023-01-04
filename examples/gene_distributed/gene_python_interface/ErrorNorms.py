'''
Created on Jan 10, 2014

@author: kowitz
'''
import numpy as np

## from https://www.rocq.inria.fr/modulef/Doc/GB/Guide6-14/node21.html

def L2Error(approx,reference,species_index=None):
    if species_index is not None:
        return np.linalg.norm(approx[species_index].flatten()-reference[species_index].flatten(),2)
    else:
        return np.linalg.norm(approx.flatten()-reference.flatten(),2)

def L2RelativeError(approx,reference):
    return L2Error(approx, reference)/ np.linalg.norm(reference.flatten(),2)

def L1Error(approx,reference):
    return np.linalg.norm(approx.flatten()-reference.flatten(),1)

def L1RelativeError(approx,reference):
    return L1Error(approx, reference)/np.linalg.norm(reference,1) / len(approx)

def L2AverageError(approx,reference):
    return L2Error(approx, reference)/np.sqrt(len(approx))

def L1AverageRelativeError(approx,reference):
    return np.sum(abs(approx-reference)/abs(reference))/len(approx)

def L1AverageError(approx,reference):
    assert approx.size == reference.size
    return L1Error(approx,reference)/approx.size

def LInfError(approx,reference):
    return np.linalg.norm(approx.flatten()-reference.flatten(),np.inf)

def LInfRelativeError(approx,reference):
    return LInfError(approx, reference)/ np.linalg.norm(reference.flatten(),np.inf)

def L1AverageErrorGene(approx,reference):
    '''
     That makes the error computation just on the inner cells and not on the artificial zero boundaries
    '''
    assert approx.size==reference.size
    # return L1Error(approx[:,:-1,1:-1,:-1,:,:],reference[:,:-1,1:-1,:-1,:,:])/approx[:,:-1,1:-1,:-1,1,0].size
    return L1Error(approx[:,:-1,1:-1,:-1,:,:],reference[:,:-1,1:-1,:-1,:,:])/approx[:,:-1,1:-1,:-1,:,:].size

def L2AverageErrorGene(approx,reference):
    assert len(approx.shape)==6
    assert len(reference.shape)==6
    ap=   approx[:,:-1,1:-1,:-1,1:-1,:]
    re=reference[:,:-1,1:-1,:-1,1:-1,:]
    return L2Error(ap,re)/np.sqrt(ap[:,:,:,:,0,0].size)

def L2AverageErrorGeneField(approx,reference):
    assert len(approx.shape)==3
    assert len(reference.shape)==3

    ap = approx[:-1,1:-1,:]
    re = reference[:-1,1:-1,:]

    return L2Error(ap,re)/np.sqrt(ap.shape[0])


if __name__=='__main__':
    def f1(x):
        return np.exp(-x)
    def f2(x):
        return np.exp(-10.*x)
    import pylab as py
    
    x = np.arange(0,1,0.1)
    
    for f in [f1,f2]:
        yRef = f(x)
        yApprox = f(x)+0.1
        plotArg = ['o-']
        py.figure()
        py.plot(yRef,*plotArg)
        py.plot(yApprox,*plotArg)
        print(L2Error(yApprox, yRef))
        print(L2RelativeError(yApprox, yRef))
        print(L1Error(yApprox, yRef))
        print(L1RelativeError(yApprox, yRef))
        print(L2AverageError(yApprox, yRef))
        print(L1AverageRelativeError(yApprox, yRef))
        print('--------------------------------------')
    py.show()
    
    
    
