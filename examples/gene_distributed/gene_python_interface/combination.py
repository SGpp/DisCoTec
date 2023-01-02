'''
Created on Nov 14, 2013

@author: kowitz
'''
import logging
from pysgpp import CombiSchemeBasis

class combinationBase(object):
    '''
    class doing the combination of various objects
    '''


    def __init__(self,scheme):
        '''
        Constructor
        '''
        self.scheme = scheme
        self.log = logging.getLogger(__name__)
        self.log.info('Created combigrid')
        
        
        
    def updateScheme(self,newScheme):
        self.scheme = newScheme
        
    def getLevels(self):
        return self.scheme.getLevels()
    
    