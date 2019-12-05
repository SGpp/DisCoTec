'''
Created on Jan 28, 2014

@author: ckow
'''
import logging,sys
reload(logging)
import string


class AllLoggers(object):
    '''
    classdocs
    '''
    standardLoggers = ['ActiveSetFactory',\
                       'combinationScheme',\
                       'combineGrids',\
                       'geneCombinationGrid',\
                       'GeneInterface',\
                       'GenePdf',\
                       'combiGrid']

    standardFormatter = logging.Formatter('%(module)20s:%(funcName)20s:%(lineno)4d:   %(message)s')

    def __init__(self,level='info',loggerNames=standardLoggers,\
                 formatter=standardFormatter,handler=logging.StreamHandler(sys.stdout), isIpython = False):
        '''
        Constructor
        '''

        self.handler = logging.NullHandler()
        if isIpython is not True:
            self.formatter= formatter
            self.handler = handler
            self.handler.setFormatter(self.formatter)

        main = 'main'
        self.setLoggerDict(loggerNames+[main])    
        self.setLevel(level)
        self.mainLog = self.loggers[main]
        
        self.info=self.mainLog.info
        self.debug=self.mainLog.debug
        self.error=self.mainLog.error
        self.warning=self.mainLog.warning



    def setLevel(self,level):
        if level=='info':
            self.level=logging.INFO
        elif level=='debug':
            self.level=logging.DEBUG
        elif level=='warn':
            self.level=logging.WARN
        elif level=='error':
            self.level=logging.ERROR
            
        for l in self.loggers.values():
            l.setLevel(self.level)
        
    def setLoggerDict(self,loggerNames):        
        self.loggers = {}
        for logName in loggerNames:
            self.loggers[logName] = logging.getLogger(logName)
            self.loggers[logName].addHandler(self.handler)
            
    def getLevel(self):
        return self.level
    

        