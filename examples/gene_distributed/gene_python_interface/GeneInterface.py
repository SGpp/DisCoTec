'''
Created on Apr 26, 2013

@author: ckow
'''
import exceptions
import shutil
import pickle
import numpy as np
import logging
import subprocess
import os
import struct
import sys
import abc
from numpy.lib.scimath import sqrt
import shlex


# all_loggers = ['parameterFileObject', \
#                'parameterFileObjectForChangingStuff', \
#                'geneEnvironment2', \
#                'checkpoint2Base', \
#                'checkpoint2Local']

def num(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            try:
                return complex(s)
            except ValueError:
                try:
                    return complex(s.replace('i', 'j'))
                except ValueError:
                    string = str(((s).strip()).strip("'")).split(' ')
                    if len(string) > 1:
                        ret = list(map(num, string))
                        return ret
                    else:
                        return str(((s).strip()).strip("'"))

class commonFileNames(object):
    '''
    class containing static variables with the names of the ouputfiles (they shouldn't change)
    '''

    parameterFile = 'parameters'
    checkpointFile = 'checkpoint'
    fieldFile = 'field.dat'
    omegaFile = 'omega.dat'
    eigenvalueFile = 'eigenvalues.dat'
    stdOutFile = 'GENEstdout'


class parameterFileObject(object):

    def __init__(self, parameterFile):
        self.log = logging.getLogger(__name__)
        self.log.debug('created strings')
        self.inputFile = parameterFile
        self.speciesName = 'name'
        self.multipleParameters = [ 'omn', 'omt', 'mass', 'temp', 'dens', 'charge']
        self.identifierKey = '#key#'
        self.parameterDictionary = self.readFile(self.inputFile)
        
        self.importantPar = ['nx0', 'nky0', 'nz0', 'nv0', 'nw0', 'comp_type', 'kymin']
        
        self.log.debug('created parameterfileobjject')


    def readFile(self, inputFile):
        self.log.debug('start reading file')
        fr = open(inputFile, 'r')
        parameterDictionary = {}
        for line in fr:
            if line.find('=') >= 0:
                # self.log.debug(line)
                beforeComment = line.split('!', 1)[0]
                if beforeComment == '':
                    continue
                varName = beforeComment.split('=', 1)[0].strip()
                value = beforeComment.split('=', 1)[1]  # .replace('i','j')
                value = value.strip()
                
                if varName.strip() == self.speciesName:
                    currentSpeciesName = value
                    varName = currentSpeciesName + self.identifierKey + varName
                    self.log.debug('Found name')
                
                if varName in self.multipleParameters:
                    varName = currentSpeciesName + self.identifierKey + varName
                    
                parameterDictionary[varName] = value
                # #setattribute create .name access to variable name --> maybe not so good
#                setattr(self, varName, num(value))
        fr.close()
        return parameterDictionary

    def getDictionaryCopy(self):
        return self.parameterDictionary.copy()

    def getPar(self, par):
        return num(self.parameterDictionary[par])

    def getMultiplePar(self, par):
        returnDict = {}
        for p in par:
            returnDict[p] = self.parameterDictionary[p]
        return returnDict
    
    def diffTo(self, parFileObject):
        setSelf = set(self.parameterDictionary.keys())
        setOther = set(parFileObject.parameterDictionary.keys())
        
        self.log.debug(str(setSelf - setOther))
        self.log.debug(str(setOther - setSelf))

    def isSameAs(self, parFileObject):
        ''' Checks if the other parameterfile has the same content for certain paramters'''
        for key in self.importantPar:
            if self.parameterDictionary[key] != parFileObject.parameterDictionary[key]:
                self.log.debug('Difference found for ' + key + ' with ' + str(self.parameterDictionary[key]) + '   ' + str(parFileObject.parameterDictionary[key]))
                return False
        return True
    
    def getGridSizesDict(self):
        return self.getMultiplePar(['n_spec','nw0','nv0','nz0','nky0','nx0'])
        
        

class parameterFileObjectForChangingStuff(parameterFileObject):

    def __init__(self, parameterFile, outputParameterFile):
        parameterFileObject.__init__(self, parameterFile)
        self.outputFile = outputParameterFile
        # self.outParDict = self.parameterDictionary.copy()
        self.saveDict = self.parameterDictionary.copy()
#         self.log = logging.getLogger(self.__class__.__name__)

    def updateParameter(self, par, value):
        ''' just for existing parameters due to namelists!'''
        outParDict = {}
        if par not in self.parameterDictionary:
            raise KeyError()
        self.parameterDictionary[par] = value
        outParDict[par] = value
        self.updateFile()
        self.log.debug('updated entry: ' + str(par) + ' with ' + str(value))

    def updateParameters(self, parDict):
        ''' just for existing parameters due to namelists!'''
        outParDict = {}
        if set(self.parameterDictionary.keys()) != set(parDict.keys()):
            raise KeyError()
        outParDict.update(parDict)
        self.parameterDictionary.update(parDict)
        self.updateFile()

    def updateFile(self):
        f_save = open(self.inputFile, 'r')
        f_output = open(self.outputFile, 'w')
        outParDict = self.parameterDictionary
        for line in f_save:
            if line.find('=') >= 0:
                par_name = line.split('=')[0].strip()
                value = line.split('=')[1].split('!')[0]
#                 self.log.debug(par_name)
                if par_name == self.speciesName:
                    currentSpeciesName = value
                
                if par_name in self.multipleParameters:
                    par_name_ext = currentSpeciesName + self.identifierKey + par_name
                else:
                    par_name_ext = par_name
                
                if par_name_ext in outParDict:
                    f_output.write(par_name + ' = ' + outParDict[par_name_ext] + '\n')
#                    if type(outParDict[par_name]) == str:
#                        f_output.write(par_name + ' = \'' + outParDict[par_name] + '\'\n')
#                    else:
#                        f_output.write(par_name + ' = ' + str(outParDict[par_name]) + '\n')
                else:
                    f_output.write(line)
            else:
                f_output.write(line)

        f_save.close()
        f_output.close()
        self.parameterDictionary.update(outParDict)

    def getChanges(self):
        # saveDict = self.readFile(self.inputFile)
        diff = {}
        for k in list(self.saveDict.keys()):
            if self.saveDict[k] != self.parameterDictionary[k]:
                diff[k] = (self.saveDict[k], self.parameterDictionary[k])
        return diff

class geneEnvironment2(object):
    '''
    Environment containing all information about GENE. Refurbished version of the original environment
    '''
    
    geneStdout='GENEstdout'


    def __init__(self, genePath, problemFolder, diagDir, chptDir=None):
        '''
        Constructor
        '''
        self.log = logging.getLogger(__name__)
#         self.log = logging.getLogger('GeneInterface')
#         print self.__class__.__name__

        # #test for existing folder
        ls = os.listdir(genePath)
        if not problemFolder.replace('/', '') in ls:
            raise IOError('problem Folder not existing')

        self.runCommand = 'mpirun -np 1 ./gene_new_machine'
        self.originalFile = 'parameters.save'
        self.parFile = 'parameters'
        self.diagDir = diagDir

        if chptDir == None:
            self.chptDir = self.diagDir
        else:
            self.chptDir = chptDir
        self.genePath = genePath
        self.problemFolder = problemFolder

        self.startFromCheckpoint = False
        
        self.outputName = 'GENEstdout'

        # os.chdir(self.getProblemPath())

        self.parFileIO = \
            parameterFileObjectForChangingStuff(self.genePath + self.problemFolder + self.originalFile, \
                                                self.genePath + self.problemFolder + self.parFile)

        self.parFileIO.updateParameter('diagdir', '"' + self.diagDir + '"')
        self.log.info('GeneEnv created with '+self.genePath+self.problemFolder)


    def update(self, newPath):
        '''
        update the diagdir
        '''
        self.parFileIO.updateParameter('diagdir', '"'+newPath+'"')
        self.diagDir = newPath
        
    def setGridSizeWithLevel(self,level):
        self.parFileIO.updateParameter('nx0', str(2**level[4]+1))
        self.parFileIO.updateParameter('nky0', str(2**level[3]-1))
        self.parFileIO.updateParameter('nz0', str(2**level[2]))
        self.parFileIO.updateParameter('nv0', str(2**level[1]))
        self.parFileIO.updateParameter('nw0', str(2**level[0]))

    def getProblemPath(self):
        return self.genePath + self.problemFolder

    def buildGENE(self):
        currentDir = os.getcwd()
        os.chdir(self.genePath + self.problemFolder)

        self.log.info('Starting compilation!')
        subprocess.call(['make', '-f', '../makefile', '-j4', \
                         'MACHINE=ubuntu_gfortran', \
                         'SLEPC_DIR=/home/kowitz/work_fast/slepc', \
                         'PETSC_DIR=/home/kowitz/work_fast/petsc', \
                         'PETSC_ARCH=arch-linux2-c-debug'])
        self.log.info('Compilation finished')

        os.chdir(currentDir)

    def runGENE(self, commandLineParm=''):
        currentDir = os.getcwd()

        f_stdout = open(self.diagDir + self.outputName, 'w')
        self.log.debug('Stdout will be written to: ' + self.diagDir + self.outputName)
        os.chdir(self.genePath + self.problemFolder)
        self.log.info('Start GENE in ' + os.getcwd()+' with grid sizes \n'+ \
                      str(self.parFileIO.getGridSizesDict()))
        runCommand = self.runCommand + ' ' + str(commandLineParm)
        self.log.debug('Starting command: ' +str(shlex.split(runCommand)))
        self.log.debug('Starting command all:'+str(runCommand))
        self.log.debug("{}".format(runCommand.split()))
        # self.parFileIO.updateFile()
        # subprocess.call(shlex.split(runCommand),stdout=f_stdout)
        subprocess.call(runCommand,stdout=f_stdout,stderr=f_stdout,shell=True)
        #subprocess.call(runCommand.split(),stdout=f_stdout,stderr=f_stdout)
        #subprocess.call(runCommand.split())

        # # write output to file
        os.chdir(currentDir)
        f_stdout.close()
        self.log.info('GENE finished, output written to: ' + str(self.diagDir))
        
    def getStdout(self,stdoutFileName=geneStdout):
        f = open(self.diagDir+stdoutFileName)
        s = f.read()
        f.close()
        return s

    def getNrProcessors(self,stdoutFileName=geneStdout):
        f_stdout = open(self.diagDir + stdoutFileName)

        for line in f_stdout:
            if 'MPI tasks.' in line:
                return float(line.split()[-3])
        raise IOError('Line giving the number of MPI tasks not found')

    def getAlgTime(self, stdoutFileName=geneStdout,normalize=False):
        # seems to be gone in current GENE
        f_stdout = open(self.diagDir + stdoutFileName)

        if normalize:
            n_procs = self.getNrProcessors()
            t=lambda x: x*n_procs
        else:
            t=lambda x:x

        for line in f_stdout:
            if 'Total time of' in line:
                return t(float(line.split()[-2]))
            elif 'time for eigenvalue computation:' in line:
                return t(float(line.split()[-2]))
        raise IOError('No time in file!')

    def readCheckpoint(self):
        try:
            x_local = self.parFileIO.getPar('x_local')
        except KeyError:
            x_local = True

        if x_local:
            self.cp = checkpoint2Local(self.diagDir, chptDir=self.chptDir)
        else:
            raise NotImplementedError
        
    def hasValidCheckpoint(self):
        cp = checkpoint2Base(self.diagDir)
        if cp.hasParFile:
            return self.parFileIO.isSameAs(cp.parFile)
        else:
            self.log.debug('no checkpoint data found at ' + self.diagDir)
            return False



class checkpoint2Base(object):
    '''
    checkpoint base class
    '''

#     __metaclass__=abc.ABCMeta

    def __init__(self, diagDir, chptDir=None):
        '''
        creates checkpoint from parameters.dat file
        '''
        self.log = logging.getLogger(__name__)

        self.diagDir = diagDir
        if chptDir == None:
            self.chptDir = self.diagDir
        else:
            self.chptDir = chptDir

        self.parameterFile = diagDir + commonFileNames.parameterFile
        self.path = self.chptDir + commonFileNames.checkpointFile
        self.fieldpath = diagDir + commonFileNames.fieldFile
        
        self.pathToOmega = diagDir + commonFileNames.omegaFile
        self.pathToEV = diagDir + commonFileNames.eigenvalueFile

        self.log.debug('Try to read parameter file')
        try:
            self.readParameterfile()
        except IOError as e:
            self.log.info('No parameterfile existing. Could not be read')
            self.hasParFile = False

        self.log.debug(self.parameterFile + '\t' + self.chptDir)

    def getPar(self, entry):
        try:
            buf = self.parameterDictionary[entry]
        except AttributeError:
            self.log.debug('Read parfile')
            self.readParameterfile()
            buf = self.parameterDictionary[entry]
        return num(buf)

    def setPar(self, entry, value):
        try:
            self.parameterDictionary[entry] = value
        except IOError:
            self.log.debug('Read parfile')
            self.readParameterfile()
            self.parameterDictionary[entry] = value

    def getNrProcessors(self,stdoutFileName=commonFileNames.stdOutFile):
        f_stdout = open(self.diagDir + stdoutFileName)

        for line in f_stdout:
            if 'MPI tasks.' in line:
                return float(line.split()[-3])
        raise IOError('Line giving the number of MPI tasks not found')
            
    def getAlgTime(self, stdoutFileName=commonFileNames.stdOutFile,normalize=False):
        # seems to be gone in current GENE
        f_stdout = open(self.diagDir + stdoutFileName)

        if normalize:
            n_procs = self.getNrProcessors()
            t=lambda x: x*n_procs
        else:
            t=lambda x:x

        for line in f_stdout:
            if 'Total time of' in line:
                return t(float(line.split()[-2]))
            elif 'time for eigenvalue computation:' in line:
                return t(float(line.split()[-2]))
        raise IOError('No time in file!')

    def isEVrun(self):
        if self.getPar('comp_type') == 'EV':
            return True
        else:
            return False

    def readParameterfile(self):
        self.parFile = parameterFileObject(self.parameterFile)
        self.parameterDictionary = self.parFile.getDictionaryCopy()
        self.hasParFile = True

    def getKx_offset(self):
        return self.getPar('nx0') / 2
    
    def getN(self):
        return  int(round(self.getPar('shat') * self.getPar('kymin') * self.getPar('lx'), 0))

    def storeCheckpoint(self, filename):
        f = open(filename, 'wb')
        self.log.debug('Storing checkpoint in ' + filename)
        pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def clearData(self):
        ''' might not delete all GEne Data, but the stuff which is used here'''
        for fName in self.files:
            os.remove(fName)

    def saveFiles(self, name, additionalFiles=[]):
        if not os.path.exists(name):
            self.log.debug('created checkpoint dir: ' + name)
            os.mkdir(name)
            
        files = [self.parameterFile, self.path, self.fieldpath]
        
        if self.isEVrun():
            files.append(self.pathToEV)
        else:
            files.append(self.pathToOmega)
            
        files.extend(additionalFiles)
        self.log.debug(str(files))
        
        for f in files:
#             try:
                self.log.debug('Create new files from:   ' + f + '   ' + name + '/' + f.split('/')[-1])
                shutil.copy(f, name + '/')  # + f.split('/')[-1])
                
#             except IOError:
#                 self.log.warning('IOError: something wrong with the filename')

        self.log.info('copied ' + str(files) + ' into ' + name)
        
    def createGenEnv(self):
        #TODO to implement
        raise NotImplementedError
    
#     @abc.abstractmethod
    def getData(self):
        raise NotImplementedError
        return 
    
    def getLevel(self):
        coords = ['nx0', 'nky0', 'nz0', 'nv0', 'nw0', 'n_species']
#         gridSizes = []
#         for c in coords[-1::-1]:
#             gridSizes.append(num(c))
            
#         l[5] += 1
#         l[4] -= 1
#         l[3] += 1
#         l[2] -= 1
#         l[1] += 1
        nx = self.getPar('nx0') - 1
        ny = self.getPar('nky0') + 1
        nz = self.getPar('nz0')
        nv = self.getPar('nv0')
        nw = self.getPar('nw0')


        check = self.checkIfPowerOf2
        hasLevel = check(nx) * \
                   check(ny) * \
                   check(nv) * \
                   check(nw)
                   
        self.log.debug(str(nx) + ' ' + str(ny) + ' ' + str(nz) + ' ' + str(nv) + ' ' + str(nw))
        self.log.debug(str(((np.log2(np.array([nw, nv, nz, ny, nx]) ).astype(int) ))))
                   
        if hasLevel:
            return tuple((np.log2(np.array([nw, nv, nz, ny, nx]))).astype(int))
        else:
            return None

    def checkIfPowerOf2(self, num):
        return ((num & (num - 1)) == 0) and num != 0





class checkpoint2Local(checkpoint2Base):
    
    def __init__(self, *args, **kwargs):
        checkpoint2Base.__init__(self, *args, **kwargs)
        # super(checkpoint2Base,self).__init__(*args,**kwargs)
        self.dataIsCreated = False
#         if self.hasParFile:
#             self.createAllData()
#         else:
#             self.log.warning('No parameterfile at hand: Data is not created')
        self.log.debug(self.__class__.__name__ + ' created')

    def getData(self):
        if self.isEVrun():
            return self.eigenvectors[0]
        else:
            return self.snapshots[-1]
            
    def createAllData(self):
        self.log.debug('start creating all data')
        self.readParameterfile()
        self.createData()
        #self.createFieldData()
        assert self.getPar('nonlinear') == 'F'
        #self.createEV()

    def checkSanityOfData(create):
        ##TODO should check for sizes etc.
        def newCreate(self):
            self.log.debug('See if eigenvalue stuff converged')
            stdoutFile = self.diagDir+commonFileNames.stdOutFile
            try:
                f = open(stdoutFile,'r')
                fr = f.read()
                f.close()
                if '**NO EIGENVALUES' in fr:
                    raise GeneNotConvergedError(self.diagDir, \
                                                'Eigenvalues were not computed')
                    ##hier koennen noch andere Kriterien hin
            except IOError:
                self.log.debug('gene stdout not found')
            create(self)
            # self.log.debug('end Decorator')
        return newCreate
        
    @checkSanityOfData
    def createData(self):
        if not hasattr(self, 'parameterDictionary'):
            self.readParameterfile()

        self.dataIsCreated = True
        self.testfield = 1.1
        self.f = open(self.path, 'rb')
        self.precision = self.f.read(6)
        databuffer = []
        stepbuffer = []

        parDict = self.parameterDictionary
        resParFile = (int(parDict['nx0']), int(parDict['nky0']), int(parDict['nz0']), int(parDict['nv0']), int(parDict['nw0']), int(parDict['n_spec']))

        while(1):
            readbuffer = self.f.read(8)
            if readbuffer == '':
                break
            step = np.array(struct.unpack("d", readbuffer))
            stepbuffer.append(step[0])
            
            self.dtcp = np.array(struct.unpack("d", self.f.read(8)))
            self.res = (struct.unpack("iiiiii", self.f.read(24)))
            self.log.debug('res: '+str(self.res))

            assert self.res == resParFile, 'resolution doesnt fit checkpoint' + str(self.res) + '  ' + str(resParFile)

            self.data = []
#             lengthOfData = np.multiply.accumulate(self.res)[-1]
            lengthOfData = self.res[0]*self.res[1]*self.res[2]*self.res[3]*self.res[4]*self.res[5]
            self.log.info('read array of length '+str(lengthOfData)+' from file '+self.path)

            stringData = self.f.read(16*lengthOfData)
            self.data=np.fromstring(stringData,dtype=complex)
            
            self.data = np.array(self.data)
            self.log.debug(str(self.data.dtype))
            self.log.debug(str(self.data.shape))
            self.log.debug(str((self.res[5], self.res[4], self.res[3], self.res[2], self.res[1], self.res[0])))

            self.data = self.data.reshape((self.res[5], self.res[4], self.res[3], self.res[2], self.res[1], self.res[0]))

            # rearrangement x to -3,-2,-1,0,1,2,3
            #truncInd = self.data.shape[5] / 2 + 1
            #self.data = np.concatenate((self.data[:, :, :, :, :, truncInd:], self.data[:, :, :, :, :, :truncInd]), axis=5)

            # # this should be shifted to geneCombigrid
            # append artificial boundary
            #            self.data=np.concatenate((self.data[:,:,:,0:1,:,:],self.data),axis=3)
            #self.data = np.concatenate((self.data, self.data[:, :, :, -1:, :, :]), axis=3)
            self.log.debug('data shape after rearrangement' + str(self.data.shape))
 
            N = self.getN()
            kx_offset = self.getKx_offset()
 
            # #phasefactor
            # # TODO hier stimmt das indexing nicht wen nky>1
            for j in range(self.getPar('ky0_ind'), self.data.shape[4] + self.getPar('ky0_ind')):
                for i in (np.arange(self.data.shape[5] - 1) - self.data.shape[5] / 2):
                    kx_ = (i + N * j)
#                     print j-1,i+kx_offset,self.data.shape[4],self.getPar('ky0_ind')
                    #self.data[:, :, :, -1, j - 1, i + kx_offset] = self.data[:, :, :, 0, j - 1, kx_ + kx_offset] * (-1) ** (N * j)

            databuffer.append(self.data)

        if self.isEVrun():
            self.eigenvectors = databuffer
        else:
            self.snapshots = databuffer
            self.snapshot_times = stepbuffer
        self.f.close()


    @checkSanityOfData
    def createEV(self):
        if not hasattr(self, 'parameterDictionary'):
            self.readParameterfile()
        self.log.debug('start reading ev file'+str(self.isEVrun()))
        if self.isEVrun():
            self.log.debug("which_ev : {}".format(self.parFile.getPar('which_ev')))
            if self.parFile.getPar('which_ev') not in  ['rhs_only','solve_system']:
                self.log.debug('read ev file '+self.pathToEV)
                f = open(self.pathToEV, 'r')
                buffer = f.readline()
                buffer = []
                for line in f:
                    buffer.append(complex(float(line.split()[0]), float(line.split()[1])))
                self.log.debug(str(buffer))
            else:
                buffer = None
        else:
            f = open(self.pathToOmega, 'r')
            buffer = []
            for line in f:
    #            print line
                buffer.append(complex(float(line.split()[1]), float(line.split()[2])))
        self.EV = buffer

    @checkSanityOfData
    def createFieldData(self):
        if not hasattr(self, 'parameterDictionary'):
            self.readParameterfile()
        f = open(self.fieldpath)
        self.fieldData = []

        while(1):
            readbuffer = f.read(4)
            if readbuffer == '':
                break
            bytes = struct.unpack("i", readbuffer)[0]
            readbuffer = f.read(bytes)
            t = struct.unpack("d", readbuffer)[0]
            readbuffer = f.read(4)
            field_t = []

            N = self.getN()
            kx_offset = self.getKx_offset()

            for i in range(self.getPar('n_fields')):
                fieldbuffer = []
                bytes = struct.unpack("i", f.read(4))[0]
                for j in range(bytes / 16):
                    readbuffer = struct.unpack("dd", f.read(16))
                    fieldbuffer.append(complex(readbuffer[0], readbuffer[1]))
                if bytes != struct.unpack("i", f.read(4))[0]:
                    raise IOError('struct size does not fit')

                fieldbuffer = np.array(fieldbuffer).reshape([self.getPar('nz0'), self.getPar('nky0'), self.getPar('nx0')])
                # ##rearrange
                truncInd = fieldbuffer.shape[2] / 2 + 1
                fieldbuffer = np.concatenate((fieldbuffer[:, :, truncInd:], fieldbuffer[:, :, :truncInd]), axis=2)

                # ##put in artificial boudnary
                fieldbuffer = np.concatenate((fieldbuffer, fieldbuffer[-1:, :, :]), axis=0)
                # ##put in boundary stuff fluxtube region ( for combination)

                for j in range(self.getPar('ky0_ind'), self.getPar('nky0') + self.getPar('ky0_ind')):
                    for i in (np.arange(self.getPar('nx0') - 1) - self.getPar('nx0') / 2):
                        kx_ = (i + N * j)

                        fieldbuffer[-1, j - 1, i + kx_offset] = fieldbuffer[0, j - 1, kx_ + kx_offset] * (-1) ** (N * j)

                field_t.append(fieldbuffer)

            self.fieldData.append(field_t)
        


class GeneNotConvergedError(Exception):
    def __init__(self,diagDir,msg):
        self.diagDir = diagDir
        self.msg = msg



def storeInCheckpoint(extrapolated_data,pathToCheckpoint,chpt_name='checkpoint_extra'):
        """
        here data is stored in a checkpointfile and the full array is cut appropriately
        :param extrapolated_data:
        :param pathToCheckpoint:
        :param chpt_name:
        """
        extrapolated_data = extrapolated_data[:, :-1, :-1, :-1, 1:-1, :]

        # #open file
        # try:
        #os.mkdir(pathToCheckpoint)
        if not os.path.exists(pathToCheckpoint):
            os.makedirs(pathToCheckpoint)

        # except OSError:
        # self.log.debug('Directory already existing')

        filename = pathToCheckpoint + chpt_name
        f = open(filename, 'wb')
        # #datatype
        f.write('DOUBLE')
        # # current count of timestep
        f.write(struct.pack('d', 1.0))
        # #write current timestep
        f.write(struct.pack('d', 1.0))
        # #write resolution
        f.write(struct.pack("iiiiii", *(list(extrapolated_data.shape)[::-1])))

        truncInd = extrapolated_data.shape[5] / 2
        extrapolated_data = np.concatenate((extrapolated_data[..., truncInd:], extrapolated_data[..., :truncInd]),
                                              axis=5)

        extrapolated_data = extrapolated_data.flatten()

        # # write data
        f.write(extrapolated_data.tostring())
        f.close()
