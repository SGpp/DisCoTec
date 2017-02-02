# -*- coding: utf-8 -*-
"""FSAmomData.py

Contains the infrastructure to handle fsamom_spec_fext files

"""
# pylint: disable=E1101
import sys
import copy
import numpy as np
from bisect import bisect_left
from ParIO import Parameters


class FSAmomFileData(object):
    """Class to handle a specific fsamom file """

    def __init__(self, filename, nx0, starttime, endtime):
        self.fsacolumns = 13  # Number of columns in fsamom file
        self.filename = filename
        self.nx0 = nx0      # TODO: Make common module
        self.timefld = []
        self.starttime = starttime
        self.endtime = endtime
        self.fsadata = []
        self.readfsamom()

    def get_minmaxtime(self):
        return self.timefld[0], self.timefld[-1]

    def readfsamom(self):
        """ Read flux surface averaged moments """
        blocks = []
        print('Reading  %s\n'%(self.filename))
        try:
            with open(self.filename) as fsafile:
                next(fsafile)  # skip header
                for line in fsafile:
                    if not line or line.startswith('\n'):
                        continue
                    if line.startswith('#'):
                        self.timefld.append(float(line.split()[1]))
                        blocks.append([])
                    else:
                        blocks[-1].append(line)
        except IOError:
            sys.exit("IOError: probably fsa file does not exist: %s"%self.filename)
        self.timefld = np.array(self.timefld)
        first_time, last_time = self.get_minmaxtime()
        if len(self.timefld) != len(set(self.timefld)):
            print("Error: %s contains 2 blocks with identical timestamp"%(self.filename))
            return
        if self.starttime == -1 or (self.starttime > 0 and
                                    self.starttime < first_time):
            print("Using first time present in fsa moments data for starttime")
            self.starttime = first_time
        if self.endtime == -1:
            print("Using first time present in fsa moments data for endtime")
            self.endtime = first_time
        if self.starttime == -2:
            print("Using last time present in fsa moments data for starttime")
            self.starttime = last_time
        if (self.endtime == -2) or (self.endtime > last_time):
            print("Using last time present in fsa moments data for endtime")
            self.endtime = last_time
        if (self.endtime < first_time) or (self.starttime > last_time):
            print("Time window not contained in fsa moments data")
            return
        print(('starttime=%f, endtime=%f, first_time=%f, last_time=%f' %
              (self.starttime, self.endtime, first_time, last_time)))
        #  For a single time, find element closest to given input time
        if self.starttime == self.endtime:
            pos = np.array([bisect_left(self.timefld, self.starttime)])
        else:
            pos = np.where((self.timefld >= float(self.starttime)) &
                           (self.timefld <= self.endtime))[0]
        self._addfsadata(blocks, pos)
        self.timefld = self.timefld[pos]

    def _addfsadata(self, blocks, pos):
        for entry in pos:
            self.fsadata.append(np.empty((self.nx0, self.fsacolumns)))
            for i in range(0, self.nx0):
                line = blocks[entry][i]
                for var in range(self.fsacolumns):
                    self.fsadata[-1][i, var] = float(line.split()[var])


class FSAmomData(object):
    def __init__(self, fext, starttime, endtime, plotset):
        self.isDataPresent = False
        self.fext = fext
        self.starttime = starttime
        self.endtime = endtime
        # each dataset has its own collection of desired plots
        self.plotset = copy.deepcopy(plotset)
        self.fsaspec = []
        self.readpar()

    def readpar(self):
        """ Get some essential parameters from par file"""
        par = Parameters()
        par.Read_Pars('parameters%s'%(self.fext))
        self.nspec = int(par.pardict['n_spec'])
        self.nx0 = int(par.pardict['nx0'])
        self.istep_fsamom = 0
        if 'plotfsa' in self.plotset and 'istep_fsa_moments' in par.pardict:
            self.istep_fsamom = int(par.pardict['istep_fsa_moments'])
        self.specname = [[], []]
        for n in range(self.nspec):
            self.specname[n] = str(par.pardict['name%1d' % (n+1)])
            self.specname[n] = self.specname[n][1:-1]

    def getfsadata(self):
        """Read the data from the fsa_moments diagnostic"""
        if self.istep_fsamom == 0:
            print("No fsamom diagnostic possible/desired")
            return
        for n in range(0, self.nspec):
            filename = 'fsamom_%s%s'%(self.specname[n], self.fext)
            self.fsaspec.append(FSAmomFileData(filename, self.nx0, self.starttime, self.endtime))
        self.isDataPresent = True
        
