"""GridBase is a base class for constractions of an optimal structures of the grid"""

import numpy as np
from scipy.optimize import fmin_powell
from scipy.integrate import quad
from math import ceil, fabs, pi
import profiles
import abc

class GridBase(object):
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, profiles, precision_pair):
        self._profiles        = profiles
        self._precision_pair  = precision_pair
        self._parameter_array = None

    @abc.abstractmethod
    def getDPrecision(self, parameter_parent, parameter_current):        
        """return difference between profile and block structure precision."""
        return

    @abc.abstractmethod
    def getDPrecisionRadialLambda(self, parameter):
        """return function that returns precision difference \
        for a fixed parameter value and given radial position"""
        return

    @abc.abstractmethod
    def getParameter(self, mark):
        """ return a parameter that is used to generate \
        a precision difference function"""
        return

    @abc.abstractmethod
    def getDPrecision(self, parameter_down, parameter_up):
        """retrun difference between ideal profile and \
        block structure precision."""
        return

    def getProbIntOutDomain(self, blk_mks):

        """return probability difference between the smooth and blocked regions"""
        
        pdiff = 0.0
        penalty = 1000.0
        # punish for being out the range by returning a reasonably big number
        for blk_mark in blk_mks:
            if (blk_mark < self._profiles.xval_a_min) or \
                   (blk_mark > self._profiles.xval_a_max):
                pdiff += penalty
        if (pdiff > 0.0): return pdiff

        mark_down = self._profiles.xval_a_min
        if blk_mks.size == 0:
            mrange = np.array([self._profiles.xval_a_max])
        else:
            mrange = np.hstack((blk_mks, self._profiles.xval_a_max))

        for mark_up in mrange:
            
            # get parameter (can be mu, v||, standard deviation, ...) from the mark
            parameter_down = self.getParameter(mark_down)
            
            # generate function for integration
            getDPrecision = self.getDPrecisionRadialLambda(parameter_down)

            # integrate precision diff for region between mark_down and mark_up
            pdiffinc, err = quad(getDPrecision, mark_down, mark_up)

            # increment global diff counter
            pdiff += fabs(pdiffinc)

            # move mark_down
            mark_down = mark_up

        return pdiff
    
    def getProbSumOutDomain(self, blk_mks):

        """return minimization parameter value. """

        badness = 0.0
        penalty = 10.0
        
        blk_mks_sorted = np.sort(blk_mks)

        max_index = blk_mks_sorted.size - 1
        max_mark = blk_mks_sorted[max_index]

        # punish for being larger than max_mark
        while max_mark > self._profiles.psize-1 :
            badness += penalty
            blk_mks_sorted[max_index] = self._profiles.psize-1
            max_index -= 1
            max_mark = blk_mks_sorted[max_index]

        min_index = 0
        min_mark = blk_mks_sorted[min_index]

        # punish for being smaller than min_mark
        while min_mark < 0:
            badness += penalty
            blk_mks_sorted[min_index] = 0
            min_index += 1
            min_mark = blk_mks_sorted[min_index]

        blk_mks_sorted = np.hstack((blk_mks_sorted, self._profiles.psize-1))

        mark_down = 0
        for mark_up in blk_mks_sorted:
            
            # get parameter (can be mu, v||, standard deviation, ...) from the mark
            parameter_down = self._parameter_array[ceil(mark_down)]

            for mark_current in range(int(ceil(mark_down))+1, \
                                          int(ceil(mark_up))+1):
                parameter_current = self._parameter_array[mark_current]
                
                badness += fabs(self.getDPrecision(
                    parameter_down,parameter_current))

                # badness += fabs(parameter_down - parameter_current)
                
            # move mark_down
            mark_down = mark_up

        return badness

    def getOptBGridSmoothStr(self, num_of_blks):

        """ method finds and optimal block grid structure by minimizing the
            probability to find particles in the region between the block and
            ideal smooth contuors """

        # initialize block division marks
        blk_mks = np.linspace(self._profiles.xval_a_min, self._profiles.xval_a_max, \
                                      num_of_blks + 1)[1:-1]

        # @TODO: remove after testing
        print "initial block division marks: " + str(blk_mks) + \
            " with minimization parameter: " + str(self.getProbIntOutDomain(blk_mks))

        # optimizing block marks positions
        res = fmin_powell(self.getProbIntOutDomain, blk_mks, \
                              xtol = 1e-8, ftol = 1e-6,  disp = 1)

        if num_of_blks > 2 :  blk_mks = np.array(res)
        else : blk_mks = np.array([res])
        
        # @TODO: remove after testing
        print "final division indices: " + str(blk_mks) + \
            " with minimization parameter: " + str(self.getProbIntOutDomain(blk_mks))

        if blk_mks.size == 0:
            blk_mks = np.hstack((self._profiles.xval_a_min,
                                     self._profiles.xval_a_max))
        else:
            blk_mks = np.hstack((self._profiles.xval_a_min,
                                     blk_mks,
                                     self._profiles.xval_a_max))

        return blk_mks
        

    def getOptBGridDiscreteStr(self, num_of_blks):
        
        """ method finds an optimal block grid structure \
        based on getProbSumOutDomain return value"""
        
        # initialization of block division marks
        blk_mks = np.array(np.int_(np.linspace(0, self._profiles.psize,\
                                      num_of_blks + 1)[1:-1]))
        if num_of_blks == 2: blk_mks = np.array([blk_mks])

        # @TODO remove after testing
        print "initial block division marks: " + str(blk_mks) + \
            " with minimization parameter: " + str(self.getProbSumOutDomain(blk_mks))

        # optimizing block marks positions
        res = fmin_powell(self.getProbSumOutDomain, blk_mks, \
                              xtol = 1e-8, ftol = 1e-6, disp = 1)

        if num_of_blks > 2 :  blk_mks = np.array(res)
        else : blk_mks = np.array([res])

        smooth_blk_mks = \
            self._profiles.profile_array[:,0][np.int_(np.ceil(np.sort(blk_mks)))]

        # @TODO: remove after testing
        print "final division indices: " + str(smooth_blk_mks) + \
            " with minimization parameter: " + str(self.getProbSumOutDomain(blk_mks))

        smooth_blk_mks = \
            np.hstack((self._profiles.xval_a_min,
                       smooth_blk_mks,
                       self._profiles.xval_a_max))

        return smooth_blk_mks
