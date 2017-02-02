"""Gridwx class constructs an optimal structure of the mu, x grid."""

from grid import GridBase
from math import log, exp
from gquad import gaulag
import numpy as np
import profiles

class Gridwx(GridBase):

    def __init__(self, profiles, precision_pair):
        super(Gridwx, self).__init__(profiles, precision_pair)
        # compute mu on the counter of the given tolerance for all profile points
        self._parameter_array = np.divide(
            np.multiply(
                precision_pair[1],
                profiles.profile_array[:,2]
                ),
            profiles.ref_temperature
            )
        self._nx_blks = None
        self._nw_blks = None
        self._blk_mks_r_w = None
        self._blk_mks_w = None
        self._iblk_mks_r_w = None
        self._iblk_mks_w = None
        self._n_wx_blks = None
        self._mu = None
        
    def getDPrecisionRadialLambda(self, mu_down):
        """return function that returns precision difference \
        for a fixed mu_parent value and given radial position"""
        # return lambda x: 1 - exp(-mu_down*self._profiles.ref_temperature/ \
        #                          self._profiles.getTemperature(x)) - \
        #                          self._precision_pair[0]
        return lambda x: mu_down - self._precision_pair[1]*self._profiles.getTemperature(x)/ \
               self._profiles.ref_temperature

    def getParameter(self, x):
        """ return mu at the mark position \
        that corresponds to the tolerance value"""
        return self._precision_pair[1]*self._profiles.getTemperature(x)/ \
               self._profiles.ref_temperature

    def getDPrecision(self, mu_down, mu_current):

        """return difference between profile and block structure precision."""

        return 1 - exp(-self._precision_pair[1]*mu_down/mu_current) \
               - self._precision_pair[0]

    def getMagneticMomentBlkMks(self, blk_mks_r):
        """returns magnetic moment maximum value for each block"""
        return self.getParameter(blk_mks_r)


    @property
    def mu(self):
        """'_parameter_array' property."""
        return self._parameter_array

    
    def getBGridVizStructure(self, blk_mks_r, blk_mks_w):

        """ returns the structure of block grid for vizualization purpose """

        sblk_mks_r = blk_mks_r[1:]
        if blk_mks_r.size > 1: sblk_mks_r.shape = (sblk_mks_r.size, -1)
        sblk_mks_r = np.hstack((sblk_mks_r, sblk_mks_r))
        sblk_mks_r = sblk_mks_r.ravel()
        sblk_mks_r = np.hstack((blk_mks_r[0], sblk_mks_r))

        sblk_mks_w = blk_mks_w[:-1]
        if blk_mks_w.size > 1: sblk_mks_w.shape = (sblk_mks_w.size, -1)
        sblk_mks_w = np.hstack((sblk_mks_w, sblk_mks_w))
        sblk_mks_w = sblk_mks_w.ravel()
        sblk_mks_w = np.hstack((sblk_mks_w, blk_mks_w[-1]))
        
        return np.vstack((sblk_mks_r, sblk_mks_w))

    def getDistrSurf(self, x, y):

        """return surface points of the distribution function on the mesh \
        represented by (x,y) points."""

        temp = self._profiles.getTemperature(y)
        temp = np.divide(temp, self._profiles.ref_temperature)
        return np.divide(np.exp(-np.divide(x,temp)),temp)
        
    ################### Postprocessing -> Estimation of Grid #####################

    def optimizeBlkStr(self, blk_mks_r_w, blk_mks_w):

        """ optimizes block structure of the x, w subgrid """

        dx = self._profiles.dx
        iblk_mks_r_w = np.add(
            np.round_(
                np.divide(
                    np.subtract(blk_mks_r_w,blk_mks_r_w[0]), dx
                    )
                ),
            1)
        iblk_mks_r_w[0] = 1

        n_wx_blks = blk_mks_r_w.size - 1 
        iblk_mks_w = np.empty(n_wx_blks + 1)
        nw0 = int(self._profiles.nw0)
        (self._mu, wmu) = gaulag(blk_mks_w[0],nw0)
        for ib in range(0, n_wx_blks):
            for i in range(0, nw0):
                if ((self._mu[i] >= blk_mks_w[ib]) or (i == nw0-1)):
                    iblk_mks_w[ib] = i + 1
                    break
        iblk_mks_w[-1] = 0
        
        # merge blocks
        idel = []
        # for i in range(1, n_wx_blks):
        #     if(iblk_mks_w[i] == iblk_mks_w[i-1]):
        #         idel.append(i)

        self._blk_mks_r_w  = np.delete(blk_mks_r_w, idel)
        self._blk_mks_w    = np.delete(blk_mks_w, idel)
        self._iblk_mks_r_w = np.delete(iblk_mks_r_w, idel)
        self._iblk_mks_w   = np.delete(iblk_mks_w, idel)
        self._n_wx_blks    = self._blk_mks_r_w.size - 1

        self._nx_blks = np.empty(self._n_wx_blks)
        self._nw_blks = np.empty(self._n_wx_blks)
        for i in range(0, self._n_wx_blks):
            self._nx_blks[i] = self._iblk_mks_r_w[i+1] - self._iblk_mks_r_w[i] + 1
            self._nw_blks[i] = self._iblk_mks_w[i]

            
    def getBlkTicks(self):

        """ returns coordinates lines ticks for each block. """

        if ((self._nx_blks == None) or (self._nw_blks == None) \
            or (self._blk_mks_w == None) or (self._mu == None)):
            print "optimize block structured grid before getting block ticks!"
            return (None, None)

        dx = self._profiles.dx
        x_blks = []
        w_blks = []
        for i in range(0, self._n_wx_blks):
            x_blk = np.arange(0,self._nx_blks[i])
            x_blk = np.multiply(x_blk,dx)
            x_blk = np.add(x_blk,self._blk_mks_r_w[0] + \
                           (self._iblk_mks_r_w[i]-1)*dx)
            x_blks.append(x_blk)
            w_blk = np.array(self._mu[0:int(self._nw_blks[i])])
            w_blks.append(w_blk)
        
        return (x_blks, w_blks)

    @property
    def blk_mks_r_w(self):
        """'blk_mks_r_w' property"""
        return self._blk_mks_r_w

    @property
    def blk_mks_w(self):
        """'blk_mks_w' property"""
        return self._blk_mks_w

    @property
    def n_wx_blks(self):
        """'n_wx_blks' property"""
        return self._n_wx_blks

    ##############################################################################
