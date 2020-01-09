"""Gridvx class constructs an optimal structure of the v||, x grid."""

from grid import GridBase
from math import sqrt, pi, fabs
from scipy.special import erf
import numpy as np
import profiles

class Gridvx(GridBase):

    def __init__(self, profiles, precision_pair):
        super(Gridvx, self).__init__(profiles, precision_pair)
        # compute scaled standard deviations for all profile points
        self._parameter_array = np.sqrt(
            np.divide(
                np.multiply(profiles.profile_array[:,2],2),
                profiles.ref_temperature
                )
            )

    def getDPrecision(self, sigma_down, sigma_current):
        
        """return difference between profile and block structure precision."""
        
        n_new = self._precision_pair[1]*sigma_down/sigma_current
        return erf(n_new/sqrt(2)) - self._precision_pair[0]
    
    def getDPrecisionLambda(self, sigma_parent):

        """return function that returns precision difference \
        for a fixed sigma_parent value and given current sigma"""

        return lambda sigma_current: self.getDPrecision(sigma_parent, sigma_current)

    def getDPrecisionRadialLambda(self, sigma_parent):
        
        """return function that returns precision difference \
        for a fixed sigma_parent value and given radial position"""
        return lambda x: self.getDPrecision(
            sigma_parent,
            sqrt(2*self._profiles.getTemperature(x)/self._profiles.ref_temperature)
            )

    def getParameter(self, mark):
        """return normal distribution funtion \
        standard deviation at the mark position """
        return sqrt(2*self._profiles.getTemperature(mark)/ \
                         self._profiles.ref_temperature)

    def getBGridVizStructure(self, blk_mks):

        """ returns the structure of block grid for vizualization purpose """

        blk_mks_r = blk_mks[1:-1]
        if blk_mks.size > 1: blk_mks_r.shape = (blk_mks_r.size, -1)
        blk_mks_r = np.hstack((blk_mks_r, blk_mks_r))
        blk_mks_r = blk_mks_r.ravel()
        blk_mks_r = np.hstack((self._profiles.xval_a_min, blk_mks_r))
        blk_mks_r = np.hstack((blk_mks_r, self._profiles.xval_a_max))
        blk_mks_r = np.hstack((blk_mks_r, blk_mks_r[::-1]))

        blk_mks_v = blk_mks[:-1]
        if blk_mks.size > 1: blk_mks_v.shape = (blk_mks_v.size, -1)
        blk_mks_v = np.hstack((blk_mks_v, blk_mks_v))
        blk_mks_v = blk_mks_v.ravel()
        
        blk_vel_range = self._profiles.getTemperature(blk_mks_v)
        blk_vel_range = np.sqrt(np.divide(np.multiply(blk_vel_range, 2), \
                                                self._profiles.ref_temperature))
        blk_vel_range = np.multiply(blk_vel_range, self._precision_pair[1])
        blk_vel_range = np.hstack((blk_vel_range, -blk_vel_range[::-1]))

        return np.vstack((blk_mks_r, blk_vel_range))
    
    def getDistrSurf(self, x, y):

        """return surface points of the distribution function on the mesh \
        represented by (x,y) points."""

        temp = self._profiles.getTemperature(y)
        temp = np.divide(temp, self._profiles.ref_temperature)
        velp2 = x**2
        F    = np.exp(-np.divide(velp2,temp))
        temp = np.sqrt(np.multiply(pi,temp))
        temp = np.divide(1, temp)
        return np.multiply(temp, F)
        
        
    def getVelocityBlkMks(self, blk_mks_r):
        """returns velocity range for each block"""
        blk_mks_v = self._profiles.getTemperature(blk_mks_r)
        blk_mks_v = np.sqrt(np.divide(np.multiply(blk_mks_v, 2), \
                                                self._profiles.ref_temperature))
        blk_mks_v = np.multiply(blk_mks_v, self._precision_pair[1])
        return blk_mks_v

    @property
    def sigma(self):
        """'_parameter_array' property."""
        return self._parameter_array

    ################### Transfinite Interpolation Grid ############################
    def getGridTfI(self, m, n):
        
        """ generates transfinite interpolation grid and return it in x, y arrays """

        # TODO:
        # 1. choose discretization (m, n)
        # 2. compute boundary arrays xb, yb, xt, yt, xl, yl, xr, yr
        # 3. transfer boundary arrays to grid arrays x and y
        # 4. transfinite interpolation on eight boundary arrays

        # TODO (figure this out, should be GL points):
        # start with: equidistant grid in radial direction

        dx = 1./(m - 1)
        dy = 1./(n - 1)
    
        radius = np.linspace(self._profiles.xval_a_min, self._profiles.xval_a_max, n)
        temper = self._profiles.getTemperature(radius)
        sigma = np.sqrt(np.divide(np.multiply(temper, 2), self._profiles.ref_temperature))
        vel_range = np.multiply(sigma, self._precision_pair[1])
        
        # bottom boundary
        xb = np.linspace(-vel_range[0], vel_range[0], m)
        yb = np.linspace(radius[0], radius[0], m)
        
        # top boundary
        xt = np.linspace(-vel_range[-1], vel_range[-1], m)
        yt = np.linspace(radius[-1], radius[-1], m)
        
        # left boundary
        xl = -vel_range
        yl = radius
        
        # right boundary
        xr = vel_range
        yr = radius
        
        x = np.empty((m, n))
        y = np.empty((m, n))

        # transfers the eight arrays containing 1D boundary data to
        # two 2D arrays (x,y)

        x[:, 0] = xb
        y[:, 0] = yb

        x[:, -1] = xt
        y[:, -1] = yt
        
        x[0, :] = xl
        y[0, :] = yl
        
        x[-1, :] = xr
        y[-1, :] = yr

        # performs 2D transfinite interpolation on the eight boundary
        # arrays to obtain the grid (x,y)
        
        for j in range(1, n - 1):
            s = j*dy
            for i in range(1, m - 1):
                r = i*dx
                
                x[i, j] = (1. - s)*xb[i] + s*xt[i] + (1. - r)*xl[j] \
                    + r*xr[j] - (r*s*xr[-1] + (1. - r)*s*xl[-1] \
                    + r*(1. - s)*xr[0] + (1. - r)*(1. - s)*xl[0])
                
                y[i, j] = (1. - s)*yb[i] + s*yt[i] + (1. - r)*yl[j] \
                    + r*yr[j] - (r*s*yr[-1] + (1. - r)*s*yl[-1] \
                    + r*(1. - s)*yr[0] + (1. - r)*(1. - s)*yl[0]) 
                
        return (x, y)

    ##############################################################################

    ################### Postprocessing -> Estimation of Grid #####################
    def getBlkTicks(self, blk_mks_r_v, blk_mks_v):

        """ returns coordinates lines ticks for each block. """

        dx = self._profiles.dx
        iblk_mks_r_v = np.add(
            np.round_(
                np.divide(
                    np.subtract(blk_mks_r_v,blk_mks_r_v[0]), dx
                    )
                ),
            1)
        iblk_mks_r_v[0] = 1

        dv = 2*blk_mks_v[0]/(self._profiles.nv0 - 1)
        iblk_mks_v = np.add(np.round_(np.divide(np.multiply(2,blk_mks_v),dv)),1)
        iblk_mks_v = np.add(
            iblk_mks_v,
            np.fabs(np.subtract(self._profiles.nv0%2,np.mod(iblk_mks_v,2)))
            )
        iblk_mks_v[-1] = 0

        n = iblk_mks_r_v.size-1
        nx_blk = np.empty(n)
        nv_blk = np.empty(n)
        x_blks = []
        v_blks = []

        for i in range(0, n):
            nx_blk[i] = iblk_mks_r_v[i+1] - iblk_mks_r_v[i] + 1
            nv_blk[i] = iblk_mks_v[i]
            x_blk = np.arange(0,nx_blk[i])
            x_blk = np.multiply(x_blk,dx)
            x_blk = np.add(x_blk,blk_mks_r_v[0] + (iblk_mks_r_v[i]-1)*dx)
            x_blks.append(x_blk)
            v_blk = np.arange(0,nv_blk[i])
            v_blk = np.multiply(v_blk,dv)
            v_blk = np.subtract(v_blk,(nv_blk[i]-1)*dv/2)
            v_blks.append(v_blk)

        return (x_blks, v_blks)

    def getBlkTicks_badptv(self, blk_mks_r_v, blk_mks_v):

        """ returns coordinates lines ticks for each block. """

        dx = self._profiles.dx
        iblk_mks_r_v = np.add(
            np.round_(
                np.divide(
                    np.subtract(blk_mks_r_v,blk_mks_r_v[0]), dx
                    )
                ),
            1)
        iblk_mks_r_v[0] = 1

        dv = 2*blk_mks_v[0]/(self._profiles.nv0 - 1)
        iblk_mks_v = np.add(np.round_(np.divide(np.multiply(2,blk_mks_v),dv)),1)
        iblk_mks_v = np.add(
            iblk_mks_v,
            np.fabs(np.subtract(self._profiles.nv0%2,np.mod(iblk_mks_v,2)))
            )
        iblk_mks_v[-1] = 0

        n = iblk_mks_r_v.size-1
        nx_blk = np.empty(n)
        x_blks = []
        v_blks = []

        print "number of points in each v|| block: ", iblk_mks_v

        # standard -> normalized number of points
        nv_blk_std = iblk_mks_v[-2]
        for i in range(0, n):
            nx_blk[i] = iblk_mks_r_v[i+1] - iblk_mks_r_v[i] + 1
            dv = 2*blk_mks_v[i]/(nv_blk_std - 1)
            x_blk = np.arange(0,nx_blk[i])
            x_blk = np.multiply(x_blk,dx)
            x_blk = np.add(x_blk,blk_mks_r_v[0] + (iblk_mks_r_v[i]-1)*dx)
            x_blks.append(x_blk)
            v_blk = np.arange(0,nv_blk_std)
            v_blk = np.multiply(v_blk,dv)
            v_blk = np.subtract(v_blk,(nv_blk_std-1)*dv/2)
            v_blks.append(v_blk)

        return (x_blks, v_blks)

    ##############################################################################
