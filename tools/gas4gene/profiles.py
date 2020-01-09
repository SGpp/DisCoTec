"""Profiles class keeps necessary data about density and temperature profiles"""

import numpy as np
from scipy.interpolate import interp1d
import math

class ProfilesData(object):
    
    def __init__(self, profile_dict, profile_array = None):
        
        self._prof_type = profile_dict['prof_type']
        self._kappa_T = profile_dict['kappa_T']
        self._LT_center = profile_dict['LT_center']
        self._LT_width = profile_dict['LT_width']
        self._kappa_n = profile_dict['kappa_n']
        self._Ln_center = profile_dict['Ln_center']
        self._Ln_width = profile_dict['Ln_width']
        self._minor_r = profile_dict['minor_r']
        self._major_R = profile_dict['major_R']
        self._xval_a_min = profile_dict['xval_a_min']
        self._xval_a_max = profile_dict['xval_a_max']
        self._x0 = profile_dict['x0']
        self._dx = profile_dict['dx']
        self._nx0 = profile_dict['nx0']
        self._nky0 = profile_dict['nky0']
        self._nz0 = profile_dict['nz0']
        self._nv0 = profile_dict['nv0']
        self._nw0 = profile_dict['nw0']
        self._n_spec = profile_dict['n_spec']
        
        if(self._prof_type == 3):
            self._delta_x_T = profile_dict['delta_x_T']
            self._delta_x_n = profile_dict['delta_x_n']

        # consider case of experimental profiles (data has to be interpolated)
        if(profile_array != None and self._prof_type == -1):
            # generate interpolation functions for temperatures and densties
            self.getTemperatureTM1 = interp1d(profile_array[:,0],  \
                                              profile_array[:,2], \
                                              kind = 'cubic')
            self.getDensityTM1 = interp1d(profile_array[:,1], \
                                          profile_array[:,3], \
                                          kind = 'cubic')
        else:
            self.getTemperatureTM1 = None
            self.getDensityTM1 = None

        # set dictionary of temperature function for different profiles
        self._temperature_dict = { -1 : self.getTemperatureTM1, \
                                   2 : self.getTemperatureT2,  \
                                   3 : self.getTemperatureT3}

        # set dictionary of density functions for different profiles
        self._density_dict = { -1 : self.getDensityTM1, \
                               2 : self.getDensityT2, \
                               3 : self.getDensityT3}

        # set profile array
        self.setProfileArray(profile_array)

        # generate a function that return radius from a given temperature
        self.getInvTemperature = interp1d(self._profile_array[::-1,2],  \
                                          self._profile_array[::-1,0], \
                                          bounds_error = False, \
                                          fill_value = self._xval_a_max, \
                                          kind = 'cubic')

        # set reference temperature
        self._ref_temperature = self.getTemperature(self._x0)

        # safe profile size
        self._profile_size = self.profile_array[:,0].size

    def setProfileArray(self, profile_array):
        """setup the profile array according to the given x coordinate range"""
        if(profile_array == None):
            # fake profile array for consistency!
            radius = np.linspace(self._xval_a_min, self._xval_a_max, 32)
            radius.shape = (radius.size, -1)
            temper = self.getTemperature(radius)
            densit = self.getDensity(radius)
            self._profile_array = np.hstack((radius, radius))
            self._profile_array = np.hstack((self._profile_array, temper))
            self._profile_array = np.hstack((self._profile_array, densit))
        else:
            min_i = 0
            max_i = len(profile_array[:, 0]) - 1
            for i in range(min_i+1, max_i+1):
                if profile_array[i,0] > self._xval_a_min :
                    min_i = i
                    break
            for i in range(min_i+1, max_i+1):
                if profile_array[i,0] > self._xval_a_max:
                    max_i = i
                    break
            radius_t = np.zeros(max_i - min_i + 2)
            radius_t[1:-1] = profile_array[min_i:max_i,0]
            radius_t[0] = self._xval_a_min
            radius_t[-1] = self._xval_a_max
            radius_d = np.zeros(max_i - min_i + 2)
            radius_d[1:-1] = profile_array[min_i:max_i,1]
            radius_d[0] = self._xval_a_min
            radius_d[-1] = self._xval_a_max
            temper = np.zeros(max_i - min_i + 2)
            temper[1:-1] = profile_array[min_i:max_i,2]
            temper[0] = self.getTemperature(self._xval_a_min)
            temper[-1] = self.getTemperature(self._xval_a_max)
            densit = np.zeros(max_i - min_i + 2)
            densit[1:-1] = profile_array[min_i:max_i,3]
            densit[0] = self.getDensity(self._xval_a_min)
            densit[-1] = self.getDensity(self._xval_a_max)
            
            radius_t.shape = (radius_t.size, -1)
            radius_d.shape = (radius_t.size, -1)
            temper.shape = (temper.size, -1)
            densit.shape = (densit.size, -1)
            self._profile_array = np.hstack((radius_t, radius_d))
            self._profile_array = np.hstack((self._profile_array, temper))
            self._profile_array = np.hstack((self._profile_array, densit))
            

    def getTemperature(self, x):
        return self._temperature_dict[self._prof_type](x)


    def getDensity(self, x):
        return self._density_dict[self._prof_type](x)

    
    def getTemperatureT2(self, x):
        """ return temperature at point x for profile type two """
        return np.exp(-self._kappa_T*self._minor_r*self._LT_width* \
                        np.tanh((x - self._LT_center)/self._LT_width))

    def getDensityT2(self, x):
        """ return density at point x for profile type two """
        return np.exp(-self._kappa_n*self._minor_r*self._Ln_width* \
                            np.tanh((x - self._Ln_center)/self._Ln_width))
    
    def getTemperatureT3(self, x):
        """ return temperature at point x for profile type three """
        return np.power( \
                np.cosh(( \
                        x - self._LT_center + self._delta_x_T)/self._LT_width)/ \
                np.cosh(( \
                        x - self._LT_center - self._delta_x_T)/self._LT_width), \
                -0.5*self._kappa_T*self._minor_r*self._LT_width)

    def getDensityT3(self, x):
        """ return density at point x for profile type three """
        return np.power( \
                np.cosh(( \
                        x - self._Ln_center + self._delta_x_n)/self._Ln_width)/ \
                np.cosh(( \
                        x - self._Ln_center - self._delta_x_n)/self._Ln_width), \
                -0.5*self._kappa_n*self._minor_r*self._Ln_width)

    # decarators that serve as getters in the pythonic way

    @property
    def ref_temperature(self):
        """'ref_temperature' property."""
        return self._ref_temperature

    @property
    def xval_a_min(self):
        """'xval_a_min' property."""
        return self._xval_a_min

    @property
    def xval_a_max(self):
        """'xval_a_max' property."""
        return self._xval_a_max

    @property
    def sigma(self):
        """'sigma' property."""
        return self._sigma

    @property
    def dx(self):
        """'dx' property."""
        return self._dx

    @property
    def nx0(self):
        """'nx0' property."""
        return self._nx0

    @property
    def nky0(self):
        """'nky0' property."""
        return self._nky0

    @property
    def nz0(self):
        """'nz0' property."""
        return self._nz0

    @property
    def nv0(self):
        """'nv0' property."""
        return self._nv0

    @property
    def nw0(self):
        """'nw0' property."""
        return self._nw0

    @property
    def n_spec(self):
        """'nx0' property."""
        return self._n_spec

    @property
    def profile_array(self):
        """'profile_array' property."""
        return self._profile_array

    @property
    def psize(self):
        """number of points in the profile"""
        return self._profile_size
