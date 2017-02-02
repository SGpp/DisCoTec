""" ProfileData.py

Contains the infrastructure to handle profile_spec_fext files
Also calculates nustar and the Chang-Hinton prediction
"""
# pylint: disable=E1101
import sys
import os
import math
import copy
import numpy as np
from bisect import bisect_left
from ParIO import Parameters


class ProfileFileData(object):
    """Class to handle the profile data for a single profile_{spec}{fext}

    """

    def __init__(self, filename, nx0, starttime, endtime):
        self.filename = filename
        self.nx0 = nx0
        self.timefld = []
        self.starttime = starttime
        self.endtime = endtime
        self.Ts = []
        self.ns = []
        self.omts = []
        self.omns = []
        self.Gammanc = []
        self.Gammaturb = []
        self.Qnc = []
        self.Qturb = []
        self.jbs = []
        self.Pinc = []
        self.Piturb = []

        self.readprof()

    def get_minmaxtime(self):
        return self.timefld[0], self.timefld[-1]

    def readprof(self):
        """ Read density, temperature, gradients and fluxes """
        blocks = []
        print('Reading  %s\n'%(self.filename))
        try:
            with open(self.filename) as prfile:
                next(prfile)  # skip header
                for line in prfile:
                    if not line or line.startswith('\n'):
                        continue
                    if line.startswith('#'):
                        self.timefld.append(float(line.split()[1]))
                        blocks.append([])
                    else:
                        blocks[-1].append(line)
        except IOError:
            sys.exit("IOError: probably profile file does not exist: %s"%self.filename)
        self.timefld = np.array(self.timefld)
        first_time, last_time = self.get_minmaxtime()
        if len(self.timefld) != len(set(self.timefld)):
            print("Error: %s contains 2 blocks with identical timestamp"%(self.filename))
            return
        if self.starttime == -1 or (self.starttime > 0 and
                                    self.starttime < first_time):
            print("Using first time present in profile data for starttime")
            self.starttime = first_time
        if self.endtime == -1:
            print("Using first time present in profile data for endtime")
            self.endtime = first_time
        if self.starttime == -2:
            print("Using last time present in profile data for starttime")
            self.starttime = last_time
        if (self.endtime == -2) or (self.endtime > last_time):
            print("Using last time present in profile data for endtime")
            self.endtime = last_time
        if (self.endtime < first_time) or (self.starttime > last_time):
            print("Time window not contained in profile data")
            return
        print(('starttime=%f, endtime=%f, first_time=%f, last_time=%f' %
              (self.starttime, self.endtime, first_time, last_time)))
        #  For a single time, find element closest to given input time
        if self.starttime == self.endtime:
            pos = np.array([bisect_left(self.timefld, self.starttime)])
        else:
            pos = np.where((self.timefld >= float(self.starttime)) &
                           (self.timefld <= self.endtime))[0]

        self._addprofdata(blocks, pos)
        self.timefld = self.timefld[pos]

    def _addprofdata(self, blocks, pos):
        """ Fill vectors with profile data """
        for entry in pos:
            self.Ts.append(np.zeros(self.nx0))
            self.ns.append(np.zeros(self.nx0))
            self.omts.append(np.zeros(self.nx0))
            self.omns.append(np.zeros(self.nx0))
            self.Gammaturb.append(np.zeros(self.nx0))
            self.Qturb.append(np.zeros(self.nx0))
            self.Piturb.append(np.zeros(self.nx0))
            self.Gammanc.append(np.zeros(self.nx0))
            self.Qnc.append(np.zeros(self.nx0))
            self.Pinc.append(np.zeros(self.nx0))
            self.jbs.append(np.zeros(self.nx0))
            for i in range(0, self.nx0):
                line = blocks[entry][i]
                self.Ts[-1][i] = float(line.split()[2])
                self.ns[-1][i] = float(line.split()[3])
                self.omts[-1][i] = float(line.split()[4])
                self.omns[-1][i] = float(line.split()[5])
                self.Gammaturb[-1][i] = float(line.split()[6])
                self.Qturb[-1][i] = float(line.split()[7])
                self.Piturb[-1][i] = float(line.split()[8])
                self.Gammanc[-1][i] = float(line.split()[9])
                self.Qnc[-1][i] = float(line.split()[10])
                self.Pinc[-1][i] = float(line.split()[11])
                self.jbs[-1][i] = float(line.split()[12])
            # Correct for convective energy flux
            self.Qnc[-1] -= 2.5 * self.Gammanc[-1] * self.Ts[-1]
            self.Qturb[-1] -= 2.5 * self.Gammaturb[-1] * self.Ts[-1]


class ProfileData(object):
    """Class handling the profile input

    Assignment to numpy vectors happens for a single given extension
    Special choices for starttime and endtime:
    -1: use first time in profile file
    -2: use last time in profile file

    """

    def __init__(self, fext, starttime, endtime, plotset):
        """ Read parameter file and create empty arrays for profile data """
        self.isDataPresent = False
        self.fext = fext
        # each dataset has its own collection of desired plots
        self.plotset = copy.deepcopy(plotset)
        # directory extension for the local equilibrium data... TODO: hardwired
        self.locncpathnr = "000"
        self.starttime = starttime
        self.endtime = endtime

        self.readpar()

        self.prspec = []
        self.xs = np.empty(self.pnt.nx0)
        self.T0s = np.empty((self.pnt.nx0, self.pnt.n_spec))
        self.n0s = np.empty((self.pnt.nx0, self.pnt.n_spec))
        self.omt0s = np.empty((self.pnt.nx0, self.pnt.n_spec))
        self.omn0s = np.empty((self.pnt.nx0, self.pnt.n_spec))

        # Chang-Hinton prediction for main ion species
        self.Qcharr = []
        # Ion-ion collisionality
        self.nustar = []

    def readpar(self):
        """ Get some essential parameters from par file"""
        par = Parameters()
        par.Read_Pars('parameters%s'%(self.fext))
        self.pnt = par.asnamedtuple()
        # self.nspec = par.pardict['n_spec']
        # self.rhostar = par.pardict['rhostar']
        # self.nz0 = par.pardict['nz0']
        #  if 'l_buffer_size' not in par.pardict:
        #      par..l_buffer_size = 0
        #  if 'u_buffer_size' not in par.pardict:
        #      self.pnt.ubuffer = 0
        # self.major_R = par.pardict['major_R']
        # self.minor_r = par.pardict['minor_r']
        q_coeffs = str(par.pardict['q_coeffs'])
        self.q_coeff_list = [float(qc) for qc in q_coeffs.split(",")]
        # self.coll = par.pardict['coll']
        self.specname = [[], []]
        for n in range(self.pnt.n_spec):
            self.specname[n] = str(par.pardict['name%1d' % (n+1)])
            self.specname[n] = self.specname[n][1:-1]
        # self.istep_prof = par.pardict['istep_prof']
        # self.comp_type = str(par.pardict['comp_type'])
        # self.nx0 = par.pardict['nx0']

        self.units = {'nref': 1.0, 'Tref': 1.0, 'Bref': 1.0, 'Lref': 1.0,
                      'mref': 1.0}
        for k in self.units:
            if k in par.pardict:
                self.units[k] = par.pardict[k]

    def q_fun(self):
        """ safety factor profile """
        x = self.xs
        qval = np.zeros(self.pnt.nx0)
        for i, qc in enumerate(self.q_coeff_list):
            qval += qc*x**i
        return qval

    def _chang_hinton(self, prspec):
        """ Calculate the Chang-Hinton estimate
        for the neoclassical heat flux """
        epsilon = np.array(self.xs*self.pnt.minor_r/self.pnt.major_R)
        mspec = 1  # since we only support one ion species so far
        # alpha = Zeff - 1
        Z = 1
        f1 = 1 + 1.5*epsilon**2
        f2 = (1-epsilon**2)**0.5
        # for alpha == 0: mustar == nustar
        for tb in range(len(prspec.Ts)):
            # k1 = 0 # relevant only for Zeff != 1
            mustar = self.nustar[tb]
            T = prspec.Ts[tb]
            n = prspec.ns[tb]
            omt = prspec.omts[tb]
            k2 = ((0.66 + 1.88*epsilon**0.5 - 1.54*epsilon) /
                  (1 + 1.03*mustar**0.5 + 0.31*mustar)*f1 +
                  1.6*0.74*mustar*epsilon /
                  (1+0.74*mustar*epsilon**1.5)*0.5*(f1-f2))
            qchval = (16.0/3/math.sqrt(math.pi) * n**2 *
                      self.q_fun()**2 / epsilon**1.5 * (1-epsilon**2) *
                      Z**4 * self.pnt.coll * np.sqrt(2*T*mspec) * k2 * omt)
            self.Qcharr.append(qchval)

    def _read_zero_prof(self):
        """Get the initial profile from profiles_spec_fext"""
        for n in range(self.pnt.n_spec):
            prof0file = open('profiles_%s%s'%(self.specname[n], self.fext))
            lines = prof0file.readlines()
            for i in range(0, self.pnt.nx0):
                l = 2+i
                self.xs[i] = lines[l].split()[0]
                self.T0s[i, n] = float(lines[l].split()[2])/self.units['Tref']
                self.n0s[i, n] = float(lines[l].split()[3])/self.units['nref']
                self.omt0s[i, n] = lines[l].split()[4]
                self.omn0s[i, n] = lines[l].split()[5]
            prof0file.close()

    def get_profiles(self):
        """ Organize the profile data, calls most other functions"""
        self._read_zero_prof()
        epsilon = np.array(self.xs*self.pnt.minor_r/self.pnt.major_R)
        for n in range(0, self.pnt.n_spec):
            filename = 'profile_%s%s'%(self.specname[n], self.fext)
            self.prspec.append(ProfileFileData(filename, self.pnt.nx0,
                                               self.starttime, self.endtime))
        # TODO: Add some consistency checks here
        for n, T in zip(self.prspec[0].ns, self.prspec[0].Ts):
            self.nustar.append(np.array(8.0/3.0/math.sqrt(math.pi) *
                               self.pnt.major_R*self.pnt.coll/epsilon**(3.0/2) *
                               self.q_fun()*n/T**2))
        self._chang_hinton(self.prspec[0])  # Only main ion species
        self._readlocalfile()
        self.isDataPresent = True

    def _readlocalfile(self):
        """Read the data from local runs set up with glob2loc"""
        locncfilename = "glob2locscan%s/neo.log"%(self.locncpathnr)
        if not os.path.isfile(locncfilename):
            self.plotset.discard('plotlocQ')
            print("No local NC equilibrium present")
            print(locncfilename)
            return
        if 'plotlocQ' in self.plotset:
            self.locxpos = []
            self.locqval = []
            self.locjbval = []
            locncfile = open(locncfilename)
            lines = locncfile.readlines()
            for lin in lines:
                li = lin.strip()
                if not li.startswith('#'):
                    self.locxpos.append(float(lin.split()[6]))
                    tmpnum = ((float(lin.split()[9]))*(float(lin.split()[2])) *
                              (float(lin.split()[4]))**2.5)
                    tmpnum = (tmpnum - 2.5*(float(lin.split()[8])) *
                              (float(lin.split()[2]))*(float(lin.split()[4]))**2.5)
                    self.locqval.append(tmpnum)
                    tmpnum = ((float(lin.split()[11])) *
                              (float(lin.split()[2]))*(float(lin.split()[4])))
                    self.locjbval.append(tmpnum)
            locncfile.close()


def glueprofiledata(prdlist):
    """ Take profdata objects with matching time intervals and combine them into a larger one
    Result is returned in the first element

    """
    # TODO: sanity checks for the parameters of prdlist entries
#    delta = 1e-6
#    if not (prd1.pnt.n_spec == prd2.pnt.n_spec and
#            abs(prd1.pnt.rhostar - prd2.pnt.rhostar) < delta and
#            prd1.pnt.nz0 == prd2.pnt.nz0 and
#            abs(prd1.pnt.major_R - prd2.pnt.major_R) < delta and
#            abs(prd1.pnt.minor_r - prd2.pnt.minor_r) < delta and
#            abs(prd1.pnt.coll - prd2.pnt.coll) < delta and
#            prd1.pnt.nx0 == prd2.pnt.nx0):
#                raise ValueError("Parameter mismatch when combining profile data")
    # Nothing to combine, nothing to do
    if len(prdlist) == 1:
        return prdlist[0]
    # Use first ProfileData as basis for the construction of the glued object
    result = prdlist.pop(0)
    for prd in prdlist:
        for n in range(result.pnt.n_spec):
            # When times overlap, give following ProfileData preference
            while result.prspec[n].timefld[-1] > prd.prspec[n].timefld[0]:
                del result.prspec[n].timefld[-1]
                del result.prspec[n].Ts[-1]
                del result.prspec[n].ns[-1]
                del result.prspec[n].omts[-1]
                del result.prspec[n].omns[-1]
                del result.prspec[n].Gammanc[-1]
                del result.prspec[n].Gammaturb[-1]
                del result.prspec[n].Qnc[-1]
                del result.prspec[n].Qturb[-1]
                del result.prspec[n].jbs[-1]
                del result.prspec[n].Pinc[-1]
                del result.prspec[n].Piturb[-1]
                if n == 0:
                    del result.Qcharr[-1]
                    del result.nustar[-1]
            result.prspec[n].timefld = np.concatenate((result.prspec[n].timefld, prd.prspec[n].timefld))
            result.prspec[n].Ts += prd.prspec[n].Ts
            result.prspec[n].ns += prd.prspec[n].ns
            result.prspec[n].omns += prd.prspec[n].omns
            result.prspec[n].omts += prd.prspec[n].omts
            result.prspec[n].Gammanc += prd.prspec[n].Gammanc
            result.prspec[n].Gammaturb += prd.prspec[n].Gammaturb
            result.prspec[n].Qnc += prd.prspec[n].Qnc
            result.prspec[n].Qturb += prd.prspec[n].Qturb
            result.prspec[n].jbs += prd.prspec[n].jbs
            result.prspec[n].Pinc += prd.prspec[n].Pinc
            result.prspec[n].Piturb += prd.prspec[n].Piturb
        result.nustar += prd.nustar
        result.Qcharr += prd.Qcharr
    result.endtime = prdlist[-1].endtime
    del prdlist
    return result


#     # TODO: Does not make sense here anymore
#     def avprof_output(self, avproffile_prefix):
#         print "Writing averaged profile to file"
#         for n in range(0, self.pnt.n_spec):
#             avprofilename = (avproffile_prefix + "_" + self.specname[n])#.strip()
#             avproffile = open(avprofilename, 'w')
#             avproffile.write('#   x/a             x/a(dummy)       T/Tref          n/nref\n')
#             avproffile.write('# Averaged profile generated by ProfileData.py: t = %4.1f - %4.1f\n'
#                               %(self.starttime, self.endtime))
#             for i in xrange(0, self.nx0):
#                 avproffile.write(str('%E'%self.xs[i]) + "\t" + str('%E'%self.xs[i]) + "\t")
#
#                avproffile.write(str('%E'%self.Ts[i, n]) + "\t" + str('%E'%self.ns[i, n]) + "\n")
#            avproffile.close()
# ==============================================================================
