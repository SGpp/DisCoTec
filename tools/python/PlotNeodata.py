"""PlotNeodata.py: Module to do the plotting """
# pylint: disable=E1101
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2.5
mpl.rc('xtick', labelsize=18)
mpl.rc('ytick', labelsize=18)
from ParIO import Parameters
import ProfileData # analysis:ignore
import FSAmomData # analysis:ignore
import read_mom_field.fieldlib as flib
from read_mom_field.averages import averages
from read_mom_field.common import common_data
from read_mom_field.geom import geometry


class PlotNeodata(object):
    """ Plot time averaged NC and FSA profiles

    Generates plots of temperature and density, particle, heat and momentum flux,
    bootstrap current, collisionality and flux surface averaged moments
    Arguments: lists of profile and fsamom class objects
    """

    def __init__(self, profdata, fsadata):
        self.leg1 = []
        self.leg2 = []
        self.leg3 = []
        self.leg4 = []
        self.leg5 = []
        self.leg6 = []
        self.leg7 = []
        self.prdlist = profdata
        self.fsalist = fsadata
        for prd in self.prdlist:
            if prd.isDataPresent is False:
                print("No profile data present in object ", prd.fext)
                prd.plotset.clear()
        for fsa in self.fsalist:
            if fsa.isDataPresent is False:
                print("No fsa moment data present in object ", prd.fext)
                fsa.plotset.clear()
        self.font = FontProperties(size=16)
        self.xyfs = 20
        self.xmin = min([prd.xs[1] for prd in self.prdlist])
        self.xmax = max([prd.xs[-1] for prd in self.prdlist])

    @classmethod
    def _get_erad(cls, prd, kspace=False):
        """Fetch the potential from the field diag and calculate E_r"""
        # TODO: This is redundant and ugly...
        par = Parameters()
        par.Read_Pars('parameters%s'%(prd.fext))
        pars = par.pardict
        phifile = flib.fieldfile("field"+prd.fext, pars)
        phifile.get_minmaxtime()
        times = []
        delta_t = int(pars['istep_field'])*float(pars['dt_max'])
        times.extend([time for time in phifile.tfld
                      if time > prd.prspec[0].timefld[0]-delta_t and
                      time < prd.prspec[0].timefld[-1]+delta_t])
        common = common_data(pars)
        common.tlen = len(times)
        common.times = np.array(times, dtype=np.float64)
        common.fileext = prd.fext
        common.rundir = ""
        geom = geometry(common, pars)
        avg = averages(common, geom)
        phi_ky0 = np.zeros((len(times), prd.pnt.nx0))
        for time in times:
            phifile.set_time(time)
            phi_ky0[times.index(time), :] = avg.z(phifile.phi())[0, :]
        phi_ky0_av = avg.t(phi_ky0)
        if kspace:
            # mirror at the inner boundary
            mirror_phi = np.concatenate((np.flipud(phi_ky0_av), phi_ky0_av))
            phi_ky0_kx = np.fft.rfft(mirror_phi)
            kx = np.fft.rfftfreq(2*prd.pnt.nx0, d=(prd.xs[1]-prd.xs[0])/prd.pnt.rhostar)
            erad = phi_ky0_kx * (-kx)
        else:
            erad = np.zeros(phi_ky0_av.shape, np.float)
            # erad[0:-1] = np.diff(phi_ky0_av)/np.diff(prd.xs)
            # erad[-1] = (phi_ky0_av[-1] - phi_ky0_av[-2])/(prd.xs[-1] - prd.xs[-2])
            erad = np.gradient(phi_ky0_av, (prd.xs[1]-prd.xs[0]))
            erad = -erad*prd.pnt.rhostar
        return erad

    def _calc_fb_param(self, prd, q, epsilon):
        """Calculate the k in the radial force balance equation"""
        erad = self._get_erad(prd)
        fb_param = []
        q = prd.q_fun()
        for prs in prd.prspec:
            fb_tmp = []
            for tb, jb in enumerate(prs.jbs):
                fb_tmp.append((-(jb/(q/epsilon*np.sqrt(1-epsilon**2) * prs.ns[tb]*prs.Ts[tb]) -
                              (prs.omns[tb] + prs.omts[tb] + erad/prs.Ts[tb]))/prs.omts[tb]))
            fb_param = self._mytrapz(fb_tmp, prs.timefld)
        return fb_param

    @classmethod
    def _mytrapz(cls, yvar, timefld):  # Trapezoid rule, adjusted for single timestep
        timefld = np.array(timefld)        
        if len(timefld) == 1:
            return yvar[0]
        else:
            tmpt = timefld / (timefld[-1] - timefld[0])
            return np.trapz(yvar, x=tmpt, axis=0)

    def _autocorrtime(self, vartx, timefld):
        corrtimes = np.zeros(len(vartx[0]))
        meanvar = self._mytrapz(vartx, timefld)
        for xind, var in enumerate(np.array(vartx).T):
            # The mean needs to be substracted to arrive at
            # the statistical definition of the autocorrelation
            var = var - meanvar[xind]
            result = np.correlate(var, var, mode="full")
            result = result[result.size/2:]
            overe = np.squeeze(np.where(result > np.exp(-1)*result[0]))
            # Find the first jump > 1 in the indices (avoiding long time correlations)
            cort_ind = np.array_split(overe, np.where(np.diff(overe) != 1)[0]+1)[0][-1]
            corrtimes[xind] = timefld[cort_ind]-timefld[0]
        return corrtimes

    def _windowerr(self, yvar, timefld):
        """ Calculate the windowed average and standard deviation at each radial position """
        if len(timefld) <= 2:
            print("Not enough times for error calculation")
            return np.zeros(len(yvar[0]))
        else:
            corrtimes = self._autocorrtime(yvar, timefld)
            yvar = np.array(yvar)
            window_std = np.empty(len(corrtimes))
            for ix, corrtx in enumerate(corrtimes):
                window_std[ix] = self._windowerr_atx(yvar[:, ix], corrtimes[ix], timefld)
            return window_std

    def _windowerr_atx(self, yvarx, corrtx, timefld):
        total_dt = timefld[-1] - timefld[0]
        # start with window width = 5 tcorr; if fewer than 10, reduce width, min 2
        # same as in IDL diag, tools.pro, stat_avg_err
        if total_dt >= 50 * corrtx:
            n_win = int(np.floor(total_dt/(5*corrtx)))
            wwidth = total_dt / n_win
        else:
            if total_dt >= 20 * corrtx:
                n_win = 10
                wwidth = total_dt / n_win
            else:
                return 0
        win_inds = []
        for n in range(n_win):
            win_inds.append(np.where(np.logical_and(
                np.greater_equal(timefld - timefld[0], n*wwidth),
                np.less(timefld - timefld[0], (n+1)*wwidth))))
        if np.any(np.diff(win_inds) <= 2):
            print('windowerr: windows with <= 2 steps, no errors computed')
            return 0
        means = np.empty(n_win)
        for iwin, window in enumerate(win_inds):
            means[iwin] = self._mytrapz(yvarx[window], timefld[window])
        std_error = np.std(means)
        return std_error

    def plot_T(self, poutput):
        """Plot the temperature and density profiles

        Profiles are temporally averaged by the trapezoid rule
        if the lists contain more than one time block
        The initial profiles are given as reference.

        """
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        for prd in self.prdlist:
            if 'plotT' not in prd.plotset:
                continue
            for prs in prd.prspec:
                ax1.plot(prd.xs, self._mytrapz(prs.Ts, prs.timefld), 'c')
                ax1.plot(prd.xs, self._mytrapz(prs.ns, prs.timefld), 'k')
                ax1.plot(prd.xs, self._mytrapz(prs.omts, prs.timefld)*prd.pnt.major_R, 'r')
                ax1.plot(prd.xs, self._mytrapz(prs.omns, prs.timefld)*prd.pnt.major_R, 'b')
                ax1.plot(prd.xs, prd.T0s[:, 0], 'c--')
                ax1.plot(prd.xs, prd.omt0s[:, 0]*prd.pnt.major_R, 'r--')
                ax1.plot(prd.xs, prd.n0s[:, 0], 'k--')
                ax1.plot(prd.xs, prd.omn0s[:, 0]*prd.pnt.major_R, 'b--')
                self.leg1.append(r'$T$')
                self.leg1.append(r'$n$')
                self.leg1.append(r'$R/L_T$')
                self.leg1.append(r'$R/L_n$')

        ax1.legend(self.leg1, loc=1, prop=self.font)
        ax1.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax1.set_xlim(self.xmin, self.xmax)
        ax1.set_ylabel(r'$T, n$ / $T_{ref}$, $n_{ref}$', fontsize=self.xyfs)
        if poutput == 1:
            fig1.frameon = False
            fig1.savefig('profiles_%s.pdf'%(self.prdlist[0].fext))

    def plot_Q(self, poutput):
        """Plot the averaged fluxes

        Plots the time averaged particle, heat and momentum fluxes.
        Standard is to plot only neoclassical fluxes and a prediction
        from the Chang-Hinton formula.

        """
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.set_color_cycle(['Crimson', 'DarkBlue', 'Red'])
        labely = ""
        for prd in self.prdlist:
            if ('plotQ' not in prd.plotset and
               'plotG' not in prd.plotset and 'plotP' not in prd.plotset):
                continue
            if 'plotch' in prd.plotset:
                ax2.plot(prd.xs, self._mytrapz(prd.Qcharr, prd.prspec[0].timefld), '--')
                self.leg2.append(r"Chang-Hinton")
            for n, prs in enumerate(prd.prspec):
                if 'plotQ' in prd.plotset:
                    meanQnc = self._mytrapz(prs.Qnc, prs.timefld)
                    errorQnc = self._windowerr(prs.Qnc, prs.timefld)
                    base_line, = ax2.plot(prd.xs, meanQnc)
                    ax2.fill_between(prd.xs, meanQnc-errorQnc, meanQnc+errorQnc,
                                     alpha=0.2, facecolor=base_line.get_color())
                    self.leg2.append(r"$Q_{nc}$ %s"%(prd.specname[n]))
                    if 'plotturb' in prd.plotset:
                        errorQt = self._windowerr(prs.Qturb, prs.timefld)
                        meanQt = self._mytrapz(prs.Qturb, prs.timefld)
                        base_line, = ax2.plot(prd.xs, meanQt)
                        ax2.fill_between(prd.xs, meanQt-errorQt, meanQt+errorQt,
                                         alpha=0.2, facecolor=base_line.get_color())
                        self.leg2.append(r"$Q_{turb}$ %s"%(prd.specname[n]))
                    if 'plotlocQ' in prd.plotset:
                        ax2.plot(prd.locxpos, prd.locqval, '+')
                        self.leg2.append(r"$Q_{nc}$ local")
                    if labely.find(r'$Q/Q_{gB}$') == -1:
                        labely += r'$Q/Q_{gB}$'
                if 'plotG' in prd.plotset:
                    ax2.plot(prd.xs, self._mytrapz(prs.Gammanc, prs.timefld))
                    self.leg2.append(r"$\Gamma_{nc}$ %s"%(prd.specname[n]))
                    if 'plotturb' in prd.plotset:
                        ax2.plot(prd.xs, self._mytrapz(prs.Gammaturb, prs.timefld))
                        self.leg2.append(r"$\Gamma_{turb} $%s"%(prd.specname[n]))
                    if labely.find(r'$\Gamma/\Gamma_{gB}$') == -1:
                        labely += r'  $\Gamma/\Gamma_{gB}$'
                if 'plotP' in prd.plotset:
                    ax2.plot(prd.xs, self._mytrapz(prs.Pinc, prs.timefld))
                    self.leg2.append(r"$\Pi_{nc}$ %s"%(prd.specname[n]))
                    if 'plotturb' in prd.plotset:
                        ax2.plot(prd.xs, self._mytrapz(prs.Piturb, prs.timefld))
                        self.leg2.append(r"$\Pi_{turb} $%s"%(prd.specname[n]))
                    if labely.find(r'$\Pi/\Pi_{gB}$') == -1:
                        labely += r'  $\Pi/\Pi_{gB}$'
        # self.leg2.append(r"$Q_{nc}$ (w. $\mathbf{v}_{D} \cdot \nabla f_1$)")
        # self.leg2.append(r"$Q_{nc}$ (wo. $\mathbf{v}_{D} \cdot \nabla f_1$)")
        ax2.legend(self.leg2, loc=1, prop=self.font)
        ax2.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax2.set_xlim(self.xmin, self.xmax)
        # ax2.set_xlim(0.04, 0.9)
        # ax2.set_ylim(0.0, 0.45)
        ax2.set_ylabel(labely, fontsize=self.xyfs)
        if poutput == 1:
            fig2.frameon = False
            fig2.savefig("radfluxes%s.pdf"%(self.prdlist[0].fext))

    def plot_chi(self, poutput):
        """Plot the heat conductivity """
        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111)
        for prd in self.prdlist:
            if 'plotchi' not in prd.plotset:
                continue
            if 'plotch' in prd.plotset:
                chitmp = []
                for tb, qchtmp in enumerate(prd.Qcharr):
                    chitmp.append(qchtmp*prd.prspec[0].Ts[tb]/prd.prspec[0].ns[tb]**2 /
                                  prd.prspec[0].omts[tb]/prd.prspec[0].omns[tb])
                ax3.plot(prd.xs, self._mytrapz(chitmp, prd.prspec[0].timefld), 'r')
                self.leg3.append(r"Chang-Hinton")
            for n, prs in enumerate(prd.prspec):
                chitmp = []
                for tb, qnctmp in enumerate(prs.Qnc):
                    chitmp.append(qnctmp/prs.Ts[tb]/prs.ns[tb]**2/prs.omts[tb]/prs.omns[tb])
                ax3.plot(prd.xs, self._mytrapz(chitmp, prs.timefld))
                self.leg3.append(r"$\chi_{nc}$  %s"%(prd.specname[n]))
                if 'plotturb' in prd.plotset:
                    chitmp = []
                    for tb, qtutmp in enumerate(prs.Qturb):
                        chitmp.append(qtutmp/prs.Ts[tb]/prs.ns[tb]**2/prs.omts[tb]/prs.omns[tb])
                    ax3.plot(prd.xs, self._mytrapz(chitmp, prs.timefld))
                    self.leg3.append(r"$\chi_{turb} $ %s"%(prd.specname[n]))
        ax3.legend(self.leg3, loc=1, prop=self.font)
        ax3.set_xlim(self.xmin, self.xmax)
        ax3.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax3.set_ylabel(r'$\chi [c_s(x_0) \rho_s^{2}(x_0)/L_n]$', fontsize=self.xyfs)
        if poutput == 1:
            fig3.frameon = False
            fig3.savefig("raddiffusivity%s.pdf"%(self.prdlist[0].fext))

    def plot_jbs(self, poutput):
        """Plot the species-wise bootstrap current

        Plots time averaged < u_\par B >_FS for every species.
        The bootstrap current would be their sum

        """
        fig4 = plt.figure()
        ax4_1 = fig4.add_subplot(121)
        ax4_2 = fig4.add_subplot(122)
        for prd in self.prdlist:
            if 'plotj' not in prd.plotset:
                continue
            epsilon = np.array(prd.xs*prd.pnt.minor_r/prd.pnt.major_R)
            q = prd.q_fun
            for n, prs in enumerate(prd.prspec):
                ax4_1.plot(prd.xs, self._mytrapz(prs.jbs, prs.timefld))
                self.leg4.append(r"$j_{b}$ %s"%(prd.specname[n]))
                # ax4_1.plot(prd.xs, np.array(q*q/epsilon*np.sqrt(1-epsilon**2)*prd.omns[:,n]*prd.ns[:,n]**2 *prd.Ts[:,n] *( 1.46*np.sqrt(epsilon) - 1.35 * np.sqrt(prd.nustar[:,n])/epsilon**0.25)))
                # self.leg4.append(r"Prediction")
                if 'plotlocQ' in prd.plotset:
                    ax4_1.plot(prd.locxpos, np.array(prd.locjbval), '+')
                    self.leg4.append(r"$j_{bs}$ local")
            fb_param = self._calc_fb_param(prd, q, epsilon)
            ax4_2.plot(prd.xs, fb_param)
            # ax4_1.plot(prd.xs, np.array(q/epsilon*np.sqrt(1-epsilon**2)*
            # prd.ns[:,n]*prd.Ts[:,n]*(prd.omns[:,n] +  prd.omts[:,n] - fb_param*prd.omts[:,n])))
            # still lacks erad!!!
        ax4_1.legend(self.leg4, loc=1, prop=self.font)
        ax4_1.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax4_2.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax4_1.set_xlim(self.xmin, self.xmax)
        ax4_2.set_xlim(self.xmin, self.xmax)
        ax4_1.set_ylabel(r'$j_B$', fontsize=self.xyfs)
        if poutput == 1:
            fig4.frameon = False
            fig4.savefig("j_botstrap%s.pdf"%(self.prdlist[0].fext))

    def plot_nustar(self, poutput):
        """ Plot the collisionality and the inverse aspect ratio for comparison """
        fig5 = plt.figure()
        ax5 = fig5.add_subplot(111)
        mspec = 1  # TODO: include species mass in prd!!!!
        for prd in self.prdlist:
            if 'plotnu' not in prd.plotset:
                continue
            eps = prd.xs*prd.pnt.minor_r/prd.pnt.major_R
            ax5.plot(prd.xs, self._mytrapz(prd.nustar, prd.prspec[0].timefld) * eps**1.5)
            self.leg5.append(r"${\nu^*}_{i} \cdot \epsilon^{3/2}$ ")
            ax5.plot(prd.xs, eps**1.5)
            self.leg5.append(r"$\epsilon^{3/2}$")
            q = prd.q_fun()
            colltmp = []
            for tb, nu in enumerate(prd.nustar):
                colltmp.append((nu/q/prd.pnt.major_R *
                               np.sqrt(2*prd.prspec[0].Ts[tb]/mspec)*eps**1.5)**(-1))
            # This is tau_ii as defined in Helander/Sigmar, note tau_i = sqrt(2) * tau_ii
            colltime = self._mytrapz(colltmp, prd.prspec[0].timefld)
            ax5.text(0.1, 0.5, "Highest collision time: %f"%np.max(colltime))
            ax5.text(0.1, 0.3, "Lowest collision time: %f"%np.min(colltime))
        ax5.set_yscale('log')
        ax5.legend(self.leg5, loc=1, prop=self.font)
        ax5.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax5.set_xlim(self.xmin, self.xmax)
        if poutput == 1:
            fig5.frameon = False
            fig5.savefig("nustar%s.pdf"%(self.prdlist[0].fext))

    def plot_erad(self, poutput):  # based on Vernay2012 eq. 45
        """Plot the radial electric field and the ExB shear """
        fig6 = plt.figure()
        ax6_1 = fig6.add_subplot(211)
        ax6_2 = fig6.add_subplot(212)
        for prd in self.prdlist:
            if 'ploterad' not in prd.plotset:
                continue
            erad = self._get_erad(prd)
            ax6_1.plot(prd.xs, erad)
            #  slope, intercept, r_value, p_value,
            #  std_err=scipy.stats.linregress(prd.omns[:, 0], erad)
            #  print slope, intercept, r_value, p_value, std_err
            #  ax7_1.plot(prd.xs, slope*prd.omns[:,0]+intercept,label="Regression")
            #  ax7_2.plot(prd.omns[:,0],erad,'+')
            #  ax7_2.plot(prd.omns[:,0],slope*prd.omns[:,0]+intercept,'-',label="Regression")
            #  temp = erad * prd.q_fun()/prd.xs
            temp = self._get_erad(prd, kspace=True)
            # wexb = np.zeros(prd.nx0)
            wexb = np.gradient(erad, (prd.xs[1]-prd.xs[0])) * prd.pnt.rhostar
            # wexb[0:-1] = np.diff(temp)/np.diff(prd.xs)
            # wexb[-1] = (temp[-1] - temp[-2])/(prd.xs[-1] - prd.xs[-2])
            # wexb = -wexb*prd.rhostar*prd.xs/prd.q_fun() #There's a B_0 missing here
            kx = np.fft.rfftfreq(2*prd.pnt.nx0, d=(prd.xs[1]-prd.xs[0])/prd.pnt.rhostar)
            # print kx, len(temp)
            # wexb_fft = temp*kx
            # ax6_2.plot(kx,kx*np.abs(temp), '-+')
            ax6_2.plot(prd.xs, 1/prd.pnt.rhostar * np.fft.irfft(-kx*temp)[prd.pnt.nx0:])
            ax6_2.plot(prd.xs, wexb)
            # ax6_2.plot(kx,np.imag(temp), '-+')
            # ax6_2.set_yscale('log')
            # ax6_2.set_xscale('log')
    #    ax6_1.set_xlim(0.25,0.7)
        ax6_1.set_xlabel(r'$x/a$', fontsize=self.xyfs)
    #    ax6_1.legend(["w. NC", "wo. NC"])
        # ax6_2.set_xlabel(r'$x/a$',fontsize=self.xyfs)
        ax6_1.set_ylabel(r'$E_r/(T/e)$', fontsize=self.xyfs)
        # ax6_2.set_ylabel(r'$\omega_{E\times B}$',fontsize=self.xyfs)
        if poutput == 1:
            fig6.frameon = False
            fig6.savefig("erad.pdf")  # %(self.prdlist[0].fext))

    def plot_fsa(self, poutput):
        fig7 = plt.figure()
        ax7 = fig7.add_subplot(111)
        for fsa in self.fsalist:
            if 'plotfsa' not in fsa.plotset:
                continue
            for fsd in fsa.fsaspec:
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 2])
                self.leg7.append(r"$dE/dt_{tot}$")
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 4])
                self.leg7.append(r"$d\chi/dxy$")
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 5])
                self.leg7.append(r"$f_1 sources$")
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 6])
                self.leg7.append(r"$f_1 buffer$")
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 7])
                self.leg7.append(r"$f_0 term$")
                ax7.plot(fsd.fsadata[0][:, 0], self._mytrapz(fsd.fsadata, fsd.timefld)[:, 8])
                self.leg7.append(r"$-v_D f_1$")
                ax7.legend(self.leg7, loc=1, prop=self.font)
        ax7.set_xlabel(r'$x/a$', fontsize=self.xyfs)
        ax7.set_xlim(self.xmin, self.xmax)
        ax7.set_ylabel(r'$\int \frac{dE}{dt} dx$', fontsize=self.xyfs)
        if poutput == 1:
            fig7.frameon = False
            fig7.savefig("fsamoments%s.pdf"%(self.prdlist[0].fext))

    @classmethod
    def show(cls):
        plt.show()


class PlotNC_xt(object):
    """Plot contours of profiles and fluxes in x-t space """

    def __init__(self, profdata):
        self.prdlist = profdata
        for prd in self.prdlist:
            if prd.isDataPresent is False:
                print("No profile data present in object ", prd.fext)
                prd.plotset.clear()
        self.font = FontProperties(size=16)
        self.xyfs = 20
        self.xmin = min([prd.xs[1] for prd in self.prdlist])
        self.xmax = max([prd.xs[-1] for prd in self.prdlist])

    def plotT_cont(self, poutput):
        for prd in self.prdlist:
            if 'plotT' not in prd.plotset:
                continue
            fig1 = plt.figure()
            for num, prs in enumerate(prd.prspec):
                ax1 = fig1.add_subplot(len(prd.prspec), 1, num)
                ax1.pcolormesh(prd.xs, prs.timefld, np.array(prs.omts))
                ax1.set_xlabel(r"$x/a$", fontsize=self.xyfs)
                ax1.set_ylabel(r"$t$", fontsize=self.xyfs)
            if poutput == 1:
                fig1.frameon = False
                fig1.savefig('profiles_xt%s.pdf'%(prd.fext))

    def plotQ_cont(self, poutput):
        for prd in self.prdlist:
            if 'plotQ' not in prd.plotset:
                continue
            if 'plotch' in prd.plotset:
                fig1 = plt.figure()
                ax1 = fig1.add_subplot(111)
                ax1.pcolormesh(prd.xs, prd.prspec[0].timefld, np.array(prd.Qcharr),
                               shading="gouraud")
                ax1.set_xlabel(r"$x/a$", fontsize=self.xyfs)
                ax1.set_ylabel(r"$t$", fontsize=self.xyfs)
            fig1 = plt.figure()
            for num, prs in enumerate(prd.prspec):
                ax1 = fig1.add_subplot(2, len(prd.prspec), num+1)
                cm1 = ax1.pcolormesh(prd.xs, prs.timefld, np.array(prs.Qnc), shading="gouraud")
                # , vmin=-3, vmax=4)
                fig1.colorbar(cm1)
#                ax2 = fig1.add_subplot(len(prd.prspec),3,2*(num+1))
#                cm2 = ax2.pcolormesh(prd.xs, prs.timefld, np.array(prs.Qnc)-np.array(prd.Qcharr),
#                              shading="gouraud", cmap=mpl.cm.Paired)
#               fig1.colorbar(cm2)
                ax3 = fig1.add_subplot(2, len(prd.prspec), 2*(num+1))
                cm3 = ax3.pcolormesh(prd.xs, prs.timefld, np.array(prs.Qturb), shading="gouraud")
                # , vmin=0, vmax=24.5)
                fig1.colorbar(cm3)
                ax1.set_title(r"$Q_{nc}$", fontsize=self.xyfs)
                ax3.set_title(r"$Q_{turb}$", fontsize=self.xyfs)
                ax1.set_xlabel(r"$x/a$", fontsize=self.xyfs)
                ax1.set_xlim(0.1, 0.9)
                ax3.set_xlim(0.1, 0.9)
                ax3.set_xlabel(r"$x/a$", fontsize=self.xyfs)
                ax1.set_ylabel(r"$t$", fontsize=self.xyfs)
                ax3.set_ylabel(r"$t$", fontsize=self.xyfs)
                plt.tight_layout()
                fig1.frameon = False
            if poutput == 1:
                fig1.frameon = False
                fig1.savefig('radfluxes_xt%s.pdf'%(prd.fext))

    def plotjbs_cont(self, poutput):
        for prd in self.prdlist:
            if 'plotj' not in prd.plotset:
                continue
            fig1 = plt.figure()
            for num, prs in enumerate(prd.prspec):
                ax1 = fig1.add_subplot(len(prd.prspec), 1, num)
                cm1 = ax1.pcolormesh(prd.xs, prs.timefld, np.array(prs.jbs))
                fig1.colorbar(cm1)
                ax1.set_xlabel(r"$x/a$", fontsize=self.xyfs)
                ax1.set_ylabel(r"$t$", fontsize=self.xyfs)
            if poutput == 1:
                fig1.frameon = False
                fig1.savefig('profiles_xt%s.pdf'%(prd.fext))

    @classmethod
    def show(cls):
        plt.show()
