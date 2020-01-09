from Tkinter import *
import machine as machine
import os
from tkMessageBox import askokcancel

class Variables:
    def __init__(self,launch):
        global launcher
        launcher=launch
        # find out machine name
        mach=machine.machine(launcher,self)
        # standard background color "#rgb" specified in hexadecimals
        self.stdbg="#ccc"
        # Tkinter variables have to be initialized like: var = StringVar() / IntVar() / DoubleVar()
        # remembers the prob directory the current job came from
        self.jobpath=StringVar()
        #jobpath.set('')
        # variables that signify first opening of the respective frame
        self.Opfirst=1
        self.Inoutfirst=1
        self.Domfirst=1
        self.Genfirst=1
        self.Collfirst=1
        self.Geofirst=1
        self.Specfirst=1
        self.Reffirst=1
        self.Nonlocalfirst=1
        self.Devfirst=1
        # Dictionary containing the namelist each parameter belongs to
        self.nmldict={}
        # Dictionary containing parameters names and values
        self.pardict={}
        # switch for expert mode
        self.expertV=IntVar()
        self.expertV.set(0)
        self.form_cleared=0
        # for multiple job management: number of current job, number of jobs to be managed (just change the number!)
        self.currentjob=0
        self.nummulti=5
        self.job_saved=[0]*self.nummulti
    
        self.allow_nonlocal=1
        
        # Define variables for Opframe
        self.nl=IntVar()
        self.ivev=IntVar()
        self.nonlocal=IntVar()
        self.nonlocal.set(0)
        self.scan=IntVar()
        self.scan.set(0)
        # Time stepping scheme
        self.timeschV=StringVar()
        # parallel derivative scheme
        self.parschV=StringVar()
        # Number of processes to distribute each dimension to
        self.sV=StringVar()
        self.vV=StringVar()
        self.wV=StringVar()
        self.xV=StringVar()
        self.yV=StringVar()
        self.zV=StringVar()
        # auto-parallelize the dimension or not
        self.sautoV=IntVar()
        self.vautoV=IntVar()
        self.wautoV=IntVar()
        self.xautoV=IntVar()
        self.yautoV=IntVar()
        self.zautoV=IntVar()
        self.npsautoV=IntVar()
        # total number of processes
        self.nprocsV=IntVar()
        #number of parallel simulations
        self.nparsimV=IntVar()
        
        # Hyperdiffusion in x,z,v directions
        self.hypzV=StringVar()
        self.hypzordV=StringVar()
        self.hypxV=StringVar()
        self.hypxordV=StringVar()
        self.hypyV=StringVar()
        self.hypyordV=StringVar()
        self.hypvV=StringVar()

        # Parallel Arakawa treatment
        self.arakawazvV=IntVar()
        
        # Define variables for Inoutframe
        self.diagdirV=StringVar()
        self.diagdirV.set("'"+self.outdir.strip("'")+"'")
        self.chptdirV=StringVar()
        self.rdchptV=IntVar()
        self.wrtchptV=IntVar()
        self.chpth5V=IntVar()
        self.wrth5V=IntVar()
        self.wrtstdV=IntVar()
        self.iterdbfileV=StringVar()
        self.iterdbtimeV=StringVar()
        self.istepnrgV=StringVar()
        self.istepfldV=StringVar()
        self.istepomegaV=StringVar()
        self.istepmomV=StringVar()
        self.istepenergyV=StringVar()
        self.istepenergy3dV=StringVar()
        self.istepprofV=StringVar()
        self.istepschptV=StringVar()
        self.istepvspV=StringVar()
        self.istepncV=StringVar()
        
        # Define variables for Domframe
        # Grid resolutions
        self.nxV=StringVar()
        self.nyV=StringVar()
        self.nzV=StringVar()
        self.nvV=StringVar()
        self.nwV=StringVar()
        # Box parameters
        self.kyminV=StringVar()
        self.lxV=StringVar()
        self.lvV=StringVar()
        self.lwV=StringVar()
        self.ky0indV=StringVar()
        self.nexcV=StringVar()
        self.kxcentV=StringVar()
        self.adaptlx=IntVar()
        self.adaptly=IntVar()
        self.x0V=StringVar()
        
        # Define variables for Genframe
        # Initial condition
        self.initcondV=StringVar()
        self.auxxV=StringVar()
        self.auxyV=StringVar()
        self.auxzV=StringVar()
        # Simulation time
        self.ntimeV=StringVar()
        self.timelimV=StringVar()
        self.simtimelimV=StringVar()
        # Timestep
        self.dtmaxV=StringVar()
        self.calcdtV=IntVar()
        self.courantV=StringVar()
        # Frequency precision
        self.omprV=StringVar()
        # Overflow limit
        self.oflV=StringVar()
        self.uflV=StringVar()
        
        # Dimensionless parameters
        self.betaV=StringVar()
        self.debye2V=StringVar()
                
        # Eigenvalue calculations
        self.whichevV=StringVar()
        self.nevV=StringVar()
        self.maxitV=StringVar()
        self.shiftV=StringVar()
        
        # Switches
        self.pressV=IntVar()
        self.delzonalV=IntVar()
#        self.quasilintmV=IntVar()        
        self.profswV=IntVar()

        #Neoclassical
        self.incf0V=IntVar()
        self.delphiV=IntVar()
        
        # Flow shear
        self.ExBrateV=StringVar()
        self.pfsrateV=StringVar()
        self.ExBtimeV=StringVar()
        self.kxExBV=StringVar()
        self.PhikxV=StringVar()
        
        # Define variables for Collframe
        self.collV=StringVar()
        self.ZeffV=StringVar()
        self.collisV=StringVar()
        self.collonoffV=IntVar()
        self.spacediffV=IntVar() 
        self.collcmV=StringVar()
        self.collffmV=IntVar()

       
        # Define variables for Geoframe
        self.geomread=0
        self.geomtype=StringVar()
        self.minorrV=StringVar()
        self.lengthV=StringVar()
        self.shatV=StringVar()
        self.trpepsV=StringVar()
        self.q0V=StringVar()
        self.geomdirV=StringVar()
        self.geomfileV=StringVar()
        self.fluxlabelV=StringVar()
        self.fluxposV=StringVar()
        self.rrefV=StringVar()
        self.amhdV=StringVar()
        self.kappaV=StringVar()
        self.s_kappaV=StringVar()
        self.deltaV=StringVar()
        self.s_deltaV=StringVar()
        self.zetaV=StringVar()
        self.s_zetaV=StringVar()
        self.drRV=StringVar()
        self.rhostarV=StringVar()
        self.magprofV=IntVar()
        
        # Define variables for Specframe
        # Number of species
        self.numspecV=StringVar()
        
        # Species currently displayed in the lhs/rhs entry fields
        self.currspec1V=IntVar()
        self.currspec2V=IntVar()
        
        # Variables for lhs entry fields
        self.specname1V=StringVar()
        self.omn1V=StringVar()
        self.omt1V=StringVar()
        self.mass1V=StringVar()
        self.charge1V=StringVar()
        self.temp1V=StringVar()
        self.dens1V=StringVar()
        self.passive1V=IntVar()
        self.fileV=StringVar()
        self.filetimeV=StringVar()
        # nonlocal
        self.proftypeV=IntVar()
        
        
        # Variables for rhs entry fields
        self.specname2V=StringVar()
        self.omn2V=StringVar()
        self.omt2V=StringVar()
        self.mass2V=StringVar()
        self.charge2V=StringVar()
        self.temp2V=StringVar()
        self.dens2V=StringVar()
        self.passive2V=IntVar()
        
        # nonlocal

        #profile propagation
        self.proftoTrefV=IntVar()
        self.proftonrefV=IntVar()
        self.proftobetaV=IntVar()
        self.proftocollV=IntVar()
        self.proftodebyeV=IntVar()
        self.proftorhostarV=IntVar()
        self.proftoExBV=IntVar()
        self.proftopfsV=IntVar()
        
        # Create species arrays
        self.specname=[]
        self.omn=[]
        self.omt=[]
        self.mass=[]
        self.charge=[]
        self.temp=[]
        self.dens=[]
        self.passive=[]
        # the following ones are for global runs
        self.proftype=[]
        self.kappa_T=[]
        self.kappa_n=[]
        self.LT_center=[]
        self.Ln_center=[]
        self.LT_width=[]
        self.Ln_width=[]
        self.delta_x_T=[]
        self.delta_x_n=[]
        self.src_amp=[]
        self.src_width=[]
        self.src_x0=[]
        self.src_proftype=[]
        
        # Define variables for Nonlocalframe
        self.radbctypeV=StringVar()
        self.lbufsizeV=StringVar()
        self.lcoefkrookV=StringVar()
        self.lpowkrookV=StringVar()
        self.ubufsizeV=StringVar()
        self.ucoefkrookV=StringVar()
        self.upowkrookV=StringVar()
        self.arakawaV=IntVar()
        self.ckheatV=StringVar()
        self.resetlimitV=StringVar()
        self.shiftedV=IntVar()
        
        # Define variables for Refframe
        self.TrefV=StringVar()
        self.nrefV=StringVar()
        self.LrefV=StringVar()
        self.BrefV=StringVar()
        self.mrefV=StringVar()
        self.betaRefV=StringVar()
        self.collRefV=StringVar()
        self.debyeRefV=StringVar()
        self.rhostarV=StringVar()
        self.rhostarRefV=StringVar()
        
        
        # Define variables for Devframe
        self.listofnamelists=['parallelization','box','in_out','general','geometry','units','external_contr','nonlocal_x','scan']
        self.addedpars=[]
        self.addedvalues=[]
        self.namelistM=[]
        self.namelistMB=[]
        self.namelist=[]
        self.EntryDevPar=[]
        self.EntryDevVal=[]
        
        # for Check_Pars, number is non-zero if errors have been found
        self.pars_okay=0
        
        # for Tooltips
        self.follow=TRUE
        self.delay=0.1

    def reset_listofnamelists(self):
        self.listofnamelists=['parallelization','box','in_out','general','geometry','units','external_contr','nonlocal_x','scan']

    def clear_button(self):
        if askokcancel(message='This will reset all entries, please confirm!'):
            self.clear()
        else:
            return
        
    def clear(self):
        self.reset_listofnamelists()
        self.geomread=0
        self.scan.set(0)
        self.numspecV.set(0)
        #self.nonlocal.set(0)
        self.omn=[]
        self.omt=[]
        self.mass=[]
        self.charge=[]
        self.temp=[]
        self.dens=[]
        self.specname=[]
        self.passive=[]
        self.kappa_T=[]
        self.kappa_n=[]
        self.LT_center=[]
        self.Ln_center=[]
        self.LT_width=[]
        self.Ln_width=[]
        self.delta_x_T=[]
        self.delta_x_n=[]
        self.proftype=[]
        self.src_proftype=[]
        self.src_amp=[]
        self.src_width=[]
        self.src_x0=[]
        self.form_cleared=1
        #self.nl.set(1)
        self.adaptlx.set(0)
        self.adaptly.set(0)
        launcher.adaptlx_sel
        self.timeschV.set("'RK4'")
        self.parschV.set("'c4th'")
        self.sV.set('0')
        self.vV.set('0')
        self.wV.set('0')
        self.xV.set('0')
        self.yV.set('0')
        self.zV.set('0')
        self.sautoV.set(1)
        self.vautoV.set(1)
        self.wautoV.set(1)
        self.xautoV.set(1)
        self.yautoV.set(1)
        self.zautoV.set(1)
        self.nprocsV.set(0)
        self.hypzV.set('')
        self.hypzordV.set(4)
        self.hypxV.set('')
        self.hypxordV.set(4)
        self.hypyV.set('')
        self.hypyordV.set(4)
        self.hypvV.set('')
        self.arakawazvV.set(1)
        self.rdchptV.set(0)
        self.wrtchptV.set(1)
        self.chpth5V.set(0)
        self.wrth5V.set(0)
        self.wrtstdV.set(1)
        self.istepnrgV.set('')
        self.istepomegaV.set('')
        self.istepfldV.set('')
        self.istepmomV.set('')
        self.istepenergyV.set('')
        self.istepenergy3dV.set('')
        self.istepprofV.set('')
        self.istepvspV.set('')
        self.istepncV.set('')
        self.istepschptV.set('')
        self.nxV.set('')
        self.nyV.set('')
        self.nzV.set('')
        self.nvV.set('')
        self.nwV.set('')
        self.kyminV.set('')
        self.lxV.set('')
        self.lvV.set('')
        self.lwV.set('')
        self.ky0indV.set('')
        self.nexcV.set('')
        self.kxcentV.set('')
        self.initcondV.set("'ppj'")
        self.auxxV.set(-100)
        self.auxyV.set(-100)
        self.auxzV.set(-100)
        self.ntimeV.set("")
        self.timelimV.set('')
        self.simtimelimV.set('')
        self.dtmaxV.set('')
        self.calcdtV.set(1)
        launcher.calcdtclick()
        self.omprV.set('')
        self.oflV.set('')
        self.uflV.set('')
        self.betaV.set('')
        self.ExBrateV.set('')
        self.pfsrateV.set('')
        self.ExBtimeV.set('')
        self.kxExBV.set('')
        self.PhikxV.set('')
        self.courantV.set('')
        self.incf0V.set(0)
        self.delphiV.set(0)
        self.collV.set('')
        self.ZeffV.set('1.0')
        self.collisV.set("'none'")
        self.collonoffV.set(0)
        self.collcmV.set("'xu_rosenbluth'")
        self.collffmV.set(0)
        self.spacediffV.set(0)
        self.debye2V.set('')
        self.whichevV.set("'none'")
        self.nevV.set('')
        self.maxitV.set('')
        self.shiftV.set('')
        self.pressV.set(0)
        self.delzonalV.set(0)
#        self.quasilintmV.set(0)
        self.geomtype.set('')
        self.minorrV.set('')
        self.lengthV.set('')
        self.shatV.set('')
        self.trpepsV.set('')
        self.q0V.set('')
        self.amhdV.set('')
        self.kappaV.set('')
        self.s_kappaV.set('')
        self.deltaV.set('')
        self.s_deltaV.set('')
        self.zetaV.set('')
        self.s_zetaV.set('')
        self.drRV.set('')
        self.geomdirV.set("''")
        self.geomfileV.set("''")
        self.fluxlabelV.set('')
        self.fluxposV.set('')
        self.rrefV.set('')
        self.specname1V.set("''")
        self.dens1V.set('')
        self.mass1V.set('')
        self.temp1V.set('')
        self.charge1V.set('')
        self.omn1V.set('')
        self.omt1V.set('')
        self.proftypeV.set(0)
        self.passive1V.set(0)
        self.specname2V.set("''")
        self.dens2V.set('')
        self.mass2V.set('')
        self.temp2V.set('')
        self.charge2V.set('')
        self.omn2V.set('')
        self.omt2V.set('')
        self.passive2V.set(0)
        self.proftoTrefV.set(1)
        self.proftonrefV.set(1)
        self.proftobetaV.set(1)
        self.proftocollV.set(1)
        self.proftodebyeV.set(1)
        #set to zero by default
        self.proftorhostarV.set(0)
        self.proftoExBV.set(0)
        self.proftopfsV.set(0)

        self.radbctypeV.set(1)
        self.lcoefkrookV.set(0.)
        self.lpowkrookV.set(4)
        self.lbufsizeV.set(0.0)
        self.ucoefkrookV.set(0.)
        self.upowkrookV.set(4)
        self.ubufsizeV.set(0.0)
        self.resetlimitV.set(1000.)
        self.ckheatV.set(0.)
        self.TrefV.set('')
        self.LrefV.set('')
        self.BrefV.set('')
        self.nrefV.set('')
        self.mrefV.set('')
        for i in range(len(self.EntryDevPar)):
            self.addedpars[i].set('')
            self.addedvalues[i].set('')
            self.namelist[i].set('')
        launcher.switch_nonlocal()
        launcher.switch_lin_nonlin()
        self.ch_spec1()
        self.job_saved[self.currentjob]=0
        launcher.autosw()
        launcher.switch_scan()
        launcher.switch_collop()

           
    def set_defaults(self):
        self.clear()
        self.timeschV.set("'RK4'")
        self.parschV.set("'c4th'")
        self.sV.set(2)
        self.vV.set(1)
        self.wV.set(4)
        self.xV.set(1)
        self.yV.set(1)
        self.zV.set(1)
        self.sautoV.set(0)
        self.vautoV.set(0)
        self.wautoV.set(0)
        self.xautoV.set(0)
        self.yautoV.set(0)
        self.zautoV.set(0)
        self.nprocsV.set(8)
        self.hypzV.set(2.0)
        self.hypzordV.set(4)
        self.hypxV.set(0.)
        self.hypxordV.set(4)
        self.hypyV.set(0.)
        self.hypyordV.set(4)
        self.hypvV.set(0.2)
        self.arakawazvV.set(1)
        self.rdchptV.set(0)
        self.wrtchptV.set(1)
        self.chpth5V.set(0)
        self.wrth5V.set(0)
        self.wrtstdV.set(1)
        self.istepnrgV.set(10)
        self.istepomegaV.set(20)
        self.istepfldV.set(100)
        self.istepmomV.set(100)
        self.istepenergyV.set(0)
        self.istepenergy3dV.set(0)
        self.istepprofV.set(100)
        self.istepvspV.set(500)
        self.istepncV.set(0)
        self.istepschptV.set(5000)
        self.lvV.set(3.0)
        self.lwV.set(9.0)
        self.x0V.set(0.5)
        self.kxcentV.set(0.)
        self.auxxV.set(-100)
        self.auxyV.set(-100)
        self.auxzV.set(-100)
        self.ntimeV.set(500000)
        self.omprV.set(1.e-3)
        #self.oflV.set()
        self.uflV.set(1.e-15)
        self.betaV.set(0.001)
        self.courantV.set(1.25)
        self.incf0V.set(0)
        self.delphiV.set(0)
        self.collV.set(0.0)
        self.ZeffV.set(1.0)
        self.collisV.set("'none'")
        self.collonoffV.set(0)
        self.collcmV.set("'xu_rosenbluth'")
        self.collffmV.set(0)
        self.spacediffV.set(0)
        self.ExBrateV.set(0.0)
        self.pfsrateV.set(0.0)
        self.ExBtimeV.set(0.0)
        self.kxExBV.set(0)
        self.PhikxV.set(0.0)
        self.debye2V.set(0.0)
        self.whichevV.set("'none'")
        self.nevV.set(1)
        self.maxitV.set(0)
        self.shiftV.set("(20,0)")
        self.pressV.set(0)
        self.delzonalV.set(0)
#        self.quasilintmV.set(0)
        self.geomtype.set("'circular'")
        self.minorrV.set(0.35)
        self.lengthV.set(1)
        if not self.nonlocal.get():
            self.shatV.set(0.796)
            self.trpepsV.set(0.18)
        if self.nonlocal.get():
            self.q0V.set('0.868, 0., 2.2')
            self.magprofV.set(1)
        else:
            self.q0V.set(1.4)
            self.magprofV.set(0)
        self.amhdV.set(0.0)
        self.kappaV.set(1.0)
        self.s_kappaV.set(0.0)
        self.deltaV.set(0.0)
        self.s_deltaV.set(0.0)
        self.zetaV.set(0.0)
        self.s_zetaV.set(0.0)
        self.drRV.set(0.0)
        if self.nonlocal.get():
            self.rhostarV.set(0.0125)
        else:
            self.rhostarV.set('')
        self.geomdirV.set("''")
        self.geomfileV.set("''")
        self.fluxlabelV.set("'arho_t'")
        self.fluxposV.set(0.0)
        self.rrefV.set(0.0)
        self.resetlimitV.set(1000.)
        self.ckheatV.set(0.)
        self.lcoefkrookV.set(1.0)
        self.ucoefkrookV.set(1.0)
        self.lbufsizeV.set(0.1)
        self.ubufsizeV.set(0.1)
        self.arakawaV.set(1)
        self.newspec('CBC')
        self.newspec('CBC')
        self.numspecV.set(2)
        self.ch_spec1(1)
        self.ch_spec2(2)
        self.TrefV.set(3)
        self.nrefV.set(5)
        self.BrefV.set(2.5)
        self.LrefV.set(1.65)
        self.mrefV.set(1.)
        self.form_cleared=0
        if self.nl.get()==1:
            self.nxV.set(64)
            self.nyV.set(16)
            self.nzV.set(16)
            self.nvV.set(32)
            self.nwV.set(8)
            self.kyminV.set(0.05)
            if self.nonlocal.get():
                self.lxV.set(60)                
            else:
                self.lxV.set(125.528)
            self.ky0indV.set(0)
            self.nexcV.set(5)
            self.adaptlx.set(0)
            self.adaptly.set(0)
            if self.nonlocal.get():
                self.initcondV.set("'db'")
            else:
                self.initcondV.set("'ppj'")
            self.timelimV.set(86000)
            self.dtmaxV.set(0.0385)
            self.calcdtV.set(1)
            self.job_saved[self.currentjob]=0
        else:
            self.nxV.set(8)	
            self.nyV.set(1)
            self.nzV.set(16)
            self.nvV.set(32)
            self.nwV.set(8)
            self.kyminV.set(0.3)
            self.lxV.set(0)
            self.ky0indV.set(1)
            self.nexcV.set(1)
            self.adaptlx.set(1)
            self.adaptlx.set(0)
            if self.nonlocal.get():
                self.lxV.set(60)
                self.nxV.set(32)
                self.initcondV.set("'db'")
            else:
                self.initcondV.set("'alm'")
            self.timelimV.set(3600)
            self.dtmaxV.set(0.0385)
            self.calcdtV.set(1)
            self.job_saved[self.currentjob]=0
        launcher.adaptlx_sel
        launcher.calcdtclick()
        launcher.switch_nonlocal(1)
        launcher.switch_lin_nonlin()
        launcher.sel_geom()
        launcher.autosw()
        launcher.switch_scan()
        launcher.switch_collop()
                
    def clearspec(self,arg):
        if arg==1:
            for item in [self.specname1V,self.omn1V,self.omt1V,self.mass1V,self.charge1V,self.temp1V,self.dens1V]:
                item.set('')
            self.passive1V.set(0)
        else: 
            for item in [self.specname2V,self.omn2V,self.omt2V,self.mass2V,self.charge2V,self.temp2V,self.dens2V]:
                item.set('')
            self.passive2V.set(0)

    def save_spec1(self):
        i=int(self.currspec1V.get())-1
        #this try-except part tests if the array element for the respective species exists (which it should per construction)
        #and saves the adjustments made; otherwise nothing is done
        try:
            #dummy which returns an error if the element does not exist
            #since the general case, where specname[i] is some string is not treated, it will also 'pass'
            if self.specname[i]==0: pass
        except:
            return
        self.specname[i]=self.specname1V.get()
        self.mass[i]=float(self.mass1V.get())
        self.charge[i]=int(self.charge1V.get())
        self.passive[i]=int(self.passive1V.get())
        self.proftype[i]=int(self.proftypeV.get())
        if int(self.proftypeV.get()) in range(1,6):
#            self.temp[i]=self.temp1V.get()
#            self.dens[i]=self.dens1V.get()
            self.omn[i]=-1
            self.omt[i]=-1
#            self.kappa_n[i]=self.omn1V.get()
#            self.kappa_T[i]=self.omt1V.get()
        else:
            if self.proftypeV.get()==0:
                self.temp[i]=self.temp1V.get()
                self.dens[i]=self.dens1V.get()
                self.omn[i]=self.omn1V.get()
                self.omt[i]=self.omt1V.get()
            else:
                self.temp[i]=-1
                self.dens[i]=-1
                self.omn[i]=-1
                self.omt[i]=-1
            #variables are set to dummy values to avoid problems with multi-job management
            self.kappa_n[i]=-1
            self.kappa_T[i]=-1
            self.LT_center[i]=-1
            self.Ln_center[i]=-1
            self.LT_width[i]=-1
            self.Ln_width[i]=-1
            self.delta_x_T[i]=-1
            self.delta_x_n[i]=-1
            self.src_proftype[i]=-1
            self.src_amp[i]=-1
            self.src_width[i]=-1
            self.src_x0[i]=-1
            
            
    def save_spec2(self):
        i=int(self.currspec2V.get())-1
        try: 
            if self.specname[i]==0: pass
        except:
            return
        self.specname[i]=self.specname2V.get()
        self.mass[i]=self.mass2V.get()
        self.charge[i]=int(self.charge2V.get())
        self.temp[i]=self.temp2V.get()
        self.dens[i]=self.dens2V.get()
        self.passive[i]=int(self.passive2V.get())
        self.proftype[i]=int(self.proftypeV.get())
        if self.proftypeV.get() in range(1,6):
#            self.temp[i]=self.temp2V.get()
#            self.dens[i]=self.dens2V.get()
            self.omn[i]=-1
            self.omt[i]=-1
#            self.kappa_n[i]=self.omn2V.get()
#            self.kappa_T[i]=self.omt2V.get()
        else:
            if int(self.proftypeV.get())==0:
                self.temp[i]=self.temp2V.get()
                self.dens[i]=self.dens2V.get()
                self.omn[i]=self.omn2V.get()
                self.omt[i]=self.omt2V.get()
            else:
                self.temp[i]=-1
                self.dens[i]=-1
                self.omn[i]=-1
                self.omt[i]=-1
            #variables are set to dummy values to avoid problems with multi-job management
            self.kappa_n[i]=-1
            self.kappa_T[i]=-1
            self.LT_center[i]=-1
            self.Ln_center[i]=-1
            self.LT_width[i]=-1
            self.Ln_width[i]=-1
            self.delta_x_T[i]=-1
            self.delta_x_n[i]=-1
            self.src_proftype[i]=-1
            self.src_amp[i]=-1
            self.src_width[i]=-1
            self.src_x0[i]=-1


    def ch_spec1(self,i=0):
        if not i: i=int(self.currspec1V.get())-1
        #transform to array indices
        else: i=i-1
        try: 
            if self.specname[i]: pass
        except: return
        try: self.specname1V.set(self.specname[i])
        except:
            launcher.Msg('Species name '+str(i+1)+' missing.')
            self.specname1V.set("'unnamed'")
        try:
            if self.nonlocal.get():
                self.omn1V.set(self.kappa_n[i])
            else:
                self.omn1V.set(self.omn[i])
        except:
            launcher.Msg('Density gradient '+str(i+1)+' missing.')
            self.omn1V.set(0)
        try:
            if self.nonlocal.get():
                self.omt1V.set(self.kappa_T[i])
            else:
                self.omt1V.set(self.omt[i])
        except:
            launcher.Msg('Temperature gradient '+str(i+1)+' missing.')
            self.omt1V.set(0)
        try: self.mass1V.set(self.mass[i])
        except:
            launcher.Msg('Mass '+str(i+1)+' missing.')
            self.mass1V.set(1)
        try: self.charge1V.set(self.charge[i])
        except:
            launcher.Msg('Charge '+str(i+1)+' missing.')
            self.charge1V.set(1)
        try: self.temp1V.set(self.temp[i])
        except:
            launcher.Msg('Temperature '+str(i+1)+' missing.')
            self.temp1V.set(1)
        try: self.dens1V.set(self.dens[i])
        except:
            launcher.Msg('Density '+str(i+1)+' missing.')
            self.dens1V.set(1)
        try: self.passive1V.set(self.passive[i])
        except:
            self.passive1V.set(0)
        try:
            self.proftypeV.set(self.proftype[i])
        except:
            self.proftypeV.set(0)
        self.currspec1V.set(i+1)
        return "TRUE"
                    
    def ch_spec2(self,i=0):
        if not i: i=int(self.currspec2V.get())-1
        #transform to array indices
        else: i=i-1
        try: 
            if self.specname[i]: pass
        except: 
            return
        try: self.specname2V.set(self.specname[i])
        except:
            launcher.Msg('Species name '+str(i+1)+' missing.')
            self.specname2V.set("'unnamed'")
        try:
            if self.nonlocal.get():
                self.omn2V.set(self.kappa_n[i])
            else:
                self.omn2V.set(self.omn[i])
        except:
            launcher.Msg('Density gradient '+str(i+1)+' missing.')
            self.omn2V.set(0)
        try:
            if self.nonlocal.get():
                self.omt2V.set(self.kappa_T[i])
            else:
                self.omt2V.set(self.omt[i])
        except:
            launcher.Msg('Temperature gradient '+str(i+1)+' missing.')
            self.omt2V.set(0)
        try: self.mass2V.set(self.mass[i])
        except:
            launcher.Msg('Mass '+str(i+1)+' missing.')
            self.mass2V.set(1)
        try: self.charge2V.set(self.charge[i])
        except:
            launcher.Msg('Charge '+str(i+1)+' missing.')
            self.charge2V.set(1)
        try: self.temp2V.set(self.temp[i])
        except:
            launcher.Msg('Temperature '+str(i+1)+' missing.')
            self.temp2V.set(1)
        try: self.dens2V.set(self.dens[i])
        except:
            launcher.Msg('Density '+str(i+1)+' missing.')
            self.dens2V.set(1)
        try: self.passive2V.set(self.passive[i])
        except:
            self.passive2V.set(0)

        self.currspec2V.set(i+1)
        return "TRUE"
                
    def newspec(self,case=None):
        if not case: self.save_spec1()
        if int(self.numspecV.get())>=2: self.save_spec2()
        if case==None: 
            self.specname.append("'unnamed'")
            self.omn.append('0.')
            self.omt.append('0.')
            self.mass.append('1.')
            self.charge.append(1)
            self.temp.append('1.')
            self.dens.append('1.')
            self.passive.append(0)
            self.kappa_n.append('0.')
            self.kappa_T.append('0.')
            self.LT_center.append('0.')
            self.Ln_center.append('0.')
            self.Ln_width.append('0.')
            self.LT_width.append('0.')
            self.delta_x_T.append('0.')
            self.delta_x_n.append('0.')
            self.proftype.append(0)
            self.src_proftype.append(0)
            self.src_amp.append('0.')
            self.src_width.append('0.')
            self.src_x0.append('0.')
            if int(self.numspecV.get())>=1:
                self.currspec2V.set(len(self.specname))
            else: self.currspec1V.set(len(self.specname))
        elif case=='CBC':
            if len(self.specname)==0:
                self.specname.append("'Ions'")
                self.omn.append('2.22')
                self.omt.append('6.92')
                self.mass.append('1.')
                self.charge.append(1)
                self.temp.append('1.')
                self.dens.append('1.')
                self.passive.append(0)
                if not self.nonlocal.get():
                    self.kappa_n.append('-1.')
                    self.kappa_T.append('-1.')
                    self.LT_center.append('-1.')
                    self.Ln_center.append('-1.')
                    self.LT_width.append('-1.')
                    self.Ln_width.append('-1.')
                    self.delta_x_T.append('-1.')
                    self.delta_x_n.append('-1.')
                    self.proftype.append(0)
                    self.src_proftype.append(0)
                    self.src_amp.append('-1.')
                    self.src_width.append('-1.')
                    self.src_x0.append('-1.')
                else:
                    self.kappa_n.append('2.22')
                    self.kappa_T.append('6.92')
                    self.LT_center.append('0.5')
                    self.Ln_center.append('0.5')
                    self.LT_width.append('0.4')
                    self.Ln_width.append('0.4')
                    self.delta_x_T.append('0.')
                    self.delta_x_n.append('0.')
                    self.proftype.append(2)
                    self.src_proftype.append(0)
                    self.src_amp.append('0.')
                    self.src_width.append('0.')
                    self.src_x0.append('0.')
            elif len(self.specname)==1:
                self.specname.append("'Electrons'")
                self.omn.append('2.22')
                self.omt.append('6.92')
                self.mass.append('2.725e-4')
                self.charge.append(-1)
                self.temp.append('1.')
                self.dens.append('1.')
                self.passive.append(0)
                if not self.nonlocal.get():
                    self.kappa_n.append('-1.')
                    self.kappa_T.append('-1.')
                    self.LT_center.append('-1.')
                    self.Ln_center.append('-1.')
                    self.LT_width.append('-1.')
                    self.Ln_width.append('-1.')
                    self.delta_x_T.append('-1.')
                    self.delta_x_n.append('-1.')
                    self.proftype.append(0)
                    self.src_proftype.append(0)
                    self.src_amp.append('-1.')
                    self.src_width.append('-1.')
                    self.src_x0.append('-1.')
                else:
                    self.kappa_n.append('2.22')
                    self.kappa_T.append('6.92')
                    self.LT_center.append('0.5')
                    self.Ln_center.append('0.5')
                    self.LT_width.append('0.4')
                    self.Ln_width.append('0.4')
                    self.delta_x_T.append('0.')
                    self.delta_x_n.append('0.')
                    self.proftype.append(2)
                    self.src_proftype.append(0)
                    self.src_amp.append('0.')
                    self.src_width.append('0.')
                    self.src_x0.append('0.')
            else:
                self.specname.append("'Spec"+str(len(self.specname)+1)+"'")
                self.omn.append('2.22')
                self.omt.append('6.92')
                self.mass.append('1.')
                self.charge.append(-1)
                self.temp.append('1.')
                self.dens.append('0.')
                self.passive.append(0)
                if not self.nonlocal.get():
                    self.kappa_n.append('-1.')
                    self.kappa_T.append('-1.')
                    self.LT_center.append('-1.')
                    self.Ln_center.append('-1.')
                    self.LT_width.append('-1.')
                    self.Ln_width.append('-1.')
                    self.delta_x_T.append('-1.')
                    self.delta_x_n.append('-1.')
                    self.proftype.append(0)
                    self.src_proftype.append(0)
                    self.src_amp.append('-1.')
                    self.src_width.append('-1.')
                    self.src_x0.append('-1.')
                else:
                    self.kappa_n.append('2.22')
                    self.kappa_T.append('6.92')
                    self.LT_center.append('0.5')
                    self.Ln_center.append('0.5')
                    self.LT_width.append('0.4')
                    self.Ln_width.append('0.4')
                    self.delta_x_T.append('0.')
                    self.delta_x_n.append('0.')
                    self.proftype.append(2)
                    self.src_proftype.append(0)
                    self.src_amp.append('0.')
                    self.src_width.append('0.')
                    self.src_x0.append('0.')
        self.numspecV.set(int(self.numspecV.get())+1)
        self.ch_spec1()
        self.ch_spec2()
