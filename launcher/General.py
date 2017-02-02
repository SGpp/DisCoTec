from Tkinter import *
from Tooltips import *
from Pmw import Group

class General(object):
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Genframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.TermG = Group(self.Genframe,tag_text='Termination conditions',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.TermG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.ntimeL=Label(self.TermG.interior(),bg=vars.stdbg,text="Timesteps:")
        self.timelimL=Label(self.TermG.interior(),bg=vars.stdbg,text="Time limit [s]:")
        self.simtimelimL=Label(self.TermG.interior(),bg=vars.stdbg,text="Simul. time limit:")
        self.omprL=Label(self.TermG.interior(),bg=vars.stdbg,text="Freq. precision:")
        self.oflL=Label(self.TermG.interior(),bg=vars.stdbg,text="Overflow limit:")
        self.uflL=Label(self.TermG.interior(),bg=vars.stdbg,text="Underflow limit:")
        self.ntimeE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.ntimeV)
        self.timelimE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.timelimV)
        self.simtimelimE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.simtimelimV)
        self.omprE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.omprV)
        self.oflE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.oflV)
        self.uflE=Entry(self.TermG.interior(),bg="white",width=6,text=vars.uflV)
        #self.exGenframe=Frame(self.Genframe,bg=vars.stdbg,height=96)
        #self.stickyframe=Frame(self.Genframe,bg=vars.stdbg,height=96)

        self.MiscG = Group(self.Genframe,tag_text='Miscellaneous',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.MiscG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.pressureC=Checkbutton(self.MiscG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,variable=vars.pressV,
                                   text="Pressure term")
        self.delzonalC=Checkbutton(self.MiscG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,variable=vars.delzonalV,
                                   text="No Zonal Flows")
#        self.quasilintmC=Checkbutton(self.MiscG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,variable=vars.quasilintmV,
#                                   text="Quasilin. transp. model")
        self.FlowshearG = Group(self.Genframe,tag_text='Equilibrium flows',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.FlowshearG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.ExBrateL=Label(self.FlowshearG.interior(),bg=vars.stdbg,text="ExB shear:")
        self.pfsrateL=Label(self.FlowshearG.interior(),bg=vars.stdbg,text="Parallel shear:")
        self.ExBtimeL=Label(self.FlowshearG.interior(),bg=vars.stdbg,text="Start time:")
        self.PhikxL=Label(self.FlowshearG.interior(),bg=vars.stdbg,text="Amplitude:")
        self.kxExBL=Label(self.FlowshearG.interior(),bg=vars.stdbg,text="kx index:")
        self.ExBrateE=Entry(self.FlowshearG.interior(),bg="white",width=6,text=vars.ExBrateV)
        self.pfsrateE=Entry(self.FlowshearG.interior(),bg="white",width=6,text=vars.pfsrateV)
        self.ExBtimeE=Entry(self.FlowshearG.interior(),bg="white",width=6,text=vars.ExBtimeV)
        self.kxExBE=Entry(self.FlowshearG.interior(),bg="white",width=6,text=vars.kxExBV)
        self.PhikxE=Entry(self.FlowshearG.interior(),bg="white",width=6,text=vars.PhikxV)

        self.DimparG = Group(self.Genframe,tag_text='Dimensionless parameters',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.DimparG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.betaL=Label(self.DimparG.interior(),bg=vars.stdbg,text="Beta:")
        self.debye2L=Label(self.DimparG.interior(),bg=vars.stdbg,text="Sq. Debye length:")
        self.betaE=Entry(self.DimparG.interior(),bg="white",width=18,text=vars.betaV)
        self.debye2E=Entry(self.DimparG.interior(),bg="white",width=18,text=vars.debye2V)
                
        self.InitG = Group(self.Genframe,tag_text='Initial value solver',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.InitG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        
        self.dtmaxL=Label(self.InitG.interior(),width=18,anchor='w',bg=vars.stdbg,text="Max. timestep:")
        self.calcdtC=Checkbutton(self.InitG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Auto",
                                 variable=vars.calcdtV,command=self.calcdtclick)
        self.courantL=Label(self.InitG.interior(),bg=vars.stdbg,text="Courant factor:")
        self.initcondL=Label(self.InitG.interior(),bg=vars.stdbg,text="Initial condition:")
        self.dtmaxE=Entry(self.InitG.interior(),bg="white",width=6,text=vars.dtmaxV)
        self.courantE=Entry(self.InitG.interior(),bg="white",width=9,text=vars.courantV)

        self.initcondMB=Menubutton(self.InitG.interior(),textvariable=vars.initcondV,bg="#eee")
        self.initcondM=Menu(self.initcondMB,bg=vars.stdbg)
        self.initcondMB.config(relief='sunken',menu=self.initcondM)
        self.initcondM.add_radiobutton(label="ppj",variable=vars.initcondV,value="'ppj'")
        self.initcondM.add_radiobutton(label="ppjmt",variable=vars.initcondV,value="'ppjmt'")
        self.initcondM.add_radiobutton(label="ppjrn",variable=vars.initcondV,value="'ppjrn'")
        self.initcondM.add_radiobutton(label="ppg",variable=vars.initcondV,value="'ppg'")
        self.initcondM.add_radiobutton(label="mmj",variable=vars.initcondV,value="'mmj'")
        self.initcondM.add_radiobutton(label="alm",variable=vars.initcondV,value="'alm'")
        self.initcondM.add_radiobutton(label="almmt",variable=vars.initcondV,value="'almmt'")
        self.initcondM.add_radiobutton(label="sw",variable=vars.initcondV,value="'sw'")
        self.initcondM.add_radiobutton(label="fb",variable=vars.initcondV,value="'fb'")
        self.initcondM.add_radiobutton(label="db",variable=vars.initcondV,value="'db'")
        self.initcondM.add_radiobutton(label="zero",variable=vars.initcondV,value="'zero'")
        self.initcondM.add_radiobutton(label="cs",variable=vars.initcondV,value="'cs'")
        self.initcondM.add_radiobutton(label="gam",variable=vars.initcondV,value="'gam'")
        self.initcondM.add_radiobutton(label="cosn",variable=vars.initcondV,value="'cosn'")
        self.auxxL=Label(self.InitG.interior(),bg=vars.stdbg,text="aux x:")
        self.auxyL=Label(self.InitG.interior(),bg=vars.stdbg,text="aux y:")
        self.auxzL=Label(self.InitG.interior(),bg=vars.stdbg,text="aux z:")
        self.auxxE=Entry(self.InitG.interior(),bg="white",width=6,text=vars.auxxV)
        self.auxyE=Entry(self.InitG.interior(),bg="white",width=6,text=vars.auxyV)
        self.auxzE=Entry(self.InitG.interior(),bg="white",width=6,text=vars.auxzV)
        
        self.EVgroup = Group(self.Genframe,tag_text='Eigenvalue computations',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.EVgroup.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.nevL=Label(self.EVgroup.interior(),bg=vars.stdbg,text="Number of EVs:")
        self.maxitL=Label(self.EVgroup.interior(),bg=vars.stdbg,text="Max. iterations:")
        self.shiftL=Label(self.EVgroup.interior(),bg=vars.stdbg,text="Target:")
        self.nevE=Entry(self.EVgroup.interior(),bg="white",width=6,text=vars.nevV)
        self.maxitE=Entry(self.EVgroup.interior(),bg="white",width=6,text=vars.maxitV,state=DISABLED)
        self.shiftE=Entry(self.EVgroup.interior(),bg="white",width=6,text=vars.shiftV,state=DISABLED)
        
        self.NCgroup = Group(self.Genframe,tag_text='Neoclassics', hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                             tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.NCgroup.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.f0contribC=Checkbutton(self.NCgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Use neoclassical term",
                                    variable=vars.incf0V)
        self.delphiC=Checkbutton(self.NCgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Remove electrostatic potential",
                                 variable=vars.delphiV)


        self.advanced=[self.initcondL,self.initcondMB,self.omprL,self.omprE,self.oflL,self.oflE,self.uflL,self.uflE,
                         self.pressureC,self.delzonalC,self.auxxL,self.auxyL,self.auxzL,self.auxxE,
                         self.auxyE,self.auxzE,self.courantL,self.courantE]



    def display(self):
        self.TermG.grid(row=0,column=0,padx=15,sticky='news')
        self.ntimeL.grid(sticky='w')
        self.timelimL.grid(sticky='w')
        self.simtimelimL.grid(sticky='w')
        self.omprL.grid(sticky='w')
        self.oflL.grid(sticky='w')
        self.uflL.grid(sticky='w')
        self.ntimeE.grid(row=0,column=1,sticky='ew')
        self.timelimE.grid(row=1,column=1,sticky='ew')
        self.simtimelimE.grid(row=2,column=1,sticky='ew')
        self.omprE.grid(row=3,column=1,sticky='ew')
        self.oflE.grid(row=4,column=1,sticky='ew')
        self.uflE.grid(row=5,column=1,sticky='ew')

        self.InitG.grid(row=0,column=1,columnspan=2,padx=15,sticky='news')
        self.dtmaxL.grid(sticky='w')
        self.dtmaxE.grid(row=0,column=1,sticky='ew')
        self.calcdtC.grid(sticky='w',row=0,column=2)
        self.courantL.grid(sticky='w')
        self.courantE.grid(row=1,column=1,sticky='ew')
        Label(self.InitG.interior(),bg=vars.stdbg).grid()
        self.initcondL.grid(row=3,column=0,sticky='w')
        self.initcondMB.grid(row=3,column=1,columnspan=2,sticky='ew')
        self.auxxL.grid(row=4,column=1,sticky='w')
        self.auxyL.grid(row=5,column=1,sticky='w')
        self.auxzL.grid(row=6,column=1,sticky='w')
        self.auxxE.grid(row=4,column=2,sticky='ew')
        self.auxyE.grid(row=5,column=2,sticky='ew')
        self.auxzE.grid(row=6,column=2,sticky='ew')



        #self.exGenframe.grid(sticky='w',row=5,columnspan=5)
        self.DimparG.grid(row=1,column=0,columnspan=2,padx=15,pady=15,sticky='news')
        self.betaL.grid(row=0,sticky='w')
        self.debye2L.grid(row=2,sticky='w')
        self.betaE.grid(row=0,column=1,sticky='ew')
        self.debye2E.grid(row=2,column=1,sticky='ew')
 
        
        #self.MiscG.grid(row=1,column=2,padx=15,pady=15,sticky='news')
        #self.pressureC.grid(sticky='w')
        #self.delzonalC.grid(row=1,sticky='w')
        #self.quasilintmC.grid(row=2,sticky='w')

        self.FlowshearG.grid(row=1,column=2,pady=15,padx=15,sticky='news')
        self.ExBrateL.grid(sticky='w')
        self.pfsrateL.grid(sticky='w')
        self.ExBtimeL.grid(sticky='w')
        self.kxExBL.grid(sticky='w')
        self.PhikxL.grid(sticky='w')
        self.ExBrateE.grid(row=0,column=1)
        self.pfsrateE.grid(row=1,column=1)
        self.ExBtimeE.grid(row=2,column=1)
        self.kxExBE.grid(row=3,column=1)
        self.PhikxE.grid(row=4,column=1)

        self.EVgroup.grid(row=2,column=0,columnspan=2,padx=15,sticky='news')
        self.nevL.grid(sticky='w')
        self.maxitL.grid(sticky='w')
        self.shiftL.grid(sticky='w')
        self.nevE.grid(row=0,column=1,sticky='ew')
        self.maxitE.grid(row=1,column=1,sticky='ew')
        self.shiftE.grid(row=2,column=1,sticky='ew')

        self.NCgroup.grid(row=2,column=2, padx=15, sticky='news')
        self.f0contribC.grid(sticky='w')
        self.delphiC.grid(sticky='w')
        #        [self.Genframe.grid_rowconfigure(i,pad=15) for i in range(y)]
        
        
    def set_tooltips(self):
        ToolTip(self.ntimeL, msg="ntimesteps: Stop the simulation after a certain number of timesteps",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.timelimL, msg="timelim: Stop the simulation after a certain wall clock time (this does not include checkpoint writing)", follow=vars.follow, delay=vars.delay)
        ToolTip(self.simtimelimL, msg="simtimelim: Time limit in normalized simulation time",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.dtmaxL, msg="dt_max: Maximum timestep for simulation", follow=vars.follow, delay=vars.delay)
        ToolTip(self.calcdtC, msg="calc_dt=.T./.F.: Calculate maximum timestep for simulation using SLEPc (may require a lot of time/memory for nonlinear runs)", follow=vars.follow, delay=vars.delay)
        ToolTip(self.courantL, msg="courant: Factor by which timestep will be reduced adaptively in nonlinear simulations",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.pressureC, msg="pressure_off=.F./.T.: Take into account pressure gradient term",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.delzonalC, msg="delzonal=.F./.T.: Suppress zonal flows in nonlinear simulations",
                follow=vars.follow, delay=vars.delay)
#        ToolTip(self.quasilintmC, msg="quasilinear_tm=.F./.T.: Apply quasilinear transport model to linear run (results are sent to standard output)", follow=vars.follow, delay=vars.delay)
        ToolTip(self.initcondL, msg="init_cond: Initial condition for distribution function",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.initcondMB, msg="""init_cond: Initial condition for distribution function
ppj = Power laws in x,y, Powers of the Jacobian in z
ppjmt = Same as ppj, but with odd parity
ppjrn = Same as ppj, but with randomized phases
ppg = Power laws in x,y, Gaussian in z
mmj = Single modes in x,y, Powers of the Jacobian in z
alm = All x,y modes + Gaussian in z
almmt = Same as alm, but with odd parity
sw = Single mode in x,y,z
fb = Gaussian blob in kx, ky and z
db = Gaussian blob in x, ky and z
zero = Zero perturbation
cs = Periodic current sheet
gam = GAM initialization
cosn = Cosine pattern in x,y plus white noise""",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.betaL, msg="beta: Ratio of kinetic electron pressure to magnetic pressure (can be calculated self-consistently in 'Reference' tab)", follow=vars.follow, delay=vars.delay)
        ToolTip(self.omprL, msg="omega_prec: Determines convergence criterion for linear runs",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.debye2L, msg="debye2: Squared Debye length (can be calculated in 'Reference' tab)",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.nevL, msg="n_ev: Number of eigenvalues to look for", follow=vars.follow, delay=vars.delay)
        ToolTip(self.ExBrateL, msg="ExBrate: ExB shearing rate in units of cref/Lref", follow=vars.follow, delay=vars.delay)
        ToolTip(self.pfsrateL, msg="pfsrate: Parallel flow shearing rate in units of cref/Lref; set to -1 for automatic calculation", follow=vars.follow, delay=vars.delay)
        ToolTip(self.ExBtimeL, msg="ExB_stime: Start ExB shearing at specific simulation time (no restart necessary at large shearing rates)", follow=vars.follow, delay=vars.delay)
        ToolTip(self.kxExBL, msg="kxind_phi_ext: kx value for sinusoidal ExB flow profile", follow=vars.follow, delay=vars.delay)
        ToolTip(self.PhikxL, msg="phi0_ext: Amplitude for sinusoidal ExB flow profile", follow=vars.follow, delay=vars.delay)
        ToolTip(self.maxitL, msg="ev_max_it: Maximum number of iterations for eigenvalue solver",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.shiftL, msg="ev_shift: Eigenvalues closest to this complex value are detected first.",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.f0contribC, msg="include_f0_contr: Add the neoclassical term v_D*grad F_0 to the gyrokinetic equation", 
        follow=vars.follow, delay=vars.delay)
        ToolTip(self.delphiC, msg="del_phi: Set the electrostatic potential to zero. Strongly recommended for neoclassical runs.",
        follow=vars.follow, delay=vars.delay)
        
    # calling this routine with weak=1 does not force the recommended adaptlx, calcdt and initial condition settings
    def switch_lin_nonlin(self,weak=0):
        if not weak: 
            vars.calcdtV.set(1)
            self.calcdtclick()
            if vars.nl.get()==1:
                if 'alm' in vars.initcondV.get(): 
                    vars.initcondV.set("'ppj'")
            else:
                if not 'alm' in vars.initcondV.get():
                    vars.initcondV.set("'alm'")

    def switch_advanced(self):
        if not vars.expertV.get():
            [item.config(state=DISABLED) for item in self.advanced]
        else:
            [item.config(state=NORMAL) for item in self.advanced]

    # here, weak does not force the setting to harmonic etc.
    def switch_ev_iv(self,weak=0):
        if vars.ivev.get()==0:
            if not weak: 
                vars.nevV.set(1)
                vars.calcdtV.set(0)
            self.maxitE.config(state=NORMAL)
            self.shiftE.config(state=NORMAL)
            self.calcdtclick()
        else:
            self.maxitE.config(state=DISABLED)
            self.shiftE.config(state=DISABLED)
            if not weak:
                if vars.ivev.get()==1: vars.nevV.set(1)
                vars.calcdtV.set(1)
                self.calcdtclick()
       
    def calcdtclick(self):
        if vars.calcdtV.get()==1:
            self.dtmaxE.config(state=DISABLED)
        if vars.calcdtV.get()==0:
            self.dtmaxE.config(state=NORMAL)
    
