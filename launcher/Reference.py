from Tkinter import *
from Tooltips import *
import math
from Pmw import Group

class Reference:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Refframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()


    def define_widgets(self):
        self.RefG = Group(self.Refframe, tag_text='Reference plasma parameters',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.RefG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.TrefL=Label(self.RefG.interior(),bg=vars.stdbg,text="Ref. Temperature (keV):")
        self.nrefL=Label(self.RefG.interior(),bg=vars.stdbg,text="Ref. Density (1e19/m^3):")
        self.LrefL=Label(self.RefG.interior(),bg=vars.stdbg,text="Ref. Length (m):")
        self.BrefL=Label(self.RefG.interior(),bg=vars.stdbg,text="Ref. Magnetic Field (T):")
        self.mrefL=Label(self.RefG.interior(),bg=vars.stdbg,text="Ref. Mass (m_proton):")
                        
        self.TrefE=Entry(self.RefG.interior(),bg="white",width=20,text=vars.TrefV)
        self.nrefE=Entry(self.RefG.interior(),bg="white",width=20,text=vars.nrefV)
        self.LrefE=Entry(self.RefG.interior(),bg="white",width=20,text=vars.LrefV)
        self.BrefE=Entry(self.RefG.interior(),bg="white",width=20,text=vars.BrefV)
        self.mrefE=Entry(self.RefG.interior(),bg="white",width=20,text=vars.mrefV)
        self.readrefsB=Button(self.RefG.interior(),bg=vars.stdbg,width=12,text='From Tracer',command=launcher.Read_Tracer_r)
        
        self.CalcdimG = Group(self.Refframe, tag_text='Compute dimensionless parameters',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.CalcdimG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.calcbetaB=Button(self.CalcdimG.interior(),bg=vars.stdbg,width=15,text='-> Beta',command=self.calcbeta)
        self.calccollB=Button(self.CalcdimG.interior(),bg=vars.stdbg,width=15,text='-> Collisionality',command=self.calccoll)
        self.calcdebyeB=Button(self.CalcdimG.interior(),bg=vars.stdbg,width=15,text='-> Debye length',command=self.calcdebye)
        self.calcrhostarB=Button(self.CalcdimG.interior(),bg=vars.stdbg,width=15,text='-> rhostar',command=self.calcrhostar)
        self.betaRefE=Entry(self.CalcdimG.interior(),bg="white",width=10,text=vars.betaRefV,state=DISABLED)
        self.collRefE=Entry(self.CalcdimG.interior(),bg="white",width=10,text=vars.collRefV,state=DISABLED)
        self.debyeRefE=Entry(self.CalcdimG.interior(),bg="white",width=10,text=vars.debyeRefV,state=DISABLED)
        self.rhostarRefE=Entry(self.CalcdimG.interior(),bg="white",width=10,text=vars.rhostarV,state=DISABLED)


        
    def display(self):
        self.RefG.grid(padx=15,ipadx=15,sticky='news')
        self.TrefL.grid(sticky='w')
        self.nrefL.grid(sticky='w')
        self.LrefL.grid(sticky='w')
        self.BrefL.grid(sticky='w')
        self.mrefL.grid(sticky='w')
        self.TrefE.grid(row=0,column=1)
        self.nrefE.grid(row=1,column=1)
        self.LrefE.grid(row=2,column=1)
        self.BrefE.grid(row=3,column=1)
        self.mrefE.grid(row=4,column=1)
        Label(self.RefG.interior(),bg=vars.stdbg).grid()
        self.readrefsB.grid(row=6,column=0,sticky='we')

        self.CalcdimG.grid(padx=15,pady=15,ipadx=15,sticky='news')
        self.calcbetaB.grid(row=0,column=0,sticky='w')
        self.calccollB.grid(row=1,column=0,sticky='w')
        self.calcdebyeB.grid(row=2,column=0,sticky='w')
        self.calcrhostarB.grid(row=3,column=0,sticky='w')
        self.betaRefE.grid(row=0,padx=15,column=1,sticky='e')
        self.collRefE.grid(row=1,padx=15,column=1,sticky='e')
        self.debyeRefE.grid(row=2,padx=15,column=1,sticky='e')
        self.rhostarRefE.grid(row=3,padx=15,column=1,sticky='e')
                        
                        

    def set_tooltips(self):
        ToolTip(self.TrefL, msg="Tref: Temperature which will correspond to T = 1.0 in the code", 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.nrefL, msg="nref: Density which will correspond to n = 1.0 in the code", 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.LrefL, msg="Lref: Macroscopic normalization length (e.g. minor radius)", 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.BrefL, msg="Bref: Magnetic field at position of the magnetic axis", 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.mrefL, msg="Reference mass (in units of proton mass)", 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.readrefsB, msg="Read only reference values from Tracer file", 
                follow=vars.follow, delay=vars.delay)

        ToolTip(self.calcbetaB, msg="Compute beta from reference temperature, density and magnetic field.\nWarning: This will overwrite any previous setting!", follow=vars.follow, delay=vars.delay)
        ToolTip(self.calccollB, msg="Compute collisionality from reference parameters.\nWarning: This will overwrite any previous setting!", follow=vars.follow, delay=vars.delay)
        ToolTip(self.calcdebyeB, msg="Compute squared Debye length from reference parameters.\nWarning: This will overwrite any previous setting!", follow=vars.follow, delay=vars.delay)
        ToolTip(self.calcrhostarB, msg="Compute ratio of gyroradius to minor radius from reference parameters.\nWarning: This will overwrite any previous setting!", follow=vars.follow, delay=vars.delay)

    def calcbeta(self):
        try:
            Bref=float(vars.BrefV.get())
            nref=float(vars.nrefV.get())
            Tref=float(vars.TrefV.get())
            Lref=float(vars.LrefV.get())
            mref=float(vars.mrefV.get())
        except:
            launcher.Msg('Please specify all reference quantities.')
            return FALSE
        vars.betaV.set('%1.4e' %(2.0*1.2566370614E-6/Bref**2*(nref*1e19*1.60217733E-19)*1e3*Tref))
        vars.betaRefV.set(vars.betaV.get())
        #self.betaRefE.config(text='= %.4g' %float(vars.betaRefV.get()))

    def calccoll(self):
        try:
            Bref=float(vars.BrefV.get())
            nref=float(vars.nrefV.get())
            Tref=float(vars.TrefV.get())
            Lref=float(vars.LrefV.get())
            mref=float(vars.mrefV.get())
        except:
            launcher.Msg('Please specify all reference quantities.')
            return FALSE
        vars.collV.set('%1.4e' %(2.3031E-5*(nref)/Tref**2*Lref*(24.-math.log(math.sqrt(nref*1e19*1.0E-6)/Tref*1e-3))))
        vars.collRefV.set(vars.collV.get())
        #self.collRefE.config(text='= %.4g' %float(vars.collRefV.get()))
        

    def calcdebye(self):
        try:
            Bref=float(vars.BrefV.get())
            nref=float(vars.nrefV.get())
            Tref=float(vars.TrefV.get())
            Lref=float(vars.LrefV.get())
            mref=float(vars.mrefV.get())
        except:
            launcher.Msg('Please specify all reference quantities.')
            return FALSE
        vars.debye2V.set('%1.4e' %(8.85419e-12*Bref**2/nref/1e19/mref/1.67262e-27))
        vars.debyeRefV.set(vars.debye2V.get())
        #self.debyeRefE.config(text='= %.4g' %float(vars.debyeRefV.get()))

    def calcrhostar(self):
        try:
            Bref=float(vars.BrefV.get())
            Tref=float(vars.TrefV.get())
            mref=float(vars.mrefV.get())
            Lref=float(vars.LrefV.get())
            minorr=float(vars.minorrV.get())
            launcher.Msg('Assuming minor radius of %s*Lref for rhostar calculation.'% (vars.minorrV.get()))
        except:
            launcher.Msg('Please specify all reference quantities, including minor radius (Geometry tab).')
            return FALSE
        vars.rhostarV.set('%1.4e' %(math.sqrt(Tref*mref*1.e3/1.60217733E-19*1.67262e-27)/Bref/minorr/Lref))
        # old calculation: %(1.02*(Tref*1e3)**0.5/(Bref*1e4)/minorr))
        vars.rhostarRefV.set(vars.rhostarV.get())
        #self.rhostarRefE.config(text='= %.4g' %float(vars.rhostarRefV.get()))
        


    def sel_geom(self):
        if vars.geomtype.get()=="'tracer'":
            self.readrefsB.config(state=NORMAL)
        else:
            self.readrefsB.config(state=DISABLED)
