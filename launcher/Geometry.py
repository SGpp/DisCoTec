from Tkinter import *
from tkFileDialog   import askopenfilename,askdirectory,asksaveasfilename
from Tooltips import *
from Pmw import Group
import re,os

class Geometry:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Geoframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.GeoG = Group(self.Geoframe,tag_text='Geometry',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.GeoG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.GeotypeL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Select geometry:")
        self.geomtypeMB=Menubutton(self.GeoG.interior(),textvariable=vars.geomtype,bg=vars.stdbg)
        self.geomtypeM=Menu(self.geomtypeMB,bg=vars.stdbg)
        self.geomtypeMB.config(relief='sunken',menu=self.geomtypeM,bg="#eee")
        self.geomtypeM.add_radiobutton(label="slab",variable=vars.geomtype,value="'slab'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="s-alpha",variable=vars.geomtype,value="'s_alpha'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="circular",variable=vars.geomtype,value="'circular'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="CHEASE",variable=vars.geomtype,value="'chease'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="GIST",variable=vars.geomtype,value="'gist'",command=self.sel_geom)
#        self.geomtypeM.add_radiobutton(label="TRACER",variable=vars.geomtype,value="'tracer'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="TRACER_EFIT",variable=vars.geomtype,value="'tracer_efit'",command=self.sel_geom)
        self.geomtypeM.add_radiobutton(label="Miller",variable=vars.geomtype,value="'miller'",command=self.sel_geom)
        
        self.minorrL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Minor radius:")
        self.lengthL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Major Radius:")
        self.shatL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Global shear (shat):")
        self.q0L=Label(self.GeoG.interior(),bg=vars.stdbg,text="Safety factor (q0):")
        self.trpepsL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Inv. aspect ratio:")
        self.amhdL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Pressure gradient (amhd):")
        self.geomdirL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Geometry directory:")
        self.geomfileL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Geometry file:")
        self.fluxlabelL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Flux label:")
        self.x0L=Label(self.GeoG.interior(),bg=vars.stdbg,text="Radial position:")
        self.rrefL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Reference length:")
        self.rhostarL=Label(self.GeoG.interior(),bg=vars.stdbg,text="rhostar:")

        self.kappaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Elongation (kappa):")
        self.s_kappaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="--> Derivative (s_kappa):")
        self.deltaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Triangularity (delta):")
        self.s_deltaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="--> Derivative (s_delta):")
        self.zetaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Squareness (zeta):")
        self.s_zetaL=Label(self.GeoG.interior(),bg=vars.stdbg,text="--> Derivative (s_zeta):")
        self.drRL=Label(self.GeoG.interior(),bg=vars.stdbg,text="Flux Surf. Shift (drR):")
        
        self.minorrE=Entry(self.GeoG.interior(),bg="white",width=10,text=vars.minorrV)
        self.lengthE=Entry(self.GeoG.interior(),bg="white",width=10,text=vars.lengthV)
        self.shatE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.shatV)
        self.trpepsE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.trpepsV)
        self.amhdE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.amhdV)
        self.q0E=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.q0V)
        self.geomdirE=Entry(self.GeoG.interior(),bg="white",width=30,
                            text=vars.geomdirV)
        self.geomfileE=Entry(self.GeoG.interior(),bg="white",width=30,
                             text=vars.geomfileV)

        self.kappaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.kappaV)
        self.s_kappaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.s_kappaV)
        self.deltaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.deltaV)
        self.s_deltaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.s_deltaV)
        self.zetaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.zetaV)
        self.s_zetaE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.s_zetaV)
        self.drRE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.drRV)


#        self.readgeoaB=Button(self.GeoG.interior(),text="Read All",bg=vars.stdbg,command=launcher.Read_Tracer_a)
        self.readgeogB=Button(self.GeoG.interior(),text="Read Geometry",bg=vars.stdbg,command=launcher.Read_Tracer_g)
        self.fluxlabelMB=Menubutton(self.GeoG.interior(),textvariable=vars.fluxlabelV,bg=vars.stdbg,)
        self.fluxlabelM=Menu(self.fluxlabelMB,bg=vars.stdbg)
        self.fluxlabelMB.config(relief='sunken',menu=self.fluxlabelM,bg="#eee")
        self.fluxlabelM.add_radiobutton(label="C_psi",variable=vars.fluxlabelV,value="'C_psi'")
        self.fluxlabelM.add_radiobutton(label="rho_toroidal",variable=vars.fluxlabelV,value="'arho_t'")
        self.fluxlabelM.add_radiobutton(label="rho_poloidal",variable=vars.fluxlabelV,value="'arho_p'")
        self.fluxlabelM.add_radiobutton(label="rho_Volume",variable=vars.fluxlabelV,value="'arho_v'")
        self.x0E=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.x0V)
        self.rrefE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.rrefV)
        self.rhostarE=Entry(self.GeoG.interior(),bg="white",width=30,text=vars.rhostarV)
        self.magprofC=Checkbutton(self.GeoG.interior(),bg=vars.stdbg,text="Global geometry",variable=vars.magprofV)
        
        
        
    def display(self):
        self.GeoG.grid(row=0,column=0,padx=15,sticky='news')
        self.GeotypeL.grid(columnspan=3,sticky='w')
        Label(self.GeoG.interior(),bg=vars.stdbg).grid(row=1)
        self.minorrL.grid(columnspan=3,sticky='w')
        self.lengthL.grid(columnspan=3,sticky='w')
        self.q0L.grid(columnspan=3,sticky='w')
        self.shatL.grid(columnspan=3,sticky='w')
        self.trpepsL.grid(columnspan=3,sticky='w')
        self.amhdL.grid(columnspan=3,sticky='w')
        self.kappaL.grid(columnspan=3,sticky='w')
        self.s_kappaL.grid(columnspan=3,sticky='w')
        self.deltaL.grid(columnspan=3,sticky='w')
        self.s_deltaL.grid(columnspan=3,sticky='w')
        self.zetaL.grid(columnspan=3,sticky='w')
        self.s_zetaL.grid(columnspan=3,sticky='w')
        self.drRL.grid(columnspan=3,sticky='w')
        self.geomdirL.grid(columnspan=3,sticky='w')
        self.geomfileL.grid(columnspan=3,sticky='w')
        self.fluxlabelL.grid(columnspan=3,sticky='w')
        self.x0L.grid(columnspan=3,sticky='w')
#        self.rrefL.grid(columnspan=3,sticky='w')
        self.rhostarL.grid(columnspan=3,sticky='w')
        self.magprofC.grid(columnspan=3,sticky='w')
        
        
        self.geomtypeMB.grid(row=0,column=3,sticky='w')
        self.minorrE.grid(row=2,column=3,sticky='w')
        self.lengthE.grid(row=3,column=3,sticky='w')
        self.q0E.grid(row=4,column=3,sticky='w')
        self.shatE.grid(row=5,column=3,sticky='w')
        self.trpepsE.grid(row=6,column=3,sticky='w')
        self.amhdE.grid(row=7,column=3,sticky='w')
        self.kappaE.grid(row=8,column=3,sticky='w')
        self.s_kappaE.grid(row=9,column=3,sticky='w')
        self.deltaE.grid(row=10,column=3,sticky='w')
        self.s_deltaE.grid(row=11,column=3,sticky='w')
        self.zetaE.grid(row=12,column=3,sticky='w')
        self.s_zetaE.grid(row=13,column=3,sticky='w')
        self.drRE.grid(row=14,column=3,sticky='w')
        self.geomdirE.grid(row=15,column=3,columnspan=3,sticky='w')
        self.geomfileE.grid(row=16,column=3,columnspan=3,sticky='w')
        self.fluxlabelMB.grid(row=17,column=3,sticky='w')
        self.x0E.grid(row=18,column=3,sticky='w')
#        self.rrefE.grid(row=13,column=3,sticky='w')
        self.rhostarE.grid(row=19,column=3,sticky='w')
        self.readgeogB.grid(row=20,column=0,sticky='ew')
#        self.readgeoaB.grid(row=21,column=0,sticky='ew')
    def set_tooltips(self):
        ToolTip(self.geomtypeMB, msg='Select geometry type', follow=FALSE, delay=0.5)
        ToolTip(self.minorrL, msg='Minor radius in units of L_ref', follow=FALSE, delay=0.5)
        ToolTip(self.lengthL, msg="""slab: parscale; Parallel box length in units of Lref\n s_alpha/circular: major_R; Major radius in units of Lref""", follow=FALSE, delay=0.5)
        ToolTip(self.shatL, msg='shat: Global magnetic shear value', follow=FALSE, delay=0.5)
        ToolTip(self.q0L, msg="""q0: Safety factor\nIn case of nonlocal simulations, a set of coefficients can be given as follows: <c0>,<c1>,...,<c5>""", follow=FALSE, delay=0.5)
        ToolTip(self.trpepsL, msg='trpeps: Inverse aspect ratio r/R', follow=FALSE, delay=0.5)
        ToolTip(self.amhdL, msg='amhd: Shafranov shift for s_alpha and Miller geometry', follow=FALSE, delay=0.5)
        ToolTip(self.kappaL, msg='kappa: Flux surface elongation', follow=FALSE, delay=0.5)
        ToolTip(self.s_kappaL, msg='s_kappa: Radial derivative parameter of the elongation, defined as r/kappa*dkappa/dr', follow=FALSE, delay=0.5)
        ToolTip(self.deltaL, msg='delta: Triangularity', follow=FALSE, delay=0.5)
        ToolTip(self.s_deltaL, msg='s_delta: Radial derivative parameter of triangularity, defined as r/sqrt(1-delta**2)*ddelta/dr', follow=FALSE, delay=0.5)
        ToolTip(self.zetaL, msg='zeta: Squareness', follow=FALSE, delay=0.5)
        ToolTip(self.s_zetaL, msg='s_zeta: Radial derivative parameter of squareness, defined as r*dzeta/dr', follow=FALSE, delay=0.5)
        ToolTip(self.drRL, msg='drR: Radial derivative of flux surface center', follow=FALSE, delay=0.5)
        ToolTip(self.rhostarL, msg='Ratio rho_ref/minor_r for nonlocal simulations', follow=FALSE, delay=0.5)
#        ToolTip(self.readgeoaB, msg='Read everything from the tracer file, including parameters and reference values',
#                follow=FALSE, delay=0.5)
        ToolTip(self.readgeogB, msg='Tracer: Read only q0, shat, and the number of gridpoints from the tracer file\nChease,Tracer_EFIT: Find/check existence of geometry file.\n', follow=FALSE, delay=0.5)
        ToolTip(self.magprofC, msg="Global runs: Use radially dependent magnetic geometry",follow=vars.follow, delay=vars.delay)
        


    def switch_nonlocal(self,weak=0):
        if vars.nonlocal.get()==1:
            self.q0L.config(text='q coefficients:')
            self.rhostarL.grid()
            self.rhostarE.grid()
            self.magprofC.grid()
            self.minorrL.grid()
            self.minorrE.grid()
            self.shatL.grid_remove()
            self.shatE.grid_remove()
            self.trpepsL.grid_remove()
            self.trpepsE.grid_remove()
        else:
            self.q0L.config(text='Safety factor:')
            self.rhostarL.grid_remove()
            self.rhostarE.grid_remove()
            self.magprofC.grid_remove()
            self.minorrL.grid_remove()
            self.minorrE.grid_remove()
            self.shatL.grid()
            self.shatE.grid()
            self.trpepsL.grid()
            self.trpepsE.grid()
        self.sel_geom()

    def sel_geom(self):
        # define lists with the fields that each geometry requires
        slablist=[self.lengthL,self.lengthE,self.trpepsL,self.trpepsE,self.shatL,self.shatE]
        salphalist=[self.lengthL,self.lengthE,self.trpepsL,self.trpepsE,self.shatL,self.shatE,self.amhdL,
                    self.amhdE,self.q0E,self.q0L]
        circularlist=[self.lengthL,self.lengthE,self.q0E,self.q0L]
        cheaselist=[self.fluxlabelL,self.x0L,self.lengthL,self.fluxlabelMB,self.x0E,
                    self.lengthE,self.geomfileL,self.geomfileE,self.geomdirL,self.geomdirE, self.readgeogB,
                    self.rhostarL,self.rhostarE]
#        tracerlist=[self.shatL,self.shatE,self.q0E,self.q0L,self.geomdirL,self.geomfileL,self.geomdirE,
                    #self.geomfileE,self.readgeoaB,self.readgeogB,self.x0L,self.x0E]
        gistlist=[self.shatL,self.shatE,self.q0E,self.q0L,self.geomdirL,self.geomfileL,self.geomdirE,
                    self.geomfileE,self.readgeogB,self.x0L,self.x0E]
        efitlist=[self.geomdirL,self.geomfileL,self.geomdirE,self.x0L,self.x0E,self.geomfileE,self.readgeogB,
                  self.rhostarL,self.rhostarE]
        millerlist=[self.trpepsL,self.trpepsE,self.shatL,self.shatE,self.amhdL,
                    self.amhdE,self.q0E,self.q0L,self.kappaL,self.kappaE,self.s_kappaL,self.s_kappaE,
                    self.deltaL,self.deltaE,self.s_deltaL,self.s_deltaE,self.zetaL,self.zetaE,
                    self.s_zetaL,self.s_zetaE,self.drRL,self.drRE,self.minorrL,self.minorrE,self.lengthL,self.lengthE]
        if vars.nonlocal.get():
            circularlist+=[self.minorrL,self.minorrE,self.rhostarL,self.x0L,self.x0E,self.rhostarE]
            cheaselist+=[self.minorrL, self.minorrE]
#            efitlist+=[self.rhostarL,self.rhostarE]
            gistlist+=[self.rhostarL,self.rhostarE]
        else: circularlist+=[self.trpepsL,self.trpepsE,self.shatL,self.shatE]
        
        # ---------------------------------------------------------------------
        # remove all elements and keep only those of the selected geometry type
        if vars.geomtype.get()=="'slab'":
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in millerlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in gistlist]
            [item.grid() for item in slablist]
            self.lengthL.config(text="Parallel length scale:")
            self.shatE.config(state=NORMAL)
            self.trpepsE.config(state=NORMAL)
        if vars.geomtype.get()=="'s_alpha'":
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in millerlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in gistlist]
            [item.grid() for item in salphalist]
            self.lengthL.config(text="Major radius:")
            self.shatE.config(state=NORMAL)
            self.trpepsE.config(state=NORMAL)
            self.q0E.config(state=NORMAL)
        if vars.geomtype.get()=="'circular'":
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in millerlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in gistlist]
            [item.grid() for item in circularlist]
            self.lengthL.config(text="Major radius:")
            self.shatE.config(state=NORMAL)
            self.trpepsE.config(state=NORMAL)
            self.q0E.config(state=NORMAL)
        if vars.geomtype.get()=="'miller'":
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in gistlist]
            [item.grid() for item in millerlist]
            self.lengthL.config(text="Major radius:")
            self.shatE.config(state=NORMAL)
            self.trpepsE.config(state=NORMAL)
            self.q0E.config(state=NORMAL)
#            self.shatE.config(state="readonly")
#            self.q0E.config(state="readonly")
        if vars.geomtype.get()=="'gist'":
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in millerlist]
            [item.grid() for item in gistlist]
            self.shatE.config(state="readonly")
            self.q0E.config(state="readonly")
        if vars.geomtype.get()=="'chease'":
            self.lengthL.config(text="Major radius:")
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in efitlist]
            [item.grid_remove() for item in gistlist]
            [item.grid_remove() for item in millerlist]
            [item.grid() for item in cheaselist]
        if vars.geomtype.get()=="'tracer_efit'":
            [item.grid_remove() for item in slablist]
            [item.grid_remove() for item in salphalist]
            [item.grid_remove() for item in circularlist]
            [item.grid_remove() for item in cheaselist]
            [item.grid_remove() for item in gistlist]
            [item.grid_remove() for item in millerlist]
            [item.grid() for item in efitlist]
            vars.minorrV.set(1.0)
            if vars.nonlocal.get(): vars.magprofV.set(1)
        launcher.sel_geom(0)

    def Read_Tracer(self,fetch):
        if vars.geomdirV.get() and vars.geomfileV.get():
            geomdir=re.search(r"'(.*)'",vars.geomdirV.get()).group(1)
            geomfile=re.search(r"'(.*)'",vars.geomfileV.get()).group(1)
        try:
            tracerfile=open(os.path.join(geomdir,geomfile),"r")
        except:
            try:
                geometryfile=askopenfilename()
                tracerfile=open(geometryfile,"r")
            except:
                return
            (tmp1,tmp2)=os.path.split(geometryfile)
            vars.geomdirV.set("'"+tmp1+"'")
            vars.geomfileV.set("'"+tmp2+"'")
        if not tracerfile or vars.geomtype in ["'chease'","'tracer_efit'"]: return
        tracerdict={}
        if fetch in ['spec','all']:
            vars.omn=[]
            vars.omt=[]
            vars.mass=[]
            vars.charge=[]
            vars.temp=[]
            vars.dens=[]
            vars.specname=[]
            vars.passive=[]
            vars.kappa_n=[]
            vars.kappa_T=[]
            vars.proftype=[]
            vars.fromfile=[]
            vars.LT_center=[]
            vars.Ln_center=[]
            vars.LT_width=[]
            vars.Ln_width=[]
            vars.delta_x_T=[]
            vars.delta_x_n=[]
            vars.src_amp=[]
            vars.src_width=[]
            vars.src_x0=[]
            vars.src_proftype=[]
        for line in tracerfile:
            # if re.search(r'\s*!?\w*\s*=.*|\s*&.*',line)==None:
            # Search lines for <parameter> = <value> patterns
            p=re.compile(r'\s*!?(.*?)\s*=\s*(.*)\s*')
            m=p.search(line)
            # Pick matching lines and build parameter dictionary
            # if m: print m.groups()
            if m:
                tracerdict[m.group(1)]=m.group(2)
        if fetch in ['geo','all']:
            if 'q0' in tracerdict: vars.q0V.set(tracerdict['q0'].strip())
            if 'shat' in tracerdict: vars.shatV.set(tracerdict['shat'].strip())
            if 'gridpoints' in tracerdict: vars.gridpoints=(tracerdict['gridpoints'].strip())
            if 'rho' in tracerdict: vars.x0V.set(tracerdict['rho'].strip())
            if 's0' in tracerdict: vars.x0V.set(tracerdict['s0'].strip())
        if fetch in ['spec','all']:
            vars.numspecV.set(2)
            vars.specname=["'Ions'","'Electrons'"]
            vars.charge=['1','-1']
            vars.mass=['1.0','0.0002725']
            vars.dens=['1','1']
            if ('omn_ions') in tracerdict: vars.omn.append(tracerdict['omn_ions'].strip())
            if ('omt_ions') in tracerdict: vars.omt.append(tracerdict['omt_ions'].strip())
            if ('temp_ions') in tracerdict: vars.temp.append(tracerdict['temp_ions'].strip())
            if ('omn_electrons') in tracerdict: vars.omn.append(tracerdict['omn_electrons'].strip())
            if ('omt_electrons') in tracerdict: vars.omt.append(tracerdict['omt_electrons'].strip())
            if ('temp_electrons') in tracerdict: vars.temp.append(tracerdict['temp_electrons'].strip())
            for ind in range(2):
                vars.kappa_n.append(-1)
                vars.kappa_T.append(-1)
                vars.proftype.append(-1)
                vars.LT_center.append(-1)
                vars.Ln_center.append(-1)
                vars.LT_width.append(-1)
                vars.Ln_width.append(-1)
                vars.delta_x_T.append(-1)
                vars.delta_x_n.append(-1)
                vars.src_amp.append(-1)
                vars.src_width.append(-1)
                vars.src_x0.append(-1)
                vars.src_proftype.append(-1)
                vars.fromfile.append(0)
                vars.passive=[0,0]
            vars.ch_spec1(1)
            vars.ch_spec2(2)
            vars.save_spec1()
            vars.save_spec2()
        if fetch in ['params','all']:
            if 'beta' in tracerdict: vars.betaV.set(tracerdict['beta'])
            if 'coll' in tracerdict: vars.collV.set(tracerdict['coll'])
            if 'debye2' in tracerdict: vars.debye2V.set(tracerdict['debye2'])
        if fetch in ['ref','all']:
            if 'Tref' in tracerdict: vars.TrefV.set(tracerdict['Tref'])
            if 'nref' in tracerdict: vars.nrefV.set(tracerdict['nref'])
            if 'Lref' in tracerdict: vars.LrefV.set(tracerdict['Lref'])
            if 'Bref' in tracerdict: vars.BrefV.set(tracerdict['Bref'])
        if fetch in ['geo','all']: vars.geomread=1
        vars.job_saved[vars.currentjob]=0


