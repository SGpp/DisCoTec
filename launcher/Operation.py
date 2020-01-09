from Tkinter import *
from Tooltips import *
from Pmw import Group

class Operation:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Opframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.modegroup = Group(self.Opframe, tag_text='Operation mode',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.modegroup.interior().config(bg=vars.stdbg)
        # Selector for nonlinear/linear mode
        self.CBnonlin=Checkbutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Nonlinear",bg=vars.stdbg,variable=vars.nl,
                                  width=9,anchor="w",command=launcher.switch_lin_nonlin)
        #self.RBlin=Radiobutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Linear",bg=vars.stdbg,variable=vars.nl,value=0,
        #                       width=9,anchor="w",command=launcher.switch_lin_nonlin)
        vars.nl.set(1)

        # Selector for (non)local mode
        self.CBnonloc_x=Checkbutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Nonlocal (x)",bg=vars.stdbg,
                                    variable=vars.nonlocal, width=13,anchor="e",command=launcher.switch_nonlocal)
        #self.RBloc=Radiobutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Local",bg=vars.stdbg,
        #                       variable=vars.nonlocal,value=0,width=24,anchor="w",command=launcher.switch_nonlocal)

        # Selector for parameter scans
        self.CBscan=Checkbutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Parameter scan",bg=vars.stdbg,
                                    variable=vars.scan, width=9,anchor="w", command=launcher.switch_scan)

        # Selectors for initial/eigenvalue solvers
        self.RBiv=Radiobutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Initial value solver",bg=vars.stdbg,
                              variable=vars.ivev, value=1,width=20,anchor="w",command=launcher.switch_ev_iv)
        self.RBev=Radiobutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Eigenvalue solver",bg=vars.stdbg,
                              variable=vars.ivev,value=0,width=20,anchor="w",state=DISABLED,command=launcher.switch_ev_iv)
        self.RBnc=Radiobutton(self.modegroup.interior(),highlightbackground=vars.stdbg,text="Neoclassical equilibrium solver",bg=vars.stdbg,
                              variable=vars.ivev,value=2,width=20,anchor="w",state=DISABLED,command=launcher.switch_ev_iv)  
        self.RBiv.select()


        self.numgroup = Group(self.Opframe, tag_text='Numerical options',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.numgroup.interior().config(bg=vars.stdbg)
#        self.numgroup=Frame(self.Opframe,bg=vars.stdbg)
#        self.Numopt=Label(self.numgroup.interior(),bg=vars.stdbg,text="Numerical options:",
#                          font=("Helvetica",12,'bold'))
            
        self.parschemeL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Par. diff. scheme:")
        self.timeschemeL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Timestepping scheme:")
        
        self.timeschMB=Menubutton(self.numgroup.interior(),textvariable=vars.timeschV,bg=vars.stdbg)
        self.timeschM=Menu(self.timeschMB,bg=vars.stdbg)
        self.timeschMB.config(relief='sunken',bg="#eee",menu=self.timeschM)
        self.timeschM.add_radiobutton(label="Implicit Euler (Scalapack)",variable=vars.timeschV,
                                      value="'IE1s'",state=DISABLED)
        self.timeschM.add_radiobutton(label="Implicit Euler (PETSc)",variable=vars.timeschV,
                                      value="'IE1p'",state=DISABLED)
        self.timeschM.add_radiobutton(label="Runge-Kutta 4th order",variable=vars.timeschV,value="'RK4'")
        self.timeschM.add_radiobutton(label="Runge-Kutta 4th mod.",variable=vars.timeschV,value="'RK4M'")
            
        self.parschMB=Menubutton(self.numgroup.interior(),textvariable=vars.parschV,bg=vars.stdbg)
        self.parschM=Menu(self.parschMB,bg=vars.stdbg)
        self.parschMB.config(relief='sunken',bg="#eee",menu=self.parschM)
        self.parschM.add_radiobutton(label="Centered 4th order",variable=vars.parschV,value="'c4th'")
        self.parschM.add_radiobutton(label="Centered 2nd order",variable=vars.parschV,value="'c2nd'")
        self.parschM.add_radiobutton(label="Upwind 3rd order",variable=vars.parschV,value="'u3rd'")
        #Expert Opframe
        self.hypzL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Par. hyperdiffusion:")
        self.hypzordL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Order:")
        self.hypxL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Radial hyperdiffusion:")
        self.hypxordL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Order:")
        self.hypyL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Binormal hyperdiffusion:")
        self.hypyordL=Label(self.numgroup.interior(),bg=vars.stdbg,text="Order:")
        self.hypvL=Label(self.numgroup.interior(),bg=vars.stdbg,text="V-space hyperdiffusion:")
        
        self.hypzE=Entry(self.numgroup.interior(),bg="white",width=20,text=vars.hypzV)
        self.hypxE=Entry(self.numgroup.interior(),bg="white",width=20,text=vars.hypxV)
        self.hypyE=Entry(self.numgroup.interior(),bg="white",width=20,text=vars.hypyV)
        self.hypvE=Entry(self.numgroup.interior(),bg="white",width=20,text=vars.hypvV)

        self.hypzordE=Entry(self.numgroup.interior(),bg="white",width=3,text=vars.hypzordV)
        self.hypxordE=Entry(self.numgroup.interior(),bg="white",width=3,text=vars.hypxordV)
        self.hypyordE=Entry(self.numgroup.interior(),bg="white",width=3,text=vars.hypyordV)

        self.arakawazvC=Checkbutton(self.numgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,
                                    text="Parallel Arakawa discretization",variable=vars.arakawazvV)


#        self.sL=Label(self.Opframe,bg=vars.stdbg,text="Species:")
#        self.vL=Label(self.Opframe,bg=vars.stdbg,text="Par. velocity (v_par):")
#        self.wL=Label(self.Opframe,bg=vars.stdbg,text="Perp. velocity (mu):")
#        self.xL=Label(self.Opframe,bg=vars.stdbg,text="Radial (x):")
#        self.yL=Label(self.Opframe,bg=vars.stdbg,text="Binormal (y):")
#        self.zL=Label(self.Opframe,bg=vars.stdbg,text="Parallel (z):")
        self.pargroup = Group(self.Opframe, tag_text='Parallelization (MPI)',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.pargroup.interior().config(bg=vars.stdbg)
#        self.pargroup.tag().config(bg=vars.stdbg)
        self.nprocsL=Label(self.pargroup.interior(),bg=vars.stdbg,text="Total:")
        self.sE=Entry(self.pargroup.interior(),bg="white",text=vars.sV,width=3)
        self.vE=Entry(self.pargroup.interior(),bg="white",text=vars.vV,width=3)
        self.wE=Entry(self.pargroup.interior(),bg="white",text=vars.wV,width=3)
        self.xE=Entry(self.pargroup.interior(),bg="white",text=vars.xV,width=3)
        self.yE=Entry(self.pargroup.interior(),bg="white",text=vars.yV,width=3)
        self.zE=Entry(self.pargroup.interior(),bg="white",text=vars.zV,width=3)
        vars.sV.trace("w",self.check_parall)
        vars.vV.trace("w",self.check_parall)
        vars.wV.trace("w",self.check_parall)
        vars.xV.trace("w",self.check_parall)
        vars.yV.trace("w",self.check_parall)
        vars.zV.trace("w",self.check_parall)
        self.nparsimE=Entry(self.pargroup.interior(),bg="white",text=vars.nparsimV,width=3)
        self.nprocsE=Entry(self.pargroup.interior(),bg="white",text=vars.nprocsV,width=5)
        self.sC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Species:",
                            variable=vars.sautoV,command=self.autosw)
        self.vC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Par. velocity (v_par):",
                            variable=vars.vautoV,command=self.autosw)
        self.wC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Perp. velocity (mu):",
                            variable=vars.wautoV,command=self.autosw)
        self.xC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Radial (x):",
                            variable=vars.xautoV,command=self.autosw)
        self.yC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Binormal (y):",
                            variable=vars.yautoV,command=self.autosw)
        self.zC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Parallel (z):",
                            variable=vars.zautoV,command=self.autosw)
        self.nparsimC=Checkbutton(self.pargroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="# par. simulations:",
                            variable=vars.npsautoV,command=self.autosw)


        self.advanced=[self.nparsimC,self.nparsimE]

    def display(self):
        self.Opframe.columnconfigure(0,weight=1)
        self.Opframe.columnconfigure(1,weight=1)
        self.Opframe.columnconfigure(2,weight=1)
        self.modegroup.grid(sticky='nws',padx=15,columnspan=2,ipadx=15)
        self.modegroup.interior().config(padx=15,pady=5)
#        Label(self.Opframe,bg=vars.stdbg,text="Mode:",
#              font=("Helvetica",12,'bold')).grid(sticky='w',columnspan=3)
        self.CBnonlin.grid(sticky='w')
        self.CBnonloc_x.grid(row=0, column=1,sticky='e')
        if not vars.allow_nonlocal:
            self.CBnonloc_x.config(state=DISABLED)
        self.CBscan.grid(columnspan=2,sticky='ew')
        #self.RBlin.grid(sticky='w')
        Label(self.modegroup.interior(),bg=vars.stdbg).grid()
        self.RBiv.grid(columnspan=2,sticky='ew')
        self.RBev.grid(columnspan=2,sticky='ew')
        self.RBnc.grid(columnspan=2,sticky='ew')
        #self.RBloc.grid(sticky='w')

        self.pargroup.grid(row=0,column=2,padx=15,sticky='nes')
        self.pargroup.interior().config(padx=15,pady=5)
        self.sC.grid(row=1,column=0,sticky='w')
        self.vC.grid(row=2,column=0,sticky='w')
        self.wC.grid(row=3,column=0,sticky='w')
        self.xC.grid(row=4,column=0,sticky='w')
        self.yC.grid(row=5,column=0,sticky='w')
        self.zC.grid(row=6,column=0,sticky='w')
        self.nparsimC.grid(row=7,column=0,sticky='w')
        Label(self.pargroup.interior(),bg=vars.stdbg).grid()
        self.nprocsL.grid(row=8,column=0,sticky='w')
        self.sE.grid(column=1,row=1)
        self.vE.grid(column=1,row=2)
        self.wE.grid(column=1,row=3)
        self.xE.grid(column=1,row=4)
        self.yE.grid(column=1,row=5)
        self.zE.grid(column=1,row=6)
        self.nparsimE.grid(column=1,row=7)
        self.nprocsE.grid(column=1,row=8)

        self.numgroup.grid(row=1,sticky='ew',padx=15,pady=15,ipadx=15,columnspan=3)
        self.numgroup.interior().config(padx=15,pady=5)
#        self.timeschemeL.grid(row=0,sticky='w')
#        self.timeschMB.grid(row=0,column=1,sticky='ew')
        #Label(self.numgroup.interior(),bg=vars.stdbg).grid(row=0)
        #self.parschemeL.grid(row=1,sticky='w')
        #self.parschMB.grid(row=1,column=1,sticky='ew')
        self.hypzL.grid(row=2,sticky='w')
        self.hypzordL.grid(sticky='w',row=2,column=3)
        self.hypxL.grid(sticky='w')
        self.hypxordL.grid(sticky='w',row=3,column=3)
        self.hypyL.grid(sticky='w')
        self.hypyordL.grid(sticky='w',row=4,column=3)
        self.hypvL.grid(sticky='w')
#        Label(self.Opframe,bg=vars.stdbg).grid(row=6)
        self.hypzE.grid(row=2,columnspan=2,column=1,sticky='w')
        self.hypxE.grid(row=3,columnspan=2,column=1,sticky='w')
        self.hypyE.grid(row=4,columnspan=2,column=1,sticky='w')
        self.hypzordE.grid(row=2,column=4,sticky='w')
        self.hypxordE.grid(row=3,column=4,sticky='w')
        self.hypyordE.grid(row=4,column=4,sticky='w')
        self.hypvE.grid(row=5,columnspan=2,column=1,sticky='w')
        self.arakawazvC.grid(columnspan=2,sticky='w')



                

    def set_tooltips(self):
        ToolTip(self.CBnonlin, msg='nonlinear=.T.: Check this button to run GENE in nonlinear mode, taking into account mode interactions.',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.CBnonloc_x, msg='x_local=.F.: Run a global simulation, taking radial profile variations into account',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.CBscan, msg='Perform a parameter scan using the scanscript. Checking this box will automatically adapt the submission routine to generate an appropriate submit script.',
                follow=vars.follow, delay=vars.delay)
        #ToolTip(self.RBlin, msg='nonlinear=.F.', follow=vars.follow, delay=vars.delay)
        ToolTip(self.RBiv, msg='Solve (non)linear gyrokinetic equations as an initial value problem', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.RBev, msg='Solve linear gyrokinetic equations as an eigenvalue problem',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.RBnc, msg='Solve linear neoclassical problem (ky=0)', 
        follow=vars.follow, delay=vars.delay)
        #ToolTip(self.RBloc, msg='Use the local approximation with periodic radial boundary conditions',
        #        follow=vars.follow, delay=vars.delay)
        ToolTip(self.parschemeL, msg='parscheme: Method for calculation of parallel derivatives',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.timeschemeL, msg='timescheme', follow=vars.follow, delay=vars.delay)
        ToolTip(self.hypzL, msg='hyp_z: Hyperdiffusion in parallel direction', follow=vars.follow,
                delay=vars.delay)
        ToolTip(self.hypzordL, msg='hyp_z_ord: Order of parallel hyperdiffusion', follow=vars.follow,
                delay=vars.delay)
        ToolTip(self.hypxL, msg='hyp_x: Hyperdiffusion in radial direction', follow=vars.follow, delay=vars.delay)
        ToolTip(self.hypxordL, msg='hyp_x_ord: Order of radial hyperdiffusion', follow=vars.follow, delay=vars.delay)
        ToolTip(self.hypyL, msg='hyp_y: Hyperdiffusion in toroidal direction', follow=vars.follow, delay=vars.delay)
        ToolTip(self.hypyordL, msg='hyp_y_ord: Order of toroidal hyperdiffusion', follow=vars.follow, delay=vars.delay)
        ToolTip(self.hypvL, msg='hyp_v: Hyperdiffusion in parallel velocity direction', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.arakawazvC, msg='arakawa_zv: Energy-conserving discretization of parallel and v_parallel derivatives', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.sC, msg='n_procs_s: Number of processes for species parallelization. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.vC, msg='n_procs_v: Number of processes for parallelization in parallel velocity. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.wC, msg='n_procs_w: Number of processes for parallelization in perpendicular velocity. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.xC, msg='n_procs_x: Number of processes for radial parallelization. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.yC, msg='n_procs_y: Number of processes for binormal parallelization. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.zC, msg='n_procs_z: Number of processes for parallel parallelization. Check this button for automatic parallelization.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.nparsimC, msg='n_parallel_sims: Number of simulations run in parallel (for parameter scans). If this button is checked -- and the parallelization is not fixed --, GENE will perform an efficiency analysis to optimize large parameter scans.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.nprocsL, msg='Total number of processes (= product of all components)', 
                follow=vars.follow, delay=vars.delay)

    def check_parall(self,var,ind,mode):
        fieldlist=[vars.sV,vars.vV,vars.wV,vars.xV,vars.yV,vars.zV]
        cblist=[vars.sautoV,vars.vautoV,vars.wautoV,vars.xautoV,vars.yautoV,vars.zautoV]
        ind=[var==item._name for item in fieldlist].index(1)
        #which variable is the one that changed?
        in_var=fieldlist[ind]
        #associated check button value
        cb_var=cblist[ind]
        try:
            if int(in_var.get())<=0:
                cb_var.set(1)
            else:
                cb_var.set(0)
        except:
            pass
        product=1
        try:
            for item in fieldlist:
                product=product*int(item.get())
        except:
            pass
        if product==0:
            self.nprocsE.config(state=NORMAL)
        else:
            self.nprocsE.config(state=DISABLED)

    def autosw(self):
        if vars.sautoV.get()==1:
            vars.sV.set(0)
            self.sE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.sE.config(state=NORMAL)            
        if vars.vautoV.get()==1:
            vars.vV.set(0)
            self.vE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.vE.config(state=NORMAL)            
        if vars.wautoV.get()==1:
            vars.wV.set(0)
            self.wE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.wE.config(state=NORMAL)            
        if vars.xautoV.get()==1:
            vars.xV.set(0)
            self.xE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.xE.config(state=NORMAL)            
        if vars.yautoV.get()==1:
            vars.yV.set(0)
            self.yE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.yE.config(state=NORMAL)            
        if vars.zautoV.get()==1:
            vars.zV.set(0)
            self.zE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.zE.config(state=NORMAL)            
        if vars.npsautoV.get()==1:
            vars.nparsimV.set(1)
            self.nparsimE.config(state=DISABLED)
            self.nprocsE.config(state=NORMAL)
        else:
            self.nparsimE.config(state=NORMAL)            
        if vars.sautoV.get()==0 and vars.vautoV.get()==0 and vars.wautoV.get()==0 and vars.xautoV.get()==0 and vars.yautoV.get()==0 and vars.zautoV.get()==0:
            self.nprocsE.config(state=DISABLED)

    def switch_lin_nonlin(self,weak=0):
        if vars.nl.get()==1:
            self.RBev.config(state=DISABLED)
            self.RBnc.config(state=DISABLED)
            vars.ivev.set(1)
            launcher.switch_ev_iv(weak)
            self.timeschM.entryconfigure(1,state=DISABLED)
            self.timeschM.entryconfigure(2,state=DISABLED)
            if not weak: 
                vars.timeschV.set("'RK4'")
        else:
            self.RBev.config(state=NORMAL)
            self.RBnc.config(state=NORMAL)
            self.timeschM.entryconfigure(1,state=NORMAL)
            self.timeschM.entryconfigure(2,state=NORMAL)
                
    def switch_ev_iv(self,weak=0):
        if vars.ivev.get()==0:
            self.timeschMB.config(state=DISABLED)
        else:
            self.timeschMB.config(state=NORMAL)
            
    def switch_scan(self):
        if vars.scan.get()==1:
            if vars.expertV.get():
                self.nparsimC.config(state=NORMAL)
            self.nparsimE.config(state=DISABLED)
            vars.npsautoV.set(1)
            #vars.nparsimV.set(1)
        else:
            self.nparsimC.config(state=DISABLED)
            self.nparsimE.config(state=DISABLED)
            vars.npsautoV.set(0)
            vars.nparsimV.set(1)
        
    def switch_advanced(self):
        if not vars.expertV.get():
            [item.config(state=DISABLED) for item in self.advanced]
        else:
            [item.config(state=NORMAL) for item in self.advanced]
            self.switch_scan()
        
