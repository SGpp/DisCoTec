from Tkinter import *
from Tooltips import *
import re
from Pmw import Group

class Domain:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Domframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.GridG = Group(self.Domframe,tag_text='Number of grid points',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                           tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.GridG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.nxL=Label(self.GridG.interior(),bg=vars.stdbg,text="Radial (x):")
        self.nyL=Label(self.GridG.interior(),bg=vars.stdbg,text="Binormal (y):")
        self.nzL=Label(self.GridG.interior(),bg=vars.stdbg,text="Parallel (z):")
        self.nvL=Label(self.GridG.interior(),bg=vars.stdbg,text="Par. velocity (v_par):")
        self.nwL=Label(self.GridG.interior(),bg=vars.stdbg,text="Perp. velocity (mu):")
        self.nxE=Entry(self.GridG.interior(),bg="white",text=vars.nxV,width=15)#,justify="right")
        self.nyE=Entry(self.GridG.interior(),bg="white",text=vars.nyV,width=15)#,justify="right")
        self.nzE=Entry(self.GridG.interior(),bg="white",text=vars.nzV,width=15)#,justify="right")
        self.nvE=Entry(self.GridG.interior(),bg="white",text=vars.nvV,width=15)#,justify="right")
        self.nwE=Entry(self.GridG.interior(),bg="white",text=vars.nwV,width=15)#,justify="right")

        self.BoxG = Group(self.Domframe,tag_text='Box properties',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                          tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.BoxG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        
        self.kyminL=Label(self.BoxG.interior(),bg=vars.stdbg,text="Minimum ky:")
        self.kyminE=Entry(self.BoxG.interior(),bg="white",width=20,text=vars.kyminV)
        self.lxL=Label(self.BoxG.interior(),bg=vars.stdbg,text="Radial extent (lx):")
        self.lxE=Entry(self.BoxG.interior(),bg="white",width=7,textvariable=vars.lxV,state="readonly")
        vars.lxV.trace("w",self.radboxsize)
        self.lxButton=Button(self.BoxG.interior(),bg=vars.stdbg,text='Calculate',command=self.calc_lx_nexc,height=1)
        self.nexcE=Entry(self.BoxG.interior(),bg="white",width=20,textvariable=vars.nexcV)
        vars.nexcV.trace("w",self.radboxsize)
        self.lxinfo=Button(self.BoxG.interior(),bitmap='question',bg=vars.stdbg,command=self.lxhelp)
                
        #Expert part of domain frame
        self.exDomframe=Frame(self.BoxG.interior(),bg=vars.stdbg)
        self.lvL=Label(self.BoxG.interior(),bg=vars.stdbg,width=22,text="v_par grid extent (lv):",anchor='w')
        self.lwL=Label(self.BoxG.interior(),bg=vars.stdbg,text="mu grid extent (lw):")
        self.ky0indL=Label(self.BoxG.interior(),bg=vars.stdbg,text="Index of lowest ky:")
        self.nexcL=Label(self.BoxG.interior(),bg=vars.stdbg,text="Rad. box factor (nexc):")
        self.adaptlxC=Checkbutton(self.BoxG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Adapt radial size",
                                  variable=vars.adaptlx,command=self.adaptlx_sel)
        self.adaptlyC=Checkbutton(self.BoxG.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Force integer toroidal mode numbers",
                                  variable=vars.adaptly)#,command=self.adaptly_sel)
        self.kxcentL=Label(self.BoxG.interior(),bg=vars.stdbg,text="kx_center:")
        
        self.lvE=Entry(self.BoxG.interior(),bg="white",width=20,text=vars.lvV)
        self.lwE=Entry(self.BoxG.interior(),bg="white",width=20,text=vars.lwV)
        self.ky0indE=Entry(self.BoxG.interior(),bg="white",width=20,text=vars.ky0indV)
        self.kxcentE=Entry(self.BoxG.interior(),bg="white",width=20,text=vars.kxcentV)

        self.advanced=[self.ky0indL,self.ky0indE,self.kxcentL,self.kxcentE,self.adaptlxC,self.adaptlyC]
        
    def display(self):
        self.GridG.grid(padx=15,sticky='news')
        self.nxL.grid(sticky='w')
        self.nyL.grid(sticky='w')
        self.nzL.grid(sticky='w')
        self.nvL.grid(sticky='w')
        self.nwL.grid(sticky='w')
        self.nxE.grid(column=1,row=0,sticky='w')
        self.nyE.grid(column=1,row=1,sticky='w')
        self.nzE.grid(column=1,row=2,sticky='w')
        self.nvE.grid(column=1,row=3,sticky='w')
        self.nwE.grid(column=1,row=4,sticky='w')

        self.BoxG.grid(row=1,padx=15,pady=15,sticky='news')
        
        self.kyminL.grid(sticky='w')
        self.nexcL.grid(sticky='w')
        self.lxL.grid(sticky='w')
        self.lxinfo.grid(sticky='e',row=1,rowspan=2,column=0)
        self.kyminE.grid(row=0,sticky='w',column=1,columnspan=2)
        self.nexcE.grid(row=1,sticky='w',column=1,columnspan=2)
        self.lxE.grid(row=2,sticky='w',column=1,columnspan=2)
        self.lxButton.grid(row=2,sticky='e',column=2)
        self.ky0indL.grid(sticky='w')
        self.kxcentL.grid(sticky='w')
        self.adaptlxC.grid(sticky='w')
        self.adaptlyC.grid(sticky='w')
        #self.exDomframe.grid(sticky='w',columnspan=4)
        
        Label(self.BoxG.interior(),bg=vars.stdbg,text='').grid(sticky='news')
        self.lvL.grid(sticky='w')
        self.lwL.grid(sticky='w')
        self.ky0indE.grid(row=3,sticky='w',column=1,columnspan=2)
        self.kxcentE.grid(row=4,sticky='w',column=1,columnspan=2)
        
        self.lvE.grid(row=8,sticky='w',column=1,columnspan=2)
        self.lwE.grid(row=9,sticky='w',column=1,columnspan=2)
        
    def switch_advanced(self):
        if not vars.expertV.get():
            [item.config(state=DISABLED) for item in self.advanced]
        else:
            [item.config(state=NORMAL) for item in self.advanced]

    def set_tooltips(self):
        ToolTip(self.nxL, msg='nx0: For local runs, this should be a power of 2 if possible.', follow=vars.follow, delay=vars.delay)
        ToolTip(self.nyL, msg='nky0: Should be a power of 2 if possible', follow=vars.follow, delay=vars.delay)
        ToolTip(self.nzL, msg='nz0', follow=vars.follow, delay=vars.delay)
        ToolTip(self.nvL, msg='nv0', follow=vars.follow, delay=vars.delay)
        ToolTip(self.nwL, msg='nw0', follow=vars.follow, delay=vars.delay)
        ToolTip(self.kyminL, msg='kymin: All ky modes in the simulation will be multiples of this value.',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.lxL, msg="""lx: In local simulations, this will be calculated from the selected nexc value (dependent on magnetic shear and minimum ky wavenumber). In global simulations, lx can be freely chosen, but must be smaller than 1/rho*.""", follow=vars.follow, delay=vars.delay)
        
        ToolTip(self.lvL, msg='lv: Gives the extent of the parallel velocity space in units of thermal velocities.',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.lwL, msg='lw: Gives the extent of the perpendicular velocity space in units of squared thermal velocities.',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.ky0indL, msg="ky0_ind: Can be used in linear simulations to save computing time by not evaluating the zero mode. Usually set to 1; ignored in nonlinear simulations.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.nexcL, msg="nexc (integer): Should be set to 1 for linear simulations; in nonlinear simulations, it must be selected to give a reasonable lx which accomodates the turbulent structures.\nThe number of necessary radial grid points increases if a larger nexc is chosen.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.adaptlxC, msg="adapt_lx=.T./.F.: Adapts radial box size so that nexc=1, which means that all kx modes are linearly connected to each other.",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.adaptlyC, msg="adapt_ly=.T./.F.: Adapts kymin such that each ky mode in the system corresponds to an integer toroidal mode number.",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.kxcentL, msg="kx_center: Shift the central kx mode to a finite value (only possible in linear simulations).",
                follow=vars.follow, delay=vars.delay)
        

    def switch_nonlocal(self,weak=0):
        if vars.nonlocal.get()==1:
            self.nexcE.config(state='disabled')
            self.lxButton.config(state='disabled')
            self.lxE.config(state='normal')
            vars.nexcV.set('')
            vars.adaptlx.set(0)
            vars.adaptlx.set(1)
        else:
            self.nexcE.config(state='normal')
            self.lxButton.config(state='normal')
            self.lxE.config(state='disabled')
        self.adaptlx_sel()
        launcher.switch_lin_nonlin()

    def adaptlx_sel(self):
        if vars.adaptlx.get()==1:
            vars.nexcV.set(1)
            self.lxE.config(state=DISABLED)
            self.nexcE.config(state=DISABLED)
        if vars.adaptlx.get()==0:
            if not vars.nonlocal.get():
                self.lxE.config(state=NORMAL)
                self.nexcE.config(state=NORMAL)
                self.lxButton.config(state=NORMAL)
            else:
                self.lxE.config(state=NORMAL)
                self.nexcE.config(state=DISABLED)
                self.lxButton.config(state=DISABLED)


    def calc_lx_nexc(self):
        def lx_err(var):
            launcher.Msg('Unable to calculate lx. Check '+var+' value.')
            return "FALSE"
        if vars.geomtype.get()=="'tracer'" and vars.geomread==0:
            launcher.Msg("Please read Tracer file before lx calculation.")
            return "FALSE"
        if vars.form_cleared==0:
            kyre=re.compile(r'(\w*.?\w*)\s*!scan:\s*\w*\s*')#,\s*\w*\s*,\s*\w*\s*')
            tmp=kyre.search(vars.kyminV.get())
            tmp1=kyre.search(vars.shatV.get())
            tmp2=kyre.search(vars.nexcV.get())
            shat=-1.
            kymin=-1
            nexc=-1
            if tmp: 
                #Check whether kymin is present
                try: kymin=float(tmp.group(1))
                #else: error
                except: lx_err('kymin')
                #The same for the case without scan
            else: 
                try: kymin=float(vars.kyminV.get())
                except: lx_err('kymin')
            if tmp1: 
                try: shat=float(tmp1.group(1))
                except: 
                    lx_err('shear')
            else: 
                try: shat=float(vars.shatV.get())
                except: 
                    lx_err('shear')
            if tmp2: 
                try: nexc=float(tmp2.group(1))
                except: lx_err('nexc')
            else: 
                try: nexc=float(vars.nexcV.get())
                except: lx_err('nexc')
            if nexc<>-1. and shat <>-1. and kymin<>-1: 
                if abs(shat)>=1e-5: 
                    vars.lxV.set(nexc/abs(shat)/kymin)
                    launcher.Msg("New lx: "+vars.lxV.get())
                else: launcher.Msg("lx can be chosen arbitrarily in unsheared geometry.")
        return "TRUE" 


    # calling this routine with weak=1 does not force the recommended adaptlx, calcdt and initial condition settings
    def switch_lin_nonlin(self,weak=0):
        if vars.nl.get()==1:
            self.ky0indE.config(state=DISABLED)
            self.kxcentE.config(state=DISABLED)
            if not weak: 
                vars.adaptlx.set(0)
                self.adaptlx_sel()
                self.adaptlxC.config(state=DISABLED)
        else:
            self.ky0indE.config(state=NORMAL)
            if not vars.nonlocal.get(): self.kxcentE.config(state=NORMAL)
            if not weak: 
                vars.ky0indV.set(1)
                if not vars.nonlocal.get(): 
                    vars.adaptlx.set(1)
                    self.adaptlxC.config(state=NORMAL)
                self.adaptlx_sel()
  
    def radboxsize(self,var,ind,mode):
        if vars.nonlocal.get():
            vars.write_lx=1
            vars.write_nexc=0
        else:
            #if lx is being set, deactivate nexc field
            if var==vars.lxV._name:
                if vars.lxV.get():
                    self.nexcE.config(state=DISABLED)
                    vars.nexcV.set('')
                    vars.write_lx=1
                    vars.write_nexc=0
                else:
                    self.nexcE.config(state=NORMAL)
            if var==vars.nexcV._name:
                if vars.nexcV.get():
                    self.lxE.config(state=DISABLED)
                    vars.lxV.set('')
                    vars.write_lx=0
                    vars.write_nexc=1
                else:
                    self.lxE.config(state=NORMAL)
        self.adaptlx_sel()
        self.sel_geom()
                    
    def lxhelp(self):
        helpwin=Toplevel(bg=vars.stdbg,width=250)
        helpwin.title('Info on radial box size settings')
        Label(helpwin,text="""In local simulations, the periodic boundary condition employed in radial 
direction imposes a quantization of the box size. The requirement is that 
nexc=shat*kymin*Lx is an integer number (shat being the magnetic shear). 

The radial box size can be specified to GENE both by directly setting 
nexc to an integer number, or by setting a desired Lx, which GENE will 
then adjust to the nearest possible value with integer nexc. 

The launcher treats both parameters in a mutually exclusive way, and will 
write only one of the two to the parameters file. 

For global simulations, nexc may take arbitrary values, and only Lx can be 
specified.""",justify=LEFT,bg=vars.stdbg).grid(padx=15,pady=15) 
          
    def sel_geom(self):
        if vars.geomtype.get() in ["'chease'","'tracer_efit'","'gist'"]:
            self.lxButton.config(state=DISABLED)
        else:
            self.lxButton.config(state=NORMAL)
