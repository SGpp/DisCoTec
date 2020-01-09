from Tkinter import *
from Tooltips import *
from Plot import SimplePlot
from Pmw import Group
import Pmw
import math
from tkFileDialog   import askopenfilename,askdirectory,asksaveasfilename
from Profiles import profile_visual

class Species:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        
        self.Specframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()
        vars.Specfirst=0

    def define_widgets(self):
        self.SpecG = Group(self.Specframe,tag_text='Species information',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.SpecG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.NewSpecB=Button(self.SpecG.interior(),bg=vars.stdbg,width=17,text='New species',command=vars.newspec)
        self.proftypeL=Label(self.SpecG.interior(),bg=vars.stdbg,text="Profile type:",anchor='w',width=100)
        
        self.proftypeMB=Menubutton(self.SpecG.interior(),width=5,textvariable=vars.proftypeV,bg="#eee")
        self.proftypeM=Menu(self.proftypeMB,bg=vars.stdbg)
        self.proftypeMB.config(relief='sunken',menu=self.proftypeM)
        self.proftypeM.add_radiobutton(label="-3: IterDB (DIII-D format)",variable=vars.proftypeV,value=-3)
        self.proftypeM.add_radiobutton(label="-2: IterDB",variable=vars.proftypeV,value=-2)
        self.proftypeM.add_radiobutton(label="-1: GENE profiles",variable=vars.proftypeV,value=-1)
        self.proftypeM.add_radiobutton(label="0: Local",variable=vars.proftypeV,value=0)
        self.proftypeM.add_radiobutton(label="1: Peaked (flat boundaries)",variable=vars.proftypeV,value=1)
        self.proftypeM.add_radiobutton(label="2: Peaked",variable=vars.proftypeV,value=2)
        self.proftypeM.add_radiobutton(label="3: Flat-top",variable=vars.proftypeV,value=3)
        self.proftypeM.add_radiobutton(label="4: Falchetto",variable=vars.proftypeV,value=4)
        self.proftypeM.add_radiobutton(label="5: Exponential",variable=vars.proftypeV,value=5)
        vars.proftypeV.trace("w",self.update_profgroup)

        self.specnameL=Label(self.SpecG.interior(),bg=vars.stdbg,text="Species name:")
        self.massL=Label(self.SpecG.interior(),bg=vars.stdbg,text="Particle mass:")
        self.chargeL=Label(self.SpecG.interior(),bg=vars.stdbg,text="Charge:")

        self.ProfG = Group(self.Specframe,tag_text='Profile information',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.ProfG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        
        self.tempL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Temperature:")
        self.omtL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Temperature gradient:")
        self.densL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Density:")
        self.omnL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Density gradient:", anchor='w')

        self.fileL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Input file:")
        self.fileE=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.fileV)
        self.fileB=Button(self.ProfG.interior(),bg=vars.stdbg,text="Find",command=self.select_input_file)
        self.timeL=Label(self.ProfG.interior(),bg=vars.stdbg,text="Time:")
        self.timeE=Entry(self.ProfG.interior(),bg="white",width=20,textvariable=vars.iterdbtimeV)

        self.propagateMB=Menubutton(self.ProfG.interior(),width=5,text='Use profile data for...',bg="#eee")
        self.propagateM=Menu(self.propagateMB,bg=vars.stdbg)
        self.propagateMB.config(relief='sunken',menu=self.propagateM)
        self.propagateM.add_checkbutton(label='Ref. Temperature',variable=vars.proftoTrefV)
        self.propagateM.add_checkbutton(label='Ref. Density',variable=vars.proftonrefV)
        self.propagateM.add_checkbutton(label='Plasma beta',variable=vars.proftobetaV)
        self.propagateM.add_checkbutton(label='Collisionality',variable=vars.proftocollV)
        self.propagateM.add_checkbutton(label='Debye length',variable=vars.proftodebyeV)
        self.propagateM.add_checkbutton(label='Relative Gyroradius (rho*)',variable=vars.proftorhostarV)
        self.propagateM.add_checkbutton(label='ExB shear',variable=vars.proftoExBV)
        self.propagateM.add_checkbutton(label='Parallel flow shear',variable=vars.proftopfsV)
        

        self.ShLeft1B=Button(self.SpecG.interior(),bg=vars.stdbg,text="<",command=self.ShLeft1)
        self.ShRight1B=Button(self.SpecG.interior(),bg=vars.stdbg,text=">",command=self.ShRight1)
        self.ShLeft2B=Button(self.SpecG.interior(),bg=vars.stdbg,text="<",command=self.ShLeft2)
        self.ShRight2B=Button(self.SpecG.interior(),bg=vars.stdbg,text=">",command=self.ShRight2)
        self.currspec1E=Label(self.SpecG.interior(),bg=vars.stdbg,width=5,textvariable=vars.currspec1V,justify=CENTER)
        self.specname1E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.specname1V)
        self.mass1E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.mass1V)
        self.charge1E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.charge1V)
        self.passive1C=Checkbutton(self.SpecG.interior(),bg=vars.stdbg,text='Passive advection',anchor=W,variable=vars.passive1V)
        self.delete1B=Button(self.SpecG.interior(),bg=vars.stdbg,text="Delete",command=self.delete1)
        self.temp1E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.temp1V)
        self.omt1E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.omt1V)
        
        self.dens1E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.dens1V)
        self.omn1E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.omn1V)
        
        self.currspec2E=Label(self.SpecG.interior(),bg=vars.stdbg,width=5,textvariable=vars.currspec2V,justify=CENTER)
        self.specname2E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.specname2V)
        self.mass2E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.mass2V)
        self.charge2E=Entry(self.SpecG.interior(),bg="white",width=20,text=vars.charge2V)
        self.passive2C=Checkbutton(self.SpecG.interior(),bg=vars.stdbg,text='Passive advection',anchor=W,variable=vars.passive2V)
        self.delete2B=Button(self.SpecG.interior(),bg=vars.stdbg,text="Delete",command=self.delete2)
        self.temp2E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.temp2V)
        self.omt2E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.omt2V)
        self.dens2E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.dens2V)
        self.omn2E=Entry(self.ProfG.interior(),bg="white",width=20,text=vars.omn2V)
        
        
        self.numspecE=Entry(self.SpecG.interior(),bg="white",width=3,text=vars.numspecV,
                            validate='focusout',vcmd=self.ch_numspec,state=DISABLED)
        self.Visualize_Profiles=Button(self.ProfG.interior(),text='View/Set profiles', bg=vars.stdbg)
        self.Visualize_Profiles.bind("<Button-1>",self.profile_visual)
        


    def display(self):
        self.SpecG.grid(column=0,columnspan=3,padx=15,sticky='news')
        Label(self.SpecG.interior(),bg=vars.stdbg,text="Number of species:").grid(sticky='w',columnspan=2)
        Label(self.SpecG.interior(),bg=vars.stdbg,text=""
              ).grid(sticky='w',columnspan=3)
        self.proftypeL.grid(row=1,columnspan=3,sticky='w')
        self.proftypeMB.grid(row=1,column=3,columnspan=6,sticky='news')
        Label(self.SpecG.interior(),bg=vars.stdbg,text="Current species:"
              ).grid(sticky='w',columnspan=3)
        self.numspecE.grid(row=0,column=3,sticky='w')
        self.NewSpecB.grid(row=0,column=6,columnspan=3,sticky='ew')
        self.specnameL.grid(sticky='w')
        self.massL.grid(sticky='w')
        self.chargeL.grid(sticky='w')

        
        self.tempL.grid(sticky='w')
        self.omtL.grid(sticky='w')
        self.densL.grid(sticky='w')
        self.omnL.grid(sticky='ew')

        self.ProfG.grid(column=0,columnspan=3,padx=15,sticky='news')
        
        self.fileL.grid(row=0,column=2,sticky='news')
        self.fileE.grid(row=0,column=3,columnspan=2,sticky='news')
        self.fileB.grid(row=0,column=5,sticky='news')
        self.timeL.grid(row=0,column=6,sticky='news')
        self.timeE.grid(row=0,column=7,columnspan=2,sticky='news')
        self.propagateMB.grid(row=1,columnspan=4,sticky='news')
        
        self.ShLeft1B.grid(row=2,column=3,sticky='w')
        self.currspec1E.grid(row=2,column=4,sticky='ew')
        self.ShRight1B.grid(row=2,column=5,sticky='e')
        self.specname1E.grid(row=3,columnspan=3,column=3,sticky='ew')
        self.mass1E.grid(row=4,columnspan=3,column=3,sticky='ew')
        self.charge1E.grid(row=5,columnspan=3,column=3,sticky='ew')
        self.passive1C.grid(row=6,columnspan=3,column=3,sticky='ew')
        self.delete1B.grid(row=7,columnspan=3,column=3, sticky='ew')
        
        
        self.temp1E.grid(row=0,columnspan=3,column=3,sticky='ew')
        self.omt1E.grid(row=1,columnspan=3,column=3,sticky='ew')
        self.dens1E.grid(row=2,columnspan=3,column=3,sticky='ew')
        self.omn1E.grid(row=3,columnspan=3,column=3,sticky='ew')

        self.ShLeft2B.grid(row=2,column=6,sticky='w')
        self.currspec2E.grid(row=2,column=7,sticky='ew')
        self.ShRight2B.grid(row=2,column=8,sticky='e')

        self.specname2E.grid(row=3,columnspan=3,column=6,sticky='ew')
        self.mass2E.grid(row=4,columnspan=3,column=6,sticky='ew')
        self.charge2E.grid(row=5,columnspan=3,column=6,sticky='ew')
        self.passive2C.grid(row=6,columnspan=3,column=6,sticky='ew')

        self.delete2B.grid(row=7,columnspan=3,column=6, sticky='ew')

        self.temp2E.grid(row=0,columnspan=3,column=6,sticky='ew')
        self.omt2E.grid(row=1,columnspan=3,column=6,sticky='ew')
        self.dens2E.grid(row=2,columnspan=3,column=6,sticky='ew')
        self.omn2E.grid(row=3,columnspan=3,column=6,sticky='ew')

        Label(self.Specframe,text='',bg=vars.stdbg).grid(row=13)
        self.Visualize_Profiles.grid(row=0,column=0,sticky='news')

    def profile_visual(self,event):
        #open first species 
        spec=vars.currspec1V.get()-1
        profwin=profile_visual(spec,self,vars)
        del profwin
        
    def set_tooltips(self):
        ToolTip(self.chargeL, msg="charge: Particle charge. Electrons must be set to -1; only integer values are allowed.",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.proftypeMB, msg="""Profile type:\n-3=IterDB (DIII-D)\n-2=IterDB\n-1=GENE profiles\n0=flat\n1-5=Analytical profiles with radial variation (global only)""",follow=vars.follow, delay=vars.delay)
        ToolTip(self.specnameL, msg="name: Species name. When setting profile type to -1, this name will also determine the filename of the input files as profiles_<species name>.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.massL, msg="mass: Particle mass in units of mref.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.omtL, msg="omt: Specify temperature gradient length scale as omt=Lref/LT.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.omnL, msg="omn: Specify density gradient length scale as omn=Lref/Ln.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.densL, msg="dens: Specify density in units of nref.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.tempL, msg="temp: Specify temperature in units of Tref.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.proftypeL, msg="""Profile type:\n-3=IterDB (DIII-D)\n-2=IterDB\n-1=GENE profiles\n0=flat\n1-5=Analytical profiles with radial variation (global only)""", follow=vars.follow, delay=vars.delay)
        ToolTip(self.passive1C, msg="Treat species as passive tracers, i.e. without influence on the fields.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.passive2C, msg="Treat species as passive tracers, i.e. without influence on the fields.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.fileL, msg="Filename for profile input from IterDB files. If profile type is set to -1, the filename will be assumed to be profiles_<species name>", follow=vars.follow, delay=vars.delay)
        ToolTip(self.fileE, msg="Filename for profile input from IterDB files. If profile type is set to -1, the filename will be assumed to be profiles_<species name>", follow=vars.follow, delay=vars.delay)
        ToolTip(self.timeL, msg="Time point for profile input from IterDB files.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.timeE, msg="Time point for profile input from IterDB files.", follow=vars.follow, delay=vars.delay)
        ToolTip(self.propagateMB, msg="Select which settings will be automatically determined from the input profile files. These values will be set to -1 upon output, and will be reset appropriately when running GENE.", follow=vars.follow, delay=vars.delay)


    def switch_nonlocal(self,weak=0):
        if vars.nonlocal.get()==1:
            self.omnL.config(text='Max. density gradient:')
            self.omtL.config(text='Max. temp. gradient:')
        else:
            self.omnL.config(text='Density gradient:')
            self.omtL.config(text='Temperature gradient:')
        self.update_profgroup()
            
    def delete1(self):
        curr=vars.currspec1V.get()
        for item in [vars.specname,vars.omn,vars.omt,vars.mass,vars.charge,vars.temp,vars.dens,vars.passive,
                     vars.kappa_T,vars.kappa_n,vars.LT_center,vars.LT_width,vars.Ln_center,vars.Ln_width,
                     vars.src_amp,vars.src_width,vars.src_x0,vars.proftype,vars.src_proftype,
                     vars.delta_x_T,vars.delta_x_n]:
            del item[curr-1]
        if len(vars.specname)==1:
            vars.numspecV.set(1)
            vars.clearspec(2)
        elif len(vars.specname)>=2:
            vars.numspecV.set(int(vars.numspecV.get())-1)
        self.ch_numspec()
        vars.ch_spec1(1)

    def delete2(self):
        curr=vars.currspec2V.get()
        for item in [vars.specname,vars.omn,vars.omt,vars.mass,vars.charge,vars.temp,vars.dens,vars.passive,
                     vars.kappa_T,vars.kappa_n,vars.LT_center,vars.LT_width,vars.Ln_center,vars.Ln_width,
                     vars.src_amp,vars.src_width,vars.src_x0,vars.proftype,vars.src_proftype,
                     vars.delta_x_T,vars.delta_x_n]:
            del item[curr-1]
        vars.numspecV.set(int(vars.numspecV.get())-1)
        if len(vars.specname)==1:
            if curr==1:
                vars.currspec1V.set(1)
            elif curr==2:
                vars.clearspec(2)
        if curr>len(vars.specname):
            vars.currspec2V.set(int(vars.currspec2V.get())-1)
        self.ch_numspec()
        vars.ch_spec1()
        vars.ch_spec2()

    def ShLeft1(self):
        curr=vars.currspec1V.get()
        if not curr>len(vars.specname): vars.save_spec1()
        if curr>=2:
            vars.ch_spec1(curr-1)

    def ShRight1(self):
        curr=vars.currspec1V.get()
        if not curr>len(vars.specname): vars.save_spec1()
        if curr<=len(vars.specname)-1:
            vars.ch_spec1(vars.currspec1V.get()+1)

    def ShLeft2(self):
        curr=vars.currspec2V.get()
        if not curr>len(vars.specname): vars.save_spec2()
        if curr>=2:
            vars.ch_spec2(vars.currspec2V.get()-1)

    def ShRight2(self):
        curr=vars.currspec2V.get()
        if not curr>len(vars.specname): vars.save_spec2()
        if curr<=len(vars.specname)-1:
            vars.ch_spec2(vars.currspec2V.get()+1)

                    
    def newspec(self):
        vars.newspec()
        self.ch_numspec()
                    
                    
    def ch_numspec(self):
        # Donothing is a workaround for a Tkinter bug on HPC-FF
        def Donothing(): pass
        n_spec=int(vars.numspecV.get())
        # ------------------------------------------------------------------
        # this part: update the 'list of namelists' in the development frame
        numofnmls=9+n_spec #parallelization,in_out,box,general,geometry,
        #external_contr,nonlocal_x,units,scan namelists+species
        vars.reset_listofnamelists()
        for n in range(n_spec):
            vars.listofnamelists.append('species'+str(n+1))
        if vars.Devfirst==0:
            #delete all menu entries
            for k in range(len(vars.namelistM)):
                vars.namelistM[k].delete(0,len(vars.namelistM)-1)
            #restore menu entries including species namelists
            for i in range(numofnmls):
                for k in range(len(vars.namelistM)):
                    vars.namelistM[k].add_radiobutton(label=vars.listofnamelists[i],
                                                      variable=vars.namelist[k],value=vars.listofnamelists[i],command=Donothing)
               
        # -------------------------------------------------
        # show or don't show the second species entry fields
        spec2items=[self.ShLeft2B,self.ShRight2B,self.currspec2E,self.specname2E,
                    self.mass2E,self.charge2E,self.passive2C,self.delete2B]
        if vars.proftypeV.get()==0:
            spec2items.extend([self.omn2E,self.omt2E,self.temp2E,self.dens2E])
        if not vars.Specfirst:
            if n_spec<=1:
                for item in spec2items:
                    item.grid_remove()
            elif n_spec>=2:
                for item in spec2items:
                    item.grid()
                vars.ch_spec2()
        return "TRUE"
                

    def sel_geom(self):
        pass

    def prof_switch(self):
        pass
            
    def update_profgroup(self,var=None,ind=0,mode=0):
        #fields for analytical profiles
        fields_analyt=[self.tempL,self.omtL,self.densL,self.omnL,self.temp1E,self.omt1E,
                       self.dens1E,self.omn1E]
        try:
            if int(vars.numspecV.get())>1:
                fields_analyt.extend([self.temp2E,self.omt2E,self.dens2E,self.omn2E])
        except:
            pass
        fields_file=[self.fileL,self.fileE,self.fileB,self.timeL,self.timeE,self.propagateMB]
        if vars.proftypeV.get()==0:
            for item in fields_analyt:
                item.grid()
            for item in fields_file:
                item.grid_remove()
        elif vars.proftypeV.get() in [-1,-2,-3]:
            for item in fields_analyt:
                item.grid_remove()
            for item in fields_file:
                item.grid()
                item.config(state=NORMAL)
            if vars.proftypeV.get()==-1:
                self.timeE.config(state=DISABLED)
                vars.iterdbtimeV.set('')
                self.fileE.config(state=DISABLED)
                vars.fileV.set('profiles_<species>')
                self.fileB.config(state=DISABLED)
            else:
                self.timeE.config(state=NORMAL)
                self.fileB.config(state=NORMAL)
                self.fileE.config(state=DISABLED)
                if vars.fileV.get()=='profiles_<species>':
                    vars.fileV.set('')
                if vars.iterdbfileV.get():
                    vars.fileV.set(vars.iterdbfileV.get())
        if not vars.nonlocal.get():
            for index in range(5,10):
                self.proftypeM.entryconfigure(index,state=DISABLED)
            self.Visualize_Profiles.grid_remove()
        else:
            for index in range(5,10):
                self.proftypeM.entryconfigure(index,state=NORMAL)
            if vars.proftypeV.get() in range(1,6):
                self.Visualize_Profiles.grid()
                for item in fields_analyt:
                    item.grid_remove()
            else:
                self.Visualize_Profiles.grid_remove()
        #try to match column sizes
        val=(float(max(self.omnL.winfo_reqwidth(),self.omtL.winfo_reqwidth()))/self.proftypeL.winfo_reqwidth()
             *self.proftypeL.cget('width'))
        self.proftypeL.config(width=int((val*10+5)/10))
        
    def select_input_file(self):
        file=askopenfilename(title='Please specify profile input file.')
        if vars.proftypeV.get() in [-3,-2]:
            vars.iterdbfileV.set(file)
        vars.fileV.set(file)

    def check_proftype(self,var,ind,mode):
        if not vars.nonlocal.get() and int(vars.proftypeV.get()) in range(1,6):
            vars.proftypeV.set(0)
            launcher.Msg('Profile types 1 to 5 are only available for global runs, setting to 0.')
            
