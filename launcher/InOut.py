from Tkinter import *
from Tooltips import *
from Pmw import Group

class Inout:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Inoutframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.IOgroup = Group(self.Inoutframe, tag_text='Input/Output parameters',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                               tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.IOgroup.interior().config(bg=vars.stdbg)

        self.diagdirL=Label(self.IOgroup.interior(),bg=vars.stdbg,text="Output directory:")
        self.chptdirL=Label(self.IOgroup.interior(),bg=vars.stdbg,text="Checkpoint directory:")
        self.diagdirE=Entry(self.IOgroup.interior(),bg="white",width=30,text=vars.diagdirV)
        self.chptdirE=Entry(self.IOgroup.interior(),bg="white",width=30,text=vars.chptdirV)
        self.readchptC=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Read checkpoint",
                                   variable=vars.rdchptV)
        self.writechptC=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Write checkpoint",variable=vars.wrtchptV)
        self.readchptC=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Read checkpoint",
                                   variable=vars.rdchptV)
        self.writeh5C=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Write HDF5 output",variable=vars.wrth5V)
        self.writestdC=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Write Standard output",variable=vars.wrtstdV)
        self.chpth5C=Checkbutton(self.IOgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="HDF5 checkpoint",variable=vars.chpth5V)
        
        self.istepgroup = Group(self.Inoutframe, tag_text='Output intervals',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                               tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.istepgroup.interior().config(bg=vars.stdbg)
        self.istepnrgL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Avg. moments (nrg):")
        self.istepomegaL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Omega diag.:")
        self.istepfldL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Fields:")
        self.istepmomL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Moments:")
        self.istepenergyL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Energy:")
        self.istepenergy3dL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Energy (3D):")
        self.istepprofL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Profiles:")
        self.istepvspL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="V-space:")
        self.istepncL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="NC fluxes:")
        self.istepschptL=Label(self.istepgroup.interior(),bg=vars.stdbg,text="Secure checkpoint:")
        
        self.istepnrgE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepnrgV)
        self.istepfldE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepfldV)
        self.istepomegaE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepomegaV)
        self.istepmomE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepmomV)
        self.istepenergyE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepenergyV)
        self.istepenergy3dE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepenergy3dV)
        self.istepprofE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepprofV)
        self.istepvspE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepvspV)
        self.istepncE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepncV)
        self.istepschptE=Entry(self.istepgroup.interior(),bg="white",width=6,text=vars.istepschptV)



    def display(self):
        self.IOgroup.grid(row=0,sticky='ew',padx=15,pady=15,ipadx=15,columnspan=3)
        self.IOgroup.interior().config(padx=15,pady=5)
        self.diagdirL.grid(sticky='w')
        self.chptdirL.grid(sticky='w')
        self.diagdirE.grid(row=0,column=1,columnspan=3,sticky='w')
        self.chptdirE.grid(row=1,column=1,columnspan=3,sticky='w')
        self.readchptC.grid(sticky='w')
        self.writechptC.grid(row=2,column=1,sticky='w')
        self.writestdC.grid(sticky='w')
        self.writeh5C.grid(row=3,column=1,sticky='w')
        self.chpth5C.grid(sticky='w')
        self.istepgroup.grid(row=1,sticky='ew',padx=15,pady=15,ipadx=15,columnspan=3)
        self.istepgroup.interior().config(padx=15,pady=5)
        self.istepnrgL.grid(sticky='w')
        self.istepfldL.grid(sticky='w')
        self.istepmomL.grid(sticky='w')
        self.istepenergyL.grid(sticky='w')
        self.istepenergy3dL.grid(sticky='w')
        self.istepprofL.grid(sticky='w')
        self.istepvspL.grid(sticky='w')
        self.istepncL.grid(sticky='w')
        self.istepschptL.grid(sticky='w')
        self.istepomegaL.grid(sticky='w')
        self.istepnrgE.grid(row=0,column=1,sticky='w')
        self.istepfldE.grid(row=1,column=1,sticky='w')
        self.istepmomE.grid(row=2,column=1,sticky='w')
        self.istepenergyE.grid(row=3,column=1,sticky='w')
        self.istepenergy3dE.grid(row=4,column=1,sticky='w')
        self.istepprofE.grid(row=5,column=1,sticky='w')
        self.istepvspE.grid(row=6,column=1,sticky='w')
        self.istepncE.grid(row=7,column=1,sticky='w')
        self.istepschptE.grid(row=8,column=1,sticky='w')
        self.istepomegaE.grid(row=9,column=1,sticky='w')







    def set_tooltips(self):
        ToolTip(self.diagdirL, msg='diagdir', follow=vars.follow, delay=vars.delay)
        ToolTip(self.chptdirL, msg='chptdir: If left empty, checkpoints will be read/written to the output directory.', 
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.readchptC, msg='read_checkpoint: Restart the simulation from a checkpoint file.', follow=vars.follow, delay=vars.delay)
        ToolTip(self.writechptC, msg='write_checkpoint: Write a checkpoint file when finishing the simulation.', follow=vars.follow, delay=vars.delay)
        ToolTip(self.chpth5C, msg='chpt_h5: Read/Write checkpoint in HDF5 format (filename "checkpoint.h5")', follow=vars.follow, delay=vars.delay)
        ToolTip(self.writeh5C, msg='write_h5: Write GENE output into HDF5 files', follow=vars.follow, delay=vars.delay)
        ToolTip(self.writestdC, msg='write_std: Write GENE output to standard binary/ASCII files', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepnrgL, msg='istep_nrg: Interval for nrg.dat output', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepomegaL, msg='istep_omega: Interval for frequency/growth rate evaluation', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepfldL, msg='istep_field: Interval for field.dat output', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepmomL, msg='istep_mom: Interval for mom_<species>.dat output', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepenergyL, msg='istep_energy: Interval for energy diagnostics', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepenergy3dL, msg='istep_energy3d: Interval for energy diagnostics (3D data)', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepprofL, msg='istep_prof: Interval for profile output in nonlocal runs', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepvspL, msg='istep_vsp: Interval for vsp.dat output', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepncL, msg='istep_neoclass: Interval for neoclass.dat output', follow=vars.follow, delay=vars.delay)
        ToolTip(self.istepschptL, msg='istep_schpt: Interval secure checkpoint writing (allows restart of aborted runs)',
                follow=vars.follow, delay=vars.delay)
        





        
    def switch_nonlocal(self,weak=0):
        if vars.nonlocal.get()==1:
            self.Nonlocalbutt.config(state='normal')
            launcher.istepprofL.grid()
            launcher.istepprofE.grid()
        else:
            self.Nonlocalbutt.config(state='disabled')
            launcher.istepprofL.grid_remove()
            launcher.istepprofE.grid_remove()

