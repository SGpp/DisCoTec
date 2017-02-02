from Tkinter import *
from Tooltips import *

class Nonlocal:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Nonlocalframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()


    def define_widgets(self):
        self.radbcL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Radial boundary condition:")
        self.radbcMB=Menubutton(self.Nonlocalframe,textvariable=vars.radbctypeV,bg="#eee")
        radbcM=Menu(self.radbcMB,bg=vars.stdbg)
        self.radbcMB.config(relief='sunken',menu=radbcM, width=17)
        radbcM.add_radiobutton(label="0 -- periodic",variable=vars.radbctypeV,value="0")
        radbcM.add_radiobutton(label="1 -- Dirichlet",variable=vars.radbctypeV,value="1")
        radbcM.add_radiobutton(label="2 -- von-Neumann",variable=vars.radbctypeV,value="2")
        self.LeftbufL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Left buffer zone:")
        self.lcoefL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Coefficient:")
        self.lcoefE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.lcoefkrookV)
        self.lwidthL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Width:")
        self.lwidthE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.lbufsizeV)
        self.lpowL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Power:")
        self.lpowE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.lpowkrookV)
        self.RightbufL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Right buffer zone:")
        self.ucoefL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Coefficient:")
        self.ucoefE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.ucoefkrookV)
        self.uwidthL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Width:")
        self.uwidthE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.ubufsizeV)
        self.upowL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Power:")
        self.upowE=Entry(self.Nonlocalframe,bg="white",width=5,text=vars.upowkrookV)
        self.ckheatL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Krook-type heating:")
        self.ckheatE=Entry(self.Nonlocalframe,bg="white",width=6,text=vars.ckheatV)
        self.resetlimitL=Label(self.Nonlocalframe,bg=vars.stdbg,text="Reset limit:")
        self.resetlimitE=Entry(self.Nonlocalframe,bg="white",width=6,text=vars.resetlimitV)
        self.arakawaC=Checkbutton(self.Nonlocalframe,highlightbackground=vars.stdbg,bg=vars.stdbg,text="Arakawa nonlinearity",
                             variable=vars.arakawaV)
        self.shiftedC=Checkbutton(self.Nonlocalframe,highlightbackground=vars.stdbg,bg=vars.stdbg,text="Shifted metric",
                             variable=vars.shiftedV)





    def display(self):
        Label(self.Nonlocalframe,bg=vars.stdbg,text="Parameters for nonlocal simulations:",
              font=("Helvetica",12,'bold')).grid(sticky='w',columnspan=3)
        self.radbcL.grid(sticky='w')
        self.radbcMB.grid(sticky='ew',row=1,column=1,columnspan=2)
        self.LeftbufL.grid(sticky='w')
        self.lcoefL.grid(sticky='w',row=2,column=1)
        self.lcoefE.grid(sticky='ew',row=2,column=2)
        self.lwidthL.grid(sticky='w',row=3,column=1)
        self.lwidthE.grid(sticky='ew',row=3,column=2)
        self.lpowL.grid(sticky='w',row=4,column=1)
        self.lpowE.grid(sticky='ew',row=4,column=2)
        self.RightbufL.grid(sticky='w')
        self.ucoefL.grid(sticky='w',row=5,column=1)
        self.ucoefE.grid(sticky='ew',row=5,column=2)
        self.uwidthL.grid(sticky='w',row=6,column=1)
        self.uwidthE.grid(sticky='ew',row=6,column=2)
        self.upowL.grid(sticky='w',row=7,column=1)
        self.upowE.grid(sticky='ew',row=7,column=2)
        self.ckheatL.grid(sticky='w')
        self.ckheatE.grid(sticky='w',row=8,column=1)
        self.resetlimitL.grid(sticky='w')
        self.resetlimitE.grid(sticky='w',row=9,column=1)
        Label(self.Nonlocalframe,bg=vars.stdbg).grid()
        self.arakawaC.grid(sticky='w')
        self.shiftedC.grid(sticky='ew', row=11,column=1,columnspan=2)


    def set_tooltips(self):
        ToolTip(self.radbcMB, msg="Radial boundary condition: \n0=periodic, 1=Dirichlet, 2=von-Neumann",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.radbcL, msg="Radial boundary condition: \n0=periodic, 1=Dirichlet, 2=von-Neumann",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.LeftbufL, msg='Adjustments for the buffer zone at the inner radial boundary',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.RightbufL, msg='Adjustments for the buffer zone at the outer radial boundary',
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.ckheatL, msg='Coefficient for Krook-type heat source -- using this option allows to perform global simulations with approximately fixed profiles. Set to values of about 1/10 of the linear growth rate.', follow=vars.follow, delay=vars.delay)
        ToolTip(self.resetlimitL, msg='If delta-f exceeds a certain relative deviation from F_0, the equilibrium is initialized to stay within the gyrokinetic delta-f ordering. Value is given in units of F_0, default=1000 (never reinitialize)',follow=vars.follow, delay=vars.delay)
        ToolTip(self.arakawaC, msg='Use the Arakawa representation to compute the nonlinearity (recommended).'
                ,follow=vars.follow, delay=vars.delay)
        ToolTip(self.shiftedC, msg='Use shifted metric geometry, resulting in an orthogonalized perpendicular grid.'
                ,follow=vars.follow, delay=vars.delay)
        
