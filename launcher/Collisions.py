from Tkinter import *
from Tooltips import *
from Pmw import Group

class Collisions(object):
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Collframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        self.set_tooltips()
        self.display()

    def define_widgets(self):
        self.collgroup = Group(self.Collframe,tag_text='Collisions',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                              tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.collgroup.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.collisnoRB=Radiobutton(self.collgroup.interior(),highlightbackground=vars.stdbg,text="No collisions",bg=vars.stdbg,
                                    variable=vars.collisV, value="'none'",width=20,anchor="w",command=self.switch_collop)
        self.collislanRB=Radiobutton(self.collgroup.interior(),highlightbackground=vars.stdbg,text="Landau-Boltzmann",bg=vars.stdbg,
                                     variable=vars.collisV, value="'landau'",width=20,anchor="w",command=self.switch_collop)
        self.collispaRB=Radiobutton(self.collgroup.interior(),highlightbackground=vars.stdbg,text="Pitch-Angle",bg=vars.stdbg,
                                    variable=vars.collisV, value="'pitch-angle'",width=20,anchor="w",command=self.switch_collop)

        self.spacediffC=Checkbutton(self.collgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Spatial diffusion",
                                    variable=vars.spacediffV)
        self.collL=Label(self.collgroup.interior(),bg=vars.stdbg,text="Collisionality:")
        self.collE=Entry(self.collgroup.interior(),bg="white",width=18,text=vars.collV)
        self.ZeffL=Label(self.collgroup.interior(),bg=vars.stdbg,text="Effective charge:")
        self.ZeffE=Entry(self.collgroup.interior(),bg="white",width=18,text=vars.ZeffV)
    
        self.consgroup = Group(self.Collframe,tag_text='Conservation model', hull_bg= vars.stdbg,ring_bg=vars.stdbg,
                               tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.consgroup.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.collcmMB=Menubutton(self.consgroup.interior(),width=15,textvariable=vars.collcmV,bg="#eee")
        self.collcmM=Menu(self.collcmMB,bg=vars.stdbg)
        self.collcmMB.config(relief='sunken',menu=self.collcmM)
        self.collcmM.add_radiobutton(label="None",variable=vars.collcmV,value="'none'")
        self.collcmM.add_radiobutton(label="Xu-Rosenbluth",variable=vars.collcmV,value="'xu_rosenbluth'")
        self.collcmM.add_radiobutton(label="Self-adjoint",variable=vars.collcmV,value="'self_adj'")

        self.collffmC=Checkbutton(self.collgroup.interior(),highlightbackground=vars.stdbg,bg=vars.stdbg,text="Use f/fM form",
                                  variable=vars.collffmV)

        self.advanced=set([self.spacediffC,self.collispaRB,self.collffmC,self.collcmMB])
        self.colloptions=set([self.spacediffC,self.collL,self.collE,self.ZeffL,self.ZeffE,self.collffmC,self.collcmMB])

    def display(self):
        self.collgroup.grid(row=0,column=0,padx=15,sticky='news')
        self.collisnoRB.grid(row=0,column=0,sticky='w')
        self.collislanRB.grid(row=1,column=0,sticky='w')
        self.collispaRB.grid(row=2,column=0,sticky='w')
        self.spacediffC.grid(row=5,column=0,sticky='w')
        self.collL.grid(row=3,sticky='w')
        self.collE.grid(row=3,column=1,sticky='ew')
        self.ZeffL.grid(row=4,sticky='w')
        self.ZeffE.grid(row=4,column=1,sticky='ew')
        self.consgroup.grid(row=1,column=0, padx=15, sticky='news')
        self.collcmMB.grid(row=0,column=0,sticky='ew')
        self.collffmC.grid(row=5,column=1,sticky='w')

    def set_tooltips(self):
        ToolTip(self.collL, msg="coll: Normalized collisionality (can be calculated self-consistently in 'Reference' tab)",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.collL, msg="Zeff: Effective ion charge, will be multiplied to the collisionality of particles scattering off ions",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.collislanRB, msg="Use the standard collision operator",
                follow=vars.follow, delay=vars.delay)
        ToolTip(self.collispaRB, msg="Collisions with pitch-angle scattering only for benchmarking",
                follow=vars.follow, delay=vars.delay)

    def switch_advanced(self):
        self.switch_collop()

    def switch_collop(self):
        [item.config(state=DISABLED) for item in self.colloptions.union(self.advanced)]
        if vars.collisV.get() =="'landau'" or vars.collisV.get()=="'pitch-angle'":
            vars.collonoffV.set(1)
        else:
            vars.collonoffV.set(0)
        if vars.collisV.get() == "'none'" and vars.expertV.get():
            [item.config(state=NORMAL) for item in self.advanced.difference(self.colloptions)]
        elif vars.collisV.get() == "'landau'" and not vars.expertV.get():
            [item.config(state=NORMAL) for item in self.colloptions.difference(self.advanced)]
        elif vars.collisV.get() == "'landau'" and vars.expertV.get():
            [item.config(state=NORMAL) for item in self.colloptions.union(self.advanced)]
        elif vars.collisV.get() == "'pitch-angle'" and not vars.expertV.get():
            [item.config(state=NORMAL) for item in self.colloptions.difference(self.advanced)]
            vars.collcmV.set("'none'")
        elif vars.collisV.get() == "'pitch-angle'" and vars.expertV.get():
            [item.config(state=NORMAL) for item in self.colloptions.union(self.advanced)]
            vars.collcmV.set("'none'")
    
