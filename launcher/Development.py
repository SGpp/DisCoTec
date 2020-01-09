from Tkinter import *
from Tooltips import *
from ScrFrame import ScrolledFrame
from Pmw import Group

class Development:
    def __init__(self,parent,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        #vars=Variables()
        
        self.Devframe=Frame(parent,width=640,height=480,bg=vars.stdbg)
        self.define_widgets()
        #self.set_tooltips()
        self.display()
        vars.Devfirst=0

    def define_widgets(self):
        self.DevG = Group(self.Devframe, tag_text='Additional parameters',hull_bg=vars.stdbg,ring_bg=vars.stdbg,
                          tag_bg=vars.stdbg,tag_font=("Helvetica",12,'bold'))
        self.DevG.interior().config(bg=vars.stdbg,padx=15,pady=5)
        self.addparframe=ScrolledFrame(self.DevG.interior(), bg=vars.stdbg,width=500,height=350)

    def display(self):
        self.DevG.grid(padx=15,ipadx=15,ipady=15,sticky='news')
        self.addparframe.grid(row=1)
        for i in range(50):
            self.AddPar()
        #self.addparframe.grid_propagate(0)


    def AddPar(self):
        def Donothing(): pass
        sv=StringVar()
        sv2=StringVar()
        sv3=StringVar()
        vars.addedpars.append(sv)
        vars.addedvalues.append(sv2)
        vars.namelist.append(sv3)
        i=len(vars.addedpars)-1
        rownum=i+2
        vars.EntryDevPar.append(Entry(self.addparframe,bg="white",width=14,textvariable=vars.addedpars[i]))
        vars.EntryDevPar[i].grid(row=rownum)
        ToolTip(vars.EntryDevPar[i], msg='Enter parameter name', follow=vars.follow, delay=0.5)
        Label(self.addparframe,bg=vars.stdbg,width=1,text='=').grid(row=rownum,column=1,padx=3)
        vars.EntryDevVal.append(Entry(self.addparframe,bg="white",width=20,textvariable=vars.addedvalues[i]))
        vars.EntryDevVal[i].grid(row=rownum,column=2)
        ToolTip(vars.EntryDevVal[i], msg='Enter value', follow=vars.follow, delay=0.5)
        Label(self.addparframe,bg=vars.stdbg,width=1,text='-->').grid(row=rownum,column=3,padx=3)
        vars.namelistMB.append(Menubutton(self.addparframe,textvariable=vars.namelist[i],width=15,bg=vars.stdbg))
        vars.namelistMB[i].config(textvariable=vars.namelist[i],bg=vars.stdbg)      
        vars.namelistM.append(Menu(vars.namelistMB[i],bg=vars.stdbg))
        vars.namelistMB[i].config(relief='sunken',menu=vars.namelistM[i],bg="#eee")
        for j in range(len(vars.listofnamelists)):
            vars.namelistM[i].add_radiobutton(label=vars.listofnamelists[j],variable=vars.namelist[i],value=vars.listofnamelists[j],command=Donothing)
            ToolTip(vars.namelistMB[i], msg='Namelist to which the parameter should be added', follow=vars.follow, delay=0.5)
            vars.namelistMB[i].grid(row=rownum,column=4)
