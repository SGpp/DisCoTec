from Tkinter import *
from Plot import SimplePlot
import math
from tkFileDialog   import askopenfilename,askdirectory,asksaveasfilename

#This routine brings up a second Launcher window which allows the user to enter and view profile data
class profile_visual:

    def __init__(self,spec_in,parent,var):
        global vars, profwin
        vars=var
        #write profile data from Species entry fields to arrays
        vars.save_spec1()
        if int(vars.numspecV.get())>=2: vars.save_spec2()
        profwin=Toplevel(bg=vars.stdbg)
        profwin.title('Profile Setup')
        self.xres=500
        #flag for profile/gradient display
        self.pgvar=0
        self.first=1
        self.plot=1
        self.title='Profile'
        prof=[]
        omprof=[]
        self.spec=IntVar()
        self.spec.set(spec_in)
        self.spec_old=IntVar()
        self.spec_old.set(spec_in)
        self.currspec=StringVar()
        self.currspec.set(vars.specname[spec_in])
        self.files=[]
        for i in range(int(vars.numspecV.get())):
            self.files.append(StringVar())
        self.prof=IntVar()
        self.prof.set(0)
        self.prof_old=IntVar()
        self.prof_old.set(0)
        self.currprof=StringVar()
        self.currprof.set('Temperature')
        self.field1L=StringVar()
        self.field2L=StringVar()
        self.field3L=StringVar()
        self.field4L=StringVar()
        self.field5L=StringVar()
        self.field1V=StringVar()
        self.field2V=StringVar()
        self.field3V=StringVar()
        self.field4V=StringVar()
        self.field5V=StringVar()
        self.field1L.set('Maximum Gradient: ')
        self.field2L.set('Center value: ')
        self.field3L.set('Center position: ')
        self.field4L.set('Width: ')
        self.field5L.set('Delta: ')
        if (vars.proftypeV.get()==0):
            self.field1V.set(vars.omt[spec_in])
        else:
            self.field1V.set(vars.kappa_T[spec_in])            
        self.field2V.set(vars.temp[spec_in])
        self.field3V.set(vars.LT_center[spec_in])
        self.field4V.set(vars.LT_width[spec_in])
        self.field5V.set(vars.delta_x_T[spec_in])
        self.proftype=IntVar()
        self.proftype.set(vars.proftype[spec_in])
        self.proftype.trace("w",self.upd_profoptions)
        
        #try: del g
        #except: pass
        xsize=400
        ysize=300
        self.log=IntVar()
        self.log.set(0)
        x0=55
        y0=25
        g = SimplePlot(profwin,width=xsize+70,height=ysize+80,bg=vars.stdbg)
        g.grid(rowspan=9)#pack( expand=1, fill='both' )
        g.update()

        def save_vals():
            if self.prof_old.get()==0:
                if vars.proftypeV.get()==0:
                    vars.omt[int(self.spec_old.get())]=self.field1V.get()
                else:
                    vars.kappa_T[int(self.spec_old.get())]=self.field1V.get()
                vars.temp[int(self.spec_old.get())]=self.field2V.get()
                vars.LT_center[int(self.spec_old.get())]=self.field3V.get()
                vars.LT_width[int(self.spec_old.get())]=self.field4V.get()
                vars.delta_x_T[int(self.spec_old.get())]=self.field5V.get()
                for i in range(int(vars.numspecV.get())):
                    vars.proftype[i]=int(self.proftype.get())
                vars.proftypeV.set(int(self.proftype.get()))
            elif self.prof_old.get()==1:
                if vars.proftypeV.get()==0:
                    vars.omn[int(self.spec_old.get())]=self.field1V.get()
                else:
                    vars.kappa_n[int(self.spec_old.get())]=self.field1V.get()
                vars.dens[int(self.spec_old.get())]=self.field2V.get()
                vars.Ln_center[int(self.spec_old.get())]=self.field3V.get()
                vars.Ln_width[int(self.spec_old.get())]=self.field4V.get()
                vars.delta_x_n[int(self.spec_old.get())]=self.field5V.get()
                for i in range(int(vars.numspecV.get())):
                    vars.proftype[i]=int(self.proftype.get())
                vars.proftypeV.set(int(self.proftype.get()))
            elif self.prof_old.get()==2:
                vars.src_amp[int(self.spec_old.get())]=self.field1V.get()
                vars.src_x0[int(self.spec_old.get())]=self.field3V.get()
                vars.src_width[int(self.spec_old.get())]=self.field4V.get()
                vars.src_proftype[int(self.spec_old.get())]=self.proftype.get()

        def redraw_profiles():
            if self.prof.get()==0:
                self.currprof.set('Temperature')
                self.field1L.set('Maximum gradient')
                self.field2L.set('Norm. center value')
                self.field3L.set('Center position')
                self.field4L.set('Width: ')
                self.field5L.set('Delta: ')
                if (vars.proftypeV.get()==0):
                    self.field1V.set(vars.omt[int(self.spec.get())])
                else:
                    self.field1V.set(vars.kappa_T[int(self.spec.get())]) 
                self.field2V.set(vars.temp[int(self.spec.get())])
                self.field3V.set(vars.LT_center[int(self.spec.get())])
                self.field4V.set(vars.LT_width[int(self.spec.get())])
                self.field5V.set(vars.delta_x_T[int(self.spec.get())])
                self.proftype.set(vars.proftype[int(self.spec.get())])
                self.typeMB.config(relief='sunken',menu=self.typeM)
            elif self.prof.get()==1:
                self.currprof.set('Density')
                self.field1L.set('Maximum gradient')
                self.field2L.set('Norm. center value')
                self.field3L.set('Center position')
                self.field4L.set('Width: ')
                self.field5L.set('Delta: ')
                if (vars.proftypeV.get()==0):
                    self.field1V.set(vars.omn[int(self.spec.get())])
                else:
                    self.field1V.set(vars.kappa_n[int(self.spec.get())])            
                self.field2V.set(vars.dens[int(self.spec.get())])
                self.field3V.set(vars.Ln_center[int(self.spec.get())])
                self.field4V.set(vars.Ln_width[int(self.spec.get())])
                self.field5V.set(vars.delta_x_n[int(self.spec.get())])
                self.proftype.set(vars.proftype[int(self.spec.get())])
                self.typeMB.config(relief='sunken',menu=self.typeM)
            elif self.prof.get()==2:
                self.currprof.set('Sources')
                self.field1L.set('Amplitude')
                self.field3L.set('Center position')
                self.field4L.set('Width: ')
                if vars.proftype[int(self.spec.get())]!=-1:
                    self.field5.config(state=DISABLED)
                    self.field5E.config(state=DISABLED)
                self.field1V.set(vars.src_amp[int(self.spec.get())])
                self.field3V.set(vars.src_x0[int(self.spec.get())])
                self.field4V.set(vars.src_width[int(self.spec.get())])
                self.proftype.set(vars.src_proftype[int(self.spec.get())])
                self.typeMB.config(relief='sunken',menu=self.typeMS)
            self.currspec.set(vars.specname[int(self.spec.get())])
            def choose_profile(n):
                if vars.proftype[int(self.spec.get())] not in (-1,-2,-3):
                    try:
                        r=float(vars.minorrV.get())
                    except:
                        r=1.0
                    k=float(self.field1V.get())
                    f=float(self.field2V.get())
                    c=float(self.field3V.get())
                    w=float(self.field4V.get())
                    if self.prof.get() in [0,1]:
                        try: d=float(self.field5V.get())
                        except: d=0.
                    elif self.prof.get()==2:
                        d=0
                else:
                    r=0
                    f=0
                    k=0
                    c=0
                    w=0
                    d=0
                return r,f,k,c,w,d

            def compute_profile(n,r,f,k,c,w,d):
                prof=[]
                self.plot=1
                if self.prof.get()!=2:
                    if int(self.proftype.get()) in range(0,6):
                        #set up x-axis
                        self.vector_x=[]
                        [self.vector_x.append(1./self.xres*i) for i in range(self.xres+1)]
                    if int(self.proftype.get())==1:
                        [prof.append(f*math.exp(k*r/(math.sinh(c/w)**2)*(i-c-w*math.cosh(c/w)**2*math.tanh((i-c)/w)))) for i in self.vector_x]
                    elif int(self.proftype.get())==2:
                        [prof.append(f*math.exp(-k*r*w*math.tanh((i-c)/w))) for i in self.vector_x]
                    elif int(self.proftype.get())==3:
                        [prof.append((f*math.cosh((i-c+d)/w)/math.cosh((i-c-d)/w))**(-0.5*k*r*w)) for i in self.vector_x]
                    elif int(self.proftype.get())==4:
                        [prof.append(f*math.exp(-k*r*(i-c-w*(math.tanh((i-c-0.5*d)/w)+math.tanh((i-c+0.5*d)/w))))) for i in self.vector_x]
                    elif int(self.proftype.get())==5:
                        [prof.append(f*math.exp(-k*r*(i-c))) for i in self.vector_x]
                    elif int(self.proftype.get())==0:
                        [prof.append(1.) for i in self.vector_x]
                    elif int(self.proftype.get()) in (-1,-2,-3):
                        prof,omprof=self.get_set_proffile()
                else:
                    if w!=0:
                        if int(self.proftype.get())==1:
                            [prof.append(math.exp(-4.0*math.log(2.)*((i-c)/w)**2)) for i in self.vector_x]
                        elif int(self.proftype.get())==2:
                            [prof.append(-0.5*(math.tanh((i-c-3*w)/w)+math.tanh((c-3*w-i)/w))) for i in self.vector_x]
                        else: [prof.append(0) for i in self.vector_x]
                    else: [prof.append(0) for i in self.vector_x]
                return prof

            def compute_gradient(n,r,k,c,w,d):
                omprof=[]
                self.plot=1
                if self.prof.get()!=2:
                    if int(self.proftype.get()) in range(0,6):
                        #set up x-axis
                        self.vector_x=[]
                        [self.vector_x.append(1./self.xres*i) for i in range(self.xres+1)]
                    if int(self.proftype.get())==1:
                        [omprof.append(k*(math.cosh((i-c)/w)**(-2)-math.cosh(c/w)**(-2))/(1.-math.cosh(c/w)**(-2))) for i in self.vector_x]
                    elif int(self.proftype.get())==2:
                        [omprof.append(k*math.cosh((i-c)/w)**(-2)) for i in self.vector_x]
                    elif int(self.proftype.get())==3:
                        [omprof.append(0.5*k*(math.tanh((i-c+d)/w)-math.tanh((i-c-d)/w))) for i in self.vector_x] 
                    elif int(self.proftype.get())==4:
                        [omprof.append(k*(1.-math.cosh((i-c-0.5*d)/w)**(-2)-math.cosh((i-c+0.5*d)/w)**(-2))) for i in self.vector_x] 
                    elif int(self.proftype.get())==5:
                        [omprof.append(k) for i in self.vector_x] 
                    elif int(self.proftype.get())==0:
                        [omprof.append(k) for i in self.vector_x]
                    elif int(self.proftype.get()) in [-1,-2,-3]:
                        prof,omprof=self.get_set_proffile()
                else: [omprof.append(0) for i in self.vector_x]
                return omprof
                    #Pmw.Blt.Vector()


            r,f,k,c,w,d=choose_profile(int(self.spec.get()))
            title_old=self.title
            if self.pgvar==0:
                prof=compute_profile(int(self.spec.get()),r,f,k,c,w,d)
                self.title='Profile'
            else:
                prof=compute_gradient(int(self.spec.get()),r,k,c,w,d)
                self.title='Logarithmic gradient'

            g.delete(ALL)
            g.create_rectangle(x0, y0, xsize+x0, ysize+y0)
            if self.log.get():
                if not True in [item<=0 for item in prof]:
                    for i in range(len(prof)):
                        prof[i]=math.log10(prof[i])
                else:
                    self.log.set(0)
            if vars.proftype[int(self.spec.get())] not in (-2,-3):
                maxval,minval=g.create_axis(x0,y0,xsize,ysize,self.vector_x,prof,self.prof,self.pgvar,self.log.get())
                g.plot(self.vector_x,prof,x0,y0,xsize,ysize,maxval,minval)
                g.update()
            else:
                g.excuse(x0,xsize,y0,ysize)                
            self.spec_old.set(self.spec.get())
            self.prof_old.set(self.prof.get())

        def redraw():
            save_vals()
            if self.plot:
                redraw_profiles()
                #update all entries in Species frame
                #[w.update_idletasks() for w in self.Specframe.winfo_children()]
            #parent.upd_profinfo()


        def profgrad():
            if self.prof.get()!=2:
                if self.pgvar==0:
                    self.pgvar=1
                    self.Profbutt.config(text='View Profile')
                else:
                    self.pgvar=0
                    self.Profbutt.config(text='View Gradient')
            redraw()

        def release():
            self.files=[]
            self.tryjobpath=0
            for i in range(int(vars.numspecV.get())):
                self.files.append(StringVar())


        #name1=StringVar(); name1.set(vars.specname[0])
        #Radiobutton(profwin,command=redraw,value=0,variable=self.spec,textvariable=name1).grid(sticky='w',row=0,column=1)
        self.specMB=Menubutton(profwin,width=15,textvariable=self.currspec,bg="#eee")
        self.specM=Menu(self.specMB,bg=vars.stdbg)
        self.specMB.config(relief='sunken',menu=self.specM)
        for i in range(len(vars.specname)):
            self.specM.add_radiobutton(label=vars.specname[i],command=redraw,variable=self.spec,value=i)
        self.specMB.grid(sticky='ew',padx=15,pady=15,row=0,column=1)
        #        if int(vars.numspecV.get())>=2:
        #            name2=StringVar(); name2.set(vars.specname[1])
        #            Radiobutton(profwin,command=redraw,value=1,variable=self.spec,textvariable=name2).grid(sticky='w',row=1,column=1)
        self.profMB=Menubutton(profwin,width=15,textvariable=self.currprof,bg="#eee")
        self.profM=Menu(self.profMB,bg=vars.stdbg)
        self.profMB.config(relief='sunken',menu=self.profM)
        self.profM.add_radiobutton(label='Temperature',command=redraw,variable=self.prof,value=0)
        self.profM.add_radiobutton(label='Density',command=redraw,variable=self.prof,value=1)
        self.profM.add_radiobutton(label='Sources',command=redraw,variable=self.prof,value=2)
        self.profMB.grid(sticky='ew',padx=15,pady=15,row=0,column=2)
        self.Profbutt=Button(profwin,bg=vars.stdbg,command=profgrad,width=15,text='View Gradient')
        self.Profbutt.grid(row=1,column=3,sticky='ew',padx=15)
        self.field1=Label(profwin,bg=vars.stdbg,textvariable=self.field1L)
        self.field1.grid(row=1,column=1,sticky='w')
        self.field1E=Entry(profwin,text=self.field1V,bg="white")
        self.field1E.grid(row=1,column=2)
        self.field2=Label(profwin,bg=vars.stdbg,textvariable=self.field2L)
        self.field2.grid(row=2,column=1,sticky='w')
        self.field2E=Entry(profwin,text=self.field2V,bg="white")
        self.field2E.grid(row=2,column=2)
        self.field3=Label(profwin,bg=vars.stdbg,textvariable=self.field3L)
        self.field3.grid(row=3,column=1,sticky='w')
        self.field3E=Entry(profwin,text=self.field3V,bg="white")
        self.field3E.grid(row=3,column=2)
        self.field4=Label(profwin,bg=vars.stdbg,textvariable=self.field4L)
        self.field4.grid(row=4,column=1,sticky='w')
        self.field4E=Entry(profwin,text=self.field4V,bg="white")
        self.field4E.grid(row=4,column=2)
        self.field5=Label(profwin,bg=vars.stdbg,textvariable=self.field5L)
        self.field5.grid(row=5,column=1,sticky='w')
        self.field5E=Entry(profwin,text=self.field5V,bg="white")
        self.field5E.grid(row=5,column=2)
        self.typeL=Label(profwin,bg=vars.stdbg,text='Profile type: ')
        self.typeL.grid(row=8,column=1,sticky='w')
        pdict={-3: 'D-IIID ITER-DB', -2: 'ITER-DB', -1: 'GENE Profile file',
                0: 'Flat', 1: 'Peaked (flat boundaries)', 2: 'Peaked', 
                3: 'Flat-top',4: 'Falchetto', 5: 'Exponential'}
        spdict={0: 'Off', 1: 'Peaked source profile', 2: 'Flat-top profile (GYSELA)'}
        self.typeMB=Menubutton(profwin,textvariable=self.proftype,bg="#eee")
        self.typeM=Menu(self.typeMB,bg=vars.stdbg)
        self.typeMS=Menu(self.typeMB,bg=vars.stdbg)
        self.typeMB.config(relief='sunken',menu=self.typeM)
        for i in range(-3,6):
            self.typeM.add_radiobutton(label='%d - %s' %(i,pdict[i]),command=redraw,variable=self.proftype,value=i)
        for i in range(0,3):
            self.typeMS.add_radiobutton(label='%d - %s' %(i,spdict[i]),command=redraw,variable=self.proftype,value=i)
        self.typeMB.grid(row=8,column=2,padx=15)
        Button(profwin,bg=vars.stdbg,command=redraw,text='Redraw!').grid(row=2,padx=15,column=3,sticky='ew')
        Button(profwin,bg=vars.stdbg,command=release,text='Clear Input files').grid(sticky='ew',padx=15,row=3,column=3)
        def destroy():
            save_vals()
#            parent.upd_profinfo()
            profwin.destroy()
        Button(profwin,bg=vars.stdbg,command=destroy,text='Close').grid(sticky='ew',padx=15,row=4,column=3)
        Checkbutton(profwin,bg=vars.stdbg,command=redraw,variable=self.log,text='Logarithmic').grid(sticky='ew',padx=15,row=5,column=3)
        redraw()
        profwin.mainloop()

    def upd_profoptions(self,var,ind,mode):
        po_dict={-3: (0,0,0,0,0,1,1), -2: (0,0,0,0,0,1,1), -1: (0,0,0,0,0,1,0),
                  0: (1,1,0,0,0,0,0), 1: (1,1,1,1,0,0,0), 2: (1,1,1,1,0,0,0),
                  3: (1,1,1,1,1,0,0), 4: (1,1,1,1,1,0,0), 5: (1,1,0,1,0,0,0)}
        entryfields=(self.field1E, self.field2E, self.field3E, self.field4E,
                     self.field5E)
        for i in range(len(entryfields)):
            if po_dict[self.proftype.get()][i]==0:
                entryfields[i].config(state=DISABLED)
            else:
                entryfields[i].config(state=NORMAL)
        if self.proftype.get() in (-1,-2,-3):
            vars.proffile[self.spec.get()]=''
            self.files[self.spec.get()].set('')
               
##     def get_set_proffile(self):
##         prof=[]
##         omprof=[]
##         if self.files[self.spec.get()].get():
##             file=self.files[self.spec.get()].get()
##             f=open(file.strip("'"'"'),'r')
##         elif vars.proffile[self.spec.get()]:
##             file=vars.proffile[self.spec.get()]
##             f=open(file.strip("'"'"'),'r')
##         else:
##             file=askopenfilename(title='Please specify input file for '+vars.specname[int(self.spec.get())]+' species.',parent=profwin)
##             try:
##                 f=open(file,'r')
##                 self.field6V.set("'"+file+"'")
##             except:
##                 self.plot=0
##                 self.field6V.set("''")
##                 return
##         if int(self.proftype.get())==-1:
##             self.files[self.spec.get()].set(file)
##             lines=f.readlines()
##             self.vector_x=[float(line.split()[0]) for line in lines[2:]]
##             self.xres=len(self.vector_x)
##             [prof.append(float(line.split()[self.prof.get()+2])) for line in lines[2:]]
##             omprof.append(-(prof[1]-prof[0])/(self.vector_x[1]-self.vector_x[0])/prof[0])
##             for i in range(1,len(prof)-1):
##                 omprof.append(-(prof[i+1]-prof[i-1])/(self.vector_x[i+1]-self.vector_x[i-1])/prof[i])
##                 omprof.append(-(prof[-1]-prof[-2])/(self.vector_x[-1]-self.vector_x[-2])/prof[-1])
##             for spec in range(0,int(vars.numspecV.get())):
##                 vars.proftype[spec]=self.proftype.get()
##         elif int(self.proftype.get()) in [-2,-3]:
##             for spec in range(0,int(vars.numspecV.get())):
##                 self.files[spec].set("'"+file.strip("'")+"'")
##                 vars.proffile[spec]="'"+file.strip("'")+"'"
##                 vars.proftype[spec]=self.proftype.get()
##         f.close()    
##         return prof,omprof
