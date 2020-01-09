#!/usr/bin/env python
# File: GENE-GUI.py

import os
import re
import math
import sys
try: import Tkinter
except: sys.exit('Need module Tkinter, please install a Python distribution which includes this module.')


#import Tix
from Tkinter import *
from time import sleep
import tkFont
#include path for modules
startdir=os.path.join(os.path.split(sys.argv[0])[0])
sys.path=sys.path+[os.path.join(startdir,'launcher')]
#import subclasses from modules
from ScrFrame import ScrolledFrame
from Tooltips import *
from Operation import Operation
from InOut import Inout
from Domain import Domain
from General import General
from Collisions import Collisions
from Geometry import Geometry
from Species import Species
from Development import Development
from Reference import Reference
from Nonlocal import Nonlocal
from vars import *
from Parameters import Parameters
from Multijob import Multijob
from Submit import Submit
import Pmw
from subprocess import Popen,PIPE,call

class Launcher:

	def __init__(self,master):
		global vars,pars,sub		
		#save launch directory
		os.chdir(os.path.split(sys.argv[0])[0])
		
		self.startdir=os.getcwd()
		vars=Variables(self)
		pars=Parameters(self,vars)

		
		self.frame=Frame(master,bg=vars.stdbg)
		self.frame.grid()
		self.createWidgets()
		self.set_tooltips()
		sub=Submit(self,vars)
		self.Expert()
		self.Grid_Frames()
		vars.clear()
		vars.newspec()
		self.Msg("""Hints: 
Many labels and buttons show tooltips when the mouse pointer hovers over them. The names of the actual parameters can also be found there.
The "Advanced options" button will activate some of the greyed-out fields. Others may be inactive because they do not apply to the present settings.
Use the "Default Values" button to fill the main fields with standard values.""",1)

	def createWidgets(self):
		title=Label(self.frame, text="GENE Launcher",bg=vars.stdbg,fg="#000",font=('Helvetica', 18,'bold'),relief="ridge",padx=10)
		title.grid(columnspan=vars.nummulti,column=1,sticky='EW')

		#left hand side buttons to choose from different parameter categories
		self.Opbutt= Button(self.frame, text="Operation",bg=vars.stdbg, width=15)
		self.Opbutt.grid(row=2,column=0)
		self.Opbutt.bind("<Button-1>",self.Open_Frame)
		self.Inoutbutt= Button(self.frame, text="Input/Output", width=15, bg=vars.stdbg)
		self.Inoutbutt.grid(row=3,column=0)
		self.Inoutbutt.bind("<Button-1>",self.Open_Frame)
		self.Dombutt= Button(self.frame, text="Domain", width=15, bg=vars.stdbg)
		self.Dombutt.grid(row=4,column=0)
		self.Dombutt.bind("<Button-1>",self.Open_Frame)
		self.Genbutt= Button(self.frame, text="General", width=15, bg=vars.stdbg)
		self.Genbutt.grid(row=5,column=0)
		self.Genbutt.bind("<Button-1>",self.Open_Frame)

		self.Collbutt= Button(self.frame, text="Collisions", width=15, bg=vars.stdbg)
		self.Collbutt.grid(row=6,column=0)
		self.Collbutt.bind("<Button-1>",self.Open_Frame)

		self.Geobutt= Button(self.frame, text="Geometry", width=15, bg=vars.stdbg)
		self.Geobutt.grid(row=7,column=0)
		self.Geobutt.bind("<Button-1>",self.Open_Frame)

		self.Specbutt= Button(self.frame, text="Species", width=15, bg=vars.stdbg)
		self.Specbutt.grid(row=8,column=0)
		self.Specbutt.bind("<Button-1>",self.Open_Frame)

		self.Nonlocalbutt= Button(self.frame, text="Nonlocal", width=15, bg=vars.stdbg)
		self.Nonlocalbutt.grid(row=9,column=0)
		if vars.allow_nonlocal:
			self.Nonlocalbutt.bind("<Button-1>",self.Open_Frame)

		self.Refbutt= Button(self.frame, text="Ref. Values", width=15, bg=vars.stdbg)
		self.Refbutt.grid(row=10,column=0)
		self.Refbutt.bind("<Button-1>",self.Open_Frame)

		self.Devbutt= Button(self.frame, text="Development", width=15, bg=vars.stdbg)
		self.Devbutt.grid(row=11,column=0)
		self.Devbutt.bind("<Button-1>",self.Open_Frame)

		Label(self.frame,text='',bg=vars.stdbg).grid()
		Label(self.frame,text='',bg=vars.stdbg).grid()
		Label(self.frame,text='',bg=vars.stdbg).grid()
		Label(self.frame,text='',bg=vars.stdbg).grid()
#		Label(self.frame,text='',bg=vars.stdbg).grid()
#		Label(self.frame,text='',bg=vars.stdbg).grid()
		#Small frame for special buttons
		self.smframe=Frame(self.frame,bg=vars.stdbg,relief='groove',bd=2)
		self.smframe.grid(row=13,column=0,rowspan=3)
		
		#Special buttons
		self.Clearbutt=Button(self.smframe,text="Clear form", width=15, bg=vars.stdbg, command=vars.clear_button)
		self.Clearbutt.grid(sticky='EW')
		self.Defbutt=Button(self.smframe,text="Default values", width=15, bg=vars.stdbg, command=self.set_defaults_button)
		self.Defbutt.grid(sticky='EW')
		self.Quickbutt=Button(self.frame,text='Quick Intro',bg=vars.stdbg,command=self.quickintro)
		self.Quickbutt.grid(row=18,ipadx=5)
		self.Helpbutt=Button(self.frame,bitmap='question',bg=vars.stdbg,command=self.help)
		self.Helpbutt.grid(row=19,ipadx=5)
		testing=0
		if testing:
			self.Testbutt=Button(self.smframe,text="Test", width=15, bg=vars.stdbg, command=self.test_output)
			self.Testbutt.grid(sticky='EW')
		
		#Multitasking buttons
		Label(self.frame,bg=vars.stdbg,text='Jobs =>',anchor='e').grid(row=1,column=0,sticky='e')

		#initialization of multi-job and submit modules
		multi=Multijob(self,vars)

		#set up buttons for job selection
		self.multibutt=[]
		for i in range(vars.nummulti):
			self.multibutt.append(Button(self.frame,width=1,bg=vars.stdbg,text=str(i+1)))
			self.multibutt[i].grid(row=1,column=i+1,sticky='ew')
			def multiswitch(event,self=self,i=i):
				return multi.multiswitch(event,i)
			self.multibutt[i].bind(sequence="<Button-1>",func=multiswitch)

		self.multibutt[0].config(bg="#eee")

		self.frame.columnconfigure(1,weight=1)

		#right hand side buttons controlling launcher actions
		self.readB= Button(self.frame, text="Read parameters",width=15,bg=vars.stdbg, command=pars.Read_Pars)
		self.readB.grid(row=2,column=vars.nummulti+1)	
		self.newprobB= Button(self.frame, text="New 'prob' dir.",width=15,bg=vars.stdbg, command=self.Newprob)
		self.newprobB.grid(row=3,column=vars.nummulti+1)	
		self.writeB= Button(self.frame, text="Write parameters",width=15,bg=vars.stdbg, command=self.save)
#		self.writeB.bind("<ButtonRelease-1>",self.save)
		self.writeB.grid(row=4,column=vars.nummulti+1)	
		self.saveasB= Button(self.frame, text="Save as...",width=15,bg=vars.stdbg, command=self.saveas)
#		self.saveasB.bind("<ButtonRelease-1>",self.saveas)
		self.saveasB.grid(row=5,column=vars.nummulti+1)	
		self.checkB= Button(self.frame, text="Check",width=15,bg=vars.stdbg, command=pars.Check_Pars)
		self.checkB.grid(row=6,column=vars.nummulti+1)	
		self.submitB= Button(self.frame, text="Submit",width=15,bg=vars.stdbg, fg="red", command=self.submit_job)
		self.submitB.grid(row=7,column=vars.nummulti+1)	
		self.quitB= Button(self.frame, text="Quit",width=15,bg=vars.stdbg, command=sys.exit)
		self.quitB.grid(row=8,column=vars.nummulti+1)
		
		self.expertC= Checkbutton(self.frame, text="Advanced options",highlightbackground=vars.stdbg,variable=vars.expertV,
					  bg=vars.stdbg, command=self.Expert)
		self.expertC.grid(row=11,column=vars.nummulti+1)


		self.checkscanstatusB= Button(self.frame, text="Check scan status",width=15,bg=vars.stdbg, command=self.check_scan)
		self.checkscanstatusB.grid(row=15,column=vars.nummulti+1)

		
		currjobpathL=LabelFrame(self.frame,font=("Helvetica",12,'bold'),height=100,width=100,bd=2,relief='groove',
					bg=vars.stdbg,text="Current path:")
		ToolTip(currjobpathL, msg='Displays the path the currently selected job has been read from/written to',
			follow=vars.follow, delay=vars.delay)
		currjobpathL.grid(row=23,column=vars.nummulti+1)
		self.jobpathL=Label(currjobpathL,font=("Helvetica",12,'bold'),bg=vars.stdbg,textvariable=vars.jobpath)
		self.jobpathL.grid(row=24,column=vars.nummulti+1)
		
		#frame for displaying messages, e.g. to report errors in the parameters
		self.Msgframe=LabelFrame(self.frame,font=("Helvetica",12,'bold'),text='Messages',height=100,bg=vars.stdbg,relief='groove',bd=2)
		for i in range(vars.nummulti):
			self.frame.columnconfigure(i+1,weight=1)
		self.Msgframe.grid(row=22,column=1,columnspan=vars.nummulti,rowspan=4,sticky='ew')
		self.MsgV=StringVar()
		self.MsgM=Text(self.Msgframe,bg=vars.stdbg,width=20,height=4,wrap='word',state='disabled')#,textvariable=self.MsgV)
		self.MsgM.grid(row=0,sticky='news',padx=5,pady=5)
		self.Msgframe.columnconfigure(0,weight=1)
		self.MsgM.grid_propagate(0)
		#set up scrollbar for message frame
		self.MsgScr=Scrollbar(self.Msgframe,bg=vars.stdbg,orient="vertical",command=self.MsgM.yview)
		self.MsgScr.grid(row=0,column=1,sticky='ns')
		self.MsgM.config(yscrollcommand=self.MsgScr.set)

		#initialize all frames
		self.Inoutframe=Inout(self.frame,self,vars)
		self.Domframe=Domain(self.frame,self,vars)
		self.Genframe=General(self.frame,self,vars)
		self.Collframe=Collisions(self.frame,self,vars)
		self.Geoframe=Geometry(self.frame,self,vars)
		self.Specframe=Species(self.frame,self,vars)
		self.Nonlocalframe=Nonlocal(self.frame,self,vars)
		self.Refframe=Reference(self.frame,self,vars)
		self.Devframe=Development(self.frame,self,vars)
		self.Opframe=Operation(self.frame,self,vars)

	def set_tooltips(self):
       		ToolTip(self.Clearbutt, msg='Clear all values', follow=vars.follow, delay=vars.delay)
       		ToolTip(self.Defbutt, msg='Set all values to defaults, depending on linear/nonlinear mode', follow=vars.follow, delay=vars.delay)
       		ToolTip(self.Helpbutt, msg='Open GENE documentation', follow=vars.follow, delay=vars.delay)
		ToolTip(self.jobpathL, msg='Displays the path the currently selected job has been read from/written to',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.readB, msg='Read an existing GENE parameter file, overwriting the current settings.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.newprobB, msg='Generate a new "prob" directory and set the path of the current job accordingly.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.writeB, msg="""Brings up a dialog to select an output directory, and writes 
a "parameters" file into that directory. If the current job already has an 
output path (see "Current path" below), the file will directly be written 
there.""",
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.saveasB, msg='Save a parameters file to a customizable directory and filename.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.checkB, msg='Perform various consistency checks on the present parameter settings.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.submitB, msg="""Submit the simulation as a batch job on the present machine. 
This will only work if a template launcher.cmd submit script is available 
for this machine!""",
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.quitB, msg='Exit the launcher tool.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.expertC, msg='Activate various advanced settings.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.checkscanstatusB, msg="""If a scan has been launched using the present parameters, 
an intermediate scan.log file will be created and displayed in the 
message window.""",
			follow=vars.follow, delay=vars.delay)
		for i in range(vars.nummulti):
			ToolTip(self.multibutt[i], msg='Switch to job no. %d' %(i+1),
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Opbutt, msg='Set up parallelization and operations modes for GENE.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Inoutbutt, msg='Adjust input/output settings.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Dombutt, msg='Set resolution and size of the computational domain.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Genbutt, msg='Adjust initial/termination conditions, dimensionless parameters, equilibrium flows, eigenvalue solver options and settings for neoclassical simulations.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Collbutt, msg='Set collisionality and options for the collision operator treatment.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Geobutt, msg='Settings for the geometry of the background magnetic field.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Specbutt, msg='Set up particle species to be treated, as well as profile settings.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Nonlocalbutt, msg='Settings for global runs.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Refbutt, msg='Set reference values, and compute derived dimensionless parameters.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Devbutt, msg='Add custom parameters, which are not hard-wired in the launcher, but recognized by GENE. When reading parameter files, unrecognized values will be added here.',
			follow=vars.follow, delay=vars.delay)
		ToolTip(self.Quickbutt, msg='How to produce a simple parameters file in three steps.',
			follow=vars.follow, delay=vars.delay)
			

	def test_output(self):
		print vars.specname,vars.omt,vars.omn,vars.dens,vars.temp,vars.charge,vars.mass
		print vars.LT_center,vars.LT_width,vars.kappa_T,vars.kappa_n,vars.Ln_center
		print vars.Ln_width,vars.src_amp,vars.src_x0,vars.src_width
		
	def check_scan(self):
		if vars.scan.get() and vars.jobpath.get():
                        #scanscript must be called from the directory it is contained in
			os.chdir(vars.jobpath.get())
			#suppress scanscript output
			fnull=open(os.devnull,'w')
			call('./scanscript --mks',stdout=fnull,shell=True)
			fnull.close()
			os.chdir(self.startdir)
			try:
				pipe=Popen("grep diagdir "+vars.jobpath.get()+"/parameters", shell=True, stdout=PIPE).stdout
				diagdir=pipe.readlines()[-1].strip('\n').split('=')[1].strip().strip("'")
				file=open(diagdir+'/scan.log')
				scanlog=file.read()
				self.MsgClear()
				self.Msg(scanlog)
				file.close()
			except:
				self.Msg('Error opening scan.log.')
		else:
			if not vars.scan.get():
				self.Msg('Please check the "Parameter scan" option.')
			if not vars.jobpath.get():
				self.Msg('Job path is not yet defined, save parameters first!')
					 
	def highlight_tab(self,widget):
		tabwidgets=[self.Opbutt,self.Inoutbutt,self.Dombutt,self.Genbutt,self.Collbutt,self.Geobutt,self.Specbutt,self.Refbutt,self.Devbutt,self.Nonlocalbutt]
		for item in tabwidgets:
			if item==widget: 
				item.config(bg='#eee')
			else: 
				item.config(bg=vars.stdbg)
	
	#Open operation frame
	def Open_Frame(self,event=None):
		if event:
			dict={self.Opbutt: self.Opframe.Opframe,
			      self.Inoutbutt: self.Inoutframe.Inoutframe,
			      self.Dombutt: self.Domframe.Domframe,
			      self.Genbutt: self.Genframe.Genframe,
			      self.Collbutt: self.Collframe.Collframe,
			      self.Geobutt: self.Geoframe.Geoframe,
			      self.Specbutt: self.Specframe.Specframe,
			      self.Refbutt: self.Refframe.Refframe,
			      self.Devbutt: self.Devframe.Devframe,
			      self.Nonlocalbutt: self.Nonlocalframe.Nonlocalframe,
			      }
			self.highlight_tab(event.widget)
			for item in dict.keys():
				if item!=event.widget:
					dict[item].lower()
			dict[event.widget].lift

	def Grid_Frames(self):
		list=[self.Inoutframe.Inoutframe,
		      self.Domframe.Domframe,
		      self.Genframe.Genframe,
		      self.Collframe.Collframe,
		      self.Geoframe.Geoframe,
		      self.Specframe.Specframe,
		      self.Refframe.Refframe,
		      self.Devframe.Devframe,
		      self.Nonlocalframe.Nonlocalframe,
		      self.Opframe.Opframe]
		for item in list:
			item.grid(column=1,columnspan=vars.nummulti,row=2,rowspan=20,sticky='news',padx=20,pady=20)


	#the following routines trigger actions in various tabs
	def adaptlx_sel(self):
		self.Domframe.adaptlx_sel()

	def switch_collop(self):
		self.Collframe.switch_collop()
	
	def ch_numspec(self):
		self.Specframe.ch_numspec()

	def sel_geom(self,withgeo=1):
		if withgeo:
			self.Geoframe.sel_geom()
		self.Refframe.sel_geom()
		self.Specframe.sel_geom()
		self.Domframe.sel_geom()

	def calcdtclick(self):
		self.Genframe.calcdtclick()

	def autosw(self):
		self.Opframe.autosw()

        #calling this routine with weak=1 does not force the recommended adaptlx, calcdt and initial condition settings
	def switch_lin_nonlin(self,weak=0):
		self.Opframe.switch_lin_nonlin(weak)
		self.Domframe.switch_lin_nonlin(weak)
		self.Genframe.switch_lin_nonlin(weak)

	#here, weak does not force the setting to harmonic etc.
	def switch_ev_iv(self,weak=0):
		self.Opframe.switch_ev_iv(weak)
		self.Genframe.switch_ev_iv(weak)
			
	def switch_nonlocal(self,weak=0):
		self.Domframe.switch_nonlocal(weak)
		self.Geoframe.switch_nonlocal(weak)
		self.Specframe.switch_nonlocal(weak)
		if vars.nonlocal.get()==1:
			self.Nonlocalbutt.config(state='normal')
			self.Nonlocalbutt.bind("<Button-1>",self.Open_Frame)
		else:
			self.Nonlocalbutt.config(state='disabled')
			self.Nonlocalbutt.unbind("<Button-1>")
		self.switch_lin_nonlin()

	def switch_scan(self):
		self.Opframe.switch_scan()
		

	def set_defaults_button(self):
		if askokcancel(message='This will set all entries to default values, please confirm!'):
			vars.set_defaults()
			self.Specframe.ch_numspec()
		else:
			return

	def Newprob(self):
		os.chdir(self.startdir)
		try:
			pipe=Popen('./newprob', shell=True, stdout=PIPE).stdout
			output=''.join(pipe.readlines())
			self.Msg(output)
			probpattern=re.compile(r'(prob[0-9999]*)')
			newprobdir=re.search(probpattern,output)
			vars.jobpath.set(newprobdir.group(0).strip())
		except:
			self.Msg("Error running 'newprob' script.")

	def save(self):
		pars.Write_Pars(0)

	def saveas(self):
		pars.Write_Pars(1)

	def submit_job(self):
		sub.Submit(pars)

	#multiple routines here to avoid bindings
	def Read_Tracer_g(self): self.Geoframe.Read_Tracer('geo')
	def Read_Tracer_a(self): self.Geoframe.Read_Tracer('all')
	def Read_Tracer_s(self): self.Geoframe.Read_Tracer('spec')
	def Read_Tracer_r(self): self.Geoframe.Read_Tracer('ref')
	def Read_Tracer_p(self): self.Geoframe.Read_Tracer('params')




	def Expert(self):
		if not vars.expertV.get():
			# Remove expert options in Opframe
			[item.config(state=DISABLED) for item in self.Opframe.numgroup.interior().winfo_children() 
			 if 'state' in item.keys()]
			self.Opframe.switch_advanced()
			# Remove expert options in Inoutframe
			self.Inoutframe.istepomegaL.grid_remove()
			self.Inoutframe.istepomegaE.grid_remove()
			# Remove expert options in Domframe
			self.Domframe.switch_advanced()
			self.Genframe.switch_advanced()
			self.Collframe.switch_advanced()
		else:
			# Restore expert options in Opframe
			[item.config(state=NORMAL) for item in self.Opframe.numgroup.interior().winfo_children() 
			 if 'state' in item.keys()]
			self.Opframe.switch_advanced()

			# Restore expert options in Inoutframe		
			self.Inoutframe.istepomegaL.grid()
			self.Inoutframe.istepomegaE.grid()
			# Restore expert options in Domframe				
			self.Domframe.switch_advanced()
			# Restore expert options in Genframe
			self.Genframe.switch_advanced()
			# Restore expert options in Collframe
			self.Collframe.switch_advanced()
		self.switch_lin_nonlin()
		self.switch_nonlocal()
		self.switch_ev_iv()
			

#compiles the GENE documentation if necessary and uses acroread to display it
	def help(self):
		#os.chdir(self.startdir)
		testpath=os.path.join(self.startdir,"doc/gene.pdf")
		if os.path.exists(testpath): os.system("acroread "+testpath+" &")
		else:
			try:
				os.chdir("doc")
				os.system("pdflatex gene.tex")
				os.chdir("..")
				os.system("acroread doc/gene.pdf &")
			except: self.Msg("GENE documentation not available.")

	def quickintro(self):
		introwin=Toplevel(bg=vars.stdbg,width=250)
		introwin.title('3 steps to a GENE parameter file')
		Label(introwin,text="""1. On the 'Operations' tab, select the type of simulation, linear/nonlinear and local/nonlocal.\n2. Hit the 'Default values' button and confirm.\n3. Click 'Write Parameters' and select an output directory.""",justify=LEFT,bg=vars.stdbg).grid(padx=15,pady=15) 

	
#display a message in the lower part of the launcher window
	def Msg(self,message,highlight=0):
		self.MsgM.config(state=NORMAL)
		self.MsgM.insert("1.0",message+'\n')
		self.MsgM.config(state=DISABLED)
		if highlight:
			self.MsgM.config(bg='green')
		else:
			self.MsgM.config(bg=vars.stdbg)

#clear all messages
	def MsgClear(self):
		self.MsgM.config(state=NORMAL)
		self.MsgM.delete(1.0,END)
		self.MsgM.config(state=DISABLED)		



#main program    
root=Tk()

def_font=tkFont.Font(family='Bitstream Vera Sans',size=10)
root.option_add("*Font", def_font)
root.option_add("*Entry*Font", def_font)
root.option_add("*Text*Font", def_font)
app=Launcher(root)

root.title("GENE Launcher")
root.mainloop()









