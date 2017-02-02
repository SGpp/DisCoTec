from Tkinter import *
from tkSimpleDialog import askstring,askinteger
from tkMessageBox import askokcancel
from subprocess import Popen,PIPE
import os, string

class Submit:

    #initialize variables for multiple-job management
    def __init__(self,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod


    def Timelimit(self,maxtime):
        #maxtime is wallclock timelimit in hours
        timelimit=askstring('Wall clock limit','(hh:mm) or (s)',initialvalue=str(maxtime)+':00')
        if timelimit==None: return 0
        t=timelimit.split(':')
        if len(t)>1:
            sec=int(t[0])*3600+int(t[1])*60
        elif len(t)==1 and int(t[0])>200:
            sec=t[0]
        elif len(t)==1 and int(t[0])<200: 
            sec=int(t[0])*3600
        else: return 0
        if sec>maxtime*3600:
            launcher.Msg('This machine allows only '+str(maxtime)+' hours of runtime, decreasing time limit to this value.')
            timelimit=str(maxtime)+':00:00'
        else: 
            timelimit=str(sec/3600)+':'+str((sec%3600)/60)+':00'
        return timelimit

    def Submit(self,pars):
        if vars.job_saved[vars.currentjob]:
            #Check for recent changes
            pars.Write_Pars(2)
            os.chdir(vars.jobpath.get())
            pipe=Popen('diff parameters temp_pars', shell=True, stdout=PIPE).stdout
            saverecent=0
            if len(pipe.readlines())!=0:
                saverecent=askokcancel(message='Save recent changes to parameters?')
            if saverecent: pars.Write_Pars(3)
            pipe=Popen('rm '+os.path.join(launcher.startdir,vars.jobpath.get()+'/temp_pars'),shell=True,stdout=PIPE).stdout
            os.chdir(launcher.startdir)
        
            #Check for scans in parameters
            parfile=open(os.path.join(vars.jobpath.get(),'parameters'))
            parlines=parfile.readlines()
            found_scan=0
            for line in parlines:
                if '!scan' in line: found_scan=1
            if found_scan and not vars.scan.get():
                vars.scan.set(askokcancel(message='Found "scan" keyword in output parameters; submit a parameter scan?'))
            elif not found_scan and vars.scan.get():
                launcher.Msg('No "scan" keyword found in output parameters; submitting single run.')
                vars.scan.set(0)
            if vars.scan.get(): pars.append_scan_nml()
            self.Read_Submitscript()
            self.Write_Submitscript()
        else:
            launcher.Msg('Please save parameters first!')


    def Check_Exec(self,dir,name):
        if os.path.exists(os.path.join(dir,'gene_'+name)): 
            return TRUE
        else: 
            return FALSE
        

    def Read_Submitscript(self):
        try:
            self.submitfile=open(os.path.join(launcher.startdir+'/makefiles/',vars.machine+'/launcher.cmd'))
        except:
            launcher.Msg('Error opening '+str(os.path.join(launcher.startdir+'/makefiles/',vars.machine+'/launcher.cmd')))
            return
        self.submitscript=self.submitfile.readlines()
        self.maxwalltime='0'
        self.submitcmd=''
        self.procspernode=0
        for line in self.submitscript:
            if 'MAXWT' in line:
                self.maxwalltime=line.split()[2]
            if 'PROCSPERNODE' in line:
                self.procspernode=int(line.split()[2])
            #timeformat not yet implemented, is supposed to allow for e.g. hh:mm:ss or s
            #may not be necessary if loadlevelers recognize both
            if 'TIMEFORMAT' in line:
                self.timeformat='s'
            if 'SUBMITCMD' in line:
                self.submitcmd=line.split()[2]
        if self.maxwalltime=='0': self.maxwalltime='24'
        if self.submitcmd=='': 
            launcher.Msg('Error: Submit command is not defined in '+str(os.path.join(launcher.startdir+'/makefiles/',vars.machine+'/launcher.cmd')))
            return

    def Write_Submitscript(self):
        try: (paroutdir,paroutfile)=os.path.split(vars.parout)
        except:
            launcher.Msg('Please write parameters before submitting.')
            return
        submitfilename='launcher.cmd'
        if not self.Check_Exec(paroutdir,vars.machine): 
            launcher.Msg('Executable does not exist, please compile before submitting!')
            return
        self.submitfile=open(os.path.join(paroutdir,submitfilename),"w")
        #total processes
        nprocs=int(vars.nprocsV.get())
        #processes per node, currently only fully used nodes are allowed
        if self.procspernode!=0:
            if nprocs%self.procspernode==0: 
                n_nodes=nprocs/self.procspernode
            else:
                if nprocs/self.procspernode==0:
                    n_nodes=1
                    self.procspernode=nprocs
                else:
                    launcher.Msg('Please use a multiple of '+str(self.procspernode)+' processes.')
                    return
        timelimit=self.Timelimit(self.maxwalltime)
        if timelimit==0: return
        try: 
            jobname=askstring('Job name','',initialvalue='GENE')
            if jobname==None: return #user clicked 'cancel'
        except: return

        replacelist=['NMPIPROCS','NPROCSDIV4','PPN','NODES','WALLTIME','JOBNAME']
        valuedict={'NMPIPROCS': str(vars.nprocsV.get()), 'NPROCSDIV4': str(vars.nprocsV.get()/4),
                   'PPN': str(self.procspernode), 'NODES': str(n_nodes), 'WALLTIME': str(timelimit),
                   'JOBNAME': jobname}

        for line in range(len(self.submitscript)):
            for item in replacelist:
                if item in self.submitscript[line]:
                    #some last minute checks on requirements for job submission
                    if item=='PPN' and self.procspernode==0:
                        launcher.Msg('Error: Processes per node is set to zero. Please check the settings in '
                                     +str(os.path.join(launcher.startdir+'/makefiles/',vars.machine+'/launcher.cmd')))
                        return
                    if item=='NMPIPROCSDIV4' and vars.nprocsV.get()%4!=0:
                        launcher.Msg('Error: Number of processes has to be a multiple of 4.')
                        return
                    #replace keyword definitions in the submit script with the appropriate values
                    self.submitscript[line]=string.replace(self.submitscript[line],item,str(valuedict[item]))

        #find line with execution command (ignore trailing empty lines)
        line=len(self.submitscript)-1
        while self.submitscript[line]=='\n':
            line-=1
        execline=line
        for line in self.submitscript[0:execline]:
            self.submitfile.write(line)
        if vars.scan.get():
            self.submitfile.write('./scanscript --n_pes='+str(vars.nprocsV.get())+' --procs_per_node='+str(self.procspernode)+'\n')
        else:
            self.submitfile.write(self.submitscript[execline].strip('\n')+'\n')
        self.submitfile.close()
        os.chdir(paroutdir)
        pipe=Popen(self.submitcmd+' launcher.cmd', shell=True, stdout=PIPE).stdout
        launcher.Msg(''.join(pipe.readlines()))
        os.chdir(launcher.startdir)
     
