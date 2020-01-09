#contains machine dependent routines
import os
from subprocess import Popen,PIPE

class machine:
    def __init__(self,launch,varmod):
        global vars,launcher
        launcher=launch
        vars=varmod
        try:
            pipe=Popen("make mach_wrapper", shell=True, stdout=PIPE).stdout
            vars.machine=pipe.readlines()[-1].strip('\n')
        except:
            vars.machine='new_machine'
        try: vars.username=os.environ['USER']
        except: vars.username=''
        try: vars.outdir=os.environ['WORK']
        except:
            file=open(os.path.join(launcher.startdir,'launcher/output-dirs.txt'),'r')
            dirlist=file.readlines()
            file.close()
            machine_dirs={}
            for i in range(len(dirlist)):
                try:
                    machine_dirs[dirlist[i].split()[0]]=dirlist[i].split()[1]
                except:
                    pass
            if vars.machine not in machine_dirs.keys():
                vars.outdir='<Please enter output directory>'
            else:
                vars.outdir=machine_dirs[vars.machine]+vars.username
