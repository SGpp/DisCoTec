#!/usr/bin/env python

import subprocess
import os
import re
import sys

# First setup the test configuration as a list
class CompileEnvironment_t(object):
    def __init__(self,name):
        self.machine_name = name
        self.host_name = ""
        self.precision=['double']
        self.futils=['yes']
        self.slepc=['no']
        self.mklversion=['10.3']
        self.modules=[]

    def clone(self):
        newCE=CompileEnvironment_t(self.machine_name)
        newCE.host_name = self.host_name
        newCE.precision = self.precision[:]
        newCE.futils = self.futils[:]
        newCE.slepc = self.slepc[:]
        newCE.modules = self.modules[:]
        return newCE

    def output(self):
        print "[%s]"%self.machine_name
        print "\thost_name = ",self.host_name
        print "\tmodules = ",self.modules
        print "\tslepc = ",self.slepc
        print "\tfutils = ",self.futils
        print "\tprecision = ",self.precision


# ----------------------- start of main program -------

from optparse import OptionParser
run_remote = True

# read commandline options
parser=OptionParser()
parser.add_option("-r","--remote",action="store_true",dest="run_remote",default=True,
                  help="run with ssh connection to the configured host.")
(options,args)=parser.parse_args()


print "args = ",args
ListOfCEs=[]
config_file=open("compiler_testsuite.cfg","rt")
for line in config_file:
    if (re.match(r"^\s*$",line) or re.match(r"^#",line)):
        continue

    #print line,
    mo=re.match(r"\[([a-zA-z0-9_-]+)\]",line)
    if mo:
        # start a new compile environment
        ListOfCEs.append(CompileEnvironment_t(mo.group(1)))
    else:
        # it is a key=value pair
        parts = line.split("=")
        key=parts[0]
        value=parts[1]
        if key=="host_name":
            ListOfCEs[-1].host_name=value.strip()
        elif key=="modules":
            ListOfCEs[-1].modules=[x.strip() for x in value.split()]
        elif key=="futils":
            ListOfCEs[-1].futils=[x.strip() for x in value.split()]
        elif key=="slepc":
            ListOfCEs[-1].slepc=[x.strip() for x in value.split()]
        elif key=="mklversion":
            ListOfCEs[-1].mklversion=[x.strip() for x in value.split()]
        elif key=="precision":
            ListOfCEs[-1].precision=[x.strip() for x in value.split()]

config_file.close()
print "ListOfCEs before expansion:"
for x in ListOfCEs:
    x.output()
    print "\n"

print "Now expanding multiple entries."
listindex = 0
for CE in ListOfCEs:
    if len(CE.slepc)!=1:
        newCE=CE.clone()
        newCE.slepc[0]=CE.slepc[1]
        del CE.slepc[1]
        ListOfCEs.insert(listindex+1,newCE)

    if len(CE.futils) != 1:
        newCE=CE.clone()
        newCE.futils=(CE.futils[1])
        CE.futils=(CE.futils[0])
        ListOfCEs.insert(listindex+1,newCE)

    if len(CE.precision) != 1:
        newCE=CE.clone()
        del newCE.precision[0]
        del CE.precision[1]
        ListOfCEs.insert(listindex+1,newCE)
        
    listindex += 1

print "ListOfCEs after expansion:"
for x in ListOfCEs:
    x.output()
    print "\n"

#sys.exit(1)

cwd = os.getcwd()
previous_machine=None
counter = 0
for CE in ListOfCEs:
    if len(args)!=0:
        if (CE.machine_name not in args):
            continue
    print "\n----------------------------"
    print "Compiling with MACHINE = ",CE.machine_name
    # build a shell script to start
    skriptname = "./to_run_%2.2u.sh"%counter
    shell_script=open(skriptname,"wt")
    shell_script.write("#!/bin/bash --login\n")
    shell_script.write("cd %s\n"%cwd)
    shell_script.write("export MACHINE=%s\n"%CE.machine_name)
    shell_script.write('echo "purging modules"\n')
    shell_script.write("module purge\n")
    shell_script.write("module load %s\n"%' '.join(CE.modules))
    shell_script.write("module list\n")
    #if (previous_machine==CE.machine_name):
    shell_script.write("make -f ../makefile distclean\n")
    shell_script.write("make ")
    # setting makefile variables
    shell_script.write("SLEPC=%s "%CE.slepc[0])
    shell_script.write("FUTILS=%s "%CE.futils[0])
    shell_script.write("PRECISION=%s "%CE.precision[0])
    shell_script.write("MKLVERSION=%s "%CE.mklversion[0])
    shell_script.write(" -f ../makefile >make_output_%2.2u 2>&1\n"%counter)
    shell_script.write('if [ $? -ne 0 ]; then echo "Error occurred.";fi\n')
    shell_script.write('mv ../bin/gene_%s ../bin/gene_%s_%2.2u'%(CE.machine_name,CE.machine_name,counter))
    shell_script.close()
    os.chmod(skriptname,0755)
    subprocess.call(["scp",skriptname,"%s:"%CE.host_name],shell=False)
    if run_remote:
        subprocess.call(["ssh",CE.host_name,skriptname],shell=False)
    else:
        subprocess.call([skriptname],shell=False)

    previous_machine=CE.machine_name
    counter += 1
