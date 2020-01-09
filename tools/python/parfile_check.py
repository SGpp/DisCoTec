#! /usr/bin/env python
from ParIO import Parameters
import os

# To do:
# check that grid points are divisible by nprocs

#### User Input ####
parfile = 'parameters'
#### User Input ####

def check_batch_script():
    if os.getenv('MACHINE') == 'stampede':
        f=open('submit.cmd')
        sfile = f.read()
        sfile_lines = sfile.split('\n')
        nprocs = -1
        nprocs2 = -1
        for line in sfile_lines:
            if '-n' in line:
                line_split = line.split()
                for i in range(len(line_split)):
                    if line_split[i] == '-n':
                        nprocs = int(line_split[i+1])
                        break
            if './scanscript' == line[:12]:
                line_split = line.split()
                for i in range(len(line_split)):
                    if '--np' in line_split[i]:
                        nprocs2 = int(line_split[i+1])
                        break
        if nprocs == -1:
            print "Error: can't determine number of processors from submit.cmd"
        elif nprocs2 != -1:
            if nprocs != nprocs2:    
                print "Error: SBATCH nprocs differs from scanscript nprocs."
        return nprocs, False
    elif 'helios' in os.getenv('HOST'):
        f=open('submit.cmd')
        sfile = f.read()
        sfile_lines = sfile.split('\n')
        nprocs = -1
        nprocs2 = -1
        for line in sfile_lines:
            if '--ntasks' in line:
                line = line.replace("=", " ")
                ls = line.split()
                line_split = []
                for l in ls:
                    if l != ' ' and l != '':
                        line_split.append(l)
                for i in range(len(line_split)):
                    if 'ntasks' in line_split[i]:
                        nprocs = int(line_split[i+1])
                        break
            elif '-A' in line:
                print "Check allocation:"
                print line
                dummy = raw_input("Press any key to continue.")
    elif 'edison' in os.getenv('HOST'):
        f=open('submit.cmd')
        sfile = f.read()
        sfile_lines = sfile.split('\n')
        nprocs = -1
        nprocs2 = -1
        for line in sfile_lines:
            if 'mppwidth' in line:
                line = line.replace("=", " ")
                ls = line.split()
                line_split = []
                for l in ls:
                    if l != ' ' and l != '':
                        line_split.append(l)
                for i in range(len(line_split)):
                    if 'mppwidth' in line_split[i]:
                        nprocs = int(line_split[i+1])
                        break
        if nprocs == -1:
            print "Error: can't determine number of processors from submit.cmd"
            return nprocs, True
        return nprocs, False

    else:
        print "Error: can't determine machine."
        return -1, True
 
par = Parameters()
par.Read_Pars(parfile)
pars = par.pardict
pars['diagdir'] = pars['diagdir'].strip()[1:-1]

# print "diagdir", pars['diagdir']
true_strings = ['t','.t.','true']
# print "read_checkpoint", pars['read_checkpoint']
# print pars['read_checkpoint'].lower()
# print true_strings
# print pars['read_checkpoint'].lower() in true_strings
pars['read_checkpoint'] = pars['read_checkpoint'].strip().lower() in true_strings
# print "read_checkpoint",pars['read_checkpoint']

has_error = False
nprocs, has_error = check_batch_script()

nprocs_par = 1
if pars['n_procs_s'] > 0:
    nprocs_par *= pars['n_procs_s']
if pars['n_procs_v'] > 0:
    nprocs_par *= pars['n_procs_v']
if pars['n_procs_w'] > 0:
    nprocs_par *= pars['n_procs_w']
if pars['n_procs_x'] > 0:
    nprocs_par *= pars['n_procs_x']
if pars['n_procs_y'] > 0:
    nprocs_par *= pars['n_procs_y']
if pars['n_procs_z'] > 0:
    nprocs_par *= pars['n_procs_z']

if nprocs_par > nprocs:
    print "Error: parameters file defines more processors than submit file."
    has_error = True
if nprocs % nprocs_par != 0:
    print "Error: nprocs is not divisible by processors from parameters file."
    has_error = True

if nprocs_par > pars['n_procs_sim']:
    print "Error: n_procs_*s inconsistent with n_procs_sim"
    has_error = True
if nprocs % pars['n_procs_sim'] != 0:
    print "Error: n_procs_*s inconsistent with n_procs_sim"
    has_error = True

if pars['n_procs_sim'] and pars['n_parallel_sims']:
    if pars['n_procs_sim']*pars['n_parallel_sims'] != nprocs:
        print "Error: n_procs_sim * n_parallel_sims does not equal nprocs from submit.cmd"
        has_error = True
elif pars['n_procs_sim']:
    if pars['n_procs_sim'] != nprocs:
        print "Error: n_procs_sim does not equal nprocs from submit.cmd"
        has_error = True

if not os.path.isdir(pars['diagdir']):
    print "Error: diagdir doesn't exist."
    print "diagdir:",pars['diagdir']
    has_error = True

if pars['read_checkpoint']:
    if not os.path.exists(pars['diagdir']+'/checkpoint'):
        print "Error: checkpoint does not exist."
        has_error = True
    else:
        chp_size = os.path.getsize(pars['diagdir']+'/checkpoint')
        if chp_size < 100:
            print "Error: checkpoint file appears to be empty."
            has_error = True

if not has_error:
    print "No errors detected."
else:
    print "Fix errors before submitting!"
