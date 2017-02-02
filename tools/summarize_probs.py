#!/usr/bin/env python
import os
from time import localtime, strftime
from sys import exit

#check true/false
def istrue(par):
    return par in ['.t.','t','T','.T.']
#strip quotation marks ' and "
def sqt(par):
    return par.strip('"'"'")

#list of all local files and subdirectories
files=os.listdir('.')
#extract list of all prob directories
probs=[item for item in files if item[0:4]=='prob']
if len(probs)==0:
    print 'No "prob" subdirectories found, please run this script from a GENE directory.'
    exit()
#extract data from each probXX/parameters
for prob in probs:
    #par=Parameters()
    parfile=os.path.join(prob,'parameters')
    parfilehandle=open(parfile,'r')
    lines=parfilehandle.readlines()
    #generate a dictionary containing all parameters
    #warning: species parameters are not treated correctly using
    #this approach and will contain only the values of the last species
    p={}
    for line in lines:
        try:
            par=line.split('=')
            if len(line.split('='))==2:
                p[par[0].strip()]=par[1].strip()
        except:
            pass

    #extract a general summary of the simulation    
    print 'Directory %s:' %(prob)
    if istrue(p['nonlinear']):
        nl='nonlinear'
    else:
        nl='linear'
    if 'x_local' in p.keys():
        if istrue(p['x_local']):
            loc='Local'
        else:
            loc='Global'
    else:
        loc='Local'
    geom=sqt(p['magn_geometry'])
    res='%s x %s x %s x %s x %s' %(p['nx0'],p['nky0'],p['nz0'],p['nv0'],p['nw0'])
    sp=p['n_spec']
    if 'collision_op' in p.keys():
        if sqt(p['collision_op'])=='none':
            coll='Collisionless'
        else:
            coll='Collisional'
    else:
        coll='Collisionless'
    outdir=sqt(p['diagdir'])
    #formatted date/time of last parameters modification
    timestamp=strftime("%a, %d %b %Y %H:%M",localtime(os.path.getmtime(parfile)))
    #print everything
    print '   %s %s, %s species, %s grid points' %(loc,nl,sp,res)
    print '   %s run, %s geometry' %(coll,geom)
    print '   Output dir: %s' %(outdir)
    print '   Last edit: %s' %(timestamp)
    print ''
