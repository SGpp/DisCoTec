#path to python interface
from sys import path
path.append('/import/home_local/oberstei/Documents/ExaHD/gene_python_interface_clean/src')
path.append('/import/home_local/oberstei/Documents/ExaHD/gene_python_interface_clean/src/tools')

from ConfigParser import SafeConfigParser
import collections
from subprocess import call

# python interface stuff
import ActiveSetFactory as aSF
import combinationScheme as cS
import numpy as np

# some static definitions
spcfile = 'spaces.dat'
        
# read parameter file
parser = SafeConfigParser()
parser.read('ctparam')

config = collections.namedtuple('Config', 'lmin lmax ntimesteps_total dt_max ntimesteps_combi basename executable mpi startscript sgpplib tasklib ngroup nprocs shat kymin lx')

config.lmin = [int(x) for x in parser.get('ct','lmin').split()]
config.lmax = [int(x) for x in parser.get('ct','lmax').split()]
config.leval = [int(x) for x in parser.get('ct','leval').split()]
config.p = [int(x) for x in parser.get('ct','p').split()]
config.ntimesteps_total = int( parser.get('application', 'nsteps') )
config.dt_max = float( parser.get('application', 'dt') )
config.basename = parser.get('preproc','basename')
config.executable = parser.get('preproc','executable')
config.mpi = parser.get('preproc','mpi')
config.startscript = parser.get('preproc','startscript')
config.sgpplib = parser.get('preproc','sgpplib')
config.tasklib = parser.get('preproc','tasklib')
config.ngroup = int( parser.get('manager','ngroup') )
config.nprocs = int( parser.get('manager','nprocs') )
config.istep_omega = 10
#config.shat = parser.get('application','shat') 
config.kymin = parser.get('application','kymin') 
config.lx = parser.get('application','lx') 
config.numspecies = parser.get('application','numspecies') 
config.local = parser.get('application','GENE_local') 
config.nonlinear = parser.get('application','GENE_nonlinear') 
if config.local == "T" :
    config.shat = parser.get('application','shat') 

# if command line options given overwrite config options
'''
import sys
if( len( sys.argv ) > 1 ):
    assert( len( sys.argv ) == 5 )
    
    config.ntimesteps_total = int( sys.argv[1] )
    config.dt_max = float( sys.argv[2] )
    config.ntimesteps_combi = int( sys.argv[3] )
    config.ntimesteps_ev_calc = int( sys.argv[4] )
    
    #update config
    parser.set( 'ct', 'ntimesteps_total', str( config.ntimesteps_total ) )
    parser.set( 'ct', 'dt_max', str( config.dt_max ) )
    parser.set( 'ct', 'ntimesteps_combi', str( config.ntimesteps_combi ) )
    parser.set( 'ct', 'ntimesteps_ev_calc', str( config.ntimesteps_ev_calc ) )
    cfgfile = open('./ctparam','w')
    parser.write(cfgfile)
    cfgfile.close()
'''

# some sanity checks
lmin = config.lmin
lmax = config.lmax
leval = config.leval

if( not (len(lmin) == len(lmax) == len(leval) ) ):
    raise ValueError('config: lmin,lmax,leval not of same size')

for i in range(0,len(lmin)):
    if lmin[i] < 0 or lmax[i] < 0 or leval[i] < 0:
        raise ValueError('config: lmin,lmax,leval may not be negative')
    
    if lmin[i] > lmax[i]:
        raise ValueError('config: lmin may not be larger than lmax') 
    
'''
if( not (config.ntimesteps_total%config.ntimesteps_combi == 0) ):
    raise ValueError('config: ntimesteps_total must be a multiple of ntimesteps_combi')
    
if( config.ntimesteps_combi > config.ntimesteps_ev_calc ):
    if( not (config.ntimesteps_combi%config.ntimesteps_ev_calc == 0) ):
        raise ValueError('config: ntimesteps_combi must be a multiple of ntimesteps_ev_calc')
else:
    if( not (config.ntimesteps_total%config.ntimesteps_combi == 0) ):
        raise ValueError('config: ntimesteps_ev_calc must be a multiple of ntimesteps_combi')
'''
            
# create base folder
call(["cp","-r","template",config.basename])

# set ct specific parameters in base folder
p = config.p
px = p[0]
py = p[1]
pz = p[2]
pv = p[3]
pw = p[4]
ps = p[5]

parfile = './' + config.basename + "/parameters"
pfilein = open(parfile,'r')
pin = pfilein.read()
pfilein.close()
pout = pin.replace('$ngroup',str(config.ngroup))
pout = pout.replace('$nprocs',str(config.nprocs))
pout = pout.replace('$istep_omega', str(config.istep_omega))
pout = pout.replace('$ps',str(ps),1)
pout = pout.replace('$pv',str(pv),1)
pout = pout.replace('$pw',str(pw),1)
pout = pout.replace('$px',str(px),1)
pout = pout.replace('$py',str(py),1)
pout = pout.replace('$pz',str(pz),1)
pfileout = open(parfile,'w')
pfileout.write(pout)
pfileout.close()

# create combischeme
# note that the level vector in the python vector is stored in reverse order
lminp = lmin[::-1]
lmaxp = lmax[::-1]
#activeSet = aSF.ClassicDiagonalActiveSet(lmaxp,lminp)
#scheme = cS.combinationSchemeArbitrary(activeSet.getActiveSet())
factory = aSF.ClassicDiagonalActiveSet(lmaxp,lminp,0)
activeSet = factory.getActiveSet()
scheme = cS.combinationSchemeFaultTolerant(factory)
# detect number of simulation steps
#nsteps = config.ntimesteps_combi if config.ntimesteps_combi <= config.ntimesteps_ev_calc else config.ntimesteps_ev_calc

# loop over scheme
id = 0
spaces = ''
print len(scheme.getCombinationDictionary())
for l in scheme.getCombinationDictionary():
    # note that the level vector in the python vector is stored in reverse order
    l0 = l[5]
    l1 = l[4]
    l2 = l[3]
    l3 = l[2]
    l4 = l[1]
    l5 = l[0]
     
    # append subspaces entry
    spaces  += str(id) + " " +  str(l0) + ' ' + str(l1) + ' ' + str(l2) \
            + ' ' + str(l3) + ' ' + str(l4) + ' ' + str(l5) + ' ' \
            + str(scheme.getCombinationDictionary()[l]) + "\n"
             
    # copy template folder
    call(["cp","-r","./template",'./' + config.basename + str(id)])
    
    # set ct specific parameters
    parfile = './' + config.basename + str(id) + "/parameters"
    pfilein = open(parfile,'r')
    pin = pfilein.read()
    pfilein.close()
    
    pout = pin.replace('$nx0',str(2**l0),1)
    if config.nonlinear == "T" :
    	pout = pout.replace('$nky0',str(2**l1),1)
    else:
        pout = pout.replace('$nky0',str(2**l1-1),1)
    pout = pout.replace('$nz0',str(2**l2),1)
    pout = pout.replace('$nv0',str(2**l3),1)
    pout = pout.replace('$nw0',str(2**l4),1)
    pout = pout.replace('$nspec',str(config.numspecies),1)
    pout = pout.replace('$GENE_local',str(config.local),1)
    pout = pout.replace('$GENE_nonlinear',str(config.nonlinear),1)

    pout = pout.replace('$ps',str(ps),1)
    pout = pout.replace('$pv',str(pv),1)
    pout = pout.replace('$pw',str(pw),1)
    pout = pout.replace('$px',str(px),1)
    pout = pout.replace('$py',str(py),1)
    pout = pout.replace('$pz',str(pz),1)

    pout = pout.replace('$ngroup',str(config.ngroup))
    pout = pout.replace('$nprocs',str(config.nprocs))
    
    pout = pout.replace('$ntimesteps_combi', str(config.ntimesteps_total))
    pout = pout.replace('$istep_omega', str(config.istep_omega))
    pout = pout.replace('$dt_max', str(config.dt_max))
    if config.local == "T" :
        pout = pout.replace('$shat', str(config.shat))
    pout = pout.replace('$kymin', str(config.kymin))
    pout = pout.replace('$lx', str(config.lx))

    pout = pout.replace('$read_cp','F')
    pout = pout.replace('$write_cp','T')
    
    pfileout = open(parfile,'w')
    pfileout.write(pout)
    pfileout.close()
    
    # update counter
    id += 1
                            
# print spaces file
sfile = open('./' + config.basename + '/' + spcfile,'w')
sfile.write(spaces)
sfile.close()

# link manager executable to base folder
call(["ln","-s","../manager",'./' + config.basename + '/manager'])

# copy param file to base folder 
call(["cp","./ctparam",'./' + config.basename + '/'])

# create start script in base folder
scmd = "export LD_LIBRARY_PATH=" + config.sgpplib + ":" + config.tasklib + ":/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH\n"
#scmd = "source xprtld.bat\n"
scmd += config.mpi
scmd += " -n " + str(config.nprocs*config.ngroup) + ' ' + config.executable
scmd += " : "
scmd += " -n 1" + " ./manager"

sfile = open( config.basename + "/" + config.startscript,'w')
sfile.write(scmd)
sfile.close()


