import os
from subprocess import call, check_output
import math
import sys
from ConfigParser import SafeConfigParser

# global parameters
#mpicmd = "aprun"
#procpernode = 24
submitJobs = False

#appdir = "./scaling5_test"
#executable = "advection-example-faults"
parfile = "ctparam"

BASENAME = 'comp'
DIM = 6
LEVAL = [3,1,7,7,7,1]
LEVAL2 = [3,1,8,8,8,1]
P = [1,1,1,16,32,1]
NGROUP = 1
NPROCS = 512
DT = 5e-3
NCOMBI = 1
NSTEPS = 6000

def setctparam( parfile, lmin, lmax, leval, leval2, p, ngroups,\
                nprocs, ncombi, nsteps, dt ):
    # read parameter file
    parser = SafeConfigParser()
    parser.read( parfile )
       
    #update config
    parser.set( 'ct', 'lmin', '%s %s %s %s %s %s' % tuple(lmin) )
    parser.set( 'ct', 'lmax', '%s %s %s %s %s %s' % tuple(lmax) )
    parser.set( 'ct', 'leval', '%s %s %s %s %s %s' % tuple(leval) )
    parser.set( 'ct', 'leval2', '%s %s %s %s %s %s' % tuple(leval2) )
    parser.set( 'ct', 'p', '%s %s %s %s %s %s' % tuple(p)  )
    parser.set( 'manager', 'ngroup', str(ngroups) )
    parser.set( 'manager', 'nprocs', str(nprocs) )
    parser.set( 'ct', 'ncombi', str(ncombi)  )
    parser.set( 'application', 'nsteps', str(nsteps)  )
    parser.set( 'application', 'dt', str(dt)  )
    
    # set paths for plot files
    workdir = os.getcwd()
    plotfile1 = workdir + '/ginstance/plot.dat'
    plotfile2 = workdir + '/ginstance/plot2.dat'
    parser.set( 'ct', 'fg_file_path', plotfile1 )
    parser.set( 'ct', 'fg_file_path2', plotfile2 )
    
    cfgfile = open( parfile,'w')
    parser.write(cfgfile)
    cfgfile.close()


def createExperiment( dirname, lmin, lmax, leval, leval2, \
                      p, ngroup, nprocs, ncombi, nsteps, dt ):   
    # make file name
    nsteps_total = ncombi*nsteps
             
    if os.path.exists(dirname):
        print 'component already exists'
        return
        
    # create dir, copy files and go into dir
    path = "./" + dirname 
    call( ['mkdir', path] )
    call( ['cp', '-r', 'template', path ] )
    call( ['cp', 'manager', path ] )
    call( ['cp', 'preproc.py', path ] )
    call( ['cp', 'run.sh', path ] )
    call( ['cp', './ctparam', path ] )
    os.chdir( path )
    
    #write ctparam
    setctparam( parfile, lmin, lmax, leval, leval2, p, ngroup,\
                nprocs, ncombi, nsteps, dt )
    
    # start job
    if submitJobs:
        # write jobscript etc
        call( ["qsub", "./job.sub"] )
    else:
        call( ["./run.sh"] )
    
    # leave dir
    os.chdir('..')


if __name__ == "__main__":
    print sys.argv 

    #read filename from args
    expdir = sys.argv[1]
    
    errorMode = sys.argv[2]
    errorPlotFile = sys.argv[3]
    errorOutFile = sys.argv[4]
    expPlotFileType = sys.argv[5]
    
    
    filename = expdir + '/ginstance/spaces.dat'
    
    #loop over lines in file
    for line in open(filename):
        par = line.split()
        
        # skip empty lines
        if len(par) == 0:
            continue
        
        # read lmin
        lmin = [ int(x) for x in par[1:7] ]
        lmax = lmin
        
        basename    = BASENAME
        dim         = DIM
    
        leval       = LEVAL
        leval2      = LEVAL2
        p           = P
        
        ngroup      = NGROUP
        nprocs      = NPROCS     
        
        ncombi      = NCOMBI
        nsteps      = NSTEPS
        dt          = DT
        
        dirname = basename 
        dirname += "_" + '%s%s%s%s%s%s' % tuple(lmin)
        dirname += "_" + str(NCOMBI*NSTEPS)
	dirname += "_" + '%s%s%s%s%s%s' % tuple(leval)
            
        print dirname, lmin, lmax, leval, leval2, p, ngroup,\
                nprocs, ncombi, nsteps, dt
        
        createExperiment( dirname, lmin, lmax, leval, leval2, \
                          p, ngroup, nprocs, ncombi, nsteps, dt )
        
        # calc error
        if expPlotFileType == '1':
            expPlotFile = dirname + '/ginstance/plot.dat'
        else:
            expPlotFile = dirname + '/ginstance/plot2.dat'
            
        print './errorCalc %s %s %s %s %s' % (errorMode, errorPlotFile, expPlotFile, errorOutFile, dirname)
        call( [ 'aprun','./errorCalc', errorMode, errorPlotFile, expPlotFile, errorOutFile, dirname] )
        
       
