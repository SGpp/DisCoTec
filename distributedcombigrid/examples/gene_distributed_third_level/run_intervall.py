import os
from subprocess import call, check_output
import math
import sys
from ConfigParser import SafeConfigParser
from macpath import basename

# global parameters
#mpicmd = "aprun"
#procpernode = 24
submitJobs = False

#appdir = "./scaling5_test"
#executable = "advection-example-faults"
parfile = "ctparam"


# walltime has to be given in minutes
def writeJobscriptHornet( filename, jobname, nodes, walltime, cmd ):
    walltime_minutes = int(walltime)
        
    walltime_hours = walltime_minutes / 60
    walltime_minutes = walltime_minutes % 60
        
    walltime_str = str(walltime_hours) + ":"
    
    if walltime_minutes > 9:     
        walltime_str += str(walltime_minutes) + ":00"
    else:
        walltime_str += "0" + str(walltime_minutes) + ":00"
    
    # header
    PPN = procpernode
    header  = "#!/bin/bash\n"
    header += "##Name of job\n"
    header += "#PBS -N " + str(jobname) + "\n"
    header += "##number of nodes and procs per node\n"
    header += "#PBS -l nodes=" + str(nodes) + ":ppn=24\n"
    header += "##Maximum execution time\n"
    header += "#PBS -l walltime=" + walltime_str + "\n"
    
    if( nodes > 2731 ):
        header += "#PBS -q nolimit \n"
    
    # body
    body = "cd $PBS_O_WORKDIR \n"
    #body += "source ../xprt.bat \n"
    body += cmd + "\n"
    
    # write to file
    content = header + "\n" + body
    file = open(filename,'w')
    file.write(content)
    file.close()


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
    '''
    if submitJobs:
        # write jobscript etc
        call( ["qsub", "./job.sub"] )
    else:
        call( ["./run.sh"] )
    '''
    
    # leave dir
    os.chdir( '..' )


if __name__ == "__main__":
    #read filename from args
    filename = sys.argv[1]
    
    errorPlotFileLarge = sys.argv[2]
    errorPlotFileSmall = sys.argv[3]
    
    #loop over lines in file
    for line in open(filename):
        par = line.split()
        
        # skip empty lines
        if len(par) == 0:
            continue
        
        basename    = par[0]
        dim         = par[1] 
        
        lmin        = [ int(x) for x in par[2:8] ]
        lmax        = [ int(x) for x in par[8:14] ]
        leval       = [ int(x) for x in par[14:20] ]
        leval2      = [ int(x) for x in par[20:26] ]
        p           = [ int(x) for x in par[26:32] ]
        
        ngroup      = int(par[32])
        nprocs      = int(par[33])     
        
        ncombi      = int(par[34])
        nsteps      = int(par[35])
        dt          = float(par[36])
        
        dirname = basename 
        dirname += "_" + '%s%s%s%s%s%s' % tuple(lmax)
        dirname += "_" + '%s%s%s%s%s%s' % tuple(lmin)
        dirname += "_" + str(nsteps*ncombi)
        dirname += "_" + str(nsteps)
        
        print basename, lmin, lmax, leval, leval2, p, ngroup,\
                nprocs, ncombi, nsteps, dt
                
        print dirname
        
        createExperiment( dirname, lmin, lmax, leval, leval2, \
                          p, ngroup, nprocs, ncombi, nsteps, dt )
        
        # calc ref31888 error
        expPlotFile = dirname + '/ginstance/plot2.dat'
        errorOutFile = basename
        errorOutFile += "_" + '%s%s%s%s%s%s' % tuple(lmax)
        errorOutFile += "_" + '%s%s%s%s%s%s' % tuple(lmin)
        errorOutFile += '_ref31888.err'
        errorMode = 'gf'
        print './errorCalc %s %s %s %s %s' % (errorMode, errorPlotFileLarge, expPlotFile, errorOutFile, nsteps)
        #call( [ 'aprun','./errorCalc', errorMode, errorPlotFileLarge, expPlotFile, errorOutFile, nsteps] )
        
        # calc error to same lmax
        expPlotFile = dirname + '/ginstance/plot.dat'
        errorOutFile = basename
        errorOutFile += "_" + '%s%s%s%s%s%s' % tuple(lmax)
        errorOutFile += "_" + '%s%s%s%s%s%s' % tuple(lmin)
        errorOutFile += '_lmax.err'
        errorMode = 'ff'
        print './errorCalc %s %s %s %s %s' % (errorMode, errorPlotFileSmall, expPlotFile, errorOutFile, nsteps)
        #call( [ 'aprun','./errorCalc',errorMode , errorPlotFileSmall, expPlotFile, errorOutFile, nsteps] )
            
        