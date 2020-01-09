##submit with: qsub submit.cmd
##in the directory containing
##the gene executable
#PBS -q regular
#PBS -l mppwidth=256
#PBS -l walltime=04:00:00
#PBS -j eo
#PBS -V
#PBS -m abe

cd $PBS_O_WORKDIR

export NPROC=`qstat -f $PBS_JOBID | awk '/mppwidth/ {print $3}'`
export NTASK=`qstat -f $PBS_JOBID | awk '/mppnppn/  {print $3}'`

export OMP_NUM_THREADS=1

aprun -n $NPROC ./gene_edison

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np $NPROC --ppn 24 --mps 4
