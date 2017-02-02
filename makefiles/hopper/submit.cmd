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
aprun -n 256 ./gene_hopper

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 256 --ppn 1 --mps 4
