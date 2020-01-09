#PBS -o output.dat
#PBS -l walltime=01:00:00,nodes=1:ppn=12
#PBS -A WWUadmin
#PBS -M <enter your mail address>
#PBS -m ae
#PBS -q default
#PBS -N GENE_bench
#PBS -j oe
cd $HOME/genetrunk/testsuite
cp $PBS_NODEFILE $HOME/$PBS_JOBID.nodefile
#module add /scHome/angenent/modules/path
#module add intel/cc/11.1.059
#module add intel/mpi/3.2.2.006
#module add intel/mkl/10.2.4.032
#module add fftw/3.2.2
#module add petsc/3.0.0-p12
#module add slepc/3.0.0-p7
mpdboot --rsh=ssh -n 1 -f ~/$PBS_JOBID.nodefile  -v
mpirun --rsh=ssh -machinefile $PBS_NODEFILE -np 12 ./gene_palma

### to submit a parameter scan, comment the previous line and uncomment
### the following (see GENE documentation or ./scanscript --help)
#./scanscript --np 12 --ppn 12 --mps 4
