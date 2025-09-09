#!/bin/bash
#SBATCH -A XXXXXXX 		# Change to your account
#SBATCH -J JobName
#SBATCH -q debug		# Switch to 'regular' for runs > 30m
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4

# Setup environment
module reset
module list

# Path to QIO/QMP libraries: Change path as needed
export LD_LIBRARY_PATH=`pwd`/build/usqcd/lib:$LD_LIBRARY_PATH

export SRUN_CPUS_PER_TASK=8
export OMP_NUM_THREADS=8
export SLURM_CPU_BIND="cores"
export OMP_PROC_BIND="spread, spread, spread"

# Print environment information
export MPICH_ENV_DISPLAY=1
module list
printenv |grep SLURM
echo $LD_LIBRARY_PATH

# Input and output for MILC: Update as needed
output=sample.out
input=sample.in

# Symlink to executable: Change path to MILC as needed
ln -s `pwd`/milc_qcd/milc_qcd/ks_spectrum/ks_spectrum_hisq ks_spectrum_hisq
executable='./ks_spectrum_hisq'

# Print some information
ls -l $input
date > $output
echo "EXECUTABLE = " $executable >> $output
echo ""
echo "EXECUTABLE = " $executable
echo "ldd $executable"
ldd $executable
ls -l $executable >> $output
echo "Job = " $SLURM_JOBID  >> $output
echo "MPI tasks = " $SLURM_NTASKS >> $output
echo "nthreads = $OMP_NUM_THREADS"  >> $output

# Run
srun -N $SLURM_NNODES -n $SLURM_NTASKS $executable $input $output

# Finish up
echo "All done!"

# Print run info
scontrol show job $SLURM_JOBID
sacct -j $SLURM_JOBID

