#!/bin/bash
#SBATCH -A XXXXXXX 		# Change to your account
#SBATCH -J JobName
#SBATCH --partition=gpuA100x4
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node 4
#SBATCH --mem=0

# Setup environment
module reset

# Path to QUDA libraries: Change as needed
export LD_LIBRARY_PATH=`pwd`/build/usqcd/lib:$LD_LIBRARY_PATH

# QUDA Tuning Directory: Change location as needed
export QUDA_RESOURCE_PATH=tunecache
mkdir -p tunecache

export CRAY_ACCEL_TARGET=nvidia80
export SRUN_CPUS_PER_TASK=16
export OMP_NUM_THREADS=16
export SLURM_CPU_BIND="cores"
export OMP_PROC_BIND="spread, spread, spread"

export QUDA_ENABLE_GDR=1
export QUDA_MILC_HISQ_RECONSTRUCT=13
export QUDA_MILC_HISQ_RECONSTRUCT_SLOPPY=9
export MPICH_ENV_DISPLAY=1
export MPICH_RDMA_ENABLED_CUDA=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_NEMESIS_ASYNC_PROGRESS=1
export MPICH_SMP_SINGLE_COPY_MODE=XPMEM

# Print environment information
module list
printenv |grep SLURM
echo $LD_LIBRARY_PATH
printenv |grep QUDA

# Print GPU monitoring information periodically
bash monitor-gpu.sh &

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
srun $executable $input $output

# Finish up
echo "All done!"

# Print run info
scontrol show job $SLURM_JOBID
sacct -j $SLURM_JOBID

