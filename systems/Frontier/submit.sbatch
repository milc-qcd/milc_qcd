#!/bin/bash
#SBATCH -A XXXXXXX 		# Change to your account
#SBATCH -J JobName
#SBATCH -q debug
#SBATCH -t 00:30:00
#SBATCH -N 1
#SBATCH -C nvme
#SBATCH --cpus-per-task=7
#SBATCH --ntasks-per-node=8

module reset
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
module load craype-accel-amd-gfx90a
module load libfabric/1.15.2.0

# Path to QUDA libraries: Change as needed
export LD_LIBRARY_PATH=`pwd`/build/usqcd/lib:$LD_LIBRARY_PATH

# QUDA Tuning Directory: Change location as needed
export QUDA_RESOURCE_PATH=tunecache
mkdir -p tunecache

export QUDA_ENABLE_GDR=1
export QUDA_MILC_HISQ_RECONSTRUCT=13
export QUDA_MILC_HISQ_RECONSTRUCT_SLOPPY=9
export GPUDIRECT=" -gpudirect "

export MPICH_ENV_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1

export OMP_NUM_THREADS=7
export OMP_PROC_BIND=spread
MASK_0="0x00fe000000000000"
MASK_1="0xfe00000000000000"
MASK_2="0x0000000000fe0000"
MASK_3="0x00000000fe000000"
MASK_4="0x00000000000000fe"
MASK_5="0x000000000000fe00"
MASK_6="0x000000fe00000000"
MASK_7="0x0000fe0000000000"
MEMBIND="--mem-bind=map_mem:3,3,1,1,0,0,2,2"
CPU_MASK="--cpu-bind=mask_cpu:${MASK_0},${MASK_1},${MASK_2},${MASK_3},${MASK_4},${MASK_5},${MASK_6},${MASK_7}"

# Print environment information
module list
printenv |grep SLURM
echo $LD_LIBRARY_PATH
printenv |grep QUDA

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
srun -N $SLURM_NNODES -n $SLURM_NTASKS ${MEMBIND} ${CPU_MASK} $executable $input $output

# Finish up
echo "All done"

# Print run info
scontrol show job $SLURM_JOBID
sacct -j $SLURM_JOBID

