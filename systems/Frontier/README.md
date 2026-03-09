Documentation: [Frontier User Guide &mdash; OLCF User Documentation](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)

Access via:

    ssh [USERNAME]@frontier.olcf.ornl.gov

Useful commands:

    squeue -u [USERNAME]

## Issues

##### libfabric:

There is a known issue with the default version of libfabric. See [Frontier User Guide: olcfdev-1811-libfabric-1-20-1-cpu-buffer-performance-regression](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1811-libfabric-1-20-1-cpu-buffer-performance-regression) for details and workarounds.

##### LLVM 17:

There is a bug in LLVM 17, which causes the MILC compilation to fail. This is the version of LLVM that is loaded by default since the Frontier software stack update in July 2024. The failure is triggered during the compilation of `milc_qcd/generic/momentum_twist.c`. The bug report on this issue can be found at https://github.com/llvm/llvm-project/issues/97949 however, they note that it will not be fixed since this version of LLVM is no longer maintained. A workaround is to use the older software stack which loads LLVM 15:

```
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
```

##### ROCm 6:

**As of July 16, 2025, this remains an issue with all versions of ROCm 6+ available on Frontier.**

There is a serious bug in the ROCm 6+ tool chains that can cause incorrect QUDA results without warning when using P2P. 

```
ROCm 5.3.0: Works
ROCm 6.0.0: Works, but only if QUDA P2P disabled
ROCm 6.2.0: Works, but only if QUDA P2P disabled
ROCm 6.2.4: Works, but only if QUDA P2P disabled
ROCm 6.3.1: Works, but only if QUDA P2P disabled
ROCm 6.4.1: Works, but only if QUDA P2P disabled
```

The workaround is to disable QUDA P2P (`export QUDA_ENABLE_P2P=0`) or to use the older software stack:

```
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
```

##### QUDA Performance Regression with Large Kernel Argument:

We have seen a significant QUDA performance regression in certain cases (e.g. eigensolve) related to large kernel arguments. See [QUDA Issue 1568](https://github.com/lattice/quda/issues/1568) and [QUDA PR 1569](https://github.com/lattice/quda/pull/1569) for details.

The performance regression can be avoided by compiling QUDA with:

```
-DQUDA_MAX_KERNEL_ARG_SIZE=0
```

## Building and Running the Sample Code

1. Start by copying the scripts in this directory to a new directory on the system

2. Run the QUDA build script:
   
   ```bash
   bash compile_quda.sh
   ```
   
   This will download the QUDA code into a `quda` directory and compile the QUDA library to a `build` directory.

3. Run the MILC build script:
   
   ```bash
   bash compile_ks_spectrum_hisq.sh
   ```
   
   This will download the MILC code into a `milc_qcd` directory and compile the `ks_spectrum_hisq` executable therein.

4. Finally, edit `submit.sbatch` to replace the SBATCH account string with your own and submit it to the queue:
   
   ```bash
   sbatch submit.sbatch
   ```
   
   This will run the executable in the current directory. The output should be a slurm output file, an output file `sample.out` from the MILC executable, and correlators saved to a file `ks_spectrum_hisq.fpi.2.corrfile_t0.test-out`.
