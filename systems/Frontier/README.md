Documentation: [Frontier User Guide &mdash; OLCF User Documentation](https://docs.olcf.ornl.gov/systems/frontier_user_guide.html)

Access via:

    ssh [USERNAME]@frontier.olcf.ornl.gov

Useful commands:

    squeue -u [USERNAME]

## Issues

### <span style="color:red">QUDA P2P and ROCM</span>

If you are running QUDA on Frontier (or any AMD MI250X system) with the ROCm 6.x or 7.x modules, the default peer-to-peer halo-exchange path is
broken at the ROCm runtime level and will produce **silent numerical errors** (rocm/6.x, all versions) or **hard crashes** (rocm/7.x, tested up to 7.2.0). The silent numerical errors from rocm/6.x occur specifically in the Dslash halo exchanges that go over P2P. The workaround for both cases is to disable QUDA P2P at runtime (`export QUDA_ENABLE_P2P=0`). The older software stack (`module load PrgEnv-amd amd/5.3.0 rocm/5.3.0`) may also work for you, but it is believed this older software will be removed soon from Frontier.

### QUDA Performance Regression with Large Kernel Argument:

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
