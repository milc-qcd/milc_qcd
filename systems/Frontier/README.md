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

There is a serious bug in the ROCm 6+ tool chains that causes incorrect QUDA results without warning. The workaround is to use the older software stack:

```
module load PrgEnv-amd amd/5.3.0 rocm/5.3.0
```

## Building and Running the Sample Code

Start by running the QUDA build script:

```
bash compile_quda.sh
```

Next, run the MILC build script:

```
bash compile_ks_spectrum_hisq.sh
```

Finally, edit `submit.sbatch` to replace the SBATCH account string with your own and submit it to the queue:

```
sbatch submit.sbatch
```
