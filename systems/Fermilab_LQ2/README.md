Documentation: [https://computing.fnal.gov/lqcd/](https://computing.fnal.gov/lqcd/)

Access via:

    ssh [USERNAME]@perlmutter.nersc.gov

Useful commands:

    lfs quota -hu [USERNAME] /lustre1 # To check your lustre file system quota
    squeue -u [USERNAME]

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
