Documentation: https://docs.ncsa.illinois.edu/systems/deltaai/en/latest/

Access via:

    ssh [USERNAME]@dtai-login.delta.ncsa.illinois.edu

Useful commands:

    accounts
    quota
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

# 
