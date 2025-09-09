Documentation: [Delta User Documentation &mdash; UIUC NCSA Delta User Guide](https://docs.ncsa.illinois.edu/systems/delta/en/latest/index.html)

Access via:

    ssh [USERNAME]@login.delta.ncsa.illinois.edu

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
