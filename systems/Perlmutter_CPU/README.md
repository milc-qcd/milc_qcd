Documentation: [Architecture - NERSC Documentation](https://docs.nersc.gov/systems/perlmutter/architecture/)

Access via:

    ssh [USERNAME]@perlmutter.nersc.gov

Useful commands:

    squeue -u [USERNAME]

## Building and Running the Sample Code

Start by running the QIO/QMP build script:

```
bash compile_qioqmp.sh
```

Next, run the MILC build script:

```
bash compile_ks_spectrum_hisq.sh
```

Finally, edit `submit.sbatch` to replace the SBATCH account string with your own and submit it to the queue:

```
sbatch submit.sbatch
```

## Performance Tuning Notes

The default MPICH environment on Perlmutter seems to be well-tuned for our purposes.
