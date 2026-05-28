Documentation: https://docs.alcf.anl.gov/aurora/getting-started-on-aurora/

Access via:

    ssh [USERNAME]@aurora.alcf.anl.gov

Useful commands:

    myprojectquotas
    sbank
    qstat -u [USERNAME]

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

4. Finally, edit `submit.qsub` to replace the account string with your own and submit it to the queue:
   
   ```bash
   qsub submit.qsub
   ```
   
   This will run the executable in the current directory. The output should be a PBS output file, an output file `sample.out` from the MILC executable, and correlators saved to a file `ks_spectrum_hisq.fpi.2.corrfile_t0.test-out`.
