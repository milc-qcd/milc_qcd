/****** imp_gauge_action_gpu.c  -- ******************/
/* MIMD version 7 */
/* gauge action stuff for improved action
* T.D. and A.H. general gauge action updating code
* D.T. modified  5/97
* D.T. modified 12/97, optimized gauge_force a little
* D.T. modified 3/99, gauge action in include file
* E.W. modified 7/22, split off gauge action, added GPU implementation */

/**#define GATIME**/ /* For timing gauge action calculation */
#include "generic_includes.h"	/* definitions files and prototypes */
#include "../include/openmp_defs.h"
#include "../include/generic_quda.h"

double imp_gauge_action_gpu() {
    register int i;
    int rep;
    register site *s;
    complex trace;
    double g_action;
    double action,act2,total_action;
    su3_matrix *tempmat1;
    su3_matrix *links;
    int length;

    /* these are for loop_table  */
    int ln,iloop;

    /* get loop variables from functions */
    const int max_length = get_max_length();
    const int nloop = get_nloop();
    const int nreps = get_nreps();
    const int *loop_length = get_loop_length();
    const int *loop_num = get_loop_num();
    int ***loop_table = get_loop_table();
    Real **loop_coeff_milc = get_loop_coeff();

    if (nreps != 1){
      printf("imp_gauge_action_gpu: Does not support nreps != 1, disable gauge action offload\n");
      terminate(1);
    }

    // Count total number of loops
    int num_paths = 0;
    for (iloop = 0; iloop < nloop; iloop++)
        for (ln = 0; ln < loop_num[iloop]; ln++)
            num_paths++;

#ifdef GATIME
    int nlinks = 0;
    for (iloop = 0; iloop < nloop; iloop++) {
        nlinks += loop_num[iloop] * loop_length[iloop];
    }
    int nflop = 198 * nlinks + 8 * num_paths; /* For any action */
    double dtime = -dclock();
#endif

    // Storage for traces
    double *traces = (double*)malloc(2 * num_paths * sizeof(double));

    // Storage for input paths
    int **input_path_buf = (int**)malloc(num_paths * sizeof(int*));
    for (i = 0; i < num_paths; i++)
        input_path_buf[i] = (int*)malloc(max_length * sizeof(int));

    // Storage for path lengths
    int *path_length = (int*)malloc(num_paths * sizeof(int));

    // Storage for loop coefficients
    double *loop_coeff = (double*)malloc(num_paths * sizeof(double));

    // Overall scaling factor
    double factor = 1;
    
    // Populate arrays
    num_paths = 0;
    for (iloop = 0; iloop < nloop; iloop++) {
        length = loop_length[iloop];
        for (ln = 0; ln < loop_num[iloop]; ln++) {
            path_length[num_paths] = length; // path length
            loop_coeff[num_paths] = 1.0; // due to the "3. - [...]" convention below, we'll wait to scale then
            for (i = 0; i < length; i++)
                input_path_buf[num_paths][i] = loop_table[iloop][ln][i];
            num_paths++;
        }
    }

    site *st;

    initialize_quda();

    QudaMILCSiteArg_t arg = newQudaMILCSiteArg();

    qudaGaugeLoopTracePhased(MILC_PRECISION, traces, input_path_buf, path_length, loop_coeff, num_paths,
                             max_length, factor, &arg, phases_in);

    g_action = 0.0;

    // traces has been populated so we now accumulate the actions
    num_paths = 0;
    for (iloop = 0; iloop < nloop; iloop++) {
        for (ln = 0; ln < loop_num[iloop]; ln++) {
            action = 3.0 * volume - traces[2 * num_paths]; // extract real part
            total_action = loop_coeff_milc[iloop][0] * action;
            g_action += total_action;
            num_paths++;
        }
    }
    free(loop_coeff);
    free(path_length);
    for (i = 0; i < num_paths; i++)
        free(input_path_buf[i]);
    free(input_path_buf);
    free(traces);

#ifdef GATIME
    dtime+=dclock();
    node0_printf("GATIME:   time = %e (GaugeAction_Symanzik1_QUDA) mflops = %e\n",dtime,
                  nflop*(double)volume/(1e6*dtime*numnodes()) );
#endif

    return( g_action );
} /* imp_gauge_action_gpu */


