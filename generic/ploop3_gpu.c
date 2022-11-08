/****************** ploop3_gpu.c ************************************/
/* MIMD version 7 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* DT 12/97 use local matrix "tmat" instead of "tempmat2" */
/* EW 7/22 Added a QUDA GPU version */

#include "generic_includes.h"
#include "../include/generic_quda.h"

complex ploop_gpu() {
    complex p_loop;

    initialize_quda();

    QudaMILCSiteArg_t arg = newQudaMILCSiteArg();

    double p_loop_[2];
    qudaPolyakovLoopPhased(MILC_PRECISION, p_loop_, 3, &arg, phases_in);
    p_loop.real = p_loop_[0];
    p_loop.imag = p_loop_[1];

    return p_loop;
} /* ploop_gpu */


