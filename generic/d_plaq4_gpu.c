/************************** d_plaq4_gpu.c *******************************/
/* MIMD version 7 */

/* Double precision version of "plaquette4.c" including optional
   Schroedinger functional - UMH - 1/27/00
   Extended with QUDA support - ESW - 6/24/22 */

/* Measure the average plaquette of the space-space and
   space-time plaquettes */

#include "generic_includes.h" /* definitions files and prototypes */
#include "../include/generic_quda.h"

void d_plaquette_gpu(double *ss_plaq, double *st_plaq) {

    initialize_quda();

    QudaMILCSiteArg_t arg = newQudaMILCSiteArg();
    double plaq[3];
    qudaPlaquettePhased(MILC_PRECISION, plaq, &arg, phases_in);

    /* plaquettes */
    *ss_plaq = 3.0 * plaq[1];
    *st_plaq = 3.0 * plaq[2];

}

