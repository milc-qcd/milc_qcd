/************************ grsource.c *****************************/
/* MIMD version 6 */

#include "arb_dirac_eig_includes.h"



void grsource() {
register int i,j,k;
register site *s;
    FORALLSITES(i,s){
        for(k=0;k<4;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            s->g_rand.d[k].c[j].real = gaussian_rand_no(&node_prn);
            s->g_rand.d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
    }

}/* grsource */


