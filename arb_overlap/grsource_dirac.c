/************************ grsource_dirac.c *****************************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

void grsource_dirac(int parity) {
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




void grsource_chiral(int c,field_offset f) {
register int i,j,k;
register site *s;
if (c==PLUS)
{
    FORALLSITES(i,s){
        for(k=0;k<2;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = gaussian_rand_no(&node_prn);
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
        for(k=2;k<4;k++)for(j=0;j<3;j++){
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = 0.;
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = 0.;
	}
		
    }
}
else /* c==MINUS */
{
    FORALLSITES(i,s){
        for(k=0;k<2;k++)for(j=0;j<3;j++){
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = 0.;
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = 0.;
	}
        for(k=2;k<4;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].real = gaussian_rand_no(&node_prn);
            (*((wilson_vector*)F_PT(s,f))).d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
		
    }




}
}/* grsource */


void grsource_chiral_field(chiral_src *f) {
register int i,j,k;
register site *s;
if (f->chiral==PLUS)
{
    FORALLSITES(i,s){
        for(k=0;k<2;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            (f->src[i]).d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            (f->src[i]).d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            (f->src[i]).d[k].c[j].real = gaussian_rand_no(&node_prn);
            (f->src[i]).d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
        for(k=2;k<4;k++)for(j=0;j<3;j++){
            (f->src[i]).d[k].c[j].real = 0.;
            (f->src[i]).d[k].c[j].imag = 0.;
	}
		
    }
}
else if (f->chiral == MINUS)/* c==MINUS */
{
    FORALLSITES(i,s){
        for(k=0;k<2;k++)for(j=0;j<3;j++){
            (f->src[i]).d[k].c[j].real = 0.;
            (f->src[i]).d[k].c[j].imag = 0.;
	}
        for(k=2;k<4;k++)for(j=0;j<3;j++){
#ifdef SITERAND
            (f->src[i]).d[k].c[j].real = gaussian_rand_no(&(s->site_prn));
            (f->src[i]).d[k].c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
            (f->src[i]).d[k].c[j].real = gaussian_rand_no(&node_prn);
            (f->src[i]).d[k].c[j].imag = gaussian_rand_no(&node_prn);
#endif
        }
		
    }




} else {
    printf("bad call to grsource_chiral_field %i\n",f->chiral);
    fflush(stdout);
    terminate(1);

}
}/* grsource */


