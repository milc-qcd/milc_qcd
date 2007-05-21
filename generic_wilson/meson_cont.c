/******** meson_cont.c *************/
/* MIMD version 7 */
/* UMH April 96 */

#include "generic_wilson_includes.h"
#include <string.h>

void meson_cont_field(spin_wilson_vector *src1, spin_wilson_vector *src2,
		      int *gamma_in,int *gamma_out,int n_in,int n_out,
		      complex *prop) 
{

register int i;
register site *s; 

int my_t;
int cf, sf, si;
int i_in,i_out;

complex g1,g2;

spin_wilson_vector localmat;       /* temporary storage for quark */
spin_wilson_vector antiquark;      /* temporary storage for antiquark */



    FORALLSITES(i,s){
	my_t = s->t;

	/*first, dirac multiplication by the source gamma matrices (on left) */

	/*  antiquark = c.c. of quark propagator */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++){
	    CONJG(src1[i].d[si].d[sf].c[cf], antiquark.d[si].d[sf].c[cf]);
         }

	/* left multiply antiquark by source gamma matrices,
	   beginning with gamma_5 for quark -> antiquark */
	mult_swv_by_gamma_l( &antiquark, &localmat, GAMMAFIVE);

	/* right dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_swv_by_gamma_r( &localmat, &antiquark, GAMMAFIVE);


	/* left multiply by the particular source dirac matrices */
	for(i_in=0;i_in<n_in;i_in++){
	    mult_swv_by_gamma_l( &antiquark, &localmat, gamma_in[i_in]);
	    antiquark = localmat;
	}

	/* right dirac multiplication by the sink gamma matrices */
	for(i_out=0;i_out<n_out;i_out++) {
	    mult_swv_by_gamma_r( &antiquark, &localmat, gamma_out[i_out]);
	    antiquark = localmat;
	}

	/* copy into quark */
	/* 2/27/98 Simplify and avoid bug in Origin compiler UMH */

	/* trace over propagators */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++)
	{
	    g1 = antiquark.d[si].d[sf].c[cf];
	    g2 = src2[i].d[si].d[sf].c[cf];

	    prop[my_t].real += (g1.real*g2.real - g1.imag*g2.imag); 
	    prop[my_t].imag += (g1.real*g2.imag + g1.imag*g2.real); 
	}

    }

} /* meson_cont_field */

void meson_cont_site(field_offset src1,field_offset src2,
		     int *gamma_in,int *gamma_out,int n_in,int n_out,
		     complex *prop) 
/* src1 and src2 are type spin_wilson_vector */
{
  spin_wilson_vector *t_src1, *t_src2;
  int i;
  site *s;

  t_src1 = (spin_wilson_vector *)
    malloc(sizeof(spin_wilson_vector)*sites_on_node);
  t_src2 = (spin_wilson_vector *)
    malloc(sizeof(spin_wilson_vector)*sites_on_node);
  if(t_src1 == NULL || t_src2 == NULL){
    printf("meson_cont_site(%d): No room for temporaries\n",this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    t_src1[i] = *((spin_wilson_vector *)(F_PT(s,src1)));
    t_src2[i] = *((spin_wilson_vector *)(F_PT(s,src2)));
  }

  meson_cont_field(t_src1, t_src2, gamma_in, gamma_out, n_in, n_out, prop);

  free(t_src1);
  free(t_src2);
}
