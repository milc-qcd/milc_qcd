/******** cvc_pauli.c *************/
/* MIMD version 6 */
#include "arb_dirac_inv_includes.h"

/*  Anomalous magnetic moment current  */

void cvc_pauli(field_offset src1,field_offset src2) 
{

register int i,x,y,z,t;
register site *s; 

int my_t,k;
int cf, sf, ci, si;
complex g1,g2;
Real prop[5][16];
int j;
Real co;


spin_wilson_vector localmat;       /* temporary storage for quark */

wilson_propagator antiquark;      /* temporary storage for antiquark */
wilson_propagator quark;      /* temporary storage for quark */

/* CVC stuff */
msg_tag *tag[2];
int mu,nu;
/* for nonzero momentum */
Real cx,cy,cz,cxy,cyz,cxz,c111;

void mult_by_gamma_l( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);
void mult_by_gamma_r( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);



for(mu=0;mu<4;mu++){

	for(j=0;j<5;j++)
	for(t=0; t<nt; t++){ /* clear meson propagators */
	prop[j][t]  = 0.0;
    }



	FORALLSITES(i,s){s->cvc[mu].real = 0.0;s->cvc[mu].imag = 0.0;}




for(nu=0;nu<4;nu++)if(mu != nu){



    FORALLSITES(i,s){
	s->cvct[nu].real=0.0;
	s->cvct[nu].imag=0.0;

	/*first, dirac multiplication by the source gamma matrices (on left) */

	/*  antiquark = c.c. of quark propagator */
	for(ci=0;ci<3;ci++){

	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++){
	    CONJG(((wilson_propagator *)F_PT(s,src1))->c[ci].d[si].d[sf].c[cf],
		antiquark.c[ci].d[si].d[sf].c[cf]);
         }

	/* left multiply antiquark by source gamma matrices,
	   beginning with gamma_5 for quark -> antiquark */
	mult_by_gamma_l( &antiquark.c[ci], &localmat, GAMMAFIVE);

	/* right dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_r( &localmat, &antiquark.c[ci], GAMMAFIVE);


	/* left multiply by the particular source dirac matrices */
	    mult_by_gamma_l( &antiquark.c[ci], &localmat, mu);

	/* right dirac multiplication by the sink gamma matrices */
	    mult_by_gamma_r( &localmat, &antiquark.c[ci], mu);
	    mult_by_gamma_r( &antiquark.c[ci], &localmat, nu);
	    antiquark.c[ci]=localmat;
	}

	/* copy into quark */
            quark = *(wilson_propagator *)F_PT(s,src2);

	/* trace over propagators */
	for(ci=0;ci<3;ci++)
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++)
	{
	    g1 = antiquark.c[ci].d[si].d[sf].c[cf];
	    g2 = quark.c[ci].d[si].d[sf].c[cf];

	    s->cvct[nu].real += (g1.real*g2.real - g1.imag*g2.imag); 
	    s->cvct[nu].imag += (g1.real*g2.imag + g1.imag*g2.real); 
	}
    }
/* now gather each nu from OPP_DIR(nu) */


        tag[0]=start_gather_site( F_OFFSET(cvct[nu]), sizeof(complex),
            OPP_DIR(nu), EVENANDODD, gen_pt[0] );

/*  gather each nu from ahead */


        tag[1]=start_gather_site( F_OFFSET(cvct[nu]), sizeof(complex),
            nu, EVENANDODD, gen_pt[1] );

	wait_gather(tag[0]);
	wait_gather(tag[1]);

        FORALLSITES(i,s){
		g1= *(complex *)(gen_pt[0][i]);
		g2= *(complex *)(gen_pt[1][i]);
	    s->cvc[mu].real += (g2.real - g1.real); 
	    s->cvc[mu].imag += (g2.imag - g1.imag); 
        } 
	cleanup_gather(tag[0]);
	cleanup_gather(tag[1]);

    } /* nu */




/* print out time sliced observable */

    FORALLSITES(i,s) {
        my_t = s->t;
                for(j=0;j<3;j++){
                cz=cos(2.0*PI/(Real)nz*(Real)(s->z)*(Real)j);
                cx=cos(2.0*PI/(Real)nx*(Real)(s->x)*(Real)j);
                cy=cos(2.0*PI/(Real)ny*(Real)(s->y)*(Real)j);


                prop[j][my_t] += (s->cvc[mu].real)*
                                (cx+cy+cz)/3.0;
                }
                cxy=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y));
                cxz=cos(2.0*PI/(Real)nz*(Real)(s->x +s->z));
                cyz=cos(2.0*PI/(Real)nz*(Real)(s->y +s->z));
                c111=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y + s->z));
                prop[3][my_t] += (s->cvc[mu].real)*(cxy+cyz+cxz)/3.0;
                prop[4][my_t] += (s->cvc[mu].real)*c111;
        }

    /* sum and print meson propagators */
            for(t=0; t<nt; t++){
                if(this_node == 0)
                    printf("PAULICV%d  %d",mu,t);
                  for(j=0;j<MAX_P;j++){
                  g_floatsum( &prop[j][t] );
                  if(this_node == 0)
                   printf(" %e",(double)prop[j][t]);
                  }
                if(this_node == 0)printf("\n");
            }

} /* mu */


} /* meson_cont */


