/******** w_spectrum_z2.c *************/
/* Screening spectrum:  pion, sigma, rho, a1, but no nucleons */
/* also includes sigma propagator */
/* MIMD version 6 */

/* requires LU */
#ifndef LU
Bomb the compile
#endif

/* screening spectrum for Wilson hadrons */
#include "wi_dyn_includes.h"
#define POINT 1
#define WALL 2

int w_spectrum_z() 
{
register site *s;
register int i;
int iters,spin,color;
Real rsq;
void w_source_z(field_offset src,int color,int spin,int source,
		int x0,int y0,int z0,int t0,Real gamma);
void w_meson_z(field_offset src1,field_offset src2);

    iters=0;
    for(spin=0;spin<4;spin++) for(color=0;color<3;color++) {
	w_source_z(F_OFFSET(chi),color,spin,WALL,0,0,0,0,(Real)0.0);
	/* multiply by L inverse */
	dslash_w_site( F_OFFSET(chi), F_OFFSET(psi), PLUS, EVEN);
	FOREVENSITES(i,s){
	    scalar_mult_add_wvec( &(s->chi), &(s->psi), kappa, &(s->chi));
	}
	/* Multiply by M_adjoint */
	dslash_w_site( F_OFFSET(chi), F_OFFSET(psi), MINUS, ODD);
	dslash_w_site( F_OFFSET(psi), F_OFFSET(psi), MINUS, EVEN);
	FOREVENSITES(i,s){
	    scalar_mult_add_wvec( &(s->chi), &(s->psi),
 		-kappa*kappa, &(s->chi) );
	}
	/* take M_adjoint * M inverse, result in psi */
	   iters += congrad_w( niter,rsqprop,&rsq);
	/* fix up odd sites in psi, which congrad doesn't compute */
	/* In fact, congrad may muck up these sites of psi!! */
	FORODDSITES(i,s){ s->psi = s->chi; }
	/* Multiply by U inverse */
	dslash_w_site( F_OFFSET(psi), F_OFFSET(chi), PLUS, ODD);
	FORODDSITES(i,s){
	    scalar_mult_add_wvec( &(s->psi), &(s->chi), kappa, &(s->psi));
	}
	/* copy result to quark_propagator - KLUDGE- FIX CONGRAD! */
	FORALLSITES(i,s){
	    s->quark_propagator.d[spin].c[color] = s->psi;
	}
    }  /* end loop over source spins and colors */
    w_meson_z(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator));
    return(iters);
} /* spectrum */


void w_source_z(field_offset src,int color,int spin,int source,
		int x0,int y0,int z0,int t0,Real gamma)
{
register int i;
register site *s; 

Real rx,ry,rt,radius2;
register Real theta;
void clear_wvec(wilson_vector *);
/*printf("WSOURCE: source = %d\n",source); */
	
	/* zero src to be safe */
	FORALLSITES(i,s) {
	    clear_wvec((wilson_vector *)F_PT(s,src)); 
	}

	if(source == POINT) {
            /* load 1.0 into source at cooordinates given by source_coord */
	    /* initialize src to be a delta function at point x0,y0,z0,t0 */

	    if(node_number(x0,y0,z0,t0) == mynode()){
		i = node_index(x0,y0,z0,t0);
		((wilson_vector *)F_PT(&(lattice[i]),src))->
		    d[spin].c[color].real = 1.0;
	    }
	}
	else if(source == WALL) {
            /* Gaussian trial source centered on  midpoint of z-slice 0 */
	    /* Modulate with Matsubara phase */

	    FORALLSITES(i,s) {
		if(s->z != 0)continue;	/* only do this if z==0 */
		rx=((Real)(s->x - nx/2));
		ry=((Real)(s->y - ny/2));
		rt=((Real)(s->t - nt/2));
		radius2= rx*rx + ry*ry + rt*rt;
	        theta =  PI*(Real)( s->t)/(Real)nt;
		((wilson_vector *)F_PT(s,src))->d[spin].c[color].real = 
		    (Real)(exp((double)(- radius2*gamma)) *
		    cos( (double)theta ) );
	    }
	}

} /* w_source_z */


void w_meson_z(field_offset src1,field_offset src2) 
/* src1 and src2 are type wilson_matrix */
{
/* pion and rho with pointlike sink, p=0 */

register int i,z;
register site *s; 

int my_z;
int cf, sf, ci, si;
complex g1,g2;
Real *prop,*prop0;

wilson_matrix localmat;       /* temporary storage for quark */
wilson_matrix localmat2;       /* temporary storage for quark */

    prop = (Real *)malloc(nz*sizeof(Real));
    /* for z slice sums of propagators */
    prop0 = (Real *)malloc(nz*sizeof(Real));
    /* for z slice sums of propagators, with extra gamma_z at sink */

    /* measure the pion propagator--gamma-5 gamma-5 correlator and
       gamma-5 gamma-5*gamma-0 correlator */
    /* this is generic code which is commented to show you how to change it */
    for(z=0; z<nz; z++){ /* clear meson propgators */
	prop[z] = prop0[z] = 0.0;
    }

    /*first, dirac multiplication by the source gamma matrices */
    FORALLSITES(i,s) {
	my_z = s->z;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2,
	    GAMMAFIVE);
	/*first, dirac multiplication by gamma-5 (antiquark will need this
	  along with complex conjugation */
	mult_by_gamma_left( &localmat2, &localmat, GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	/** Commented out for this example, because gamma_5^2 = 1 **
	mult_by_gamma_right( &localmat, &localmat2, GAMMAFIVE);
	mult_by_gamma_right( &localmat2, &localmat, GAMMAFIVE);
	**/
	/* Multiply by gamma_z for gamma_5 - gamma_5*gamma_z correlator */
	mult_by_gamma_right( &localmat, &localmat2, ZUP);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {
                  g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
                  g1 = localmat.d[si].c[ci].d[sf].c[cf];
                  prop[my_z] += g1.real*g2.real + g1.imag*g2.imag;
                  g1 = localmat2.d[si].c[ci].d[sf].c[cf];
                  prop0[my_z] += g1.real*g2.real + g1.imag*g2.imag;
	}

    }
    /* sum and print meson propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	g_floatsum( prop0+z );
	if(mynode()==0)printf("PSEUDO2_Z  %d  %e   %e\n",z,
	    (double)prop[z], (double)prop0[z] );
    } 


    /* SIGMA propagator-- 1 - 1 correlator */
    for(z=0; z<nz; z++){ /* clear meson propgators */
	prop[z]=0.0;
    }
    /*first, dirac multiplication by the source gamma matrices */
    FORALLSITES(i,s) {
	my_z = s->z;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	/* dirac multiplication by gamma-5 (antiquark will need this
	   along with complex conjugation */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2,
	    GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	/* dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_right( &localmat2, &localmat, GAMMAFIVE);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {

              g1 = localmat.d[si].c[ci].d[sf].c[cf];
              g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
              prop[my_z] += g1.real*g2.real + g1.imag*g2.imag;
	}

    }
    /* print meson propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	if(mynode()==0)printf("SIGMA  %d  %e\n",z,(double)prop[z]);
    } 


    /* RHO propagator-- gamma3gamma-1 - gamma3gamma-1 correlator */
    /* this is generic code which is commented to show you how to change it */
    for(z=0; z<nz; z++){ /* clear meson propgators */
	prop[z]=0.0;
    }
    /*first, dirac multiplication by the source gamma matrices */
    FORALLSITES(i,s) {
	my_z = s->z;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2, XUP);
	mult_by_gamma_left( &localmat2, &localmat, ZUP);
	/* dirac multiplication by gamma-5 (antiquark will need this
	   along with complex conjugation */
	mult_by_gamma_left( &localmat, &localmat2, GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	mult_by_gamma_right( &localmat2, &localmat, ZUP);
	mult_by_gamma_right( &localmat, &localmat2, XUP);
	/* dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_right( &localmat2, &localmat, GAMMAFIVE);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {

              g1 = localmat.d[si].c[ci].d[sf].c[cf];
              g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
	      /*minus sign since propagator seems to be negative without it */
	      /* This is because we really should be using adjoint of
		hadron operator at one end */
              prop[my_z] -= g1.real*g2.real + g1.imag*g2.imag;
	}

    }
    /* print meson propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	if(mynode()==0)printf("RHO31312  %d  %e\n",z,(double)prop[z]);
    } 

    /* RHO propagator--gamma_1 gamma_1 correlator */
    /* this is generic code which is commented to show you how to change it */
    for(z=0; z<nz; z++){ /* clear meson propgators */
	prop[z]=0.0;
    }
    FORALLSITES(i,s) {
	my_z = s->z;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2, XUP);
	/*first, dirac multiplication by gamma-5 (antiquark will need this
	  along with complex conjugation */
	mult_by_gamma_left( &localmat2, &localmat, GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	mult_by_gamma_right( &localmat, &localmat2, XUP);
	/* dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_right( &localmat2, &localmat, GAMMAFIVE);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {

              g1 = localmat.d[si].c[ci].d[sf].c[cf];
              g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
              prop[my_z] +=  (g1.real*g2.real + g1.imag*g2.imag);
	}
    }
    /* print meson propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	if(mynode()==0)printf("RHO112  %d  %e\n",z,(double)(*(prop+z)));
    }

    /* A1 propagator--gamma_1-gamma_5 gamma_1-gamma_5 correlator */
    /* this is generic code which is commented to show you how to change it */
    for(z=0; z<nz; z++){ /* clear meson propgators */
	prop[z]=0.0;
    }
    FORALLSITES(i,s) {
	my_z = s->z;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2, XUP);
	/* then multiplication by gamma-5, then another gamma-5 for
	  anitquark, so actually do nothing. */

	/* dirac multiplication by the sink gamma matrices (on right) */
	mult_by_gamma_right( &localmat2, &localmat, XUP);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {

              g1 = localmat.d[si].c[ci].d[sf].c[cf];
              g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
              prop[my_z] +=  (g1.real*g2.real + g1.imag*g2.imag);
	}
    }
    /* print meson propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	if(mynode()==0)printf("A1_112  %d  %e\n",z,(double)(*(prop+z)));

    } 

free(prop); free(prop0);
} /* meson */
