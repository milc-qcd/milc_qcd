/******** w_spectrum_z.c *************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE. Consider using w_spectrum_z2.c */

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
    w_baryon_z(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator),
        F_OFFSET(quark_propagator));
    return(iters);
} /* spectrum */


void w_source_z(field_offset src,int color,int spin,int source,
		int x0,int y0,int z0,int t0,Real gamma)
{
register int i;
register site *s; 

short my_x,my_y,my_z;
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

register int i,x,y,z,t;
register site *s; 

int my_z,k;
int cf, sf, ci, si;
int cf1, cf2, sf1, sf2;
int ci1, ci2, si1, si2;
complex g1,g2,g3;
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


    /* measure the rho propagator-- gamma3gamma-1 - gamma3gamma-1 correlator */

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

    /* measure the rho propagator--gamma_1 gamma_1 correlator */

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

free(prop); free(prop0);
} /* meson */


/* Baryon propagator:
 propagator            point to plane baryon correlation (write)

 We are using a chiral basis for Dirac matrices (gamma-5 diagonal)
 In that basis relativistic wave functions are much easier to implement.
 For the proton we have 
         |p,ispinp>=(u C gamma_5 d)u_ispinp     (1)
                   =(u_1 d_2 - u_2 d_1 + u_3 d_4 -u_4 d_3)u_ispinp
 for a proton of Dirac index (=helicity in a chiral basis).
  The wave function is encoded in the terms chi_in(i,j,k) and chi_out(i,j,k)
 where the first index labels the spinor of the d quark and the 
 last two indices are for the two u quarks.
 The Delta is similar, only all plus signs

  We use the unsymmetrized w. f. and
  include direct and exchange terms in the propagator explicitly
         ___
         \
 b(t) =   >  < b(t , 0) b(t + t, x) >
         /        1        1
         ---
          x
*/

#define Nc 3
#define Ns 4

#define PROTON 0
#define DELTA 1
void w_baryon_z(field_offset src1,field_offset src2,field_offset src3) 
/* src1's  are type wilson_matrix */
{

register int isite,z;
register site *s; 

int my_z;
int i,j,k;

int ispinp1, ispinp2;
int si_1,si_2,si_3,sf_1,sf_2,sf_3,ci_1,ci_2,ci_3,cf_1,cf_2,cf_3;
int  chi_i, chi_f, eps_f, eps_i;
Real factor;

int  baryons;

int chi_in[4][4][4],chi_out[4][4][4];
int eps[3][3][3];


static char *bar_kind[2] = {"PROTON","DELTA"};

Real *prop;
complex psi_sq, diquark, diquark_temp;

prop = (Real *)malloc(nz*sizeof(Real));

/* ispinp1 and ispinp2 are the spins of the source and sink */
ispinp1 = 0; ispinp2 = 0;


  /* initialize epsilon factors */
	for( i=0; i< 3; i++ )
	for( j=0; j< 3; j++ )
	for( k=0; k< 3; k++ )
	{
	eps[i][j][k]=0;
	}
	eps[0][1][2]=1;eps[1][2][0]=1;eps[2][0][1]=1;
	eps[0][2][1]= -1;eps[1][0][2]= -1;eps[2][1][0]= -1;



for(baryons=PROTON;baryons<=DELTA;baryons++)
{
  /* initialize chi factors */
	for( i=0; i< 4; i++ )
	for( j=0; j< 4; j++ )
	for( k=0; k< 4; k++ )
	{
	chi_in[i][j][k] = 0;
	chi_out[i][j][k] = 0;
	}
switch(baryons)
	{
	case PROTON:
	chi_in[1][0][ispinp1] = 1;
	chi_in[0][1][ispinp1] =  -1;
	chi_in[3][2][ispinp1] = 1;
	chi_in[2][3][ispinp1] =  -1; 

	chi_out[1][0][ispinp2] = 1;
	chi_out[0][1][ispinp2] =  -1;
	chi_out[3][2][ispinp2] = 1;
	chi_out[2][3][ispinp2] =  -1;
	break;

	case DELTA:
	chi_in[1][0][ispinp1] = 1;
	chi_in[0][1][ispinp1] = 1;
	chi_in[3][2][ispinp1] = 1;
	chi_in[2][3][ispinp1] = 1; 

	chi_out[1][0][ispinp2] = 1;
	chi_out[0][1][ispinp2] = 1;
	chi_out[3][2][ispinp2] = 1;
	chi_out[2][3][ispinp2] = 1;
	break;
	}

  /* set Baryon propagator to zero */
    for(z=0; z<nz; z++)
	{
	*(prop+z)=0.0;
	}

  /* begin sum over source quark spin and color (labelled by 'f' suffix) */
          FORALLSITES(isite,s)
	{
	my_z = s->z;

	sf_3 =  ispinp2;
	for(cf_3=0; cf_3<Nc; cf_3++)
	{




 /* Sum over source diquark components */
	for(sf_1=0; sf_1<Ns; sf_1++)
	for(sf_2=0; sf_2<Ns; sf_2++) 
	{
	chi_f = chi_out[sf_1][ sf_2][ sf_3];
	if(chi_f != 0)
	{
	for(cf_1=0; cf_1<Nc; cf_1++)
	for(cf_2=0; cf_2<Nc; cf_2++)
	{
	eps_f = eps[cf_1][cf_2][cf_3];
	if(eps_f != 0)
	{


/*  combine the sink colors and spins of these three propagators
 to compose the outgoing baryon */

 /* Sum over sink quark color and spin (labelled by 'i' suffix)*/
		for(si_3=0; si_3<Ns; si_3++)
		for(ci_3=0; ci_3<Nc; ci_3++)
		{
    /* begin sum over sink diquark components */
		for(si_1=0; si_1<Ns; si_1++)
		for(si_2=0; si_2<Ns; si_2++)
		{
		chi_i = chi_in[si_1][ si_2][ si_3];
		if(chi_i != 0)
		{
		for(ci_1=0; ci_1<Nc; ci_1++)
		for(ci_2=0; ci_2<Nc; ci_2++)
		{
		eps_i = eps[ci_1][ci_2][ci_3];
		if(eps_i != 0)
		{
			factor = (Real)(eps_f*eps_i*chi_i*chi_f);

/*  build the 2-3 diquark direct term */
		CMUL(
((wilson_matrix *)F_PT(s,src2))->d[sf_2].c[cf_2].d[si_2].c[ci_2],
((wilson_matrix *)F_PT(s,src3))->d[sf_3].c[cf_3].d[si_3].c[ci_3],
diquark_temp);
		CMULREAL(diquark_temp,factor,diquark);


/*  next, build the 2-3 diquark exchange term and subtract it from the 
diquark */
		CMUL(
((wilson_matrix *)F_PT(s,src2))->d[sf_2].c[cf_2].d[si_3].c[ci_3],
((wilson_matrix *)F_PT(s,src3))->d[sf_3].c[cf_3].d[si_2].c[ci_2],
diquark_temp);
		CMULREAL(diquark_temp,factor,diquark_temp);
		CSUB(diquark,diquark_temp,diquark);



/*  finally tie the diquark to quark 1 to form the baryon */
		CMUL(
((wilson_matrix *)F_PT(s,src1))->d[sf_1].c[cf_1].d[si_1].c[ci_1],
diquark,psi_sq);
/* and accumulate into the baryon propagator */
	*(prop  + my_z) += psi_sq.real;

		}/* end if eps_i */
		} /* Close ci_1, ci_2 */
		} /* end if chi_i */
		} /* Close si_1, si_2 */
		} /* Close  si_3,ci_3 */
/*  this ends all loops over the sink */
          
	} /* end if eps_f */
	}  /* end do cf_1, cf_2 */
	}  /*end if chi_f */
	}  /*end do sf_1,sf_2 */
	}  /*end do cf_3   (and sf_3,if used) */
	}  /*end loop for sites */
  
	/* print baryon propagators */
    for(z=0; z<nz; z++){
	g_floatsum( prop+z );
	if(mynode()==0)printf("%s  %d  %e\n",
	    bar_kind[baryons],z,(double)(*(prop+z)));
    }
}      /* end loop on baryon kind */

    free(prop);
} /* spectrum */
#undef Nc
#undef Ns
#undef PROTON
#undef DELTA
