/******** w_spectrum.c *************/
/* MIMD version 6 */

/* requires LU */
#ifndef LU
Bomb the compile
#endif

/* spectrum for Wilson hadrons */
#include "wi_dyn_includes.h"
#define POINT 1
#define WALL 2

int w_spectrum() 
{
register site *s;
register int i;
int iters,spin,color;
Real rsq;
/* NOTE: these three local routines should probably be replaced
   with the generic routines in generic_wilson - since they are
   more recent - CD 11/09/97 */
void w_source2(field_offset src,int color,int spin,int source,
	 int x0,int y0,int z0,int t0,Real gamma);
void w_meson2(field_offset src1,field_offset src2);
void w_baryon2(field_offset src1,field_offset src2,field_offset src3);

    iters=0;
    for(spin=0;spin<4;spin++) for(color=0;color<3;color++) {
	w_source2(F_OFFSET(chi),color,spin,WALL,0,0,0,0,(Real)0.0);
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
/**
dslash_w_site( F_OFFSET(psi), F_OFFSET(mp), PLUS, EVENANDODD);
FORALLSITES(i,s)scalar_mult_add_wvec( &(s->psi), &(s->mp), -kappa, &(s->mp) );
printf("DUMP!!\n");
FORALLSITES(i,s){
    printf("i = %d,  coords %d %d %d %d\n",i,s->x,s->y,s->z,s->t);
    dump_wvec( &(s->mp) );
}
**/
	/* copy result to quark_propagator - KLUDGE- FIX CONGRAD! */
	FORALLSITES(i,s){
	    s->quark_propagator.d[spin].c[color] = s->psi;
	}
    }  /* end loop over source spins and colors */
    w_meson2(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator));
    w_baryon2(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator),
        F_OFFSET(quark_propagator));
    return(iters);
} /* spectrum */


void w_source2(field_offset src,int color,int spin,int source,
	 int x0,int y0,int z0,int t0,Real gamma)
{
register int i;
register site *s; 

Real rx,ry,rz,radius2;
	void clear_wvec();
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
            /* Gaussian trial source centered on  midpoint of timeslice t0 */

	    FORALLSITES(i,s) {
		if(s->t != 0)continue;	/* only do this if t==0 */
		rx=((Real)(s->x - nx/2));
		ry=((Real)(s->y - ny/2));
		rz=((Real)(s->z - nz/2));
		radius2= rx*rx + ry*ry + rz*rz;
		    ((wilson_vector *)F_PT(s,src))->d[spin].c[color].real = 
		(Real)exp((double)(- radius2*gamma));
	    }
	}

} /* w_source2 */


void w_meson2(field_offset src1,field_offset src2) 
/* src1 and src2 are type wilson_matrix */
{
/* pion and rho with pointlike sink, p=0 */

register int i,t;
register site *s; 

int my_t;
int cf, sf, ci, si;
complex g1,g2;
Real *prop,*prop0;

wilson_matrix localmat;       /* temporary storage for quark */
wilson_matrix localmat2;       /* temporary storage for quark */

    prop = (Real *)malloc(nt*sizeof(Real));
    /* for time slice sums of propagators */
    prop0 = (Real *)malloc(nt*sizeof(Real));
    /* for time slice sums of propagators, with extra gamma_0 at sink */

    /* measure the pion propagator--gamma-5 gamma-5 correlator and
       gamma-5 gamma-5*gamma-0 correlator */
    /* this is generic code which is commented to show you how to change it */
    for(t=0; t<nt; t++){ /* clear meson propgators */
	prop[t] = prop0[t] = 0.0;
    }

    /*first, dirac multiplication by the source gamma matrices */
    FORALLSITES(i,s) {
	my_t = s->t;
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
	/* Multiply by gamma_0 for gamma_5 - gamma_5*gamma_0 correlator */
	mult_by_gamma_right( &localmat, &localmat2, TUP);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {
                  g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
                  g1 = localmat.d[si].c[ci].d[sf].c[cf];
                  prop[my_t] += g1.real*g2.real + g1.imag*g2.imag;
                  g1 = localmat2.d[si].c[ci].d[sf].c[cf];
                  prop0[my_t] += g1.real*g2.real + g1.imag*g2.imag;
                  /**g3 = localmat2.d[si].c[ci].d[sf].c[cf];
                  prop0[my_t] += g3.real*g2.real + g3.imag*g2.imag;**/
	}

    }
    /* sum and print meson propagators */
    for(t=0; t<nt; t++){
	g_floatsum( prop+t );
	g_floatsum( prop0+t );
	if(mynode()==0)printf("PSEUDO2  %d  %e   %e\n",t,
	    (double)prop[t], (double)prop0[t] );
    } 


    /* measure the rho propagator-- gamma0gamma-3 - gamma0gamma-3 correlator */

    /* this is generic code which is commented to show you how to change it */
    for(t=0; t<nt; t++){ /* clear meson propgators */
	prop[t]=0.0;
    }

    /*first, dirac multiplication by the source gamma matrices */
    FORALLSITES(i,s) {
	my_t = s->t;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2, ZUP);
	mult_by_gamma_left( &localmat2, &localmat, TUP);
	/* dirac multiplication by gamma-5 (antiquark will need this
	   along with complex conjugation */
	mult_by_gamma_left( &localmat, &localmat2, GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	mult_by_gamma_right( &localmat2, &localmat, TUP);
	mult_by_gamma_right( &localmat, &localmat2, ZUP);
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
              prop[my_t] -= g1.real*g2.real + g1.imag*g2.imag;
	}

    }
    /* print meson propagators */
    for(t=0; t<nt; t++){
	g_floatsum( prop+t );
	if(mynode()==0)printf("RHO03032  %d  %e\n",t,(double)prop[t]);
    } 

    /* measure the rho propagator--gamma_3 gamma_3 correlator */

    /* this is generic code which is commented to show you how to change it */
    for(t=0; t<nt; t++){ /* clear meson propgators */
	prop[t]=0.0;
    }

    FORALLSITES(i,s) {
	my_t = s->t;
	/*first, dirac multiplication by the source gamma matrices (on left) */
	mult_by_gamma_left( ((wilson_matrix *)F_PT(s,src1)), &localmat2, ZUP);
	/*first, dirac multiplication by gamma-5 (antiquark will need this
	  along with complex conjugation */
	mult_by_gamma_left( &localmat2, &localmat, GAMMAFIVE);

	/* dirac multiplication by the sink gamma matrices (on right) */
	mult_by_gamma_right( &localmat, &localmat2, ZUP);
	/* dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_right( &localmat2, &localmat, GAMMAFIVE);
	/* trace over propagators */
	for(si=0;si<4;si++)
	  for(ci=0;ci<3;ci++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {

              g1 = localmat.d[si].c[ci].d[sf].c[cf];
              g2 = ((wilson_matrix *)F_PT(s,src2))->d[si].c[ci].d[sf].c[cf];
              prop[my_t] +=  (g1.real*g2.real + g1.imag*g2.imag);
	}
    }
    /* print meson propagators */
    for(t=0; t<nt; t++){
	g_floatsum( prop+t );
	if(mynode()==0)printf("RHO332  %d  %e\n",t,(double)(*(prop+t)));
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
#define DELTA_N 2
void w_baryon2(field_offset src1,field_offset src2,field_offset src3) 
/* src1's  are type wilson_matrix */
{

register int isite,t;
register site *s; 

int my_t;
int i,j,k;

int ispinp1, ispinp2;
int si_1,si_2,si_3,sf_1,sf_2,sf_3,ci_1,ci_2,ci_3,cf_1,cf_2,cf_3;
int  chi_i, chi_f, eps_f, eps_i;
Real factor, iso_fac;

int  baryons;

int chi_in[4][4][4],chi_out[4][4][4];
int eps[3][3][3];


/* static char *bar_kind[2] = {"PROTON","DELTA"}; */
static char *bar_kind[3] = {"PROTON","DELTA","DELTA_N"};

Real *prop;
complex psi_sq, diquark, diquark_temp;

prop = (Real *)malloc(nt*sizeof(Real));

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



/* for(baryons=PROTON;baryons<=DELTA;baryons++) */
for(baryons=PROTON;baryons<=DELTA_N;baryons++)
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

	iso_fac = 1.0;
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

	iso_fac = 1.0;
	break;

	case DELTA_N:
	chi_in[0][0][ispinp1] = 1;
	chi_in[2][2][ispinp1] = 1;

	chi_out[0][0][ispinp2] = 1;
	chi_out[2][2][ispinp2] = 1;

	iso_fac = 2.0;
	break;
	}

  /* set Baryon propagator to zero */
    for(t=0; t<nt; t++)
	{
	*(prop+t)=0.0;
	}

  /* begin sum over source quark spin and color (labelled by 'f' suffix) */
          FORALLSITES(isite,s)
	{
	my_t = s->t;

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
		factor = factor * iso_fac;
		CMULREAL(diquark_temp,factor,diquark_temp);
		CSUB(diquark,diquark_temp,diquark);



/*  finally tie the diquark to quark 1 to form the baryon */
		CMUL(
((wilson_matrix *)F_PT(s,src1))->d[sf_1].c[cf_1].d[si_1].c[ci_1],
diquark,psi_sq);
/* and accumulate into the baryon propagator */
	*(prop  + my_t) += psi_sq.real;

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
    for(t=0; t<nt; t++){
	g_floatsum( prop+t );
	if(mynode()==0)printf("%s  %d  %e\n",
	    bar_kind[baryons],t,(double)(*(prop+t)));
    }
}      /* end loop on baryon kind */

    free(prop);
} /* spectrum */
#undef Nc
#undef Ns
#undef PROTON
#undef DELTA
#undef DELTA_N
