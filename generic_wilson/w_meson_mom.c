/************************ w_meson_mom.c *****************************/
/* MIMD version 7  ***/

/* The optimizations in this version assume that the gamma_corr table
   is sorted so all entries with the same gout are contiguous */

/* Modifications
   Original version UMH April 96  w_meson.c
   cmn Nov 1996: added momentum phase factors, added ANSI protoptypes.
   C. DeTar 4/29/97 Extended source and sink meson types.
                    Switched to more general version of 
		    gamma matrix multiplication 
   C. DeTar 4/30/97 Takes a table of source and sink gamma matrices
   C. DeTar 5/24/97 Operator indexing now based on corrlist table
   C. DeTar 4/21/98 (lean version) Provision for storing
                    values only for time slices this node owns.
   C. DeTar 4/28/98 Optimized.
   C. DeTar 5/24/98 Further optimization.  
                    Essentially same behavior as w_meson_finite_mom_lean2.c
   C. DeTar 10/13/01 Prefetching		    
   C. DeTar 5/31/07  Modified to use with generalized clover_invert 
*/
#include "generic_wilson_includes.h"
#include <string.h>
#include "../include/prefetch.h"
#define FETCH_UP 1
#define INLINE
#define loopend sites_on_node   /* No checkerboarding here */

#define WMTIME

/* Local types */

/* Plan Dirac matrix */

typedef struct {
  complex d[4];
} dirac_spinor;

typedef struct {
  dirac_spinor d[4];
} dirac_matrix;
  
/* Dirac matrix organized so each element is a vector
   indexed by the possible q values.  Composite structure
   is foo.d[si].d[sf].e[p], where si and sf label
   rows and columns on the Dirac spin basis and p
   labels momentum indices */

#define MAXQ 100
typedef struct {
  complex e[MAXQ];
} dirac_element_v;

typedef struct {
  dirac_element_v d[4];
} dirac_spinor_v;

typedef struct {
  dirac_spinor_v d[4];
} dirac_matrix_v;

/*****************************************

This is a general function that can contract two quark
propagators together to form a meson correlator.

\sum_{x} exp( i p .x } 
    trace( \Gamma_in  S_2 \Gamma_out  \gamma_5  S_1^{\dagger} \gamma_5  )

   where S_1 and S_2 are quark propagators

   Note: the first spin index on the spin_wilson_vector refers to the
   source and the second to the sink.  Thus the trace above is
   unconventional.  Gamma_in operates on the source and Gamma_out on
   the sink.

  Function arguments

     On input 
        src1 :: spin_wilson_vector 
        src2 :: spin_wilson_vector 
        no_q_momenta :: number of momentum values p in table to use
	q_momstore[p] :: table of momentum values p to use
        q_parity[p] :: reflection parity for each momentum component
        no_gamma_corr :: number of source/sink gamma pairs G in table
        num_corr_mom[g] :: number of momentum/parity sets for each g
        corr_table[g] :: list of correlator indices c for gamma pair g
        p_index[c] :: p = p_index[c] is the momentum index for correlator c 
        gout[c] :: sink gamma matrix for correlator c
        gin[c] :: source gamma matrix for correlator c
        meson_phase[c] :: phase factor to multiply the correlator before
                          accumulating.  Encoded as in gammatypes.h
        meson_factor[c] :: normalization factor for each correlator c
        corr_index[c] :: correlator index m where like propagators are summed

     On output
         prop :: complex vector to the data correlators with indexing
                 prop[m][time] where m is the correlator index.


*******************************************/

#if 0

/******************  mult_swv_an.c ******************************
*									*
* dirac_matrix mult_swv_an(a,b) spin_wilson_vector *a,*b;		*
* return matrix product c = adj a * b of two spin wilson_vectors       	*
* traced over the color indices                                         *
*/

static dirac_matrix mult_swv_an( spin_wilson_vector *a, 
				 spin_wilson_vector *b )
{
  
  dirac_matrix c;
  
#ifndef INLINE
  complex temp1,temp2;
  register int si,sf,s;
  for(si=0;si<4;si++)for(sf=0;sf<4;sf++) {
    temp1.real = temp1.imag = 0.0;
    for(s=0;s<4;s++){
      CMULJ_(a->d[s].d[si].c[0],b->d[s].d[sf].c[0],temp2); CSUM(temp1,temp2);
      CMULJ_(a->d[s].d[si].c[1],b->d[s].d[sf].c[1],temp2); CSUM(temp1,temp2);
      CMULJ_(a->d[s].d[si].c[2],b->d[s].d[sf].c[2],temp2); CSUM(temp1,temp2);
    }
    c.d[si].d[sf] = temp1;
  }
  
  return c;

#else
  
#ifdef NATIVEDOUBLE
  register double ar,ai,br,bi,cr,ci;
#else
  register Real ar,ai,br,bi,cr,ci;
#endif

  register int si,sf;
  for(si=0;si<4;si++)for(sf=0;sf<4;sf++) {
    ar=a->d[0].d[si].c[0].real;  ai=a->d[0].d[si].c[0].imag;
    br=b->d[0].d[sf].c[0].real;  bi=b->d[0].d[sf].c[0].imag;
    cr = ar*br + ai*bi;
    ci = ar*bi - ai*br;
    ar=a->d[0].d[si].c[1].real;  ai=a->d[0].d[si].c[1].imag;
    br=b->d[0].d[sf].c[1].real;  bi=b->d[0].d[sf].c[1].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[0].d[si].c[2].real;  ai=a->d[0].d[si].c[2].imag;
    br=b->d[0].d[sf].c[2].real;  bi=b->d[0].d[sf].c[2].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    
    ar=a->d[1].d[si].c[0].real;  ai=a->d[1].d[si].c[0].imag;
    br=b->d[1].d[sf].c[0].real;  bi=b->d[1].d[sf].c[0].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[1].d[si].c[1].real;  ai=a->d[1].d[si].c[1].imag;
    br=b->d[1].d[sf].c[1].real;  bi=b->d[1].d[sf].c[1].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[1].d[si].c[2].real;  ai=a->d[1].d[si].c[2].imag;
    br=b->d[1].d[sf].c[2].real;  bi=b->d[1].d[sf].c[2].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    
    ar=a->d[2].d[si].c[0].real;  ai=a->d[2].d[si].c[0].imag;
    br=b->d[2].d[sf].c[0].real;  bi=b->d[2].d[sf].c[0].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[2].d[si].c[1].real;  ai=a->d[2].d[si].c[1].imag;
    br=b->d[2].d[sf].c[1].real;  bi=b->d[2].d[sf].c[1].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[2].d[si].c[2].real;  ai=a->d[2].d[si].c[2].imag;
    br=b->d[2].d[sf].c[2].real;  bi=b->d[2].d[sf].c[2].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    
    ar=a->d[3].d[si].c[0].real;  ai=a->d[3].d[si].c[0].imag;
    br=b->d[3].d[sf].c[0].real;  bi=b->d[3].d[sf].c[0].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[3].d[si].c[1].real;  ai=a->d[3].d[si].c[1].imag;
    br=b->d[3].d[sf].c[1].real;  bi=b->d[3].d[sf].c[1].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    ar=a->d[3].d[si].c[2].real;  ai=a->d[3].d[si].c[2].imag;
    br=b->d[3].d[sf].c[2].real;  bi=b->d[3].d[sf].c[2].imag;
    cr += ar*br + ai*bi;
    ci += ar*bi - ai*br;
    
    c.d[si].d[sf].real = cr;
    c.d[si].d[sf].imag = ci;
  }

  return c;
  
#endif
}

#endif

static dirac_matrix mult_swv_na( spin_wilson_vector *a, 
				 spin_wilson_vector *b )
{
  
  dirac_matrix c;
  
#ifndef INLINE
  complex temp1,temp2;
  register int si,sf,s;
  // flops 8*3*4*4*4 = 1536
  for(si=0;si<4;si++)for(sf=0;sf<4;sf++) {
    temp1.real = temp1.imag = 0.0;
    for(s=0;s<4;s++){
      CMUL_J(a->d[si].d[s].c[0],b->d[sf].d[s].c[0],temp2); CSUM(temp1,temp2);
      CMUL_J(a->d[si].d[s].c[1],b->d[sf].d[s].c[1],temp2); CSUM(temp1,temp2);
      CMUL_J(a->d[si].d[s].c[2],b->d[sf].d[s].c[2],temp2); CSUM(temp1,temp2);
    }
    c.d[si].d[sf] = temp1;
  }
  
  return c;

#else
  
#ifdef NATIVEDOUBLE
  register double ar,ai,br,bi,cr,ci;
#else
  register Real ar,ai,br,bi,cr,ci;
#endif

  register int si,sf;
  for(si=0;si<4;si++)for(sf=0;sf<4;sf++) {
    ar=a->d[si].d[0].c[0].real;  ai=a->d[si].d[0].c[0].imag;
    br=b->d[sf].d[0].c[0].real;  bi=b->d[sf].d[0].c[0].imag;
    cr = br*ar + bi*ai;
    ci = br*ai - bi*ar;
    ar=a->d[si].d[0].c[1].real;  ai=a->d[si].d[0].c[1].imag;
    br=b->d[sf].d[0].c[1].real;  bi=b->d[sf].d[0].c[1].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[0].c[2].real;  ai=a->d[si].d[0].c[2].imag;
    br=b->d[sf].d[0].c[2].real;  bi=b->d[sf].d[0].c[2].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    
    ar=a->d[si].d[1].c[0].real;  ai=a->d[si].d[1].c[0].imag;
    br=b->d[sf].d[1].c[0].real;  bi=b->d[sf].d[1].c[0].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[1].c[1].real;  ai=a->d[si].d[1].c[1].imag;
    br=b->d[sf].d[1].c[1].real;  bi=b->d[sf].d[1].c[1].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[1].c[2].real;  ai=a->d[si].d[1].c[2].imag;
    br=b->d[sf].d[1].c[2].real;  bi=b->d[sf].d[1].c[2].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    
    ar=a->d[si].d[2].c[0].real;  ai=a->d[si].d[2].c[0].imag;
    br=b->d[sf].d[2].c[0].real;  bi=b->d[sf].d[2].c[0].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[2].c[1].real;  ai=a->d[si].d[2].c[1].imag;
    br=b->d[sf].d[2].c[1].real;  bi=b->d[sf].d[2].c[1].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[2].c[2].real;  ai=a->d[si].d[2].c[2].imag;
    br=b->d[sf].d[2].c[2].real;  bi=b->d[sf].d[2].c[2].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    
    ar=a->d[si].d[3].c[0].real;  ai=a->d[si].d[3].c[0].imag;
    br=b->d[sf].d[3].c[0].real;  bi=b->d[sf].d[3].c[0].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[3].c[1].real;  ai=a->d[si].d[3].c[1].imag;
    br=b->d[sf].d[3].c[1].real;  bi=b->d[sf].d[3].c[1].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    ar=a->d[si].d[3].c[2].real;  ai=a->d[si].d[3].c[2].imag;
    br=b->d[sf].d[3].c[2].real;  bi=b->d[sf].d[3].c[2].imag;
    cr += br*ar + bi*ai;
    ci += br*ai - bi*ar;
    
    c.d[si].d[sf].real = cr;
    c.d[si].d[sf].imag = ci;
  }
  
  return c;
  
#endif
}


/* Multiply a vector Dirac matrix by a Dirac gamma on the left
   and take the trace */

static double dirac_v_tr_gamma(complex *tr, dirac_matrix_v *src, 
			       int gamma, int phase[], Real factor[],
			       int ct[], int nc, int p_index[])
{
  int s2,s; /* spin indices */
  int k, c, p, ph;
  complex z = {0.,0.};
  Real fact;
  double flops = 0;
  
  /* One value for each momentum */
  for(k=0; k<nc; k++)
    tr[k].real = tr[k].imag = 0;

  /* For each momentum in list, multiply by gamma and take trace */

  for(s=0;s<4;s++){
    s2 = gamma_mat(gamma).row[s].column;
    switch (gamma_mat(gamma).row[s].phase){
    case 0:
      for(k=0; k<nc; k++){
	c = ct[k];
	p = p_index[c];
	z =            src->d[s2].d[s].e[p];
	CSUM(tr[k],z);
      }
      break;
    case 1:
      for(k=0; k<nc; k++){
	c = ct[k];
	p = p_index[c];
	TIMESPLUSI(    src->d[s2].d[s].e[p], z);
	CSUM(tr[k],z);
      }
      break;
    case 2:
      for(k=0; k<nc; k++){
	c = ct[k];
	p = p_index[c];
	TIMESMINUSONE( src->d[s2].d[s].e[p], z);
	CSUM(tr[k],z);
      }
      break;
    case 3:
      for(k=0; k<nc; k++){
	c = ct[k];
	p = p_index[c];
	TIMESMINUSI(   src->d[s2].d[s].e[p], z);
	CSUM(tr[k],z);
      }
    }
  }

  /* Normalization and phase */

  for(k=0; k<nc; k++){
    c = ct[k];
    ph = phase[c];
    fact = factor[c];
    switch(ph){
    case 0:
      z =            tr[k];
      break;
    case 1:
      TIMESPLUSI(    tr[k], z);
      break;
    case 2:
      TIMESMINUSONE( tr[k], z);
      break;
    case 3:
      TIMESMINUSI(   tr[k], z);
    }
    CMULREAL(z,fact,tr[k]);
  }

  flops = 4*nc;
  
  return flops;
  
} /* dirac_v_tr_gamma */

/* Calculate FT weight factor */

static complex ff(Real theta, char parity, complex tmp)
{
  complex z = {0.,0.};
  
  if(parity == EVEN){
    z.real = tmp.real*cos(theta);
    z.imag = tmp.imag*cos(theta);
  }
  else if(parity == ODD){
    z.real = -tmp.imag*sin(theta);
    z.imag =  tmp.real*sin(theta);
  }
  else if(parity == EVENANDODD){
    z.real = tmp.real*cos(theta)-tmp.imag*sin(theta);
    z.imag = tmp.imag*cos(theta)+tmp.real*sin(theta);
  }
  else{
    printf("ff(%d): bad parity %d\n", this_node, parity);
    terminate(1);
  }
  return z;
} /* ff */

/********************* meson_cont_mom **************************/

void meson_cont_mom(
  complex **prop,           /* prop[m][t] is where result is accumulated */
  spin_wilson_vector *src1, /* quark propagator (to become antiquark) */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_gamma_corr,        /* # of gamma src/snk combinations (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values for each corr (gt k) */
  int **corr_table,         /* c = corr_table[g][k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  int gout[],               /* gout[c] is the sink gamma */
  int gin[],                /* gin[c] is the source gamma */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_factor[],      /* meson_factor[c] scales the correlator */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* spatial origin for defining FT phases */
		    )
{
  char myname[] = "meson_cont_mom";
  int i,k,c,m,gsnk,gsrc;
  site *s; 
  
  int sf, si;
  int g,p,t;
  int old_gamma_out;
  
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 
  int px,py,pz;
  char ex, ey, ez;
  complex fourier_fact ; 
  
  spin_wilson_vector localmat;  /* temporary storage */
  spin_wilson_vector antiquark; /* temporary storage for antiquark */
  
  dirac_matrix *meson;
  dirac_matrix_v *meson_q;
  int *nonzero;
  complex tr[MAXQ];
  complex tmp;
  complex *ftfact;
  gamma_matrix_t gm5, gmout, gmoutadj, gm;
  
  /* performance */
  double dtime;
  double flops,mflops;

  dtime = -dclock();
  flops = 0;

  /* Allocate temporary storage for meson propagator */
  
  if(no_q_momenta > MAXQ)
    {
      printf("%s(%d): no_q_momenta %d exceeds max %d\n",
	     myname, this_node, no_q_momenta,MAXQ);
      terminate(1);
    }
  
  meson = (dirac_matrix *)malloc(sites_on_node*sizeof(dirac_matrix));
  if(meson == NULL){
    printf("%s(%d): No room for meson\n",myname,this_node);
    terminate(1);
  }
  
  meson_q = (dirac_matrix_v *)malloc(nt*sizeof(dirac_matrix_v));
  if(meson_q == NULL){
    printf("%s(%d): No room for meson_q\n",myname,this_node);
    terminate(1);
  }
  
  nonzero = (int *)malloc(nt*sizeof(int));
  if(nonzero == NULL){
    printf("%s(%d): No room for nonzero array\n",myname,this_node);
    terminate(1);
  }
  
  for(t = 0; t < nt; t++)nonzero[t] = 0;
  
  ftfact = (complex *)malloc(no_q_momenta*sites_on_node*sizeof(complex));
  if(ftfact == NULL)
    {
      printf("%s(%d): No room for FFT phases\n",myname,this_node);
      terminate(1);
    }

  /* ftfact contains factors such as cos(kx*x)*sin(ky*y)*exp(ikz*z)
     with factors of cos, sin, and exp selected according to the
     requested component parity */
  
  FORALLSITES(i,s) {
    for(p=0; p<no_q_momenta; p++)
      {
	px = q_momstore[p][0];
	py = q_momstore[p][1];
	pz = q_momstore[p][2];
	
	ex = q_parity[p][0];
	ey = q_parity[p][1];
	ez = q_parity[p][2];
	
	tmp.real = 1.;
	tmp.imag = 0.;
	
	tmp = ff(factx*(s->x-r0[0])*px, ex, tmp);
	tmp = ff(facty*(s->y-r0[1])*py, ey, tmp);
	tmp = ff(factz*(s->z-r0[2])*pz, ez, tmp);
	
	ftfact[p+no_q_momenta*i] = tmp;
      }
  }      
  
  flops += (double)sites_on_node*18*no_q_momenta;
  
  
  /* Run through the table of unique source-sink gamma matrix pairs
     and for each, do the Fourier transform according to the list
     of requested momenta */
  
  /* To help suppress repetition */
  old_gamma_out = -999;
  
  for(g = 0; g < no_gamma_corr; g++)
    {

      /*  All gammas with the same index g must be the same */
      c = corr_table[g][0];  
      gsrc = gin[c];
      gsnk = gout[c];

      /* For compatibility */
      if(gsrc == GAMMAFIVE)gin[c] = G5;
      if(gsnk == GAMMAFIVE)gout[c] = G5;
      
      if(gsrc >= MAXGAMMA || gsnk >= MAXGAMMA)
	{
	  printf("%s(%d): Illegal gamma index %d or %d\n",
		 myname, this_node, gsrc, gsnk);
	  terminate(1);
	}

      for(k=0; k<num_corr_mom[g]; k++)
	{
	  if(gin[corr_table[g][k]] != gsrc ||
	     gout[corr_table[g][k]] != gsnk )
	    {
	      printf("meson_cont_mom(%d): bad gamma list\n", this_node);
	      terminate(1);
	    }
	}
      
      /* Skip reconstruction of "meson" if "out" gamma matrix hasn't changed */
      if(gsnk != old_gamma_out)
	{
	  old_gamma_out = gsnk;

	  /* Compute gm = \gamma_5 Gamma_out */
	  gm5  = gamma_mat(G5);
	  gmout = gamma_mat(gsnk);
	  gamma_adj(&gmoutadj, &gmout);
	  mult_gamma_by_gamma(&gm5, &gmoutadj, &gm);
	    
	  FORALLSITES(i,s) {
	    
	    /* antiquark = gamma_5 adjoint of quark propagator gamma_5 */
	    /* But we use a complex matrix dot product below, so don't
	       take the c.c. here and do the transpose later 
	       - so just do a gamma5 transformation now */

	    if( i < loopend-FETCH_UP){
	      prefetch_WWWW(
		    &(src1[i].d[0]), &(src1[i].d[1]),
		    &(src1[i].d[2]), &(src1[i].d[3]));
	      prefetch_WWWW(
		    &(src2[i].d[0]), &(src2[i].d[1]),
		    &(src2[i].d[2]), &(src2[i].d[3]));
	    }

	    /* antiquark = \gamma_5 S_1 \gamma_5 \Gamma_out^\dagger */

	    mult_sw_by_gamma_l( src1+i, &localmat, G5);    
	    mult_sw_by_gamma_mat_r( &localmat, &antiquark, &gm);     
	    
	    
	    /* combine with src2 quark and sew together sink colors 
	       and spins to make meson Dirac matrix 
	       "meson" is then [S_2 \Gamma_out \gamma_5 S_1^\dagger \gamma_5]
	       and is a Dirac matrix */
	    
	    meson[i] = mult_swv_na( src2+i, &antiquark);
	    
	  } /* end FORALLSITES */
	  flops += (double)sites_on_node*1536;
	} /* end if not same Gamma_out */
      
      /* Do FT on "meson" for momentum projection - 
	 Result in meson_q.  We use a dumb FT because there 
         are so few momenta needed. */
      
      FORALLSITES(i,s) {
	for(si = 0; si < 4; si++)for(sf = 0; sf < 4; sf++)
	   for(k=0; k<num_corr_mom[g]; k++)
	     {
	       c = corr_table[g][k];
	       p = p_index[c];
	       
	       meson_q[s->t].d[si].d[sf].e[p].real = 0.;
	       meson_q[s->t].d[si].d[sf].e[p].imag = 0;
	     }
      }
      
      FORALLSITES(i,s) {
	nonzero[s->t] = 1;  /* To save steps below */
	for(si = 0; si < 4; si++)
	  for(sf = 0; sf < 4; sf++)
	    {
	      for(k=0; k<num_corr_mom[g]; k++)
		{
		  c = corr_table[g][k];
		  p = p_index[c];
		  fourier_fact = ftfact[p+no_q_momenta*i];

		  meson_q[s->t].d[si].d[sf].e[p].real += 
		    meson[i].d[si].d[sf].real*fourier_fact.real -  
		    meson[i].d[si].d[sf].imag*fourier_fact.imag;
		  meson_q[s->t].d[si].d[sf].e[p].imag += 
		    meson[i].d[si].d[sf].real*fourier_fact.imag +  
		    meson[i].d[si].d[sf].imag*fourier_fact.real;
		}
	    }
      }
	
      flops += (double)sites_on_node*128*num_corr_mom[g];
      
      /* Complete the propagator by tying in the sink gamma.
         Then store it */
      
      for(t=0; t < nt; t++)if(nonzero[t]) {
	  /* Do final Dirac trace for all sink momenta q */
	  flops += dirac_v_tr_gamma(tr, &meson_q[t], gsrc,
				    meson_phase, meson_factor, corr_table[g], 
				    num_corr_mom[g], p_index);
	  /* Accumulate in corr_index location */
	  for(k=0; k<num_corr_mom[g]; k++)
	    {
	      c = corr_table[g][k];
	      m = corr_index[c];
	      prop[m][t].real += tr[k].real;
	      prop[m][t].imag += tr[k].imag;
	    }
	}
    }  /**** end of the loop over gamma table ******/
  
  free(meson);  free(meson_q);  free(nonzero);  free(ftfact);
  
  dtime += dclock();

/**  flops = sites_on_node*(nmeson_evals*1536 + nmeson_q_evals*128) +
     68*no_q_momenta) + 14*nprop_evals*ntslices*no_q_momenta; **/

#ifdef WMTIME
  if(dtime > 0)mflops = flops/(dtime*1e6);
  else mflops = 0;
  
  node0_printf("WMTIME: time %.1e sec %g flops %.1f MF\n",
	       dtime,flops,mflops);fflush(stdout);
#endif
  
} /* end of meson_cont_mom function  */

