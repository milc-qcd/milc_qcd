/*********************** w_meson_open_mom.c *****************************/
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
#include "../include/gammatypes.h"
#define FETCH_UP 1
#define INLINE
#define loopend sites_on_node   /* No checkerboarding here */

#define WMTIME

/* Local types */

/* Plan Dirac matrix */

/* Wilson propagator organized so each element is a vector
   indexed by the possible q values.  Composite structure
   is foo.c[ci].d[si].d[sf].c[cf].e[p], where ci,si and cf,sf label
   rows and columns on the Dirac color-spin basis and p
   labels momentum indices */

#define MAXQ 100
typedef struct {
  complex e[MAXQ];
} element_v;

typedef struct {
  element_v c[3];
} color_vector_v;

typedef struct {
  color_vector_v d[4];
} wilson_vector_v;

typedef struct {
  wilson_vector_v d[4];
} spin_wilson_vector_v;

typedef struct {
  spin_wilson_vector_v c[3];
} wilson_propagator_v;

/*****************************************

This is a general function that can contract two quark propagators
together to form an "open meson" correlator.  i.e. the sink color and
spin indices are left open and the sink gamma matrix is ignored.

\sum_{x} exp( i p .x } 
       \gamma_5  S_1^{\dagger} \gamma_5 \Gamma_in S_2

   where S_1 and S_2 are quark propagators

   Note: the first spin index on the spin_wilson_vector refers to the
   source and the second to the sink.  Thus the trace above is
   unconventional.  Gamma_in operates on the source.

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
        gin[c] :: source gamma matrix for correlator c
        gout[c] :: sink gamma matrix for correlator c (ignored)
        meson_phase[c] :: phase factor to multiply the correlator before
                          accumulating.  Encoded as in gammatypes.h
        meson_factor[c] :: normalization factor for each correlator c
        corr_index[c] :: correlator index m where like propagators are summed

     On output
         prop :: wilson_propagator vector of open correlators with indexing
                 prop[m][time] where m is the correlator index.


*******************************************/

static wilson_propagator open_mult_swv_an( spin_wilson_vector *a, 
					   spin_wilson_vector *b )
{
  
  wilson_propagator c;
  
#ifndef INLINE
  int si,sf,ci,cf,s;
  complex temp1,temp2;
  for(si=0;si<4;si++)
    for(ci=0;ci<3;ci++)
      for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++){
	  temp1.real = temp1.imag = 0.0;
	  for(s=0;s<4;s++){
	    CMULJ_(a->d[s].d[si].c[ci],b->d[s].d[sf].c[cf],temp2); 
	    CSUM(temp1,temp2);
	  }
	  c.c[ci].d[si].d[sf].c[cf] = temp1;
	}
  return c;

#else
  
#ifdef NATIVEDOUBLE
  register double ar,ai,br,bi,dr,di;
#else
  register Real ar,ai,br,bi,dr,di;
#endif

  register int si,sf,ci,cf;
  for(si=0;si<4;si++)
    for(ci=0;ci<3;ci++)
      for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++) {
	  ar=a->d[0].d[si].c[ci].real;  ai=a->d[0].d[si].c[ci].imag;
	  br=b->d[0].d[sf].c[cf].real;  bi=b->d[0].d[sf].c[cf].imag;
	  dr = br*ar + bi*ai;
	  di = bi*ar - br*ai;
	  
	  ar=a->d[1].d[si].c[ci].real;  ai=a->d[1].d[si].c[ci].imag;
	  br=b->d[1].d[sf].c[cf].real;  bi=b->d[1].d[sf].c[cf].imag;
	  dr += br*ar + bi*ai;
	  di += bi*ar - br*ai;
	  
	  ar=a->d[2].d[si].c[ci].real;  ai=a->d[2].d[si].c[ci].imag;
	  br=b->d[2].d[sf].c[cf].real;  bi=b->d[2].d[sf].c[cf].imag;
	  dr += br*ar + bi*ai;
	  di += bi*ar - br*ai;
	  
	  ar=a->d[3].d[si].c[ci].real;  ai=a->d[3].d[si].c[ci].imag;
	  br=b->d[3].d[sf].c[cf].real;  bi=b->d[3].d[sf].c[cf].imag;
	  dr += br*ar + bi*ai;
	  di += bi*ar - br*ai;
	  
	  c.c[ci].d[si].d[sf].c[cf].real = dr;
	  c.c[ci].d[si].d[sf].c[cf].imag = di;
	}
  
  return c;
  
#endif
}

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

/* Calculate normalization factor with phase */

static complex phase_and_factor(int ph, Real fact){
  complex z = {0.,0.};
  complex f = cmplx(fact,0);

  switch(ph){
  case 0:
    z =            f;
    break;
  case 1:
    TIMESPLUSI(    f, z);
    break;
  case 2:
    TIMESMINUSONE( f, z);
    break;
  case 3:
    TIMESMINUSI(   f, z);
  }
  return z;
}


/********************* meson_open_mom **************************/

void meson_open_mom(
  wilson_propagator **prop, /* prop[m][t] is where result is accumulated */
  spin_wilson_vector *src1, /* quark propagator (to become antiquark) */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_gamma_corr,        /* # of gamma src/snk combinations (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values for each corr */
  int **corr_table,         /* c = corr_table[g][k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  int gout[],               /* gout[c] is the sink gamma (ignored) */
  int gin[],                /* gin[c] is the source gamma */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_factor[],      /* meson_factor[c] scales the correlator */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* spatial origin for defining FT phases */

		    )
{
  char myname[] = "meson_open_mom";
  int i,k,c,m,gsrc;
  site *s; 
  
  int ci, cf, sf, si;
  int g,p,t;
  int old_gamma_in;
  
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 
  int px,py,pz;
  char ex, ey, ez;
  complex fourier_fact ; 
  
  spin_wilson_vector localmat;  /* temporary storage */
  spin_wilson_vector antiquark; /* temporary storage for antiquark */
  
  wilson_propagator *meson;
  wilson_propagator_v *meson_q;
  int *nonzero;
  complex tmp;
  complex *ftfact;
  gamma_matrix_t gm5, gmin, gm;
  
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
  
  meson = (wilson_propagator *)malloc(sites_on_node*sizeof(wilson_propagator));
  if(meson == NULL){
    printf("%s(%d): No room for meson\n",myname,this_node);
    terminate(1);
  }
  
  meson_q = (wilson_propagator_v *)malloc(nt*sizeof(wilson_propagator_v));
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
  
  flops += (double)sites_on_node*68*no_q_momenta;
  
  
  /* Run through the table of unique source-sink gamma matrix pairs
     and for each, do the Fourier transform according to the list
     of requested momenta */
  
  /* To help suppress repetition */
  old_gamma_in = -999;
  
  for(g = 0; g < no_gamma_corr; g++)
    {

      c = corr_table[g][0];  
      gsrc = gin[c];

      /* For compatibility */
      if(gsrc == GAMMAFIVE)gsrc = G5;
      
      if(gsrc >= MAXGAMMA)
	{
	  printf("%s(%d): Illegal gamma index %d\n",
		 myname,this_node,gsrc);
	  terminate(1);
	}

      /*  All gammas with the same index g must be the same */
      for(k=0; k<num_corr_mom[g]; k++)
	{
	  if(gin[corr_table[g][k]] != gsrc )
	    {
	      printf("%s(%d): bad gamma list\n", myname, this_node);
	      terminate(1);
	    }
	}
      
      /* Skip reconstruction of "meson" if "in" gamma matrix hasn't changed */
      if(gsrc != old_gamma_in)
	{
	  old_gamma_in = gsrc;

	  /* Compute gm = \gamma_5 Gamma_in */
	  gm5  = gamma_mat(G5);
	  gmin = gamma_mat(gsrc);
	  mult_gamma_by_gamma(&gm5, &gmin, &gm);
	    
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

	    /* right multiply src1 by gamma5 to make antiquark */

	    mult_sw_by_gamma_r( src1+i, &antiquark, G5);    
	    
	    /* left multiply src2 by gamma5 Gamma_in to make localmat */
	    mult_sw_by_gamma_mat_l( src2+i, &localmat, &gm);
	    
	    /* trace over source spins and colors to make open meson
	       propagator "meson" is then \gamma_5 S_1^\dagger
	       \gamma_5 \Gamma_in S_2 */
	    
	    meson[i] = open_mult_swv_an(&antiquark, &localmat);
	    
	  } /* end FORALLSITES */
	  flops += (double)sites_on_node*1536;  /* Needs correcting */
	} /* end if not same Gamma_in */
      
      /* Do FT on "meson" for momentum projection - 
	 Result in meson_q.  We use a dumb FT because there 
         are so few momenta needed. */
      
      FORALLSITES(i,s) {
	for(ci = 0; ci < 3; ci++)
	  for(si = 0; si < 4; si++)
	    for(sf = 0; sf < 4; sf++)
	      for(cf = 0; cf < 3; cf++)
		for(k=0; k<num_corr_mom[g]; k++)
		  {
		    c = corr_table[g][k];
		    p = p_index[c];
		    
		    meson_q[s->t].c[ci].d[si].d[sf].c[cf].e[p].real = 0.;
		    meson_q[s->t].c[ci].d[si].d[sf].c[cf].e[p].imag = 0;
		  }
      }
      
      FORALLSITES(i,s) {
	nonzero[s->t] = 1;  /* To save steps below in case we have
			       few time slices on this node */
	for(ci=0;ci<3;ci++)
	  for(si=0;si<4;si++)
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {
		{
		  for(k=0; k<num_corr_mom[g]; k++)
		    {
		      c = corr_table[g][k];
		      p = p_index[c];
		      fourier_fact = ftfact[p+no_q_momenta*i];
		      
		      CMUL(meson[i].c[ci].d[si].d[sf].c[cf], fourier_fact, tmp);
		      CSUM(meson_q[s->t].c[ci].d[si].d[sf].c[cf].e[p], tmp);
		    }
		}
	      }
      }
	
      flops += (double)sites_on_node*1152*num_corr_mom[g];
      
      /* Store the result */
      
      for(t=0; t < nt; t++)if(nonzero[t]) {
	  /* Accumulate in corr_index location */
	  for(k=0; k<num_corr_mom[g]; k++)
	    {
	      complex y,z;

	      c = corr_table[g][k];
	      m = corr_index[c];
	      p = p_index[c];
	      z = phase_and_factor(meson_phase[c],meson_factor[c]);
	      for(ci=0;ci<3;ci++)
		for(si=0;si<4;si++)
		  for(sf=0;sf<4;sf++)
		    for(cf=0;cf<3;cf++) {
		      /* Adjust phase and normalization */
		      CMUL(z,meson_q[t].c[ci].d[si].d[sf].c[cf].e[p],y);
		      /* Accumulate resultant open meson propagator */
		      CSUM(prop[m][t].c[ci].d[si].d[sf].c[cf],y);
		    }
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
  
} /* end of meson_open_mom function  */



