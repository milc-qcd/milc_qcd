/************* w_meson_finite_mom_lean3.c **************************/
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
   is foo.d[si].d[sf].e[q_pt], where si and sf label
   rows and columns on the Dirac spin basis and q_pt
   labels momentum transfer values */

#define MAXQ 48
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
    trace( \Gamma_in  \gamma_5  S_1^{\dagger} \gamma_5 \Gamma_out  S_2 )

   where S_1 and S_2 are quark propagators

   Note: the first spin index on the spin_wilson_vector refers to the
   source and the second to the sink.  Thus the trace above is
   unconventional.  Gamma_out operates on the source and Gamma_in on
   the sink.  In more conventional notation, where S = < \psi(t) \psibar(0)>
   we should replace S_1 and S_2 by their transposes.  Then taking the
   transpose inside the trace gives

    trace( \Gamma_out^T \gamma_5 S_1^\dagger \gamma_5 \Gamma_in^T S_2 )

   where now the first index on S is the sink index and the second the
   source.

  Function arguments

     On input 
        src1 :: spin_wilson_vector 
        src2 :: spin_wilson_vector 
        no_q_momenta :: number of momentum values in table to use
	q_momstore :: table of momentum values to use
        no_gamma_corr :: number of entries in table to use
        meson_index :: where to accumulate the propagator
        gamma_table :: table of vertex and source gamma matrices to use
        meson_phase :: phase factor to multiply the correlator before
                       accumulating
     On output
         prop :: complex vector to the data correlators with indexing
                 prop[meson_type][momentum][time]


*******************************************/

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

/* Multiply a vector Dirac matrix by a Dirac gamma 
   and take the trace */

static void dirac_v_tr_gamma(complex *tr, dirac_matrix_v *src, 
			     int gamma, int nq, complex phase)
{
  register int s2,s,iq;	/* spin indices */
  register complex z;
  
  /* For compatibility */
  if(gamma == GAMMAFIVE)gamma = G5;
  if(gamma >= MAXGAMMA)
    {
      printf("dirac_v_tr_gamma: Illegal gamma index %d\n",gamma);
      exit(1);
    }

  for(iq = 0; iq < nq; iq++){
    tr[iq].real = tr[iq].imag = 0;
  }

  for(s=0;s<4;s++){
    s2 = gamma_mat(gamma).row[s].column;
    switch (gamma_mat(gamma).row[s].phase){
    case 0:
      for(iq=0;iq<nq;iq++){
	z =            src->d[s2].d[s].e[iq];
	CSUM(tr[iq],z);
      }
      break;
    case 1:
      for(iq=0;iq<nq;iq++){
	TIMESPLUSI(    src->d[s2].d[s].e[iq], z);
	CSUM(tr[iq],z);
      }
      break;
    case 2:
      for(iq=0;iq<nq;iq++){
	TIMESMINUSONE( src->d[s2].d[s].e[iq], z);
	CSUM(tr[iq],z);
      }
      break;
    case 3:
      for(iq=0;iq<nq;iq++){
	TIMESMINUSI(   src->d[s2].d[s].e[iq], z);
	CSUM(tr[iq],z);
      }
    }
  }

  for(iq = 0; iq < nq; iq++){
    z = tr[iq];
    CMUL(z,phase,tr[iq]);
  }

} /* dirac_v_tr_gamma.c */

/********************* meson_cont_mom **************************/

void meson_cont_mom(
  complex ***prop,          /* where result is stored */
  spin_wilson_vector *src1, /* quark propagator */
  spin_wilson_vector *src2, /* quark propagator */
  int no_q_momenta,         /* number of q values */
  int q_momstore[][3],      /* q values themselves */
  int mom_index[],          /* momentum hash */
  int no_gamma_corr,        /* number of meson types */
  int meson_index[],        /* meson hash */
  int gin[],                /* Gamma matrix type for source */
  int gout[],               /* Gamma matrix type for sink */
  complex meson_phase[]     /* phase factor to apply to correlator */
		    )
{
  char myname[] = "meson_cont_mom";
  int i, j;
  site *s; 
  
  int sf, si;
  int i_gamma_corr,q_pt,t;
  int old_gamma_out;
  
  double theta ; 
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 
  int px,py,pz;
  complex phase_fact ; 
  
  spin_wilson_vector localmat;  /* temporary storage */
  spin_wilson_vector antiquark; /* temporary storage for antiquark */
  
  dirac_matrix *meson;
  dirac_matrix_v *meson_q;
  int *nonzero;
  complex tr[MAXQ];
  complex *tmp;

  /* performance */
  double dtime;
  int nmeson_evals = 0;
  int nprop_evals = 0;
  int ntslices;
  double flops,mflops;

  dtime = -dclock();

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
  
  tmp = (complex *)malloc(no_q_momenta*sites_on_node*sizeof(complex));
  if(tmp == NULL)
    {
      printf("%s(%d): No room for FFT phases\n",myname,this_node);
      terminate(1);
    }

  FORALLSITES(i,s) {
    for(q_pt=0; q_pt<no_q_momenta; q_pt++)
      {
	px = q_momstore[q_pt][0];
	py = q_momstore[q_pt][1];
	pz = q_momstore[q_pt][2];
	
	theta = factx*(s->x)*px + facty*(s->y)*py + factz*(s->z)*pz; 
	tmp[q_pt+no_q_momenta*i].real = cos(theta);
	tmp[q_pt+no_q_momenta*i].imag = sin(theta); 
    }
  }      

  
  /* Run through the table of source-sink gamma matrices */

  for(i_gamma_corr = 0; i_gamma_corr < no_gamma_corr; i_gamma_corr++)
    {
      /* Skip reconstruction of "meson" if "out" gamma matrix hasn't changed */
      if(i_gamma_corr == 0 || gout[i_gamma_corr] != old_gamma_out)
	{
	  old_gamma_out = gout[i_gamma_corr];

	  nmeson_evals++;    /* performance */

	  FORALLSITES(i,s) {
	    
	    /* antiquark = gamma_5 adjoint of quark propagator gamma_5 */
	    /* But we use a complex dot product below, so don't
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

	    /* left multiply antiquark by source gamma matrices,
	       beginning with gamma_5 for quark -> antiquark */

	    mult_sw_by_gamma_l( src1+i, &localmat, G5);    
	    
	    /* right dirac multiplication by gamma-5 
	       (finishing up antiquark) */
	    mult_sw_by_gamma_r( &localmat, &antiquark, G5);     
	    
	    /* left multiply src2 by Gamma_out.  Result in localmat */
	    mult_sw_by_gamma_l( src2+i,	&localmat, gout[i_gamma_corr]);

	    /* combine with src2 quark and sew together colors 
	       and source spins to make meson Dirac matrix 
	    "meson" is then [\gamma_5 S_1 \gamma_5]^\dagger \Gamma_out S_ 2 
	    and is a Dirac matrix */

	    meson[i] = mult_swv_an(&antiquark,&localmat);
	    
	  } /* end FORALLSITES */
	  
	  /* Do FT on "meson" for momentum projection - 
	     Result in meson_q */
	  FORALLSITES(i,s) {
	    for(si = 0; si < 4; si++)for(sf = 0; sf < 4; sf++)
	      for(q_pt=0; q_pt<no_q_momenta; q_pt++)
		{
		  meson_q[s->t].d[si].d[sf].e[q_pt].real = 0.;
		  meson_q[s->t].d[si].d[sf].e[q_pt].imag = 0;
		}
	  }
	  
	  ntslices = 0;  /* For determining performance */
	  FORALLSITES(i,s) {
	    nonzero[s->t] = 1;  /* To save steps below */
	    ntslices++;
	    for(si = 0; si < 4; si++)for(sf = 0; sf < 4; sf++)
	      for(q_pt=0; q_pt<no_q_momenta; q_pt++)
		{
		  phase_fact = tmp[q_pt+no_q_momenta*i];
		  meson_q[s->t].d[si].d[sf].e[q_pt].real += 
		    meson[i].d[si].d[sf].real*phase_fact.real -  
		    meson[i].d[si].d[sf].imag*phase_fact.imag;
		  meson_q[s->t].d[si].d[sf].e[q_pt].imag += 
		    meson[i].d[si].d[sf].real*phase_fact.imag +  
		    meson[i].d[si].d[sf].imag*phase_fact.real;
		}
	  }
	} /* end if not same Gamma_out */
      
      nprop_evals++;   /* performance */
      
      /* Complete the propagator */
      
      for(t=0; t < nt; t++)if(nonzero[t]) {
	/* Do final Dirac trace for all sink momenta q */
	dirac_v_tr_gamma(tr, &meson_q[t], gin[i_gamma_corr],
			 no_q_momenta, meson_phase[i_gamma_corr]);
	/* Store values for all sink momenta */	
	for(q_pt=0; q_pt<no_q_momenta; q_pt++)
	  {
	    /* Accumulate in meson_index location */
	    i = meson_index[i_gamma_corr];
	    j = mom_index[q_pt];
	    prop[i][j][t].real += tr[q_pt].real;
	    prop[i][j][t].imag += tr[q_pt].imag;
	  }
      }
    }  /**** end of the loop over gamma table ******/
  
  free(meson);  free(meson_q);  free(nonzero);  free(tmp);

  dtime += dclock();

  flops = sites_on_node*(nmeson_evals*(1536 + 128*no_q_momenta) +
	       68*no_q_momenta) + 14*nprop_evals*ntslices*no_q_momenta;
  if(dtime > 0)mflops = flops/(dtime*1e6);
  else mflops = 0;

  node0_printf("%s(%d): time %.1f sec %.1f MF\n",myname, this_node,
	       dtime,mflops);fflush(stdout);
  
} /* end of meson_cont_mom function  */

static char *phaselabel[4] = { "1", "i", "-1", "-i" };
static complex phase[4] = { {1,0},  {0,1},  {-1,0},  {0,-1} };

/* Map a label to the corresponding phase */

complex decode_phase(char *label){
  int i;
  for(i = 0; i < 4; i++){
    if(strcmp(label,phaselabel[i]) == 0)return phase[i];
  }
  return cmplx(0,0);  /* Error condition */
}


