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
*/

#include "generic_form_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

/* Local types */


typedef struct {
  int gin;         /* Source gamma matrix */
  int gout;        /* Sink gamma matrix */
  int oper;        /* Operator identification */
} gamma_corr;

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
        src1 :: site structure pouinter to spin_wilson_vector 
        src2 :: site structure pouinter to spin_wilson_vector 
	base_pt:: partially computed offset for storage of result in prop array
	q_stride, op_stride:: stride for q index and operator index in prop array
	q_momstore :: table of momentum values to use
        no_q_momenta :: number of momentum values in table to use
        gamma_corr :: table of vertex and source gamma matrices to use
        no_gamma_corr :: number of entries in table to use

     On output
         prop :: complex vector to the data correlators


*******************************************/

/******************  mult_swv_an.c ******************************
*									*
* dirac_matrix mult_swv_an(a,b) spin_wilson_vector *a,*b;		*
* return matrix product c = adj a * b of two spin wilson_vectors       	*
* traced over the color indices                                         *
*/

dirac_matrix mult_swv_an( spin_wilson_vector *a, spin_wilson_vector *b ){

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

void dirac_v_tr_gamma(dirac_matrix_v *src, int gamma, 
		      complex *tr, int nq)
{
  register int s2,s,iq;	/* spin indices */
  register complex z;
  
  if(gamma_initialized==0)make_gammas(gamma_mat);
  
  /* For compatibility */
  if(gamma == GAMMAFIVE)gamma = G5;
  if(gamma >= MAXGAMMA)
    {
      printf("dirac_tr_gamma: Illegal gamma index %d\n",gamma);
      exit(1);
    }

  for(iq = 0; iq < nq; iq++){
    tr[iq].real = tr[iq].imag = 0;
  }

  for(s=0;s<4;s++){
    s2 = gamma_mat[gamma].row[s].column;
    switch (gamma_mat[gamma].row[s].phase){
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
} /* dv_tr_gamma.c */

/********************* meson_cont_mom **************************/

void meson_cont_mom_lean2(
  complex prop[],        /* where result is stored */
  field_offset src1,     /* quark propagator (spin_wilson_vector) */
  field_offset src2,     /* quark propagator (spin_wilson_vector) */
  int base_pt,           /* used for indexing prop (see prop_pt below) */
  int q_stride,          /* used for indexing prop (see prop_pt below) */
  int op_stride,         /* used for indexing prop (see prop_pt below) */
  int w_meson_store_t[], /* for compact storage */
  int w_meson_my_t[],    /* for compact storage */
  int w_meson_nstore,    /* for compact storage */
  int no_q_momenta,       /* number of q values */
  int q_momstore[][3],   /* q values themselves */
  int no_gamma_corr,     /* number of pairs */
  gamma_corr gamma_table[], /* table of gamma matrix pairs */
  field_offset tmp,      /* vector of complex numbers */
  int dimtmp             /* maximum number of complex values in tmp */
  )
{
  register int i;
  register site *s; 
  
  int store_t;
  int sf, si;
  int i_gamma_corr,q_pt,prop_pt;
  int old_gamma_out;
  
  double theta ; 
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 
  int px,py,pz;
  complex phase_fact ; 
  
  complex z;
  
  spin_wilson_vector localmat;  /* temporary storage */
  spin_wilson_vector antiquark;           /* temporary storage for antiquark */
  
  dirac_matrix *meson;
  dirac_matrix_v *meson_q;
  complex tr[MAXQ];

  /* performance */
  double dtime;
  int nmeson_evals = 0;
  int nprop_evals = 0;
  double flops,mflops;

  dtime = -dclock();

  /* Allocate temporary storage for meson propagator */

  if(no_q_momenta > MAXQ)
    {
      printf("meson_cont_mom_lean2: no_q_momenta %d exceeds max %d\n",
	     no_q_momenta,MAXQ);
      terminate(1);
    }
  
  meson = (dirac_matrix *)malloc(sites_on_node*sizeof(dirac_matrix));
  if(meson == NULL){
    printf("meson_cont_mom_lean2: No room for meson\n");
    terminate(1);
  }
  
  meson_q = (dirac_matrix_v *)malloc(w_meson_nstore*sizeof(dirac_matrix_v));
  if(meson_q == NULL){
    printf("meson_cont_mom_lean2: No room for meson_q\n");
    terminate(1);
  }
  
  /* Precompute phase factors for Fourier transform - save in tmp */
  if(no_q_momenta > dimtmp)
    {
      printf("meson_cont_mom_lean2: No room for FFT phases\n");
      terminate(1);
    }

  FORALLSITES(i,s) {
    for(q_pt=0; q_pt<no_q_momenta; q_pt++)
      {
	px = q_momstore[q_pt][0];
	py = q_momstore[q_pt][1];
	pz = q_momstore[q_pt][2];
	
	theta = factx*(s->x)*px + facty*(s->y)*py + factz*(s->z)*pz; 
	((complex *)F_PT(s,tmp))[q_pt] = 
	  cmplx((Real) cos(theta)  , (Real) sin(theta)) ; 
    }
  }      

  
  /* Run through the table of source-sink gamma matrices */

  for(i_gamma_corr = 0; i_gamma_corr < no_gamma_corr; i_gamma_corr++)
    {
      /* Skip reconstruction of "meson" if "out" gamma matrix hasn't changed */
      if(i_gamma_corr == 0 || gamma_table[i_gamma_corr].gout != old_gamma_out)
	{
	  old_gamma_out = gamma_table[i_gamma_corr].gout;

	  nmeson_evals++;    /* performance */

	  FORALLSITES(i,s) {
	    
	    /* antiquark = gamma_5 adjoint of quark propagator gamma_5 */
	    /* But we use a complex dot product below, so don't
	       take the c.c. here and do the transpose later 
	       - so just do a gamma5 transformation now */
	    
	    /* left multiply antiquark by source gamma matrices,
	       beginning with gamma_5 for quark -> antiquark */

	    mult_sw_by_gamma_l((spin_wilson_vector *)F_PT(s,src1),
			       &localmat, G5);    
	    
	    /* right dirac multiplication by gamma-5 
	       (finishing up antiquark) */
	    mult_sw_by_gamma_r( &localmat, &antiquark, G5);     
	    
	    /* left multiply src2 by Gamma_out.  Result in localmat */
	    mult_sw_by_gamma_l( (spin_wilson_vector *)F_PT(s,src2), 
				&localmat,
				gamma_table[i_gamma_corr].gout);

	    /* combine with src2 quark and sew together colors 
	       and source spins to make meson Dirac matrix 
	    "meson" is then [\gamma_5 S_1 \gamma_5]^\dagger \Gamma_out S_ 2 
	    and is a Dirac matrix */

	    meson[i] = mult_swv_an(&antiquark,&localmat);
	    
	  } /* end FORALLSITES */
	  
	  /* Do FT on "meson" for momentum projection - 
	     Result in meson_q */
	  FORALLSITES(i,s) {
	    store_t = w_meson_store_t[s->t];
	    for(si = 0; si < 4; si++)for(sf = 0; sf < 4; sf++)
	      for(q_pt=0; q_pt<no_q_momenta; q_pt++)
		{
		  meson_q[store_t].d[si].d[sf].e[q_pt].real = 0.;
		  meson_q[store_t].d[si].d[sf].e[q_pt].imag = 0;
		}
	  }
	  FORALLSITES(i,s) {
	    store_t = w_meson_store_t[s->t];
	    for(si = 0; si < 4; si++)for(sf = 0; sf < 4; sf++)
	      for(q_pt=0; q_pt<no_q_momenta; q_pt++)
		{
		  phase_fact = ((complex *)F_PT(s,tmp))[q_pt];
		  meson_q[store_t].d[si].d[sf].e[q_pt].real += 
		    meson[i].d[si].d[sf].real*phase_fact.real -  
		    meson[i].d[si].d[sf].imag*phase_fact.imag;
		  meson_q[store_t].d[si].d[sf].e[q_pt].imag += 
		    meson[i].d[si].d[sf].real*phase_fact.imag +  
		    meson[i].d[si].d[sf].imag*phase_fact.real;
		}
	  }
	} /* end if not same Gamma_out */
      
      nprop_evals++;   /* performance */
      
      /* Complete the propagator */
      
      for(store_t=0; store_t < w_meson_nstore; store_t++) {
	/* Do final Dirac trace for all sink momenta q */
	dirac_v_tr_gamma(&meson_q[store_t],
		       gamma_table[i_gamma_corr].gin,tr,no_q_momenta);
	/* Store values for all sink momenta */	
	for(q_pt=0; q_pt<no_q_momenta; q_pt++)
	  {
	    prop_pt = store_t + base_pt + q_pt * q_stride + 
	      i_gamma_corr * op_stride;
	    
	    prop[prop_pt ].real += tr[q_pt].real;
	    prop[prop_pt ].imag += tr[q_pt].imag;
	  }
      }
    }  /**** end of the loop over gamma table ******/
  
  free(meson);
  free(meson_q);

  dtime += dclock();

  flops = sites_on_node*(nmeson_evals*(1536 + 128*no_q_momenta) +
	       68*no_q_momenta) + 10*nprop_evals*w_meson_nstore*no_q_momenta;
  if(dtime > 0)mflops = flops/(dtime*1e6);
  else mflops = 0;

  node0_printf("meson_cont_mom_lean2: time %.1f sec %.1f MF\n",
	       dtime,mflops);fflush(stdout);
  
} /* end of meson_cont_mom_lean3 function  */


