/******* mrilu_w_or.c - MR-ILU for  Wilson fermions ****/
/* MIMD version 6 */

/* Modifications:
   4/26/98 Moved parameters to structures CD
   8/29/97 ANSI prototyping and added comments C. D.
   Code created by ??
   */

/* double precision accumlulations */
/* Memory stingy version
  "r" overwrites src on even sites
  "p" overwrites src on odd sites
   3/29/00 EVENFIRST is the rule now. CD.
*/

/* Requires qic wilson vector temporary wv2 */

/* The source vector is in "src", and the initial guess and answer
   in "dest".  "r" is the residual vector, which is a pointer to src since
the source  is overwritten to save space 
 and "p" and "mmp" are
   working vectors for the minimal residue. 
   MinMR = minimum number of iterations.
   MaxMR = maximum number of iterations.
   size_r = desired residual, quit when we reach it.
(Square root def for residue size_r = sqrt(r*r))
   ILU resides on parity=EVEN so do only even sites
*/

#include "generic_wilson_includes.h"

/**#define MRTIME **/  /* Define for timing */

/*

  Work space in site structure specified in qic:
    wilson_vector wv2

*/


int mrilu_w_or(          /* Return value is number of iterations taken */
    field_offset src,    /* type wilson_vector (source vector - OVERWRITTEN!)*/
    field_offset dest,   /* type wilson_vector (answer and initial guess )*/
    quark_invert_control *qic, /* parameters controlling inversion */
    void *dmp            /* parameters defining the Dirac matrix */
    )
{
  /* Unpack required members of the structures */
  int MinMR = qic->min;      /* minimum number of iterations */
  int MaxMR = qic->max;      /* maximum number of iterations */
  Real RsdMR = qic->resid;  /* desired residual - 
				 normalized as sqrt(r*r)/sqrt(src_e*src_e */
  int flag = qic->start_flag;   /* 0: use a zero initial guess; 1: use dest */
  field_offset mmp = qic->wv2;    /* size of wilson_vector */

  dirac_wilson_param *dwp 
    = (dirac_wilson_param *)dmp; /* Cast pass-through pointer */
  Real Kappa = dwp->Kappa;     /* hopping */
  /* End of unpacking required members of structures */

    int N_iter;
    register int i;
    register site *s;
    Real magsq_wvec();
    Real size_src;
    double dsize_r,dsize_src;
    register complex b;
    complex a,na;
    register Real MKsq = -Kappa*Kappa;
    register field_offset r,p;
    double c[3];
    double_complex d;
#ifdef MRTIME
    double dtime;
#endif

    if(even_sites_on_node!=odd_sites_on_node){
        printf("Need same number of even and odd sites on each node\n");
        terminate(1);
    }
    r = src;
    p = src + even_sites_on_node*sizeof(site);
    /* This disgusting trick makes p for each even site actually be
        src on some corresponding odd site */
/**if(this_node==0)printf("MRILU: p=%d\n",p);**/

    /* MR_ILU: */

    /* Start Inversion */

    /* src = U^(-1)*src */
#ifdef MRTIME
    dtime = -dclock();
#endif

    dslash_w(src,mmp,PLUS,EVEN);

    /* Normalisation  */
    dsize_src=0.0;
    FOREVENSITES(i,s) {
        scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,src)), 
			      (wilson_vector *)F_PT(s,mmp),
            Kappa, ((wilson_vector *)F_PT(s,src)) );
        dsize_src += magsq_wvec( ((wilson_vector *)F_PT(s,src)) );
    }
    g_doublesum(&dsize_src);
    size_src = (Real)sqrt(dsize_src);

    /**if(this_node==0)printf("beginning inversion--size_src=%e\n",
        (double)size_src);**/

    /* r overwrites the source (bad for
        dynamical fermions but good for quenched calcs) */
    /* set r = src --- nothing to do */

    /* Initial guess */
    /* set dest = src on odd sites, whatever else you do
(even if we restart with a nonzero solution vector, the end of the
subroutine rebuilds the odd component from the even one. The (trivial)
solution of the odd component of the equation is dest = src, before
we rotate back to the basis in which  M is not checkerboard-diagonal) */
    FORODDSITES(i,s) {
        copy_wvec( ((wilson_vector *)F_PT(s,src)),
                   ((wilson_vector *)F_PT(s,dest)));
    }


    /* code if you want to start with dest=0... */
    if(flag == 0) {
        /**if(this_node==0)printf("dest_0=0\n");**/
        FOREVENSITES(i,s) {
            clear_wvec( ((wilson_vector *)F_PT(s,dest)) );
        }
        dsize_r = 1.00;
        qic->size_r = dsize_r;
        /**if(this_node==0)printf("size_r=%e\n",(double)qic->size_r));**/

    }
    /* code if you want to start dest with some particular starting value... */
    /* r=src[1]-[U^(-1)*M*L^(-1)]*dest */
    if(flag != 0) {
        /**if(this_node==0)printf("dest_0  !=0\n");**/
        /* we use mmp temporarily to construct r */
        dslash_w(dest,mmp,PLUS,ODD);
        dslash_w(mmp,mmp,PLUS,EVEN);
        FOREVENSITES(i,s) {
            scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,dest)),
                (wilson_vector *)F_PT(s,mmp),MKsq, 
				  (wilson_vector *)F_PT(s,mmp) );
            scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,r)),
                (wilson_vector *)F_PT(s,mmp),-1.0, 
				  ((wilson_vector *)F_PT(s,r)) );
        }

        dsize_r=0.0;
        FOREVENSITES(i,s) {
            dsize_r += magsq_wvec( ((wilson_vector *)F_PT(s,r)) );
        }
        g_doublesum(&dsize_r);
        qic->size_r = (Real)sqrt(dsize_r)/size_src;
        /**if(this_node==0)    printf("beginning inversion--size_r=%e\n",
            (double)(qic->size_r));**/

    }

    for( N_iter = 0; N_iter < MinMR || (N_iter < MaxMR && RsdMR  < qic->size_r); 
        N_iter = N_iter + 1) {

    	/*  mmp = [U^(-1)*M*L^(-1)]*r  */
    	dslash_w(r,mmp,PLUS,ODD);
    	dslash_w(mmp,mmp,PLUS,EVEN);

        /* c = |mmp|^2,  d = (mmp,r)  */
    	c[0]=0.0;
    	d=dcmplx(0.0,0.0);
    	FOREVENSITES(i,s) {
        scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,r)),
			      (wilson_vector *)F_PT(s,mmp),MKsq, 
			      (wilson_vector *)F_PT(s,mmp) );
        c[0] += magsq_wvec( (wilson_vector *)F_PT(s,mmp) );
        b = wvec2_dot( ((wilson_vector *)F_PT(s,mmp)), 
		       ((wilson_vector *)F_PT(s,r)) );
        CSUM(d,b);
    	}
	c[1] = d.real;
	c[2] = d.imag;
    	g_vecdoublesum(c,3);

    	a.real = (1.3) * (c[1]/c[0]);
	a.imag = (1.3) * (c[2]/c[0]);

        if( fabs((double)(a.real)) < 0.3 && fabs((double)(a.imag)) < 0.3 ) {
		a.real = 0.3;
		a.imag = 0.3;
	}


	/**if(this_node==0) printf("a-real =%e,a-imag =%e\n", a.real,a.imag);**/


	CNEGATE(a,na);

   
        /* dest = dest + a*r  */
        /* r = r - a*mmp */
 	dsize_r=0.0;
        FOREVENSITES(i,s) {
            c_scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,dest)),
                (wilson_vector *)F_PT(s,r), &a,
                ((wilson_vector *)F_PT(s,dest)) );
            c_scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,r)),
                (wilson_vector *)F_PT(s,mmp),&na, 
                ((wilson_vector *)F_PT(s,r)) );
            dsize_r += magsq_wvec( ((wilson_vector *)F_PT(s,r)) );
        }

        g_doublesum(&dsize_r);

        qic->size_r = (Real)sqrt(dsize_r)/size_src;
        /**if(this_node==0)printf("iteration= %d, residue= %e\n",N_iter,
	   (double)(qic->size_r));**/
 
    }

    /*** End of mrilu iterations ***/

#ifdef MRTIME
dtime += dclock();
#endif
if(this_node==0){
  if( N_iter == 0 )
    printf("MRILU_OR: NO iterations taken\n") ;
#ifdef MRTIME
  else
    printf("MRILU_OR: time = %.2e size_r= %.2e iters= %d MF = %.1f\n",
	   (double)3170*N_iter*even_sites_on_node/(dtime*1e6));
#endif
  fflush(stdout);
}


/** if( (qic->size_r) > RsdMR ) {
        if(this_node==0)printf(" MR_ILU Not Converged\n");
	} **/

    /* dest = L^(-1)*dest  */
    dslash_w(dest,mmp,PLUS,ODD);
    FORODDSITES(i,s) {
        scalar_mult_add_wvec( ((wilson_vector *)F_PT(s,dest)), 
			      (wilson_vector *)F_PT(s,mmp), Kappa, 
			      ((wilson_vector *)F_PT(s,dest)) );
    }

    return(N_iter);
}
/* mrilu_w_or.c */
