/*************** d_action.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
   When this routine is called the conjugate gradient should already
   have been run on the even sites, so that the vector xxx contains
   (M_adjoint*M)^(-1) * phi.
*/

#include "schroed_ks_includes.h"

static double hmom_action();
static double fermion_action();
static Real ahmat_mag_sq(anti_hermitmat *pt);

double d_action(){
double ssplaq,stplaq,g_action,h_action,f_action;
#ifdef REWEIGH
double ds_deta,bd_plaq;
#endif
Real fs,ft;
    fs = beta*(Real)(nx*ny*nz*(nt-1));
    ft = beta*(Real)(nx*ny*nz*nt);
    d_plaquette(&ssplaq,&stplaq);
    ssplaq *= -1.0; stplaq *= -1.0;	/* KS phase factors change sign */
    g_action = -fs*ssplaq-ft*stplaq;
#ifdef REWEIGH
    if(bc_flag > 0){
	coupling(&ds_deta,&bd_plaq);
	ds_deta *= -1.0;	/* KS phases change sign */
	g_action += gamma_rv*ds_deta;
    }
#endif
    h_action = hmom_action();
    f_action = fermion_action();
    if(this_node==0)printf("ACTION: g,h,f = %e  %e  %e  %e\n",
	g_action, h_action, f_action, (g_action+h_action+f_action));
    return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
static double fermion_action() {
register int i;
register site *s;
register complex cc;
double sum;
    sum=0.0;
    FOREVENSITES(i,s) if(s->t > 0){
	/* phi is defined on even sites only */
        cc = su3_dot( &(s->phi), &(s->xxx) );
        sum += (double)cc.real;
    }
    g_doublesum( &sum );
    return(sum);
}

/* gauge momentum contribution to the action */
static double hmom_action() {
register int i,dir;
register site *s;
double sum;
Real ahmat_mag_sq(anti_hermitmat *pt);

    sum=0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++) if(dir==TUP || s->t>0){
            sum += (double)ahmat_mag_sq( &(s->mom[dir]) );
	}
    }
    g_doublesum( &sum );
    return(sum);
}

/* magnitude squared of an antihermition matrix */
static Real ahmat_mag_sq(anti_hermitmat *pt) {
register Real x,sum;
    x = pt->m00im; sum  = 0.5*x*x;
    x = pt->m11im; sum += 0.5*x*x;
    x = pt->m22im; sum += 0.5*x*x;
    x = pt->m01.real; sum += x*x;
    x = pt->m01.imag; sum += x*x;
    x = pt->m02.real; sum += x*x;
    x = pt->m02.imag; sum += x*x;
    x = pt->m12.real; sum += x*x;
    x = pt->m12.imag; sum += x*x;
    return(sum);
}
