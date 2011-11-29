/*************** action.c ****************************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE. Consider using d_action.c */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
   When this routine is called the conjugate gradient should already
   have been run on the even sites, so that the vector xxx contains
   (M_adjoint*M)^(-1) * phi.
*/

#include "wi_dyn_includes.h"

Real action(){
Real hmom_action(),fermion_action();
 double ssplaq,stplaq;
 Real g_action,h_action,f_action;
    d_plaquette(&ssplaq,&stplaq);
    g_action = -beta*volume*(ssplaq+stplaq);
    h_action = hmom_action();
    f_action = fermion_action();
/**if(this_node==0)printf("ACTION: g,h,f = %e  %e  %e  %e\n",
(double)g_action,(double)h_action,(double)f_action,
(double)(g_action+h_action+f_action));**/
    return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
Real fermion_action() {
register int i;
register site *s;
register complex cc;
Real sum;
    sum=0.0;
#ifndef LU
    FORALLSITES(i,s){
#else
    FOREVENSITES(i,s){
#endif
        sum += wvec_rdot( &(s->psi), &(s->chi) );
    }
    g_floatsum( &sum );
    return(sum);
}

/* gauge momentum contribution to the action */
Real hmom_action() {
register int i,dir;
register site *s;
Real sum;
Real ahmat_mag_sq();

    sum=0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
            sum += ahmat_mag_sq( &(s->mom[dir]) );
	}
    }
    g_floatsum( &sum );
    return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(pt) anti_hermitmat *pt; {
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
