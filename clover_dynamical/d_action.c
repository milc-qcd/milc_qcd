/*************** d_action.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
   When this routine is called the conjugate gradient should already
   have been run on the even sites, so that the vector psi contains
   (M_adjoint*M)^(-1) * chi.
*/

#include "cl_dyn_includes.h"
Real ahmat_mag_sq(anti_hermitmat *pt);

double d_action(){
double d_hmom_action(),d_fermion_action();
double ssplaq,g_action,h_action,f_action;

    gauge_action(&ssplaq);
/* printf("in d_action, gauge_action=%e\n",&ssplaq); */

    g_action = beta*ssplaq/3.0;
    h_action = d_hmom_action();
    f_action = d_fermion_action();
/**if(this_node==0)printf("D_ACTION: g,h,f = %e  %e  %e  %e\n",
g_action,h_action,f_action,
(g_action+h_action+f_action));**/
    return(g_action+h_action+f_action);
}

/* fermion contribution to the action */
double d_fermion_action() {
register int i;
register site *s;
double sum;
    sum= (double)0.0;
#ifndef LU
    FORALLSITES(i,s){
#else
    FOREVENSITES(i,s){
#endif
        sum += (double)wvec_rdot( &(s->psi), &(s->chi) );
    }
    g_doublesum( &sum );
/**{Real xxx ; xxx = sum; g_floatsum( &xxx ); sum = xxx;}**/
    return(sum);
}

/* gauge momentum contribution to the action */
double d_hmom_action() {
register int i,dir;
register site *s;
double sum;

    sum= (double)0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
            sum += (double)ahmat_mag_sq( &(s->mom[dir]) );
	}
    }
    g_doublesum( &sum );
/**{Real xxx ; xxx = sum; g_floatsum( &xxx ); sum = xxx;}**/
    return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt){
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

