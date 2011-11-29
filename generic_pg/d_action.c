/*************** d_action.c ****************************************/
/* MIMD version 6 */
/* UMH: Combined with Schroedinger functional version, Jan 2000 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
*/

#include "generic_pg_includes.h"

double d_action(){
double ssplaq,stplaq,g_action,h_action;
double hmom_action();
#ifdef SCHROED_FUN
Real fs,ft;
#endif
    d_plaquette(&ssplaq,&stplaq);
#ifdef SCHROED_FUN
    fs = beta*(Real)(nx*ny*nz*(nt-1));
    ft = beta*(Real)(volume);
    g_action = -fs*ssplaq-ft*stplaq;
#else
    g_action = -beta*volume*(ssplaq+stplaq);
#endif
    h_action = hmom_action();
    if(this_node==0)printf("ACTION: g,h = %e  %e  %e\n",
	g_action,h_action, (g_action+h_action));
    return(g_action+h_action);
}

/* gauge momentum contribution to the action */
double hmom_action() {
register int i,dir;
register site *s;
double sum;
Real ahmat_mag_sq(anti_hermitmat *pt);

    sum=0.0;
    FORALLSITES(i,s){
#ifdef SCHROED_FUN
	for(dir=XUP;dir<=TUP;dir++) if(dir==TUP || s->t>0){
#else
	for(dir=XUP;dir<=TUP;dir++){
#endif
            sum += (double)ahmat_mag_sq( &(s->mom[dir]) );
	}
    }
    g_doublesum( &sum );
    return(sum);
}

/* magnitude squared of an antihermition matrix */
Real ahmat_mag_sq(anti_hermitmat *pt) {
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
