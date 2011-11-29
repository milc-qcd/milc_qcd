/*************** d_action.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
*/

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#include "symanzik_sl32_includes.h"

double d_action(){
double ssplaq,stplaq,g_action,h_action;

    d_plaquette(&ssplaq,&stplaq);
    g_action = -beta*volume*(ssplaq+stplaq);
    node0_printf("PLAQUETTE ACTION: %e\n", g_action);

    g_action = (beta/3.0)*imp_gauge_action();
    h_action = hmom_action();

    node0_printf("ACTION: g,h = %e  %e  %e\n",
	g_action, h_action, (g_action+h_action));

    return(g_action+h_action);
}

/* gauge momentum contribution to the action */
double hmom_action() {
register int i,dir;
register site *s;
double sum;
Real ahmat_mag_sq(anti_hermitmat *pt);

    sum= (double)0.0;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
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
