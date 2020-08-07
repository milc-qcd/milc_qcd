/*************** d_action.c ****************************************/
/* MIMD version 7 */

/* Measure total action, as needed by the hybrid Monte Carlo algorithm.
*/

/* Modifications:
   2/17/98  ANSI prototyping U.M.H.
   */

#include "ani_non_prod_tests_includes.h"

double d_action(){
double g_action,h_action;
#ifdef SCHROED_FUN
Real fs = beta*(Real)(nx*ny*nz*(nt-1)), ft = beta*(Real)(volume);
#endif

#ifndef ANISOTROPY
double ssplaq,stplaq;
    d_plaquette(&ssplaq,&stplaq);
#ifdef SCHROED_FUN
    g_action = -fs*ssplaq-ft*stplaq;
#else
    g_action = -beta*volume*(ssplaq+stplaq);
#endif
#else // #ifndef ANISOTROPY
int mu,nu,imunu;
double plaq[6];
    d_plaquette6(plaq);
#ifdef DEBUG
    node0_printf("Pmunu %.6e %.6e %.6e %.6e %.6e %.6e\n",plaq[0], plaq[1], plaq[2], plaq[3], plaq[4], plaq[5]);
#endif
    for ( g_action=0., imunu=0, mu = YUP; mu<=TUP; mu++ ) { 
      for ( nu = XUP; nu < mu; nu++, imunu++ ) { 
        int is_anisotropic = ( mu == ani_dir || nu == ani_dir ? 1 : 0 );
#ifdef SCHROED_FUN
        int volfac = ( mu < TUP ? fs : ft );
        g_action -= beta[is_anisotropic]*plaq[imunu]*volfac;
#else
        g_action -= beta[is_anisotropic]*plaq[imunu]*volume;
#endif
#ifdef DEBUG
        node0_printf("imunu %d is_a %d (mu=%d,nu=%d) bet=%.3f p=%.6f\n",imunu,is_anisotropic, mu,nu,beta[is_anisotropic], plaq[imunu]);
#endif
      }   
    }   
    g_action /= 3.;            
#endif
    node0_printf("PLAQUETTE ACTION: %e\n", g_action);

#ifndef ANISOTROPY
    g_action = (beta/3.0)*imp_gauge_action();
#else
/*  It would be favorable to keep the beta[0]/3. factor outside of 
    the function double imp_gauge_action() to retain the
    logic of the code. However, from the physics perspective it
    should go into the function, which is the logic we follow here 
    in the anisotropic case. JHW, 01/28/2019 */
    g_action = imp_gauge_action()/3.;
#endif
    h_action = hmom_action();

    node0_printf("ACTION: g,h = %.16e  %.16e  %.16e\n",
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
