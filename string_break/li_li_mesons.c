/************************** li_li_mesons.c *****************************/
/* MIMD version 6 */

/* Compute the light-light meson propagotors from the light
   quark propagators qprop[j].
   Here qprop[j] can either be a point-sink or smeared-sink
   quark propagator. */

#include "string_break_includes.h"

void li_li_mesons(int tot_smear, int t0, char *sink)
{
register int i;
register site *s; 
int x, y, z, t, t_off, j;
Real piprop, pi2prop, rhoprop, rho2prop;
Real *piprop_t, *pi2prop_t, *rhoprop_t, *rho2prop_t;
complex cc;

    piprop_t = (Real *)malloc(nt*sizeof(Real));
    pi2prop_t = (Real *)malloc(nt*sizeof(Real));
    rhoprop_t = (Real *)malloc(nt*sizeof(Real));
    rho2prop_t = (Real *)malloc(nt*sizeof(Real));
    for(t=0; t<nt; t++){
	piprop_t[t] = 0.0;
	pi2prop_t[t] = 0.0;
	rhoprop_t[t] = 0.0;
	rho2prop_t[t] = 0.0;
    }

    for(j=0; j<num_src; j++){

	for(t=0; t<nt; t++){
	    piprop = pi2prop = rhoprop = rho2prop = 0.0;

	    /* define the time value offset t from t0 */
	    t_off = (t+t0)%nt;

	    for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++){

		if( node_number(x,y,z,t_off) != mynode() )continue;
		i = node_index(x,y,z,t_off);
		cc = su3_dot( &lattice[i].qprop[j], &lattice[i].qprop[j]);

		piprop += cc.real;

		if( (x+y)%2==0)rhoprop += cc.real;
		else           rhoprop -= cc.real;
		if( (y+z)%2==0)rhoprop += cc.real;
		else           rhoprop -= cc.real;
		if( (z+x)%2==0)rhoprop += cc.real;
		else           rhoprop -= cc.real;

		if( x%2==0)rho2prop += cc.real;
		else       rho2prop -= cc.real;
		if( y%2==0)rho2prop += cc.real;
		else       rho2prop -= cc.real;
		if( z%2==0)rho2prop += cc.real;
		else       rho2prop -= cc.real;

		if( (x+y+z)%2==0)pi2prop += cc.real;
		else             pi2prop -= cc.real;
	    }

	    g_sync();

	    g_floatsum( &piprop );
	    g_floatsum( &pi2prop );
	    g_floatsum( &rhoprop );
	    g_floatsum( &rho2prop );
	    if(this_node == 0)
		printf("%s_MESONS_%d %d %d %e %e %e %e\n", sink, tot_smear,
		    j, t, (double)piprop, (double) rhoprop,
		    (double)pi2prop, (double) rho2prop);

	    piprop_t[t] += piprop;
	    pi2prop_t[t] += pi2prop;
	    rhoprop_t[t] += rhoprop;
	    rho2prop_t[t] += rho2prop;
	}

    } /* j < num_src */

    for(t=0; t<nt; t++){
	piprop_t[t] /= (Real)num_src;
	pi2prop_t[t] /= (Real)num_src;
	rhoprop_t[t] /= (Real)num_src;
	rho2prop_t[t] /= (Real)num_src;
	if(this_node == 0)
	    printf("%s_MESONS_T%d %d %e %e %e %e\n", sink, tot_smear, t,
		(double)piprop_t[t], (double)rhoprop_t[t],
		(double)pi2prop_t[t], (double)rho2prop_t[t]);
    }

    free(piprop_t);
    free(pi2prop_t);
    free(rhoprop_t);
    free(rho2prop_t);

} /* li_li_mesons */

