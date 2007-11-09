/******** spectrum.c *************/
/* MIMD version 7 with SG additions*/
/* Last updated 2/21/95, by SG. */
/* 12/11/95 DT change looping over source time slices */
/* 7/29/96 DT fixed names of rho and rho2 propagators.
   rho is VT (rho-b_1), rho2 is PV (rho-a_1)
   also pi and pi2 -> pi_ps_prop, pi_sc_prop */
/* Gaugefixing tolerance via macro now.  Calls ks_congrad. CD */
/* For measuring propagation in the t direction */

/* Spectrum for Kogut-Susskind pointlike hadrons, wall source, point
  source if t_source >= nt */

#include "generic_ks_includes.h"

/* return the C.G. iteration number */
int spectrum(ferm_links_t *fn) 
{
  Real pi_ps_prop,pi_sc_prop,rho_pv_prop,rho_vt_prop,barprop;
  Real mass_x2;
  register complex cc;
  Real finalrsq;
  int isrc;
  register int i,x,y,z,t,icol,cgn;
  register int t_source,t_off;
  /* Preset values for now */
  int source_start = 0;
  int source_inc = nt/2;
  int n_sources = 2;
  
  /* Fix TUP Coulomb gauge - gauge links only*/
  
  rephase( OFF );
  gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  rephase( ON );
  
  mass_x2 = 2.*mass;
  cgn=0;
  for(t_source=source_start, isrc=0; t_source<2*nt && isrc < n_sources;
       ++isrc, t_source += source_inc ) {

	  for(icol=0; icol<3; icol++) {
		      
	      /* initialize ttt and xxx */
	      clear_latvec( F_OFFSET(phi), EVENANDODD);
	      clear_latvec( F_OFFSET(xxx), EVENANDODD);

	      if(t_source < nt)  /* wall source */
	         for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
		     if( node_number(x,y,z,t_source) != mynode() )continue;
		     i=node_index(x,y,z,t_source);
		     lattice[i].phi.c[icol].real = -1.0;
		   }
	      else{  /* point source at origin */
		     if( node_number(0,0,0,t_source%nt) == mynode() )
			{ 
			   i=node_index(0,0,0,t_source%nt);
		     	   lattice[i].phi.c[icol].real = (-nx*ny*nz/8);
			}
		}
		

	      /* do a C.G. (source in phi, result in xxx) */
	      cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				niter, nrestart, rsqprop, PRECISION, 
				EVEN, &finalrsq, fn);
	      /* Multiply by -Madjoint */
	      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, fn);
	      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), EVEN);

	      
	      /* fill the hadron matrix */
	      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);
	    } /* end loop on icol */
	  
	  /* measure the meson propagator */
	  for(t=0; t<nt; t++)
	    {
	      /* clear meson propgators */
	      pi_ps_prop=rho_pv_prop=pi_sc_prop=rho_vt_prop=0.;
	      /* define the time value offset t from t_source */
	      t_off = (t+t_source)%nt;
	      
	      for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)
		for(icol=0;icol<3;icol++)
		  {
		    if( node_number(x,y,z,t_off) != mynode() )continue;
		    i=node_index(x,y,z,t_off);
		    cc = su3_dot( &lattice[i].propmat[icol],
				 &lattice[i].propmat[icol] );
		    
		    pi_ps_prop += cc.real;
		    
		    if( (x+y)%2==0)rho_pv_prop += cc.real;
		    else	   rho_pv_prop -= cc.real;
		    if( (y+z)%2==0)rho_pv_prop += cc.real;
		    else	   rho_pv_prop -= cc.real;
		    if( (z+x)%2==0)rho_pv_prop += cc.real;
		    else	   rho_pv_prop -= cc.real;
		    
		    if( x%2==0)rho_vt_prop += cc.real;
		    else       rho_vt_prop -= cc.real;
		    if( y%2==0)rho_vt_prop += cc.real;
		    else       rho_vt_prop -= cc.real;
		    if( z%2==0)rho_vt_prop += cc.real;
		    else       rho_vt_prop -= cc.real;
		    
		    if( (x+y+z)%2==0)pi_sc_prop += cc.real;
		    else	     pi_sc_prop -= cc.real;
		    
		  }
	      
	      g_sync();
	      
	      /* dump meson propagators */
	      g_floatsum( &pi_ps_prop );
	      g_floatsum( &rho_pv_prop );
	      g_floatsum( &rho_vt_prop );
	      g_floatsum( &pi_sc_prop );
	      if(mynode()==0)printf("MES_PRO%d  %d  %e  %e  %e  %e\n",t_source,
			t,(double)pi_ps_prop,(double)rho_pv_prop,
			(double)pi_sc_prop,(double)rho_vt_prop);
	    } /* nt-loop */
		fflush(stdout);
	  
	  /* measure the baryon propagator */
	  for(t=0; t<nt; t++)
	    {
	      /* clear baryon propgators */
	      barprop=0.0;
	      /* define the time value offset t from t_source */
	      t_off = (t+t_source)%nt;
	      
	      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2)
		{
		  if( node_number(x,y,z,t_off) != mynode() )continue;
		  i=node_index(x,y,z,t_off);
		  cc = det_su3( (su3_matrix *)lattice[i].propmat );
		  barprop += cc.real;
		}
		      
	      g_sync();
	      
	      /* dump baryon propagators */
	      g_floatsum( &barprop );
	      if(mynode()==0) {
		  /* must get sign right.  This looks to see if we have
			wrapped around the lattice.  "t" is the distance
			from the source to the measurement, so we are
			trying to find out if t_source+t is greater than
			or equal to nt.  the "-tsource/nt" is in there
			so that it will work correctly with tsource=nt,
			which we use for the point source and which is
			logically the same as t_source=0. */
		  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
		     printf("NUC_PRO%d  %d  %e\n",t_source,t,(double)barprop);
		  else  /* change sign because antiperiodic b.c.  sink point
			should really be in a copy of the lattice */
		     printf("NUC_PRO%d  %d  %e\n",t_source,t,-(double)barprop);
	      }
	    } /* nt-loop */
		fflush(stdout);
  
    } /* end loop on t_source */
  return(cgn);
} /* spectrum */

