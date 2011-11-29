/******** spectrum_s.c *************/
/* MIMD version 6 */
/* For measuring propagation IN THE Z DIRECTION (screening spectrum) */

/* Spectrum for Kogut-Susskind pointlike hadrons, wall source */
/* Wall source is modulated by a cosine at the lowest Matsubara frequency */
#include "ks_dyn_includes.h"

int spectrum() /* return the C.G. iteration number */
{
  Real piprop,pi2prop,rhoprop0,rhoprop1,rho2prop0,rho2prop1,barprop;
  Real mass_x2;
  register complex cc;
  Real finalrsq, th;
  register int i,x,y,z,t,icol,cgn;

  /* Fix ZUP Coulomb gauge - gauge links only*/
  rephase( OFF );
  gaugefix(ZUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  rephase( ON );
  
  mass_x2 = 2.*mass;
  cgn=0;
  /* Phase increment for minimum Matsubara frequency */
  th = PI/nt;

  for(icol=0; icol<3; icol++)
    {
      
      /* initialize ttt and xxx */
      clear_latvec( F_OFFSET(ttt), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(t=0;t<nt;t+=2)
/**for(x=0;x<1;x+=2)for(y=0;y<1;y+=2)for(t=0;t<1;t+=2)**/
	{
	  if( node_number(x,y,0,t) != mynode() )continue;
	  i=node_index(x,y,0,t);
	  /* Modulate source with Matsubara phase */
	  lattice[i].ttt.c[icol].real = -cos((double)t*th);
	}

      /* Multiply by -Madjoint */
      load_ferm_links(&fn_links);
      dslash_site( F_OFFSET(ttt), F_OFFSET(phi), ODD, &fn_links);
      scalar_mult_latvec( F_OFFSET(ttt), -mass_x2, F_OFFSET(phi), EVEN);
      /* do a C.G. */
      load_ferm_links(&fn_links);
      cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			niter, rsqprop, PRECISION, EVENANDODD, &finalrsq,
			&fn_links);
      
      /* fill the hadron matrix */
      copy_latvec( F_OFFSET(xxx), F_OFFSET(propmat[icol]), EVENANDODD);
    } /* end loop on icol */
  
  /* measure the meson propagator */
  for(z=0; z<nz; z++)
    {
      /* clear meson propgators */
      piprop=rhoprop0=rhoprop1=pi2prop=rho2prop0=rho2prop1=0.;
      
      for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(t=0;t<nt;t++)
	for(icol=0;icol<3;icol++)
	  {
	    if( node_number(x,y,z,t) != mynode() )continue;
	    i=node_index(x,y,z,t);
	    cc = su3_dot( &lattice[i].propmat[icol],
			 &lattice[i].propmat[icol] );
	    
	    piprop += cc.real;
	    
	    if( (x+y)%2==0)rhoprop0 += cc.real;
	    else	   rhoprop0 -= cc.real;
	    if( (y+t)%2==0)rhoprop1 += cc.real;
	    else	   rhoprop1 -= cc.real;
	    if( (t+x)%2==0)rhoprop1 += cc.real;
	    else	   rhoprop1 -= cc.real;
	    
	    if( x%2==0)rho2prop1 += cc.real;
	    else       rho2prop1 -= cc.real;
	    if( y%2==0)rho2prop1 += cc.real;
	    else       rho2prop1 -= cc.real;
	    if( t%2==0)rho2prop0 += cc.real;
	    else       rho2prop0 -= cc.real;
	    
	    if( (x+y+t)%2==0)pi2prop += cc.real;
	    else	     pi2prop -= cc.real;
	    
	  }
      
      g_sync();
      
      /* dump meson propagators */
      g_floatsum( &piprop );
      g_floatsum( &rhoprop0 );
      g_floatsum( &rhoprop1 );
      g_floatsum( &rho2prop0 );
      g_floatsum( &rho2prop1 );
      g_floatsum( &pi2prop );
      if(mynode()==0)printf("MES_SCREEN  %d  %e  %e  %e  %e  %e  %e\n",z,
			    (double)piprop,(double)rhoprop0,
			    (double)rhoprop1,
			    (double)pi2prop,(double)rho2prop0,
                            (double)rho2prop1);
    } /* nz-loop */
  
  /* measure the baryon propagator */
  for(z=0; z<nz; z++)
    {
      /* clear baryon propgators */
      barprop=0.0;
      
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(t=0;t<nt;t+=2)
	{
	  if( node_number(x,y,z,t) != mynode() )continue;
	  i=node_index(x,y,z,t);
	  cc = det_su3((su3_matrix *)lattice[i].propmat );
	  /* Include phase for Fourier component at minimum Matsubara freq.*/
	  barprop += cc.real*cos((double)t*th);
	}
      
      g_sync();
      
      /* dump baryon propagators */
      g_floatsum( &barprop );
      if(mynode()==0)printf("NUC_SCREEN  %d  %e\n",z,(double)barprop);
    } /* nz-loop */
  
  return(cgn);
} /* spectrum */
