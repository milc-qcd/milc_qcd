/******** spectrum_s2.c *************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE.  Consider using spectrum2.c */
/* For measuring propagation IN THE Z DIRECTION (screening spectrum) */
/* ALSO measures the nonlocal propagators for the current algebra quark mass */

/* Spectrum for Kogut-Susskind pointlike hadrons, WALL SOURCE */
/* Wall source is modulated by a cosine at the lowest Matsubara frequency */
#include "ks_dyn_includes.h"

int spectrum() /* return the C.G. iteration number */
{
  Real piprop,pi2prop,rhoprop0,rhoprop1,rho2prop0,rho2prop1,barprop;
  complex *even_forw, *even_back, *odd_forw, *odd_back;
  Real mass_x2;
  register complex cc,cc2;
  Real finalrsq, th;
  register int i,x,y,z,t,icol,cgn;
  site *st;
  msg_tag *tag0,*tag1;

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
      
      /* initialize phi and xxx */
      clear_latvec( F_OFFSET(phi), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(t=0;t<nt;t+=2)
/**for(x=0;x<1;x+=2)for(y=0;y<1;y+=2)for(t=0;t<1;t+=2)**/
	{
	  if( node_number(x,y,0,t) != mynode() )continue;
	  i=node_index(x,y,0,t);
	  /* Modulate source with Matsubara phase */
	  lattice[i].phi.c[icol].real = -cos((double)t*th);
	}

      /* do a C.G. (source in phi, result in xxx) */
      load_ferm_links(&fn_links);
      cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			niter, rsqprop, PRECISION, EVEN, &finalrsq,
			&fn_links);
      /* Multiply by -Madjoint */
      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links);
      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), EVEN);
      
      /* fill the hadron matrix */
      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);
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


  /* measure the nonlocal propagators for the current algebra quark mass */
  /* parallel transport the propagator forward and backward in time.  The
     phase factor embedded in the link matrix is just the one we
     require. */
  even_forw = (complex *)malloc( nz*sizeof(complex) );
  even_back = (complex *)malloc( nz*sizeof(complex) );
  odd_forw = (complex *)malloc( nz*sizeof(complex) );
  odd_back = (complex *)malloc( nz*sizeof(complex) );
  for(z=0;z<nz;z++){
     even_forw[z] = even_back[z] = odd_forw[z] = odd_back[z] = cmplx(0.0,0.0);
  }
  for(icol=0; icol<3;icol++){
     tag0 = start_gather_site( F_OFFSET(propmat[icol]), sizeof(su3_vector), ZUP,
	EVENANDODD, gen_pt[0] );
     FORALLSITES(i,st){
	    mult_adj_su3_mat_vec( &(st->link[ZUP]),
                &(st->propmat[icol]),  &(st->tempvec[icol]) );
     }
     tag1 = start_gather_site( F_OFFSET(tempvec[icol]), sizeof(su3_vector),
		OPP_DIR(ZUP), EVENANDODD, gen_pt[1] );
     wait_gather(tag0);
     FORALLSITES(i,st){
	    mult_su3_mat_vec( &(st->link[ZUP]),
                (su3_vector *)gen_pt[0][i],  &(st->ttt) );
     }
     /* now ttt contains the propagator parallel transported from
	displacement +z_hat, and *gen_pt[1] the propagator parallel
	transported from -z_hat. */
     cleanup_gather(tag0);
     wait_gather(tag1);

     /* for each z, put the "forward" and "backward" contributions
	for even and odd sites separately. Construct the propagator
	summed over hypercubes in the analysis code. */
     FORALLSITES(i,st){
	cc  = su3_dot( &(st->propmat[icol]), &(st->ttt) );
	cc2 = su3_dot( &(st->propmat[icol]), (su3_vector *)gen_pt[1][i] );
	if( st->parity==EVEN ){
	    CSUM( even_forw[ st->z ], cc );
	    CSUM( even_back[ st->z ], cc2);
	}
	else {
	    CSUM( odd_forw[ st->z ], cc );
	    CSUM( odd_back[ st->z ], cc2);
	}
     } /* loop over sites */
     cleanup_gather(tag1);
  } /* loop on colors */

  for(z=0; z<nz; z++){
      g_complexsum( &even_forw[z] );
      g_complexsum( &even_back[z] );
      g_complexsum( &odd_forw[z] );
      g_complexsum( &odd_back[z] );
      /**
      if(mynode()==0)printf("XXXXXXX  %d  %e  %e  %e  %e\n",z,
	    (double)even_forw[z].real,(double)even_back[z].real,
	    (double)odd_forw[z].real,(double)odd_back[z].real );
      if(mynode()==0)printf("YYYYYYY  %d  %e  %e  %e  %e\n",z,
	    (double)even_forw[z].imag,(double)even_back[z].imag,
	    (double)odd_forw[z].imag,(double)odd_back[z].imag );
      **/
    } /* nz-loop */
   /* Now compute the NLT propagator */
   for(z=0; z<nz-1; z++){
      if(mynode()==0)printf("QM_PROP  %d  %e\n",z,
	    (double)( even_forw[z].real + odd_forw[z].real
	    + even_back[z+1].real + odd_back[z+1].real) );
   }
  
  return(cgn);
} /* spectrum */
