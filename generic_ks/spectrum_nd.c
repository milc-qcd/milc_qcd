/******** spectrum_nd.c *************/
/* MILC version 7 */
/* DT 9/20/00 started, from spectrum_nlpi2.c
   KS hadron masses with nondegenerate quarks ("nd")
   compute corner source propagators with two masses, assemble
   local mesons and nucleons.

   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie phases_in=ON )

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop_SINK_[lh] */
/* Name by KS flavor content e.g. Goldstone pion = pion5 */
/* operators:
   (all mesons = heavy+light)
   pion5_ll:	local 0-+:  (flavor)gamma_5     partner=0+-  phase=(1), two light quarks
   pion5_lh:	local 0-+:  (flavor)gamma_5     partner=0+-  phase=(1)  one heavy, one light
   pion5_hh:	local 0-+:  (flavor)gamma_5     partner=0+-  phase=(1)  two heavy
   pion05_ll:	local 0-+:  gamma_0 gamma_5     partner=0++  phase=(-1)^(x+y+z+t)
   pion05_lh:	local 0-+:  gamma_0 gamma_5     partner=0++  phase=(-1)^(x+y+z+t)
   pion05_hh:	local 0-+:  gamma_0 gamma_5     partner=0++  phase=(-1)^(x+y+z+t)
   rhoi_ll:	local 1--: gamma_i              partner=1+-  phase=(-1)^(dir) (VT)
   rhoi_lh:	local 1--: gamma_i              partner=1+-  phase=(-1)^(dir) (VT)
   rhoi_hh:	local 1--: gamma_i              partner=1+-  phase=(-1)^(dir) (VT)
   rhoi0_ll:	local 1--: gamma_i gamma_0      partner=1++  phase=(-1)^(x+y+z+t+dir) (PV)
   rhoi0_lh:	local 1--: gamma_i gamma_0      partner=1++  phase=(-1)^(x+y+z+t+dir) (PV)
   rhoi0_hh:	local 1--: gamma_i gamma_0      partner=1++  phase=(-1)^(x+y+z+t+dir) (PV)
   nuc_lll:	local nucleon, three light,
   nuc_llh:	local nucleon, two light, one heavy: 
   nuc_lhh:	local nucleon, one light, two heavy: 
   nuc_hhh:	local nucleon, three heavy: 
*/
enum prop_name { 
    prop_pion5_ll,
    prop_pion5_lh,
    prop_pion5_hh,
    prop_pion05_ll,
    prop_pion05_lh,
    prop_pion05_hh,
    prop_rhoi_ll,
    prop_rhoi_lh,
    prop_rhoi_hh,
    prop_rhoi0_ll,
    prop_rhoi0_lh,
    prop_rhoi0_hh,
    prop_nuc_lll,
    prop_nuc_llh,
    prop_nuc_lhh,
    prop_nuc_hhh,
    nprops		/* nprops = number of propagators */
};

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"

int test_converge( int t_source );
su3_vector * lightprop[3];
su3_vector * heavyprop[3];

int spectrum_nd( Real mass1, Real mass2, Real tol, ferm_links_t *fn ){
  /* arguments are light and heavy quark masses, return C.G. iteration number */

  int cgn;
  register int i,j,x,y,z,t,t_off;
  register site* s;
  register complex cc;
  register int t_source;
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  complex **props;	/* arrays of propagators */
  su3_matrix tmat;
  Real finalrsq;
  
  cgn=0; /* number of CG iterations */

  /* allocate arrays to accumulate propagators */
  props = (complex **)malloc(nprops*sizeof(complex *));
  props[0] = (complex *)malloc(nprops*nt*sizeof(complex));
  for(i=1;i<nprops;i++)props[i]=props[i-1]+nt;

  /* set propagators to zero */
  for(cc.real=cc.imag=0.0,i=0;i<nprops;i++)for(j=0;j<nt;j++){
    props[i][j]=cc;
  }

  /* allocate light and heavy quark propagators for each color */
  for( color=0; color<3; color++){
    lightprop[color] = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
    heavyprop[color] = (su3_vector *)malloc( sizeof(su3_vector)*sites_on_node );
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Corner wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_nd(): source time = %d\n",t_source);

    for(color=0;color<3;color++){
      clear_latvec( F_OFFSET(quark_source), EVENANDODD );
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
        if( node_number(x,y,z,t_source) != mynode() )continue;
        i=node_index(x,y,z,t_source);
        lattice[i].quark_source.c[color].real = -1.0;
      }

      /* do a C.G. (source in quark_source, result in g_rand) */
      if(t_source%2 == 0) {
         cgn += ks_congrad( F_OFFSET(quark_source), F_OFFSET(g_rand),
			    mass1, niter, nrestart, rsqprop, PRECISION, 
			    EVEN, &finalrsq, fn);
         /* Multiply by -Madjoint */
         dslash_site( F_OFFSET(g_rand), F_OFFSET(quark_prop), ODD, fn);
         scalar_mult_latvec( F_OFFSET(g_rand), -2.0*mass1, F_OFFSET(quark_prop),EVEN);
      }
      else {
        cgn += ks_congrad( F_OFFSET(quark_source), F_OFFSET(g_rand),
			   mass1, niter, nrestart, rsqprop, PRECISION, 
			   ODD, &finalrsq, fn);
          /* Multiply by -Madjoint */
          dslash_site( F_OFFSET(g_rand), F_OFFSET(quark_prop), EVEN, fn);
          scalar_mult_latvec( F_OFFSET(g_rand), -2.0*mass1, F_OFFSET(quark_prop),ODD);
      }
      FORALLSITES(i,s){ lightprop[color][i] = lattice[i].quark_prop; }
scalar_mult_latvec( F_OFFSET(quark_prop), -1.0, F_OFFSET(g_rand), EVENANDODD );
 check_invert( F_OFFSET(g_rand), F_OFFSET(quark_source), mass1, tol, fn);

      /* repeat for heavy quark */
      if(t_source%2 == 0) {
         cgn += ks_congrad( F_OFFSET(quark_source), F_OFFSET(g_rand),
			    mass2, niter, nrestart, rsqprop, PRECISION, 
			    EVEN, &finalrsq, fn);
         /* Multiply by -Madjoint */
         dslash_site( F_OFFSET(g_rand), F_OFFSET(quark_prop), ODD, fn);
         scalar_mult_latvec( F_OFFSET(g_rand), -2.0*mass2, F_OFFSET(quark_prop),EVEN);
      }
      else {
        cgn += ks_congrad( F_OFFSET(quark_source), F_OFFSET(g_rand),
			   mass2, niter, nrestart, rsqprop, PRECISION, 
			   ODD, &finalrsq, fn);
          /* Multiply by -Madjoint */
          dslash_site( F_OFFSET(g_rand), F_OFFSET(quark_prop), EVEN, fn);
          scalar_mult_latvec( F_OFFSET(g_rand), -2.0*mass2, F_OFFSET(quark_prop),ODD);
      }
      FORALLSITES(i,s){ heavyprop[color][i] = lattice[i].quark_prop; }

      /* TEMP: test inversion, */
scalar_mult_latvec( F_OFFSET(quark_prop), -1.0, F_OFFSET(g_rand), EVENANDODD );
 check_invert( F_OFFSET(g_rand), F_OFFSET(quark_source), mass2, tol, fn);
    } /* end color loop*/

    /* add contributions into propagators */
    /* measure the meson propagator */
    for(t=0; t<nt; t++){
        /* define the time value offset t from t_source */
        t_off = (t+t_source)%nt;
        
        for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)
  	for(color=0;color<3;color++) {
  	    if( node_number(x,y,z,t_off) != mynode() )continue;
  	    i=node_index(x,y,z,t_off);

	    /* light-light mesons */
  	    cc = su3_dot( &lightprop[color][i],
  			 &lightprop[color][i] );
  	    
  	    CSUM( props[prop_pion5_ll][t], cc )
  	    
  	    if( (x+y)%2==0)CSUM( props[prop_rhoi0_ll][t], cc )
  	    else CSUB( props[prop_rhoi0_ll][t], cc, props[prop_rhoi0_ll][t] )
  	    if( (y+z)%2==0)CSUM( props[prop_rhoi0_ll][t], cc )
  	    else CSUB( props[prop_rhoi0_ll][t], cc,  props[prop_rhoi0_ll][t] )
  	    if( (z+x)%2==0)CSUM( props[prop_rhoi0_ll][t], cc )
  	    else CSUB( props[prop_rhoi0_ll][t], cc, props[prop_rhoi0_ll][t] )
  	    
  	    if( x%2==0)CSUM( props[prop_rhoi_ll][t], cc )
  	    else CSUB( props[prop_rhoi_ll][t], cc,  props[prop_rhoi_ll][t] )
  	    if( y%2==0)CSUM( props[prop_rhoi_ll][t], cc )
  	    else CSUB( props[prop_rhoi_ll][t], cc, props[prop_rhoi_ll][t] )
  	    if( z%2==0)CSUM( props[prop_rhoi_ll][t], cc )
  	    else CSUB( props[prop_rhoi_ll][t], cc, props[prop_rhoi_ll][t] )
  	    
  	    if( (x+y+z)%2==0)CSUM( props[prop_pion05_ll][t], cc )
  	    else CSUB( props[prop_pion05_ll][t], cc, props[prop_pion05_ll][t] )
  	    
	    /* light-heavy mesons */
  	    cc = su3_dot( &lightprop[color][i],
  			 &heavyprop[color][i] );
  	    
  	    CSUM( props[prop_pion5_lh][t], cc )
  	    
  	    if( (x+y)%2==0)CSUM( props[prop_rhoi0_lh][t], cc )
  	    else CSUB( props[prop_rhoi0_lh][t], cc, props[prop_rhoi0_lh][t] )
  	    if( (y+z)%2==0)CSUM( props[prop_rhoi0_lh][t], cc )
  	    else CSUB( props[prop_rhoi0_lh][t], cc,  props[prop_rhoi0_lh][t] )
  	    if( (z+x)%2==0)CSUM( props[prop_rhoi0_lh][t], cc )
  	    else CSUB( props[prop_rhoi0_lh][t], cc, props[prop_rhoi0_lh][t] )
  	    
  	    if( x%2==0)CSUM( props[prop_rhoi_lh][t], cc )
  	    else CSUB( props[prop_rhoi_lh][t], cc,  props[prop_rhoi_lh][t] )
  	    if( y%2==0)CSUM( props[prop_rhoi_lh][t], cc )
  	    else CSUB( props[prop_rhoi_lh][t], cc, props[prop_rhoi_lh][t] )
  	    if( z%2==0)CSUM( props[prop_rhoi_lh][t], cc )
  	    else CSUB( props[prop_rhoi_lh][t], cc, props[prop_rhoi_lh][t] )
  	    
  	    if( (x+y+z)%2==0)CSUM( props[prop_pion05_lh][t], cc )
  	    else CSUB( props[prop_pion05_lh][t], cc, props[prop_pion05_lh][t] )
  	    
	    /* heavy-heavy mesons */
  	    cc = su3_dot( &heavyprop[color][i],
  			 &heavyprop[color][i] );
  	    
  	    CSUM( props[prop_pion5_hh][t], cc )
  	    
  	    if( (x+y)%2==0)CSUM( props[prop_rhoi0_hh][t], cc )
  	    else CSUB( props[prop_rhoi0_hh][t], cc, props[prop_rhoi0_hh][t] )
  	    if( (y+z)%2==0)CSUM( props[prop_rhoi0_hh][t], cc )
  	    else CSUB( props[prop_rhoi0_hh][t], cc,  props[prop_rhoi0_hh][t] )
  	    if( (z+x)%2==0)CSUM( props[prop_rhoi0_hh][t], cc )
  	    else CSUB( props[prop_rhoi0_hh][t], cc, props[prop_rhoi0_hh][t] )
  	    
  	    if( x%2==0)CSUM( props[prop_rhoi_hh][t], cc )
  	    else CSUB( props[prop_rhoi_hh][t], cc,  props[prop_rhoi_hh][t] )
  	    if( y%2==0)CSUM( props[prop_rhoi_hh][t], cc )
  	    else CSUB( props[prop_rhoi_hh][t], cc, props[prop_rhoi_hh][t] )
  	    if( z%2==0)CSUM( props[prop_rhoi_hh][t], cc )
  	    else CSUB( props[prop_rhoi_hh][t], cc, props[prop_rhoi_hh][t] )
  	    
  	    if( (x+y+z)%2==0)CSUM( props[prop_pion05_hh][t], cc )
  	    else CSUB( props[prop_pion05_hh][t], cc, props[prop_pion05_hh][t] )
  	    
  	  } /* color */
        
      } /* nt-loop */
    
    /* measure the baryon propagator */
    for(t=0; t<nt; t++) {
        /* define the time value offset t from t_source */
        t_off = (t+t_source)%nt;
        
        for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
  	  if( node_number(x,y,z,t_off) != mynode() )continue;
  	  i=node_index(x,y,z,t_off);

	  /* three light quarks  */
	  for(color=0;color<3;color++){
	    (tmat.e[0][color]) = lightprop[0][i].c[color];
	    (tmat.e[1][color]) = lightprop[1][i].c[color];
	    (tmat.e[2][color]) = lightprop[2][i].c[color];
	  }
  	  cc = det_su3( &tmat );
  	      
  	  /* must get sign right.  This looks to see if we have
  		wrapped around the lattice.  "t" is the distance
  		from the source to the measurement, so we are
  		trying to find out if t_source+t is greater than
  		or equal to nt.  the "-tsource/nt" is in there
  		so that it will work correctly with tsource=nt.  */
  	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
  	     CSUM( props[prop_nuc_lll][t], cc )
  	  else  /* change sign because antiperiodic b.c.  sink point
  		should really be in a copy of the lattice */
  	     CSUB( props[prop_nuc_lll][t], cc, props[prop_nuc_lll][t] )

	  /* two lights and one heavy */
	  for(color=0;color<3;color++){
	    (tmat.e[0][color]) = lightprop[0][i].c[color];
	    (tmat.e[1][color]) = lightprop[1][i].c[color];
	    (tmat.e[2][color]) = heavyprop[2][i].c[color];
	  }
  	  cc = det_su3( &tmat );
  	      
  	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
  	     CSUM( props[prop_nuc_llh][t], cc )
  	  else CSUB( props[prop_nuc_llh][t], cc, props[prop_nuc_llh][t] )

	  /* now repeat for hhl baryon */
	  for(color=0;color<3;color++){
	    (tmat.e[0][color]) = lightprop[0][i].c[color];
	    (tmat.e[1][color]) = heavyprop[1][i].c[color];
	    (tmat.e[2][color]) = heavyprop[2][i].c[color];
	  }
  	  cc = det_su3( &tmat );
  	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
  	     CSUM( props[prop_nuc_lhh][t], cc )
  	  else  CSUB( props[prop_nuc_lhh][t], cc, props[prop_nuc_lhh][t] )

	  /* and hhh baryon (nonexistent particle!!) */
	  for(color=0;color<3;color++){
	    (tmat.e[0][color]) = heavyprop[0][i].c[color];
	    (tmat.e[1][color]) = heavyprop[1][i].c[color];
	    (tmat.e[2][color]) = heavyprop[2][i].c[color];
	  }
  	  cc = det_su3( &tmat );
  	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
  	     CSUM( props[prop_nuc_hhh][t], cc )
  	  else  CSUB( props[prop_nuc_hhh][t], cc, props[prop_nuc_hhh][t] )

        }
      } /* nt-loop */

  } /* end loop on t_source */

  /* Sum propagator arrays over nodes */
  /* print out propagators */
  g_veccomplexsum( props[0], nprops*nt );
  for(i=0;i<nprops;i++)for(j=0;j<nt;j++){
    CDIVREAL(props[i][j],n_sources,props[i][j]);
  }
  if(this_node==0){
  
    /*meson propagators*/
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass1,mass1 );
    printf("SOURCE: CORNER\n");
    printf("SINKS: PION_5 PION_05 RHO_i RHO_i0 \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e %e %e\n",j,
        props[prop_pion5_ll][j].real, props[prop_pion5_ll][j].imag,
        props[prop_pion05_ll][j].real, props[prop_pion05_ll][j].imag,
        props[prop_rhoi_ll][j].real, props[prop_rhoi_ll][j].imag,
        props[prop_rhoi0_ll][j].real, props[prop_rhoi0_ll][j].imag);
    }
    printf("ENDPROP\n");
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass1,mass2 );
    printf("SOURCE: CORNER\n");
    printf("SINKS: PION_5 PION_05 RHO_i RHO_i0 \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e %e %e\n",j,
        props[prop_pion5_lh][j].real, props[prop_pion5_lh][j].imag,
        props[prop_pion05_lh][j].real, props[prop_pion05_lh][j].imag,
        props[prop_rhoi_lh][j].real, props[prop_rhoi_lh][j].imag,
        props[prop_rhoi0_lh][j].real, props[prop_rhoi0_lh][j].imag);
    }
    printf("ENDPROP\n");
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n",mass2,mass2 );
    printf("SOURCE: CORNER\n");
    printf("SINKS: PION_5 PION_05 RHO_i RHO_i0 \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e %e %e\n",j,
        props[prop_pion5_hh][j].real, props[prop_pion5_hh][j].imag,
        props[prop_pion05_hh][j].real, props[prop_pion05_hh][j].imag,
        props[prop_rhoi_hh][j].real, props[prop_rhoi_hh][j].imag,
        props[prop_rhoi0_hh][j].real, props[prop_rhoi0_hh][j].imag);
    }
    printf("ENDPROP\n");
  
    /* Baryon propagators */
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e   %e\n",mass1,mass1,mass1);
    printf("SOURCE: CORNER\n"); printf("SINKS: NUCLEON \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e\n",j,
        props[prop_nuc_lll][j].real, props[prop_nuc_lll][j].imag);
    }
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e   %e\n",mass1,mass1,mass2);
    printf("SOURCE: CORNER\n"); printf("SINKS: NUCLEON \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e\n",j,
        props[prop_nuc_llh][j].real, props[prop_nuc_llh][j].imag);
    }
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e   %e\n",mass1,mass2,mass2);
    printf("SOURCE: CORNER\n"); printf("SINKS: NUCLEON \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e\n",j,
        props[prop_nuc_lhh][j].real, props[prop_nuc_lhh][j].imag);
    }
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e   %e\n",mass2,mass2,mass2);
    printf("SOURCE: CORNER\n"); printf("SINKS: NUCLEON \n");
    for(j=0;j<nt;j++){
      printf("%d %e %e\n",j,
        props[prop_nuc_hhh][j].real, props[prop_nuc_hhh][j].imag);
    }
    printf("ENDPROP\n");

    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  for(color=0;color<3;color++){
    free( lightprop[color] );
    free( heavyprop[color] );
  }
  
  return(cgn);
} /* spectrum_nd */

