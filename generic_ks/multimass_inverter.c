/******** multimass_inverter.c *************/
/* MIMD version 6*/
/* SG 9/03   based on fpi_2.c
   point sources only for now

   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie rephase(ON) )

    Coulomb gauge should be fixed before this routine for wall source

   Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.  Assumes "dslash" has
   been defined to be the appropriate "dslash_fn" or "dslash_eo"

    Normalize sources to 2, since ks_congrad inverts matrix with 2m on diagonal,
    which is twice the usual convention.
*/
enum prop_name {
    POINTPOINT,	/* point source point sink */
    POINTWALL, 	/* point source wall sink */
    WALLPOINT, 	/* wall source point sink */
    WALLWALL, 	/* wall source wall sink */
    NPROPS	/* number of propagators */
};

#include "generic_ks_includes.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

#ifdef FN
#define dslash dslash_fn
#define dslash_on_temp dslash_fn_on_temp
#endif
#ifdef EO
#define dslash dslash_eo
#define dslash_on_temp dslash_eo_on_temp /* actually isn't written yet (4/01) */
#endif

void f2d_vector(fsu3_vector *, su3_vector *);
void d2f_vector(su3_vector *, fsu3_vector *);
int test_converge(int t_source);
#define MAX_RECXML 65

int multimass_inverter( Real *masses, int nmasses, Real tol){
  /* arguments are array of masses, number of masses,
     tolerance for inverter check.
     return C.G. iteration number */

  int cgn;
  register int i,j,m1,m2;
  register site* s;
  register complex cc,czero;
  register int t_source;
  int sourcevec,color;	/* color for source */
  int src_count; /* number of source time slices used */
  Real finalrsq;
  su3_vector **quark_props, **quark_props_color;
  su3_vector *temp_prop;
  complex **props;	/* arrays of propagators */
  su3_vector *wall_sink_m1,*wall_sink_m2; 
  char kssavefile_tmp[MAXFILENAME],extension[24];
  int len1,len2;
  ks_prop_file *kspf; /* structure for opening ks prop file */
  char recxml[MAX_RECXML];

  cgn=0; /* number of CG iterations */
  czero.real = czero.imag = 0.0;

  /* allocate space for quark propagators and zero them */
  /* include a factor of 3 for source colors */
  quark_props = (su3_vector **)malloc(3*nmasses*sizeof(su3_vector *));
  for(j=0; j<3*nmasses; j++){
      quark_props[j] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      FORALLSITES(i,s)  clearvec( &(quark_props[j][i])) ; 
  }
   
  temp_prop = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  wall_sink_m1 = (su3_vector *)malloc(nt*sizeof(su3_vector));
  wall_sink_m2 = (su3_vector *)malloc(nt*sizeof(su3_vector));

  /* allocate space for meson propagators (NPROPS numbers per timeslice) */
  /* for meson with masses number m1 and m2, use props[NPROPS*(nmasses*m1+m2)+TYPE] */
  props = (complex **)malloc( NPROPS*nmasses*nmasses*sizeof(complex *) );
  props[0] = (complex *)malloc( NPROPS*nmasses*nmasses*nt*sizeof(complex) );
  for(i=1; i<NPROPS*nmasses*nmasses; i++) props[i] = props[i-1] + nt;

  /* set propagators to zero */
  for(i=0;i<NPROPS*nmasses*nmasses;i++)for(j=0;j<nt;j++){
    props[i][j]=czero;
  }

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    if(this_node==0)printf("multimass_inverter(): source time = %d\n",
      t_source);

    /* point source at spatial origin */
    for(color=0;color<3;color++){
	    clear_latvec( F_OFFSET(quark_source), EVENANDODD );

		if( node_number(0,0,0,t_source) == mynode() )
		   {
		     i=node_index(0,0,0,t_source);
		     lattice[i].quark_source.c[color].real
		                         = 1.0; 
		     /* this should have the same normalization as the
		      * the propagators generated earlier in the decay constant
		      * program with Fermilab */
		   }

	/* compute M^-1 * quark_source */

		/* must use the correct block of quark_props */
		quark_props_color = &(quark_props[nmasses*color]);
	cgn += ks_multicg( F_OFFSET(quark_source), quark_props_color, masses, nmasses,
	   niter, rsqprop, EVEN, &finalrsq);
	/* Multiply by Madjoint */
	for(j=0;j<nmasses;j++){
	    dslash_on_temp( quark_props_color[j], temp_prop, ODD);
	    FOREVENSITES(i,s)
		scalar_mult_su3_vector( &(quark_props_color[j][i]), 
				2.0*masses[j], &(quark_props_color[j][i]) );
	    FORODDSITES(i,s)
		scalar_mult_su3_vector( &(temp_prop[i]), -1.0,
		    &(quark_props_color[j][i]) );

	    FORALLSITES(i,s) s->ttt = quark_props_color[j][i];
	    check_invert( F_OFFSET(ttt), F_OFFSET(quark_source), masses[j], tol );
	} /* j=masses */

	/* 0-+ (kaon) propagators */
	for(m1=0;m1<nmasses;m1++)for(m2=m1;m2<nmasses;m2++){
	    for(i=0;i<nt;i++){ 
		clearvec( &(wall_sink_m1[i]) );
		clearvec( &(wall_sink_m2[i]) );
	    }
	    FORALLSITES(i,s){
		add_su3_vector( &(wall_sink_m1[s->t]),&(quark_props_color[m1][i]), &(wall_sink_m1[s->t]) );
		add_su3_vector( &(wall_sink_m2[s->t]),&(quark_props_color[m2][i]), &(wall_sink_m2[s->t]) );
	        cc = su3_dot( &(quark_props_color[m1][i]), &(quark_props_color[m2][i]) );
	        CSUM( props[NPROPS*(nmasses*m1+m2)+POINTPOINT][(s->t+nt-t_source)%nt], cc );
	    }
	    g_veccomplexsum( (complex *)wall_sink_m1, 3*nt);
	    g_veccomplexsum( (complex *)wall_sink_m2, 3*nt);
	    for(i=0;i<nt;i++){ 
	        cc = su3_dot( &(wall_sink_m1[i]), &(wall_sink_m2[i]) );
	        CSUM( props[NPROPS*(nmasses*m1+m2)+POINTWALL][(i+nt-t_source)%nt], cc );
	    }
	} /* m1,m2=masses */

    } /* source color */
    /* Save quark propagators here */
    /* create name for each mass and copy to propmat */
    len1=strlen(kssavefile);
    
     for(j=0;j<nmasses;j++){
	     /* create a file name for this propagator */
    	     strcpy(kssavefile_tmp, kssavefile); 
	     sprintf(extension,"_m%.4f_nt%d",masses[j],t_source);
	     len2=strlen(extension);
	     len2 += len1;   /* total string length */
		     if (len2 >= MAXFILENAME)  {
		     	node0_printf(
			"multimass_inverter: WARNING MAXFILENAME exceeded\n");
			node0_printf("%s and %s\n",kssavefile_tmp,extension);
		     }
	     strncat(kssavefile_tmp,extension,MAXFILENAME-1);
	     /* define propmass for this mass value for ks_info file */
    	     propmass = masses[j];
   		
        /* copy each color to propmat[color] */
	FORALLSITES(i,s){
	  s->propmat[0] = quark_props[j][i];
	  s->propmat[1] = quark_props[j+nmasses][i];
	  s->propmat[2] = quark_props[j+2*nmasses][i];
        }
        /* save KS propagator if requested */
	/* Create record XML.  Substitute string until we have XML support */
	snprintf(recxml,MAX_RECXML,"\nmass %g\nt_source %d\n",masses[j],
		 t_source);
	save_ksprop_from_site( kssaveflag, kssavefile_tmp, recxml, 
			       F_OFFSET(propmat), 0);
   } /* end loop on j (masses) */
  } /* end loop on t_source */


  /* Sum propagator arrays over nodes */
  /* print out propagators */
  for(i=0;i<NPROPS*nmasses*nmasses;i+=NPROPS){
      g_veccomplexsum( props[i+POINTPOINT] , nt );
      g_veccomplexsum( props[i+WALLPOINT] , nt );
  }
  for(i=0;i<nmasses*nmasses;i++){
      for(j=0;j<nt;j++){
          CDIVREAL(props[NPROPS*i+POINTPOINT][j],(double)n_sources*nx*ny*nz*nx*ny*nz,
	    props[NPROPS*i+POINTPOINT][j]);
          CDIVREAL(props[NPROPS*i+POINTWALL ][j],(double)n_sources*nx*ny*nz*nx*ny*nz,
	    props[NPROPS*i+POINTWALL ][j]);
          CDIVREAL(props[NPROPS*i+WALLPOINT ][j],3.0*n_sources*nx*ny*nz*nx*ny*nz,
	    props[NPROPS*i+WALLPOINT ][j]);
          CDIVREAL(props[NPROPS*i+WALLWALL  ][j],3.0*n_sources*nx*ny*nz*nx*ny*nz,
	    props[NPROPS*i+WALLWALL  ][j]);
      }
  }
  if(this_node==0){

    for(m1=0;m1<nmasses;m1++)for(m2=m1;m2<nmasses;m2++){
        printf("STARTPROP\n");
        printf("MASSES:  %.5e   %.5e\n",masses[m1],masses[m2]);
        printf("SOURCE: POINT\n");
        printf("SINKS: POINT_KAON_5 WALL_KAON_5\n");
        for(j=0;j<nt;j++){
	    printf("%d %e %e %e %e\n",j,
	      props[NPROPS*(nmasses*m1+m2)+POINTPOINT][j].real,
	      props[NPROPS*(nmasses*m1+m2)+POINTPOINT][j].imag,
	      props[NPROPS*(nmasses*m1+m2)+POINTWALL][j].real,
	      props[NPROPS*(nmasses*m1+m2)+POINTWALL][j].imag );
        }
        printf("ENDPROP\n");

    } /* masses */


    fflush(stdout);
  } /* end if(this_node==0) */

  /* free arrays */
  free(props[0]); free(props);
  for(i=0; i<nmasses; i++){
      free(quark_props[i]);
  }
  free(quark_props); free(temp_prop);
  free(wall_sink_m1); free(wall_sink_m2);
  
  return(cgn);
} /* multimass_inverter */

void f2d_vector(fsu3_vector *a, su3_vector *b){
  int i;

    for(i = 0; i < 3; i++){
      b->c[i].real = a->c[i].real;
      b->c[i].imag = a->c[i].imag;
    }
}
void d2f_vector(su3_vector *a, fsu3_vector *b){
  int i;

    for(i = 0; i < 3; i++){
      b->c[i].real = a->c[i].real;
      b->c[i].imag = a->c[i].imag;
    }
}
