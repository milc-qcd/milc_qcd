/******** mom_spec5.c *************/
/* MIMD version 6 */
/* Measure pi and rho propagators at momentum zero and one.
   We also need the flavor singlet (nonlocal) rho. */
/* The Kogut-Susskind phase factor -1^{sum_(nu<mu) eta_nu} is absorbed
   in the gauge links, so phase factors for nonlocal mesons implicitly
   contain this factor. */
/* average over origins of hypercubes */
/* The gauge should already be fixed when this routine is called. */

#include "ks_dyn_includes.h"
#define PI2 (2.0*PI)

int mom_spec() /* return the C.G. iteration number */
{
  Real mass_x2;
  Real finalrsq, scale;
  register site *st;
  register int i,x,y,z,t,dir,momdir;
  register complex cc,cc1,cc2,cc3;
  int icol, cgn, source_t, distance, n_source;
  msg_tag * tag[NDIRS];
  int coords[4],dims[4];

  /* We need propagators for local pi, local rho, nonlocal pi and
	nonlocal rho at momenta 0, 1 and -1. For all but the local
	pi we need a direction index too.
	Put these in an array of propagators indexed by species, direction
	and momentum; define a macro to remember which propagator is
	which. */
#define N_MOM 8		/* number of momentum values */
#define MOM_ZERO TUP	/* values for momentum - use previously defined 
			macros for other values: XUP = 1,0,0   
			XDOWN = -1,0,0 etc.  Steal TUP for 0,0,0  */

#define N_SPECIES 4	/* Number of hadron species */
#define PI_5_5 0	/* local pion */
#define PI_5_MU5 1	/* nonlocal pion */
#define RHO_MU_MU 2	/* local rho */
#define RHO_MU_1 3	/* nonlocal rho */

#define PROP_INDEX(species,dir,momentum,time) \
	(((species*4+dir)*N_MOM+momentum)*nt+time)
  complex *meson_prop;

  /* Fix TUP Coulomb gauge - gauge links only*/
  rephase( OFF );
  gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  rephase( ON );

  /* Allocate array for hadron propagators. */
  meson_prop = (complex *)malloc(
	PROP_INDEX(N_SPECIES,0,0,0)*sizeof(complex) );
	/* This macro gives one more than the highest possible index */

  /* initialize the hadron propagators */
  for(i=0; i< PROP_INDEX(N_SPECIES,0,0,0); i++)meson_prop[i] = cmplx(0.0,0.0);

  /* need dimensions in an array */
  dims[XUP] = nx; dims[YUP] = ny; dims[ZUP] = nz; dims[TUP] = nt;

  mass_x2 = 2.*mass;
  cgn=0;

  /* loop over source time slices */
  for ( n_source=0, source_t = 0; source_t < nt;
     n_source++, source_t += nt/4 ){
     /* For now, do all source colors */
     for(icol=0; icol<3; icol++){
         
         /* Put the source in vector "phi" */
	 /* 
		Take (M*M_adjoint)^(-1) times source,
		then multiply by -M_adjoint
		to get 1/M times source, separately for the two sources
		(For KS matrix, M_adjoint*M = M*M_adjoint).
	 */
         /* initialize phi and xxx */
         clear_latvec( F_OFFSET(phi), EVENANDODD);
         clear_latvec( F_OFFSET(xxx), EVENANDODD);

	 /* even site source is on 0,0,0,0 site of hypercube, we use a
		zero momentum source for the zero momentum mesons and a
		source with components with momentum +-1 in each direction
		for the nonzero momentum mesons. This will be used as the
	        antiquark source. */
	 /* zero momentum source */
         for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	   if( node_number(x,y,z,source_t) == this_node ){
	      i=node_index(x,y,z,source_t);
	      lattice[i].phi.c[icol] = cmplx(-1.0,0.0);
	   }
         }
         /* do a C.G. */
	 load_ferm_links(&fn_links);
         cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			   niter, rsqprop, PRECISION, EVEN, &finalrsq,
			   &fn_links);
         /* Now solution vector is in xxx */
         /* Multiply by -Madjoint, even site source -> propmat[0] */
         dslash_site( F_OFFSET(xxx), F_OFFSET(propmat[0]), ODD, &fn_links);
         FOREVENSITES(i,st){
            scalar_mult_su3_vector(&(st->xxx),-mass_x2, &(st->propmat[0]));
         }

         /* reinitialize phi and xxx even sites */
         clear_latvec( F_OFFSET(phi), EVEN);
         clear_latvec( F_OFFSET(xxx), EVEN);
	 /* nonzero momentum source */
         for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	   /* even site source is on 0,0,0,0 site of hypercube, we use a
		zero momentum source for the zero momentum mesons and a
		source with components with momentum +-1 in each direction
		for the nonzero momentum mesons. This will be used as the
	        antiquark source. */
	   if( node_number(x,y,z,source_t) == this_node ){
	      i=node_index(x,y,z,source_t);
	      lattice[i].phi.c[icol].real  = 
		- (Real)cos( (double)(PI2*(Real)(x & 0xffff)/(Real)nx) )
	        - (Real)cos( (double)(PI2*(Real)(y & 0xffff)/(Real)ny) )
	        - (Real)cos( (double)(PI2*(Real)(z & 0xffff)/(Real)nz) );
	   }
         }
         /* do a C.G. */
	 load_ferm_links(&fn_links);
         cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			   niter, rsqprop, PRECISION, EVEN, &finalrsq,
			   &fn_links);
         /* Now solution vector is in xxx */
         /* Multiply by -Madjoint, even site nonzero mom. source -> propmat[1]*/
         dslash_site( F_OFFSET(xxx), F_OFFSET(propmat[1]), ODD, &fn_links);
         FOREVENSITES(i,st){
            scalar_mult_su3_vector(&(st->xxx),-mass_x2, &(st->propmat[1]));
         }

         for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	   /* Odd site source is even site source displaced in three
		directions, to couple to nonlocal rho.  This will be
	        used as the quark source. */
	   if( node_number(x+1,y,z,source_t) == this_node ){
	      i=node_index(x+1,y,z,source_t);
	      lattice[i].phi.c[icol].real = -1.0;
	   }
	   if( node_number(x,y+1,z,source_t) == this_node ){
	      i=node_index(x,y+1,z,source_t);
	      lattice[i].phi.c[icol].real = -1.0;
	   }
	   if( node_number(x,y,z+1,source_t) == this_node ){
	      i=node_index(x,y,z+1,source_t);
	      lattice[i].phi.c[icol].real = -1.0;
	   }
         }
         /* do a C.G. */
	 load_ferm_links(&fn_links);
         cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			   niter, rsqprop, PRECISION, ODD, &finalrsq,
			   &fn_links);
         /* Now solution vector is in xxx */
         /* Multiply by -Madjoint, odd site source -> ttt */
         dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), EVEN, &fn_links);
         FORODDSITES(i,st){
            scalar_mult_su3_vector(&(st->xxx),-mass_x2, &(st->ttt));
         }
   
        /* measure the meson propagators */
        /* parallel transport the propagator forward and backward in all
           directions. */
        for(dir=XUP;dir<=TUP;dir++){
             tag[dir] = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector),
	       dir, EVENANDODD, gen_pt[dir] );
           FORALLSITES(i,st){
	     mult_adj_su3_mat_vec( &(st->link[dir]),
                   &(st->ttt),  &(st->tempvec[dir]) );
           }
           tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	      sizeof(su3_vector), OPP_DIR(dir), EVENANDODD,
	      gen_pt[OPP_DIR(dir)] );
        }
   
        for(dir=XUP;dir<=TUP;dir++){
            wait_gather(tag[dir]);
            wait_gather(tag[OPP_DIR(dir)]);
	    FORALLSITES(i,st){
	       coords[XUP] = st->x; coords[YUP] = st->y;
	       coords[ZUP] = st->z; coords[TUP] = st->t;
	       distance = (st->t - source_t + nt)%nt;
	       mult_su3_mat_vec( &(st->link[dir]),
                   (su3_vector *)gen_pt[dir][i],  &(st->cg_p) );
               /* now cg_p contains the propagator parallel transported from
	       displacement +dir, and *gen_pt[OPP_DIR(dir)] the propagator
	       parallel transported from -dir. */
   

	        /* local pion and rho */
	        /* local rho:  sign = (-1)^coords[dir] */
		/* "cc" is for pion, "cc1" for rho */

		/* for zero momentum, use both antiquark and quark source
		   at zero momentum */
	        cc  = su3_dot( &(st->propmat[0]), &(st->propmat[0]) );
	        if( coords[dir]%2 != 0){ CMULREAL(cc,-1.0,cc1); }
		else cc1 = cc;
	        CSUM( meson_prop[ PROP_INDEX(PI_5_5,dir,MOM_ZERO,distance) ],
		    cc);
	        CSUM(meson_prop[ PROP_INDEX(RHO_MU_MU,dir,MOM_ZERO,distance) ],
		    cc1);

		/* for nonzero momentum, use one source at nonzero momentum */
	        cc  = su3_dot( &(st->propmat[1]), &(st->propmat[0]) );
	        if( coords[dir]%2 != 0){ CMULREAL(cc,-1.0,cc1); }
		else cc1 = cc;
	        for( momdir=XUP; momdir <= ZUP; momdir++){
	            cc2 = ce_itheta(
			(Real)(PI2*(Real)(coords[momdir] & 0xffff)
			/(Real)dims[momdir]) );
	            CMUL(cc,cc2,cc3);
	            CSUM( meson_prop[ PROP_INDEX(PI_5_5,dir,momdir,distance)],
		       cc3);
	            CMUL(cc1,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(RHO_MU_MU,dir,momdir,distance)],
			cc3);
	            cc2 = ce_itheta(
			(Real)(-PI2*(Real)(coords[momdir] & 0xffff)
		       /(Real)dims[momdir]) );
	            CMUL(cc,cc2,cc3);
	            CSUM(meson_prop[
		       PROP_INDEX(PI_5_5,dir,OPP_DIR(momdir),distance)], cc3);
	            CMUL(cc1,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(RHO_MU_MU,dir,OPP_DIR(momdir),
		       distance) ], cc3);
	        }
   

	        /* nonlocal pi:  sign = (-1)^sum(coords[ !=dir ]) */
	        /* nonlocal rho:  sign = (-1)^sum(coords[ dir ]) */
		/* "cc" is for pion, "cc1" for rho */

		/* for zero momentum, use both antiquark and quark source
		   at zero momentum */
	        cc1  = su3_dot( &(st->propmat[0]), &(st->cg_p) );
	        cc  = su3_dot( &(st->propmat[0]), 
			       (su3_vector *)gen_pt[OPP_DIR(dir)][i]);
		CSUM(cc,cc1);
	        if( (coords[XUP]+coords[YUP]+coords[ZUP]+coords[TUP])%2 != 0)
		    { CMULREAL(cc,-1.0,cc1); }
		else {cc1 = cc;}		/* rho */
	        if( (coords[XUP]+coords[YUP]+coords[ZUP]+coords[TUP]
		    -coords[dir]) %2 != 0) { CMULREAL(cc,-1.0,cc); }
		else {}  /* pion */
	        CSUM(meson_prop[ PROP_INDEX(PI_5_MU5,dir,MOM_ZERO,distance) ],
		    cc);
	        CSUM(meson_prop[ PROP_INDEX(RHO_MU_1,dir,MOM_ZERO,distance) ],
		    cc1);

		/* for nonzero momentum, use one source at nonzero momentum */
	        cc1  = su3_dot( &(st->propmat[1]), &(st->cg_p) );
	        cc  = su3_dot( &(st->propmat[1]), 
			       (su3_vector *)gen_pt[OPP_DIR(dir)][i]);
		CSUM(cc,cc1);
	        if( (coords[XUP]+coords[YUP]+coords[ZUP]+coords[TUP])%2 != 0)
		    { CMULREAL(cc,-1.0,cc1); }
		else {cc1 = cc;}		/* rho */
	        if( (coords[XUP]+coords[YUP]+coords[ZUP]+coords[TUP]
		    -coords[dir]) %2 != 0) { CMULREAL(cc,-1.0,cc); }
		else {}  /* pion */
	        for( momdir=XUP; momdir <= ZUP; momdir++){
		    cc2 = ce_itheta(
			(Real)(PI2*(Real)(coords[momdir] & 0xffff)
			/(Real)dims[momdir]) );
	            CMUL(cc,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(PI_5_MU5,dir,momdir,distance) ],
			cc3);
	            CMUL(cc1,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(RHO_MU_1,dir,momdir,distance) ],
			cc3);
	            cc2 = ce_itheta(
		       (Real)(-PI2*(Real)(coords[momdir] & 0xffff)
			/(Real)dims[momdir]) );
	            CMUL(cc,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(PI_5_MU5,dir,OPP_DIR(momdir),
			distance) ], cc3);
	            CMUL(cc1,cc2,cc3);
	            CSUM(meson_prop[ PROP_INDEX(RHO_MU_1,dir,OPP_DIR(momdir),
			distance) ], cc3);
	        }
   
	     }/* loop over sites */
             cleanup_gather(tag[dir]);
             cleanup_gather(tag[OPP_DIR(dir)]);
         }/* direction loop */
     } /* loop on colors */
  } /* loop on source time slices */

  /* sum over nodes and normalize by number of sources */
  for(i=0; i< PROP_INDEX(N_SPECIES,0,0,0); i++){
      g_complexsum( &meson_prop[i] );
      scale = (Real)(1.0/(Real)n_source);
      CMULREAL( meson_prop[i], scale, meson_prop[i]);
  }

  /* Print results */
  if(this_node==0){
    /* Header line so I can remember format */
    printf("MESON_PROP:\tSPECIES\tPOLARIZATION\tMOMENTUM\n");
    printf("distance\treal_part\timag_part\n");
    printf("END_PROP\n");
    printf("END_SET at end of all for this lattice\n");
    /* Header line so I can remember format */
    for( i=0; i<N_SPECIES; i++){	/* loop over species */

      for(dir=XUP;dir<=ZUP;dir++){
	if( i == PI_5_5 && dir != XUP )continue;

        for(momdir=0; momdir<N_MOM; momdir++){
	  if(momdir == TDOWN)continue;
          switch(i){
	    case PI_5_5: printf("PI_5_5\t"); break;
	    case PI_5_MU5: printf("PI_5_MU5\t"); break;
	    case RHO_MU_MU: printf("RHO_MU_MU\t"); break;
	    case RHO_MU_1: printf("RHO_MU_1\t"); break;
          }
	  printf("%d\t",dir);
	  switch(momdir){
	    case(MOM_ZERO): printf("0,0,0\n"); break;
	    case(XUP): printf("1,0,0\n"); break;
	    case(YUP): printf("0,1,0\n"); break;
	    case(ZUP): printf("0,0,1\n"); break;
	    case(XDOWN): printf("-1,0,0\n"); break;
	    case(YDOWN): printf("0,-1,0\n"); break;
	    case(ZDOWN): printf("0,0,-1\n"); break;
	  }

	  for( t=0; t<nt; t++){
	    printf("%d\t%e\t%e\n",t,
		(double)meson_prop[ PROP_INDEX(i,dir,momdir,t) ].real,
		(double)meson_prop[ PROP_INDEX(i,dir,momdir,t) ].imag);
	  }
          printf("END_PROP\n");
	}
      }
    }
    printf("END_SET\n");
  }/* if this_node==0 */
  
  free( meson_prop );
  return(cgn);
} /* spectrum */
