/************************* fmunu_region.c *******************************/
/* Calculates the temporal and spatial field strength components */
/* and the topological charge for the sphaleron rate.            */

#include "wilson_flow_includes.h"
#include "../include/field_strength.h"
#include <string.h>


/* Computes the real trace of the su3 matrix product: ReTr(A*B) */
static Real
real_trace_nn( su3_matrix *a, su3_matrix *b ) {
  register int i,j;
  register complex x;
  register Real sum;

  sum = 0.0;
  for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) {
    CMUL( a->e[i][j], b->e[j][i], x );
    sum += x.real;
  }

  return sum;
}

/* Computes the field strength components and topological charge,
 * in the active bulk including the two boundaries */
void
fmunu_fmunu_bulk( double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;
	msg_tag *tag[6];

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fsi = NULL;
	su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 4 * sizeof(double) );

  // make_field_strength_bulk( F_OFFSET(link), F_OFFSET(fieldstrength) );
  /* Compute 8*F_mu,nu at each site for half-integer time */

  /* Compute 8*F_mu,nu at each site for half-integer time */
  #define NFS 9
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_bulk( fstrength );

	tag[2] = start_gather_field( fstrength[FS_XY] , sizeof(su3_matrix),
			                         TUP, EVENANDODD, gen_pt[2]);
	tag[1] = start_gather_field( fstrength[FS_XZ] , sizeof(su3_matrix),
                               TUP, EVENANDODD, gen_pt[1]);
	tag[0] = start_gather_field( fstrength[FS_YZ] , sizeof(su3_matrix),
                               TUP, EVENANDODD, gen_pt[0]);

  tag[5] = start_gather_field( fstrength[FS_XY+6] , sizeof(su3_matrix),
                               TUP, EVENANDODD, gen_pt[5]);
  tag[4] = start_gather_field( fstrength[FS_XZ+6] , sizeof(su3_matrix),
                               TUP, EVENANDODD, gen_pt[4]);
  tag[3] = start_gather_field( fstrength[FS_YZ+6] , sizeof(su3_matrix),
                               TUP, EVENANDODD, gen_pt[3]);

  for ( i=0; i<6; i++ ) wait_gather(tag[i]);

  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
 	IF_ACTIVE(s) 
  {
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 6
      fsi = &(fstrength[FS_XY+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    IF_LOWER_BULK(s) {
	    ft = &(fstrength[FS_ZT][i]);
	    time[0] -= real_trace_nn(ft, ft);
    	add_su3_matrix( fs, (su3_matrix*)(gen_pt[2][i]), &tmp );
    	fs = &tmp;
 	    charge[0] -= real_trace_nn(fs, ft);
      #if NFS > 6
        add_su3_matrix( fsi, (su3_matrix*)(gen_pt[5][i]), &tmp );
        fsi = &tmp;
        charge[1] -= real_trace_nn(fsi, ft);
      #endif
    }

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 6
      fsi = &(fstrength[FS_XZ+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    IF_LOWER_BULK(s) {
	    ft = &(fstrength[FS_YT][i]);
  	  time[0] -= real_trace_nn(ft, ft);
    	add_su3_matrix( fs, (su3_matrix*)(gen_pt[1][i]), &tmp );
    	fs = &tmp;
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
 	    charge[0] -= realtrace_su3(fs, ft);
      #if NFS > 6
        add_su3_matrix( fsi, (su3_matrix*)(gen_pt[4][i]), &tmp );
        fsi = &tmp;
        charge[1] -= realtrace_su3(fsi, ft);
      #endif
    }

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 6
      fsi = &(fstrength[FS_YZ+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    IF_LOWER_BULK(s) {
  	  ft = &(fstrength[FS_XT][i]);
	    time[0] -= real_trace_nn(ft, ft);
    	add_su3_matrix( fs, (su3_matrix*)(gen_pt[0][i]), &tmp );
    	fs = &tmp;
 	    charge[0] -= real_trace_nn(fs, ft);
      #if NFS > 6
        add_su3_matrix( fsi, (su3_matrix*)(gen_pt[3][i]), &tmp );
        fsi = &tmp;
        charge[1] -= real_trace_nn(fsi, ft);
      #endif
    }
  }

  for ( i=0; i<6; i++ ) cleanup_gather(tag[i]);

  /* Sum over all nodes */
  g_vecdoublesum( time,   1 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 2 );

  /* Norma(lizations */
  for ( i = 0; i < 1 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    time[i] /= 0.25; /* extended only forward in time, factor 1/2, then squared */
    time[i] /= LWR_BULK_LEN * 1. / nt; /* computed only on nt/2 integer time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    space[i] /= ACTV_VOL_LEN * 1. / nt; /* computed only on nt/2+1 integer time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    // charge[i] /= 0.5; /* extended only forward in time */
  }
  time[1] = charge[2] = charge[3] = 0.;
  destroy_field( &fstrength );
  #undef NFS
}

/* Computes the field strength components and topological charge,
 * in the active bulk including the two boundaries */
void
fmunu_fmunu_full( double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fti = NULL, *fsi = NULL;
  su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 4 * sizeof(double) );

  /* Compute 8*F_mu,nu at each site for half-integer time */
  #define NFS 12
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_full( fstrength );
    
  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  {
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    ft = &(fstrength[FS_ZT][i]);
    time[0] -= real_trace_nn(ft, ft);
    charge[0] -= real_trace_nn(fs, ft);
    #if NFS > 6
      fsi = &(fstrength[FS_XY+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
      charge[1] -= real_trace_nn(fsi, ft);
    #endif
    #if NFS > 9
      fti = &(fstrength[FS_ZT+6][i]);
      time[1] -= real_trace_nn(fti, fti);
      charge[2] -= real_trace_nn(fs, fti);
      charge[3] -= real_trace_nn(fsi, fti);
    #endif

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    ft = &(fstrength[FS_YT][i]);
    time[0] -= real_trace_nn(ft, ft);
    /* ANTICYCLIC! ReTr{ fs.dag * ft } */
    charge[0] -= realtrace_su3(fs, ft);
    #if NFS > 6
      fsi = &(fstrength[FS_XZ+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
      charge[1] -= realtrace_su3(fsi, ft);
    #endif
    #if NFS > 9
      fti = &(fstrength[FS_YT+6][i]);
      time[1] -= real_trace_nn(fti, fti);
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
      charge[2] -= realtrace_su3(fs, fti);
      charge[3] -= realtrace_su3(fsi, fti);
    #endif

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    ft = &(fstrength[FS_XT][i]);
    time[0] -= real_trace_nn(ft, ft);
    charge[0] -= real_trace_nn(fs, ft);
    #if NFS > 6
      fsi = &(fstrength[FS_YZ+6][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
      charge[1] -= real_trace_nn(fsi, ft);
    #endif
    #if NFS > 9
      fti = &(fstrength[FS_XT+6][i]);
      time[1] -= real_trace_nn(fti, fti);
      charge[2] -= real_trace_nn(fs, fti);
      charge[3] -= real_trace_nn(fsi, fti);
    #endif

  }

  /* Sum over all nodes */
  g_vecdoublesum( time,   2 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 4 );

  /* Norma(lizations */
  for ( i = 0; i < 2 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    // time[i] /= 0.25; /* extended only forward in time, factor 1/2, then squared */
    // time[i] /= LWR_BULK_LEN * 1. / nt; /* computed only on nt/2 integer time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    // space /= ACTV_VOL_LEN * 1. / nt; /* computed only on nt/2+1 integer time slices */
  }
  for ( i = 0; i < 4 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    // charge[i] /= 0.5; /* extended only forward in time */
  }
  destroy_field( &fstrength );
  #undef NFS
}


/* Compute loops, here only 3D spatial part,
	 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   in various regions */
static void
spatial_gauge_action_w_s_region( int region_flag, 
	double *wl1x1s, double *wl1x2s, su3_matrix **link ) {

  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
			dirb[2] = {NODIR,NODIR}, 
			diro[2] = {NODIR,NODIR};
  register int ig;
	#define NGATHER 4
  register site *s;
  msg_tag *tag[NGATHER];
  for( ig=0; ig<NGATHER; ig++ ) { tag[ig] = NULL; }
  // #define NTEMP_STORAGE 6
  // su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;
  double time_frac = 1.;
  switch ( region_flag ) {
  	case (FULLVOL):
  		time_frac = 1.;
  		break;
  	case (BOUNDARY):
  		time_frac = nt / 2.;
  		break;
  	case (BULK):
  		time_frac = nt * 1. / INR_BULK_LEN;
  		break;
  	case (LOWER_BULK):
  		time_frac = nt * 1. / LWR_BULK_LEN;
  		break;
  	case (ACTIVE):
  		time_frac = nt * 1. / ACTV_VOL_LEN;
  		break;
    case (LOWER_BOUNDARY):
      time_frac = nt / 1.;
      break;
    case (UPPER_BOUNDARY):
      time_frac = nt / 1.;
      break;
  	default:
  		node0_printf("Undefined REGION: %d\n",region_flag);
  		terminate(1);
  		break;
  }

  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
  //   su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  //   if(su3mat[ig] == NULL) {
  //       printf( "spatial_gauge_action_w_s_region: can't malloc su3mat[%d]\n", ig );
  //       fflush(stdout); terminate(1);
  //   }
  // }
  #define NTEMP 6
  su3_matrix **su3mat = new_field( NTEMP );

  // prepare accumulators
  *wl1x1s = 0;
  *wl1x2s = 0;

  for( dir[0]=YUP; dir[0]<TUP; dir[0]++ ) {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      setup_blocked_dirs( dir, dirb ,diro );
			#define LINK0 link[dir[0]]
			#define LINK1 link[dir[1]]

			ig = 0;
			// request LINK1 from direction dirb[0]
			tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
																		dirb[ig], EVENANDODD, gen_pt[ig]);

			ig = 1;
			// request LINK0 from direction dirb[1]
			tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
																		dirb[ig], EVENANDODD, gen_pt[ig]);

      // while waiting for gathers multiply two two links on the site
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
  			mult_su3_an( &(LINK1[i]), &(LINK0[i]), &(su3mat[0][i]) );

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );

      wait_gather(tag[1]);

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
 			{
        // mult_su3_nn( &(s->link[dir[1]]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        mult_su3_na( &tempmat, (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

			ig = 2;
			// request staple for "left" = su3mat[5] from dirb[1]
			tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
																		dirb[1], EVENANDODD, gen_pt[ig]);

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "bottom" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
 			{
        mult_su3_na( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        // mult_su3_na( &(su3mat[1][i]), &(s->link[dir1]), &(su3mat[4][i]) );
        mult_su3_na( &(su3mat[1][i]), &(LINK0[i]), &(su3mat[4][i]) );
      }

			ig = 3;
			// request staple for "bottom" = su3mat[4] from dirb[0]
			tag[ig] = start_gather_field( su3mat[4] , sizeof(su3_matrix),
																		dirb[0], EVENANDODD, gen_pt[ig]);
      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
 			{
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        *wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, region_flag) 
 			{
        // get the contribution of 1x2 rectangle extended in dir1
        // and of 1x1 plaquette to the accumulators
        *wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
        *wl1x1s += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      if ( tag[ig] != NULL ) 
          cleanup_gather(tag[ig]);

	    #undef LINK0
	    #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1s );
  g_doublesum( wl1x2s );

  // get densities
  *wl1x1s /= volume;
  *wl1x1s *= block_stride * block_stride * block_stride;
  *wl1x1s *= time_frac;
  *wl1x1s /= 3;
  *wl1x2s /= volume;
  *wl1x2s *= block_stride * block_stride * block_stride;
  *wl1x2s *= time_frac;
  *wl1x2s /= 6;
  
  // deallocate temporary storage
  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
  //   free( su3mat[ig] );
	// #undef NTEMP_STORAGE
  destroy_field( &su3mat );
  #undef NTEMP
	#undef NGATHER
}

/* Compute loops, here only temporal part,
	 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   restricted to the lower bulk */
static void
temporal_gauge_action_w_s_lwr_bulk( double *wl1x1t, 
  double *wl1x2t, double *wl2x1t, su3_matrix **link ) {

	// #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
			dirb[2] = {NODIR,NODIR}, 
			diro[2] = {NODIR,NODIR};
  register int ig;
	#define NGATHER 4
  register site *s;
  msg_tag *tag[NGATHER];
  for( ig=0; ig<NGATHER; ig++ ) { tag[ig] = NULL; }
  // #define NTEMP_STORAGE 6
  // su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;

  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
  //   su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  //   if(su3mat[ig] == NULL) {
  //       printf( "temporal_gauge_action_w_s_lwr_bulk: can't malloc su3mat[%d]\n", ig );
  //       fflush(stdout); terminate(1);
  //   }
  // }
  #define NTEMP 6
  su3_matrix **su3mat = new_field( NTEMP );

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

	dir[0]=TUP; {
		for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      setup_blocked_dirs( dir, dirb ,diro );
			#define LINK0 link[dir[0]]
			#define LINK1 link[dir[1]]

			ig = 0;
			// request link[dir[1]] from direction dir[0], not blocked in time direction
			tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
																		dir[0], EVENANDODD, gen_pt[ig]);

			ig = 1;
			// request link[dir[0]] from direction dirb[1]
			tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
																		dirb[1], EVENANDODD, gen_pt[ig]);

      // while waiting for gathers multiply two two links on the site
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
  			mult_su3_an( &(LINK1[i]), &(LINK0[i]), &(su3mat[0][i]) );

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
      {
				// su3_adjoint( &(LINK1[i]), &(su3mat[0][i]) );
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[1]);

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
      {
        mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        mult_su3_na( &tempmat, (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
        // mult_su3_na( &(s->link[dir[1]]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

			ig = 2;
			// request staple for "left" = su3mat[5] from dirb[1]
			tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
																		dirb[1], EVENANDODD, gen_pt[ig]);

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "bottom" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
      // must include here one time step more, since su3mat[4] is gathered down one time step
 			IF_ACTIVE(s) 
      { 
        mult_su3_na( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        mult_su3_na( &(su3mat[1][i]), &(LINK0[i]), &(su3mat[4][i]) );
        // su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        // su3mat_copy( &(su3mat[1][i]), &(su3mat[4][i]) );
      }

			ig = 3;
			// request staple for "bottom" = su3mat[4] from dir[0], not blocked in time direction
			tag[ig] = start_gather_field( su3mat[4], sizeof(su3_matrix),
																		dir[0], EVENANDODD, gen_pt[ig]);
      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
      {
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // su3_adjoint( &(LINK1[i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dirb[1]
        // to the accumulator
				*wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
      {
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
				if ( s->t < UPR_BDRY-1 )
          *wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
				*wl1x1t += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      if ( tag[ig] != NULL ) 
          cleanup_gather(tag[ig]);

	    #undef LINK0
	    #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );
  g_doublesum( wl1x2t );
  g_doublesum( wl2x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1. / LWR_BULK_LEN;
  *wl1x1t /= 3;
  *wl1x2t /= volume;
  *wl1x2t *= block_stride * block_stride * block_stride;
  *wl1x2t *= nt * 1. / ( LWR_BULK_LEN - 1.);
  *wl1x2t /= 3;
	*wl2x1t /= volume;
  *wl2x1t *= block_stride * block_stride * block_stride;
  *wl2x1t *= nt * 1. / LWR_BULK_LEN;
  *wl2x1t /= 3;

  // deallocate temporary storage
  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
  //   free( su3mat[ig] );
	// #undef NTEMP_STORAGE
  destroy_field( &su3mat );
  #undef NTEMP
	#undef NGATHER
}

/* Compute loops, here only temporal part,
   1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   restricted to the lower bulk */
static void
dropped_temporal_gauge_action_w_s_lwr_bulk( double *wl1x1t, 
  double *wl1x2t, double *wl2x1t, su3_matrix **link ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  #define NGATHER 4
  register site *s;
  msg_tag *tag[NGATHER];
  for( ig=0; ig<NGATHER; ig++ ) { tag[ig] = NULL; }
  // #define NTEMP_STORAGE 6
  // su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;

  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
  //   su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  //   if(su3mat[ig] == NULL) {
  //       printf( "dropped_temporal_gauge_action_w_s_lwr_bulk: can't malloc su3mat[%d]\n", ig );
  //       fflush(stdout); terminate(1);
  //   }
  // }

  #define NTEMP 6
  su3_matrix **su3mat = new_field( NTEMP );

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;


  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      setup_blocked_dirs( dir, dirb ,diro );
      #define LINK0 link[dir[0]]
      #define LINK1 link[dir[1]]

      ig = 0;
      // request link[dir[1]] from direction dir[0], not blocked in time direction
      tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);

      // ig = 1;
      // // request link[dir[0]] from direction dirb[1]
      // tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
      //                              dirb[1], EVENANDODD, gen_pt[ig]);

      // while waiting for gathers multiply two two links on the site
      // FORALLSITES(i, s) 
      // IF_BLOCKED(s, block_stride)  
      // IF_LOWER_BULK(s) 
      //  mult_su3_an( &(s->link[dir[1]]), &(s->link[dir[0]]), &(su3mat[0][i]) );

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        su3_adjoint( &(LINK1[i]), &(su3mat[0][i]) );
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );
      }

      // wait_gather(tag[1]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        // mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        // mult_su3_na( &tempmat, (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
        mult_su3_na( &(LINK1[i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

      ig = 2;
      // request staple for "left" = su3mat[5] from dirb[1]
      tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
                                    dirb[1], EVENANDODD, gen_pt[ig]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "bottom" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      // must include here one time step more, since su3mat[4] is gathered down one time step
      IF_ACTIVE(s) 
      // if ( s->x == 0 && s->y == 0 && s->z == 0 )
      { 
        // mult_su3_na( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        // mult_su3_na( &(su3mat[1][i]), &(s->link[dir[0]]), &(su3mat[4][i]) );
        su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        su3mat_copy( &(su3mat[1][i]), &(su3mat[4][i]) );
      }

      ig = 3;
      // request staple for "bottom" = su3mat[4] from dir[0], not blocked in time direction
      tag[ig] = start_gather_field( su3mat[4], sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);
      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        // form a staple for "top" link in the plaquette
        // mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        su3_adjoint( &(LINK1[i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dirb[1]
        // to the accumulator
        *wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
        if ( s->t < UPR_BDRY-1 )
          *wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
        *wl1x1t += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
        if ( tag[ig] != NULL ) 
          cleanup_gather(tag[ig]);

      #undef LINK0
      #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );
  g_doublesum( wl1x2t );
  g_doublesum( wl2x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1. / LWR_BULK_LEN;
  *wl1x1t /= 3;
  *wl1x2t /= volume;
  *wl1x2t *= block_stride * block_stride * block_stride;
  *wl1x2t *= nt * 1. / ( LWR_BULK_LEN - 1.);
  *wl1x2t /= 3;
  *wl2x1t /= volume;
  *wl2x1t *= block_stride * block_stride * block_stride;
  *wl2x1t *= nt * 1. / LWR_BULK_LEN;
  *wl2x1t /= 3;

  // deallocate temporary storage
  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
  //   free( su3mat[ig] );
  // #undef NTEMP_STORAGE
  destroy_field( &su3mat );
  #undef NTEMP
  #undef NGATHER
}

/* Compute loops, here only temporal part,
   1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action */
static void
temporal_gauge_action_w_s_full( double *wl1x1t, 
  double *wl1x2t, double *wl2x1t, su3_matrix **link ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  #define NGATHER 4
  register site *s;
  msg_tag *tag[NGATHER];
  for( ig=0; ig<NGATHER; ig++ ) { tag[ig] = NULL; }
  // #define NTEMP_STORAGE 6
  // su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;

  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
  //   su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  //   if(su3mat[ig] == NULL) {
  //       printf( "temporal_gauge_action_w_s_full: can't malloc su3mat[%d]\n", ig );
  //       fflush(stdout); terminate(1);
  //   }
  // }

  #define NTEMP 6
  su3_matrix **su3mat = new_field( NTEMP );

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      setup_blocked_dirs( dir, dirb ,diro );
      #define LINK0 link[dir[0]]
      #define LINK1 link[dir[1]]

      ig = 0;
      // request link[dir[1]] from direction dir[0], not blocked in time direction
      tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);

      ig = 1;
      // request link[dir[0]] from direction dirb[1]
      tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                    dirb[1], EVENANDODD, gen_pt[ig]);

      // while waiting for gathers multiply two two links on the site
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
        mult_su3_an( &(LINK1[i]), &(LINK0[i]), &(su3mat[0][i]) );

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      {
        // su3_adjoint( &(LINK1[i]), &(su3mat[0][i]) );
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[1]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      {
        mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[1][i]), &tempmat );
        mult_su3_na( &tempmat, (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
        // mult_su3_na( &(s->link[dir[1]]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

      ig = 2;
      // request staple for "left" = su3mat[5] from dirb[1]
      tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
                                    dirb[1], EVENANDODD, gen_pt[ig]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "bottom" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      { 
        mult_su3_na( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        mult_su3_na( &(su3mat[1][i]), &(LINK0[i]), &(su3mat[4][i]) );
        // su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &(su3mat[1][i]) );
        // su3mat_copy( &(su3mat[1][i]), &(su3mat[4][i]) );
      }

      ig = 3;
      // request staple for "bottom" = su3mat[4] from dir[0], not blocked in time direction
      tag[ig] = start_gather_field( su3mat[4], sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);
      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      {
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // su3_adjoint( &(LINK1[i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dirb[1]
        // to the accumulator
        *wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      {
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
        *wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
        *wl1x1t += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
        if ( tag[ig] != NULL ) 
          cleanup_gather(tag[ig]);

      #undef LINK0
      #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );
  g_doublesum( wl1x2t );
  g_doublesum( wl2x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t /= 3;
  *wl1x2t /= volume;
  *wl1x2t *= block_stride * block_stride * block_stride;
  *wl1x2t /= 3;
  *wl2x1t /= volume;
  *wl2x1t *= block_stride * block_stride * block_stride;
  *wl2x1t /= 3;

  // deallocate temporary storage
  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
  //   free( su3mat[ig] );
  // #undef NTEMP_STORAGE
  destroy_field( &su3mat );
  #undef NTEMP
  #undef NGATHER
}


/* Compute loops: 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   temporal and spatial part separately,
   in the active bulk including the two boundaries */
void
gauge_action_w_s_bulk( double wlt[2], double wls[2] ) {
	register int i, dir;
	register site *s;
  // #undef DROP_TIME_LINKS
  #ifndef DROP_TIME_LINKS
	#define NLINK 4
  #else
  #define NLINK 3
  #endif
	su3_matrix **link = new_links_from_site( ACTIVE, NLINK );

	spatial_gauge_action_w_s_region( ACTIVE, &(wls[0]), &(wls[1]), link );
	// spatial_gauge_action_w_s_region( LOWER_BULK, wl1x1s, wl1x2s, link );

	double wl2x1t[0];
  #ifndef DROP_TIME_LINKS
    temporal_gauge_action_w_s_lwr_bulk( &(wlt[0]), &(wlt[1]), wl2x1t, link );
  #else
    dropped_temporal_gauge_action_w_s_lwr_bulk( &(wlt[0]), &(wlt[1]), wl2x1t, link );
  #endif
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}


/* Compute loops: 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   temporal and spatial part separately,
   in the full 4D volume */
void
gauge_action_w_s_full( double wlt[2], double wls[2] ) {

  register int i, dir;
  register site *s;
  #define NLINK 4
  su3_matrix **link = new_links_from_site( FULLVOL, NLINK );

  spatial_gauge_action_w_s_region( FULLVOL, &(wls[0]), &(wls[1]), link );

  double wl2x1t[0];
  temporal_gauge_action_w_s_full( &(wlt[0]), &(wlt[1]), wl2x1t, link );
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}


#ifdef SPHALERON

void
fmunu_fmunu_bdry( su3_matrix **link_last_flow, 
                  double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;
  // msg_tag *tag0 = NULL, *tag1 = NULL, *tag2 = NULL;

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fsi = NULL;
  su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 2 * sizeof(double) );

  /* Compute 8*F_mu,nu at each site on the boundary slices */  
  #define NFS (5 * N_LAST_FLOW)
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_bdry( F_OFFSET(link), link_last_flow, fstrength );
  
  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_BOUNDARY(s) 
  {
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XY + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_BOUNDARY(s) 
    {
      ft = &(fstrength[FS_ZT + 1 * N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XY + 1 * N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XY + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_BOUNDARY(s) 
    {
      ft = &(fstrength[FS_YT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
      IF_LWR_BDRY(s) {
        charge[0] -= realtrace_su3(fs, ft); 
      } else {
        charge[0] += realtrace_su3(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= realtrace_su3(fsi, ft);
        } else {
          charge[1] += realtrace_su3(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_YZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_BOUNDARY(s) 
    {
      ft = &(fstrength[FS_XT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_YZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_YZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }
  }

  // cleanup_gather( tag2 );
  // cleanup_gather( tag1 );
  // cleanup_gather( tag0 );

  /* Sum over all nodes */
  g_vecdoublesum( time,   1 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 2 );

  /* Norma(lizations */
  for ( i = 0; i < 1 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    time[i] /= 0.25; /* extended only forward in time, factor 1/2, then squared */
    time[i] /= 2. / nt; /* computed only on two time slice */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    space[i] /= 2. / nt; /* computed only on two time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    // charge[i] /= 0.5; /* extended only forward in time */
  }
  time[1] = 0.;
  charge[2] += charge[0];
  charge[3] += charge[1];

  destroy_field( &fstrength );
  #undef NFS
}

void
fmunu_fmunu_lwr_bdry( su3_matrix **link_last_flow, 
                  double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;
  // msg_tag *tag0 = NULL, *tag1 = NULL, *tag2 = NULL;

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fsi = NULL;
  su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 2 * sizeof(double) );

  /* Compute 8*F_mu,nu at each site on the boundary slices */  
  #define NFS (5 * N_LAST_FLOW)
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_lwr_bdry( F_OFFSET(link), link_last_flow, fstrength );
  
  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_LWR_BDRY(s) 
  {
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XY + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_LWR_BDRY(s) 
    {
      ft = &(fstrength[FS_ZT + 1 * N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XY + 1 * N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XY + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_LWR_BDRY(s) 
    {
      ft = &(fstrength[FS_YT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
      IF_LWR_BDRY(s) {
        charge[0] -= realtrace_su3(fs, ft); 
      } else {
        charge[0] += realtrace_su3(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= realtrace_su3(fsi, ft);
        } else {
          charge[1] += realtrace_su3(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_YZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_LWR_BDRY(s) 
    {
      ft = &(fstrength[FS_XT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_YZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_YZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }
  }

  // cleanup_gather( tag2 );
  // cleanup_gather( tag1 );
  // cleanup_gather( tag0 );

  /* Sum over all nodes */
  g_vecdoublesum( time,   1 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 2 );

  /* Norma(lizations */
  for ( i = 0; i < 1 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    time[i] /= 0.25; /* extended only forward in time, factor 1/2, then squared */
    time[i] /= 1. / nt; /* computed only on one time slice */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    space[i] /= 1. / nt; /* computed only on one time slice */
  }
  for ( i = 0; i < 2 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    // charge[i] /= 0.5; /* extended only forward in time */
  }
  time[1] = 0.;
  charge[2] += charge[0];
  charge[3] += charge[1];

  destroy_field( &fstrength );
  #undef NFS
}

void
fmunu_fmunu_upr_bdry( su3_matrix **link_last_flow, 
                  double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;
  // msg_tag *tag0 = NULL, *tag1 = NULL, *tag2 = NULL;

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fsi = NULL;
  su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 2 * sizeof(double) );

  /* Compute 8*F_mu,nu at each site on the boundary slices */  
  #define NFS (5 * N_LAST_FLOW)
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_upr_bdry( F_OFFSET(link), link_last_flow, fstrength );
  
  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_UPR_BDRY(s) 
  {
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XY + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_UPR_BDRY(s) 
    {
      ft = &(fstrength[FS_ZT + 1 * N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XY + 1 * N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XY + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_UPR_BDRY(s) 
    {
      ft = &(fstrength[FS_YT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_XZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      /* ANTICYCLIC! ReTr{ fs.dag * ft } */
      IF_LWR_BDRY(s) {
        charge[0] -= realtrace_su3(fs, ft); 
      } else {
        charge[0] += realtrace_su3(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_XZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= realtrace_su3(fsi, ft);
        } else {
          charge[1] += realtrace_su3(fsi, ft);
        }
      #endif
    }

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_YZ + 3 * N_LAST_FLOW][i]);
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    // IF_UPR_BDRY(s) 
    {
      ft = &(fstrength[FS_XT+N_LAST_FLOW][i]);
      time[0] -= real_trace_nn(ft, ft);

      add_su3_matrix( fs, &(fstrength[FS_YZ+N_LAST_FLOW][i]), &tmp );
      fs = &tmp;
      IF_LWR_BDRY(s) {
        charge[0] -= real_trace_nn(fs, ft);
      } else {
        charge[0] += real_trace_nn(fs, ft);
      }
      #if NFS > 12
        add_su3_matrix( fsi, &(fstrength[FS_YZ + 4 * N_LAST_FLOW][i]), &tmp );
        fsi = &tmp;
        IF_LWR_BDRY(s) {
          charge[1] -= real_trace_nn(fsi, ft);
        } else {
          charge[1] += real_trace_nn(fsi, ft);
        }
      #endif
    }
  }

  // cleanup_gather( tag2 );
  // cleanup_gather( tag1 );
  // cleanup_gather( tag0 );

  /* Sum over all nodes */
  g_vecdoublesum( time,   1 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 2 );

  /* Norma(lizations */
  for ( i = 0; i < 1 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    time[i] /= 0.25; /* extended only forward in time, factor 1/2, then squared */
    time[i] /= 1. / nt; /* computed only on one time slice */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    space[i] /= 1. / nt; /* computed only on one time slice */
  }
  for ( i = 0; i < 2 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    // charge[i] /= 0.5; /* extended only forward in time */
  }
  time[1] = 0.;
  charge[2] += charge[0];
  charge[3] += charge[1];

  destroy_field( &fstrength );
  #undef NFS
}

/* Computes the field strength components and topological charge,
 * at half-integer time steps within the bulk */
void
fmunu_fmunu_half( double *time, double *space, double *charge ) {
  /* Site variables */
  register int i;
  register site *s;

  /* Temporary component storage */
  // su3_matrix *ft, *fs;
  su3_matrix *ft = NULL, *fs = NULL, *fsi = NULL;
  su3_matrix tmp;

  /* Initialize sums */
  memset( time,   '\0', 2 * sizeof(double) );
  memset( space,  '\0', 2 * sizeof(double) );
  memset( charge, '\0', 4 * sizeof(double) );

  // make_field_strength_bulk( F_OFFSET(link), F_OFFSET(fieldstrength) );
  /* Compute 8*F_mu,nu at each site for half-integer time */

  /* Compute 8*F_mu,nu at each site for half-integer time */
  #define NFS 12
  su3_matrix **fstrength = new_field( NFS );
  make_fieldstrength_half( F_OFFSET(link), fstrength );

  /* Loop over each site to sum F_mu,nu components */
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
 	IF_LOWER_BULK(s) 
 	{
    fs = &(fstrength[FS_XY][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XY+9][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    ft = &(fstrength[FS_ZT+3][i]);
    add_su3_matrix( ft, &(fstrength[FS_ZT+0][i]), &tmp );
    ft = &tmp;
    time[0] -= real_trace_nn(ft, ft);
    charge[0] -= real_trace_nn(fs, ft);
    #if NFS > 9
      charge[1] -= real_trace_nn(fsi, ft);
    #endif

    fs = &(fstrength[FS_XZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_XZ+9][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
    ft = &(fstrength[FS_YT+3][i]);
    add_su3_matrix( ft, &(fstrength[FS_YT+0][i]), &tmp );
    ft = &tmp;
    time[0] -= real_trace_nn(ft, ft);
    /* ANTICYCLIC! ReTr{ fs.dag * ft } */
    charge[0] -= realtrace_su3(fs, ft);
    #if NFS > 9
      charge[1] -= realtrace_su3(fsi, ft);
    #endif

    fs = &(fstrength[FS_YZ][i]);
    space[0] -= real_trace_nn(fs, fs);
    #if NFS > 9
      fsi = &(fstrength[FS_YZ+9][i]); // improved
      space[1] -= real_trace_nn(fsi, fsi);
    #endif
	  ft = &(fstrength[FS_XT+3][i]);
    add_su3_matrix( ft, &(fstrength[FS_XT+0][i]), &tmp );
    ft = &tmp;
    time[0] -= real_trace_nn(ft, ft);
    charge[0] -= real_trace_nn(fs, ft);
    #if NFS > 9
      charge[1] -= real_trace_nn(fsi, ft);
    #endif
  }

  /* Sum over all nodes */
  g_vecdoublesum( time,   1 );
  g_vecdoublesum( space,  2 );
  g_vecdoublesum( charge, 2 );

  /* Norma(lizations */
  for ( i = 0; i < 1 ; i++ ) {
    time[i] /= ( volume * 64.0 );
    time[i] *= block_stride * block_stride * block_stride;
    time[i] /= 0.25;  /* extended only forward in time, factor 1/2, then squared */
    time[i] /= LWR_BULK_LEN * 1. / nt; /* computed only on nt/2 integer time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    space[i] /= ( volume * 64.0 );
    space[i] *= block_stride * block_stride * block_stride;
    space[i] /= LWR_BULK_LEN * 1. / nt; /* computed only on nt/2+1 integer time slices */
  }
  for ( i = 0; i < 2 ; i++ ) {
    charge[i] *= 0.0003957858736028819197; /* normalization of 1/(8^2 * 4 * PI^2) */
    charge[i] /= 0.5; /* extended only forward in time */
  }
  time[1] = charge[2] = charge[3] = 0.;

  destroy_field( &fstrength );
  #undef NFS
}

/* Compute loops, here only temporal part,
   1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   using half-integer time steps */
static void
dropped_temporal_gauge_action_w_s_half( double *wl1x1t, 
  double *wl1x2t, double *wl2x1t, su3_matrix **link, su3_matrix **link_half ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  #define NGATHER 4
  register site *s;
  msg_tag *tag[NGATHER];
  for( ig=0; ig<NGATHER; ig++ ) { tag[ig] = NULL; }
  // #define NTEMP_STORAGE 6
  // su3_matrix *su3mat[NTEMP_STORAGE];
  su3_matrix tempmat;
  double tt;

  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
  //   su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  //   if(su3mat[ig] == NULL) {
  //       printf( "gauge_action_w_s_bulk: can't malloc su3mat[%d]\n", ig );
  //       fflush(stdout); terminate(1);
  //   }
  // }
  #define NTEMP 6
  su3_matrix **su3mat = new_field( NTEMP );

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      setup_blocked_dirs( dir, dirb ,diro );
      #define LINKF link[dir[1]]
      #define LINKH link_half[dir[1]]

      ig = 0;
      // request link[dir[1]] at integer time from direction dir[0], not blocked in time direction
      tag[ig] = start_gather_field( LINKF , sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);

      wait_gather(tag[0]);

      ig = 1;
      // request link[dir[1]] at half-int time from direction dir[0], not blocked in time direction
      tag[ig] = start_gather_field( LINKH , sizeof(su3_matrix),
                                    dir[0], EVENANDODD, gen_pt[ig]);

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        mult_su3_an( &(LINKH[i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[1]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        mult_su3_na( &(LINKH[i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[5][i]) );
      }

      ig = 2;
      // request staple at half-int to int time for "left" = su3mat[5] from dirb[1]
      tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
                                    dirb[1], EVENANDODD, gen_pt[ig]);

      //fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      /* reinterpret here what is being done with su3mat[0] and su3mat[4] 
         in order to average the rectangles of int to half-int and 
         half-int to int */
      IF_ACTIVE(s) 
      { 
        mult_su3_an( &(LINKF[i]), &(LINKH[i]), &(su3mat[0][i]) );
        mult_su3_na( &(LINKF[i]), &(LINKH[i]), &(su3mat[4][i]) );
      }

      ig = 3;
      // request staple for "bottom" = su3mat[4] from dir[0], not blocked in time direction
      // tag[ig] = start_gather_field( su3mat[4], sizeof(su3_matrix),
      //                               dir[0], EVENANDODD, gen_pt[ig]);

      // request staple at int to half-int time for "left" = su3mat[5] from dirb[1],
      // but that combination yields exactly the same as the previous one...
      // tag[ig] = start_gather_field( su3mat[4], sizeof(su3_matrix),
      //                               dirb[1], EVENANDODD, gen_pt[ig]);
      wait_gather(tag[2]);

      // wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        // form a staple for "top" link in the plaquette
        su3_adjoint( &(LINKF[i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dirb[1]
        // to the accumulator
        /* The 2x1t loops from int to half-int yields the same as the 
           the ones from half-int to int, but is not defined on right 
           beneath the upper boundary. Weights included at the end. */
        *wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
        // *wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[0][i]) );
      }

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LOWER_BULK(s) 
      {
        // get the contribution of 1x2 rectangle extended in dir[0]
        *wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[0][i]), &(LINKF[i]) );
        if ( s->t < UPR_BDRY-1 )
          *wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[1][i]), &(LINKH[i]) );

        // and of 1x1 plaquette to the accumulators
        /* The 1x1t loop from int to half-int yields the same as the 
           the one from half-int to int, but is not defined on right 
           beneath the upper boundary. The weights account for this. */
        *wl1x1t += ( s->t < UPR_BDRY-1 ? 2. : 1. ) * 
          realtrace_su3( (su3_matrix *)(gen_pt[0][i]), &(LINKH[i]) );
        // if ( s->t < UPR_BDRY-1 )
        //   *wl1x1t += realtrace_su3( &(LINKF[i]), &(LINKH[i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ ) 
        if ( tag[ig] != NULL ) 
          cleanup_gather(tag[ig]); 

      #undef LINKF
      #undef LINKH
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );
  g_doublesum( wl1x2t );
  g_doublesum( wl2x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1. / ( LWR_BULK_LEN * 2. - 1.);
  *wl1x1t /= 3;
  *wl1x2t /= volume;
  *wl1x2t *= block_stride * block_stride * block_stride;
  *wl1x2t *= nt * 1. / ( LWR_BULK_LEN * 2. - 1.);
  *wl1x2t /= 3;
  *wl2x1t /= volume;
  *wl2x1t *= block_stride * block_stride * block_stride;
  *wl2x1t *= nt * 1. / ( LWR_BULK_LEN * 1. );
  *wl2x1t /= 3;

  // deallocate temporary storage
  // for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
  //   free( su3mat[ig] );
  // #undef NTEMP_STORAGE
  destroy_field( &su3mat );
  #undef NTEMP
  #undef NGATHER
}

static void
dropped_temporal_gauge_action_w_s_bdry( double *wl1x1t, 
                                double *wl1x2t, double *wl2x1t, 
                                su3_matrix **link, 
                                su3_matrix **link_lf ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR};
  register int ig;
  register site *s;

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      #define LINK0 link[dir[1]]
      #define LINK1 link_lf[dir[1]]

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_BOUNDARY(s) {
        // get the contribution of 1x1 plaquette to the accumulators
        *wl1x1t += 
          realtrace_su3( &(LINK1[i]), &(LINK0[i]) );
      }

      #undef LINK0
      #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1. / 2.;
  *wl1x1t /= 3;

}

static void
dropped_temporal_gauge_action_w_s_lwr_bdry( double *wl1x1t, 
                                double *wl1x2t, double *wl2x1t, 
                                su3_matrix **link, 
                                su3_matrix **link_lf ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR};
  register int ig;
  register site *s;

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      #define LINK0 link[dir[1]]
      #define LINK1 link_lf[dir[1]]

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_LWR_BDRY(s) {
        // get the contribution of 1x1 plaquette to the accumulators
        *wl1x1t += 
          realtrace_su3( &(LINK1[i]), &(LINK0[i]) );
      }

      #undef LINK0
      #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1.;
  *wl1x1t /= 3;

}

static void
dropped_temporal_gauge_action_w_s_upr_bdry( double *wl1x1t, 
                                double *wl1x2t, double *wl2x1t, 
                                su3_matrix **link, 
                                su3_matrix **link_lf ) {

  // #undef DROP_TIME_LINKS
  register int i;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR};
  register int ig;
  register site *s;

  // prepare accumulators
  *wl1x1t = 0;
  *wl1x2t = 0;
  *wl2x1t = 0;

  dir[0]=TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
      #define LINK0 link[dir[1]]
      #define LINK1 link_lf[dir[1]]

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride) 
      IF_UPR_BDRY(s) {
        // get the contribution of 1x1 plaquette to the accumulators
        *wl1x1t += 
          realtrace_su3( &(LINK1[i]), &(LINK0[i]) );
      }

      #undef LINK0
      #undef LINK1
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( wl1x1t );

  // get densities
  *wl1x1t /= volume;
  *wl1x1t *= block_stride * block_stride * block_stride;
  *wl1x1t *= nt * 1.;
  *wl1x1t /= 3;

}

void
gauge_action_w_s_bdry ( su3_matrix **link_last_flow, 
                        double wlt[2], double wls[2] ) {

  register int i, dir;
  register site *s;
  #define NSPATIAL 3
  su3_matrix **link = new_links_from_site( BOUNDARY, NSPATIAL );

  spatial_gauge_action_w_s_region( BOUNDARY, &(wls[0]), &(wls[1]), link_last_flow );

  double wl2x1t[0];
  dropped_temporal_gauge_action_w_s_bdry( &(wlt[0]), &(wlt[1]), wl2x1t, link, link_last_flow );
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NSPATIAL
}

void
gauge_action_w_s_lwr_bdry ( su3_matrix **link_last_flow, 
                        double wlt[2], double wls[2] ) {

  register int i, dir;
  register site *s;
  #define NSPATIAL 3
  su3_matrix **link = new_links_from_site( LOWER_BOUNDARY, NSPATIAL );

  spatial_gauge_action_w_s_region( LOWER_BOUNDARY, &(wls[0]), &(wls[1]), link_last_flow );

  double wl2x1t[0];
  dropped_temporal_gauge_action_w_s_lwr_bdry( &(wlt[0]), &(wlt[1]), wl2x1t, link, link_last_flow );
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NSPATIAL
}

void
gauge_action_w_s_upr_bdry ( su3_matrix **link_last_flow, 
                        double wlt[2], double wls[2] ) {

  register int i, dir;
  register site *s;
  #define NSPATIAL 3
  su3_matrix **link = new_links_from_site( UPPER_BOUNDARY, NSPATIAL );

  spatial_gauge_action_w_s_region( UPPER_BOUNDARY, &(wls[0]), &(wls[1]), link_last_flow );

  double wl2x1t[0];
  dropped_temporal_gauge_action_w_s_upr_bdry( &(wlt[0]), &(wlt[1]), wl2x1t, link, link_last_flow );
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NSPATIAL
}

/* Compute loops: 1x1 -- plaquette, 1x2 + 2x1 -- rectangle
   for Wilson (one-plaquette) and
   Symanzik tree-level (plaquette and rectangle) action,
   temporal and spatial part separately,
   in the lower bulk including the two boundaries,
   at half-integer time steps */
void
gauge_action_w_s_half( double wlt[2], double wls[2] ) {

  #define NSPATIAL 3
  su3_matrix **link = new_links_from_site( ACTIVE, NSPATIAL);
	su3_matrix **link_half = new_half_links();

	spatial_gauge_action_w_s_region( LOWER_BULK, &(wls[0]), &(wls[1]), link_half );

  // temporal action still not corrected
  double wl2x1t[0];
  dropped_temporal_gauge_action_w_s_half( &(wlt[0]), &(wlt[1]), wl2x1t, link, link_half );
  wlt[1] = ( wlt[1] + *wl2x1t ) * 0.5;

  destroy_field( &link );
	destroy_field( &link_half );
  #undef NSPATIAL
}

#endif