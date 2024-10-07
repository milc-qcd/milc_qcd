/************************* blocking.c **********************/
/* Utilities associated with blocking the spatial lattice */

#include "wilson_flow_includes.h"
#include "defines.h"
#include <string.h>


void setup_blocked_dirs( int *dir, int *dirb, int *diro ) {

	#ifdef BLOCKING
	  switch( block_stride ) {
	    case (1):
	      dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
	      dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
	      break;
	    case (2):
	      dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
	      dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
	      break;
	    case (4):
	      dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
	      dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
	      break;
	    case (8):
	      dirb[0] = DIR8(dir[0]); diro[0] = DIR8(OPP_DIR(dir[0]));
	      dirb[1] = DIR8(dir[1]); diro[1] = DIR8(OPP_DIR(dir[1]));
	      break;
	    default:
	      node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
	      terminate(1);
	      break;
	  }
	#else
	  dirb[0] = dir[0]; diro[0] = OPP_DIR(dir[0]);
    dirb[1] = dir[1]; diro[1] = OPP_DIR(dir[1]);
	#endif
  if ( dir[0] == TUP ) {
    dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
  }
  if ( dir[1] == TUP ) {
    dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
  }
}

#ifdef BLOCKING

void spatial_blocking( void ) {
	register int i;
	int dir, dirb = NODIR, icur;
	int stride, new_stride = 2 * block_stride;
	register site *s = NULL;

	// #define NBLOCK_STORAGE 2
	// su3_matrix *temp[NBLOCK_STORAGE];
	msg_tag *tag = NULL;
  #define NBLOCK 2
  su3_matrix **temp = new_field( NBLOCK );	

	if ( nx % new_stride !=0 || ny % new_stride !=0 || nz % new_stride !=0 
		|| nx  < new_stride || ny  < new_stride || nz  < new_stride ) {
		node0_printf("Illegal BLOCKING option rejected: new block_stride %d not an even integer divisor\n",
			new_stride);
		fflush(stdout); return;
	}
	if ( new_stride > MAX_BLOCK_STRIDE ) {
		node0_printf("Requested BLOCKING option rejected: new block_stride %d exceeds MAX_BLOCK_STRIDE %d\n",
			new_stride, MAX_BLOCK_STRIDE);
		fflush(stdout); return;
	}

	node0_printf("BLOCK spatial links from STRIDE %d to STRIDE %d\n",block_stride,new_stride);

	// for( icur=0; icur<NBLOCK_STORAGE; icur++ ) {
  // 	temp[icur] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  // 	if( temp[icur] == NULL ) {
  //       printf( "spatial_blocking: can't malloc temp[%d]\n", icur );
  //       fflush(stdout); terminate(1);
  //   }
  // }

	// FORALLSITES(i,s) {
	// 	if ( s-> x % 2 == 0 && s-> y == 0 && s-> t==0 ) {
	// 	// if ( s-> t == 1 ) {
	// 		node0_printf("i %04d: s %02d %02d %02d %02d\n",i,s->x,s->y,s->z,s->t);
	// 		dumpmat( &(s->link[ZUP]) );
	// 	}
	// }

	FORALLUPDIRBUT(TUP,dir) {
		FORALLSITES(i,s)
		IF_BLOCKED(s, block_stride ) 
		{
			su3mat_copy( &(s->link[dir]), &(temp[0][i]) );
			// su3mat_copy( &(s->link[dir]), &(temp[2+dir][i]) );
		}

		for (stride = new_stride, icur = 0; stride > block_stride; stride--, icur = 1-icur ) 
		{
			dirb = ( stride > block_stride && block_stride > 1 ? DIR2(dir) : dir);
			// node0_printf("dir %d stride %d -> dirb %d icur %d\n", dir, stride, dirb, icur);
			tag = start_gather_field( temp[icur] , sizeof(su3_matrix),
																dirb, EVENANDODD, gen_pt[0]);

    	wait_gather( tag );

			if ( stride <= 2 + block_stride ) 
			{ // multiplication after the last gather
				FORALLSITES(i,s) 
				{
					mult_su3_nn( &(s->link[dir]), (su3_matrix*)(gen_pt[0][i]), &(temp[1-icur][i]) );
					// mult_su3_nn( &(s->link[dir]), (su3_matrix*)(gen_pt[0][i]), &(temp[5+dir][i]) );
				}
				stride--;
			} else {// move through this stage without multiplication
				if ( stride % 2 == 0 ) 
					stride--;
				FORALLSITES(i,s) 
					su3mat_copy( (su3_matrix*)(gen_pt[0][i]), &(temp[1-icur][i]) );
			}

			cleanup_gather( tag ); 
		}
		// if (dir == XUP) {
		// 	node0_printf("stride %d dirb-dir %d icur %d post\n",stride, dirb-dir,icur);
		// 	dumpmat( &(temp[1-icur][0]) );
		// }
		FORALLSITES(i,s) {
			IF_BLOCKED(s, new_stride ) {
				su3mat_copy( &(temp[icur][i]), &(s->link[dir]) );
			} else {
				set_identity( &(s->link[dir]) );
			}
		}
	}

	// FORALLSITES(i,s) {
	// 	if ( s-> x % 2 == 0 && s-> y == 0 && s-> t==0 ) {
	// 	// if ( s-> t == 1 ) {
	// 		node0_printf("i %04d: s %02d %02d %02d %02d\n",i,s->x,s->y,s->z,s->t);
	// 		dumpmat( &(temp[2+ZUP][i]) );
	// 		dumpmat( &(s->link[ZUP]) );
	// 		dumpmat( &(temp[5+ZUP][i]) );
	// 	}
	// }

	// deallocate temporary storage
  // for( icur=icur; icur<NBLOCK_STORAGE; icur++ ) {
  //   free( temp[icur] );
  // }
	// #undef NBLOCK_STORAGE
	 destroy_field( &temp );
  #undef NBLOCK
  block_stride = new_stride;
}

#ifdef DEBUG_BLOCKING
void test_blocking ( void ) {

	#undef DROP_TIME_LINKS

  #define NLINK 4
	#define NTEMP_STORAGE 6
	#define NGATHER 4
	#define NBLOCK_STORAGE (2*NGATHER)
  #define NLOOP_AXB 18

  register int i, idx;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
			dirb[2] = {NODIR,NODIR}, 
			diro[2] = {NODIR,NODIR};
  register int ig;
  int disp[4];
  int icur[NGATHER], ncur[NGATHER];
  double maxdev;
  register site *s;
  msg_tag *tag[2*NGATHER];
  su3_matrix tempmat;

  node0_printf("\nTest the Blocking starting from unblocked links: form plaquettes\n\
and rectangles with two steps instead of one for spatial links\n\
then block, and use the standard code; compare the untraced loops\n\n");
  int my_i=0;

  su3_matrix *su3mat[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    su3mat[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    memset( su3mat[ig], '\0', sizeof(su3_matrix)*sites_on_node );
    if(su3mat[ig] == NULL) {
        printf( "test_blocking: can't malloc su3mat[%d]\n", ig );
        fflush(stdout); terminate(1);
    }
  }

  su3_matrix *loop_axb[NLOOP_AXB];
  for( ig=0; ig<NLOOP_AXB; ig++ ) {
    loop_axb[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    memset( loop_axb[ig], '\0', sizeof(su3_matrix)*sites_on_node );
    if(loop_axb[ig] == NULL) {
        printf( "test_blocking: can't malloc loop_axb[%d]\n", ig );
        fflush(stdout); terminate(1);
    }
  }

  su3_matrix *temp[NBLOCK_STORAGE];
  msg_tag *btag;
  if ( block_stride >= 1 ) {
		for( ig=0; ig<NBLOCK_STORAGE; ig++ ) {
  		temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
  		memset( temp[ig], '\0', sizeof(su3_matrix)*sites_on_node );
  		if( temp[ig] == NULL ) {
        printf( "test_blocking: can't malloc temp[%d]\n", ig );
        fflush(stdout); terminate(1);
    	}
  	}
	}
	su3_matrix *dlink_dir1_at_s00 = &(temp[ 0][0]);
	su3_matrix *dlink_dir1_at_s20 = &(temp[ 1][0]);
	su3_matrix *dlink_dir0_at_s00 = &(temp[ 2][0]);
	su3_matrix *dlink_dir0_at_s02 = &(temp[ 3][0]);

	if ( block_stride > 1 ) {
		node0_printf("test_blocking: must start from unblocked links, instead block_stride %d\n",block_stride);
		terminate(1);
	}

	//WHICH ONE IS COUNTED AS LEFT/RIGHT, AND WHICH ONE AS TOP/BOTTOM?
	// dir[0] THOUGHT AS VERTICAL, wrt to dir LEFT/RIGHT
	// dir[1]

	su3_matrix **link = new_links_from_site( FULLVOL, NLINK );

	int this_region = FULLVOL;
	double time_frac = 1.;
  switch ( this_region ) {
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
  	default:
  		node0_printf("Undefined REGION: %d\n",this_region);
  		terminate(1);
  		break;
  }

  // accumulators
  double wl1x1s = 0;
  double wl1x2s = 0;

	setup_blocked_dirs( dir, dirb, diro );
  // for( dir[0]=YUP; dir[0]<TUP; dir[0]++ ) {
  //   for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {

  //  		idx = ( dir[0] == ZUP ? dir[1] : ZUP) ;
	// 		switch( 2 * block_stride ) {
	// 			case (1):
	// 				dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
	// 				dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
	// 				break;
	// 			case (2):
	// 				dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
	// 				dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
	// 				break;
	// 			case (4):
	// 				dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
	// 				dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
	// 				break;
	// 			default:
	// 				node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
	// 				terminate(1);
	// 			break;
	// 		}
			#define LINK0 link[dir[0]]
			#define LINK1 link[dir[1]]

   	  for (ig = 0; ig < NGATHER; ig++) {
		  	icur[ig] = 2 * ig;
  			ncur[ig] = 4 * ig + 1;
  		}

  		ig = 0;
  		for ( i = 0; i < 4; i++) 
  			disp[i] = 0;
  		disp[dir[0]] = 2;
  		disp[dir[1]] = 1;

     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_REGION(s, this_region) 
     		su3mat_copy( &(LINK1[i]) , &(temp[icur[ig]][i]) );

      /* request link[dir[1]] from direction dirb[0], 
       * then extend via general_gather to make double dir[1]-link -> &(temp[icur[0]]),
       * finally gather link[dir[1]] from direction dir[1] 
       * to locally make a another double dir[1]-link later -> gen_pt[8] */
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig] - icur[ig] ) 
			if ( block_stride <=2 ) {
     		tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
       		     												dirb[0], EVENANDODD, gen_pt[ig]);
				btag = start_general_gather_field( LINK1 , sizeof(su3_matrix),
																			disp, EVENANDODD, gen_pt[8-ig]);
				wait_gather(tag[ig]);
				wait_general_gather( btag );
				FORALLSITES(i,s) 
				IF_BLOCKED(s, block_stride)	
				IF_REGION(s, this_region) 
					mult_su3_nn( (su3_matrix*)(gen_pt[ig][i]), (su3_matrix*)(gen_pt[8-ig][i]), 
											&(dlink_dir1_at_s20[i]) );
											// &(temp[ncur[ig]-icur[ig]][i]) );
				cleanup_general_gather( btag );
				tag[8-ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
																			dir[1], EVENANDODD, gen_pt[8-ig]);
			}

			ig = 1;
  		for ( i = 0; i < 4; i++) 
  			disp[i] = 0;
  		disp[dir[0]] = 1;
  		disp[dir[1]] = 2;

      FORALLSITES(i,s)
      IF_BLOCKED(s, block_stride)
      IF_REGION(s, this_region) 
      	su3mat_copy( &(LINK0[i]) , &(temp[icur[ig]][i]) );

      /* request link[dir[0]] from direction dirb[1], 
       * then extend via general_gather to make double dir[0]-link -> &(temp[icur[1]]),
       * finally gather link[dir[0]] from direction dir[0] 
       * to locally make a another double dir[0]-link later -> gen_pt[8-1]*/
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig]-icur[ig] )
			if ( block_stride <=2 ) {
     		tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
											       		    	dirb[1], EVENANDODD, gen_pt[ig]);
     		btag = start_general_gather_field( LINK0 , sizeof(su3_matrix),
																			disp, EVENANDODD, gen_pt[8-ig]);
				wait_gather(tag[ig]);
				wait_general_gather( btag );
				FORALLSITES(i,s) 
				IF_BLOCKED(s, block_stride)	
				IF_REGION(s, this_region)
					mult_su3_nn( (su3_matrix*)(gen_pt[ig][i]), (su3_matrix*)(gen_pt[8-ig][i]), 
											&(dlink_dir0_at_s02[i]) );
											// &(temp[ncur[ig]-icur[ig]][i]) );
				cleanup_gather( btag );
				tag[8-ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
																			dir[0], EVENANDODD, gen_pt[8-ig]);
			}

			wait_gather(tag[8-0]);

			// extend to double-link in dir[1] direction -> dlink_dir1_at_s00
			FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
  			mult_su3_nn( &(LINK1[i]), (su3_matrix*)(gen_pt[8-0][i]), &(dlink_dir1_at_s00[i]) );
			// cleanup_gather(tag[8-0]); 

			wait_gather(tag[8-1]);
			// extend to double-link in dir[0] direction -> dlink_dir0_at_s00
			FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
  			mult_su3_nn( &(LINK0[i]), (su3_matrix*)(gen_pt[8-1][i]), &(dlink_dir0_at_s00[i]) );
  		// cleanup_gather(tag[8-1]); 

      // while waiting for gathers multiply the two double-links on the site
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
  			mult_su3_an( &(dlink_dir1_at_s00[i]), &(dlink_dir0_at_s00[i]), &(su3mat[0][i]) );

      // form a double-staple for "right" double-link in the plaquette
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
        mult_su3_nn( &(su3mat[0][i]), &(dlink_dir1_at_s20[i]), &(su3mat[3][i]) );

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a double-staple for "left" double-link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        mult_su3_nn( &(dlink_dir1_at_s00[i]), &(dlink_dir0_at_s02[i]), &tempmat );
        mult_su3_na( &tempmat, &(dlink_dir1_at_s20[i]), &(su3mat[5][i]) );
      }

      ig = 2;
     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_REGION(s, this_region) 
     		su3mat_copy( &(su3mat[5][i]) , &(temp[icur[ig]][i]) );

      // request double-staple for "left" = su3mat[5] from dirb[1]
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig]-icur[ig] )
      if ( block_stride <=2 ) {
      	tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
      										  		    	dirb[1], EVENANDODD, gen_pt[ig]);
			}

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a double-staple for "bottom" double-link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        mult_su3_na( &(dlink_dir0_at_s02[i]), &(dlink_dir1_at_s20[i]), &(su3mat[1][i]) );
        mult_su3_na( &(su3mat[1][i]), &(dlink_dir0_at_s00[i]), &(su3mat[4][i]) );
      }

      ig = 3;
     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_REGION(s, this_region) 
     		su3mat_copy( &(su3mat[4][i]) , &(temp[icur[ig]][i]) );

      // request double-staple for "bottom" = su3mat[4] from dirb[0]
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig]-icur[ig] )
      if ( block_stride <=2 ) {
      	tag[ig] = start_gather_field( su3mat[4] , sizeof(su3_matrix),
      										  		     dirb[0], EVENANDODD, gen_pt[ig]);
			}

      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        mult_su3_an( &(dlink_dir0_at_s02[i]), &(su3mat[0][i]), &(su3mat[2][i]) );
			}

			wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        // mult_su3_an( &(dlink_dir0_at_s02[i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
        mult_su3_an( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]), &(loop_axb[idx+6][i]) );
      }

      // wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
				mult_su3_an( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]), &(loop_axb[idx+3][i]) );
				mult_su3_an( &(su3mat[1][i]), &(su3mat[0][i]), &(loop_axb[idx+0][i]) );
      }

      FORALLSITES(i, s) 
      IF_BLOCKED(s, 2 * block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        // mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
        wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      FORALLSITES(i, s) 
      IF_BLOCKED(s, 2 * block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
        wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
        wl1x1s += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      cleanup_gather(tag[ig]);
    } // dir[1]
  } // dir[0]

	// global sum
  g_doublesum( &wl1x1s );
  g_doublesum( &wl1x2s );

  // get densities
  wl1x1s /= volume;
  wl1x1s *= 2 * block_stride * 2 * block_stride * 2 * block_stride;
  wl1x1s *= time_frac;
  wl1x1s /= 3;
  wl1x2s /= volume;
  wl1x2s *= 2 * block_stride * 2 * block_stride * 2 * block_stride;
  wl1x2s *= time_frac;
  wl1x2s /= 6;

	su3_matrix *dlink_dir1_at_s10 = &(temp[ 1][0]);
	su3_matrix *slink_dir0_at_s02 = &(temp[ 3][0]);

  // accumulators
  double wl1x1t = 0;
  double wl1x2t = 0;
  double wl2x1t = 0;

  dir[0] = TUP; {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
   		idx = dir[1] + 9;

			setup_blocked_dirs( dir, dirb, diro );
      // dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
      // switch( 2 * block_stride ) {
      //   case (1):
      //     dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
      //     break;
      //   case (2):
      //     dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
      //     break;
      //   case (4):
      //     dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
      //     break;
      //   default:
      //     node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
      //     terminate(1);
      //   break;
      // }
      #define LINK0 link[dir[0]]
      #define LINK1 link[dir[1]]

   	  for (ig = 0; ig < NGATHER; ig++) {
		  	icur[ig] = 2 * ig;
  			ncur[ig] = 4 * ig + 1;
  		}

  		ig = 0;
  		for ( i = 0; i < 4; i++) 
  			disp[i] = 0;
  		disp[dir[0]] = 1;
  		disp[dir[1]] = 1;

     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_ACTIVE(s) 
     		su3mat_copy( &(LINK1[i]) , &(temp[icur[ig]][i]) );

			/* request link[dir[1]] from direction dir[0], not blocked in time direction
       * then extend via general_gather to make double dir[1]-link -> &(temp[icur[0]]),
       * finally gather link[dir[1]] from direction dir[1] 
       * to locally make a another double dir[1]-link later -> gen_pt[8] */
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig] - icur[ig] ) 
			if ( block_stride <=2 ) {
     		tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
       		     											dir[0], EVENANDODD, gen_pt[ig]);
				btag = start_general_gather_field( LINK1 , sizeof(su3_matrix),
																		disp, EVENANDODD, gen_pt[8-ig]);
				wait_gather(tag[ig]);
				wait_general_gather( btag );
				FORALLSITES(i,s) 
				IF_BLOCKED(s, block_stride)	
				IF_ACTIVE(s) 
					mult_su3_nn( (su3_matrix*)(gen_pt[ig][i]), (su3_matrix*)(gen_pt[8-ig][i]), 
											&(dlink_dir1_at_s10[i]) );
											// &(temp[ncur[ig]-icur[ig]][i]) );
				cleanup_general_gather( btag );
				tag[8-ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
																			dir[1], EVENANDODD, gen_pt[8-ig]);
			}

			ig = 1;
  		for ( i = 0; i < 4; i++) 
  			disp[i] = 0;
  		disp[dir[0]] = 1;
  		disp[dir[1]] = 2;

      FORALLSITES(i,s)
      IF_BLOCKED(s, block_stride)
      IF_ACTIVE(s) 
      	su3mat_copy( &(LINK0[i]) , &(temp[icur[ig]][i]) );

      /* request link[dir[0]] from direction dirb[1], 
       * then extend via general_gather to make double dir[0]-link -> &(temp[icur[1]]),
       * finally gather link[dir[0]] from direction dir[0] 
       * to locally make a another double dir[0]-link later -> gen_pt[8-1]*/
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig]-icur[ig] )
			if ( block_stride <=2 ) {
     		tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
											       		    dirb[1], EVENANDODD, gen_pt[ig]);
     		// btag = start_general_gather_field( LINK0 , sizeof(su3_matrix),
				// 														disp, EVENANDODD, gen_pt[8-ig]);
				wait_gather(tag[ig]);
				// wait_general_gather( btag );
				FORALLSITES(i,s) 
				IF_BLOCKED(s, block_stride)	
				IF_ACTIVE(s)
					su3mat_copy( (su3_matrix*)(gen_pt[ig][i]), &(dlink_dir0_at_s02[i]) );
					// mult_su3_nn( (su3_matrix*)(gen_pt[ig][i]), (su3_matrix*)(gen_pt[8-ig][i]), 
					// 						&(dlink_dir0_at_s02[i]) );
					// 						// &(temp[ncur[ig]-icur[ig]][i]) );
				// cleanup_gather( btag );
				// tag[8-ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
				// 														dir[0], EVENANDODD, gen_pt[8-ig]);
			}

			wait_gather(tag[8-0]);

			// extend to double-link in dir[1] direction -> dlink_dir1_at_s00
			FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
  			mult_su3_nn( &(LINK1[i]), (su3_matrix*)(gen_pt[8-0][i]), &(dlink_dir1_at_s00[i]) );
			// cleanup_gather(tag[8-0]); 

      // while waiting for gathers multiply the double- and single-links on the site
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
  			mult_su3_an( &(dlink_dir1_at_s00[i]), &(LINK0[i]), &(su3mat[0][i]) );

      // form a double-single-double-staple for "right" single-link in the plaquette
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
        mult_su3_nn( &(su3mat[0][i]), &(dlink_dir1_at_s10[i]), &(su3mat[3][i]) );

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a double-sinlge-double-staple for "left" single-link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        mult_su3_nn( &(dlink_dir1_at_s00[i]), &(slink_dir0_at_s02[i]), &tempmat );
        mult_su3_na( &tempmat, &(dlink_dir1_at_s10[i]), &(su3mat[5][i]) );
      }

      ig = 2;
     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_ACTIVE(s) 
     		su3mat_copy( &(su3mat[5][i]) , &(temp[icur[ig]][i]) );

      // request double-single-staple for "left" = su3mat[5] from dirb[1]
      for (stride = block_stride; stride >= 1; stride -= 2, icur[ig] = ncur[ig]-icur[ig] )
      if ( block_stride <=2 ) {
      	tag[ig] = start_gather_field( su3mat[5] , sizeof(su3_matrix),
      										  		     dirb[1], EVENANDODD, gen_pt[ig]);
			}

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a single-double-single-staple for "bottom" double-link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        mult_su3_na( &(slink_dir0_at_s02[i]), &(dlink_dir1_at_s10[i]), &(su3mat[1][i]) );
        mult_su3_na( &(su3mat[1][i]), &(LINK0[i]), &(su3mat[4][i]) );
      }

      ig = 3;
     	FORALLSITES(i,s)
     	IF_BLOCKED(s, block_stride)
     	IF_ACTIVE(s) 
     		su3mat_copy( &(su3mat[4][i]) , &(temp[icur[ig]][i]) );

      // request single-double-single-staple for "bottom" = su3mat[4] from dirb[0]
     	tag[ig] = start_gather_field( su3mat[4] , sizeof(su3_matrix),
     										  		     dir[0], EVENANDODD, gen_pt[ig]);

      wait_gather(tag[2]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        // form a single-double-single-staple for "top" double-link in the plaquette
 				mult_su3_an( &(slink_dir0_at_s02[i]), &(su3mat[0][i]), &(su3mat[2][i]) );
 			}

      FORALLSITES(i, s) 
      IF_BLOCKED(s, 2 * block_stride)	
 			IF_LOWER_BULK(s) 
 			{
        // form a staple for "top" link in the plaquette
        // mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
				wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, 2 * block_stride)	
 			IF_LOWER_BULK(s) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
        if ( s->t < UPR_BDRY-1 )
        	wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
				wl1x1t += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        // form a single-double-single-staple for "top" double-link in the plaquette
        mult_su3_an( &(slink_dir0_at_s02[i]), &(su3mat[0][i]), &(su3mat[2][i]) );
			}

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        // form a single-double-single-staple for "top" double-link in the plaquette
        // mult_su3_an( &(slink_dir0_at_s02[i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
        mult_su3_an( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]), &(loop_axb[idx+6][i]) );
      }

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_ACTIVE(s) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
				mult_su3_an( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]), &(loop_axb[idx+3][i]) );
				mult_su3_an( &(su3mat[1][i]), &(su3mat[0][i]), &(loop_axb[idx+0][i]) );
      }

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      cleanup_gather(tag[ig]);
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( &wl1x1t );
  g_doublesum( &wl1x2t );
  g_doublesum( &wl2x1t );

  // get densities
  wl1x1t /= volume;
  wl1x1t *= 2 * block_stride * 2 * block_stride * 2 * block_stride;
  wl1x1t *= nt * 1. / LWR_BULK_LEN;
  wl1x1t /= 3;
  wl1x2t /= volume;
  wl1x2t *= 2 * block_stride * 2 * block_stride * 2 * block_stride;
  wl1x2t *= nt * 1. / ( LWR_BULK_LEN - 1.);
  wl1x2t /= 3;
	wl2x1t /= volume;
  wl2x1t *= 2 * block_stride * 2 * block_stride * 2 * block_stride;
  wl2x1t *= nt * 1. / LWR_BULK_LEN;
  wl2x1t /= 3;

  // node0_printf("pss %g rs2s %g \n",wl1x1s,wl1x2s );
  node0_printf("pss %g pst %g rs2s %g rs2t %g r2st %g rst2 %g\n",
  						wl1x1s,wl1x1t,wl1x2s,wl1x2t,wl2x1t,0.5*(wl1x2t+wl2x1t) );

  // deallocate temporary storage
  destroy_field( &link );

	spatial_blocking();

  link = new_links_from_site( FULLVOL, NLINK );

  // accumulators
  wl1x1s = 0;
  wl1x2s = 0;

  for( dir[0]=YUP; dir[0]<TUP; dir[0]++ ) {
    for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {

    	idx = ( dir[0] == ZUP ? dir[1] : ZUP) ;
			setup_blocked_dirs( dir, dirb, diro );
			// switch( block_stride ) {
			// 	case (1):
			// 		dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
			// 		dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
			// 		break;
			// 	case (2):
			// 		dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
			// 		dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
			// 		break;
			// 	case (4):
			// 		dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
			// 		dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
			// 		break;
			// 	default:
			// 		node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
			// 		terminate(1);
			// 	break;
			// }
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
 			IF_REGION(s, this_region) 
 				mult_su3_an( &(LINK1[i]), &(LINK0[i]), &(su3mat[0][i]) );

      wait_gather(tag[0]);

      // form a staple for "right" link in the plaquette
      FORALLSITES(i, s)
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
        mult_su3_nn( &(su3mat[0][i]), (su3_matrix *)(gen_pt[0][i]), &(su3mat[3][i]) );

      wait_gather(tag[1]);

			//fflush(stdout);printf("dirb[0]=%d dirb[1]=%d\n",dirb[0],dirb[1]);fflush(stdout);
      // form a staple for "left" link in the plaquette
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
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
 			IF_REGION(s, this_region) 
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
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // get the contribution of 1x2 rectangle extended in dir1
        // and of 1x1 plaquette to the accumulators
        wl1x2s += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
        wl1x1s += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }


      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
      }

      maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // form a staple for "top" link in the plaquette
        // mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
        mult_su3_an( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]), &tempmat );
        sub_su3_matrix( &tempmat, &(loop_axb[idx+6][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) 
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) ) 
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
      }
      node0_printf("maxdev %.16g in 1x2 rectangle extended in dir[1] %d\n",maxdev,dir[1]);

      // wait_gather(tag[3]);

			maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
				mult_su3_an( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]), &tempmat );
				sub_su3_matrix( &tempmat, &(loop_axb[idx+3][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) 
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) ) 
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
 			}
      node0_printf("maxdev %.16g in 1x2 rectangle extended in dir[0] %d\n",maxdev,dir[0]);

			maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_REGION(s, this_region) 
 			{
        // and of 1x1 plaquette to the accumulators
				mult_su3_an( &(su3mat[1][i]), &(su3mat[0][i]), &tempmat );
				sub_su3_matrix( &tempmat, &(loop_axb[idx][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) {
        	if (maxdev < cabs( ((complex*)(&tempmat) + ig) )) {
        		node0_printf("i %d: s %d %d %d %d; ig %d\n",i,s->x,s->y,s->z,s->t, ig);
        		dumpmat( &tempmat );
        	}
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) )
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
        }
      }
      node0_printf("maxdev %.16g in 1x1 plaquette \n",maxdev);

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      cleanup_gather(tag[ig]);

	    #undef LINK0
	    #undef LINK1
    } // dir[1]
  } // dir[0]

	// global sum
  g_doublesum( &wl1x1s );
  g_doublesum( &wl1x2s );

  // get densities
  wl1x1s /= volume;
  wl1x1s *= block_stride * block_stride * block_stride;
  wl1x1s *= time_frac;
  wl1x1s /= 3;
  wl1x2s /= volume;
  wl1x2s *= block_stride * block_stride * block_stride;
  wl1x2s *= time_frac;
  wl1x2s /= 6;

  // accumulators
  wl1x1t = 0;
  wl1x2t = 0;
  wl2x1t = 0;

	dir[0]=TUP; {
		for( dir[1]=XUP; dir[1]<dir[0]; dir[1]++ ) {
   		idx = dir[1] + 9;

			setup_blocked_dirs( dir, dirb, diro );
			// dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
			// switch( block_stride ) {
			// 	case (1):
			// 		dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
			// 		break;
			// 	case (2):
			// 		dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
			// 		break;
			// 	case (4):
			// 		dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
			// 		break;
			// 	default:
			// 		node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
			// 		terminate(1);
			// 	break;
			// }
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
				wl2x1t += realtrace_su3( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]) );
      }

      wait_gather(tag[3]);

      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
      {
        // get the contribution of 1x2 rectangle extended in dir[0]
        // and of 1x1 plaquette to the accumulators
				if ( s->t < UPR_BDRY-1 )
          wl1x2t += realtrace_su3( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]) );
				wl1x1t += realtrace_su3( &(su3mat[1][i]), &(su3mat[0][i]) );
      }

      maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
 			{
        // form a staple for "top" link in the plaquette
        // mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(su3mat[0][i]), &(su3mat[2][i]) );
        // get the contribution of 1x2 rectangle extended in dir[1]
        // to the accumulator
        mult_su3_an( (su3_matrix *)(gen_pt[2][i]), &(su3mat[3][i]), &tempmat );
        sub_su3_matrix( &tempmat, &(loop_axb[idx+6][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) {
        	// if (maxdev < cabs( ((complex*)(&tempmat) + ig) )) {
        	// 	node0_printf("i %d: s %d %d %d %d; ig %d\n",i,s->x,s->y,s->z,s->t, ig);
        	// 	dumpmat( &tempmat );
        	// }
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) ) 
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
        }
      }
      node0_printf("maxdev %.16g in 2x1 rectangle extended in dir[1] %d\n",maxdev,dir[1]);

      // wait_gather(tag[3]);

			maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
 			{
        // get the contribution of 1x2 rectangle extended in dir[0]
				mult_su3_an( (su3_matrix *)(gen_pt[3][i]), &(su3mat[2][i]), &tempmat );
				sub_su3_matrix( &tempmat, &(loop_axb[idx+3][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) {
        	// if (maxdev < cabs( ((complex*)(&tempmat) + ig) )) {
        	// 	node0_printf("i %d: s %d %d %d %d; ig %d\n",i,s->x,s->y,s->z,s->t, ig);
        	// 	dumpmat( &tempmat );
        	// }
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) ) 
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
        }
 			}
      node0_printf("maxdev %.16g in 1x2 rectangle extended in dir[0] %d\n",maxdev,dir[0]);

			maxdev = 0.0;
      FORALLSITES(i, s) 
      IF_BLOCKED(s, block_stride)	
 			IF_LOWER_BULK(s) 
 			{
        // and of 1x1 plaquette to the accumulators
				mult_su3_an( &(su3mat[1][i]), &(su3mat[0][i]), &tempmat );
				sub_su3_matrix( &tempmat, &(loop_axb[idx][i]), &tempmat );
        for ( ig = 0; ig < 9; ig ++ ) {
        	maxdev = ( maxdev > cabs( ((complex*)(&tempmat) + ig) )
        					 ? maxdev : cabs( ((complex*)(&tempmat) + ig) ) );
        }
      }
      node0_printf("maxdev %.16g in 1x1 plaquette \n",maxdev);

      // clean up all gathers
      for ( ig = 0; ig < NGATHER; ig++ )
	      cleanup_gather(tag[ig]);

	    #undef LINK0
	    #undef LINK1
	    // break;
    } // dir[1]
  } // dir[0]

  // global sum
  g_doublesum( &wl1x1t );
  g_doublesum( &wl1x2t );
  g_doublesum( &wl2x1t );

  // get densities
  wl1x1t /= volume;
  wl1x1t *= block_stride * block_stride * block_stride;
  wl1x1t *= nt * 1. / LWR_BULK_LEN;
  wl1x1t /= 3;
  wl1x2t /= volume;
  wl1x2t *= block_stride * block_stride * block_stride;
  wl1x2t *= nt * 1. / ( LWR_BULK_LEN - 1.);
  wl1x2t /= 3;
	wl2x1t /= volume;
  wl2x1t *= block_stride * block_stride * block_stride;
  wl2x1t *= nt * 1. / LWR_BULK_LEN;
  wl2x1t /= 3;

  // node0_printf("pss %g rs2s %g \n",wl1x1s,wl1x2s );
  node0_printf("pss %g pst %g rs2s %g rs2t %g r2st %g rst2 %g\n",
  						wl1x1s,wl1x1t,wl1x2s,wl1x2t,wl2x1t,0.5*(wl1x2t+wl2x1t) );

  // deallocate temporary storage
  destroy_field( &link );

  for( ig=0; ig<NLOOP_AXB; ig++ ) 
    free( loop_axb[ig] );
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) 
    free( su3mat[ig] );
	if ( block_stride >= 1 ) 
  	for( ig=0; ig < NBLOCK_STORAGE; ig++ ) 
    	free( temp[ig] );
	#undef NBLOCK_STORAGE
	#undef NGATHER
	#undef NTEMP_STORAGE
	#undef NLINK
  // normal_exit(0);
}

void make_test_lat (void) {
	register int i,dir;
	register site *s;
	FORALLUPDIRBUT(TUP,dir) 
	FORALLSITES(i, s) {
		// if (dir < ZUP)
		set_identity( &(s -> link[dir]) );
		// if (dir == XUP)
		// 	s -> link[dir].e[dir][dir].real = s->x+1;
		// if (dir == YUP)
		// 	s -> link[dir].e[dir][dir].real = s->y+1;
		// if (dir == ZUP)
		// 	s -> link[dir].e[dir][dir].real = s->z+1;
			s -> link[dir].e[XUP][XUP].real = s->x+1;
			s -> link[dir].e[YUP][YUP].real = s->y+1;
			s -> link[dir].e[ZUP][ZUP].real = s->z+1;
	}
}
#endif 
#endif