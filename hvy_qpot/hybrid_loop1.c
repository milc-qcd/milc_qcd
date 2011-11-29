/********************** hybrid_loop1.c *******************************/
/* MIMD version 7 */
/* Computes on-axis Wilson loops appropriate for excited flux-tube
   potentials (including those relevant to hybrids) with symmetry
   \Pi_u, and Delta_g 

   ALSO computes the conventional Wilson loop with symmetry \Sigma_g^+

   For the diatomic molecular notation applied to this problem, see
   K.J. Juge, J. Kuti, and C.J. Morningstar, Lattice '97 Nucl Phys
   (Proc Supp) B 63, 326 (1998)

   In terms of representations of D_4h these are
   \Sigma_g^+ = A_1g
   \Pi_u = E_u
   \Delta_g = B_1g

 */

/* We generate the hybrid loops by doing parallel displacements of the
   linear space-like Wilson line segment.  The displacement is in the
   plane perpendicular to the space-like segment as follows (Q denotes
   the static quark and o the locus of the Wilson line):

   \Sigma_g^+    Qo    (no displacement - conventional Wilson loop)

                                          o
                                          |
                                          |
   \Pi_u         Q----o - o----Q    or    Q - Q
                                              |
                                              |
                                              o

                                    o
                                    |
                                    |
   \Delta_g      Q----o + o----Q -  Q - Q
                                        |
                                        |
                                        o

*/


/* Assumes temporal axial gauge with the nontrivial link copied to
   all time links as in generic/ax_gauge.c */

/* based on w_loop1 by UMH */

/* 4/16/02 CD */

#define MAX_PATH_LENGTH 4

typedef struct
{
  char name[12];                  /* name of path */
  int length;                     /* length of path */
  int disp[4];                    /* net displacement of path */
  int dir[MAX_PATH_LENGTH];       /* list of directions in path */
} link_path;

#include "hvy_qpot_includes.h"
#include "../include/dirs.h"


#define WILS_LOOP1(i,t,r)   wils_loop1[ (i)*nxh*nth + (t)*nxh + r ]
#define TRANS_LINKS(i,j)    (trans_links + (i)*NTRANS + (j))
#define TRANS_LINKS_F(i,j)  (trans_links_f + (i)*NTRANS_F + (j))
/* flux_links must be laid out with site index varying most rapidly since
   we want to gather j's separately */
#define FLUX_LINKS(i,j)     (flux_links + (i) + (j)*sites_on_node)
#define FLUX_LINKS_F(i,j)   (flux_links_f + (i)*NFLUX_F + (j))

void hybrid_loop1(int tot_smear) {

  char myname[] = "hybrid_loop1";
  register int i,j,dir,r,t;
  int dir1=0,dir2=0,trans_path1=0,trans_path2=0;
  int disp[4];
  int nth,nxh;
  register site *s;
  su3_matrix tmat1,tmat2;
  su3_matrix *tmatp;
  Real *wils_loop1;
  su3_matrix *flux_links_f;
  su3_matrix *trans_links, *trans_links_f; 
  su3_matrix *flux_links;
  su3_matrix *tempmat1;
  
  /* Names of messages used to index "mtag" and "gen_pt" */
  /* NMSGS must not exceed N_POINTERS */
  enum{ M_S_LINK, M_F_LINKS_F, M_T_LINKS_F, M_STAP_POS1, 
	  M_STAP_NEG1, M_STAP_POS2, M_STAP_NEG2, NMSGS };
  
  msg_tag *mtag[NMSGS];

  /* Names of flux tube shapes built from transverse links 
     and on-axis links used to index "flux_links" */
  
  enum{ S_LINK, STAP_POS1, STAP_NEG1, STAP_POS2, STAP_NEG2, NFLUX };

  /* Names of transverse links 
     used to index "trans_path", "trans_links" and "trans_links_f" */
  enum{ XX, YY, ZZ, NTRANS };

  /* Names of space-shifted transverse links */
  /* used to index "trans_links_f" */
  enum{ TRANS_PATH1_F, TRANS_PATH2_F, T_LINK_F, NTRANS_F };
  
  /* Paths for transverse links */
  
  const link_path trans_path[NTRANS] =
    {
      { "XX",   2, {2,0,0,0}, {XUP, XUP, NODIR, NODIR} },
      { "YY",   2, {0,2,0,0}, {YUP, YUP, NODIR, NODIR} },
      { "ZZ",   2, {0,0,2,0}, {ZUP, ZUP, NODIR, NODIR} }
    };

  /* Names of flux tube shapes transported forward in time */
  /* Used to index "flux_links_f" */
  enum{ S_LINK_F, STAP_POS1_F, NFLUX_F };

  /* Names of loop observables */
  /* Used to index "wils_loop1" */
  enum{ W_LOOP1, STAP_SIG_GP1, STAP_PI_U1, STAP_DELTA_G1, NWLOOP1 };

  /* Check the rules */
  if(NMSGS > N_POINTERS){
    if(this_node == 0)fprintf(stderr,"%s: Aborted. gen_pt array too small.",myname);
    terminate(1);
  }

  if( nx != ny || nx != ny){
    if(this_node == 0)fprintf(stderr,"%s: Aborted. requires nx=ny=nz",myname);
    terminate(1);
  }
  
  /* Allocate space for observables */
  nth = nt/2;  nxh = nx/2;
  wils_loop1 = (Real *)malloc(nth*nxh*sizeof(Real)*NWLOOP1);
  if(wils_loop1 == NULL){
    fprintf(stderr,"%s: CAN'T MALLOC wils_loop1\n",myname);
    fflush(stderr);
    terminate(1);
  }
  
  for(i=0;i<NWLOOP1;i++) for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    WILS_LOOP1(i,t,r) = 0.0;

  /* Allocate space for timeward shifted flux tube shapes */
  flux_links_f = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix)*NFLUX_F);
  if(flux_links_f == NULL){
    fprintf(stderr,"%s: CAN'T MALLOC flux_links_f\n",myname);
    fflush(stderr);
    terminate(1);
  }
  
  /* Allocate space for transverse link products */
  trans_links = 
    (su3_matrix *)malloc(NTRANS*sites_on_node*sizeof(su3_matrix));
  if(trans_links == NULL){
    fprintf(stderr,"%s: CAN'T MALLOC trans_links\n",myname);
    fflush(stderr);
    terminate(1);
  }
  
  /* Allocate space for shifted auxiliary link products */
  trans_links_f = 
    (su3_matrix *)malloc(NTRANS_F*sites_on_node*sizeof(su3_matrix));
  if(trans_links_f == NULL){
    fprintf(stderr,"%s: CAN'T MALLOC trans_links_f\n",myname);
    fflush(stderr);
    terminate(1);
  }
     
  tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s(%d): Can't malloc temporary\n",myname,this_node);
    terminate(1);
  }

  /* Compute and store products of transverse links */
  /* trans_links[i][j] is set to the link product for the
     jth path type that ends at site i */
  for(j = 0; j < NTRANS; j++){
    path_product(trans_path[j].dir, trans_path[j].length, tempmat1);
    FORALLSITES(i,s){
      su3mat_copy(tempmat1+i, TRANS_LINKS(i,j));
    }
  }

  free(tempmat1);

  /* Allocate space for flux-tube shapes */
  flux_links = 
    (su3_matrix *)malloc(NFLUX*sites_on_node*sizeof(su3_matrix));
  if(flux_links == NULL){
    fprintf(stderr,"%s: CAN'T MALLOC flux_links\n",myname);
    fflush(stderr);
    terminate(1);
  }
  
  
  /* Main loop over spatial directions */
  for(dir=XUP;dir<=ZUP;dir++){
    
      /* Build giant staples */
      switch(dir){
      case XUP:
	dir1 = YUP; dir2 = ZUP; 
	trans_path1 = YY; trans_path2 = ZZ;
	break;
      case YUP:
	dir1 = ZUP; dir2 = XUP; 
	trans_path1 = ZZ; trans_path2 = XX;
	break;
      case ZUP:
	dir1 = XUP; dir2 = YUP; 
	trans_path1 = XX; trans_path2 = YY;
	break;
      default:
	if(this_node == 0)fprintf(stderr,"%s unknown direction %d\n", myname, dir);
	break;
      }
      
      /* Initialization for loop operators */
      
      /* Start building s_link */
      /* Prepare to shift trans_links */
      FORALLSITES(i,s){
	su3mat_copy( &(s->link[dir]), FLUX_LINKS(i,S_LINK) );
	su3mat_copy(TRANS_LINKS(i,trans_path1), 
		    TRANS_LINKS_F(i,TRANS_PATH1_F));
	su3mat_copy(TRANS_LINKS(i,trans_path2), 
		    TRANS_LINKS_F(i,TRANS_PATH2_F));
	su3mat_copy( &(s->link[TUP]), TRANS_LINKS_F(i,T_LINK_F));
      }
      
      /* Start spatial (dir) transport of all transverse links */
      /* gen_pt[M_T_LINKS_F][i] will point to trans_links_f[i+dir] 
	 for site i */
      mtag[M_T_LINKS_F] = start_gather_field( 
		       (void *)(trans_links_f), 
		       NTRANS_F*sizeof(su3_matrix),
		       dir, EVENANDODD, gen_pt[M_T_LINKS_F] );
      
      /* Recursively construct the space-like segments and compute
	 the Wilson loops with that segment */
      
      for(r=0;r<nxh;r++){
	
	if( r>0 ){
	  /* Collect the space-like segment and extend it by one
	     link in the direction dir. */
	  /* s_link[i] then contains the product of links BEGINNING at
	     site i */
	  wait_gather( mtag[M_S_LINK]);
	  FORALLSITES(i,s){
	    su3mat_copy( (su3_matrix *)(gen_pt[M_S_LINK][i]), &(s->staple));
	  }
	  FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir]), &(s->staple), 
			 FLUX_LINKS(i,S_LINK) );
	  }
	}
	
	/* Shift the space-like segment for next r, if still needed. */
	/* gen_pt[M_S_LINK][i] will point to s_link[i+dir] */
	if( r==0 ){
	  mtag[M_S_LINK] = start_gather_field( 
			   (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix),
			   dir, EVENANDODD, gen_pt[M_S_LINK] );
	}
	else if( r<(nxh-1) ){
	  restart_gather_field( 
			   (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix),
			   dir, EVENANDODD, gen_pt[M_S_LINK], mtag[M_S_LINK] );
	}
	else{
	  cleanup_gather( mtag[M_S_LINK]);
	}
	
	/* Collect the transverse links shifted from dir. */
	/* trans_links_f[i][j] <- trans_links_f[i+dir][j] */
	wait_gather( mtag[M_T_LINKS_F]);
	for(j = 0; j < NTRANS_F; j++){ 
	  FORALLSITES(i,s){
	    su3mat_copy( (su3_matrix *)(gen_pt[M_T_LINKS_F][i]) + j, 
			 &(s->staple));
	  }
	  FORALLSITES(i,s){
	    su3mat_copy( &(s->staple), TRANS_LINKS_F(i,j) );
	  }
	}
	
	/* Start making flux_links[i][STAP_NEG1] <- backward s_link
	   at -disp */
	FORALLUPDIR(j){disp[j] = -trans_path[dir1].disp[j];}
	mtag[M_STAP_NEG1] = start_general_gather_field(
			 (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix),
			 disp, EVENANDODD, gen_pt[M_STAP_NEG1] );

	/* flux_links[i][STAP_POS1] <- forward staple at disp */
	FORALLSITES(i,s){
	  mult_su3_nn( TRANS_LINKS(i,trans_path1),
		       FLUX_LINKS(i,S_LINK), &tmat1 );
	  mult_su3_na( &tmat1, 
		       TRANS_LINKS_F(i,TRANS_PATH1_F), 
		       FLUX_LINKS(i,STAP_POS1) );
	}
	
	/* Can't overlap gen_gathers!! */
	wait_general_gather(mtag[M_STAP_NEG1]);
	
	
	/* Start making flux_links[i][STAP_NEG2] <- backward s_link
	   at -disp */
	FORALLUPDIR(j){disp[j] = -trans_path[dir2].disp[j];}
	mtag[M_STAP_NEG2] = start_general_gather_field(
			 (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix),
			 disp, EVENANDODD, gen_pt[M_STAP_NEG2] );
	
	
	/* flux_links[i][STAP_POS2] <- forward staple at disp */
	FORALLSITES(i,s){
	  mult_su3_nn( TRANS_LINKS(i,trans_path2), 
		       FLUX_LINKS(i,S_LINK), &tmat1 );
	  mult_su3_na( &tmat1, 
		       TRANS_LINKS_F(i,TRANS_PATH2_F), 
		       FLUX_LINKS(i,STAP_POS2) );
	}
	
	wait_general_gather(mtag[M_STAP_NEG2]);

	/* Shift staples to the sites where the t_links join them */
	FORALLUPDIR(j){disp[j] = trans_path[dir1].disp[j];}
	mtag[M_STAP_POS1] = start_general_gather_field(
			 (void *)FLUX_LINKS(0,STAP_POS1), sizeof(su3_matrix),
			 disp, EVENANDODD, gen_pt[M_STAP_POS1] );
	

	/* Finish flux_links[i][STAP_NEG1] <- backward staple in dir1 */
	FORALLSITES(i,s){
	  mult_su3_an( TRANS_LINKS(i,trans_path1), 
		       (su3_matrix *)gen_pt[M_STAP_NEG1][i], 
		       &tmat1 );
	  mult_su3_nn( &tmat1, 
		       TRANS_LINKS_F(i,TRANS_PATH1_F), 
		       FLUX_LINKS(i,STAP_NEG1) );
	}
	
	wait_general_gather(mtag[M_STAP_POS1]);


	/* Collect results of the shift of forward staples */
	FORALLSITES(i,s){
	  su3mat_copy( (su3_matrix *)gen_pt[M_STAP_POS1][i], &(s->staple));
	}
	FORALLSITES(i,s){
	  su3mat_copy( &(s->staple), FLUX_LINKS(i,STAP_POS1) );
	}

	/* Shift staples to the sites where the t_links join them */
	FORALLUPDIR(j){disp[j] = trans_path[dir2].disp[j];}
	mtag[M_STAP_POS2] = start_general_gather_field(
			 (void *)FLUX_LINKS(0,STAP_POS2), sizeof(su3_matrix),
			 disp, EVENANDODD, gen_pt[M_STAP_POS2] );
	
	/* Finish flux_links[i][STAP_NEG2] <- forward staple in dir2 */
	FORALLSITES(i,s){
	  mult_su3_an( TRANS_LINKS(i,trans_path2), 
		       (su3_matrix *)gen_pt[M_STAP_NEG2][i], 
		       &tmat1 );
	  mult_su3_nn( &tmat1, 
		       TRANS_LINKS_F(i,TRANS_PATH2_F), 
		       FLUX_LINKS(i,STAP_NEG2) );
	}
	
	wait_general_gather(mtag[M_STAP_POS2]);

	/* Collect results of the shift of forward staples */
	FORALLSITES(i,s){
	  su3mat_copy( (su3_matrix *)gen_pt[M_STAP_POS2][i], &(s->staple));
	}
	FORALLSITES(i,s){
	  su3mat_copy( &(s->staple), FLUX_LINKS(i,STAP_POS2) );
	}

	/* Prepare the forward space-like segment for parallel
	   transport in time */
	/* Note: we time shift only one giant staple */
	/* s_link_f[i] <- s_link[i] */
	/* stap_pos1_f[i] <- stap_pos[i] */

	FORALLSITES(i,s){
	  su3mat_copy( FLUX_LINKS(i,S_LINK), FLUX_LINKS_F(i,S_LINK_F) );
	  su3mat_copy( FLUX_LINKS(i,STAP_POS1), FLUX_LINKS_F(i,STAP_POS1_F) );
	}
	
	cleanup_general_gather(mtag[M_STAP_NEG1]);
	cleanup_general_gather(mtag[M_STAP_POS1]);
	cleanup_general_gather(mtag[M_STAP_NEG2]);
	cleanup_general_gather(mtag[M_STAP_POS2]);

	/* Start gather of forward space-like segments for next t 
	   gen_pt[M_F_LINKS_F][i] will point to flux_links_f[i+TUP] */
	mtag[M_F_LINKS_F] = start_gather_field( 
			   (void *)flux_links_f, 
			   sizeof(su3_matrix)*NFLUX_F,
			   TUP, EVENANDODD, gen_pt[M_F_LINKS_F] );
	
	
	/* Recursively compute the Wilson loops of different time
	   extent for fixed spatial extent */
	for(t=0;t<nth;t++){
	  
	  /* Collect forward space-like segments */
	  wait_gather( mtag[M_F_LINKS_F]);
	  FORALLSITES(i,s){
	    su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+S_LINK_F, 
			 &(s->staple));
	    su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+STAP_POS1_F, 
			 &(s->diag));
	  }
	  FORALLSITES(i,s){
	    su3mat_copy( &(s->staple), FLUX_LINKS_F(i,S_LINK_F) );
	    su3mat_copy( &(s->diag), FLUX_LINKS_F(i,STAP_POS1_F) );
	  }
	  
	  /* Start gather for next t, if still needed. */
	  if( t<(nth-1) ){
	    restart_gather_field( 
		     (void *)flux_links_f, 
		     sizeof(su3_matrix)*NFLUX_F,
		     TUP, EVENANDODD, gen_pt[M_F_LINKS_F], 
		     mtag[M_F_LINKS_F] );
	  }
	  else{
	    cleanup_gather( mtag[M_F_LINKS_F]);
	  }
	  
	  /* Finally, compute the Wilson loops. */
	  
	  /* Compute naive Wilson loop term */
	  FORALLSITES(i,s){
	    /* If the loop extends past t = Nt - 1 the temporal
	       axial gauge link is nontrivial */
	    if( ((s->t)+t+1)>=nt ){
	      mult_su3_nn( &(s->link[TUP]), 
			   FLUX_LINKS_F(i,S_LINK_F), &tmat1);
	      mult_su3_na( &tmat1, 
			   TRANS_LINKS_F(i,T_LINK_F), &tmat2);
	      tmatp = &tmat2;
	    }
	    else
	      tmatp = FLUX_LINKS_F(i,S_LINK_F);
	    
	    WILS_LOOP1(W_LOOP1,t,r) += 
	      realtrace_su3( tmatp, FLUX_LINKS(i,S_LINK) );
	    
	  }
	  
	  /* Compute big staple terms */
	  
	  /* Do projection for source and use dir1 for sink */
	  
	  FORALLSITES(i,s){
	    /* If the loop extends past t = Nt - 1 the temporal
	       axial gauge link is nontrivial */
	    if( ((s->t)+t+1)>=nt ){
	      mult_su3_nn( &(s->link[TUP]), 
			   FLUX_LINKS_F(i,STAP_POS1_F), &tmat1);
	      mult_su3_na( &tmat1, 
			   TRANS_LINKS_F(i,T_LINK_F), &tmat2);
	      tmatp = &tmat2;  /* Use tmat2 to close loop */
	    }
	    else
	      /* Use stap_pos1_f to close loop */
	      tmatp = FLUX_LINKS_F(i,STAP_POS1_F);
	    
	    /* Projection of source onto irreps */
	    
	    WILS_LOOP1(STAP_SIG_GP1,t,r) +=
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) +
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) +
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) );

	    WILS_LOOP1(STAP_PI_U1,t,r) +=
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) -
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) );


	    WILS_LOOP1(STAP_DELTA_G1,t,r) +=
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) -
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) -
	      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) );
	  }
	  
	} /* end loop over t */
	
	/* Start gather of forward time-like links for next r. */
	if( r<(nxh-1) ){
	  restart_gather_field( 
		   (void *)trans_links_f, 
		   NTRANS_F*sizeof(su3_matrix),
		   dir, EVENANDODD, gen_pt[M_T_LINKS_F], 
		   mtag[M_T_LINKS_F] );
	}
	else{
	  cleanup_gather( mtag[M_T_LINKS_F]);
	}
	
      } /* end loop over r */
      
  } /* end loop over dir */
  
  /* Normalize and print the loops */
  g_vecfloatsum(wils_loop1,nxh*nth*NWLOOP1);
  
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    node0_printf("WILS_LOOP1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)WILS_LOOP1(W_LOOP1,t,r) / (double)(9*volume) );
  
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    node0_printf("STAP_SIG_GP1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)WILS_LOOP1(STAP_SIG_GP1,t,r) / (double)(9*volume) );
    
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    node0_printf("STAP_PI_U1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)WILS_LOOP1(STAP_PI_U1,t,r) / (double)(9*volume) );
    
  for(t=0;t<nth;t++) for(r=0;r<nxh;r++)
    node0_printf("STAP_DELTA_G1_%d  %d  %d  %e\n", tot_smear, r, t, 
		 (double)WILS_LOOP1(STAP_DELTA_G1,t,r) / (double)(9*volume) );
  
  free( trans_links ); 
  free( trans_links_f ); 
  free( wils_loop1 );
  free( flux_links_f );
  free( flux_links );
  
} /* hybrid_loop1 */

