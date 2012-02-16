/****** fermion_links_fn_qop_milc_P.c  -- ******************/
/* MIMD version 7 */
/* This is a substitute MILC implementation of the QOP fermion links routine */
/* It is used for testing the Level 3 interface */

/* Link fattening routines for varions staggered actions
   CD 9/8/06 separated from quark_stuff.c 
*/

#if ( QOP_Precision == 1)
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix
#else
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix
#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_qopmilc.h"
#include "../include/generic_ks_qop.h"
#include "../include/qop_milc.h"
#define IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#include <quark_action.h>
#include "../include/prefetch.h"
#define FETCH_UP 1

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static void path_product_qop_milc( const int *dir, const int length, 
			    su3_matrix *tempmat1, QOP_GaugeField *gauge);
static void compute_gen_staple_field_qop_milc(su3_matrix *staple, 
					      int mu, int nu, 
					      su3_matrix *linkmu, 
					      su3_matrix *links, 
					      su3_matrix *tfl, Real cf);

/********************************************************************/
/* Sum over paths connecting to nearest neighbor point (fat link) and to third
   nearest neighbor (longlinks) */
/********************************************************************/

/* Taken from fermion_links_fn.c */

static su3_matrix *create_longlinks_qop_milc(QOP_info_t *info, 
					     QOP_asqtad_coeffs_t *coeffs,
					     QOP_GaugeField *gauge) {
  register int i;
  register site *s;
  int ipath,dir;
  int disp[4];
  int num_q_paths = ks_act_paths.num_q_paths;
  Q_path *q_paths = ks_act_paths.q_paths;
  register su3_matrix *long1;
  su3_matrix *staple, *tempmat1;
  int nflop = 1804;
  double dtime;
  su3_matrix *t_ll;
  char myname[] = "create_longlinks_qop_milc";
  
  dtime=-dclock();

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(0);
  }
  /* Allocate space for t_longlink if NULL */
  t_ll = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_ll==NULL){
    printf("NODE %d: no room for t_ll\n",this_node);
    terminate(1);
  }
  
  staple = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(staple == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n",myname);
    terminate(1);
  }

  for (dir=XUP; dir<=TUP; dir++){ /* loop over longlink directions */
    /* set longlink to zero */
    FORALLSITES(i,s){
      long1 = &(t_ll[4*i+dir]);
      clear_su3mat( long1 );
    }

    /* loop over paths, checking for ones with total displacement 3*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=3,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("ipath = %d, found a path:  ",ipath);
for(j=0;j<q_paths[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

	path_product_qop_milc( q_paths[ipath].dir, q_paths[ipath].length, 
			       tempmat1, gauge );
	FORALLSITES(i,s){
	  su3_adjoint( &tempmat1[i], &staple[i] );
	  long1 = &(t_ll[4*i+dir]);
          scalar_mult_add_su3_matrix( long1,
	    &staple[i], -q_paths[ipath].coeff, long1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */

  } /* loop over directions */

  free(staple);
  free(tempmat1);

dtime += dclock();
#ifdef FLTIME
node0_printf("FLTIME(long): time =  %e (Naik) mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

 info->final_sec = dtime;
 info->final_flop = ((double)nflop*volume)/numnodes();
 info->status = QOP_SUCCESS;

 return t_ll;

}  /* load_longlinks_qop_milc() */


/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fatlinks:       61632 for load_fatlinks
 */

/* Taken from fermion_links_fn.c */

/* KS phases and APBC must be in the links. See long comment at 
   end of fermion_forcee_eo_milc.c */
static su3_matrix *create_fatlinks_qop_milc(QOP_info_t *info, 
					  QOP_asqtad_coeffs_t *coeffs,
					  QOP_GaugeField *gauge) {
  register int i;
  register site *s;
  int dir;
  register su3_matrix *fat1;
  su3_matrix *staple, *tempmat1;
  su3_matrix *t_fl;
  char myname[] = "create_fatlinks_qop_milc";

  asqtad_coeffs_t act_path_coeff;

  int  nu,rho,sig ;
  Real one_link; /* needed to fix the problem with the Lepage
		       term */
  /* Convert specific QOP precision to prevailing MILC precision */
#if ( QOP_Precision == 1 )
  su3_matrix *links = create_links_from_qop_milc_F(gauge->g);
#else
  su3_matrix *links = create_links_from_qop_milc_D(gauge->g);
#endif
  Real cf;

  int nflop = 61632;
  double dtime;
  dtime=-dclock();

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(0);
  }

  /* Unload coeff */
  act_path_coeff.one_link =  coeffs->one_link;
  act_path_coeff.naik =  coeffs->naik;
  act_path_coeff.three_staple =  coeffs->three_staple;
  act_path_coeff.five_staple =  coeffs->five_staple;
  act_path_coeff.seven_staple =  coeffs->seven_staple;
  act_path_coeff.lepage =  coeffs->lepage;

  /* Allocate space for t_fl if NULL */
  t_fl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_fl==NULL){
    printf("NODE %d: no room for t_fl\n",this_node);
    terminate(1);
  }
  
  staple = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(staple == NULL){
    printf("%s: Can't malloc temporary\n", myname);
    terminate(1);
  }

  tempmat1 = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("%s: Can't malloc temporary\n", myname);
    terminate(1);
  }

/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
 
 /* to fix up the Lepage term, included by a trick below */
 one_link = (act_path_coeff.one_link - 6.0*act_path_coeff.lepage);

 for (dir=XUP; dir<=TUP; dir++){
   FORALLSITES(i,s) /* Intialize fat links with c_1*U_\mu(x) */
     {
       fat1 = &(t_fl[4*i+dir]);
       //scalar_mult_su3_matrix(&(s->link[dir]), one_link ,
       //		      fat1 );
       scalar_mult_su3_matrix(links + dir*sites_on_node + i, one_link ,
			      fat1 );

     }
   for(nu=XUP; nu<=TUP; nu++) if(nu!=dir)
     {
       cf = act_path_coeff.thre_staple;
       compute_gen_staple_field_qop_milc(staple,dir,nu,links+dir*sites_on_node,
					 links, t_fl, cf);
       /* The Lepage term */
       /* Note this also involves modifying c_1 (above) */
       cf = act_path_coeff.lepage;
       compute_gen_staple_field_qop_milc(NULL,dir,nu,staple,links,t_fl,cf);
       for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   cf = act_path_coeff.five_staple;
	   compute_gen_staple_field_qop_milc( tempmat1, dir, rho, staple, 
					      links, t_fl, cf);
	   for(sig=XUP; sig<=TUP; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 cf = act_path_coeff.seven_staple;
		 compute_gen_staple_field_qop_milc(NULL,dir,sig,
						   tempmat1, links, t_fl, cf);
	       } /* sig */
	 } /* rho */
     } /* nu */
 }/* dir */  

  free(staple);
  free(tempmat1);

#ifdef FLTIME
dtime += dclock();
 node0_printf("FLTIME(Fat): time = %e (Asqtad opt) mflops = %e\n",dtime,
	      (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

 info->final_sec = dtime;
 info->final_flop = ((double)nflop*volume)/numnodes();
 info->status = QOP_SUCCESS;

#if ( QOP_Precision == 1 )
  destroy_links_from_qop_milc_F(links);
#else
  destroy_links_from_qop_milc_D(links);
#endif

 return t_fl;
}  /* load_fatlinks_qop_milc() */


QOP_FermionLinksAsqtad *
  QOP_asqtad_create_L_from_G(QOP_info_t *info,
			     QOP_asqtad_coeffs_t *coeffs,
			     QOP_GaugeField *gauge){
  QOP_info_t flinfo = {0., 0., 0, 0, 0};
  QOP_info_t llinfo = {0., 0., 0, 0, 0};
  QOP_FermionLinksAsqtad *lqop;
  MYSU3_MATRIX **raw_fat_links, **raw_long_links;
  su3_matrix *t_fl, *t_ll;

  t_fl = create_fatlinks_qop_milc( &flinfo, coeffs, gauge);
  raw_fat_links  = create_raw4_G_from_field(t_fl,EVENANDODD);
  free(t_fl);
  if(raw_fat_links == NULL)terminate(1);

  t_ll = create_longlinks_qop_milc( &llinfo, coeffs, gauge);
  raw_long_links = create_raw4_G_from_field(t_ll,EVENANDODD);
  free(t_ll);
  if(raw_long_links == NULL)terminate(1);
  
  lqop = QOP_asqtad_create_L_from_raw((MYREAL **)raw_fat_links, 
				      (MYREAL **)raw_long_links,
				      QOP_EVENODD);

  destroy_raw4_G(raw_fat_links);   raw_fat_links = NULL;
  destroy_raw4_G(raw_long_links);  raw_long_links = NULL;

  info->final_sec = flinfo.final_sec + llinfo.final_sec;
  info->final_flop = flinfo.final_flop + llinfo.final_flop;
  if(flinfo.status == QOP_SUCCESS && llinfo.status == QOP_SUCCESS)
    info->status = QOP_SUCCESS;
  else
    info->status = QOP_FAIL;

  return lqop;
}

#ifndef FN
BOMB THE COMPILE
#endif
static void compute_gen_staple_field_qop_milc(su3_matrix *staple, 
					      int mu, int nu, 
					      su3_matrix *linkmu, 
					      su3_matrix *links,
					      su3_matrix *t_fl,
					      Real cf) {

  su3_matrix tmat1,tmat2;
  su3_matrix *linknu = links + nu*sites_on_node;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat ;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;
  Real coef = cf;

  /* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
    Where the mu link can be any su3_matrix. The result is saved in staple.
    if staple==NULL then the result is not saved.
    It also adds the computed staple to the fatlink[mu] with weight coef.
  */

  /* Upper staple */
  mtag0 = start_gather_field( linkmu, sizeof(su3_matrix), nu, EVENANDODD, gen_pt[0] );
  mtag1 = start_gather_field( linknu, sizeof(su3_matrix), mu, 
			EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( linknu+i, &tmat1, &staple[i] );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_nn( linknu+i, &tmat1, &tmat2 );
      fat1 = &(t_fl[4*i+mu]);
      scalar_mult_add_su3_matrix(fat1, &tmat2, coef,
				 fat1) ;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
  mtag0 = start_gather_field( linknu,
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    mult_su3_an( linknu+i, linkmu+i, &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
				  OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( &staple[i],(su3_matrix *)gen_pt[0][i], 
		      &staple[i] );
      fat1 = &(t_fl[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 &staple[i], coef, 
				 fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(t_fl[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
    }
  }

  free(tempmat);
  cleanup_gather(mtag0);
} /* compute_gen_staple_field_qop_milc */

/* LOOPEND is required now -CD */
#undef FORALLSITES
#define FORALLSITES(i,s) \
{ register int loopend; loopend=sites_on_node; \
for( i=0,  s=lattice ; i<loopend; i++,s++ )
#define END_LOOP }

/* Taken from path_product.c */

void path_product_qop_milc( const int *dir, const int length, 
			    su3_matrix *tempmat1, QOP_GaugeField *gauge) {
    register int i;
    register site *s;
    msg_tag *mtag0 = NULL;
    su3_matrix *tempmat2t, *tempmat3t;
#if ( QOP_Precision == 1 )
    su3_matrix *links = create_links_from_qop_milc_F(gauge->g);
#else
    su3_matrix *links = create_links_from_qop_milc_D(gauge->g);
#endif
    int j;
    /* a forward step leaves the answer in gen_pt[0], which points into
	link, tempmat1 or tempmat2, and backwards step in tempmat1 or tempmat2,
	After a forwards step, need to wait and clean a gather.
	  STEP	leaves answer in
	  even # forward	gen_pt[0]->tempmat1  (gen_pt[0]->link for step 0
	  even # backward	tempmat1
	  odd  # forward	gen_pt[0]->tempmat2
	  odd  # backward	tempmat2
	*/

    /* Trivial path case */
    if(length == 0){
      FORALLSITES(i,s){
	clear_su3mat(&tempmat1[i]);
	tempmat1[i].e[0][0].real = tempmat1[i].e[1][1].real 
	  = tempmat1[i].e[2][2].real = 1.;
      } END_LOOP
      return;
    }

    /* allocate temporary space */
    tempmat3t = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
    tempmat2t = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );

    /* j=0 */
    if( GOES_FORWARDS(dir[0]) )  {
      //mtag0 = start_gather_site( F_OFFSET(link[dir[0]]), sizeof(su3_matrix),
      mtag0 = start_gather_field( 
            links + dir[0]*sites_on_node, sizeof(su3_matrix),
	    OPP_DIR(dir[0]), EVENANDODD, gen_pt[0] );
    }
    else{  /* if GOES_BACKWARDS(dir[0]) */
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &tempmat1[i+FETCHUP] );
	  }
	  //su3_adjoint(&(s->link[OPP_DIR(dir[0])]),&tempmat1[i] );
	  su3_adjoint(links + OPP_DIR(dir[0])*sites_on_node + i ,&tempmat1[i] );
	} END_LOOP
    }

    for(j=1;j<length;j++) {
	if( j%2==1 ){
	    if( GOES_FORWARDS(dir[j]) ) {
	      if( GOES_FORWARDS(dir[j-1]) ){
	        wait_gather(mtag0);
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		    prefetch_M(  &(tempmat2t[i+FETCH_UP]) );
		  }
		  //mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir[j]]),
		  //&(tempmat2t[i]) );
		  mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), 
			       links + dir[j]*sites_on_node + i,
			       &(tempmat2t[i]) );
	        } END_LOOP
	        cleanup_gather(mtag0);
	      }
	      else{ /* last link was backwards */
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( &(tempmat1[i+FETCH_UP]) );
		    prefetch_M(  &(tempmat2t[i+FETCH_UP]) );
		  }
		  //mult_su3_nn( &tempmat1[i],&(s->link[dir[j]]),
		  //&(tempmat2t[i]) );
		  mult_su3_nn( &tempmat1[i],
			       links + dir[j]*sites_on_node + i,
			       &(tempmat2t[i]) );
	        } END_LOOP
	      }
	      mtag0 = start_gather_field( tempmat2t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }  /* for GOES_FORWARDS */

	    else{ /* GOES_BACKWARDS(dir[j]), which is an odd numbered step */
	      if( GOES_FORWARDS(dir[j-1]) ){
	        wait_gather(mtag0);
	        FORALLSITES(i,s){
		  if( i < loopend-FETCH_UP ){
		    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  }
	          su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(tempmat3t[i]) );
	        } END_LOOP
	        cleanup_gather(mtag0);
	        mtag0 = start_gather_field( tempmat3t, sizeof(su3_matrix),
		  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	      }
	      else{ /*last step was backwards */
	        mtag0 = start_gather_field( tempmat1, sizeof(su3_matrix),
		  OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	      }
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  //prefetch_M(  &((s+FETCH_UP)->link[OPP_DIR(dir[j])]) );
		  prefetch_M(  links + OPP_DIR(dir[j])*sites_on_node + i + FETCH_UP );
		}
		// mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		// &(s->link[OPP_DIR(dir[j])]), &(tempmat2t[i]) );
		  mult_su3_na((su3_matrix *)(gen_pt[0][i]),
			      links + OPP_DIR(dir[j])*sites_on_node + i, 
			      &(tempmat2t[i]) );
	      } END_LOOP
	      cleanup_gather(mtag0);
	    } /* end for GOES_BACKWARDS */
	} /* end for j=odd */

	else{	/* j=even */
	  if( GOES_FORWARDS(dir[j]) ) {
	    if( GOES_FORWARDS(dir[j-1]) ){
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  //prefetch_M(  &((s+FETCH_UP)->link[dir[j]]) );
		  prefetch_M( links + dir[j]*sites_on_node + i + FETCH_UP );
		}
		//mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir[j]]),
		//  &tempmat1[i] );
		mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), 
			     links + dir[j]*sites_on_node + i,
			     &tempmat1[i] );
	      } END_LOOP
	      cleanup_gather(mtag0);
	    }
	    else{ /* last link goes backwards */
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( &(tempmat2t[i+FETCH_UP]) );
		  //prefetch_M(  &((s+FETCH_UP)->link[dir[j]]) );
		  prefetch_M(links + dir[j]*sites_on_node + i + FETCH_UP);
		}
		//mult_su3_nn( &(tempmat2t[i]),&(s->link[dir[j]]),
		//  &tempmat1[i] );
		mult_su3_nn( &(tempmat2t[i]),
			     links + dir[j]*sites_on_node + i,
			     &tempmat1[i] );
	      } END_LOOP
	    }
	    mtag0 = start_gather_field( tempmat1, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	  }  /* for GOES_FORWARDS */

	  else{ /* GOES_BACKWARDS(dir[j]) */
	    if( GOES_FORWARDS(dir[j-1]) ){
	      wait_gather(mtag0);
	      FORALLSITES(i,s){
		if( i < loopend-FETCH_UP ){
		  prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		  prefetch_M(  &(tempmat3t[i+FETCH_UP]) );
		}
	        su3mat_copy((su3_matrix *)(gen_pt[0][i]),&(tempmat3t[i]) ); 
	      } END_LOOP
	      cleanup_gather(mtag0);
	      mtag0 = start_gather_field( tempmat3t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }
	    else{ /* last step was backwards */
	      mtag0 = start_gather_field( tempmat2t, sizeof(su3_matrix),
		OPP_DIR(dir[j]), EVENANDODD, gen_pt[0] );
	    }
	    wait_gather(mtag0);
	    FORALLSITES(i,s){
	      if( i < loopend-FETCH_UP ){
		prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
		//prefetch_M(  &((s+FETCH_UP)->link[OPP_DIR(dir[j])]) );
		prefetch_M( links + OPP_DIR(dir[j])*sites_on_node + i + FETCH_UP );
	      }
	      //mult_su3_na((su3_matrix *)(gen_pt[0][i]),
	      //    &(s->link[OPP_DIR(dir[j])]), &tempmat1[i] );
	      mult_su3_na((su3_matrix *)(gen_pt[0][i]),
		    links + OPP_DIR(dir[j])*sites_on_node + i, &tempmat1[i] );
	    } END_LOOP
	    cleanup_gather(mtag0);
	  } /* for GOES_BACKWARDS */
	} /* for j=even */

    }  /* j=link in loop */

    /* Want to end in tempmat1 */
    if( length%2==0 ){  /* last step was odd */
      if( GOES_FORWARDS(dir[length-1]) ){
	wait_gather(mtag0);
	  FORALLSITES(i,s){
	    if( i < loopend-FETCH_UP ){
	      prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
	    }
	    su3mat_copy((su3_matrix *)(gen_pt[0][i]),&tempmat1[i] ); 
	} END_LOOP
	cleanup_gather(mtag0);
      }
      else{
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &(tempmat2t[i+FETCH_UP]) );
	  }
	  su3mat_copy(&(tempmat2t[i]),&tempmat1[i] );
	} END_LOOP
      }
    }
    else{ /* odd length path */
      if( GOES_FORWARDS(dir[length-1]) ){
	wait_gather(mtag0);
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( (su3_matrix *)(gen_pt[0][i+FETCH_UP]) );
	  }
	  su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(tempmat3t[i]) );
	} END_LOOP
	cleanup_gather(mtag0);
	FORALLSITES(i,s){
	  if( i < loopend-FETCH_UP ){
	    prefetch_M( &(tempmat3t[i+FETCH_UP]) );
	  }
	  su3mat_copy( &(tempmat3t[i]), &tempmat1[i] );
	} END_LOOP
      }
      else{
      }
    }
    free(tempmat3t);
    free(tempmat2t);

#if ( QOP_Precision == 1 )
    destroy_links_from_qop_milc_F(links);
#else
    destroy_links_from_qop_milc_D(links);
#endif

} /* path_product_qop_milc */


