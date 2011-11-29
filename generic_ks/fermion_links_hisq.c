OBSOLETE

/****** fermion_links_hisq.c  -- ******************/
/* MIMD version 7 */
/* Link fattening routines for various hisq actions
*/

/* External entry points

   init_ferm_links
   load_ferm_links
   load_ferm_links_dmdu0 (ifdef DM_DU0)
   invalidate_all_ferm_links

 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/umethod.h"
#ifdef HAVE_QOP
#include "../include/generic_ks_qop.h"
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

// high precison matrix dump
void dumpmat_hp( su3_matrix *m ){
int i,j;
    for(i=0;i<3;i++){
	for(j=0;j<3;j++)printf("(%.15f,%.15f)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}

// Some notation:
// U -- original link, SU(3), copied to "field" from "site"
// V -- after 1st level of smearing, non-SU(3)
// Y -- unitarized, U(3)
// W -- special unitarized, SU(3)
// X -- after 2nd level of smearing, non-SU(3)

// Forward declarations.
static void  load_U_from_site(ferm_links_t *fn);
static void  load_V_from_U(ferm_links_t *fn, ks_component_paths *ap1);
static void  load_Y_from_V(ferm_links_t *fn, int umethod);
void  load_W_from_Y(ferm_links_t *fn, int umethod, int ugroup);
void  load_X_from_W(ferm_links_t *fn, ks_component_paths *ap2,
                    su3_matrix **X_fatlinks, su3_matrix **X_longlink);
static void  reorder_links(su3_matrix *link4, su3_matrix *link[]);


/* Make the various fat links.
	U = copy of links in site structure
	V = once smeared
	W = projected onto U3
	Y = projected onto SU3
	Xfat = twice smeared (depends on mass used in Naik term)
	Xlong = three link product (depends on mass used in Naik term)
        X_fatbacklink = adjoint of links coming in from backwards, for DBLSTORE
        X_longbacklink = adjoint of links coming in from backwards, for DBLSTORE
*/
void 
load_ferm_links(ferm_links_t *fn){
  ks_action_paths *ap = fn->ap;
  int i,dir,inaik;
  su3_matrix **matfield;
  hisq_links_t *hl = &fn->hl;
  site *s;
//   su3_matrix *Xt_fatbacklink  = fn->fatback;
//   su3_matrix *Xt_longbacklink = fn->lngback;

  int retwist = fn->twist_in;  /* So we remember to restore the boundary twist if needed. */
  char myname[] = "load_ferm_links";

  if(fn->ap == NULL){
    printf("load_ferm_links(%d): KS action paths not initialized\n",this_node);
    terminate(1);
  }

  if( !(hl->valid_all_links == 1) ) {
    // Make sure all the required space is allocated
    int ii, jj, ishift=4, imax=ishift+2*fn->hl.n_naiks;
    if( phases_in != 1){
      node0_printf("BOTCH: %s needs phases in\n",myname); terminate(0);
    }
    for( i=0; i<imax; i++){
      if(i<4) {
        switch(i){
        case 0: matfield = hl->U_link; break;
        case 1: matfield = hl->V_link; break;
        case 2: matfield = hl->Y_unitlink; break;
        // 3 has to be handled in specail way, see below
        case 3: matfield = hl->W_unitlink; break;
	default: matfield = NULL;
        }
      }
      if(i>3) {
        jj = i - ishift;
        ii = jj/2;
        if(0==jj%2) {
          matfield = hl->XX_fatlink[ii];
        }
        else {
          matfield = hl->XX_longlink[ii];
        }
      }
      if(i>=imax) {
        node0_printf("%s BOTCH\n",myname); terminate(0);
      }
      for(dir=XUP;dir<=TUP;dir++){
        // In U(3) case hl->W_unitlink array is not malloc'ed
        // but set as a pointer to hl->Y_unitlink array
        if( !(ap->ugroup==UNITARIZE_U3 && 3==i ) ) { // 3 corresponds to W_unitlink
          if(matfield[dir] == NULL){
            matfield[dir] = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));
            if(matfield[dir]==NULL){
              printf("%s(%d): no room for matfield %d\n", myname,this_node,i); 
              terminate(1);
            }
          }
        }
        else {
          hl->W_unitlink[dir] = hl->Y_unitlink[dir];
        }
      } // dir loop
    } // i loop

    if(fn->fl.fat == NULL){
      fn->fl.fat = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
      if(fn->fl.fat==NULL){
        printf("%s(%d): no room for fn->fl.fat\n", myname,this_node); 
        terminate(1);
      }
    }

    if(fn->fl.lng == NULL){
      fn->fl.lng = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
      if(fn->fl.lng==NULL){
        printf("%s(%d): no room for fn->fl.lng\n", myname,this_node); 
        terminate(1);
      }
    }

    load_U_from_site(fn);
    load_V_from_U(fn, &ap->p1);
    load_Y_from_V(fn,ap->umethod);
    load_W_from_Y(fn,ap->umethod,ap->ugroup);

    // building different sets of X links SKETCH:
    // if n_naiks > 1, say, n_naiks = 3 in this example
    // a) calculate 3rd path table set in XX_fat/long[0]
    // b) multiply by eps_naik[i] and store in XX_fat/long[i], i.e.
    //    XX_fat/long[1] = eps_naik[1]*XX_fat/long[0],
    //    XX_fat/long[2] = eps_naik[2]*XX_fat/long[0]
    // c) calculate 2nd path table set in XX_fat/long[0],
    //    this is the set with 0 correction
    // d) add XX_fat/long[0] to all other sets, i.e.
    //    XX_fat/long[1] += XX_fat/long[0],
    //    XX_fat/long[2] += XX_fat/long[0]
    if(fn->hl.n_naiks>1 ) {
      // 3rd path table set
      hl->valid_X_links = 0;
      load_X_from_W(fn, &ap->p3, hl->XX_fatlink[0], hl->XX_longlink[0]);
      for( inaik=1;inaik<fn->hl.n_naiks;inaik++ ) {
        for(dir=XUP;dir<=TUP;dir++) {
          FORALLSITES(i,s) {
          	// multiply by epsilon correction
            scalar_mult_su3_matrix( &( hl->XX_fatlink[0][dir][i] ),
              fn->hl.eps_naik[inaik], &( hl->XX_fatlink[inaik][dir][i] ) );
            scalar_mult_su3_matrix( &( hl->XX_longlink[0][dir][i] ),
              fn->hl.eps_naik[inaik], &( hl->XX_longlink[inaik][dir][i] ) );
          }
        }
      }
      // 2nd path table set
      hl->valid_X_links = 0;
      load_X_from_W(fn, &ap->p2, hl->XX_fatlink[0], hl->XX_longlink[0]);
      for( inaik=1;inaik<fn->hl.n_naiks;inaik++ ) {
        for(dir=XUP;dir<=TUP;dir++) {
          FORALLSITES(i,s) {
            // add 2nd and 3rd path table sets
            add_su3_matrix( &( hl->XX_fatlink[inaik][dir][i] ),
                            &( hl->XX_fatlink[    0][dir][i] ),
                            &( hl->XX_fatlink[inaik][dir][i] ) );
            add_su3_matrix( &( hl->XX_longlink[inaik][dir][i] ),
                            &( hl->XX_longlink[    0][dir][i] ),
                            &( hl->XX_longlink[inaik][dir][i] ) );
          }
        }
      }
      hl->last_used_X_set = hl->current_X_set;
    }
    else {
      // 2nd path table set only, no other terms with Naik corrections
      load_X_from_W(fn, &ap->p2, hl->XX_fatlink[0], hl->XX_longlink[0]);
      hl->last_used_X_set = 0;
    }


    // set pointers to actual arrays
    for(dir=XUP;dir<=TUP;dir++) {
      hl->X_fatlink[dir]=hl->XX_fatlink[hl->last_used_X_set][dir];
      hl->X_longlink[dir]=hl->XX_longlink[hl->last_used_X_set][dir];
    }

    reorder_links(fn->fl.fat, hl->X_fatlink);
    reorder_links(fn->fl.lng, hl->X_longlink);
    fn->twist_in = OFF;

#ifdef DBLSTORE_FN
    load_fatbacklinks(fn);
    load_longbacklinks(fn);
#endif

    hl->valid_all_links = 1;

    /* Reapply the boundary twist */
    if(retwist)
      boundary_twist_fn(fn, ON);

    /* In case we are using this code in place of QOP code to manage
       the links, we have to refresh the QOP links now that the MILC
       fermion links may have changed. */
#if HAVE_QOP
    invalidate_qop_ferm_links_F(&fn->ql);
    invalidate_qop_ferm_links_D(&fn->ql);
#endif  
  }
  else {
//CHECK current_X_set AND SET POINTERS DEPENDING ON IT
//CALL reorder_links(), THEN load_...backlinks()
    if( hl->last_used_X_set != hl->current_X_set ) {
#ifdef MILC_GLOBAL_DEBUG
      node0_printf("load_ferm_links: choosing set %d, previous set %d\n",
		   hl->current_X_set, hl->last_used_X_set );
#endif
      
      // set pointers to actual arrays
      hl->last_used_X_set = hl->current_X_set;
      for(dir=XUP;dir<=TUP;dir++) {
        hl->X_fatlink[dir]=hl->XX_fatlink[hl->last_used_X_set][dir];
        hl->X_longlink[dir]=hl->XX_longlink[hl->last_used_X_set][dir];
      }

      reorder_links(fn->fl.fat, hl->X_fatlink);
      reorder_links(fn->fl.lng, hl->X_longlink);
      fn->twist_in = OFF;

#ifdef DBLSTORE_FN
      load_fatbacklinks(fn);
      load_longbacklinks(fn);
#endif

      /* Reapply the boundary twist */
      if(retwist)
	boundary_twist_fn(fn, ON);

    /* In case we are using this code in place of QOP code to manage
       the links, we have to refresh the QOP links now that the MILC
       fermion links may have changed. */
#if HAVE_QOP
      invalidate_qop_ferm_links_F(&fn->ql);
      invalidate_qop_ferm_links_D(&fn->ql);
#endif  
    }
    else {
#ifdef MILC_GLOBAL_DEBUG
    	node0_printf("load_ferm_links: nothing is done, sets are the same: %d\n",
    	  hl->current_X_set );
#endif
      ;
    }
  }

  return;

}

void 
load_ferm_links_dmdu0(ferm_links_t *fn){
  printf("load_ferm_links_dmdu0: NOT SUPPORTED YET FOR HISQ!\n");
  terminate(1);
}

// New routines for HISQ
static void  
load_U_from_site(ferm_links_t *fn){
  int dir,i; site *s;
  hisq_links_t *hl = &fn->hl;
  su3_matrix **U_link = hl->U_link;

  if(hl->valid_U_links)return;
  FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++) U_link[dir][i] = s->link[dir];
  hl->valid_U_links = 1;
  hl->phases_in_U = phases_in;
}

static void  
load_V_from_U(ferm_links_t *fn, ks_component_paths *ap1){
  hisq_links_t *hl = &fn->hl;
  su3_matrix **U_link = hl->U_link;
  su3_matrix **V_link = hl->V_link;

  if(  hl->valid_V_links )return;
  if( !hl->valid_U_links ){node0_printf("load_V_from_U: Link validity botched\n"); terminate(0);}
  load_fatlinks_hisq(U_link, ap1, V_link );
  hl->valid_V_links = hl->valid_U_links;
  hl->phases_in_V   = hl->phases_in_U;
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  int i, idir;
  site *s;
  complex cdetV;
  FORALLSITES(i,s) {
    for(idir=XUP;idir<=TUP;idir++) {
      if( lattice[i].on_step_V[idir] < global_current_time_step ) {
        lattice[i].on_step_V[idir] = global_current_time_step;
        cdetV = det_su3( &(V_link[idir][i]) );
        lattice[i].Vdet[idir] = cabs( &cdetV );
      }
    }
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
}

static void  
load_Y_from_V(ferm_links_t *fn, int umethod){
  int dir,i; site *s; su3_matrix tmat;
  hisq_links_t *hl = &fn->hl;
  su3_matrix **U_link = hl->U_link;
  su3_matrix **V_link = hl->V_link;
  su3_matrix **Y_unitlink = hl->Y_unitlink;

  if(  hl->valid_Y_links )return;
  if( !hl->valid_V_links ){node0_printf("load_Y_from_V: Link validity botched\n"); terminate(0);}
  switch(umethod){
  case UNITARIZE_NONE:
    //node0_printf("WARNING: UNITARIZE_NONE\n");
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++) 
      Y_unitlink[dir][i] = V_link[dir][i];
    
    hl->valid_Y_links = hl->valid_V_links;
    hl->phases_in_Y   = hl->phases_in_V;
    break;
  case UNITARIZE_APE:
    node0_printf("UNITARIZE_APE: derivative is not ready for this method\n"); 
    terminate(0);
    
    int nhits = 100; // this is just a guess
    Real tol = 1.0e-5; // this is just a guess
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      /* Use partially reunitarized link for guess */
      tmat = V_link[dir][i];
      reunit_su3(&tmat);
      project_su3(&tmat, &(V_link[dir][i]), nhits, tol);
      Y_unitlink[dir][i] = tmat;
    }
    break;
  case UNITARIZE_ROOT:
    //node0_printf("WARNING: UNITARIZE_ROOT is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS
    /* rephase (out) V_link array */
    rephase_field( V_link, OFF, &hl->phases_in_V );
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
      su3_unitarize( &( V_link[dir][i] ), &tmat );
      Y_unitlink[dir][i] = tmat;
    }
    hl->valid_Y_links = hl->valid_V_links;
    hl->phases_in_Y   = hl->phases_in_V;
    /* rephase (in) V_link array */
    rephase_field( V_link, ON, &hl->phases_in_V );
    /* initialize status and rephase (in) Y_unitlink array */
    rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;
  case UNITARIZE_RATIONAL:
    //node0_printf("WARNING: UNITARIZE_RATIONAL is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS
    /* rephase (out) V_link array */
    rephase_field( V_link, OFF, &hl->phases_in_V );
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
      su3_unitarize_rational( &( V_link[dir][i] ), &tmat );
      Y_unitlink[dir][i] = tmat;
    }
    hl->valid_Y_links = hl->valid_V_links;
    hl->phases_in_Y   = hl->phases_in_V;
    /* rephase (in) V_link array */
    rephase_field( V_link, ON, &hl->phases_in_V );
    /* rephase (in) Y_unitlink array */
    rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;
  case UNITARIZE_ANALYTIC:
    //node0_printf("WARNING: UNITARIZE_ANALYTIC is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS
    /* rephase (out) V_link array */
    rephase_field( V_link, OFF, &hl->phases_in_V );
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
      su3_unitarize_analytic_index( &( V_link[dir][i] ), &tmat, i, dir );
#else  /* HISQ_REUNITARIZATION_DEBUG */
      su3_unitarize_analytic( &( V_link[dir][i] ), &tmat );
#endif /* HISQ_REUNITARIZATION_DEBUG */
#else /* MILC_GLOBAL_DEBUG */
      su3_unitarize_analytic( &( V_link[dir][i] ), &tmat );
#endif /* MILC_GLOBAL_DEBUG */
      Y_unitlink[dir][i] = tmat;
    }
    hl->valid_Y_links = hl->valid_V_links;
    hl->phases_in_Y = hl->phases_in_V;
    /* rephase (in) V_link array */
    rephase_field( V_link, ON, &hl->phases_in_V );
    /* rephase (in) Y_unitlink array */
    rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;
  case UNITARIZE_STOUT:
    rephase_field( U_link, OFF, &hl->phases_in_U );
    rephase_field( V_link, OFF, &hl->phases_in_V );
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on SU(3) with "stout" procedure */
      stout_smear( &( Y_unitlink[dir][i] ), 
                   &( V_link[dir][i] ),
                   &( U_link[dir][i] ) );
    }
    hl->valid_Y_links = hl->valid_V_links;
    hl->phases_in_Y = hl->phases_in_V;
    /* rephase (in) U_link array */
    rephase_field( U_link, ON, &hl->phases_in_U );
    /* rephase (in) V_link array */
    rephase_field( V_link, ON, &hl->phases_in_V );
    /* rephase (in) Y_unitlink array */
    rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;
  case UNITARIZE_HISQ:
    node0_printf("UNITARIZE_HISQ: not ready!\n"); 
    terminate(0);
    break;
  default:
    node0_printf("Unknown unitarization method\n"); terminate(0);
  } /* umethod */
}

void  load_W_from_Y(ferm_links_t *fn, int umethod, int ugroup){
  int dir,i; site *s; su3_matrix tmat;
  hisq_links_t *hl = &fn->hl;
  su3_matrix **W_unitlink = hl->W_unitlink;
  su3_matrix **Y_unitlink = hl->Y_unitlink;

  if(  hl->valid_W_links )return;
  if( !hl->valid_Y_links ){node0_printf("load_W_from_Y: Link validity botched\n"); terminate(0);}

if( ugroup==UNITARIZE_SU3 ) {
  switch(umethod){
  case UNITARIZE_NONE:
    //node0_printf("WARNING: UNITARIZE_NONE\n");
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++) 
      W_unitlink[dir][i] = Y_unitlink[dir][i];
    hl->valid_W_links = hl->valid_Y_links;
    hl->phases_in_W = hl->phases_in_Y;
    break;
  case UNITARIZE_APE:
    node0_printf("UNITARIZE_APE: not ready\n"); 
    terminate(0);
    break;
  case UNITARIZE_ROOT:
    {
      complex cdet;
      //node0_printf("WARNING: SPECIAL UNITARIZE_ROOT is performed\n");
      /* rephase (out) Y_unitlink array */
      rephase_field( Y_unitlink, OFF, &hl->phases_in_Y );
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
	su3_spec_unitarize( &( Y_unitlink[dir][i] ), &tmat, &cdet );
	W_unitlink[dir][i] = tmat;
      }
      hl->valid_W_links = hl->valid_Y_links;
      hl->phases_in_W = hl->phases_in_Y;
      /* rephase (in) Y_unitlink array */
      rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
      /* rephase (in) W_unitlink array */
      rephase_field( W_unitlink, ON, &hl->phases_in_W );
      //printf("SPECIAL UNITARIZATION RESULT (ROOT)\n");
      //dumpmat( &( Y_unitlink[TUP][3] ) );
      //dumpmat( &( W_unitlink[TUP][3] ) );
    }
    break;
  case UNITARIZE_RATIONAL:
    {
      complex cdet;
      //node0_printf("WARNING: SPECIAL UNITARIZE_ROOT is performed\n");
      /* rephase (out) Y_unitlink array */
      rephase_field( Y_unitlink, OFF, &hl->phases_in_Y );
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
	su3_spec_unitarize( &( Y_unitlink[dir][i] ), &tmat, &cdet );
	W_unitlink[dir][i] = tmat;
      }
      hl->valid_W_links = hl->valid_Y_links;
      hl->phases_in_W = hl->phases_in_Y;
      /* rephase (in) Y_unitlink array */
      rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
      /* rephase (in) W_unitlink array */
      rephase_field( W_unitlink, ON, &hl->phases_in_W );
      //printf("SPECIAL UNITARIZATION RESULT (RATIONAL)\n");
      //dumpmat( &( Y_unitlink[TUP][3] ) );
      //dumpmat( &( W_unitlink[TUP][3] ) );
    }
    break;
  case UNITARIZE_ANALYTIC:
    {
      complex cdet;
      //node0_printf("WARNING: SPECIAL UNITARIZE_ROOT is performed\n");
      /* rephase (out) Y_unitlink array */
      rephase_field( Y_unitlink, OFF, &hl->phases_in_Y );
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
	su3_spec_unitarize_index( &( Y_unitlink[dir][i] ), &tmat, &cdet, i, dir );
#else  /* HISQ_REUNITARIZATION_DEBUG */
	su3_spec_unitarize( &( Y_unitlink[dir][i] ), &tmat, &cdet );
#endif /* HISQ_REUNITARIZATION_DEBUG */
#else  /* MILC_GLOBAL_DEBUG */
	su3_spec_unitarize( &( Y_unitlink[dir][i] ), &tmat, &cdet );
#endif /* MILC_GLOBAL_DEBUG */
	W_unitlink[dir][i] = tmat;
      }
      hl->valid_W_links = hl->valid_Y_links;
      hl->phases_in_W = hl->phases_in_Y;
      /* rephase (in) Y_unitlink array */
      rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
      /* rephase (in) W_unitlink array */
      rephase_field( W_unitlink, ON, &hl->phases_in_W );
//TEMP TEST
//{int ii,jj; complex cc; su3_matrix tmat;
//cc = ce_itheta ( 2.0*3.1415926/3.0 );
//c_scalar_mult_su3mat( &(hl->W_unitlink[XUP][0]), &cc, &tmat );
 //hl->W_unitlink[YUP][0] = tmat;
//printf("ACTION: TEMPTEST  %e  %e\n", cc.real, cc.imag);
//} //END TEMPTEST
      //printf("SPECIAL UNITARIZATION RESULT (RATIONAL)\n");
      //dumpmat( &( Y_unitlink[TUP][3] ) );
      //dumpmat( &( W_unitlink[TUP][3] ) );
    }
    break;
  case UNITARIZE_STOUT:
    {
      //node0_printf("WARNING: SPECIAL UNITARIZE_ROOT is performed\n");
      /* rephase (out) Y_unitlink array */
      rephase_field( Y_unitlink, OFF, &hl->phases_in_Y );
      FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
	/* QUICK FIX: copy Y to W */
        su3mat_copy( &( Y_unitlink[dir][i] ), &( W_unitlink[dir][i] ) );
      }
      hl->valid_W_links = hl->valid_Y_links;
      hl->phases_in_W = hl->phases_in_Y;
      /* rephase (in) Y_unitlink array */
      rephase_field( Y_unitlink, ON, &hl->phases_in_Y );
      /* rephase (in) W_unitlink array */
      rephase_field( W_unitlink, ON, &hl->phases_in_W );
    }
    break;
  case UNITARIZE_HISQ:
    node0_printf("UNITARIZE_HISQ: not ready!\n"); 
    terminate(0);
    break;
  default:
    node0_printf("Unknown unitarization method\n"); terminate(0);
  }
}
//  FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
//    /* simply copy Y links to W links */
//    su3mat_copy( &( Y_unitlink[dir][i] ), &( W_unitlink[dir][i] ) );
//  }
if ( ugroup==UNITARIZE_U3 ) {
  hl->valid_W_links = hl->valid_Y_links;
  hl->phases_in_W = hl->phases_in_Y;
}
}

/* Put links in site-major order for dslash */
static void 
reorder_links(su3_matrix *link4, su3_matrix *link[]){
  int i, dir;
  site *s;

  for(dir = 0; dir < 4; dir++){
    FORALLSITES(i,s){
      link4[4*i+dir] = link[dir][i];
    }
  }
}

void  load_X_from_W(ferm_links_t *fn, ks_component_paths *ap2,
  su3_matrix **X_fatlink, su3_matrix **X_longlink){
  hisq_links_t *hl = &fn->hl;
  su3_matrix **W_unitlink = hl->W_unitlink;
//  su3_matrix **X_fatlink  = hl->X_fatlink;
//  su3_matrix **X_longlink = hl->X_longlink;

  if(  hl->valid_X_links )return;
  if( !hl->valid_W_links ){node0_printf("load_X_from_W: Link validity botched\n"); terminate(0);}

  load_fatlinks_hisq(W_unitlink, ap2, X_fatlink );
  load_longlinks_hisq(W_unitlink, ap2, X_longlink );
//TEMP TEST
//{int ii,jj; complex cc; su3_matrix tmat;
//cc.real = 2.10;  cc.imag=0.00;
//c_scalar_mult_su3mat( &(hl->X_fatlink[XUP][0]), &cc, &tmat );
 //hl->X_fatlink[XUP][0] = tmat;
//printf("ACTION: TEMPTEST  %e  %e\n", cc.real, cc.imag);
//} //END TEMPTEST

  hl->valid_X_links = hl->valid_W_links;
  hl->phases_in_Xfat = hl->phases_in_W;
//  hl->valid_Xfat_mass = ap2->naik_mass;
  hl->phases_in_Xlong = hl->phases_in_W;
//  hl->valid_Xlong_mass = ap2->naik_mass;

  fn->fl.valid = hl->valid_X_links;
}

static void
invalidate_hisq_links(hisq_links_t *hl)
{
  hl->valid_U_links = 0;
  hl->valid_V_links = 0;
  hl->valid_W_links = 0;
  hl->valid_Y_links = 0;
  hl->valid_X_links = 0;
  hl->valid_all_links = 0;
  hl->valid_Xfat_mass = 0.;
  hl->valid_Xlong_mass = 0.;
}

// Invalidate all the links
void
invalidate_all_ferm_links(ferm_links_t *fn)
{
  fn->fl.valid = 0;

  invalidate_hisq_links(&fn->hl);
#if HAVE_QOP
  invalidate_qop_ferm_links_F(&fn->ql);
  invalidate_qop_ferm_links_D(&fn->ql);
#endif  
}

// Invalidate only the fat and long links
//void
//invalidate_fn_links(ferm_links_t *fn){
//  fn->fl.valid = 0;
//  fn->hl.valid_X_links = 0;
//  fn->hl.valid_all_links = 0;
//}

static void 
init_hisq_links(hisq_links_t *hl){
  int i,dir;

  invalidate_hisq_links(hl);

  hl->phases_in_U = OFF;
  hl->phases_in_V = OFF;
  hl->phases_in_W = OFF;
  hl->phases_in_Y = OFF;
  hl->phases_in_Xfat = OFF;
  hl->phases_in_Xlong = OFF;

  hl->current_X_set = 0;
  hl->last_used_X_set = 0;

  for(dir = 0; dir < 4; dir++){
    hl->U_link[dir] = NULL;
    hl->V_link[dir] = NULL;
    hl->Y_unitlink[dir] = NULL;
    hl->W_unitlink[dir] = NULL;
  }

  hl->n_naiks = 0;
  for(i = 0; i < MAX_NAIK; i++){
    hl->eps_naik[i] = 0.;
  }
}

void 
init_ferm_links(ferm_links_t *fn, ks_action_paths *ap){

  invalidate_all_ferm_links(fn);

  fn->mass  = 0.;
  fn->fl.fat = NULL;
  fn->fl.lng = NULL;
  fn->fl.fatback = NULL;
  fn->fl.lngback = NULL;
  fn->ap = ap;
  init_path_table(fn->ap);
  init_hisq_links(&fn->hl);
}

// print a value of link on site i in direction dir
void
look_at_link(ferm_links_t *fn, int *x, int dir){
  int i, idir;
  site *s;
  su3_matrix Usum, Vsum, Ysum, Xfatsum, Xlongsum;
  clear_su3mat( &Usum );
  clear_su3mat( &Vsum );
  clear_su3mat( &Ysum );
  clear_su3mat( &Xfatsum );
  clear_su3mat( &Xlongsum );
  // NEED TO convert site coordinates to site number
  i=0;
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Enter look_at_link fermion_links_hisq.c\n");
#endif
  printf("*** Looking at links:\n");
  printf("U links\n");
  dumpmat( &(fn->hl.U_link[dir][i]) );
  printf("V links\n");
  dumpmat( &(fn->hl.V_link[dir][i]) );
  printf("Y links\n");
  dumpmat( &(fn->hl.Y_unitlink[dir][i]) );
  printf("Xfat links\n");
  dumpmat( &(fn->hl.X_fatlink[dir][i]) );
  printf("Xlong links\n");
  dumpmat( &(fn->hl.X_longlink[dir][i]) );
  // global sum of links
  for(idir=XUP;idir<=TUP;idir++) {
    FORALLSITES(i,s) {
      add_su3_matrix( &Usum, &(fn->hl.U_link[idir][i]), &Usum );
      add_su3_matrix( &Vsum, &(fn->hl.V_link[idir][i]), &Vsum );
      add_su3_matrix( &Ysum, &(fn->hl.Y_unitlink[idir][i]), &Ysum );
      add_su3_matrix( &Xfatsum, &(fn->hl.X_fatlink[idir][i]), &Xfatsum );
      add_su3_matrix( &Xlongsum, &(fn->hl.X_longlink[idir][i]), &Xlongsum );
    }
  }
  printf("*** Looking at global sums of links:\n");
  printf("U links global sum\n");
  dumpmat( &Usum );
  printf("V links global sum\n");
  dumpmat( &Vsum );
  printf("Y links global sum\n");
  dumpmat( &Ysum );
  printf("Xfat links global sum\n");
  dumpmat( &Xfatsum );
  printf("Xlong links global sum\n");
  dumpmat( &Xlongsum );
#ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
  printf("Exit  look_at_link fermion_links_hisq.c\n");
#endif
}

// end new routines

