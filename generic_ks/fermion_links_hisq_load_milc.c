/**************** fermion_links_hisq_load_milc.c *****************************/
/* (Formerly fermion_links_hisq.c AND fermion_links_helpers_hisq.c) */
/* MILC Version 7 */
/* Routines used by all of the fermion_links*.c files */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include <quark_action.h>
#include "../include/fermion_links.h"
#include "../include/info.h"
#include <string.h>

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

/*-------------------------------------------------------------------*/
/* Special memory allocations for field with one su3_matrix per site */
/*-------------------------------------------------------------------*/

#if 0
static su3_matrix *
create_m_special(void){
  char myname[] = "create_m_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*sizeof(su3_matrix));

  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  return m;
}

static void
destroy_m_special(su3_matrix *m){
  special_free(m);
}
#endif

/*--------------------------------------------------------------------*/

static su3_matrix *create_m4_special(void){
  char myname[] = "create_m4_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));

  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  return m;
}

static void destroy_m4_special(su3_matrix *m){
  special_free(m);
}

/*--------------------------------------------------------------------*/
// higher precison su3_matrix dump

void dumpmat_hp( su3_matrix *m ){
int i,j;
    for(i=0;i<3;i++){
	for(j=0;j<3;j++)printf("(%.15f,%.15f)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}

/*--------------------------------------------------------------------*/
// Some notation:
// U -- original link, SU(3), copied to "field" from "site"
// V -- after 1st level of smearing, non-SU(3)
// Y -- unitarized, U(3)
// W -- special unitarized, SU(3)
// X -- after 2nd level of smearing, non-SU(3)

/*--------------------------------------------------------------------*/

static void  
load_U_from_field(info_t *info, hisq_auxiliary_t *aux, su3_matrix *links){
  //  int dir,i;
  su3_matrix * U_link = aux->U_link;

//  FORALLFIELDSITES(i){
//    for(dir=XUP;dir<=TUP;dir++) 
//      U_link[dir][i] = links[4*i + dir];
//  }
  memcpy(U_link, links, 4*sizeof(su3_matrix)*sites_on_node);
}

/*--------------------------------------------------------------------*/

static void  
load_V_from_U(info_t *info, hisq_auxiliary_t *aux, ks_component_paths *ap1){
  su3_matrix *U_link = aux->U_link;
  su3_matrix *V_link = aux->V_link;

  load_fatlinks(info, V_link, ap1, U_link);

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  int i, idir;
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

/*--------------------------------------------------------------------*/

#ifdef USE_REUNIT_GPU

static void 
load_Y_from_V(info_t *info, hisq_auxiliary_t *aux, int umethod){
  int dir, i;
  su3_matrix *V_link = aux->V_link;
  su3_matrix *Y_unitlink = aux->Y_unitlink;
  const int num_its = 20;
  double dtime = -dclock();
  int svd_calls = 0;

  /* Until we support performance statistics here */
  info->final_flop = 0.0;

  switch(umethod){
    case UNITARIZE_NONE:
      FORALLFIELDSITES(i) for(dir=XUP; dir<=TUP; ++dir)
	      Y_unitlink[4*i+dir] = V_link[4*i+dir];
	    break;

    case UNITARIZE_ANALYTIC:
	    unitarize_field_analytic_gpu(Y_unitlink, V_link);  // this needs to be defined elsewhere
      break;

    case UNITARIZE_ITERATIVE:
            unitarize_field_iterative_gpu(Y_unitlink, V_link, num_its);
      break;

    default:
      node0_printf("Unknown unitarization method\n"); terminate(0);
  }

  dtime += dclock();
  info->final_sec = dtime;
  INFO_HISQ_SVD_COUNTER(info) = svd_calls;
  return;
}

#else

/*--------------------------------------------------------------------*/

static void  
load_Y_from_V(info_t *info, hisq_auxiliary_t *aux, int umethod){
  int dir,i; 
  su3_matrix tmat;
  su3_matrix *U_link = aux->U_link;
  su3_matrix *V_link = aux->V_link;
  su3_matrix *Y_unitlink = aux->Y_unitlink;
  /* Temporary.  FIX THIS! */
  int r0[4] = {0, 0, 0, 0};

  double dtime;
  dtime=-dclock();

  int nhits = 100; // this is just a guess
  Real tol = 1.0e-5; // this is just a guess

  info->final_flop = 0.;

  switch(umethod){

  case UNITARIZE_NONE:
    //node0_printf("WARNING: UNITARIZE_NONE\n");
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++) 
      Y_unitlink[4*i+dir] = V_link[4*i+dir];
    break;

  case UNITARIZE_APE:
    node0_printf("UNITARIZE_APE: derivative is not ready for this method\n"); 
    terminate(0);
    
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      /* Use partially reunitarized link for guess */
      tmat = V_link[4*i+dir];
      reunit_su3(&tmat);
      project_su3(&tmat, &(V_link[4*i+dir]), nhits, tol);
      Y_unitlink[4*i+dir] = tmat;
    }
    break;

  case UNITARIZE_ROOT:
    //node0_printf("WARNING: UNITARIZE_ROOT is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS

    /* rephase (out) V_link array */
    rephase_field_offset( V_link, OFF, NULL , r0);
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
      u3_unitarize( &( V_link[4*i+dir] ), &tmat );
      Y_unitlink[4*i+dir] = tmat;
    }

    /* rephase (in) Y_unitlink array */
    rephase_field_offset( Y_unitlink, ON, NULL , r0);
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;

  case UNITARIZE_RATIONAL:
    //node0_printf("WARNING: UNITARIZE_RATIONAL is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS

    /* rephase (out) V_link array */
    rephase_field_offset( V_link, OFF, NULL , r0);
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
      u3_unitarize_rational( &( V_link[4*i+dir] ), &tmat );
      Y_unitlink[4*i+dir] = tmat;
    }
    rephase_field_offset( V_link, ON, NULL , r0);

    /* rephase (in) Y_unitlink array */
    rephase_field_offset( Y_unitlink, ON, NULL , r0);
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );
    break;

  case UNITARIZE_ANALYTIC:
    //node0_printf("WARNING: UNITARIZE_ANALYTIC is performed\n");
    // REPHASING IS NOT NEEDED IN THIS ROUTINE BUT THIS
    // HAS TO BE CHECKED FIRST FOR NONTRIVIAL V_link ARRAYS
    /* rephase (out) V_link array */
    rephase_field_offset( V_link, OFF, NULL , r0);
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on U(3) */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
      u3_unitarize_analytic_index( &( V_link[4*i+dir] ), &tmat, i, dir );
#else  /* HISQ_REUNITARIZATION_DEBUG */
      u3_unitarize_analytic( info, &( V_link[4*i+dir] ), &tmat );
#endif /* HISQ_REUNITARIZATION_DEBUG */
#else /* MILC_GLOBAL_DEBUG */
      u3_unitarize_analytic( info, &( V_link[4*i+dir] ), &tmat );
#endif /* MILC_GLOBAL_DEBUG */
      Y_unitlink[4*i+dir] = tmat;
    }
    /* rephase (in) V_link array */
    rephase_field_offset( V_link, ON, NULL , r0);

    /* rephase (in) Y_unitlink array */
    rephase_field_offset( Y_unitlink, ON, NULL , r0);
    //printf("UNITARIZATION RESULT\n");
    //dumpmat( &( V_link[TUP][3] ) );
    //dumpmat( &( Y_unitlink[TUP][3] ) );

    break;

  case UNITARIZE_STOUT:

    rephase_field_offset( U_link, OFF, NULL , r0);
    rephase_field_offset( V_link, OFF, NULL , r0);
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      /* unitarize - project on SU(3) with "stout" procedure */
      stout_smear( &( Y_unitlink[4*i+dir] ), 
                   &( V_link[4*i+dir] ),
                   &( U_link[4*i+dir] ) );
    }
    rephase_field_offset( U_link, ON, NULL , r0);
    rephase_field_offset( V_link, ON, NULL , r0);

    /* rephase (in) Y_unitlink array */
    rephase_field_offset( Y_unitlink, ON, NULL , r0);
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

  dtime += dclock();
  info->final_sec = dtime;
}

#endif

/*--------------------------------------------------------------------*/

void 
load_W_from_Y(info_t *info, hisq_auxiliary_t *aux, int umethod, int ugroup){
  su3_matrix *W_unitlink = aux->W_unitlink;
  su3_matrix *Y_unitlink = aux->Y_unitlink;
  int dir,i; su3_matrix tmat;
  /* Temporary.  FIX THIS! */
  int r0[4] = {0, 0, 0, 0};
  char myname[] = "load_W_from_Y";

  info->final_flop = 0; /* Currently not counted if we do use SU(3) */

  /* If we have already set the W pointers equal to the V
     pointers,  no copy is needed. */

  if(aux->WeqY)return;

  if( ugroup != UNITARIZE_SU3 ) {
    node0_printf("%s: Unrecognized ugroup value\n", myname);
    terminate(1);
  }

  switch(umethod){
    
  case UNITARIZE_NONE:
    //node0_printf("WARNING: UNITARIZE_NONE\n");
    FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++) 
      W_unitlink[4*i+dir] = Y_unitlink[4*i+dir];
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
      rephase_field_offset( Y_unitlink, OFF, NULL , r0);
      FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
	su3_spec_unitarize( &( Y_unitlink[4*i+dir] ), &tmat, &cdet );
	W_unitlink[4*i+dir] = tmat;
      }
      rephase_field_offset( Y_unitlink, ON, NULL , r0);
      
      /* rephase (in) W_unitlink array */
      rephase_field_offset( W_unitlink, ON, NULL , r0);
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
      rephase_field_offset( Y_unitlink, OFF, NULL , r0);
      FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
	su3_spec_unitarize( &( Y_unitlink[4*i+dir] ), &tmat, &cdet );
	W_unitlink[4*i+dir] = tmat;
      }
      rephase_field_offset( Y_unitlink, ON, NULL , r0);
      
      /* rephase (in) W_unitlink array */
      rephase_field_offset( W_unitlink, ON, NULL , r0);
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
      rephase_field_offset( Y_unitlink, OFF, NULL , r0);
      FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
	/* special unitarize - project U(3) on SU(3)
	   CAREFUL WITH FERMION PHASES! */
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
	su3_spec_unitarize_index( &( Y_unitlink[dir][i] ), &tmat, 
				  &cdet, i, dir );
#else  /* HISQ_REUNITARIZATION_DEBUG */
	su3_spec_unitarize( &( Y_unitlink[dir][i] ), &tmat, &cdet );
#endif /* HISQ_REUNITARIZATION_DEBUG */
#else  /* MILC_GLOBAL_DEBUG */
	su3_spec_unitarize( &( Y_unitlink[4*i+dir] ), &tmat, &cdet );
#endif /* MILC_GLOBAL_DEBUG */
	W_unitlink[4*i+dir] = tmat;
      }
      rephase_field_offset( Y_unitlink, ON, NULL , r0);
      
      /* rephase (in) W_unitlink array */
      rephase_field_offset( W_unitlink, ON, NULL , r0);
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
      rephase_field_offset( Y_unitlink, OFF, NULL , r0);
      FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
	/* QUICK FIX: copy Y to W */
	su3mat_copy( &( Y_unitlink[4*i+dir] ), &( W_unitlink[4*i+dir] ) );
      }
      rephase_field_offset( Y_unitlink, ON, NULL , r0);
      
      /* rephase (in) W_unitlink array */
      rephase_field_offset( W_unitlink, ON, NULL , r0);
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

#if 0
/*--------------------------------------------------------------------*/

/* Put links in site-major order for dslash */
static void 
reorder_links(su3_matrix *link4, su3_matrix *link[]){
  int i, dir;

  for(dir = 0; dir < 4; dir++){
    FORALLFIELDSITES(i){
      link4[4*i+dir] = link[4*dir][i];
    }
  }
}
#endif
/*--------------------------------------------------------------------*/
/* Smear the W links to make the HISQ links.  Result in fn            */
/*--------------------------------------------------------------------*/

static void  
load_X_from_W(info_t *info, fn_links_t *fn, hisq_auxiliary_t *aux,
	      ks_component_paths *ap){

  su3_matrix *fat = get_fatlinks(fn);
  su3_matrix *lng = get_lnglinks(fn);
  double final_flop = 0.0;

  load_fatlinks(info, fat, ap, aux->W_unitlink );
  final_flop += info->final_flop;
  load_lnglinks(info, lng, ap, aux->W_unitlink );
  final_flop += info->final_flop;

  info->final_flop = final_flop;

}

/*--------------------------------------------------------------------*/
/* Make the various auxiliary fat links.
	U = copy of links in site structure
	V = once smeared
	Y = projected onto U3
	W = projected onto SU3
*/

static void 
load_hisq_aux_links_cpu(info_t *info, ks_action_paths_hisq *ap, hisq_auxiliary_t *aux,
			su3_matrix *links){
  char myname[] = "load_hisq_aux_links";
  double final_flop = 0.0;

  if(ap == NULL){
    printf("%s(%d): KS action paths not initialized\n", myname, this_node);
    terminate(1);
  }

  load_U_from_field(info, aux, links );

  load_V_from_U(info, aux, &ap->p1);
  final_flop += info->final_flop;

  load_Y_from_V(info, aux, ap->umethod );
  final_flop += info->final_flop;

  load_W_from_Y(info, aux, ap->umethod, ap->ugroup);
  final_flop += info->final_flop;

  info->final_flop = final_flop;
}

/*--------------------------------------------------------------------*/

static void 
load_hisq_fn_links(info_t *info, fn_links_t **fn, fn_links_t *fn_deps,
		   hisq_auxiliary_t *aux, ks_action_paths_hisq *ap, 
		   su3_matrix *links, int want_deps, int want_back){

  //char myname[] = "load_hisq_fn_links";
  int inaik;
  int n_naiks = ap->n_naiks;
  double *eps_naik = ap->eps_naik;
  double final_flop = 0.0;
  
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


  if(n_naiks > 1 ) {
    // 3rd path table set
    load_X_from_W(info, fn[0], aux, &ap->p3);
    final_flop += info->final_flop;
    if(want_deps)
      copy_fn(fn[0], fn_deps);
    for( inaik = 1; inaik < n_naiks; inaik++ ){
      scalar_mult_fn( fn[0], eps_naik[inaik], fn[inaik] );
      final_flop += 18.*volume/numnodes();
    }

    // 2nd path table set
    load_X_from_W(info, fn[0], aux, &ap->p2);
    final_flop += info->final_flop;
    for( inaik = 1; inaik < n_naiks; inaik++ )
      add_fn( fn[inaik], fn[0], fn[inaik] );
      final_flop += 18.*volume/numnodes();
  }
  else {
    // 2nd path table set only, no other terms with Naik corrections
    load_X_from_W(info, fn[0], aux, &ap->p2);
    final_flop += info->final_flop;
    if(want_deps){
      load_X_from_W(info, fn_deps, aux, &ap->p3);
      final_flop += info->final_flop;
    }
  }

  /* Move up the back links */
  if(want_back){
    for( inaik = 0; inaik < n_naiks; inaik++)
      load_fn_backlinks(fn[inaik]);
    if(want_deps)
      load_fn_backlinks(fn_deps);
  }

  info->final_flop = final_flop;
}

/*-------------------------------------------------------------------*/
/* Alexei's API follows...                                             */
/*-------------------------------------------------------------------*/

// // print a value of link on site i in direction dir
// void
// look_at_link(ferm_links_t *fn, int *x, int dir){
//   int i, idir;
//   site *s;
//   su3_matrix Usum, Vsum, Ysum, Xfatsum, Xlongsum;
//   clear_su3mat( &Usum );
//   clear_su3mat( &Vsum );
//   clear_su3mat( &Ysum );
//   clear_su3mat( &Xfatsum );
//   clear_su3mat( &Xlongsum );
//   // NEED TO convert site coordinates to site number
//   i=0;
// #ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//   printf("Enter look_at_link fermion_links_hisq.c\n");
// #endif
//   printf("*** Looking at links:\n");
//   printf("U links\n");
//   dumpmat( &(fn->hl.U_link[dir][i]) );
//   printf("V links\n");
//   dumpmat( &(fn->hl.V_link[dir][i]) );
//   printf("Y links\n");
//   dumpmat( &(fn->hl.Y_unitlink[dir][i]) );
//   printf("Xfat links\n");
//   dumpmat( &(fn->hl.X_fatlink[dir][i]) );
//   printf("Xlong links\n");
//   dumpmat( &(fn->hl.X_longlink[dir][i]) );
//   // global sum of links
//   for(idir=XUP;idir<=TUP;idir++) {
//     FORALLSITES(i,s) {
//       add_su3_matrix( &Usum, &(fn->hl.U_link[idir][i]), &Usum );
//       add_su3_matrix( &Vsum, &(fn->hl.V_link[idir][i]), &Vsum );
//       add_su3_matrix( &Ysum, &(fn->hl.Y_unitlink[idir][i]), &Ysum );
//       add_su3_matrix( &Xfatsum, &(fn->hl.X_fatlink[idir][i]), &Xfatsum );
//       add_su3_matrix( &Xlongsum, &(fn->hl.X_longlink[idir][i]), &Xlongsum );
//     }
//   }
//   printf("*** Looking at global sums of links:\n");
//   printf("U links global sum\n");
//   dumpmat( &Usum );
//   printf("V links global sum\n");
//   dumpmat( &Vsum );
//   printf("Y links global sum\n");
//   dumpmat( &Ysum );
//   printf("Xfat links global sum\n");
//   dumpmat( &Xfatsum );
//   printf("Xlong links global sum\n");
//   dumpmat( &Xlongsum );
// #ifdef AB_DEBUG_ENTRY_EXIT_ROUTINES
//   printf("Exit  look_at_link fermion_links_hisq.c\n");
// #endif
// }

/*-------------------------------------------------------------------*/
/* Public API follows...                                             */
/*-------------------------------------------------------------------*/

/*-----------------------------------------------------*/
/* Create/destroy the hisq_auxiliary_t structure      */
/*-----------------------------------------------------*/

/* Create all but the last asqtad smearing */
static hisq_auxiliary_t *
create_hisq_auxiliary_t(ks_action_paths_hisq *ap, su3_matrix *links){

  hisq_auxiliary_t *aux;
  char myname[] = "create_hisq_auxiliary_t";

  if(ap == NULL){
    printf("%s(%d): KS action paths not initialized\n", myname, this_node);
    terminate(1);
  }

  aux = (hisq_auxiliary_t *)malloc(sizeof(hisq_auxiliary_t));
  if(aux == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  
  aux->U_link = create_m4_special();
  aux->V_link = create_m4_special();
  aux->Y_unitlink = create_m4_special();
  if(ap->ugroup == UNITARIZE_U3){
    /* If we are doing only U(3), we equate the W and Y links */
    aux->W_unitlink = aux->Y_unitlink;
    aux->WeqY = 1;
  } else {
    aux->W_unitlink = create_m4_special();
    aux->WeqY = 0;
  }

//  /* Copy the gauge field */
//  load_U_from_field(aux, links );
//
//  /* First level of smearing */
//  load_V_from_U(aux, &ap->p1);
//
//  /* Project onto U(3) */
//  load_Y_from_V(aux, ap->umethod);
//
//  /* Project onto SU(3) if requested */
//  load_W_from_Y(aux, ap->umethod, ap->ugroup);
  
  return aux;
}

void
destroy_hisq_auxiliary_t(hisq_auxiliary_t *aux){

  if(aux == NULL)return;
  destroy_m4_special(aux->W_unitlink);
  if(aux->WeqY)
    aux->Y_unitlink = NULL;
  else
    destroy_m4_special(aux->Y_unitlink);
  destroy_m4_special(aux->V_link);
  destroy_m4_special(aux->U_link);

  free(aux);
}

/*-------------------------------------------------------------------*/
/* Create/destroy HISQ fn links                                      */
/*-------------------------------------------------------------------*/

void 
create_hisq_links_milc(info_t *info, fn_links_t **fn, fn_links_t **fn_deps,
		       hisq_auxiliary_t **aux, ks_action_paths_hisq *ap, 
		       su3_matrix *links, int want_deps, int want_back){
  //char myname[] = "create_hisq_links_milc";

  int n_naiks = ap->n_naiks;
  int i;
  double final_flop = 0.;
  double dtime = -dclock();

  *aux = create_hisq_auxiliary_t(ap, links);
  
  load_hisq_aux_links(info, ap, *aux, links);
  final_flop += info->final_flop;
  
  for(i = 0; i < n_naiks; i++)
    fn[i] = create_fn_links();

  if(want_deps)
    *fn_deps = create_fn_links();
  else
    *fn_deps = NULL;

  load_hisq_fn_links(info, fn, *fn_deps, *aux, ap, links, 
		     want_deps, want_back);
  final_flop += info->final_flop;

  dtime += dclock();
  info->final_sec = dtime;
}

void
destroy_hisq_links_milc(ks_action_paths_hisq *ap, hisq_auxiliary_t *aux, 
			fn_links_t **fn, fn_links_t *fn_deps){

  //char myname[] = "destroy_hisq_links_milc";

  int n_naiks = ap->n_naiks;
  int i;

  destroy_hisq_auxiliary_t(aux);
  
  /* Destroy and free the structures */
  for(i = 0; i < n_naiks; i++){
    destroy_fn_links(fn[i]);
    fn[i] = NULL;
  }

  destroy_fn_links(fn_deps);
}

