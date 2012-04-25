/**************** fermion_links_fn_load_milc.c *****************************/
/* MILC Version 7 */

/* Methods for the fn_links_t structure  */
/* Compute and load the fat and long links. */

#include "generic_ks_includes.h"
#define IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#include <quark_action.h>
#include "../include/info.h"

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

/*-------------------------------------------------------------------*/
/* Create/destroy Naik links                                         */
/*-------------------------------------------------------------------*/
void
load_lnglinks(info_t *info, su3_matrix *lng, ks_component_paths *p,
	      su3_matrix *links ) {
  register int i;
  
  int ipath,dir;
  int disp[4];
  int num_q_paths = p->num_q_paths;
  Q_path *q_paths = p->q_paths;
  register su3_matrix *long1;
  su3_matrix *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_lnglinks";
  double dtime = -dclock();

  if( phases_in != 1){
    node0_printf("BOTCH: %s needs phases in\n",myname);
    terminate(0);
  }

  staple = create_m_special();
  tempmat1 = create_m_special();

  for (dir=XUP; dir<=TUP; dir++){ /* loop over longlink directions */
    /* set longlink to zero */
    FORALLFIELDSITES(i){
      long1 = lng + 4*i +dir;
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

//	path_product_field( q_paths[ipath].dir, q_paths[ipath].length, 
//			    tempmat1, links );
	path_product_fields( links, q_paths[ipath].dir, q_paths[ipath].length, 
			     tempmat1 );
	FORALLFIELDSITES(i){
	  su3_adjoint( &tempmat1[i], &staple[i] );
	  long1 = lng + 4*i + dir;
          scalar_mult_add_su3_matrix( long1,
	    &staple[i], -q_paths[ipath].coeff, long1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */

  } /* loop over directions */


  destroy_m_special(staple); staple = NULL;
  destroy_m_special(tempmat1); tempmat1 = NULL;


  dtime += dclock();
  info->final_sec = dtime;
  info->final_flop = 1728.*volume/numnodes();  /* (formerly 1804) */

}  /* load_lnglinks() */

/*-------------------------------------------------------------------*/
/* Create/destroy fat links                                          */
/*-------------------------------------------------------------------*/
/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED_FATTENING - C. DeTar
 * Fatlinks:       61632 for load_fatlinks
 */

/* KS phases and APBC must be in the links. See long comment at 
   end of fermion_forcee_eo_milc.c */

void
load_fatlinks_cpu(info_t *info, su3_matrix *fat, ks_component_paths *p, 
		  su3_matrix *links){
  register int i;
  int dir;
  register su3_matrix *fat1;
  su3_matrix *staple = NULL, *tempmat1 = NULL;
  char myname[] = "load_fatlinks_cpu";

#ifdef ASQ_OPTIMIZED_FATTENING
  int  nu,rho,sig ;
  Real one_link;
#else
  int ipath;
  int disp[4];
  int num_q_paths = p->num_q_paths;
  Q_path *q_paths = p->q_paths;
#endif

  double dtime = -dclock();

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

#ifndef  ASQ_OPTIMIZED_FATTENING   /* general case code */
  for (dir=XUP; dir<=TUP; dir++){ /* loop over fatlink directions */
    /* set fatlink to zero */
    FORALLFIELDSITES(i){
      fat1 = fat + 4*i + dir;
      clear_su3mat( fat1 );
    }
    
    /* loop over paths, checking for ones with total displacement 1*dir */
    for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
	/* compute total displacement of path */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	for( i=0; i<q_paths[ipath].length; i++){
	  if( GOES_FORWARDS(q_paths[ipath].dir[i]) )
	    disp[        q_paths[ipath].dir[i]  ]++;
	  else
	    disp[OPP_DIR(q_paths[ipath].dir[i]) ]--;
	}
	for( disp[dir]+=1,i=XUP; i<=TUP; i++)if(disp[i]!=0)break;
	if( i<=TUP )continue;  /* skip if path doesn't go to right place */
/**printf("dir = %d, found a path:  ",dir);
for(j=0;j<q_paths.[ipath].length;j++)printf("\t%d", q_paths[ipath].dir[j]);
printf("\n");**/

//	path_product( q_paths[ipath].dir, q_paths[ipath].length, tempmat1 );
//	path_product_field( q_paths[ipath].dir, q_paths[ipath].length, 
//			    tempmat1, links );
	path_product_fields( links, q_paths[ipath].dir, q_paths[ipath].length, 
			     tempmat1 );
	FORALLFIELDSITES(i){
	  su3_adjoint( &tempmat1[i], &staple[i] );
	  fat1 = fat +  4*i + dir;
          scalar_mult_add_su3_matrix( fat1,
	    &staple[i], -q_paths[ipath].coeff, fat1 );
		/* minus sign in coeff. because we used backward path*/
	}
    } /* ipath */
  } /* loop over directions */
#else	/* ASQ_OPTIMIZED_FATTENING, for Asq and Asqtad actions */
/*  Optimized fattening code for the Asq and Asqtad actions.           *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lepage term. *
 * Path 1 is the Naik term.                                            */
 
 /* to fix up the Lepage term, included by a trick below */
 one_link = (p->act_path_coeff.one_link - 6.0*p->act_path_coeff.lepage);
 
 for (dir=XUP; dir<=TUP; dir++){
   FORALLFIELDSITES(i) /* Intialize fat links with c_1*U_\mu(x) */
     {
       fat1 = fat +  4*i + dir;
       scalar_mult_su3_matrix(links + 4*i + dir, one_link,
			      fat1 );
     }
   /* Skip the rest of the calculation if the remaining coefficients vanish */
   if( p->act_path_coeff.three_staple == 0.0 &&
       p->act_path_coeff.lepage == 0.0 &&
       p->act_path_coeff.five_staple == 0.0)continue;

   for(nu=XUP; nu<=TUP; nu++) if(nu!=dir)
     {
//       compute_gen_staple_site(staple,dir,nu,F_OFFSET(link[dir]),
//			       *t_fl, act_path_coeff.three_staple);

       compute_gen_staple_field(staple, dir, nu, links + dir, 4,
				fat, p->act_path_coeff.three_staple, links);
       /* The Lepage term */
       /* Note this also involves modifying c_1 (above) */
       compute_gen_staple_field(NULL, dir, nu, staple, 1,
				fat, p->act_path_coeff.lepage, links);
       for(rho=XUP; rho<=TUP; rho++) if((rho!=dir)&&(rho!=nu))
	 {
	   compute_gen_staple_field( tempmat1, dir, rho, staple, 1,
				     fat, p->act_path_coeff.five_staple, links);
	   for(sig=XUP; sig<=TUP; sig++)
	     if((sig!=dir)&&(sig!=nu)&&(sig!=rho))
	       {
		 compute_gen_staple_field(NULL,dir,sig,tempmat1, 1,
				  fat, p->act_path_coeff.seven_staple, links);
	       } /* sig */
	 } /* rho */
     } /* nu */
 }/* dir */  
#endif

 special_free(staple);  staple = NULL;
 special_free(tempmat1); tempmat1 = NULL;

 dtime += dclock();
 info->final_sec += dtime;
 info->final_flop = 61632.*volume/numnodes();
 if( p->act_path_coeff.three_staple == 0.0 &&
     p->act_path_coeff.lepage == 0.0 &&
     p->act_path_coeff.five_staple == 0.0)
   info->final_flop = 72.*volume/numnodes();

}  /* load_fatlinks_cpu */

/*-------------------------------------------------------------------*/
/* Fill in fn_links_t structure                                      */
/*-------------------------------------------------------------------*/

/* This procedure is meant to be called only by fermion_links*.c
   routines. The public API is create_fermion_links and
   restore_fermion_links.  */

void
load_fn_links(info_t *info, fn_links_t *fn, ks_action_paths *ap, 
	      su3_matrix *links, int want_back){
  ks_component_paths *p = &ap->p;
  double final_flop = 0;
  double dtime = -dclock();

  //  fn->fat = create_fatlinks();
  load_fatlinks(info, fn->fat, p, links);
  final_flop += info->final_flop;

  //  fn->lng = create_lnglinks();
  load_lnglinks(info, fn->lng, p, links);
  final_flop += info->final_flop;

  if(want_back)
    load_fn_backlinks(fn);
  else
    destroy_fn_backlinks(fn);

  dtime += dclock();
  info->final_sec = dtime;
  info->final_flop = final_flop;
}
