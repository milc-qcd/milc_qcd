/********************** general_staple.c *****************************/
/* MILC Version 7 */

#include "generic_includes.h"	/* definitions files and prototypes */

/* Computes the staple :
                   mu
                +-------+
          nu	|	|
		|	|
		X	X
   Where the mu link can be any su3_matrix. The result is saved in staple.
   if staple==NULL then the result is not saved.
   It also adds the computed staple to fatlink[mu] with weight coef.
*/

/* Here "link" represents the "mu" link and it is stored as "stride"
   matrices per site.  The nu "links" are the original gauge field
   links, stored as four matrices per site. The result is accumulated
   in "fatlink", stored as four matrices per site */

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static msg_tag *
start_gather_field_strided(void *field, int stride, int size,
			   int index, int parity, char **dest){
  msg_tag *mt;
  mt = declare_strided_gather( field, stride, size, index, parity, dest );
  prepare_gather(mt);
  do_gather(mt);
  return mt;
}

void 
compute_gen_staple_field(su3_matrix *staple, int mu, int nu, 
			 su3_matrix *link, int stride,
			 su3_matrix *fatlink, Real coef,
			 su3_matrix *links) {
  su3_matrix tmat1,tmat2;
  msg_tag *mtag0,*mtag1;
  su3_matrix *tempmat = NULL;
  register site *s ;
  register int i ;
  register su3_matrix *fat1;


  /* Upper staple */
  mtag0 = start_gather_field_strided( link, stride*sizeof(su3_matrix),
				      sizeof(su3_matrix), nu, EVENANDODD, 
				      gen_pt[0] );
//  mtag1 = start_gather_site( F_OFFSET(link[nu]), sizeof(su3_matrix), mu, 
//			EVENANDODD, gen_pt[1] );
  mtag1 = start_gather_field_strided( links + nu, 4*sizeof(su3_matrix),
				      sizeof(su3_matrix), mu, 
				      EVENANDODD, gen_pt[1] );
  wait_gather(mtag0);
  wait_gather(mtag1);
  
  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      //      mult_su3_nn( &(s->link[nu]), &tmat1, &staple[i] );
      mult_su3_nn( links + 4*i + nu, &tmat1, staple + i );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)gen_pt[0][i],
		   (su3_matrix *)gen_pt[1][i], &tmat1 );
      //      mult_su3_nn( &(s->link[nu]), &tmat1, &tmat2 );
      mult_su3_nn( links + 4*i + nu, &tmat1, &tmat2 );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix(fat1, &tmat2, coef, fat1) ;
    }
  }
  cleanup_gather(mtag0);
  cleanup_gather(mtag1);

  /* lower staple */
  tempmat = (su3_matrix *)special_alloc( sites_on_node*sizeof(su3_matrix) );
//  mtag0 = start_gather_site( F_OFFSET(link[nu]),
//			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  mtag0 = start_gather_field_strided( links + nu, 4*sizeof(su3_matrix),
			sizeof(su3_matrix), mu, EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);
  FORALLSITES(i,s){
    //    mult_su3_an( &(s->link[nu]),&link[i], &tmat1 );
    mult_su3_an( links + 4*i + nu, link + stride*i, &tmat1 );
    mult_su3_nn( &(tmat1),(su3_matrix *)gen_pt[0][i], &(tempmat[i]) );
  }
  cleanup_gather(mtag0);
  mtag0 = start_gather_field( tempmat, sizeof(su3_matrix),
			      OPP_DIR(nu), EVENANDODD, gen_pt[0] );
  wait_gather(mtag0);

  if(staple!=NULL){/* Save the staple */
    FORALLSITES(i,s){
      add_su3_matrix( staple + i,(su3_matrix *)gen_pt[0][i], 
		      staple + i );
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1, staple + i, coef, fat1 );
    }
  }
  else{ /* No need to save the staple. Add it to the fatlinks */
    FORALLSITES(i,s){
      fat1 = &(fatlink[4*i+mu]);
      scalar_mult_add_su3_matrix( fat1,
				 (su3_matrix *)gen_pt[0][i], coef, 
				 fat1 );
    }
  }

  special_free(tempmat); tempmat = NULL;
  cleanup_gather(mtag0);
} /* compute_gen_staple_field */
