/******* dslash_eo.c - dslash for improved KS fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- improved.  

   D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions.

   Stupid dslash routine that follows all the paths.  Should optimize
   by precomputing sums of paths to all the displacements.

   Use dslash_fn for actions that involve only +X and +X+X+X
   couplings.
*/
/* C. DeTar 3/05 Code split from quark_stuff?.c */
/* UMH 2/11 Fix not-allocated flag */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fermion_links.h"

/* Temporary work space for dslash_eo_field */ 
static su3_vector *temp ;
/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

void cleanup_dslash_temps(){
  if(!temp_not_allocated)
      free(temp) ; 
  temp_not_allocated=1 ;
}

void dslash_eo_site( field_offset src, field_offset dest, int parity,
		     eo_links_t *eo )
{
  register int i;
  register site *s;
  register int ipath,otherparity;
  register Real x; /* coefficient of path */
  ks_action_paths *ap = eo->ap;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;
  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Parallel transport by all the paths in the action.  
     Multiply by coefficient in table
  */
  FORSOMEPARITY(i,s,parity){ clearvec( (su3_vector *)F_PT(s,dest) ); }
  
  for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
    path_transport_site( src, F_OFFSET(tempvec[0]), parity,
		    q_paths[ipath].dir, q_paths[ipath].length );
    x=q_paths[ipath].coeff;
    FORSOMEPARITY(i,s,parity){
      scalar_mult_add_su3_vector(  (su3_vector *)F_PT(s,dest),
		   &(s->tempvec[0]), x, (su3_vector *)F_PT(s,dest) );
    }
  }   /* ipath */
} /* dslash_eo_site */

void dslash_eo_field( su3_vector *src, su3_vector *dest, int parity,
		      eo_links_t *eo )
{
  register int i;
  register site *s;
  register int ipath,otherparity;
  register Real x; /* coefficient of path */
  ks_action_paths *ap = eo->ap;
  int num_q_paths = ap->p.num_q_paths;
  Q_path *q_paths = ap->p.q_paths;
  
  /* allocate temporary work space only if not already allocated */
  if(temp_not_allocated)
    {
      temp =(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
      temp_not_allocated = 0 ;
    }
  
  switch(parity){
  case EVEN:	otherparity=ODD; break;
  case ODD:	otherparity=EVEN; break;
  case EVENANDODD:	otherparity=EVENANDODD; break;
  }
  
  /* Parallel transport by all the paths in the action.  
     Multiply by coefficient in table
  */
  FORSOMEPARITY(i,s,parity){ clearvec( &(dest[i]) ); }
  
  for( ipath=0; ipath<num_q_paths; ipath++ ){  /* loop over paths */
    path_transport_field( src, temp, parity,
		    q_paths[ipath].dir, q_paths[ipath].length );
    x=q_paths[ipath].coeff;
    FORSOMEPARITY(i,s,parity){
      scalar_mult_add_su3_vector(  &(dest[i]), &(temp[i]), x, &(dest[i]) );
    }
  }   /* ipath */
} /* dslash_eo_field */


