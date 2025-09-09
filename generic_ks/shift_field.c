/******************************** shift_field.c *************************/
/* MIMD version 7                                                       */

/* Carry out a symmetric shift or a plain forward or backward shift  
   of a color field with parallel transport defined by a set
   of provided links */

/*-------------------------------------------------------------*/
/* Apply the shift operator in direction "dir" with type fb    *
 * The KS phases MUST BE in the links                          */

#include "generic_ks_includes.h"
#include "../include/imp_ferm_links.h"

void 
shift_field_cpu(int dir, enum shift_dir fb, su3_vector *dest, const su3_vector *const src,
		const su3_matrix *const links)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2] = {NULL, NULL};
  su3_vector *tvec = create_v_field();

  //node0_printf("Using CPU shift\n");

  if(fb == SHIFT_FORWARD || fb == SHIFT_SYMMETRIC)
    tag[0] = start_gather_field( (void *)src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
  
  if(fb == SHIFT_BACKWARD || fb == SHIFT_SYMMETRIC){
    FORALLFIELDSITES(i)
    {
      mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec+i ) ;
    }
    tag[1] = start_gather_field((void *)tvec, sizeof(su3_vector), OPP_DIR(dir), 
				EVENANDODD, gen_pt[1] );
  }
  
  if(fb == SHIFT_FORWARD || fb == SHIFT_SYMMETRIC){
    wait_gather(tag[0]);
    FORALLSITES(i,s)
      {
	mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i );
      }
    cleanup_gather(tag[0]);

  } else {

    clear_v_field(dest);

  }

  if(fb == SHIFT_BACKWARD || fb == SHIFT_SYMMETRIC){
    wait_gather(tag[1]);
    FORALLSITES(i,s)
      {
	add_su3_vector( dest+i, (su3_vector *)gen_pt[1][i], dest+i ) ;    
      }
    cleanup_gather(tag[1]);
  }

  if(fb == SHIFT_SYMMETRIC){
    /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
    FORALLSITES(i,s)
      {
	scalar_mult_su3_vector( dest+i, 0.5, dest+i );
      }
  }

  destroy_v_field(tvec);
}

#if defined(HAVE_QUDA) && defined(USE_SHIFT_GPU)
#include <quda_milc_interface.h>
void 
shift_field(int dir, enum shift_dir fb, su3_vector *dest, const su3_vector *const src,
	    const su3_matrix *const links, int *refresh_links)
{
  int quda_precision = MILC_PRECISION;
  int sym = 0;

  node0_printf("Using GPU shift\n");

  switch(fb){
  case SHIFT_FORWARD:
    sym = 1;
    break;
  case SHIFT_BACKWARD:
    sym = 2;
    break;


  case SHIFT_SYMMETRIC:
    sym = 3;
    break;
  }

  qudaShift(MILC_PRECISION, quda_precision, links, src, dest, dir, sym, *refresh_links);
  *refresh_links = 0;
}

#else

/* CPU version */
void 
shift_field(int dir, enum shift_dir fb, su3_vector *dest, const su3_vector *const src,
	    const su3_matrix *const links, int *refresh_links)
{
  shift_field_cpu(dir, fb, dest, src, links);
}

#endif
