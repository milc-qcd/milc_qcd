/******************* smearing.c ****************************************/

/* MIMD version 7 */
/* original version of 5/23/94 by UMH */
/* 2/19/98 Version 5 port CD */
/* 11/25/01 Removed redundant code. Changed to call ape_smear_dir CD */

/* Perform one iteration of "APE smearing" on the space-like links.
   For the projection back to SU(3) use multiples of 3 hits over the
   3 different SU(2) subgroups until the procedure converges.
   Note: only this complete projection makes the smearing truly
   gauge invariant. */

#include "hvy_qpot_includes.h"
#define TOL 1e-5
#define MAXCOUNT 50

void smearing( void ) {
  register int dir,i;
  register site *s;
  su3_matrix *x_link, *y_link, *z_link;
  
  /* Allocate space for temporary storage of results */
  x_link = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(x_link == NULL){
    fprintf(stderr,"smearing: CAN'T MALLOC x_link\n");
    fflush(stderr);
    terminate(1);
  }

  /* Allocate space for temporary storage of results */
  y_link = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(y_link == NULL){
    fprintf(stderr,"smearing: CAN'T MALLOC y_link\n");
    fflush(stderr);
    terminate(1);
  }

  /* Allocate space for temporary storage of results */
  z_link = 
    (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(z_link == NULL){
    fprintf(stderr,"smearing: CAN'T MALLOC z_link\n");
    fflush(stderr);
    terminate(1);
  }

  /* Loop over the space directions and APE-smear the links.
     The results will temporarily be stored in x_link, y_link
     and z_link for dir=XUP, YUP and ZUP and copied do link[dir]
     at the end. */
  
  for(dir=XUP;dir<=ZUP;dir++){

    /* Do the APE smearing and SU(3) projection for all links in dir */
    ape_smear_dir(F_OFFSET(link[0]),dir,F_OFFSET(diag),
		  1./smear_fac,1.,1,3*MAXCOUNT,TOL);
    
    /* Storage management: Temporarily store the new link */
    FORALLSITES(i,s){
      switch (dir){
      case XUP:
	su3mat_copy( &(s->diag), x_link+i );
	break;
      case YUP:
	su3mat_copy( &(s->diag), y_link+i );
	break;
      case ZUP:
	su3mat_copy( &(s->diag), z_link+i );
	break;
      default:
	if(this_node == 0)printf("unknown direction %d\n", dir);
	break;
      }
      
    }
    
  } /* end loop over dir */
  
  /* Now copy the smeared links, overwriting the previous ones */
  FORALLSITES(i,s){
    su3mat_copy( x_link+i, &(s->link[XUP]) );
    su3mat_copy( y_link+i, &(s->link[YUP]) );
    su3mat_copy( z_link+i, &(s->link[ZUP]) );
  }

  free(x_link);
  free(y_link);
  free(z_link);
  
} /* smearing */
