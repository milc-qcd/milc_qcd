/******************* smearing.c ****************************************/

/* MIMD version 7 */
/* original version of 5/23/94 by UMH */
/* 2/19/98 Version 5 port CD */
/* 11/25/01 Removed redundant code. Changed to call ape_smear_dir CD */
/* See string_break/smearing.c and hvy_pot/smearing.c */

/* Perform one iteration of "APE smearing" on the space-like links.
   For the projection back to SU(3) use multiples of 3 hits over the
   3 different SU(2) subgroups until the procedure converges.
   Note: only this complete projection makes the smearing truly
   gauge invariant. */
#include <qdp.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"
#define TOL 1e-5
#define MAXCOUNT 50

void smearing( void ) {
  register int dir,i;
  register site *s;

  /* Loop over the space directions and APE-smear the links.
     The results will temporarily be stored in s_link, s_link_f
     and t_link_f for dir=XUP, YUP and ZUP and copied do link[dir]
     at the end. */
  
  for(dir=XUP;dir<=ZUP;dir++){

    /* Do the APE smearing and SU(3) projection for all links in dir */
    ape_smear_dir(F_OFFSET(link[0]),dir,F_OFFSET(diag),
		  1./smear_fac,1.,1,3*MAXCOUNT,TOL);
    /* Storage management: Temporarily store the new link */
    FORALLSITES(i,s){
      switch (dir){
      case XUP:
	su3mat_copy( &(s->diag), &(s->s_link) );
	break;
      case YUP:
	su3mat_copy( &(s->diag), &(s->s_link_f) );
	break;
      case ZUP:
	su3mat_copy( &(s->diag), &(s->t_link_f) );
	break;
      default:
	if(this_node == 0)printf("unknown direction %d\n", dir);
	break;
      }
      
    }
    
  } /* end loop over dir */

#ifdef CHECK_SMEAR_MILC_2  
  FORALLSITES(i,s){
    su3mat_copy( &(s->s_link), &(s->sm_link[XUP]) );
    su3mat_copy( &(s->s_link_f), &(s->sm_link[YUP]) );
    su3mat_copy( &(s->t_link_f), &(s->sm_link[ZUP]) );
  }
#else
  /* Now copy the smeared links, overwriting the previous ones */
  FORALLSITES(i,s){
    su3mat_copy( &(s->s_link), &(s->link[XUP]) );
    su3mat_copy( &(s->s_link_f), &(s->link[YUP]) );
    su3mat_copy( &(s->t_link_f), &(s->link[ZUP]) );
  }
#endif
  
} /* smearing */
