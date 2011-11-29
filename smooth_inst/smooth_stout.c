/********** smooth_stout.c ****************************************************/
/* MIMD version 7 */

/*
 Make stout links

 Here, the input "ape_weight", is the parameter rho of Morningstar-Peardon,
 ie. it does not have a factor 1/6 implied.
*/

#include "smooth_inst_includes.h"

void smooth(void)
{
   register int i,dir;
   register site *s;
   Real alpha, w_link, staple_weight;
   su3_matrix C, W;

   alpha = 6. * ape_weight;
   staple_weight = alpha / (6.*(1.0-alpha));
   w_link = (1.0-alpha);

   /* compute staple links */
   /* initial matrices in link.  staple links in fatlink */
   ape_smear(F_OFFSET(link[0]),F_OFFSET(fatlink[0]), 
	     staple_weight, 1., 0, 0, 0.);

   FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++)
   {
     /* subtract link, to get staple */
     scalar_mult_sub_su3_matrix( &(s->fatlink[dir]), &(s->link[dir]),
				 w_link, &C);
     /* do the stout smearing */
     stout_smear( &W, &C, &(s->link[dir]));
     /* copy result */
     su3mat_copy( &W, &(s->link[dir]) );
   }
}
