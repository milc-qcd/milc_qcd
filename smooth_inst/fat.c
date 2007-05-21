/********** fat.c ****************************************************/
/* MIMD version 7 */
/* Based on update.c */

/*
 Make fat links only
*/

#include "smooth_inst_includes.h"

void smooth(void)
{
   register int i,dir;
   register site *s;

   /* compute fat links */
   /* initial matrices in link.  fat links in fatlink */
   ape_smear(F_OFFSET(link[0]),F_OFFSET(fatlink[0]), 
	     ape_weight/(6*(1.0-ape_weight)), 1., 0, 0, 0.);

   /* copy result */
   FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++)
   {
     su3mat_copy( &(s->fatlink[dir]), &(s->link[dir]) );
   }
}
