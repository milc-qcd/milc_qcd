/****** fermion_links_fn.c  -- ******************/
/* MIMD version 7 */
/* Link fattening routines for varions staggered actions
   CD 9/8/06 separated from quark_stuff.c 
   CD 10/15/06 Moved dm_du0 stuff to fermion_links_fn_dmdu0.c
*/

/* External entry points

   init_ferm_links
   load_ferm_links
   load_ferm_links_dmdu0 (ifdef DM_DU0)
   invalidate_all_ferm_links

 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

/********************************************************************/
/* Sum over paths connecting to nearest neighbor point (fat link) and to third
   nearest neighbor (longlinks) */
/********************************************************************/
/* Doug Toussaint 2/4/98 */
/* modified to use t_longlinks, S. Gottlieb 7/13/01 */
/* long link calculating routine */
/* path_product() follows the path starting at step 0, and
   leaves the answer at the end of the path.  We want the answer
   at the site where the path begins.  So we look for paths with
   the opposite displacement from the displacement of the point
   that we want to transport to this site, and take the adjoint
   of the matrix at the end. clear? */
/* KS phases and APBC must be in the links. See long comment at 
   end of fermion_force_general.c */

/* Load fat links into t_fatlink, t_longlink and t_fatbacklink,
   t_longbacklink */

void load_ferm_links(ferm_links_t *fn, ks_action_paths *ap){

  if(fn->valid == 1)return;

#ifdef FN
  load_fatlinks(fn, ap);
  load_longlinks(fn, ap);
#endif

  fn->ap = ap;

#ifdef DBLSTORE_FN
  load_fatbacklinks(fn);
  load_longbacklinks(fn);
#endif

  fn->valid = 1;
}

#ifdef DM_DU0
void load_ferm_links_dmdu0(ferm_links_t *fn, ks_action_paths *ap){
  if(fn->valid == 1)return;

#ifdef FN
  load_fatlinks(fn, ap);
#endif

  fn->ap = ap;
  fn->valid = 1;
}
#endif

void
invalidate_all_ferm_links(ferm_links_t *fn)
{
  fn->valid = 0;
}

/* For compatibility with fermion_links_hisq.c */
void
invalidate_fn_links(ferm_links_t *fn)
{
  invalidate_all_ferm_links(fn);
}

void 
init_ferm_links(ferm_links_t *fn){
  fn->valid = 0;
  fn->fat = NULL;
  fn->lng = NULL;
  fn->fatback = NULL;
  fn->lngback = NULL;
  fn->ap = NULL;
}
