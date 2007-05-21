/****** fermion_links_fn.c  -- ******************/
/* MIMD version 7 */
/* Link fattening routines for varions staggered actions
   CD 9/8/06 separated from quark_stuff.c 
   CD 10/15/06 Moved dm_du0 stuff to fermion_links_fn_dmdu0.c
*/

/* External entry points

   load_fn_links
   load_fn_links_dmdu0 (ifdef DM_DU0)
   invalidate_fn_links

 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static int valid_fn_links = 0;
static int valid_fn_links_dmdu0 = 0;

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
void load_fn_links(){

  if(valid_fn_links == 1)return;

  load_fatlinks(&t_fatlink, get_quark_path_coeff(), get_q_paths());
  load_longlinks(&t_longlink);

#ifdef DBLSTORE_FN
  load_fatbacklinks(&t_fatbacklink, t_fatlink);
  load_longbacklinks(&t_longbacklink, t_longlink);
#endif

  valid_fn_links = 1;
}

#ifdef DM_DU0
void load_fn_links_dmdu0(){
  if(valid_fn_links_dmdu0 == 1)return;

  load_fatlinks(&t_dfatlink_du0, get_quark_path_coeff_dmdu0(), 
		get_q_paths_dmdu0());
  valid_fn_links_dmdu0 = 1;
}
#endif

void
invalidate_fn_links( void )
{
  valid_fn_links = 0;
  valid_fn_links_dmdu0 = 0;
}


