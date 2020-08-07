/*************************** check_ani.c ************************/
/* MIMD version 7 */
/* Check implementation of anisotropy
   a) of anisotropic factors in fermion links 
      with point source on a unit configuration 
   b) pure gauge HMC updating with a single step 
      using a funny link setting
   c) full QCD HMC fermion force on a unit gauge 
      configuration with one special link 
      (still to be done)
*/

/* JHW 05/08/2020 */

/* Usage ...

   check_ani
   */

#define CONTROL

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "ani_non_prod_tests_includes.h"

#ifdef HAVE_QIO
#include <qio.h>
#endif

#define MAXERRCOUNT 10
#define TOLERANCE (0.0001)

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/


int main(int argc, char *argv[])
{

  int prompt, status;
  int todo;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  char myname="check_ani";

  ks_prop_field *source[MAX_SOURCE];
  ks_prop_field *prop[MAX_PROP];
  
  initialize_machine(&argc,&argv);

  for(int i = 0; i < MAX_SOURCE; i++)source[i] = NULL;
  for(int i = 0; i < MAX_PROP; i++)prop[i] = NULL;

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  this_node = mynode();
  number_of_nodes = numnodes();


  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  if( readin(prompt) == 0){

#ifdef FREE_KS_ANI_TEST
    node0_printf("Checking the anisotropic KS Dirac operator in free field theory \n");
    quark_source *qs = &param.src_qs[0];
    int mytype = qs->type;
    source[0] = create_ksp_field(qs->ncolor);;
    for(int color = 0; color < qs->ncolor; color++){
      if(v_source_field(source[0]->v[color], qs)){
        printf("%s(%d): error getting source\n",myname,this_node);
        terminate(1);
      }
    } 
    qs->type = mytype;
    prop[0] = create_ksp_field(qs->ncolor);

    for(int color = 0; color < qs->ncolor; color++){
      free_KS_ani_test(&(source[0]->v[color]),&(param.src_qs[0]),&(param.ksp[0]),1,&(prop[0]->v[color]));
    }

    /* Now destroy all sourcess */
    for(int i = 0; i < param.num_base_source+param.num_modified_source; i++){
      if(source[i] != NULL){
        node0_printf("destroy source[%d]\n",i);
        destroy_ksp_field(source[i]);
        source[i] = NULL;
      }
    }

    /* Now destroy all propagators */
    for(int i = 0; i <= param.end_prop[param.num_set-1]; i++){
      if(prop[i] != NULL){
        node0_printf("destroy prop[%d]\n",i);
        destroy_ksp_field(prop[i]);
        prop[i] = NULL;
      }
    }
#endif

    /* Destroy fermion links (created in readin() */
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
    destroy_fermion_links_hypisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
    starttime = endtime;

#ifdef FUNNYLINKS
    rephase(OFF);
    node0_printf("Checking the anisotropic Symanzik gauge action with funny link setup\n");
    funnylinks( param.umu );
    for(todo=param.trajecs; todo > 0; --todo ){  
      /* do the trajectories */
      update();
//#ifdef G_MEASURE
      /* do "extensive local" measurements */
      g_measure();
//#endif
    }

#endif
  }

#ifdef HAVE_QUDA
  qudaFinalize();
#endif

#ifdef HAVE_QPHIX
  finalize_qphix();
#endif

#ifdef HAVE_GRID
  finalize_grid();
#endif

  normal_exit(0);

  return 0;
}
