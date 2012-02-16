/***************** control.c *****************************************/

/* Main procedure for measuring the staggered fermion chiral
   condensate and related quantities */

/* MIMD version 7 */

/* The chiral condensate is computed on a supplied background field */

/* Modifications ... */
   
//  $Log: control.c,v $
//  Revision 1.2  2012/02/16 16:54:33  detar
//  Print hisq_svd_counter diagnostics only from node 0
//
//  Revision 1.1  2011/12/02 04:38:14  detar
//  Add
//
   

#define CONTROL
#include "ks_measure_includes.h"
#include <string.h>

int main(int argc, char *argv[])
{
  int prompt;
  int k;
  double starttime, endtime;
  
  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();
  /* set up */
  prompt = setup();
  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;

    starttime=dclock();
    
    total_iters=0;
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif
    
    /**************************************************************/
    /* Compute chiral condensate and related quantities           */

    /* Make valid fermion links if not already done */


    for(k = 0; k < param.num_set; k++){
      int num_pbp_masses = param.num_pbp_masses[k];
      int i0 = param.begin_pbp_masses[k];
      restore_fermion_links_from_site(fn_links, param.qic_pbp[i0].prec);

      if(num_pbp_masses == 1){
	f_meas_imp_field( param.npbp_reps[k], &param.qic_pbp[i0], param.ksp_pbp[i0].mass, 
			  param.ksp_pbp[i0].naik_term_epsilon_index, fn_links);
#ifdef D_CHEM_POT
	Deriv_O6_field(param.npbp_reps[k], &param.qic_pbp[i0], param.ksp_pbp[i0].mass,
		       fn_links, param.ksp_pbp[i0].naik_term_epsilon_index, 
		       param.ksp_pbp[i0].naik_term_epsilon);
#endif
      } else {
	f_meas_imp_multi( param.num_pbp_masses[k], param.npbp_reps[k], &param.qic_pbp[i0], 
			  &param.ksp_pbp[i0], fn_links);
#ifdef D_CHEM_POT
	Deriv_O6_multi( param.num_pbp_masses[k], param.npbp_reps[k], &param.qic_pbp[i0],
			&param.ksp_pbp[i0], fn_links);
#endif
      }
    }

    node0_printf("RUNNING COMPLETED\n");
    endtime=dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
      node0_printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
    fflush(stdout);

    /* save lattice if requested */
    if( param.saveflag != FORGET ){
      rephase( OFF );
      save_lattice( param.saveflag, param.savefile, param.stringLFN );
      rephase( ON );
    }

    destroy_ape_links_3D(ape_links);

    /* Destroy fermion links (created in readin() */

#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
  } /* readin(prompt) */

  return 0;
}
