/***************** control.c *****************************************/

/* Main procedure for staggered fermion spectroscopy     	     */
/* MIMD version 7 */

/* This version computes propagators for staggered fermions on a supplied
   background field config and ties them together according to the
   input parameters. */

/* Modifications ... */
   
//  $Log: control.c,v $
//  Revision 1.8  2013/12/28 20:57:27  detar
//  Fix improper type for creating ksprop quark field
//
//  Revision 1.7  2013/12/24 05:39:00  detar
//  Add combo operation
//
//  Revision 1.6  2012/11/24 05:14:20  detar
//  Add support for U(1) fields and for future HYPISQ action
//
//  Revision 1.5  2012/05/08 20:39:54  detar
//  Call qudaFinalize to allow writing optimization file.
//
//  Revision 1.4  2012/04/25 03:23:21  detar
//  Fix rephase flag
//
//  Revision 1.3  2012/01/21 21:34:36  detar
//  Move start time to beginning.  Remake APE links after gauge fixing.
//
//  Revision 1.2  2011/12/03 03:43:39  detar
//  Cosmetic
//
//  Revision 1.1  2011/11/30 22:11:38  detar
//  Add
//
   

#define CONTROL
#include "ks_deflation_includes.h"
#include <string.h>
#ifdef HAVE_QUDA
#include <quda_milc_interface.h>
#endif


int main(int argc, char *argv[])
{
  const int m = 144;//max search dimension
  const int k = 8; //Nev

  printf("\nEigCG test, m = %d, k = %d, precision %d.\n", m, k, sizeof(Real) );

  int prompt;
  int i, j, s;

  int naik_index, naik_index0, naik_index1;
  double mass;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif
  
  initialize_machine(&argc,&argv);

  if(this_node==0)printf("\nEnd initializing.\n");

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  if(this_node==0)printf("\nPassed.\n");  

  g_sync();

  starttime=dclock();
    
  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("setup");

  /* loop over input sets */

  while( readin(prompt) == 0){
    
    if(prompt == 2)continue;

    total_iters=0;
#ifdef HISQ_SVD_COUNTER
    hisq_svd_counter = 0;
#endif
#ifdef HYPISQ_SVD_COUNTER
    hypisq_svd_counter = 0;
#endif

    if(this_node==0)printf("\nStarting %d mass loop\n", param.num_pbp_masses);
    
    /**************************************************************/
    /* Compute chiral condensate and related quantities           */

    /* Make fermion links if not already done */

    for(i = 0; i < 1 /* param.num_pbp_masses*/; i++){

        restore_fermion_links_from_site(fn_links, param.qic_pbp[i].prec);
        naik_index = param.ksp_pbp[i].naik_term_epsilon_index;
        mass = param.ksp_pbp[i].mass;

        printf("\nMass = %le\n", mass);
        //f_meas_imp_field( param.npbp_reps, &param.qic_pbp[i], mass, naik_index, fn_links);


        /*no gauge fixing for tests*/
    
        if(this_node==0)printf("END OF HEADER\n");

        /**************************************************************/

        /* Loop over rhs (random sources here) */

        const int deflation_grid = param.Nvecs / k;

        const int num_rhs = param.nrhs;
        const int num_evs = param.Nvecs;

        if(this_node==0)printf("\nRHS number = %d, Evecs number = %d deflation grid = %d\n", num_rhs, num_evs, deflation_grid);

        su3_vector **dst, **rtz;

        Real *evals = (Real*)malloc(deflation_grid*num_evs*sizeof(Real));
        //
        imp_ferm_links_t *fn = NULL;

        fn = get_fm_links(fn_links)[naik_index];

        dst = (su3_vector **)malloc(num_rhs*sizeof(su3_vector *));

        for(j = 0; j < num_rhs; j++) dst[j] = create_v_field();

        rtz = (su3_vector **)malloc(num_evs*sizeof(su3_vector *));

        for(j = 0; j < num_evs; j++) rtz[j] = create_v_field();

        //param.qic.parity = EVENANDODD;

        STARTTIME;
        for(s = 0; s < num_rhs; s++){

           su3_vector *src = NULL;//current rhs

           su3_vector *aux = NULL; 

           src = create_v_field();

           aux = create_v_field();

           //generate random source
#ifndef Z2RSOURCE
           grsource_plain_field( src, EVENANDODD);
#else
           z2rsource_plain_field( src, EVENANDODD);
#endif
           
           ks_dirac_adj_op( src, aux, mass, EVENANDODD, fn);//defined in ../generic_ks/mat_invert.c

           total_iters = ks_eigcg_uml_field(s, num_rhs, src, dst[s], rtz, evals, &param.qic_pbp[i] , mass, m, k, deflation_grid, param.restart_tol, fn);

           destroy_v_field(aux);
           destroy_v_field(src);
      
        } /* rhs */
        ENDTIME("compute rhs");

        for(j = 0; j < num_rhs; j++) destroy_v_field(dst[j]);
        free(dst) ;

        for(j = 0; j < num_evs; j++) destroy_v_field(rtz[j]);
        free(rtz) ;

        free(evals); 
    
    }//end of mass loop
   
    node0_printf("RUNNING COMPLETED\n");
 
    endtime=dclock();
    
    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
#ifdef HISQ_SVD_COUNTER
    printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
#ifdef HYPISQ_SVD_COUNTER
    printf("hypisq_svd_counter = %d\n",hypisq_svd_counter);
#endif
    fflush(stdout);
    
    /* Destroy fermion links (created in readin() */
    
#if FERM_ACTION == HISQ
    destroy_fermion_links_hisq(fn_links);
#elif FERM_ACTION == HYPISQ
    destroy_fermion_links_hypisq(fn_links);
#else
    destroy_fermion_links(fn_links);
#endif
    fn_links = NULL;
  } /* readin(prompt) */
  

#ifdef HAVE_QUDA
  qudaFinalize();
#endif
  
  return 0;
}
