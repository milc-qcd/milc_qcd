/******************** control.c *****************************/
/* MIMD version 7 */
/* Main procedure for SU3 eigenvalues with improved dynamical fermions */

#define CONTROL
#include "ks_eig_includes_qdp.h"	/* definitions files and prototypes */
#include "lattice_qdp.h"

EXTERN gauge_header start_lat_hdr;     /* Input gauge field header */

int
main(int argc, char *argv[])
{
  double chirality;
  double dtime;
  double *eigVal;
  su3_vector **eigVec;
  su3_vector *tmp;
  site *s;
  int i, si;
  int prompt;
  int total_R_iters=0;

  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#ifndef QDP_PROFILE
  QDP_profcontrol(0);
#endif
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();

  /* set up */
  prompt = setup();

  dtime = -dclock();
  /* loop over input sets */
  while( readin(prompt) == 0) {

    /* call fermion_variable measuring routines */
    /* results are printed in output file */
    /* NOTE: How do we initialize fat and long links in Dslash? */
    f_meas_imp(F_OFFSET(phi), F_OFFSET(xxx), mass);
    eigVal = malloc(Nvecs*sizeof(double));
    eigVec = malloc(Nvecs*sizeof(su3_vector*));
    for(i=0; i<Nvecs; i++)
      eigVec[i] = (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));

    total_R_iters += Kalkreuter(eigVec, eigVal, eigenval_tol,
				error_decr, Nvecs, MaxIter, Restart,
				Kiters, EVEN, &fn_links, &ks_act_paths);
    tmp = (su3_vector*)malloc(sites_on_node*sizeof(su3_vector));
    load_ferm_links(&fn_links, &ks_act_paths);
    for(i=0; i<Nvecs; i++) {
      /* Construct to odd part of the vector.                *
       * Note that the true odd part of the eigenvector is   *
       * i/sqrt(eigVal) Dslash Psi. But since I only compute *
       * the chirality the i factor is irrelevant (-i)*i=1!! */
      dslash_fn_field(eigVec[i], tmp, ODD, &fn_links);
      FORSOMEPARITY(si,s,ODD) {
	scalar_mult_su3_vector( &(tmp[si]),
				1.0/sqrt(eigVal[i]),
				&(eigVec[i][si]) );
      }

      measure_chirality(eigVec[i], &chirality, EVENANDODD);
      /* Here I divide by 2 since the EVEN vector is normalized            *
       * to 1. The EVENANDODD vector is normalized to 2. I could have      *
       * normalized the EVENANDODD vector to 1 and then not devide by to.  *
       * The measure_chirality routine assumes vectors normalized to 1.    */
      node0_printf("Chirality(%i): %g\n",i,chirality/2) ;
    }
    free(tmp);
    /**
       for(i=0;i<Nvecs;i++)
       {
       sprintf(label,"DENSITY(%i)",i) ;
       print_densities(eigVec[i], label, ny/2,nz/2,nt/2, EVEN) ;
       }
    **/
    for(i=0; i<Nvecs; i++) free(eigVec[i]);
    free(eigVec);
    free(eigVal);
#ifdef FN
    invalidate_all_ferm_links(&fn_links);
#endif
    fflush(stdout);

    node0_printf("RUNNING COMPLETED\n"); fflush(stdout);

    dtime += dclock();
    if(this_node==0) {
      printf("Time = %e seconds\n", dtime);
      printf("total_iters = %d\n", total_iters);
      printf("total Rayleigh iters = %d\n", total_R_iters);
    }
    fflush(stdout);
  }
  QDP_finalize();
  return 0;
}
