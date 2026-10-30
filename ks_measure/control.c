/***************** control.c *****************************************/

/* Main procedure for measuring the staggered fermion chiral
   condensate and related quantities */

/* MIMD version 7 */

/* The chiral condensate is computed on a supplied background field */

/* Modifications ... */

// CD 2/5/16 Merged Hiroshi Ohno's support for deflated and eigcg solves */
   
//  $Log: control.c,v $
//  Revision 1.3  2012/05/08 20:40:16  detar
//  Call qudaFinalize to allow writing optimization file
//
//  Revision 1.2  2012/02/16 16:54:33  detar
//  Print hisq_svd_counter diagnostics only from node 0
//
//  Revision 1.1  2011/12/02 04:38:14  detar
//  Add
//
   

#define CONTROL
#include "ks_measure_includes.h"
#include <string.h>
#ifdef HAVE_QUDA
#include "../include/generic_quda.h"
#include <quda_milc_interface.h>
#endif
#if defined(U1_VALENCE) || defined(U1_SEA)
#include "../include/generic_u1.h"
#include "../include/io_u1lat.h"
#endif

extern int gethostname (char *__name, size_t __len); // Should get this from unistd.h

/********************************************************************/
int main(int argc, char *argv[])
{
  int prompt;
  int k;
  double starttime, endtime;
#ifdef PRTIME
  double dtime;
#endif

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);
  
  g_sync();

  /* set up */
  STARTTIME;
  prompt = setup();
  ENDTIME("set up");

#if 0
  gethostname(hostname, 128);
  
#endif

  /* loop over input sets */

  while( 1 ){

    starttime = dclock();
    if(readin(prompt) != 0)break;
    
    /* Skip calculations if just proofreading input parameters */
    if(prompt == 2)continue;

    total_iters=0;


    /**************************************************************/
    /* Compute or reread eigenpairs */

    if(param.eigen_param.Nvecs > 0){

      /* Compute or reread eigenpairs if requested and check them */

      STARTTIME;

      restore_eigenvectors_from_file(get_fm_links(fn_links, 0), &param.eigen_param,
				     param.ks_eigen_startflag, param.ks_eigen_startfile, 
				     param.ks_eigen_saveflag, param.ks_eigen_savefile);

      ENDTIME("calculate/reload Dirac eigenpairs"); fflush(stdout);
    }

    /**************************************************************/
    /* Compute chiral condensate and other observables            */
    
    STARTTIME;
    
    for(k = 0; k < param.num_set; k++){

#if 0 /* If needed for debugging */
      /* (The Kalkreuter routine uses the random number generator to
	 initialize the eigenvector search, so, if you want to compare
	 results with and without deflation, you need to re-initialize
	 here.) */
      initialize_site_prn_from_seed(iseed);
#endif

      int num_pbp_masses = param.num_pbp_masses[k];
      Real masses[num_pbp_masses];
      Real charges[num_pbp_masses];
      char *eigstartfiles[num_pbp_masses];
      int eigstartflags[num_pbp_masses];
      imp_ferm_links_t *fn_mass[num_pbp_masses];
      int i0 = param.begin_pbp_masses[k];
      Real bdry_phase[4] = {0.,0.,0.,param.time_bc};
      
      /* Make table of FN links and masses and set boundary phases if
	 requested */
      imp_ferm_links_t *fn_pt[MAX_NAIK];
      for(int naik_index = 0; naik_index < MAX_NAIK; naik_index++)
	fn_pt[naik_index] = NULL;

#ifdef U1_SEA
      fn_mass[0] = get_fm_links(fn_links, 0);
      set_boundary_twist_fn(fn_mass[0], bdry_phase, param.coord_origin);
      /* Apply the operation if not already done */
      if(twist_status(fn_mass[0]) == OFF)boundary_twist_fn(fn_mass[0], ON);
#else
      for(int j = 0; j < num_pbp_masses; j++){
	int naik_index = param.ksp_pbp[i0+j].naik_term_epsilon_index;
	int charge_index = param.ksp_pbp[i0+j].charge_index;
	if(fn_pt[naik_index] == NULL)
	  fn_pt[naik_index] = get_fm_links(fn_links_charge[charge_index], naik_index);
	fn_mass[j] = fn_pt[naik_index];
	/* Set values in the structure fn */
	set_boundary_twist_fn(fn_mass[j], bdry_phase, param.coord_origin);
	/* Apply the operation if not already done */
	if(twist_status(fn_mass[j]) == OFF)boundary_twist_fn(fn_mass[j], ON);
      }
#endif
      for(int j = 0; j < num_pbp_masses; j++){
	eigstartfiles[j] = param.ks_eigen_startQfile[i0+j];
	eigstartflags[j] = param.ks_eigen_startQflag[i0+j];
	masses[j] = param.ksp_pbp[i0+j].mass;
	charges[j] = param.ksp_pbp[i0+j].charge;
      }
      

#ifdef CURRENT_DISC

#ifdef U1_SEA

      /* Read one U(1) field from a file or Generate 
	 an array of U(1) fields if requested */

      int nU1_A = param.nU1A_generate[k];

      if(param.start_u1flag != GENERATE){
	u1_A_array = create_u1_A_field_array(1);
	/* reload_u1_lattice reads into u1_A. The pointer assignment
	   makes it go into the u1_A_array[0] implicitly */
	u1_A = u1_A_array[0];
	reload_u1_lattice(param.start_u1flag, param.start_u1file);
	nU1_A = 0;
      } else {
	/* Define gauge and scheme */
	initialize_u1_gauge_field(U1_COULOMB, U1_QED_L);
	u1_A_array = create_u1_A_field_array(nU1_A);
	generate_u1_gauge_field_array(nU1_A, u1_A_array);
      }

      /* (Re)read a U(1) field if requested */
      if(param.start_u1flag != GENERATE){
	nU1_A = 0;
      }

      /* We support only zero Naik epsilon here, so fn_mass[0]  */
      if(param.truncate_diff[k])
	f_meas_sea_current_diff( num_pbp_masses, param.npbp_reps[k],
				 param.thinning[k], &param.qic_pbp[i0], &param.qic_pbp_sloppy[i0],
				 masses, charges, nU1_A, u1_A_array, fn_mass[0], 
				 &param.pbp_filenames[i0] );
      else
	f_meas_sea_current( num_pbp_masses, param.npbp_reps[k],
			    param.thinning[k], &param.qic_pbp[i0],
			    masses, charges, nU1_A, u1_A_array, fn_mass[0],
			    &param.pbp_filenames[i0] );
#else

      if(param.truncate_diff[k])
	f_meas_valence_current_diff( num_pbp_masses, param.npbp_reps[k],
				     param.thinning[k],
				     &param.qic_pbp[i0], &param.qic_pbp_sloppy[i0],
				     masses, charges,
				     eigstartflags, eigstartfiles,
				     &param.eigen_param,
				     fn_mass, u1_A,
				     &param.pbp_filenames[i0] );
      else
	f_meas_valence_current( num_pbp_masses, param.npbp_reps[k],
				param.thinning[k],
				&param.qic_pbp[i0], masses, charges,
				eigstartflags, eigstartfiles,
				&param.eigen_param,
				fn_mass, u1_A,
				&param.pbp_filenames[i0] );
#endif

#else
      // THESE NEED FIXING NOW
      if(num_pbp_masses == 1)
	f_meas_imp_field( param.npbp_reps[k], &param.qic_pbp[i0],
			  param.ksp_pbp[i0].mass, 
			  param.ksp_pbp[i0].naik_term_epsilon_index,
			  fn_links);
      else
	f_meas_imp_multi( param.num_pbp_masses[k], param.npbp_reps[k],
			  &param.qic_pbp[i0], &param.ksp_pbp[i0], fn_links);
      
#ifdef D_CHEM_POT
      if(num_pbp_masses == 1)
	Deriv_O6_field(param.npbp_reps[k], &param.qic_pbp[i0],
		       param.ksp_pbp[i0].mass,
		       fn_links, param.ksp_pbp[i0].naik_term_epsilon_index, 
		       param.ksp_pbp[i0].naik_term_epsilon);
      else
	Deriv_O6_multi( param.num_pbp_masses[k], param.npbp_reps[k],
			&param.qic_pbp[i0], &param.ksp_pbp[i0], fn_links);
#endif
#endif

      
      /* Clean up */
      for(int naik_index = 0; naik_index < MAX_NAIK; naik_index++)
	if(fn_pt[naik_index] != NULL)
	  destroy_fn_links(fn_pt[naik_index]);

#ifdef U1_SEA

      if(nU1_A == 0)
	destroy_u1_A_field_array(1, u1_A_array);
      else
	destroy_u1_A_field_array(nU1_A, u1_A_array);
      u1_A_array = NULL;
      u1_A = NULL;
	
#endif

    } /* k num_set */
 
#ifdef HISQ_SVD_COUNTER
    node0_printf("hisq_svd_counter = %d\n",hisq_svd_counter);
#endif
    
    ENDTIME("calculate observables");
    
    /* save lattice if requested */
    if( param.saveflag != FORGET ){
      STARTTIME;
      rephase( OFF );
      save_lattice( param.saveflag, param.savefile, param.stringLFN );
      rephase( ON );
      ENDTIME("save lattice");
    }

#ifdef U1_VALENCE
    if( param.save_u1flag != FORGET ){
      save_u1_lattice( param.save_u1flag, param.save_u1file );
    }
#endif

    /* Free host eigenvector storage if MILC allocated it.
       When QUDA owns the deflation space, these stayed NULL.
       free(NULL) is a no-op, so this is safe and self-consistent */

    destroy_local_eigenpairs();
    
    node0_printf("RUNNING COMPLETED\n");
    endtime = dclock();

    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));
    node0_printf("total_iters = %d\n",total_iters);
    fflush(stdout);

    destroy_ape_links_4D(ape_links);

    /* Destroy fermion links (created in readin() */


    for(int j = 0; j < n_charges; j++){
      if(fn_links_charge[j] != NULL){
#if FERM_ACTION == HISQ
	destroy_fermion_links_hisq(fn_links_charge[j]);
#else
	destroy_fermion_links(fn_links_charge[j]);
#endif
	fn_links_charge[j] = NULL;
      }
    }
  } /* readin(prompt) */

  free(fn_links_charge);
  free_lattice();

#ifdef HAVE_QUDA
  finalize_quda();
#endif

  normal_exit(0);
  return 0;
}
