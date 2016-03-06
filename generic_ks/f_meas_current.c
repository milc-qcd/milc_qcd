/**************** f_meas_current.c ***************************************/
/* MIMD version 7 */
/* CD 1/15 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure the fermionic observable:

    J_\mu
    
    Write the result to a file as a real four-vector field 

    Entry points

    f_meas_current
    f_meas_current_multi
*/

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fn_links.h"
#include "../include/io_scidac.h"
#include "../include/imp_ferm_links.h"
#include <qio.h>
#include <string.h>

#define NMU 4
#define NRECINFO 128

/* Create spin-taste indices for current */
static int *
get_spin_taste(void){
  
  /* Current spin-taste list */
  char *spin_taste_list[NMU] = {"rhoxsfape", "rhoysfape", "rhozsfape", "rhotsfape"};
  static int spin_taste[NMU];
  int mu;
  
  /* Decode spin-taste label */
  for(mu = 0; mu < NMU; mu++){
    char dummy[32];
    strncpy(dummy, spin_taste_list[mu], 32);
    spin_taste[mu] = spin_taste_index(dummy);
  }
  
  return spin_taste;
}

/* Open file for writing */
static QIO_Writer *
open_vector_current_file(char *filename){
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG disconnected vector current</title>";
  int volfmt = QIO_SINGLEFILE;
  int serpar = QIO_SERIAL;
  QIO_String *filexml = QIO_string_create();
  QIO_string_set(filexml, default_file_xml);
  QIO_Layout layout;
  QIO_Filesystem fs;
  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  QIO_Writer *outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO, NULL,
					   &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  return outfile;
}

/* Close file */
static void
close_vector_current_file(QIO_Writer *outfile){
  QIO_close_write(outfile);
}

/* Thin the random source */
static void
thin_random_source(su3_vector *src, int thinning, int ex, int ey, int ez, int et){
  site *s;
  int i;

  FORALLSITES(i,s) {
    if(s->x % thinning != ex || s->y % thinning != ey || s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
  }
}

/*Write current record for the accumulated average over random sources */
static int
write_vector_current_record(QIO_Writer *outfile, int jrand, int nwrite, Real mass, Real *j_mu){
  int status = QIO_SUCCESS;
  QIO_String *recxml = QIO_string_create();
  char recinfo[NRECINFO];
  Real *j_mu_avg;
  int i, mu;

  /* Compute average */
  j_mu_avg = create_r_array_field(NMU);
  FORALLFIELDSITES(i){
    for(mu = 0; mu < NMU; mu++)
      j_mu_avg[NMU*i + mu] = j_mu[NMU*i + mu]/nwrite;
  }

  snprintf(recinfo, NRECINFO, "source index %d mass %g", jrand, mass);
  QIO_string_set(recxml, recinfo);
  if(PRECISION == 1)
    status = write_F_R_from_field(outfile, recxml, j_mu_avg, NMU);
  else
    status = write_D_R_from_field(outfile, recxml, j_mu_avg, NMU);
  QIO_string_destroy(recxml);

  destroy_r_array_field(j_mu_avg, NMU);

  node0_printf("Wrote current density for source %d and mass %g\n", jrand, mass);
  
  return status;
}

void 
f_meas_current_diff( int nrand, int nwrite, int thinning, 
		     quark_invert_control *qic_precise, quark_invert_control *qic_sloppy,
		     Real mass, int naik_term_epsilon_index, fermion_links_t *fl, 
		     char *filename){
  
  char myname[] = "f_meas_current";

  imp_ferm_links_t* fn = get_fm_links(fl)[naik_term_epsilon_index];

  /* local variables for accumulators */
  int i, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  Real *j_mu = create_r_array_field(NMU);
  int ex, ey, ez, et, d = thinning;

  /* Open file for writing */
  QIO_Writer *outfile = open_vector_current_file(filename);
  if(outfile == NULL){
    node0_printf("f_meas_current: Failed to open %s\n", filename);
    exit(1);
  }
  
  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){

	/* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );
	    
	    /* First, sloppy solution */
	    /* M_inv_gr = M^{-1} gr */
	    mat_invert_uml_field( gr, M_inv_gr, qic_sloppy, mass, fn );
	    
	    /* Loop over directions for the current */
	    for(mu = 0; mu < NMU; mu++){
	      
	      /* Apply the appropriate spin_taste operator for
		 a nearly conserved current.  */
	      spin_taste_op_fn(fn, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	      spin_taste_op_fn(fn, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	      
	      /* J_mu = imag[gr.M_inv_gr] */
	      /* SUBTRACT the sloppy result */
	      FORALLFIELDSITES(i){
		complex cc = su3_dot( gr+i, gr_mu+i );
		j_mu[NMU*i + mu] -= cc.imag;
	      }
	    } /* mu */
	    
	    /* Next, continue to a "precise" solution from the same source */
	    /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	    mat_invert_uml_field( gr, M_inv_gr, qic_precise, mass, fn );
	    
	    /* Loop over directions for the current */
	    for(mu = 0; mu < NMU; mu++){
	      
	      /* Apply the appropriate spin_taste operator for
		 a nearly conserved current.  */
	      spin_taste_op_fn(fn, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	      spin_taste_op_fn(fn, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	      
	      /* J_mu = imag[gr.M_inv_gr] */
	      /* ADD the precise result, which then gives the difference */
	      FORALLFIELDSITES(i){
		complex cc = su3_dot( gr+i, gr_mu+i );
		j_mu[NMU*i + mu] += cc.imag;
	      }
	    } /* mu */
	    
#if 0      
	    /* DEBUG */
	    FORALLFIELDSITES(i){
	      printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	      for(mu = 0; mu < NMU; mu++)
		printf("%g ",j_mu[NMU*i + mu]);
	      printf("\n");
	    }
#endif
	    
	  } /* ex, ey, ez, et */
    
    /* Write at intervals of nwrite random values */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      int status = write_vector_current_record(outfile, jrand, nwrite, mass, j_mu);
      clear_r_array_field(j_mu, NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filename);
      }
    }
  } /* jrand */

  close_vector_current_file(outfile);
  node0_printf("Time to write %d records = %e\n", nrand, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
  destroy_r_array_field(j_mu, NMU);
}

void 
f_meas_current( int nrand, int nwrite, int thinning, quark_invert_control *qic, 
		Real mass, int naik_term_epsilon_index, fermion_links_t *fl, 
		char *filename){

  char myname[] = "f_meas_current";

  imp_ferm_links_t* fn = get_fm_links(fl)[naik_term_epsilon_index];

  /* local variables for accumulators */
  int i, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  Real *j_mu = create_r_array_field(NMU);
  int ex, ey, ez, et, d = thinning;

  /* Open file for writing */
  QIO_Writer *outfile = open_vector_current_file(filename);
  if(outfile == NULL){
    node0_printf("f_meas_current: Failed to open %s\n", filename);
    exit(1);
  }
  
  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){

	/* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif

    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );
	    
	    /* M_inv_gr = M^{-1} gr */
	    mat_invert_uml_field( gr, M_inv_gr, qic, mass, fn );
	    
	    /* Loop over directions for the current */
	    for(mu = 0; mu < NMU; mu++){
	      
	      /* Apply the appropriate spin_taste operator for
		 a nearly conserved current.  */
	      spin_taste_op_fn(fn, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	      spin_taste_op_fn(fn, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	      
	      /* J_mu = imag[gr.M_inv_gr] */
	      FORALLFIELDSITES(i){
		complex cc = su3_dot( gr+i, gr_mu+i );
		j_mu[NMU*i + mu] += cc.imag;
	      }
	    } /* mu */
	    
#if 0      
	    /* DEBUG */
	    FORALLFIELDSITES(i){
	      printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	      for(mu = 0; mu < NMU; mu++)
		printf("%g ",j_mu[NMU*i + mu]);
	      printf("\n");
	    }
#endif

	  } /* ex, ey, ez, et */
    
    /* Write at intervals of nwrite random values */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      int status = write_vector_current_record(outfile, jrand, nwrite, mass, j_mu);
      clear_r_array_field(j_mu, NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filename);
      }
    } /* if write */
    
  } /* jrand */

  close_vector_current_file(outfile);
  node0_printf("Time to write %d records = %e\n", nrand, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
  destroy_r_array_field(j_mu, NMU);
}

/*#if EIGMODE == EIGCG || EIGMODE == DEFLATE*/
#ifdef EIGMODE

/* Entry point for multiple masses with iterated single-mass inverter.
   This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.
   Designed for use with deflation or eigcg */

void 
f_meas_current_multi_diff( int n_masses, int nrand, int nwrite, int thinning,
			   quark_invert_control *qic_precise,
			   quark_invert_control *qic_sloppy, 
			   ks_param *ksp, fermion_links_t *fl, 
			   char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_diff";

  int i, j, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
 

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;
  
  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  
  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );
	    
	    for(j = 0; j < n_masses; j++){
	      
	      /* First, sloppy solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      mat_invert_uml_field( gr, M_inv_gr, qic_sloppy, mass[j], fn_multi[j]);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		/* SUBTRACT the sloppy result */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] -= cc.imag;
		}
	      } /* mu */
	      
	      /* Next, continue to a "precise" solution from the same source */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      mat_invert_uml_field( gr, M_inv_gr, qic_precise, mass[j], fn_multi[j]);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		/* ADD the precise result, which then gives the difference */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
	      
#if 0      
	      /* DEBUG */
	      FORALLFIELDSITES(i){
		printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
		for(mu = 0; mu < NMU; mu++)
		  printf("%g ",j_mu[j][NMU*i + mu]);
		printf("\n");
	      }
#endif
	    } /* j */
	  } /* ex, ey, ez, et */
	    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	} 
      }
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}


/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg */

void 
f_meas_current_multi( int n_masses, int nrand, int nwrite, int thinning,
		      quark_invert_control *qic, ks_param *ksp,
		      fermion_links_t *fl, 
		      char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi";

  int i, j, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
 

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;
  
  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  
  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );

	    for(j = 0; j < n_masses; j++){
	      
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      mat_invert_uml_field( gr, M_inv_gr, qic, mass[j], fn_multi[j]);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
	      
#if 0      
	      /* DEBUG */
	      FORALLFIELDSITES(i){
		printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
		for(mu = 0; mu < NMU; mu++)
		  printf("%g ",j_mu[j][NMU*i + mu]);
		printf("\n");
	      }
#endif
	    } /* j */
	  } /* ex, ey, ez, et */
      
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	} 
	clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
      } /* j */
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

#else

    /* Entry point for multiple masses.  Uses the multimass inverter */

void 
f_meas_current_multi_diff( int n_masses, int nrand, int nwrite, int thinning,
			   quark_invert_control *qic_precise,
			   quark_invert_control *qic_sloppy, 
			   ks_param *ksp, fermion_links_t *fl, 
			   char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_diff";

  int i, j, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr[n_masses];
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

  /* Create vector fields */
  for(j = 0; j < n_masses; j++){
    M_inv_gr[j] = create_v_field();
  }

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;

  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];

  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
      
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );
	    
	    /* First, sloppy solution */
	    /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	    total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic_sloppy, fn_multi );
	    
	    for(j = 0; j < n_masses; j++){
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[j]);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		/* SUBTRACT the sloppy result */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] -= cc.imag;
		}
	      } /* mu */
	    }
	    
	    /* Next, get a "precise" solution from the same source */
	    /* This won't be a continuation of the sloppy solution, because multimass always starts again from 0 */
	    /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	    total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic_precise, fn_multi );
	    
	    for(j = 0; j < n_masses; j++){
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[j]);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		/* ADD the precise result, which then gives the difference */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
	    } /* j */
	  } /* ex, ey, ez, et */
    
	    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	} 
	clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
      } /* j */
    } /* if write */
  } /* jrand */
    
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }
  
  node0_printf("Time to write %d records for %d masses = %e\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

/* Entry point for multiple masses.  Uses the multimass inverter */

void 
f_meas_current_multi( int n_masses, int nrand, int nwrite, int thinning,
		      quark_invert_control *qic, ks_param *ksp, 
		      fermion_links_t *fl, 
		      char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi";

  int i, j, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr[n_masses];
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

  /* Create vector fields */
  for(j = 0; j < n_masses; j++){
    M_inv_gr[j] = create_v_field();
  }

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;

  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];

  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
      
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){

	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_random_source( gr, d, ex, ey, ez, et );
	    
	    /* M_inv_gr = M^{-1} gr */
	    total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic, fn_multi );
	    
	    for(j = 0; j < n_masses; j++){
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[j]);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = imag[gr.M_inv_gr] */
		FORALLFIELDSITES(i){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
	    } /* j */
	  } /* ex, ey, ez, et */
	      
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	}
	clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
      } /* j */
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }

  node0_printf("Time to write %d records for %d masses = %e\n", nrand/nwrite, n_masses, wtime);

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

#endif

