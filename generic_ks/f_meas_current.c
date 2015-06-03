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

/*Write current record for one random source */
static int
write_vector_current_record(QIO_Writer *outfile, int jrand, Real mass, Real *j_mu){
  int status;
  QIO_String *recxml = QIO_string_create();
  char recinfo[NRECINFO];
  snprintf(recinfo, NRECINFO, "source index %d mass %g", jrand, mass);
  QIO_string_set(recxml, recinfo);
  if(PRECISION == 1)
    status = write_F_R_from_field(outfile, recxml, j_mu, NMU);
  else
    status = write_D_R_from_field(outfile, recxml, j_mu, NMU);
  QIO_string_destroy(recxml);
  return status;
}

void 
f_meas_current( int nrand, quark_invert_control *qic, Real mass,
		int naik_term_epsilon_index, fermion_links_t *fl, 
		char *filename){

  imp_ferm_links_t* fn = get_fm_links(fl)[naik_term_epsilon_index];

  /* local variables for accumulators */
  register int i;

  int jrand;
  int mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  Real *j_mu = create_r_array_field(NMU);

  /* Open file for writing */
  QIO_Writer *outfile = open_vector_current_file(filename);
  if(outfile == NULL){
    node0_printf("f_meas_current: Failed to open %s\n", filename);
    exit(1);
  }
  
  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    /* Loop over directions for the current */
	/* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr, EVENANDODD );
#else
    z2rsource_plain_field( gr, EVENANDODD );
#endif
    
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
	  j_mu[NMU*i + mu] = cc.imag;
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

      wtime -= dclock();
      int status = write_vector_current_record(outfile, jrand, mass, j_mu);
      if(status != QIO_SUCCESS){
	node0_printf("f_meas_curent: Failed to write record to %s\n", filename);
      } else {
	node0_printf("f_meas_current: Wrote current density for source %d and mass %g on file %s\n", 
		     jrand, mass, filename);
      }
      wtime += dclock();

  } /* jrand */

  close_vector_current_file(outfile);
  node0_printf("Time to write %d records = %e\n", nrand, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_r_array_field(j_mu, NMU);
}

static Real *
create_real_array(int n){
  Real *a;
  int i;

  a = (Real *)malloc(n*sizeof(Real));
  if(a == NULL){
    printf("f_meas: No room for array\n");
    terminate(1);
  }
  for(i = 0; i < n; i++)a[i] = 0.;
  return a;
}

static void 
destroy_real_array(Real *a){
  if(a == NULL)return;
  free(a);
}

/* Entry point for multiple masses.  Saves a few cycles because one
   inversion can be done with the multimass inverter */

void 
f_meas_current_multi( int n_masses, int nrand, quark_invert_control *qic, 
		      ks_param *ksp, fermion_links_t *fl, 
		      char filenames[][MAXFILENAME]){
  
  Real *mass = create_real_array(n_masses);
  imp_ferm_links_t **fn = get_fm_links(fl);

  int i, j;
  int jrand;
  imp_ferm_links_t **fn_multi;
  int mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr[n_masses];
  Real *j_mu = create_r_array_field(NMU);

  QIO_Writer *outfile[n_masses];

  /* Create vector fields */
  for(j = 0; j < n_masses; j++){
    M_inv_gr[j] = create_v_field();
  }

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("f_meas_current_multi: Failed to open %s\n", filenames[j]);
      exit(1);
    }
  }

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;

  /* Load pointers for fermion links, based on Naik epsilon indices */
  fn_multi = (imp_ferm_links_t **)malloc(sizeof(imp_ferm_links_t *)*n_masses);
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];

  double wtime = 0.;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
      
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr, EVENANDODD );
#else
    z2rsource_plain_field( gr, EVENANDODD );
#endif

    /* M_inv_gr = M^{-1} gr */
    total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic, fn_multi );

    for(j = 0; j < n_masses; j++){

      /* Loop over directions for the current */
      for(mu = 0; mu < NMU; mu++){

	/* Apply the appropriate spin_taste operator for
	   a nearly conserved current.  */
	spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[j]);
	spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
	/* J_mu = imag[gr.M_inv_gr] */
	FORALLFIELDSITES(i){
	  complex cc = su3_dot( gr+i, gr_mu+i );
	  j_mu[NMU*i + mu] = cc.imag;
	}
      } /* mu */

      wtime -= dclock();
      int status = write_vector_current_record(outfile[j], jrand, mass[j], j_mu);
      if(status != QIO_SUCCESS){
	node0_printf("f_meas_current_multi: Failed to write record to %s\n", filenames[j]);
      } else {
	node0_printf("f_meas_current_multi: Wrote current density for source %d and mass %g on file %s\n", 
		     jrand, mass[j], filenames[j]);
      }
      wtime += dclock();

    } /* j */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }

  node0_printf("Time to write %d records for %d masses = %e\n", nrand, n_masses, wtime);

  destroy_r_array_field(j_mu, NMU); j_mu = NULL;
  
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  free(fn_multi);
  destroy_real_array(mass); mass = NULL;
}

