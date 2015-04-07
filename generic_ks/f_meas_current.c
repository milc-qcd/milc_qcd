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
  char *spin_taste_label[NMU] = {"GX-G1", "GY-G1", "GZ-G1", "GT-G1"};
  static int spin_taste[NMU];
  int mu;
  
  /* Decode spin-taste label */
  for(mu = 0; mu < NMU; mu++){
    char dummy[6];
    strncpy(dummy, spin_taste_label[mu], 6);
    spin_taste[mu] = spin_taste_index(dummy);
  }
  
  return spin_taste;
}

/* Open file for writing */
static QIO_Writer *
open_vector_current_file(char *filename){
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG disconnected vector current</title>";
  int volfmt = QIO_SINGLEFILE;
  int serpar = QIO_PARALLEL;
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
write_vector_current_record(QIO_Writer *outfile, int jrand, Real mass, complex *j_mu){
  int status;
  QIO_String *recxml = QIO_string_create();
  char recinfo[NRECINFO];
  snprintf(recinfo, NRECINFO, "source index %d mass %g", jrand, mass);
  QIO_string_set(recxml, recinfo);
  if(PRECISION == 1)
    status = write_F_C_from_field(outfile, recxml, j_mu, NMU);
  else
    status = write_D_C_from_field(outfile, recxml, j_mu, NMU);
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
  su3_vector *M_inv_gr_mu = create_v_field();
  complex *j_mu = create_c_array_field(NMU);

  /* Open file for writing */
  QIO_Writer *outfile = open_vector_current_file(filename);
  if(outfile == NULL){
    node0_printf("f_meas_current: Failed to open %s\n", filename);
    exit(1);
  }
  
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    /* Loop over directions for the current */
      for(mu = 0; mu < NMU; mu++){
	
	/* Make random source, and do inversion */
#ifndef Z2RSOURCE
	grsource_plain_field( gr, EVENANDODD );
#else
	z2rsource_plain_field( gr, EVENANDODD );
#endif
	/* Apply the appropriate spin_taste operator for
	   a local current */
	spin_taste_op(spin_taste[mu], r_offset, gr_mu, gr);
	
	/* M_inv_gr_mu = M^{-1} gr_mu */
	
	mat_invert_uml_field( gr_mu, M_inv_gr_mu, qic, mass, fn );
	
	/* J_mu = gr.M_inv_gr_mu */
	double norm = 1./(double)volume;
	FORALLFIELDSITES(i){
	  complex cc = su3_dot( gr+i, M_inv_gr_mu+i );
	  CMULREAL(cc, norm, j_mu[NMU*i + mu]);
	}
      } /* mu */

#if 0      
      /* DEBUG */
      FORALLFIELDSITES(i){
	printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	for(mu = 0; mu < NMU; mu++)
	  printf("(%g, %g) ",j_mu[NMU*i + mu].real,j_mu[NMU*i + mu].imag);
	printf("\n");
      }
#endif

      int status = write_vector_current_record(outfile, jrand, mass, j_mu);
      if(status != QIO_SUCCESS){
	node0_printf("f_meas_curent: Failed to write record to %s\n", filename);
      } else {
	node0_printf("f_meas_current: Wrote current density for source %d and mass %g on file %s\n", 
		     jrand, mass, filename);
      }

  } /* jrand */
  
  close_vector_current_file(outfile);
  destroy_v_field(M_inv_gr_mu); M_inv_gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_c_array_field(j_mu, NMU);
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
  su3_vector *M_inv_gr_mu[n_masses];
  complex *j_mu[n_masses];

  QIO_Writer *outfile[n_masses];

  /* Create vector fields */
  for(j = 0; j < n_masses; j++){
    M_inv_gr_mu[j] = create_v_field();
    j_mu[j] = create_c_array_field(NMU);
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

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    /* Loop over directions for the current */
    for(mu = 0; mu < NMU; mu++){
      
      /* Make random source, and do inversion */
#ifndef Z2RSOURCE
      grsource_plain_field( gr, EVENANDODD );
#else
      z2rsource_plain_field( gr, EVENANDODD );
#endif

      /* Apply the appropriate spin_taste operator for
	 a local current */
      spin_taste_op(spin_taste[mu], r_offset, gr_mu, gr);

      /* M_inv_gr_mu = M^{-1} gr_mu */

      total_iters += mat_invert_multi( gr_mu, M_inv_gr_mu, ksp, n_masses, qic, fn_multi );
      
      /* J_mu = gr.M_inv_gr_mu */
      double norm = 1./(double)volume;

      for(j = 0; j < n_masses; j++){
	
	/* psi-bar-psi on even sites = gr.M_inv_gr */
	FORALLFIELDSITES(i){
	  complex cc = su3_dot( gr+i, M_inv_gr_mu[j]+i );
	  CMULREAL(cc, norm, j_mu[j][NMU*i + mu]);
	}
      } /* j */
    } /* mu */
    
    for(j = 0; j < n_masses; j++){
      int status = write_vector_current_record(outfile[j], jrand, mass[j], j_mu[j]);
      if(status != QIO_SUCCESS){
	node0_printf("f_meas_current_multi: Failed to write record to %s\n", filenames[j]);
      } else {
	node0_printf("f_meas_current_multi: Wrote current density for source %d and mass %g on file %s\n", 
		     jrand, mass[j], filenames[j]);
      }
    }
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_c_array_field(j_mu[j], NMU); j_mu[j] = NULL;
    destroy_v_field(M_inv_gr_mu[j]); M_inv_gr_mu[j] = NULL;
  }
  
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  free(fn_multi);
  destroy_real_array(mass); mass = NULL;
}

