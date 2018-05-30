/**************** f_meas_current.c ***************************************/
/* MIMD version 7 */
/* CD 1/15 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure the fermionic observable:

    J_mu
    
    Write the result to a file as a real four-vector field 

    Entry points

    f_meas_current_diff
    f_meas_current

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
thin_source(su3_vector *src, int thinning, int ex, int ey, int ez, int et){
  site *s;
  int i;

  FORALLSITES(i,s) {
    if(s->x % thinning != ex || s->y % thinning != ey ||
       s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
  }
}

/* Average the vector current density */
static void
average_vector_current(int nwrite, Real *j_mu){
  int i, mu;
  /* Compute average in place */
  FORALLFIELDSITES(i){
    for(mu = 0; mu < NMU; mu++)
      j_mu[NMU*i + mu] = j_mu[NMU*i + mu]/nwrite;
  }
}

/* Average the vector current density */
static void
average_vector_current_and_sum(int nwrite, Real *j_mu, Real *jadd_mu){
  int i, mu;
  /* Compute average in place */
  FORALLFIELDSITES(i){
    for(mu = 0; mu < NMU; mu++){
      j_mu[NMU*i + mu] = j_mu[NMU*i + mu]/nwrite + jadd_mu[NMU*i + mu];
      //      printf("j_mu  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, j_mu[NMU*i+mu]);
    }
  }
}

/*Write current record for the accumulated average over random sources */
static int
write_vector_current_record(QIO_Writer *outfile, int jrand, int nwrite, Real mass, Real *j_mu){
  int status = QIO_SUCCESS;
  QIO_String *recxml = QIO_string_create();
  char recinfo[NRECINFO];

  snprintf(recinfo, NRECINFO, "source index %d mass %g", jrand/nwrite, mass);
  QIO_string_set(recxml, recinfo);
  if(MILC_PRECISION == 1)
    status = write_F_R_from_field(outfile, recxml, j_mu, NMU);
  else
    status = write_D_R_from_field(outfile, recxml, j_mu, NMU);
  QIO_string_destroy(recxml);

  node0_printf("Wrote current density for source %d and mass %g\n", jrand, mass);
  
  return status;
}

/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
		   double_complex *dot, int parity) {
  register double re,im ;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORSOMEFIELDPARITY(i,parity){
    cc = su3_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  }
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
			  su3_vector *vec2, int parity){

  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEFIELDPARITY(i,parity){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  }
}

/************************************************************************/
/*  Projects out the *vectors from the  vec. Num is the Number of vectors  *
 * and parity is the parity on which we work on.                           *
 * The vectors are assumed to be orthonormal.                              */
static void project_out(su3_vector *vec, su3_vector *vector[], int Num, int parity){
  register int i ;
  double_complex cc ;
  double ptime = -dclock();

  if(Num == 0)return;
  
  for(i=Num-1;i>-1;i--){
    dot_product(vector[i], vec, &cc, parity) ;
    complex_vec_mult_sub(&cc, vector[i], vec, parity);
  }

  ptime += dclock();
  node0_printf("Time to project out low modes %g sec\n", ptime);
}

/************************************************************************/
static void
collect_evenodd_sources(su3_vector *gr[], int ns, int parity, int thinning,
			su3_vector *gr0){
  /* Create thinned sources of the specified parity */
  /* Result in gr */
  
  /* Iterate over displacements within a d^4 cube for this parity. */
  int ex, ey, ez, et;
  int is = 0;
  int d = thinning;
  for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
      if ( ((ex+ey+ez+et)%2==0 && parity == EVEN) ||
	   ((ex+ey+ez+et)%2==1 && parity == ODD) ){

	node0_printf("Source %d is %d %d %d %d\n", is, ex, ey, ez, et);

	/* Apply source thinning */
	copy_v_field( gr[is], gr0 );
	thin_source( gr[is], d, ex, ey, ez, et );
	
	/* Project out the low mode part, based on the given eigenvectors */
	int Nvecs = param.eigen_param.Nvecs;
	if(Nvecs > 0)
	  project_out(gr[is], eigVec, Nvecs, parity);
#if 0
	/* DEBUG */
	/* Check the norm of the reduced source */
	double_complex dd;
	dot_product(gr[is], gr[is], &dd, parity);
	node0_printf("Deflated source norm %g\n", dd.real);
#endif
	is++;
	if(is > ns){
	  node0_printf("collect_evenodd_sources: Internal error: too many sources\n");
	  terminate(1);
	}
      } /* ex, ey, ez, et */
}

/************************************************************************/
/* Collect diluted random sources for mrhs */

static void
collect_sources(su3_vector *gr_even[], su3_vector *gr_odd[], int nr, int thinning, int evol){

  su3_vector *gr0 = create_v_field();

  for(int jr = 0; jr < nr; jr++){
    
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    node0_printf("EVEN sources\n");
    collect_evenodd_sources(gr_even + jr*evol, nr*evol, EVEN, thinning, gr0);
    node0_printf("ODD sources\n");
    collect_evenodd_sources(gr_odd  + jr*evol, nr*evol, ODD,  thinning, gr0);
  } /* jr */

  destroy_v_field(gr0);
}

/************************************************************************/
/* Calculate current densities for a given quark mass and source
   parity using a list of thinned stochastic estimators */
/* nsrc sources in gr.  Results in j_mu */

static void
block_current_stochastic( Real *j_mu_mass, Real mass, int nsrc, int sign, int parity, 
			  quark_invert_control *qic, imp_ferm_links_t *fn_mass,
			  su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic";

  int otherparity;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  su3_vector **M_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    M_inv_gr[is] = create_v_field();
  }

  switch(parity){
  case(EVEN): otherparity=ODD; break;
  case(ODD):  otherparity=EVEN; break;
  }

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  qic->parity = parity;
  ks_congrad_block_field(nsrc, gr, M_inv_gr, qic, mass, fn_mass);

  /* Might be better to use a block dslash here?? */
  for(int is = 0; is < nsrc; is++)
    dslash_fn_field( M_inv_gr[is], M_inv_gr[is], otherparity, fn_mass);
  
  /* For each source, apply current in various directions at the sink */
  su3_vector *gr_mu = create_v_field();
  for(int is = 0; is < nsrc; is++)
    for(int mu = 0; mu < NMU; mu++){
      
      /* Apply the appropriate spin_taste operator for
	 a nearly conserved current.  */
      spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
      spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      /* J_mu = -imag[gr * M_inv_gr] */
      /* If sign = +1, add the result.  If sign = -1, subtract it */
      int i;
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_mu+i );
	j_mu_mass[NMU*i + mu] += -sign*cc.imag;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu_mass[NMU*i + mu]);
#endif
      }
    } /* is, mu */
  destroy_v_field(gr_mu);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M_inv_gr[is]);

}
      
/************************************************************************/

static void
block_current_diff(Real *j_mu[], int nwrite, int thinning, int n_masses,
		   quark_invert_control *qic_precise, quark_invert_control *qic_sloppy,
		   ks_param *ksp, fermion_links_t *fl){

  char myname[] = "block_current_diff";
  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  /* Block solver parameters -- temporary */
  int nr = 2;  /* Number of random sources to block */
  if(nwrite < nr)nr = nwrite;
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;

  su3_vector **gr_even = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **gr_odd = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    gr_even[is] = create_v_field();
    gr_odd[is] = create_v_field();
  }

  /* Loop over random sources in groups of nr */
  su3_vector *gr_mu = create_v_field();

  for(int jrand = 0; jrand < nwrite; jrand += nr){
    
    /* Create sources in gr_even and gr_odd */
    collect_sources(gr_even, gr_odd, nr, d, evol);

    /* Construct current density from the list of sources */
    imp_ferm_links_t **fn = get_fm_links(fl);
    for(int j = 0; j < n_masses; j++){
      Real mass = ksp[j].mass;
      imp_ferm_links_t *fn_mass = fn[j];
      
      /* First, the sloppy high-mode solution */
      node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", mass);
      
      block_current_stochastic( j_mu[j], mass, nsrc, -1, EVEN, qic_sloppy + j, fn_mass,
				gr_even);
      node0_printf("Solving sloppily for all ODD displacements for mass %g\n", mass);
      block_current_stochastic( j_mu[j], mass, nsrc, -1, ODD, qic_sloppy+ j, fn_mass,
				gr_odd);
      
      /* Next, continue to a "precise" solution from the same sources */
      node0_printf("Solving precisely for all EVEN displacements for mass %g\n", mass);

      block_current_stochastic( j_mu[j], mass, nsrc, +1, EVEN, qic_precise + j, fn_mass,
				gr_even);
      node0_printf("Solving precisely for all ODD displacements for mass %g\n", mass);
      block_current_stochastic( j_mu[j], mass, nsrc, +1, ODD, qic_precise + j, fn_mass,
				gr_odd);
    } /* j */
  } /* jrand */
  
  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  free(gr_even); gr_even = NULL;
  free(gr_odd); gr_odd = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
}

/************************************************************************/
/* Entry point for multiple masses with deflation and iterated single-mass inverter.
   This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  Designed for use with deflation or eigcg. 
   Requires a set of accurate low-mode eigenpairs */

void 
f_meas_current_diff( int n_masses, int nrand, int nwrite, int thinning,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy,
		     ks_param *ksp, fermion_links_t *fl, 
		     char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_diff";
  node0_printf("Entered %s\n", myname); fflush(stdout);

  double wtime = 0.;
  
  /* Create fields for current densities, one for each mass */
  Real *j_mu[n_masses];
  for(int j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open file(s) for writing */
  QIO_Writer *outfile[n_masses];
  for(int j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Calculate and write at intervals of nwrite random values */
  for(int jrand = 0; jrand < nrand; jrand += nwrite){

    block_current_diff(j_mu, nwrite, thinning, n_masses,
		       qic_precise, qic_sloppy, ksp, fl);
   
    wtime -= dclock();
    for(int j = 0; j < n_masses; j++){
#if 0      
      /* DEBUG */
      FORALLFIELDSITES(i){
	printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	for(mu = 0; mu < NMU; mu++)
	  printf("%d %g ",mu, j_mu[j][NMU*i + mu]);
	printf("\n");
      }
#endif
      node0_printf("For rand %d and mass %g\n", jrand, mass);
      average_vector_current(nwrite, j_mu[j]);
      int status = write_vector_current_record(outfile[j], jrand, nwrite, ksp[j].mass, j_mu[j]);
      clear_r_array_field(j_mu[j], NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } 
    } /* j */
  } /* jrand */

  for(int j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  fflush(stdout);
  
} /* f_meas_current_diff */

/*********************************************************************/

static void
exact_current(Real *jlow_mu[], int n_masses, ks_param *ksp, fermion_links_t *fl){

  /* Compute exact low-mode current density */

  int r_offset[4] = {0, 0, 0, 0};
  int *spin_taste = get_spin_taste();
  double dtime = -dclock();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  int i;
  imp_ferm_links_t **fn = get_fm_links(fl);

  for(int n = 0; n < Nvecs; n++){
    for(int j = 0; j < n_masses; j++){
      Real mass = ksp[j].mass;
      imp_ferm_links_t *fn_mass = fn[j];
      dslash_fn_field(eigVec[n], gr0, ODD, fn_mass);
      for(int mu = 0; mu < NMU; mu++){

        spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, gr0);
        spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
        FOREVENFIELDSITES(i){
          complex z;
          z = su3_dot( eigVec[n] + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass*mass);
        } /* i */

        spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
        spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
        FORODDFIELDSITES(i){
          complex z;
          z = su3_dot( gr0 + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += z.imag/(eigVal[n]+4.0*mass*mass);
        } /* i */

      } /* mu */
    } /* j */
  } /* n */

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;

#if 0
  for(int j = 0; j < n_masses; j++){
    Real mass = ksp[j].mass;
    for(int mu = 0; mu < NMU; mu++){
      node0_printf("Exact low modes For mass %g\n", mass);
      FORALLFIELDSITES(i){
	node0_printf("j_mu_low  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, jlow_mu[j][NMU*i+mu]);
      }
    }
  }
#endif

  dtime += dclock();
  node0_printf("Time for exact low modes %g sec\n", dtime);
}

/*********************************************************************/
static void 
block_current( Real *j_mu[], int nwrite, int n_masses, int thinning,
	       quark_invert_control *qic, ks_param *ksp, fermion_links_t *fl){

  char myname[] = "block_current";
  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  /* Block solver parameters -- temporary */
  int nr = 2;  /* Number of random sources to block */
  if(nwrite < nr)nr = nwrite;
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;

  su3_vector **gr_even = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **gr_odd = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    gr_even[is] = create_v_field();
    gr_odd[is] = create_v_field();
  }

  /* Loop over random sources in groups of nr */
  su3_vector *gr_mu = create_v_field();

  for(int jrand = 0; jrand < nwrite; jrand += nr){
    
    /* Create sources in gr_even and gr_odd */
    collect_sources(gr_even, gr_odd, nr, d, evol);

    /* Construct current density from the list of sources */
    imp_ferm_links_t **fn = get_fm_links(fl);
    for(int j = 0; j < n_masses; j++){
      Real mass = ksp[j].mass;
      imp_ferm_links_t *fn_mass = fn[j];
      
      node0_printf("Solving for all EVEN displacements for mass %g\n", mass);
      block_current_stochastic( j_mu[j], mass, nsrc, +1, EVEN, qic + j, fn_mass,
				gr_even);
      node0_printf("Solving for all ODD displacements for mass %g\n", mass);
      block_current_stochastic( j_mu[j], mass, nsrc, +1, ODD, qic + j, fn_mass,
				gr_odd);
    } /* j */
  } /* jrand */

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  destroy_v_field(gr_mu); gr_mu = NULL;
}

/*********************************************************************/
/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */

void 
f_meas_current( int n_masses, int nrand, int nwrite, int thinning,
		quark_invert_control *qic, ks_param *ksp,
		fermion_links_t *fl, char filenames[][MAXFILENAME]){

  char myname[] = "f_meas_current";
  node0_printf("Entered %s\n", myname); fflush(stdout);

  int i;
  double wtime = 0.;
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;

#if 0
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(int j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(((i == j) && (fabs(cc.real - 1) > 1e-8)) || ((i != j && fabs(cc.real) > 1e-8)))
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  Real *j_mu[n_masses];
  Real *jlow_mu[n_masses];
  for(int j = 0; j < n_masses; j++){
    j_mu[j] = create_r_array_field(NMU);
    if(Nvecs > 0)
      jlow_mu[j] = create_r_array_field(NMU);
    else
      jlow_mu[j] = NULL;
  }

  /* Open file(s) for writing */
  QIO_Writer *outfile[n_masses];
  for(int j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Compute exact low-mode current density if we have eigenvectors to do it */
  if(Nvecs > 0)
    exact_current(jlow_mu, n_masses, ksp, fl);
  
  /* Calculate and write at intervals of nwrite random values */
  for(int jrand = 0; jrand < nrand; jrand += nwrite){

    block_current( j_mu, nwrite, n_masses, thinning, qic, ksp, fl );
      
    wtime -= dclock();
    for(int j = 0; j < n_masses; j++){
      Real mass = ksp[j].mass;
#if 0
      /* DEBUG */
      node0_printf("For rand %d and mass %g\n", jrand, mass);
      FORSOMEFIELDPARITY(i, parity){
	printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	for(mu = 0; mu < NMU; mu++)
	  printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	printf("\n");
      }
#endif
      if(Nvecs > 0)
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j]);
      else
	average_vector_current(nwrite, j_mu[j]);
      int status = write_vector_current_record(outfile[j], jrand, nwrite, mass, j_mu[j]);
      if(Nvecs > 0)
	clear_r_array_field(j_mu[j], NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } 
    } /* j */
  } /* jrand */
  
  for(int j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
} /* f_meas_current */


#if 0

/* Traditional versions */

/************************************************************************/
/* Entry point for multiple masses with deflation and iterated single-mass inverter.
   This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  Designed for use with deflation or eigcg. 
   Requires a set of accurate low-mode eigenpairs */

void 
f_meas_current_diff( int n_masses, int nrand, int nwrite, int thinning,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy,
		     ks_param *ksp, fermion_links_t *fl, 
		     char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_diff";

  int i, j, jrand, mu;
  int parity, otherparity;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  double wtime = 0.;

  node0_printf("Entered %s\n", myname); fflush(stdout);
 
  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

  /* Open file(s) for writing */
  for(j = 0; j < n_masses; j++){
    outfile[j] = open_vector_current_file(filenames[j]);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
      exit(1);
    }
  }

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){

    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    
    /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){
	    parity = (ex+ey+ez+et)%2==0?EVEN:ODD;
	    switch(parity){
	    case(EVEN): otherparity=ODD; break;
	    case(ODD):  otherparity=EVEN; break;
	    }
	    
	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    if(param.eigen_param.Nvecs > 0)
	      project_out(gr, eigVec, Nvecs_tot, parity);

	    for(j = 0; j < n_masses; j++){
	      Real mass = ksp[j].mass;
	      imp_ferm_links_t *fn_mass = fn[j];

	      /* First, the sloppy high-mode solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
	      clear_v_field(M_inv_gr);
	      qic_sloppy[j].parity = parity;
	      ks_congrad_field( gr, M_inv_gr, qic_sloppy + j, mass, fn_mass);
	      dslash_fn_field( M_inv_gr, M_inv_gr, otherparity, fn_mass);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = -imag[gr * M_inv_gr] */
		/* SUBTRACT the sloppy result */
		FORSOMEFIELDPARITY(i, parity){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] -= -cc.imag;
		}
	      } /* mu */
	      
	      /* Next, continue to a "precise" solution from the same source */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	      qic_precise[j].parity = parity;
	      ks_congrad_field( gr, M_inv_gr, qic_precise + j, mass, fn_mass);
	      dslash_fn_field( M_inv_gr, M_inv_gr, otherparity, fn_mass);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = -imag[gr * gr_mu] */
		/* ADD the precise result, which then gives the difference */
		FORSOMEFIELDPARITY(i, parity){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += -cc.imag;
		  //printf("j_mu %d %d %d %d %d %g %g\n",
		  //	 lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, -cc.real, -cc.imag);
		}
	      } /* mu */
#if 0
	      /* DEBUG */
	      FORSOMEFIELDPARITY(i,parity){
		printf("diff j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
		for(mu = 0; mu < NMU; mu++)
		  printf("%g ",j_mu[j][NMU*i + mu]);
		printf("\n");
	      }
#endif
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write at intervals of nwrite random values */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	Real mass = ksp[j].mass;
#if 0
	/* DEBUG */
	FORALLFIELDSITES(i){
	  printf("write diff j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g ",mu, j_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	node0_printf("For rand %d and mass %g\n", jrand, mass);
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass, j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	}
      } /* j */
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_diff -- non-BLOCKCG version */

/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */

void 
f_meas_current( int n_masses, int nrand, int nwrite, int thinning,
		quark_invert_control *qic, ks_param *ksp,
		fermion_links_t *fl, char filenames[][MAXFILENAME]){

  char myname[] = "f_meas_current";

  int i, j, is, jr, jrand, mu, n;
  int parity, otherparity;

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
  Real *jlow_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  double wtime = 0.;
  int Nvecs = param.eigen_param.Nvecs;

  node0_printf("Entered %s\n", myname); fflush(stdout);

#if 0
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, ODD) ;
      if(((i == j) && (fabs(cc.real - 1) > 1e-8)) || ((i != j) && (fabs(cc.real) > 1e-8)))
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++){
    j_mu[j] = create_r_array_field(NMU);
    if(Nvecs > 0)
      jlow_mu[j] = create_r_array_field(NMU);
    else
      jlow_mu[j] = NULL;
  }

  /* Open file(s) for writing */
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
  
  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    for(j = 0; j < n_masses; j++){
      dslash_fn_field(eigVec[n], gr, ODD, fn_multi[j]);
      for(mu = 0; mu < NMU; mu++){

        spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr);
        spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
        FOREVENFIELDSITES(i){
          complex z;
          z = su3_dot( eigVec[n] + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
        } /* i */

        spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, eigVec[n]);
        spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
        FORODDFIELDSITES(i){
          complex z;
          z = su3_dot( gr + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
        } /* i */

      } /* mu */
    } /* j */
  } /* n */
  
#if 0
  for(j = 0; j < n_masses; j++){
    for(mu = 0; mu < NMU; mu++){
      node0_printf("For mass %g and mu %d\n", mass[j], mu);
      FORALLFIELDSITES(i){
	node0_printf("j_mu_low  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, jlow_mu[j][NMU*i+mu]);
      }
    }
  }
#endif

  dtime += dclock();
  node0_printf("Time for exact low modes %g sec\n", dtime);

  /* HACK to get only result from low modes  */
  if(nrand == 0){
    for(j = 0; j < n_masses; j++){
      average_vector_current_and_sum(1, j_mu[j], jlow_mu[j]);
      int status = write_vector_current_record(outfile[j], 0, 1, mass[j], j_mu[j]);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } 
      clear_r_array_field(j_mu[j], NMU);
    }
  }

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
    /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif

    /* Iterate over displacements within a d^4 cube.*/
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){
	    parity = (ex+ey+et+ez)%2==0?EVEN:ODD;
	    switch(parity){
	    case(EVEN): otherparity=ODD; break;
	    case(ODD):  otherparity=EVEN; break;
	    }
	    
	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    if(Nvecs > 0)
	      project_out(gr, eigVec, Nvecs, parity);
	    
#if 0
	    /* DEBUG */
	    /* Check the norm of the reduced source */
	    double_complex dd;
	    dot_product(gr, gr, &dd, parity);
	    node0_printf("Deflated source norm %g\n", dd.real);
#endif

	    for(j = 0; j < n_masses; j++){

	      /* First, the sloppy high-mode solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j]);
	      clear_v_field(M_inv_gr);
	      qic[j].parity = parity;
	      ks_congrad_field( gr, M_inv_gr, qic + j, mass[j], fn_multi[j]);
	      dslash_fn_field( M_inv_gr, M_inv_gr, otherparity, fn_multi[j]);
	      
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		
		/* J_mu = -imag[gr * gr_mu] */
		/* Add the result */
		FORSOMEFIELDPARITY(i, parity){
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += -cc.imag;
		  //printf("j_mu %d %d %d %d %d %g %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
		}
	      } /* mu */
	      
#if 0
	      /* DEBUG */
	      FORSOMEFIELDPARITY(i,parity){
		printf("j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
		for(mu = 0; mu < NMU; mu++)
		  printf("%g ",j_mu[j][NMU*i + mu]);
		printf("\n");
	      }
#endif
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write at intervals of nwrite random values */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
#if 0
	/* DEBUG */
	FORALLFIELDSITES(i){
	  printf("write j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);
	if(Nvecs > 0)
	  average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j]);
	else
	  average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(Nvecs > 0)
	  clear_r_array_field(j_mu[j], NMU);
	wtime += dclock();
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	}
	wtime += dclock();
      } /* j */
    } /* if write */
  } /* jrand */

  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    if(Nvecs > 0)
      destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current -- non-BLOCKCG version */

#endif  /* End of traditional version */

