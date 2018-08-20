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

#define DIFF

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
static void
write_tslice_values(char *tag, int jrand, int nwrite, Real mass, Real *j_mu){
  node0_printf("JTMU%s source %d mass %g\n", tag, jrand, mass);
  double *jtmu = (double *)malloc(sizeof(double)*4*param.nt);
  for(int tmu = 0; tmu < 4*param.nt; tmu++)
    jtmu[tmu] = 0.;

  int i;
  FORALLFIELDSITES(i){
    for(int mu = 0; mu < 4; mu++)
      jtmu[4*lattice[i].t + mu] += j_mu[4*i + mu];
  }
  for(int t = 0; t < param.nt; t++){
    node0_printf("JTMU%s %d ", tag, t);
    for(int mu = 0; mu < 4; mu++){
      g_doublesum(&jtmu[4*t+mu]);
      node0_printf("%g ", jtmu[4*t+mu]);
    }
    node0_printf("\n");
  }
}

/*Write current record for the accumulated average over random sources */
static int
write_vector_current_record(QIO_Writer *outfile, char *filename, int jrand,
			    int nwrite, Real mass, Real *j_mu){
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

  node0_printf("Wrote current density for source %d and mass %g to file %s\n",
	       jrand, mass, filename);
  
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
#ifdef CGTIME
  node0_printf("Time to project out low modes %g sec\n", ptime);
#endif
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

#if 0
	node0_printf("Source %d is %d %d %d %d\n", is, ex, ey, ez, et);
#endif

	/* Apply source thinning */
	copy_v_field( gr[is], gr0 );
	thin_source( gr[is], d, ex, ey, ez, et );
	
	/* Project out (remove) the low mode part, based on the given eigenvectors */
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
collect_sources(su3_vector *gr_even[], su3_vector *gr_odd[], int nr,
		int thinning, int evol){

  su3_vector *gr0 = create_v_field();

  for(int jr = 0; jr < nr; jr++){
    
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    //    node0_printf("EVEN sources\n");
    collect_evenodd_sources(gr_even + jr*evol, nr*evol, EVEN, thinning, gr0);
    //    node0_printf("ODD sources\n");
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

  int otherparity = ODD;  /* Initialized to humor the compiler */

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
  free(M_inv_gr);

} /* block_current_stochastic */
      
/************************************************************************/
/* Calculate the difference in current densities for a given pair of quark masses
   and source parity using a list of thinned stochastic estimators.
   Calculates mass0 result minus mass1 result. */
/* nsrc sources in gr.  Results in j_mu01 */

static void
block_current_stochastic_deltam( Real *j_mu01, Real mass0, Real mass1,
				 imp_ferm_links_t *fn01, int nsrc,
				 int sign, int parity, 
				 quark_invert_control *qic, su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic_deltam";

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  Real deltamsq4 = 4*(mass1*mass1 - mass0*mass0);

  su3_vector *M_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    M_inv_gr[is] = create_v_field();
  }

  int otherparity;
  if(parity == EVEN)
    otherparity=ODD;
  else
    otherparity=EVEN;

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *M0_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    M0_inv_gr[is] = create_v_field();
  }

  qic->parity = parity;
  /* M0_inv_gr = 1/[D^2 + 4*mass0^2] gr */
  ks_congrad_block_field(nsrc, gr, M0_inv_gr, qic, mass0, fn01);
  /* M_inv_gr = 1/[D^2 + 4*mass0^2] M0_inv_gr */
  ks_congrad_block_field(nsrc, M0_inv_gr, M_inv_gr, qic, mass1, fn01);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M0_inv_gr[is]);

  /* Might be better to use a block dslash here?? */
  for(int is = 0; is < nsrc; is++)
    dslash_fn_field( M_inv_gr[is], M_inv_gr[is], otherparity, fn01);
  
  /* For each source, apply current in various directions at the sink */
  su3_vector *gr_mu = create_v_field();
  for(int is = 0; is < nsrc; is++)
    for(int mu = 0; mu < NMU; mu++){
      
      /* Apply the appropriate spin_taste operator for
	 a nearly conserved current.  */
      spin_taste_op_fn(fn01, spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
      spin_taste_op_fn(fn01, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      /* J_mu = -imag[gr * M_inv_gr] */
      /* Multiply result by the mass difference squared times four */
      int i;
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_mu+i );
	j_mu01[NMU*i + mu] += -sign*cc.imag*deltamsq4;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu01[NMU*i + mu]);
#endif
      }
    } /* is, mu */
  destroy_v_field(gr_mu);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M_inv_gr[is]);

} /* block_current_stochastic_deltam */
      
/************************************************************************/

/* Calculate the current density difference between fine and sloppy
   solves separately for all masses */

static void
block_current_diff(int n_masses, Real *j_mu[], Real masses[], imp_ferm_links_t *fn_mass[],
		   int nwrite, int thinning, 
		   quark_invert_control *qic_precise, quark_invert_control *qic_sloppy){

  char myname[] = "block_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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
    for(int j = 0; j < n_masses; j++){
      
      /* First, the sloppy high-mode solution */
      node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", masses[j]);
      
      block_current_stochastic( j_mu[j], masses[j], nsrc, -1, EVEN, qic_sloppy + j,
				fn_mass[j], gr_even);
      //      node0_printf("j_mu[%d][0] = %g\n",j,j_mu[j][0]);
      node0_printf("Solving sloppily for all ODD displacements for mass %g\n", masses[j]);
      block_current_stochastic( j_mu[j], masses[j], nsrc, -1, ODD, qic_sloppy + j,
				fn_mass[j], gr_odd);
      //      node0_printf("j_mu[%d][4*node_index(1,0,0,0)] = %g\n",j,j_mu[j][4*node_index(1,0,0,0)]);
      
      /* Next, continue to a "precise" solution from the same sources */
      node0_printf("Solving precisely for all EVEN displacements for mass %g\n", masses[j]);

      block_current_stochastic( j_mu[j], masses[j], nsrc, +1, EVEN, qic_precise + j,
				fn_mass[j], gr_even);
      node0_printf("Solving precisely for all ODD displacements for mass %g\n", masses[j]);
      block_current_stochastic( j_mu[j], masses[j], nsrc, +1, ODD, qic_precise + j,
				fn_mass[j], gr_odd);
    } /* j */
  } /* jrand */
  
  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  free(gr_even); gr_even = NULL;
  free(gr_odd); gr_odd = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
} /* block_current_diff */

/************************************************************************/

/* Calculate the current density difference between fine and sloppy
   solves separately for three masses and take the difference: the
   mass0 result minus the mass1 result. */

static void
block_current_diff_01diff(Real *j_mu01, Real *j_mu2,
			  Real mass0, Real mass1, Real mass2,
			  imp_ferm_links_t *fn01, imp_ferm_links_t *fn2,
			  int nwrite, int thinning, 
			  quark_invert_control *qic_precise01,
			  quark_invert_control *qic_precise2,
			  quark_invert_control *qic_sloppy01,
			  quark_invert_control *qic_sloppy2){

  char myname[] = "block_current_diff_01diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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

    /* First, the sloppy high-mode solution */
    node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
		 mass0, mass1);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, -1, EVEN,
				     qic_sloppy01, gr_even);
    //    node0_printf("j_mu01[0] = %g\n",j_mu01[0]);
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", mass2);
    block_current_stochastic( j_mu2, mass2, nsrc, -1, EVEN, qic_sloppy2, fn2, gr_even);
    //    node0_printf("j_mu2[0] = %g\n",j_mu2[0]);
    node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
		 mass0, mass1);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, -1, ODD,
				     qic_sloppy01, gr_odd);
    //    node0_printf("j_mu01[4*node_index(1,0,0,0)] = %g\n",j_mu01[4*node_index(1,0,0,0)]);
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", mass2);
    block_current_stochastic( j_mu2, mass2, nsrc, -1, ODD, qic_sloppy2, fn2, gr_odd);
    //    node0_printf("j_mu2[4*node_index(1,0,0,0)] = %g\n",j_mu2[4*node_index(1,0,0,0)]);

      
    /* Next, continue to a "precise" solution from the same sources */
    node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
		 mass0, mass1);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, +1, EVEN,
				     qic_precise01, gr_even);
    node0_printf("Solving precisely for all EVEN displacements for mass %g\n", mass2);
    block_current_stochastic( j_mu2, mass2, nsrc, +1, EVEN, qic_precise2, fn2, gr_even);
    node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
		 mass0, mass1);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, +1, ODD,
				     qic_precise01, gr_odd);
    node0_printf("Solving precisely for all ODD displacements for mass %g\n", mass2);
    block_current_stochastic( j_mu2, mass2, nsrc, +1, ODD, qic_precise2, fn2, gr_odd);

  } /* jrand */
  
  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  free(gr_even); gr_even = NULL;
  free(gr_odd); gr_odd = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
} /* block_current_diff_01diff */

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
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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

  /* Load arrays with masses and the HISQ link structure for each */
  imp_ferm_links_t *fn_mass[n_masses];
  Real masses[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  for(int j = 0; j < n_masses; j++){
    masses[j] = ksp[j].mass;
    fn_mass[j] = fn[ksp[j].naik_term_epsilon_index];
  }

  /* Calculate and write at intervals of nwrite random values */
  for(int jrand = 0; jrand < nrand; jrand += nwrite){

#ifndef DIFF
    block_current_diff(n_masses, j_mu, masses, fn_mass, nwrite, thinning, 
		       qic_precise, qic_sloppy);
#else
    block_current_diff_01diff(j_mu[0], j_mu[2],
			      masses[0], masses[1], masses[2],
			      fn_mass[0], fn_mass[2], nwrite, thinning, 
			      &qic_precise[0], &qic_precise[2],
			      &qic_sloppy[0], &qic_sloppy[2]);
#endif
    
   
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
      int status = write_vector_current_record(outfile[j], filenames[j], jrand,
					       nwrite, masses[j], j_mu[j]);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	terminate(1);
      } 
      write_tslice_values("DIFF", jrand, nwrite, masses[j], j_mu[j]);
      clear_r_array_field(j_mu[j], NMU);
    } /* j */
  } /* jrand */

  for(int j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  

#ifdef CGTIME
  node0_printf("Time to write %d records for %d masses = %e sec\n",
	       nrand/nwrite, n_masses, wtime);
#endif
  fflush(stdout);
  
} /* f_meas_current_diff */

/*********************************************************************/

static void
exact_current(Real *jlow_mu, Real mass, imp_ferm_links_t *fn_mass){

  /* Compute exact low-mode current density for a single mass */

  int r_offset[4] = {0, 0, 0, 0};
  int *spin_taste = get_spin_taste();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  int i;

  for(int n = 0; n < Nvecs; n++){
    dslash_fn_field(eigVec[n], gr0, ODD, fn_mass);
    for(int mu = 0; mu < NMU; mu++){
      
      spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	jlow_mu[NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass*mass);
      } /* i */
      
      spin_taste_op_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FORODDFIELDSITES(i){
	complex z;
	z = su3_dot( gr0 + i, gr_mu + i);
	jlow_mu[NMU*i + mu] += z.imag/(eigVal[n]+4.0*mass*mass);
      } /* i */
      
    } /* mu */
  } /* n */

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

/*********************************************************************/

static void
exact_current_deltam(Real *jlow_mu01, Real mass0, Real mass1, imp_ferm_links_t *fn01){

  /* Compute the difference in exact low-mode current densities for two masses */
  /* Takes densities for mass0 minus densities for mass1 */

  int r_offset[4] = {0, 0, 0, 0};
  int *spin_taste = get_spin_taste();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  int i;
  Real deltamsq4 = 4*(mass1*mass1 - mass0*mass0);

  for(int n = 0; n < Nvecs; n++){
    dslash_fn_field(eigVec[n], gr0, ODD, fn01);
    for(int mu = 0; mu < NMU; mu++){
      
      spin_taste_op_fn(fn01, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_fn(fn01, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	Real d1 = eigVal[n]+4.0*mass0*mass0;
	Real d2 = eigVal[n]+4.0*mass1*mass1;
	jlow_mu01[NMU*i + mu] += -z.imag*deltamsq4/(d1*d2);
      } /* i */
      
      spin_taste_op_fn(fn01, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_fn(fn01, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FORODDFIELDSITES(i){
	complex z;
	z = su3_dot( gr0 + i, gr_mu + i);
	Real d1 = eigVal[n]+4.0*mass0*mass0;
	Real d2 = eigVal[n]+4.0*mass1*mass1;
	jlow_mu01[NMU*i + mu] += z.imag*deltamsq4/(d1*d2);
      } /* i */
    } /* mu */
  } /* n */

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

/*********************************************************************/

static void
exact_currents(int n_masses, Real *jlow_mu[], Real masses[],
	       imp_ferm_links_t *fn_mass[]){

  /* Compute exact low-mode current density for all masses */

  double dtime = -dclock();

  for(int j = 0; j < n_masses; j++){
    exact_current(jlow_mu[j], masses[j], fn_mass[j]);
  } /* j */

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
#ifdef CGTIME
  node0_printf("Time for exact low modes %g sec\n", dtime);
#endif
}

/*********************************************************************/

static void
exact_currents_deltam(Real jlow_mu01[], Real jlow_mu2[],
		      Real mass0, Real mass1, Real mass2,
		      imp_ferm_links_t *fn01, imp_ferm_links_t *fn2){

  /* Compute the difference in exact low-mode current densities for 
     mass0 and mass1 and compute the current density for mass2 */

  double dtime = -dclock();

  exact_current_deltam(jlow_mu01, mass0, mass1, fn01);
  exact_current(jlow_mu2, mass2, fn2);

#if 0
  for(int mu = 0; mu < NMU; mu++){
    node0_printf("Exact low modes For mass %g minus mass %g\n", mass0, mass1);
    FORALLFIELDSITES(i){
      node0_printf("j_mu_low  %d %d %d %d %d %g\n",
		   lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		   mu, jlow_mu01[NMU*i+mu]);
    }
    node0_printf("Exact low modes For mass %g\n", mass2);
    FORALLFIELDSITES(i){
      node0_printf("j_mu_low  %d %d %d %d %d %g\n",
		   lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		   mu, jlow_mu2[NMU*i+mu]);
    }
  }
#endif

  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time for exact low modes %g sec\n", dtime);
#endif
}

/*********************************************************************/
// This is the original version

static void 
block_current( int n_masses, Real *j_mu[], Real masses[],
	       imp_ferm_links_t *fn_mass[], int nwrite, int thinning,
	       quark_invert_control *qic){

  char myname[] = "block_current";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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
    for(int j = 0; j < n_masses; j++){
      node0_printf("Solving for all EVEN displacements for mass %g\n", masses[j]);
      block_current_stochastic( j_mu[j], masses[j], nsrc, +1, EVEN, qic + j, fn_mass[j],
				gr_even);
      node0_printf("Solving for all ODD displacements for mass %g\n", masses[j]);
      block_current_stochastic( j_mu[j], masses[j], nsrc, +1, ODD, qic + j, fn_mass[j],
				gr_odd);
    } /* j */
  } /* jrand */

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  destroy_v_field(gr_mu); gr_mu = NULL;
}

// This version takes explicit differences of mass0 and mass1 but not mass2
/*********************************************************************/
static void 
block_current_01diff( Real *j_mu01, Real *j_mu2,
		      Real mass0, Real mass1, Real mass2,
		      imp_ferm_links_t *fn01, imp_ferm_links_t *fn2,
		      int nwrite, int thinning,
		      quark_invert_control *qic01, quark_invert_control *qic2 ){

  char myname[] = "block_current_01diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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
    node0_printf("Solving for all EVEN displacements for mass diff  %g %g\n", mass0, mass1);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, +1, EVEN,
				     qic01, gr_even);
    node0_printf("Solving for all EVEN displacements for mass %g\n", mass2);
    block_current_stochastic( j_mu2, mass2, nsrc, +1, EVEN, qic2, fn2, gr_even);

    node0_printf("Solving for all ODD displacements for mass %g\n", mass);
    block_current_stochastic_deltam( j_mu01, mass0, mass1, fn01, nsrc, +1, ODD,
				     qic01, gr_odd);
    node0_printf("Solving for all ODD displacements for mass %g\n", mass);
    block_current_stochastic( j_mu2, mass2, nsrc, +1, ODD, qic2, fn2, gr_odd);

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
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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

  /* Load arrays with masses and the HISQ link structure for each */
  imp_ferm_links_t *fn_mass[n_masses];
  Real masses[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  for(int j = 0; j < n_masses; j++){
    masses[j] = ksp[j].mass;
    fn_mass[j] = fn[ksp[j].naik_term_epsilon_index];
  }

  /* Compute exact low-mode current density if we have eigenvectors to do it */
  if(Nvecs > 0){
#ifndef DIFF
    exact_currents(n_masses, jlow_mu, masses, fn_mass);
#else
    exact_currents_deltam(jlow_mu[0], jlow_mu[2],
      masses[0], masses[1], masses[2], fn_mass[0], fn_mass[2]);
#endif
  }
    
  /* Construct current density from the list of sources */
  /* Calculate and write at intervals of nwrite random values */
  for(int jrand = 0; jrand < nrand; jrand += nwrite){

#ifndef DIFF
    block_current( n_masses, j_mu, masses, fn_mass, nwrite, thinning, qic );
#else
    block_current_01diff( j_mu[0], j_mu[2], masses[0], masses[1], masses[2],
			  fn_mass[0], fn_mass[2], nwrite, thinning,
			  &qic[0], &qic[2] );
#endif
      
    wtime -= dclock();
    for(int j = 0; j < n_masses; j++){
#if 0
      /* DEBUG */
      node0_printf("For rand %d and mass %g\n", jrand, masses[j]);
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
      int status = write_vector_current_record(outfile[j], filenames[j], jrand,
					       nwrite, masses[j], j_mu[j]);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	terminate(1);
      }
      write_tslice_values("", jrand, nwrite, masses[j], j_mu[j]);
      clear_r_array_field(j_mu[j], NMU);
    } /* j */
  } /* jrand */
  
  for(int j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_r_array_field(jlow_mu[j], NMU);
  }
  
#ifdef CGTIME
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
#endif
  
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

  //  node0_printf("Entered %s\n", myname); fflush(stdout);
 
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
	      imp_ferm_links_t *fn_mass = fn[ksp[j].naik_term_epsilon_index];

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
	int status = write_vector_current_record(outfile[j], filenames[j], jrand,
						 nwrite, ksp[j].mass, j_mu[j]);
	wtime += dclock();
	if(status != QIO_SUCCESS){
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	  terminate(1);
	}
	write_tslice_values("DIFF", jrand, nwrite, ksp[j].mass, j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
      } /* j */
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
#ifdef CGTIME
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
#endif
  
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

  //  node0_printf("Entered %s\n", myname); fflush(stdout);

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
#ifdef CGTIME
  node0_printf("Time for exact low modes %g sec\n", dtime);
#endif

  /* HACK to get only result from low modes  */
  if(nrand == 0){
    for(j = 0; j < n_masses; j++){
      average_vector_current_and_sum(1, j_mu[j], jlow_mu[j]);
      int status = write_vector_current_record(outfile[j], filenames[j], 0,
					       1, ksp[j].mass, j_mu[j]);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	terminate(1);
      } 
      write_tslice_values("", jrand, nwrite, ksp[j].mass, j_mu[j]);
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
	int status = write_vector_current_record(outfile[j], filenames[j], jrand,
						 nwrite, mass[j], j_mu[j]);
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
  
#ifdef CGTIME
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
#endif
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current -- non-BLOCKCG version */

#endif  /* End of traditional version */

