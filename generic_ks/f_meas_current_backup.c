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

    f_meas_current_diff
    f_meas_current
    f_meas_current_pt2pt
    f_meas_current_multi_diff
    f_meas_current_multi

    With EIGMODE defined:

    f_meas_current_multi_diff_eig
    f_meas_current_multi_diff
    f_meas_current_multi_eig
    f_meas_current_multi
    f_meas_current__multi_pt2pt <- under construction
    f_meas_current_multi_diff_eig_corr
    f_meas_current_multi_eig_multi_solve

    Control Structure:

    # For single mass:
    f_meas_current_diff
    f_meas_current 
    f_meas_current_pt2pt

    # For multi-mass:
    if EIGMODE
      if BLOCKCG
        f_meas_current_multi_diff_eig
	f_meas_current_multi_eig   
      else
        f_meas_current_multi_diff_eig
	f_meas_current_multi_eig
	f_meas_current_multi_diff_eig_eo
	f_meas_current_multi_eig_eo
	f_meas_current_multi_diff_eig_corr
	f_meas_current_multi_eig_multi_solve
      f_meas_current_multi_diff <- use iterated single-mass inverter
      f_meas_current_multi      <- use iterated single-mass inverter  
      f_meas_current_multi_pt2pt<- use iterated single-mass inverter: under construction
    else
      f_meas_current_multi_diff <- use multi-mass inverter 
      f_meas_current_multi      <- use multi-mass inverter
      f_meas_current_multi_pt2pt<- use iterated single-mass inverter: need to be tested
*/

/* definitions files and prototypes */
#include "generic_ks_includes.h"	
#include "../include/fn_links.h"
#include "../include/io_scidac.h"
#include "../include/imp_ferm_links.h"
#include <qio.h>
#include <string.h>
#include <stdio.h>

#define NMU 4
#define NRECINFO 128
#define FORALLFIELDSITESLIN(x,y,z,t,a,b,c,d) for(x=0;x<a;x++)for(y=0;y<b;y++)for(z=0;z<c;z++)for(t=0;t<d;t++)

/************************************ HELPER FUNCTIONS *****************************/
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
    if(s->x % thinning != ex || s->y % thinning != ey || s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
  }
}

/* Thin the random source w/r/t color too*/
static void
thin_source_color(su3_vector *src, int thinning, int ex, int ey, int ez, int et, int color){
  site *s;
  int i;

  FORALLSITES(i,s) {
    if(s->x % thinning != ex || s->y % thinning != ey || s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
    else{
      for(int shift=1; shift < 3; shift++)
	src[i].c[(color + shift)%3] = (complex) {0.0,0.0};
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
      //printf("j_mu  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, j_mu[NMU*i+mu]);
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
  if(PRECISION == 1)
    status = write_F_R_from_field(outfile, recxml, j_mu, NMU);
  else
    status = write_D_R_from_field(outfile, recxml, j_mu, NMU);
  QIO_string_destroy(recxml);

  node0_printf("Wrote current density for source %d and mass %g\n", jrand, mass);
  
  return status;
}


/************************************ PUBLIC FUNCTIONS *****************************/
void 
f_meas_current_diff( int nrand, int nwrite, int thinning, 
		     quark_invert_control *qic_precise, quark_invert_control *qic_sloppy,
		     Real mass, int naik_term_epsilon_index, fermion_links_t *fl, 
		     char *filename){
  
  char myname[] = "f_meas_current_diff";

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
    node0_printf("%s: Failed to open %s\n", myname, filename);
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
    /* Iterate over 8 even displacements within a d^4 cube */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      //r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );
	    
	    /* First, sloppy solution */
	    /* M_inv_gr = M^{-1} gr */
	    node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
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
	    node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
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
	    FOREVENFIELDSITES(i){
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
      average_vector_current(nwrite, j_mu);
      int status = write_vector_current_record(outfile, jrand, nwrite, mass, j_mu);
      clear_r_array_field(j_mu, NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filename);
      }
    }
  } /* jrand */

  close_vector_current_file(outfile);
  node0_printf("Time to write %d records = %e sec\n", nrand, wtime);
  
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
    node0_printf("%s: Failed to open %s\n", myname, filename);
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
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      //r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );
	    
	    /* M_inv_gr = M^{-1} gr */
	    node0_printf("Solving for %d %d %d %d\n", ex, ey, ez, et);
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
	    FOREVENFIELDSITES(i){
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
      average_vector_current(nwrite, j_mu);
      int status = write_vector_current_record(outfile, jrand, nwrite, mass, j_mu);
      clear_r_array_field(j_mu, NMU);
      wtime += dclock();
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filename);
      }
    } /* if write */
  } /* jrand */

  close_vector_current_file(outfile);
  node0_printf("Time to write %d records = %e sec\n", nrand, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
  destroy_r_array_field(j_mu, NMU);
}

/************************************************************************
 * Entry point for single masse.  Uses the multimass inverter.
 * This function computes current explicitly by working with pt. sources.
 ************************************************************************/
void 
f_meas_current_pt2pt( quark_invert_control *qic, Real mass, 
		      int naik_term_epsilon_index, fermion_links_t *fl,
		      char *filename){
  
  char myname[] = "f_meas_current_pt2pt";

  int i, j, mu, color;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr;
  int x, y, z, t;

  Real j_mu[NMU];
  //QIO_Writer *outfile;
  FILE *outfile;

  /* Create vector fields */
  M_inv_gr = create_v_field();

  /* Create fields for current densities, one for each mass */
  //j_mu = create_r_array_field(NMU);

  /* Open files for writing */
  //outfile = open_vector_current_file(filename);
  if(this_node==0){
  outfile = fopen(filename,"w");

  if(outfile == NULL){
    node0_printf("%s: Failed to open %s\n", myname, filename);
    exit(1);
  }
  }

  /* Load pointers for fermion links, based on Naik epsilon indices */
  imp_ferm_links_t* fn = get_fm_links(fl)[naik_term_epsilon_index];
  
  double wtime = 0.;
  /* Loop over pt sources */
  FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
    for(mu = 0; mu < NMU; mu++) j_mu[mu]=0;
    /* Loop over colors to take the trace */
    for(color=0;color<3;color++){
      clear_v_field(gr);
      if(this_node==node_number(x,y,z,t)) gr[node_index(x,y,z,t)].c[color]=(complex){1,0};
      
      /* M_inv_gr = M^{-1} gr  */
      clear_v_field(M_inv_gr);
      mat_invert_uml_field_projected( gr, M_inv_gr, qic, (x+y+z+t)%2==0?EVEN:ODD, mass, fn );
      /* Apply current in various directions at the sink */
      for(mu = 0; mu < NMU; mu++){
	/* Apply the appropriate spin_taste operator for
	   a nearly conserved current.  */
	spin_taste_op_fn(fn, spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	spin_taste_op_fn(fn, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	/* J_mu = imag[gr.M_inv_gr] */
	if(this_node==node_number(x,y,z,t)){
	  complex cc = su3_dot( gr+node_index(x,y,z,t), gr_mu+node_index(x,y,z,t) );
          j_mu[mu] += cc.imag;
	}
      } /* mu */
    } /* color */
    g_vecfloatsum(j_mu,NMU);

    wtime -= dclock();
    for(mu=0;mu<NMU;mu++)
      if(this_node==0) fprintf(outfile,"J(%d,%d,%d,%d)_%d= %g\n",x,y,z,t,mu,j_mu[mu]); 
    wtime += dclock();

  } /* FORALLFIELDSITESLIN */
  
  //close_vector_current_file(outfile);
  fclose(outfile);
  
  node0_printf("Time to write %d records = %e sec\n", 1, wtime);

  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
} /* f_meas_current_pt2pt: singlemass inverter */

//#if EIGMODE == EIGCG || EIGMODE == DEFLATE
#ifdef EIGMODE
/*********************** HELPER FUNCTIONS ************************************/
/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
		   double_complex *dot, int parity) {
  register double re,im ;
  register site *s;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORSOMEPARITY(i,s,parity){
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

  register site *s;
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEPARITY(i,s,parity){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  }
}

/************************************************************************/
/* Projects out the vectors from the  vec. Num is the Number of vectors  
   and parity is the parity on which we work on.                           
 * The vectors are assumed to be orthonormal.                              */
static void project_out(su3_vector *vec, su3_vector **vector, int Num, int parity){
  register int i ;
  double_complex cc ;
  double ptime = -dclock();

  for(i=Num-1;i>-1;i--){
    dot_product(vector[i], vec, &cc, parity) ;
    complex_vec_mult_sub(&cc, vector[i], vec, parity);
  }

  ptime += dclock();
  node0_printf("Time to project out low modes %g sec\n", ptime);
}

/*********************** ENTRY POINTS ***********************************/
#ifdef BLOCKCG //EIGMODE && BLOCKCG
/************************************************************************
 * Entry point for multiple masses with deflation and iterated single-mass inverter.
 * This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  
 * Designed for use with deflation or eigcg. 
 * Requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
*************************************************************************/
void 
f_meas_current_multi_diff_eig( int n_masses, int nrand, int nwrite, int thinning,
			       quark_invert_control *qic_precise,
			       quark_invert_control *qic_sloppy, 
			       su3_vector **eigVec, double *eigVal, int Nvecs,
			       ks_param *ksp, fermion_links_t *fl, 
			       char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_diff_eig";

  int i, j, is, jr, jrand, mu;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  double wtime = 0.;
  
  /* Block solver parameters -- temporary */
  int nr = 2;  /* Number of random sources to block */
  su3_vector **gr;  /* Storage for sources */
  su3_vector **M_inv_gr; /* Storage for solutions */

  node0_printf("Entered %s\n", myname); fflush(stdout);
 
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

  /* Allocate source vectors */
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  M_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(is = 0; is < nsrc; is++){
    gr[is] = create_v_field();
    M_inv_gr[is] = create_v_field();
  }
  
  /* Loop over random sources in groups of nr */
  for(jrand = 0; jrand < nrand; jrand += nr){
    /* Block of random sources */
    for(jr = 0; jr < nr; jr++){
      /* Make random source, and do inversion */
#ifndef Z2RSOURCE
      grsource_plain_field( gr0, EVENANDODD );
#else
      z2rsource_plain_field( gr0, EVENANDODD );
#endif
      /* Iterate over displacements within a d^4 cube. Use even displacements only */
      for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
	if((ex+ey+ez+et)%2==0){
	  is = (et + d*(ez + d*(ey + d*ex)))/2 + evol*jr;
	  /* Apply source thinning */
	  copy_v_field( gr[is], gr0 );
	  thin_source( gr[is], d, ex, ey, ez, et );
	  /* Project out the low mode part, based on the given eigenvectors */
	  project_out(gr[is], eigVec, Nvecs, EVEN);
	}
    } /* jr */
      
    for(j = 0; j < n_masses; j++){
      /* First, the sloppy high-mode solution */
      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
      node0_printf("Solving sloppily\n");
      for(is = 0; is < nsrc; is++)
	clear_v_field(M_inv_gr[is]);
      mat_invert_block_uml_field( nsrc, gr, M_inv_gr, qic_sloppy + j, mass[j], fn_multi[j]);
      /* For each source, apply current in various directions at the sink */
      for(is = 0; is < nsrc; is++)
	for(mu = 0; mu < NMU; mu++){
	  /* Apply the appropriate spin_taste operator for
	     a nearly conserved current.  */
	  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  /* SUBTRACT the sloppy result */
	  FORALLFIELDSITES(i){
	    complex cc = su3_dot( gr[is]+i, gr_mu+i );
	    j_mu[j][NMU*i + mu] -= cc.imag;
	  } /* mu */
	} /* is */
      
      /* Next, continue to a "precise" solution from the same source */
      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
      node0_printf("Solving precisely\n");
      mat_invert_block_uml_field( nsrc, gr, M_inv_gr, qic_precise + j, mass[j], fn_multi[j]);
      /* For each source, apply current in various directions at the sink */
      for(is = 0; is < nsrc; is++)
	for(mu = 0; mu < NMU; mu++){
	  /* Apply the appropriate spin_taste operator for
	     a nearly conserved current.  */
	  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  /* ADD the precise result, which then gives the difference */
	  FORALLFIELDSITES(i){
	    complex cc = su3_dot( gr[is]+i, gr_mu+i );
	    j_mu[j][NMU*i + mu] += cc.imag;
	  } /* mu */
	} /* is */
    } /* j */
    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
#if 0      
	/* DEBUG */
	FOREVENFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g ",mu, j_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif

	/* Write record */
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  for(is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(gr[is]);
  }
  free(M_inv_gr); M_inv_gr = NULL; 
  free(gr); gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_diff_eig: BLOCKCG version */

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
*************************************************************************/
void 
f_meas_current_multi_eig( int n_masses, int nrand, int nwrite, int thinning,
			  quark_invert_control *qic,
			  su3_vector **eigVec, double *eigVal, int Nvecs,
			  ks_param *ksp, fermion_links_t *fl, 
			  char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_eig";

  int i, j, is, jr, jrand, mu, n;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_masses];
  Real *jlow_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  double wtime = 0.;

  /* Block solver parameters -- temporary */
  int nr = 2;  /* Number of random sources to block */
  su3_vector **gr;  /* Storage for sources */
  su3_vector **M_inv_gr; /* Storage for solutions */

  node0_printf("Entered %s\n", myname); fflush(stdout);

#if 1
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++){
    j_mu[j] = create_r_array_field(NMU);
    jlow_mu[j] = create_r_array_field(NMU);
  }

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
  
  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    for(j = 0; j < n_masses; j++){
      dslash_fn_field(eigVec[n], gr0, ODD, fn_multi[j]);
      for(mu = 0; mu < NMU; mu++){
	/* Add in the exact low-mode solution */
        spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr0);
        spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	
        FOREVENFIELDSITES(i){
          complex z;
          z = su3_dot( eigVec[n] + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
        } /* i */
      } /* mu */
    } /* j */
  } /* n */

#if 0
  /* DEBUG */
  for(j = 0; j < n_masses; j++){
    for(mu = 0; mu < NMU; mu++){
      node0_printf("For mass %g\n", mass[j]);
      FOREVENFIELDSITES(i){
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


  /* Allocate source vectors */
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  M_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(is = 0; is < nsrc; is++){
    gr[is] = create_v_field();
    M_inv_gr[is] = create_v_field();
  }
  
  /* Loop over random sources in groups of nr */
  for(jrand = 0; jrand < nrand; jrand += nr){
    /* Block of random sources */
    for(jr = 0; jr < nr; jr++){
      /* Make random source, and do inversion */
#ifndef Z2RSOURCE
      grsource_plain_field( gr0, EVENANDODD );
#else
      z2rsource_plain_field( gr0, EVENANDODD );
#endif
      
      /* DEBUG */
      //    clear_v_field(gr0);
      //    FOREVENFIELDSITES(i){
      //      gr0[i].c[jrand].real = 1.0;
      //    }
      
      /* Iterate over displacements within a d^4 cube. Use even displacements only */
      for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
	if((ex+ey+ez+et)%2==0){
	  is = (et + d*(ez + d*(ey + d*ex)))/2 + evol*jr;
	  /* Apply source thinning */
	  copy_v_field( gr[is], gr0 );
	  thin_source( gr[is], d, ex, ey, ez, et );
	  
	  /* Project out the low mode part, based on the given eigenvectors */
	  project_out(gr[is], eigVec, Nvecs, EVEN);
	  
#if 1
	  /* DEBUG */
	  /* Check the norm of the reduced source */
	  double_complex dd;
	  dot_product(gr[is], gr[is], &dd, EVEN);
	  node0_printf("Deflated source norm %g\n", dd.real);
#endif
	}
    } /* jr */
    
    for(j = 0; j < n_masses; j++){
      
      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
      node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j]);
      
      for(is = 0; is < nsrc; is++)
	clear_v_field(M_inv_gr[is]);
      
      mat_invert_block_uml_field( nsrc, gr, M_inv_gr, qic + j, mass[j], fn_multi[j]);
      
      /* For each source, apply current in various directions at the sink */
      for(is = 0; is < nsrc; is++){
#if 0
	/* DEBUG */
	su3_vector *M_inv_gr_test = create_v_field();
	for(n = 0; n < Nvecs; n++){
	  complex cc = {0., 0.};
	  FOREVENFIELDSITES(i){
	    complex z;
	    z = su3_dot(eigVec[n] + i, gr[is] + i);
	    CSUM(cc, z);
	  }
	  CMULREAL(cc, 2*mass[j]/(eigVal[n] + 4*mass[j]*mass[j]), cc);
	  FOREVENFIELDSITES(i){
	    for(int c = 0; c < 3; c++){
	      complex z;
	      CMUL(cc, eigVec[n][i].c[c], z);
	      CSUM(M_inv_gr_test[i].c[c], z);
	    }
	  }
	}
	destroy_v_field(M_inv_gr_test);
#endif
	
	/* Apply current in various directions at the sink */
	for(mu = 0; mu < NMU; mu++){
	  
	  /* Apply the appropriate spin_taste operator for
	     a nearly conserved current.  */
	  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  
	  /* J_mu = imag[gr.M_inv_gr] */
	  FORALLFIELDSITES(i){
	    complex cc = su3_dot( gr[is]+i, gr_mu+i );
	    j_mu[j][NMU*i + mu] += cc.imag;
	    //		  printf("j_mu %d %d %d %d %d %g %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
	  }
	  
#if 0
	  /* DEBUG */
	  su3_vector *gr_mu_test = create_v_field();
	  su3_vector *grp = create_v_field();
	  su3_vector *grpp = create_v_field();
	  
	  for(n = 0; n < Nvecs; n++){
	    dslash_fn_field(eigVec[n], grp, ODD, fn_multi[j]);
	    spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, grpp, grp);
	    spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, grpp, grpp);
	    
	    complex cc = {0., 0.};
	    FOREVENFIELDSITES(i){
	      complex z;
	      z = su3_dot(eigVec[n] + i, gr[is] + i);
	      CSUM(cc, z);
	    }
	    CMULREAL(cc, -1./(eigVal[n] + 4*mass[j]*mass[j]), cc);
	    FOREVENFIELDSITES(i){
	      for(int c = 0; c < 3; c++){
		complex z;
		CMUL(cc, grpp[i].c[c], z);
		CSUM(gr_mu_test[i].c[c], z);
	      } /* c */
	    } /* i */
	  } /* n */
	  
	  destroy_v_field(grpp);
	  destroy_v_field(grp);
	  destroy_v_field(gr_mu_test);
#endif
	} /* mu */
      } /* is */	    
    } /* j */

    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
#if 0
	/* DEBUG */
	FOREVENFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  for(is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(gr[is]);
  }
  free(M_inv_gr); M_inv_gr = NULL; 
  free(gr); gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig: BLOCKCG version */


#else //EIGMODE BUT NOT BLOCKCG

/************************************************************************
 * Entry point for multiple masses with deflation and iterated single-mass inverter.
 * This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  
 * Designed for use with deflation or eigcg. 
 * Requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
 ************************************************************************/
void 
f_meas_current_multi_diff_eig( int n_masses, int nrand, int nwrite, int thinning,
			       quark_invert_control *qic_precise,
			       quark_invert_control *qic_sloppy, 
			       su3_vector **eigVec, double *eigVal, int Nvecs,
			       ks_param *ksp, fermion_links_t *fl, 
			       char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_diff_eig";

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
  QIO_Writer * outfile[n_masses]; // = malloc( n_masses * sizeof(QIO_Writer*) );
  imp_ferm_links_t **fn = get_fm_links(fl);

  node0_printf("Entered %s\n", myname); fflush(stdout);
 
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
  
  double wtime = 0.0;

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
     /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
     /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){
	      
      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, EVEN);

	    for(j = 0; j < n_masses; j++){
	      /* First, the sloppy high-mode solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic_sloppy + j, EVEN, mass[j], fn_multi[j]);
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
	      node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic_precise + j, EVEN, mass[j], fn_multi[j]);
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		/* J_mu = imag[gr.M_inv_gr] */
		/* ADD the precise result, which then gives the difference */
		FORALLFIELDSITES(i){// if i corresponds to an odd site, j_mu[k][NMU*i + mu]=0
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write record */	    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);fflush(stdout);
	average_vector_current(nwrite, j_mu[j]);node0_printf("hhh\n");fflush(stdout);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);fflush(stdout);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current_multi_diff_eig: iterated single-mass inverter */


/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
 ************************************************************************/
void 
f_meas_current_multi_eig( int n_masses, int nrand, int nwrite, int thinning,
			  quark_invert_control *qic,
			  su3_vector **eigVec, double *eigVal, int Nvecs,
			  ks_param *ksp, fermion_links_t *fl, 
			  char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_eig";

  int i, j, n, jrand, mu;

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
 
  node0_printf("Entered %s\n", myname); fflush(stdout);

#if 1
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
      dot_product(eigVec[i], eigVec[j], &cc, ODD) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
        node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++){
    j_mu[j] = create_r_array_field(NMU);
    jlow_mu[j] = create_r_array_field(NMU);
  }

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
  
  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    for(j = 0; j < n_masses; j++){
      dslash_fn_field(eigVec[n], gr, ODD, fn_multi[j]);
      for(mu = 0; mu < NMU; mu++){
	/* Add in the exact low-mode solution */
        spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr);
        spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
        FOREVENFIELDSITES(i){
          complex z;
          z = su3_dot( eigVec[n] + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
        } /* i */
      } /* mu */
    } /* j */
  } /* n */

#if 0
  for(j = 0; j < n_masses; j++){
    for(mu = 0; mu < NMU; mu++){
      node0_printf("For mass %g\n", mass[j]);
      FOREVENFIELDSITES(i){
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

    /* DEBUG */
    //    clear_v_field(gr0);
    //    FOREVENFIELDSITES(i){
    //      gr0[i].c[jrand].real = 1.0;
    //    }
    
    /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, EVEN);
	    
#if 1
	    /* DEBUG */
	    /* Check the norm of the reduced source */
	    double_complex dd;
	    dot_product(gr, gr, &dd, EVEN);
	    node0_printf("Deflated source norm %g\n", dd.real);
#endif
	    
	    for(j = 0; j < n_masses; j++){
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j]);
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic + j, EVEN, mass[j], fn_multi[j]);
	      
#if 0
	      /* DEBUG */
	      su3_vector *M_inv_gr_test = create_v_field();
	      for(n = 0; n < Nvecs; n++){
		complex cc = {0., 0.};
		FOREVENFIELDSITES(i){
		  complex z;
		  z = su3_dot(eigVec[n] + i, gr + i);
		  CSUM(cc, z);
		}
		CMULREAL(cc, 2*mass[j]/(eigVal[n] + 4*mass[j]*mass[j]), cc);
		FOREVENFIELDSITES(i){
		  for(int c = 0; c < 3; c++){
		    complex z;
		    CMUL(cc, eigVec[n][i].c[c], z);
		    CSUM(M_inv_gr_test[i].c[c], z);
		  }
		}
	      }
	      destroy_v_field(M_inv_gr_test);
#endif

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
		  //		  printf("j_mu %d %d %d %d %d %g %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
		}

#if 0
		/* DEBUG */
		su3_vector *gr_mu_test = create_v_field();
		su3_vector *grp = create_v_field();
		su3_vector *grpp = create_v_field();

		for(n = 0; n < Nvecs; n++){
		  dslash_fn_field(eigVec[n], grp, ODD, fn_multi[j]);
		  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, grpp, grp);
		  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, grpp, grpp);
		  
		  complex cc = {0., 0.};
		  FOREVENFIELDSITES(i){
		    complex z;
		    z = su3_dot(eigVec[n] + i, gr + i);
		    CSUM(cc, z);
		  }
		  CMULREAL(cc, -1./(eigVal[n] + 4*mass[j]*mass[j]), cc);
		  FOREVENFIELDSITES(i){
		    for(int c = 0; c < 3; c++){
		      complex z;
		      CMUL(cc, grpp[i].c[c], z);
		      CSUM(gr_mu_test[i].c[c], z);
		    } /* c */
		  } /* i */
		} /* n */
		
		destroy_v_field(grpp);
		destroy_v_field(grp);
		destroy_v_field(gr_mu_test);
#endif
	      } /* mu */
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
#if 0
	/* DEBUG */
	FOREVENFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with deflation and iterated single-mass inverter.
   This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  Designed for use with deflation or eigcg. 
 * Requires a set of accurate low-mode eigenpairs
 * Computes current density at sites of both parities
 ************************************************************************/
void 
f_meas_current_multi_diff_eig_eo( int n_masses, int nrand, int nwrite, int thinning,
				  quark_invert_control *qic_precise,
				  quark_invert_control *qic_sloppy, 
				  su3_vector **eigVec, double *eigVal, int Nvecs,
				  ks_param *ksp, fermion_links_t *fl,
				  char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_diff_eig_eo";

  int i, j, jrand, mu, parity;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  int charges[3] = {1,-1,2}; //better if we take this as an input
  Real mass[n_masses];
  Real *j_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer * outfile[n_masses]; // = malloc( n_masses * sizeof(QIO_Writer*) );
  imp_ferm_links_t **fn = get_fm_links(fl);

  node0_printf("Entered %s\n", myname); fflush(stdout);
 
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
  
  double wtime = 0.0;

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
     /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
     /* Iterate over displacements within a d^4 cube. Use displacements of the given parity only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){
	    parity = (ex+ey+ez+et)%2==0?EVEN:ODD;
      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, parity);

	    for(j = 0; j < n_masses; j++){
	      /* First, the sloppy high-mode solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic_sloppy + j, parity, mass[j], fn_multi[j]);
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		/* J_mu = imag[gr.M_inv_gr] */
		/* SUBTRACT the sloppy result */
		FORALLFIELDSITES(i){// if parity(i)!=parity(ex,ey,ez,et), j_mu[k][NMU*i + mu]=0
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] -= cc.imag;
		}
	      } /* mu */
	      
	      /* Next, continue to a "precise" solution from the same source */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic_precise + j, parity, mass[j], fn_multi[j]);
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		/* J_mu = imag[gr.M_inv_gr] */
		/* ADD the precise result, which then gives the difference */
		FORALLFIELDSITES(i){// if parity(i)!=parity(ex,ey,ez,et), j_mu[k][NMU*i + mu]=0
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  j_mu[j][NMU*i + mu] += cc.imag;
		}
	      } /* mu */
/****** ADD DEBUG BEGIN ************************/
#if 0
              /* DEBUG */
              FORSOMEFIELDPARITY(i,parity){
                printf("diff j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
                for(mu = 0; mu < NMU; mu++)
                  printf("%g ",j_mu[j][NMU*i + mu]);
                printf("\n");
              }
#endif
/****** ADD DEBUG END ************************/
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write record */	    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);fflush(stdout);
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);fflush(stdout);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current_multi_diff_eig_eo: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Computes current density at sites of both parities
 ************************************************************************/
//#define TESTEO // tests if the solving on both even and odd sites is working.
#define FINDVAR
//#define FINDCOV
void 
f_meas_current_multi_eig_eo( int n_masses, int nrand, int nwrite, int thinning,
			     quark_invert_control *qic,
			     su3_vector **eigVec, double *eigVal, int Nvecs,
			     ks_param *ksp, fermion_links_t *fl, 
			     char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_eig_eo";

  int i, j, n, jrand, mu, parity;

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
 
#if defined(FINDVAR) || defined(FINDCOV)
  Real *j_mur[n_masses];
  Real jt_mu[nt*NMU];
#endif
#ifdef FINDVAR
  Real *jt2_mu[n_masses];
#endif
#ifdef FINDCOV
  //Real *cov[n_masses][NMU][sites_on_node];
  int svolume=nx*ny*nz;
  Real cov[n_masses][NMU][nt][svolume][svolume];
#endif

  node0_printf("Entered %s\n", myname); fflush(stdout);

#if 1
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
      dot_product(eigVec[i], eigVec[j], &cc, ODD) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
        node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++){
    j_mu[j] = create_r_array_field(NMU);
    jlow_mu[j] = create_r_array_field(NMU);
#if defined(FINDVAR) || defined(FINDCOV)
    j_mur[j] = create_r_array_field(NMU);
    if(j==0)for(i=0;i<nt*NMU;i++) jt_mu[i]=0;
#endif
#ifdef FINDVAR
    jt2_mu[j] = (Real *) malloc(nt*NMU*sizeof(Real));
    for(i=0;i<nt*NMU;i++) jt2_mu[j][i]=0;
#endif
#ifdef FINDCOV
    for(int k1=0;k1<svolume;k1++)
      for(int k2=0;k2<svolume;k2++)
	for(mu=0;mu<NMU;mu++)
	  for(int t=0;t<nt;t++)
	    cov[j][mu][t][k1][k2]=0;
#if 0
    FORALLFIELDSITES(i){
      cov[j][mu][i] = (Real *) malloc(nx*ny*nz*sizeof(Real));
    }
#endif
#if 0
    if(this_node==0) 
      for(int i=0;i<volume;i++)for(mu=0;mu<NMU;mu++){
	  cov[j][mu][i/(nx*ny*nz)][i%(nx*ny*nz)] = (Real *) malloc((i%nx*ny*nz+1)*sizeof(Real));
	  for(int k=0;k<i%nx*ny*nz+1;k++) cov[j][mu][i/(nx*ny*nz)][i%(nx*ny*nz)][i]=0;
	}
#endif
#endif
  }

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
  
  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    for(j = 0; j < n_masses; j++){
#if 0 // shorter version
      dslash_fn_field(eigVec[n], gr, EVENANDODD, fn_multi[j]);
      for(mu = 0; mu < NMU; mu++){
	/* Add in the exact low-mode solution */
	spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr);
	spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	FORALLFIELDSITES(i){
	    complex z;
	    z = su3_dot( eigVec[n] + i, gr_mu + i);
	    jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
	} /* i */
      } /* mu */
#endif // less affected by inaccuracy of odd EV's
      dslash_fn_field(eigVec[n], gr, ODD, fn_multi[j]);
      for(int parity=1;parity<3;parity++){
	for(mu = 0; mu < NMU; mu++){
	  if(parity == EVEN) spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr);
	  else               spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, eigVec[n]);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  FORSOMEFIELDPARITY(i,parity){
	    complex z;
	    if(parity == EVEN) z = su3_dot( eigVec[n] + i, gr_mu + i);
	    else               z = su3_dot( gr + i, gr_mu + i);
	    jlow_mu[j][NMU*i + mu] += pow(-1,parity-1)*z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
	  } /* i */
	} /* mu */
      } /* parity */
    } /* j */
  } /* n */

#if 0
  for(j = 0; j < n_masses; j++){
    for(mu = 0; mu < NMU; mu++){
      node0_printf("For mass %g\n", mass[j]);
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

  //Real ran=0;
  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
     /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif

#if defined(FINDCOV) || defined(FINDVAR)
    for(j=0;j<n_masses;j++) clear_r_array_field(j_mur[j],NMU);
#endif
#ifdef TESTEO /* DEBUG */
    for(int color=0;color<3;color++){
      clear_v_field(gr0);
      if(jrand==0){
	parity=EVEN;
	if(this_node==node_number(0,0,0,0)) gr0[node_index(0,0,0,0)].c[color]=(complex){1,0};
      }
      else if(jrand==1){
	parity=ODD;
	if(this_node==node_number(1,0,0,0)) gr0[node_index(1,0,0,0)].c[color]=(complex){1,0};
      }
      copy_v_field(gr, gr0);
#else
    //if(this_node==0) ran+=gr0[0].c[0].real*gr0[0].c[1].real+gr0[0].c[0].imag*gr0[0].c[1].imag;
    //for(int color=0;color<3;color++) node0_printf("gr_%d= (%e,%e)\n",color,gr0[0].c[color].real,gr0[0].c[color].imag);
    /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)for(int color=0;color<3;color++){
	      //for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++){
	    parity = (ex+ey+et+ez)%2==0?EVEN:ODD;
      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    //thin_source( gr, d, ex, ey, ez, et );
	    thin_source_color(gr, d, ex, ey, ez, et, color);
#endif
	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, parity);
	    
#if 1	    /* DEBUG */
	    /* Check the norm of the reduced source */
	    double_complex dd;
	    dot_product(gr, gr, &dd, EVEN);
	    node0_printf("Deflated source norm %g\n", dd.real);
#endif
	    
	    for(j = 0; j < n_masses; j++){
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j]);
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic + j, parity, mass[j], fn_multi[j]);
	      
#if 0 	      /* DEBUG */
	      su3_vector *M_inv_gr_test = create_v_field();
	      for(n = 0; n < Nvecs; n++){
		complex cc = {0., 0.};
		FOREVENFIELDSITES(i){
		  complex z;
		  z = su3_dot(eigVec[n] + i, gr + i);
		  CSUM(cc, z);
		}
		CMULREAL(cc, 2*mass[j]/(eigVal[n] + 4*mass[j]*mass[j]), cc);
		FOREVENFIELDSITES(i){
		  for(int c = 0; c < 3; c++){
		    complex z;
		    CMUL(cc, eigVec[n][i].c[c], z);
		    CSUM(M_inv_gr_test[i].c[c], z);
		  }
		}
	      }
	      destroy_v_field(M_inv_gr_test);
#endif

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
#ifdef TESTEO
		  //if(j_mu[j][NMU*i + mu]!=0) printf("j_mu %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, j_mu[j][NMU*i + mu]);
		  if(cc.real!=0 || cc.imag!=0) printf("j_mu %d %d %d %d %d %g %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
#endif
#if defined(FINDCOV) || defined(FINDVAR)
		  j_mur[j][NMU*i + mu] += cc.imag;
#endif
		}

#if 0		/* DEBUG */
		su3_vector *gr_mu_test = create_v_field();
		su3_vector *grp = create_v_field();
		su3_vector *grpp = create_v_field();

		for(n = 0; n < Nvecs; n++){
		  dslash_fn_field(eigVec[n], grp, ODD, fn_multi[j]);
		  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, grpp, grp);
		  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, grpp, grpp);
		  
		  complex cc = {0., 0.};
		  FOREVENFIELDSITES(i){
		    complex z;
		    z = su3_dot(eigVec[n] + i, gr + i);
		    CSUM(cc, z);
		  }
		  CMULREAL(cc, -1./(eigVal[n] + 4*mass[j]*mass[j]), cc);
		  FOREVENFIELDSITES(i){
		    for(int c = 0; c < 3; c++){
		      complex z;
		      CMUL(cc, grpp[i].c[c], z);
		      CSUM(gr_mu_test[i].c[c], z);
		    } /* c */
		  } /* i */
		} /* n */
		
		destroy_v_field(grpp);
		destroy_v_field(grp);
		destroy_v_field(gr_mu_test);
#endif

	      } /* mu */

/****** ADD DEBUG BEGIN ************************/
#if 0
              /* DEBUG */
              FORSOMEFIELDPARITY(i,parity){
                printf("j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
                for(mu = 0; mu < NMU; mu++)
                  printf("%g ",j_mu[j][NMU*i + mu]);
                printf("\n");
              }
#endif
/****** ADD DEBUG END ************************/
	    } /* j */
	  } /* ex, ey, ez, et */

#if defined(FINDVAR) || defined(FINDCOV)
    for(j=0;j<n_masses;j++){
      for(i=0;i<NMU*nt;i++) jt_mu[i]=0;
      for(mu=0;mu<NMU;mu++){
#ifdef FINDVAR
	for(int t=0;t<nt;t++)
	  FORALLFIELDSITES(i){
	    if(lattice[i].t==t) jt_mu[NMU*t+mu]+=j_mur[j][NMU*i+mu];
	  }
      }
      g_vecfloatsum(jt_mu,NMU*nt);
      for(i=0;i<NMU*nt;i++) jt2_mu[j][i]+=jt_mu[i]*jt_mu[i]/nwrite;
#endif
#ifdef FINDCOV
#if 0
      int x2,y2,z2,t2,x,y,z,t;
      Real jtmp;
      FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	FORALLFIELDSITESLIN(x2,y2,z2,t2,nx,ny,nz,1){
	  jtmp=0;
	  if(this_node==node_number(x2,y2,z2,t)) jtmp+=j_mur[j][NMU*node_index(x2,y2,z2,t)+mu];
	  g_floatsum(&jtmp);
	  if(this_node==node_number(x,y,z,t))
	    cov[j][mu][node_index(x,y,z,t)][x2+y2*nx+z2*nx*ny]+=j_mur[j][NMU*node_index(x,y,z,t)+mu]*jtmp/nwrite;
	  g_sync();
	}
      }
#endif
#if 0
      int x,y,z,t;
      Real jtmp1=0,jtmp2=0;
      FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	if(node_number(x,y,z,t)==0) jtmp1=j_mur[j][NMU*node_index(x,y,z,t)+mu]; 
	else{
	  if(this_node==node_number(x,y,z,t)) send_field((char*) &j_mur[j][NMU*node_index(x,y,z,t)+mu],sizeof(Real),0);
	  if(this_node==0) get_field((char *) &jtmp1,sizeof(Real),node_number(x,y,z,t));
	  }
	for(i=0;i<x+y*nx+z*nx*ny+1;i++){
	  int x2=i%nx,y2=(i/nx)%ny,z2=i/(nx*ny);
	  if(node_number(x,y,z,t)==0) jtmp2=j_mur[j][NMU*node_index(x2,y2,z2,t)+mu];
	  else{
	    if(this_node==node_number(x2,y2,z2,t)) send_field((char*) &j_mur[j][NMU*node_index(x2,y2,z2,t)+mu],sizeof(Real),0);
	      if(this_node==0) get_field((char *) &jtmp2,sizeof(Real),node_number(x2,y2,z2,t));
	  }
	  if(this_node==0) cov[j][mu][t][x+y*nx+z*nx*ny][i]+=jtmp1*jtmp2/nwrite;
	}
	g_sync();
      }
#endif
    } /* mu */
#endif
    } /* j:n_masses */
#endif

    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
#ifdef TESTEO   /* DEBUG */
	Real *msg;
	for(mu=0;mu<NMU;mu++){
	  if(jrand==0 && this_node==node_number(0,0,0,0)) msg=&j_mu[j][NMU*node_index(0,0,0,0)+mu];
	  if(jrand==1){
	    if(node_number(0,0,0,0)!=node_number(1,0,0,0)){
	      if(this_node==node_number(1,0,0,0)){
		msg=&j_mu[j][NMU*node_index(1,0,0,0)+mu];
		send_field((char *)msg,sizeof(Real),0);
	      }
	      if(this_node==node_number(0,0,0,0)) get_field((char *)msg,sizeof(Real),node_number(1,0,0,0));
	    }
	    else if(this_node==node_number(0,0,0,0)) msg=&j_mu[j][NMU*node_index(1,0,0,0)+mu];
	  }
	  g_sync();
	  node0_printf("J_%d= %g for mass %d\n",mu,*msg,j);fflush(stdout);
	}
#endif
#if defined(FINDVAR) || defined(FINDCOV)  /* DEBUG */
	int x,y,z,t;
	double msg[NMU];
	for(int k=0;k<nt*NMU;k++){
	  jt_mu[k]=0;
	  FORALLFIELDSITES(i){
	    if(lattice[i].t==k/NMU) jt_mu[k]+=j_mu[j][NMU*i+k%NMU];
	  }
	}
	g_vecfloatsum(jt_mu,nt*NMU);
#ifdef FINDVAR
	for(i=0;i<NMU*nt;i++){
	  jt2_mu[j][i]-=jt_mu[i]*jt_mu[i]/nwrite/nwrite;
	  node0_printf("j_t(%d)_%d= %e %e\n",i/NMU,i%NMU,jt_mu[i]/nwrite,sqrt(jt2_mu[j][i]/nwrite));
	}
#endif
	// print j_mu at all sites
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  for(mu=0;mu<NMU;mu++) msg[mu]=0;
	  if(this_node==node_number(x,y,z,t))for(mu=0;mu<NMU;mu++) msg[mu]=j_mu[j][NMU*node_index(x,y,z,t)+mu];
	  g_vecfloatsum(msg,NMU);
	  for(mu=0;mu<NMU;mu++) 
	    node0_printf("j(%d,%d,%d,%d)_%d= %.16e\n",x,y,z,t,mu,msg[mu]/nwrite);
	}
#ifdef FINDCOV //computing the covariance matrix
	node0_printf("start2\n");fflush(stdout);
	int  x2,y2,z2,t2,x,y,z,t;
	Real jtmp=0;
#if 0
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  FORALLFIELDSITESLIN(x2,y2,z2,t2,nx,ny,nz,1){
	    jtmp=0;
	    if(this_node==node_number(x2,y2,z2,t)) jtmp+=j_mur[j][NMU*node_index(x2,y2,z2,t)+mu];
	    g_floatsum(&jtmp);
	    if(this_node==node_number(x,y,z,t))//find the node the matrix entry is located
	      cov[j][mu][node_index(x,y,z,t)][x2+y2*nx+z2*nx*ny]-=j_mu[j][NMU*node_index(x,y,z,t)+mu]*jtmp/nwrite/nwrite;
	  }
	}
#endif
#if 0 // the version with some bugs leading to "MPI_Barrier(comm=0xc4000002) failed"
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  FORALLFIELDSITESLIN(x2,y2,z2,t2,nx,ny,nz,1){
	    if(node_number(x,y,z,t)!=node_number(x2,y2,z2,t)){// the current is somewhere else.
	      if(this_node==node_number(x2,y2,z2,t)) //Are you the sender?
		send_field((char*) &j_mu[j][NMU*node_index(x2,y2,z2,t)+mu],sizeof(Real),node_number(x,y,z,t));
	      else if(this_node==node_number(x,y,z,t))// If not, are you the receiver?
		get_field((char *) &jtmp,sizeof(Real),node_number(x2,y2,z2,t));
	    }
	    else
	      jtmp=j_mu[j][NMU*node_index(x2,y2,z2,t)+mu];
	    if(this_node==node_number(x,y,z,t))//find the node the matrix entry is located
	      cov[j][mu][node_index(x,y,z,t)][x2+y2*nx+z2*nx*ny]-=j_mu[j][NMU*node_index(x,y,z,t)+mu]*jtmp/nwrite/nwrite;
	  }
	}
#endif
#if 0
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  FORALLFIELDSITESLIN(x2,y2,z2,t2,nx,ny,nz,1){
	    Real jtmp2=0;
	    for(i=0;i<x+y*nx+z*nx*ny+1;i++){
	      int x2=i%nx,y2=((i-i%nx)/nx)%ny,z2=i/(nx*ny);
	      if(node_number(x,y,z,t)==0) jtmp2=j_mu[j][NMU*node_index(x2,y2,z2,t)+mu];
	      else{
		if(this_node==node_number(x2,y2,z2,t)) send_field((char*) &j_mu[j][NMU*node_index(x2,y2,z2,t)+mu],sizeof(Real),0);
		if(this_node==0) get_field((char *) &jtmp2,sizeof(Real),node_number(x2,y2,z2,t));
	      }
	      for(mu=0;mu<NMU;mu++) if(this_node==0) cov[j][mu][t][x+y*nx+z*nx*ny][i]-=msg[mu]*jtmp2/nwrite/nwrite;
	    }
	  }
	}
#endif
	//Find std
#if 0
	Real std2=0;
	for(t=0;t<nt;t++)
	  for(mu=0;mu<NMU;mu++){
	    FORALLFIELDSITES(i){
	      if(lattice[i].t==t) 
		FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,1){
		  std2+=cov[j][mu][i][x+y*nx+z*nx*ny];
		}
	    }
	    g_floatsum(&std2);
	    node0_printf("j_t(%d)_%d= %e %e\n",t,mu,j_t[t*NMU+mu],sqrt(std2));
	  }
#endif
#if 0
	if(this_node==0){
	  Real std2=0;
	  for(t=0;t<nt;t++)
	    for(mu=0;mu<NMU;mu++){
	      for(i=0;i<nx*ny*nz;i++)
		for(int k=0;k<i+1;k++)
		  std2+=pow(2,(int) i!=k)*cov[j][mu][t][i][j];
	      node0_printf("j_t(%d)_%d= %e %e\n",t,mu,j_t[i/NMU+i%NMU],sqrt(std2));
	    }
	}
#endif
#endif
	g_sync();
#endif
#if 0	/* DEBUG */
	FOREVENFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	node0_printf("For rand %d and mass %g\n", jrand, mass[j]);
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  //node0_printf("grxgr= %e\n",ran/nrand);
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_r_array_field(jlow_mu[j], NMU);
#ifdef FINDVAR
    destroy_r_array_field(j_mur[j],NMU);
    free(jt2_mu[j]);
  }
   jt2_mu[0]=NULL;
#else
  }
#endif
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig_eo: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with deflation and iterated single-mass inverter.
 * This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities with multiple sloppy precisions on sites with the spcified parity.
 * It also computes correlation of sloppy and precise current densities for
   multiple  sloppy precisions.  
 * If numMassPair != 0, the difference of the two current densities for a pair of
   masses (m1,m2) is estimated using the mass-difference trick
   * The difference is J_m2-J_m1
   * In this case, the first numMassPair entries of ksp are assumed to be the 
     pairs for which the trick is applied.
 * Designed for use with deflation or eigcg. 
 * Requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
 ************************************************************************/
void 
f_meas_current_multi_diff_eig_corr( int n_masses, int const nrand, int nwrite, int thinning,
				    quark_invert_control *qic_precise,
				    quark_invert_control *qic_sloppy, 
				    su3_vector **eigVec, double *eigVal, int Nvecs, ks_param *ksp, 
				    fermion_links_t *fl, Real *precs_sloppy, int n_slp,
				    char filenames[][MAXFILENAME],
				    int parity, int numMassPair){
  
  char myname[] = "f_meas_current_multi_diff_eig_corr";

  int i, j, jrand, mu, ns;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  int charge[n_masses];
  Real mass[n_masses];
  Real *jd_mu[n_slp][n_masses];
  Real *cov_mu[n_slp];
  Real *js_mu[n_slp];
  Real *jsa_mu[n_slp];
  Real *js_mu_var[n_slp];
  Real *jp_mu;
  Real *jpa_mu;
  Real *jp_mu_var;
  Real jsa[n_slp],jpa=0;
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer * outfile[n_slp*n_masses]; // = malloc( n_slp * n_masses * sizeof(QIO_Writer*) );
  FILE **cfps = (FILE **) malloc(n_slp*n_masses * sizeof(FILE*));  
  imp_ferm_links_t **fn = get_fm_links(fl);


  node0_printf("Entered %s\n", myname); fflush(stdout);
 
  /* Create fields for current densities, one for each mass */
  for(ns=0;ns<n_slp;ns++){
    for(j = 0; j < n_masses; j++) jd_mu[ns][j] = create_r_array_field(NMU);
    cov_mu[ns]=create_r_array_field(NMU);
    js_mu[ns] = create_r_array_field(NMU);
    jsa_mu[ns]=create_r_array_field(NMU);
    js_mu_var[ns]=create_r_array_field(NMU);
    jsa[ns]=0;
  }
  jp_mu = create_r_array_field(NMU);
  jpa_mu=create_r_array_field(NMU);
  jp_mu_var=create_r_array_field(NMU);

  /* Open files for writing */ 
  char *p;
  char es[50];
  char fname[MAXFILENAME];
  char ftmp[MAXFILENAME];
  char fnames[n_slp*n_masses][MAXFILENAME];
  for(j = 0; j < n_slp*n_masses; j++){
    strcpy(fname,filenames[j%n_masses]);
    p=strstr(fname,"esls"); //to be replaced by esls
    *p='\0';    
    p=strstr(p+1,"epls");
    strcpy(ftmp,p);
    sprintf(es,"%s%.2e%s%.2e","esls",precs_sloppy[(j/n_masses)*n_masses],"esc",precs_sloppy[(j/n_masses)*n_masses+2]);
    strcat(fname,es);
    strcat(fname,ftmp);
    outfile[j] = open_vector_current_file(fname);
    strcpy(fnames[j],fname);
    node0_printf("%s\n",fname);fflush(stdout);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, fname);
      exit(1);
    }
#if 0
    if(this_node == 0){
      while((p=strstr(fname,"vcd")) != NULL){
        *p='\0';
        p+=3;
        strcpy(ftmp,p);
        strcat(fname,"vcc");
        strcat(fname,ftmp);
        cfps[j]=fopen(fname,"w");
        if((cfps[j]) == NULL){
          node0_printf("%s: Failed to open %s\n", myname, fname);
          exit(1);
        }
      }
    }
#endif
  }

  /* Load parameters */
  for(j = 0; j < n_masses; j++){
    /* Load charges from ks_param */
    charge[j] = ksp[j].charge;
    /* Load masses from ks_param */
    mass[j] = ksp[j].mass;
    /* Load pointers for fermion links, based on Naik epsilon indices */
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  }

  double wtime = 0.0;

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
     /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
     /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){
	      
      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, EVEN);

	    for(j = 0; j < n_masses; j++){
	      /* First, the sloppy high-mode solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
	      clear_v_field(M_inv_gr);
	      for(ns=0;ns < n_slp; ns++){
		(qic_sloppy+j)->resid = precs_sloppy[j+n_masses*ns];
		mat_invert_uml_field_projected( gr, M_inv_gr, qic_sloppy + j, EVEN, mass[j], fn_multi[j]);
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
		    jd_mu[ns][j][NMU*i + mu] -= cc.imag;
		    js_mu[ns][NMU*i + mu] += 2*charge[j]*cc.imag/3;// The factor of 2 comes from correcting for MILC normalization convention of M
		    jsa[ns] += 2*charge[j]*cc.imag/3;// The factor of 2 comes from correcting for MILC normalization convention of M
		  }/* i:field sites*/
		}/* mu */
	      }/* ns:n_slp */

	      /* Next, continue to a "precise" solution from the same source */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	      mat_invert_uml_field_projected( gr, M_inv_gr, qic_precise + j, EVEN, mass[j], fn_multi[j]);
	      /* Apply current in various directions at the sink */
	      for(mu = 0; mu < NMU; mu++){
		/* Apply the appropriate spin_taste operator for
		   a nearly conserved current.  */
		spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		/* J_mu = imag[gr.M_inv_gr] */
		/* ADD the precise result, which then gives the difference */
		FORALLFIELDSITES(i){// if i corresponds to an odd site, j_mu[k][NMU*i + mu]=0
		  complex cc = su3_dot( gr+i, gr_mu+i );
		  for(ns=0;ns<n_slp;ns++) jd_mu[ns][j][NMU*i + mu] += cc.imag;
		  jp_mu[NMU*i + mu] += 2*charge[j]*cc.imag/3;// The factor of 2 comes from correcting for MILC normalization convention of M
		  jpa += 2*charge[j]*cc.imag/3;// The factor of 2 comes from correcting for MILC normalization convention of M
		}/* i:field sites*/
	      } /* mu */
	    } /* j: n_masses */
	    
	    /* compute covariance and variance at each site */
	    for(mu=0;mu<NMU;mu++){
	      node0_printf("j(jrand=%d)_%d: ",jrand,mu);
	      FORALLFIELDSITES(i){
		for(ns=0;ns<n_slp;ns++){
		  cov_mu[ns][NMU*i+mu] += js_mu[ns][NMU*i + mu]*jp_mu[NMU*i + mu];
		  jsa_mu[ns][NMU*i + mu] += js_mu[ns][NMU*i + mu];
		  js_mu_var[ns][NMU*i + mu] += js_mu[ns][NMU*i + mu]*js_mu[ns][NMU*i + mu];
		  if(i==node_index(0,0,0,0)) node0_printf("%.16e  ",js_mu[ns][NMU*i + mu]);
		  js_mu[ns][NMU*i + mu]=0;
		}
		jpa_mu[NMU*i + mu] += jp_mu[NMU*i + mu];
		jp_mu_var[NMU*i + mu] += jp_mu[NMU*i + mu]*jp_mu[NMU*i + mu];
		if(i==node_index(0,0,0,0)) node0_printf("%.16e at the origin\n",jp_mu[NMU*i + mu]);
		jp_mu[NMU*i + mu]=0;
	      }
	    }
	  } /* ex, ey, ez, et *

    /* Write record */
    /* No division by #thinned sublattices:
       This is b/c current density value at a site recieves a contribution only from a single sublattice.
       If we divide by the factor, the current density  will not be properly normalized.
     */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_masses; j++){
	node0_printf("For rand %d and mass %g\n", jrand, mass[j%n_masses]);fflush(stdout);
	average_vector_current(nwrite, jd_mu[j/n_masses][j%n_masses]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j%n_masses], jd_mu[j/n_masses][j%n_masses]);
	clear_r_array_field(jd_mu[j/n_masses][j%n_masses], NMU);//node0_printf("hhh3\n");fflush(stdout);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, fnames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */

  /* Compute Statistics */
  /* get the average of current density values over all sites */
  g_vecfloatsum(jsa,n_slp);
  g_floatsum(& jpa);
  for(ns=0;ns<n_slp;ns++) jsa[ns] /= (nrand*volume/2*NMU); //factor of 2 comes from avg over only even sites
  jpa /= (nrand*volume/2*NMU); //factor of 2 comes from avg over only even sites
  /* get the variance and covariance field */
  for(mu=0;mu<NMU;mu++){
    FORALLFIELDSITES(i){
      for(ns=0;ns<n_slp;ns++){
        cov_mu[ns][NMU*i+mu] = (cov_mu[ns][NMU*i+mu]-jsa_mu[ns][NMU*i+mu]*jpa_mu[NMU*i+mu]/nrand)/nrand;
	js_mu_var[ns][NMU*i + mu]=(js_mu_var[ns][NMU*i + mu]-jsa_mu[ns][NMU*i + mu]*jsa_mu[ns][NMU*i + mu]/nrand)/nrand;
      }
      jp_mu_var[NMU*i + mu]=(jp_mu_var[NMU*i + mu]-jpa_mu[NMU*i + mu]*jpa_mu[NMU*i + mu]/nrand)/nrand;
    }
  }
  /* get the correlation of the current densities over sites 
     along w/ avg. of correlation, js_mu, jp_mu, js_mu_sig, and jp_mu sig over sites and comp's */
  Real js[n_slp],jp=0,js_var[n_slp],jp_var=0,nume[n_slp],r=0,r_mu[n_slp][NMU],js_mu_sig[n_slp],js_mu_siga[n_slp],jp_mu_sig=0,jp_mu_siga=0,ra[n_slp],rvar[n_slp];
  for(ns=0;ns<n_slp;ns++){
    js[ns]=0;
    js_var[ns]=0;
    js_mu_sig[ns]=0;
    js_mu_siga[ns]=0;
    nume[ns]=0;
    ra[ns]=0;
    rvar[ns]=0;
    for(mu=0;mu<NMU;mu++) r_mu[ns][mu]=0;
  }
  for(int x=0;x<nx;x++){
    for(int y=0;y<ny;y++){
      for(int z=0;z<nz;z++){
	for(int t=0;t<nt;t++){
	  if((x+y+z+t)%2==0){
	    for(mu=0;mu<NMU;mu++){
	      if(node_number(x,y,z,t) == this_node){
		i = node_index(x,y,z,t);
		for(ns=0;ns<n_slp;ns++){
		  js[ns] = jsa_mu[ns][NMU*i + mu]/nrand;
		  js_mu_sig[ns] = sqrt(js_mu_var[ns][NMU*i + mu]);
		  js_mu_siga[ns]+=js_mu_sig[ns];
		}
		jp = jpa_mu[NMU*i + mu]/nrand;
		jp_mu_sig = sqrt(jp_mu_var[NMU*i + mu]);
		jp_mu_siga+=jp_mu_sig;
		for(ns=0;ns<n_slp;ns++){
		  r_mu[ns][mu] =  cov_mu[ns][NMU*i+mu]/sqrt(js_mu_var[ns][NMU*i + mu]*jp_mu_var[NMU*i + mu]);
		  ra[ns]+=r_mu[ns][mu];
		  rvar[ns]+=r_mu[ns][mu]*r_mu[ns][mu];
		  nume[ns]   += (js[ns]-jsa[ns])*(jp-jpa);
		  js_var[ns] += (js[ns]-jsa[ns])*(js[ns]-jsa[ns]);
		}
		jp_var += (jp-jpa)*(jp-jpa);
	      }
	      g_vecfloatsum(js,n_slp);
	      g_vecfloatsum(js_mu_sig,n_slp);
	      g_floatsum(&jp);
	      g_floatsum(&jp_mu_sig);
	      for(ns=0;ns<n_slp;ns++){
		g_vecfloatsum(r_mu[ns],NMU);
		node0_printf("(ns,mu)=(%d,%d): js= %.16e js_sig= %.16e jp= %.16e jp_sig= %.16e r= %.16e at (%d,%d,%d,%d)\n",
			     ns,mu,js[ns],js_mu_sig[ns],jp,jp_mu_sig,r_mu[ns][mu],x,y,z,t);fflush(stdout);
		r_mu[ns][mu]=0;
		js[ns]=0;
		js_mu_sig[ns]=0;
	      }
	      jp=0;
	      jp_mu_sig=0;
	      g_sync();
	    }
	  }
	}
      }
    }
  }
  g_vecfloatsum(nume,n_slp);
  g_vecfloatsum(ra,n_slp);
  g_vecfloatsum(rvar,n_slp);
  g_vecfloatsum(js_var,n_slp);
  g_vecfloatsum(js_mu_siga,n_slp);
  g_floatsum(&jp_var);
  g_floatsum(&jp_mu_siga);
  /* Report Statistics */
  for(ns=0;ns<n_slp;ns++){
    r = nume[ns]/sqrt(js_var[ns]*jp_var);
    node0_printf("sloppy precision: "); for(j=0;j<n_masses;j++) node0_printf("%.2e ",precs_sloppy[j+n_masses*ns]); node0_printf("\n");
    node0_printf("statistics over all %d components of current densities on %d even sites\n",NMU,volume/2);
    node0_printf("stats: js= %.16e pm %.16e jp= %.16e pm %.16e\n",
		 jsa[ns],sqrt(js_var[ns]/NMU/volume*2/NMU/volume*2),jpa,sqrt(jp_var/NMU/volume*2/NMU/volume*2)); fflush(stdout);
    node0_printf("r: %lf\n",r);fflush(stdout);
    node0_printf("average of statistics on each sites for each components\n");
    node0_printf("avg_sigma(js,jp): %.16e %.16e\n",js_mu_siga[ns]/NMU/volume*2,jp_mu_siga/NMU/volume*2);
    node0_printf("r_avg: %lf pm sig_bar= %e\n",ra[ns]/NMU/volume*2,sqrt((rvar[ns]-ra[ns]*ra[ns]/NMU/volume*2)/NMU/volume*2/NMU/volume*2)); // stdm=std/sqrt(N)
  }

  node0_printf("\nTime to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);fflush(stdout);

  /* clean up */  
  for(j = 0; j < n_slp*n_masses; j++){
    close_vector_current_file(outfile[j]);//node0_printf("hi222hi\n");fflush(stdout);
    destroy_r_array_field(jd_mu[j/n_masses][j%n_masses], NMU);
#if 0
    if(this_node==0 && fclose(cfps[j]) == EOF) node0_printf("closing %dth vcc file failed\n",j);
#endif
  }
  for(ns=0;ns<n_slp;ns++){
    destroy_r_array_field(cov_mu[ns],NMU);
    destroy_r_array_field(js_mu[ns],NMU);
    destroy_r_array_field(jsa_mu[ns],NMU);
    destroy_r_array_field(js_mu_var[ns],NMU);
  }
  destroy_r_array_field(jp_mu,NMU);
  destroy_r_array_field(jpa_mu,NMU);
  destroy_r_array_field(jp_mu_var,NMU);

  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
  //free(outfile); outfile=NULL;
  free(cfps); cfps=NULL;

} /* f_meas_current_multi_diff_eig_corr: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * This variant does solves for the same source at multiple precisions.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
 ************************************************************************/
void 
f_meas_current_multi_eig_multi_solve( int n_masses, int nrand, int nwrite, int thinning,
				      quark_invert_control *qic,
				      su3_vector **eigVec, double *eigVal, int Nvecs,
				      ks_param *ksp, fermion_links_t *fl, 
				      Real *precs_sloppy, int n_slp,
				      char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_eig_multi_solve";

  int i, j, n, jrand, mu, ns;

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
  Real *j_mu[n_slp][n_masses];
  Real *jlow_mu[n_masses];
  imp_ferm_links_t *fn_multi[n_masses];
  QIO_Writer *outfile[n_slp*n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  double wtime = 0.;
 
  node0_printf("Entered %s\n", myname); fflush(stdout);

#if 0
  /* DEBUG */
  /* Check orthonormality of a few eigenvectors */
  for(j = 0; j < Nvecs; j += 8)
    for(i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
      dot_product(eigVec[i], eigVec[j], &cc, ODD) ;
      if(i == j && fabs(cc.real - 1) > 1e-8 || i != j && fabs(cc.real) > 1e-8)
        node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
#endif

  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++){
    for(ns = 0; ns < n_slp; ns++) j_mu[ns][j] = create_r_array_field(NMU);
    jlow_mu[j] = create_r_array_field(NMU);
  }

  /* Open files for writing */
  char *p;
  char es[50];
  char fname[MAXFILENAME];
  char ftmp[MAXFILENAME];
  char fnames[n_slp*n_masses][MAXFILENAME];
  for(j = 0; j < n_slp*n_masses; j++){
    strcpy(fname,filenames[j%n_masses]);
    p=strstr(fname,"esl");
    *p='\0';
    p=strstr(p+1,"srcsp");
    strcpy(ftmp,p);
    sprintf(es,"%s%.2e%s%.2e","esls",precs_sloppy[(j/n_masses)*n_masses],"esc",precs_sloppy[(j/n_masses)*n_masses+2]);
    strcat(fname,es);
    strcat(fname,ftmp);
    outfile[j] = open_vector_current_file(fname);
    strcpy(fnames[j],fname);
    if(outfile[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, fname);
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
	/* Add in the exact low-mode solution */
        spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, gr);
        spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
        FOREVENFIELDSITES(i){
          complex z;
          z = su3_dot( eigVec[n] + i, gr_mu + i);
	  jlow_mu[j][NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass[j]*mass[j]);
        } /* i */
      } /* mu */
    } /* j */
  } /* n */

  dtime += dclock();
  node0_printf("Time for exact low modes %g sec\n", dtime);

  /* HACK to get only result from low modes  */
  if(nrand == 0){
    for(j = 0; j < n_slp*n_masses; j++){
      average_vector_current_and_sum(1, j_mu[j/n_masses][j%n_masses], jlow_mu[j%n_masses]);
      int status = write_vector_current_record(outfile[j], 0, 1, mass[j%n_masses], j_mu[j/n_masses][j%n_masses]);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, fnames[j]);
      } 
      clear_r_array_field(j_mu[j/n_masses][j%n_masses], NMU);
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

    /* DEBUG */
    //    clear_v_field(gr0);
    //    FOREVENFIELDSITES(i){
    //      gr0[i].c[jrand].real = 1.0;
    //    }
    
    /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    /* Project out the low mode part, based on the given eigenvectors */
	    project_out(gr, eigVec, Nvecs, EVEN);
	    
#if 1
	    /* DEBUG */
	    /* Check the norm of the reduced source */
	    double_complex dd;
	    dot_product(gr, gr, &dd, EVEN);
	    node0_printf("Deflated source norm %g\n", dd.real);
#endif
	    
	    for(j = 0; j < n_masses; j++){
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j]);
	      clear_v_field(M_inv_gr);
	      for(ns=0;ns<n_slp;ns++){
		(qic+j)->resid = precs_sloppy[j+ns*n_masses];
		mat_invert_uml_field_projected( gr, M_inv_gr, qic + j, EVEN, mass[j], fn_multi[j]);
	      
#if 0
		/* DEBUG */
		su3_vector *M_inv_gr_test = create_v_field();
		for(n = 0; n < Nvecs; n++){
		  complex cc = {0., 0.};
		  FOREVENFIELDSITES(i){
		    complex z;
		    z = su3_dot(eigVec[n] + i, gr + i);
		    CSUM(cc, z);
		  }
		  CMULREAL(cc, 2*mass[j]/(eigVal[n] + 4*mass[j]*mass[j]), cc);
		  FOREVENFIELDSITES(i){
		    for(int c = 0; c < 3; c++){
		      complex z;
		      CMUL(cc, eigVec[n][i].c[c], z);
		      CSUM(M_inv_gr_test[i].c[c], z);
		    }
		  }
		}
		destroy_v_field(M_inv_gr_test);
#endif
		
		/* Apply current in various directions at the sink */
		for(mu = 0; mu < NMU; mu++){
		  /* Apply the appropriate spin_taste operator for
		     a nearly conserved current.  */
		  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
		  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
		  /* J_mu = imag[gr.M_inv_gr] */
		  FORALLFIELDSITES(i){
		    complex cc = su3_dot( gr+i, gr_mu+i );
		    j_mu[ns][j][NMU*i + mu] += cc.imag;
		  }
		
#if 0
		/* DEBUG */
		su3_vector *gr_mu_test = create_v_field();
		su3_vector *grp = create_v_field();
		su3_vector *grpp = create_v_field();

		for(n = 0; n < Nvecs; n++){
		  dslash_fn_field(eigVec[n], grp, ODD, fn_multi[j]);
		  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, grpp, grp);
		  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, grpp, grpp);
		  
		  complex cc = {0., 0.};
		  FOREVENFIELDSITES(i){
		    complex z;
		    z = su3_dot(eigVec[n] + i, gr + i);
		    CSUM(cc, z);
		  }
		  CMULREAL(cc, -1./(eigVal[n] + 4*mass[j]*mass[j]), cc);
		  FOREVENFIELDSITES(i){
		    for(int c = 0; c < 3; c++){
		      complex z;
		      CMUL(cc, grpp[i].c[c], z);
		      CSUM(gr_mu_test[i].c[c], z);
		    } /* c */
		  } /* i */
		} /* n */
		
		destroy_v_field(grpp);
		destroy_v_field(grp);
		destroy_v_field(gr_mu_test);
#endif
		} /* mu */
	      } /* ns:n_slp */
	    } /* j:mass */
	  } /* ex, ey, ez, et */

    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_masses; j++){
	node0_printf("For rand %d and mass %g\n", jrand, mass[j%n_masses]);
	average_vector_current_and_sum(nwrite, j_mu[j/n_masses][j%n_masses], jlow_mu[j%n_masses]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j/n_masses][j%n_masses]);
	clear_r_array_field(j_mu[j/n_masses][j%n_masses], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, fnames[j]);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_slp*n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j/n_masses][j%n_masses], NMU);
    if(j/n_masses==0) destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig_multi_solve: iterated single-mass inverter */

#endif //END: if-else BLOCKCG clause; still within EIGMODE

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.
 * Return: current densities at EVEN sites 
 ************************************************************************/
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
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      //r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );
	    
	    for(j = 0; j < n_masses; j++){
	      
	      /* First, sloppy solution */
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      node0_printf("Solving sloppily for %d %d %d %d\n", ex, ey, ez, et);
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field( gr, M_inv_gr, qic_sloppy + j, mass[j], fn_multi[j]);
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
	      node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	      mat_invert_uml_field( gr, M_inv_gr, qic_precise + j, mass[j], fn_multi[j]);
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
	      FOREVENFIELDSITES(i){
		printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
		for(mu = 0; mu < NMU; mu++)
		  printf("%g ",j_mu[j][NMU*i + mu]);
		printf("\n");
	      }
#endif
	    } /* j */
	  } /* ex, ey, ez, et */

    /* Write record */	    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* n_masses */
      wtime += dclock();
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
} /* f_meas_current_multi_diff: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg
 * Return: current densities at EVEN sites
 ************************************************************************/
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
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      //r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );

	    for(j = 0; j < n_masses; j++){
	      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
	      clear_v_field(M_inv_gr);
	      mat_invert_uml_field( gr, M_inv_gr, qic + j, mass[j], fn_multi[j]);
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
      
    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	clear_r_array_field(j_mu[j], NMU);
      } /* j */
      wtime += dclock();
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
} /* f_meas_current_multi: iterated single-mass inverter */

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg
need to be tested
 ************************************************************************/
void 
f_meas_current_multi_pt2pt( int n_masses, quark_invert_control *qic, 
			    ks_param *ksp, fermion_links_t *fl, 
			    char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_pt2pt";

  int i, j, jrand, mu, color;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int x, y, z, t;

  Real mass[n_masses];
  Real j_mu[n_masses*NMU];
  imp_ferm_links_t *fn_multi[n_masses];
  //QIO_Writer *outfile[n_masses];
  FILE **outfile = (FILE**) malloc(sizeof(FILE*)*n_masses);
  imp_ferm_links_t **fn = get_fm_links(fl);

 
  /* Open files for writing */
  if(this_node==0){
    for(j = 0; j < n_masses; j++){
      //outfile[j] = open_vector_current_file(filenames[j]);
      outfile[j] = fopen(filenames[j],"w");
      if(outfile[j] == NULL){
	node0_printf("%s: Failed to open %s\n", myname, filenames[j]);
	exit(1);
      }
    }
  }
  
  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;
  
  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  
  double wtime = 0.;
  /* Loop over pt sources */
  FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
    for(j = 0; j < NMU*n_masses; j++) j_mu[j]=0;
    /* Loop over colors to take the trace */
    for(color=0;color<3;color++){
      clear_v_field(gr);
      if(this_node==node_number(x,y,z,t)) gr[node_index(x,y,z,t)].c[color]=(complex){1.0,0};
      for(j = 0; j < n_masses; j++){
	/* M_inv_gr = M^{-1} gr  */
	clear_v_field(M_inv_gr);
	mat_invert_uml_field_projected( gr, M_inv_gr, qic + j, mass[j], (x+y+z+t)%2==0?EVEN:ODD, fn_multi[j]);
	/* Apply current in various directions at the sink */
	for(mu = 0; mu < NMU; mu++){
	  /* Apply the appropriate spin_taste operator for
	     a nearly conserved current.  */
	  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  if(this_node==node_number(x,y,z,t)){
	    i = node_index(x,y,z,t);
	    complex cc = su3_dot( gr+i, gr_mu+i );
	    j_mu[j*n_masses+mu] += cc.imag;
	  }
	} /* mu */
      } /* j:n_masses */
    } /* color */
    g_vecfloatsum(j_mu,n_masses*NMU);
    
    wtime -= dclock();
    for(j = 0; j < n_masses; j++)
      for(mu=0;mu<NMU;mu++)
        //node0_printf("J(%d,%d,%d,%d)_%d= %g\n",x,y,z,t,mu,j_mu[j*n_masses+mu]);
        if(this_node==0) fprintf(outfile[j],"J(%d,%d,%d,%d)_%d= %g for mass %d\n",x,y,z,t,mu,j_mu[j*n_masses+mu],j);
    wtime += dclock();

  } /* FORALLFIELDSITESLIN */
  
  for(j = 0; j < n_masses; j++)
    //close_vector_current_file(outfile[j]);
    if(this_node==0) fclose(outfile[j]);
  
  node0_printf("Time to write 1 records for %d masses = %e sec\n", n_masses, wtime);
  
  free(outfile); outfile = NULL;
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
} /* f_meas_current_multi_pt2pt: iterated single-mass inverter */

#else //NOT EIGMODE

/************************************************************************/
/* Entry point for multiple masses.  Uses the multimass inverter        */
/* Return: current densities at EVEN sites                              */
/************************************************************************/
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
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      //r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );
	    
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
	    } /* n_masses */
	    
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

    /* Write record */    
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	clear_r_array_field(j_mu[j], NMU);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
    
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);
  
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current_multi_diff: multimass inverter */

/************************************************************************/
/* Entry point for multiple masses.  Uses the multimass inverter        */
/* Return: current densities at EVEN sites                              */
/************************************************************************/
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
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if((ex+ey+ez+et)%2==0){

	      // r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

	    /* Apply source thinning */
	    copy_v_field(gr, gr0);
	    thin_source( gr, d, ex, ey, ez, et );
	    
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
	      
    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_masses; j++){
	average_vector_current(nwrite, j_mu[j]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, mass[j], j_mu[j]);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
	clear_r_array_field(j_mu[j], NMU);
      } /* j */
      wtime += dclock();
    } /* if write */
  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(j_mu[j], NMU);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }

  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current_multi: multimass inverter */

/************************************************************************
   Entry point for multiple masses.  Uses the multimass inverter.
   This function computes current explicitly by working with pt. sources.
need to be tested.
 ************************************************************************/
void 
f_meas_current_multi_pt2pt( int n_masses, quark_invert_control *qic, 
			    ks_param *ksp, fermion_links_t *fl, 
			    char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_multi_pt2pt";

  int i, j, mu, color;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr[n_masses];
  int x, y, z, t;

  Real mass[n_masses];
  Real j_mu[n_masses*NMU];
  imp_ferm_links_t *fn_multi[n_masses];
  //QIO_Writer *outfile[n_masses];
  FILE ** outfile = (FILE **) malloc(sizeof(FILE*)*n_masses);
  imp_ferm_links_t **fn = get_fm_links(fl);
 
  /* Create vector fields */
  for(j = 0; j < n_masses; j++)
    M_inv_gr[j] = create_v_field();

  /* Open files for writing */
  for(j = 0; j < n_masses; j++){
    //outfile[j] = open_vector_current_file(filenames[j]);
    outfile[j] = fopen(filenames[j],"w");
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
  /* Loop over pt sources */
  FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
    for(j = 0; j < NMU; j++) j_mu[j]=0;
    /* Loop over colors to take the trace */
    for(color=0;color<3;color++){
      clear_v_field(gr);
      if(this_node==node_number(x,y,z,t)) gr[node_index(x,y,z,t)].c[color]=(complex){1,0};

      /* M_inv_gr = M^{-1} gr  */
      total_iters += mat_invert_multi( gr, M_inv_gr, ksp, n_masses, qic, fn_multi );

      for(j = 0; j < n_masses; j++){
	/* Apply current in various directions at the sink */
	for(mu = 0; mu < NMU; mu++){
	  /* Apply the appropriate spin_taste operator for
	     a nearly conserved current.  */
	  spin_taste_op_fn(fn_multi[j], spin_taste[mu], r_offset, gr_mu, M_inv_gr[j]);
	  spin_taste_op_fn(fn_multi[j], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  if(this_node==node_number(x,y,z,t)){
	    complex cc = su3_dot( gr+node_index(x,y,z,t), gr_mu+node_index(x,y,z,t) );
	    j_mu[j*n_masses+mu] += cc.imag;
	  }
	} /* mu */
      } /* j:n_masses */
    } /* color */
    g_vecfloatsum(j_mu,NMU);

    wtime -= dclock();
    for(j = 0; j < n_masses; j++)
      for(mu=0;mu<NMU;mu++)
	//node0_printf("J(%d,%d,%d,%d)_%d= %g\n",x,y,z,t,mu,j_mu[j*n_masses+mu]);
	if(this_node==0) fprintf(outfile[j],"J(%d,%d,%d,%d)_%d= %g for mass %d\n",x,y,z,t,mu,j_mu[j*n_masses+mu],j);
    wtime += dclock();

  } /* FORALLFIELDSITESLIN */
  
  for(j = 0; j < n_masses; j++){
    //close_vector_current_file(outfile[j]);
    fclose(outfile[j]);
    destroy_v_field(M_inv_gr[j]); M_inv_gr[j] = NULL;
  }
  
  node0_printf("Time to write %d records for %d masses = %e sec\n", 1, n_masses, wtime);
  
  free(outfile); outfile=NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
} /* f_meas_current_multi_pt2pt: multimass inverter */

#endif //END EIGMODE

