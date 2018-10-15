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

    Control Structure:

    if EIGMODE
      if BLOCKCG
        f_meas_current_diff <- use iterated single-mass inverter 
	f_meas_current      <- use iterated single-mass inverter
      else
        f_meas_current_diff <- use iterated single-mass inverter
	f_meas_current      <- use iterated single-mass inverter
      f_meas_current_pt2pt  <- use iterated single-mass inverter: under construction (deflation not impremented)
    else
      f_meas_current_diff   <- use multi-mass inverter 
      f_meas_current        <- use multi-mass inverter
      f_meas_current_pt2pt  <- use iterated single-mass inverter: need to be tested
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

/* Thin the random source w/r/t color too if isColorThin != 0 */
static void
thin_source_color(su3_vector *src, int thinning, int ex, int ey, int ez, int et, int color, int isColorThin){
  site *s;
  int i;

  FORALLSITES(i,s) {
    if(s->x % thinning != ex || s->y % thinning != ey || s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
    else{
      for(int shift=1; shift < 3*isColorThin; shift++)
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

static int
write_total_vector_current_record_t(int jrand, Real j_mu[], FILE* outfile[], int j){
  int i,k,mu,status;
  Real jt_mu[nt*NMU];

  for(k=0;k<nt*NMU;k++) jt_mu[k]=0;
  for(mu=0;mu<NMU;mu++){
    FORALLFIELDSITES(i){
      jt_mu[lattice[i].t*NMU+mu]+=j_mu[NMU*i+mu];
    }
  }
  g_vecfloatsum(jt_mu,NMU*nt);
  if(this_node==0) for(k=0;k<nt*NMU;k++) status+=fprintf(outfile[j],"%d %d %d %.16e\n", jrand, k/NMU,k%NMU,jt_mu[k]);

  return status!=0;
}

/* Open up multiple files for multiple precision.
   mass pair corresponds to a single file.
   Assume: if 2*n_MassPair < n_masses, the last mass corresponds to charm
   Assume: original names contain the string esls#.##e+/-##esc#.##e+/-##
   Assume: all masses except for the charm are inverted w/ the same residual
           i.e., the names do not contain separate field for each mass.
 */
static void 
#ifdef TIME_CURR
open_files(char myname[], char filenames[][MAXFILENAME], FILE *outfiles[], char fnames[][MAXFILENAME],
           int n_slp, Real *precs_sloppy, int n_masses, Real masses[], int n_MassPair){
#else
open_files(char myname[], char filenames[][MAXFILENAME], QIO_Writer *outfiles[], char fnames[][MAXFILENAME],
	   int n_slp, Real *precs_sloppy, int n_masses, Real masses[], int n_MassPair){
#endif
  char *p;
  char es[50];
  char fname[MAXFILENAME];
  char ftmp[MAXFILENAME];
  int n_jmu=n_masses-n_MassPair;
  char mass[50];
  int size;
  int shift;
  for(int j = 0; j < n_slp*n_jmu; j++){
    if(j%n_jmu==0) shift=0;
    strcpy(fname,filenames[j%n_jmu+shift]);
    if(n_masses!=1){
      if(j%n_jmu<n_MassPair){
	snprintf(mass, sizeof(mass),"%f",masses[(j%n_jmu)+shift]);
	size=strlen(mass);
	p=strstr(fname,mass);
	p=p+size;
	strcpy(ftmp,p);
	*p='\0';
	strcat(fname,"-");
	snprintf(mass, sizeof(mass),"%f",masses[(j%n_jmu)+shift+1]);
	strcat(fname,mass);
	strcat(fname,ftmp);
	shift++;
      }
      p=strstr(fname,"esl");
      *p='\0';
      if(strstr(p+1,"epl")==NULL) p=strstr(p+1,"srcsp");
      else                        p=strstr(p+1,"epl");
      strcpy(ftmp,p);
      if(j%n_jmu<n_MassPair) sprintf(es,"%s%.2e","esmd",precs_sloppy[(j/n_jmu)*n_masses+j%n_jmu+shift]);
      else sprintf(es,"%s%.2e%s%.2e",
		   "esls",precs_sloppy[(j/n_jmu)*n_masses+2*shift],
		   "esc",precs_sloppy[(j/n_jmu)*n_masses+n_masses-1]); 
      strcat(fname,es);
      strcat(fname,ftmp);
    }
#ifdef TIME_CURR
    if(this_node==0) outfiles[j] = fopen(filenames[j],"w");
#else
    outfiles[j] = open_vector_current_file(fname);
#endif
    strcpy(fnames[j],fname);
    if(outfiles[j] == NULL){
      node0_printf("%s: Failed to open %s\n", myname, fname);
      exit(1);
    }
  }

}

//#if EIGMODE == EIGCG || EIGMODE == DEFLATE
#ifdef EIGMODE
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

/************************************ PUBLIC FUNCTIONS *****************************/
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
#define FINDVAR
void 
f_meas_current_diff( int n_masses, int const nrand, int nwrite, int thinning,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy,
		     su3_vector **eigVec, double *eigVal, int Nvecs, ks_param *ksp,
		     fermion_links_t *fl, Real *precs_sloppy, int n_slp,
		     char filenames[][MAXFILENAME],
		     int parity, int isCorr, int numMassPair, int isColorThin){
  
  char myname[] = "f_meas_current_multi_diff_eig";

  int i, is, j, jr, jrand, mu, ns, ir, thin_parity, shift, n_jmu=n_masses-numMassPair;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int ex, ey, ez, et, d = thinning;

  /* Block solver parameters -- temporary */
  int nr = 2;            /* Number of random sources to block */
  int evol = d*d*d*d*(isColorThin*3+(int) pow(0,isColorThin));
  int nsrc = evol*nr;
  su3_vector **gr;       /* Storage for sources */
  su3_vector **M_inv_gr; /* Storage for solutions */
  su3_vector **M_inv_gr0;

  /* Variables for computation of correlation */
  int charge[n_masses]; 
  Real *cov_mu[n_slp];     // covariance btwn slp and fine at each site for each component
  Real *js_mur[nr][n_slp]; // temp container for sloppy current
  Real *jsa_mu[n_slp];     // estimate of sloppy current density
  Real *js_mu_var[n_slp];  // variance of slp current over 4-volume and comp's for each residual
  Real *jp_mur[nr];        // temp container for fine current
  Real *jpa_mu;            // estimate of fine current density
  Real *jp_mu_var;         // variance of fine current over 4-volume and comp's for each residual
#if defined(FINDVAR) || defined(TIME_CURR)
  Real *jd_mur[nr][n_slp];
#endif
#ifdef FINDVAR
  Real jt_mu[nt*NMU];
  Real *jt2_mu[n_slp][n_jmu];
#endif
  
  /* Variables for computation of current densities at multiple precisions*/
  Real mass[n_masses];
#ifndef TIME_CURR
  Real *jd_mu[n_slp][n_jmu];
#endif
  imp_ferm_links_t *fn_multi[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

  /* Variables for files */
#ifdef TIME_CURR
  FILE **outfile = (FILE**) malloc(sizeof(FILE*)*n_slp*n_masses);
#else
  QIO_Writer * outfile[n_slp*n_jmu]; // = malloc( n_slp * n_masses * sizeof(QIO_Writer*) );
#endif
  char outfilenames[n_slp*n_jmu][MAXFILENAME];

  double wtime = 0.;

  /* Initialize fields for current densities, one for each mass */
  for(ns=0;ns<n_slp;ns++){
#ifndef TIME_CURR
    for(j = 0; j < n_jmu; j++)
      jd_mu[ns][j] = create_r_array_field(NMU);
#endif
    if(isCorr){
      cov_mu[ns]=create_r_array_field(NMU);
      for(ir=0;ir<nr;ir++) js_mur[ir][ns] = create_r_array_field(NMU);
      jsa_mu[ns]=create_r_array_field(NMU);
      js_mu_var[ns]=create_r_array_field(NMU);
    }
#if defined(FINDVAR) || defined(TIME_CURR)
    for(ir=0;ir<nr;ir++) jd_mur[ir][ns] = create_r_array_field(NMU);
#elif defined(FINDVAR)
    for(j = 0; j < n_jmu; j++){
      jt2_mu[ns][j] = (Real *) malloc(nt*NMU*sizeof(Real));
      for(i=0;i<nt*NMU;i++) jt2_mu[ns][j][i]=0;
    }
#endif  
  }
  if(isCorr){
    for(ir=0;ir<nr;ir++) jp_mur[ir] = create_r_array_field(NMU);
    jpa_mu=create_r_array_field(NMU);
    jp_mu_var=create_r_array_field(NMU);
  }

  /* Allocate source vectors for block inversion */
  gr  = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  M_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(is = 0; is < nsrc; is++){
    gr[is] = create_v_field();
    M_inv_gr[is] = create_v_field();
    if(numMassPair!=0) M_inv_gr0[is] = create_v_field();
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

  /* Open files for writing for multiple vcd files */ 
  open_files(myname, filenames, outfile, outfilenames, n_slp, precs_sloppy, n_masses, mass, numMassPair);

  node0_printf("Entered %s\n", myname); fflush(stdout);
 
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
      /* Iterate over displacements within a d^4 cube. Use displacements of specified parity */
      /* gr is indexed [(even sites,odd sites) x size of blocked random sources ] 
	 in the order of following for-loop expansion. */
      is = 0;
      for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
	if(parity==EVENANDODD || ex+ey+ez+et%2==parity%2)for(int color=0;color<isColorThin*3+(int) pow(0,isColorThin);color++){
		  
          thin_parity = (ex+ey+et+ez)%2==0?EVEN:ODD;
	
	  /* Apply source thinning */
	  copy_v_field(gr[is+evol/2*(thin_parity%2)+jr*evol], gr0);
	  thin_source_color( gr[is+evol/2*(thin_parity%2)+jr*evol], d, ex, ey, ez, et, color, isColorThin );
	  /* Project out the low mode part, based on the given eigenvectors */
	  project_out(gr[is+evol/2*(thin_parity%2)+jr*evol], eigVec, Nvecs, thin_parity);

	  if((ex+ey+ez+et)%2==0) is++;
	  if(is > nr*evol){
	    node0_printf("collect_evenodd_sources: Internal error: too many sources\n");
	    terminate(1);
	  }
      } /* ex, ey, ez, et, color */
    } /* jr */
    
    /* Compute the differenvce of sloppy and fine solve */
    /* M_inv_gr = M^{-1} gr (same random source for each mass) */
    /* Assume: The residual of inversion of mass difference is the same. */
    shift=0;
    for(j = 0; j < n_jmu; j++){
      
      /* Solve sloppily at multiple precisions */
      if(j<numMassPair){
	shift++;
	node0_printf("Solving sloppily for mass %g-%g\n", mass[j+shift-1],mass[j+shift]);
      }else
	node0_printf("Solving sloppily for mass %g\n", mass[j+shift]);
      
      for(is = 0; is < nsrc; is++){
	clear_v_field(M_inv_gr[is]);
	clear_v_field(M_inv_gr0[is]);
      }
      for(ns=0;ns < n_slp; ns++){
	/* invert */
	(qic_sloppy+j+shift)->resid = precs_sloppy[j+shift+n_masses*ns];
	node0_printf("Goal Sloppy Residual: %g\n",(qic_sloppy+j+shift)->resid);
	if(j<numMassPair){
	  mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr0, qic_sloppy + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);
	  mat_invert_block_uml_field_projected( nsrc/2, M_inv_gr0, M_inv_gr, qic_sloppy + j + shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
	}
	else
	  mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr, qic_sloppy + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);

	/* Apply the appropriate spin_taste operator for a nearly conserved current.  */
	for(is = 0; is < nsrc; is++){
	  for(mu = 0; mu < NMU; mu++){
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	    /* J_mu = imag[gr.M_inv_gr] */
	    /* SUBTRACT the sloppy result */
	    FORALLFIELDSITES(i){
	      complex cc = su3_dot( gr[is]+i, gr_mu+i );
	      if(j<numMassPair) cc.imag*=4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1]);
#ifndef TIME_CURR
	      jd_mu[ns][j][NMU*i + mu] -= cc.imag;
#endif
#if defined(FINDVAR) || defined(TIME_CURR)
	      jd_mur[is/evol][ns][NMU*i + mu] -= cc.imag;
#endif	      
	      if(isCorr)
		/* The factor of 2 comes from correcting for MILC normalization convention of M */
		js_mur[is/evol][ns][NMU*i + mu] += 2*charge[j+shift]*cc.imag/3;
	    }/* i:field sites*/
	  }/* mu */
	  if(j<numMassPair) copy_v_field(M_inv_gr[is],M_inv_gr0[is]);
	} /* is:nsrc  */
      }/* ns:n_slp */
	
      /* Next, continue to a fine solution from the same source */
      node0_printf("Solving precisely\n");
      if(j<numMassPair){
	mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr0, qic_precise + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);
	mat_invert_block_uml_field_projected( nsrc/2, M_inv_gr0, M_inv_gr, qic_precise + j+shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
      }
      else
	mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr, qic_precise + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);

      /* Apply the appropriate spin_taste operator for a nearly conserved current. */
      for(is = 0; is < nsrc; is++){
	for(mu = 0; mu < NMU; mu++){
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  /* ADD the precise result, which then gives the difference */
	  FORALLFIELDSITES(i){
	    complex cc = su3_dot( gr[is]+i, gr_mu+i );
	    if(j<numMassPair) cc.imag*=4*(mass[j+1]*mass[j+1]-mass[j]*mass[j]);
#ifndef TIME_CURR
	    for(ns=0;ns<n_slp;ns++) jd_mu[ns][j][NMU*i + mu] += cc.imag;
#endif
#if defined(FINDVAR) || defined(TIME_CURR)
	    for(ns=0;ns<n_slp;ns++) jd_mur[is/evol][ns][NMU*i + mu] += cc.imag;
#endif	      
	    if(isCorr)
	      /* The factor of 2 comes from correcting for MILC normalization convention of M */
	      jp_mur[is/evol][NMU*i + mu] += 2*charge[j+shift]*cc.imag/3;
	  }/* i:field sites*/
	} /* mu */
#ifdef FINDVAR
	if((is+1)%evol==0){
	  for(ns=0;ns<n_slp;ns++){
	    for(i=0;i<NMU*nt;i++) jt_mu[i]=0;
	    for(mu=0;mu<NMU;mu++){
	      FORALLFIELDSITES(i){
		jt_mu[NMU*lattice[i].t+mu]+=jd_mur[is/evol][ns][NMU*i+mu];
	      }
	    }
	    g_vecfloatsum(jt_mu,NMU*nt);
	    for(i=0;i<NMU*nt;i++) jt2_mu[ns][j][i]+=jt_mu[i]*jt_mu[i]/nwrite;
	  }
	}
#endif
#ifdef TIME_CURR
	if((is+1)%evol==0){
	  wtime -= dclock();
	  for(ns=0;ns<n_slp;ns++)
	    write_total_vector_current_record_t(jrand+is/evol, jd_mur[is/evol][ns], outfile, ns*n_jmu+j);
	  wtime += dclock();
	}
#endif
      } /* is: nrsc */
#if defined(FINDVAR) || defined(TIME_CURR)
      for(ir=0;ir<nr;ir++)for(ns=0;ns<n_slp;ns++) clear_r_array_field(jd_mur[ir][ns],NMU);
#endif
    } /* j: n_jmu */
    
    /* compute covariance and variance at each site */
    if(isCorr){
      for(ir=0;ir<nr;ir++){
	for(mu=0;mu<NMU;mu++){
	  //node0_printf("j(jrand=%d)_%d: ",jrand,mu);
	  FORALLFIELDSITES(i){
	    for(ns=0;ns<n_slp;ns++){
	      cov_mu[ns][NMU*i+mu] += js_mur[ir][ns][NMU*i + mu]*jp_mur[ir][NMU*i + mu];
	      jsa_mu[ns][NMU*i + mu] += js_mur[ir][ns][NMU*i + mu];
	      js_mu_var[ns][NMU*i + mu] += js_mur[ir][ns][NMU*i + mu]*js_mur[ir][ns][NMU*i + mu];
	      //if(i==node_index(0,0,0,0)) node0_printf("%.16e  ",js_mur[ir][ns][NMU*i + mu]);
	      js_mur[ir][ns][NMU*i + mu]=0;
	    }
	    jpa_mu[NMU*i + mu] += jp_mur[ir][NMU*i + mu];
	    jp_mu_var[NMU*i + mu] += jp_mur[ir][NMU*i + mu]*jp_mur[ir][NMU*i + mu];
	    //if(i==node_index(0,0,0,0)) node0_printf("%.16e at the origin\n",jp_mur[ir][NMU*i + mu]);
	    jp_mur[ir][NMU*i + mu]=0;
	  }
	}
      }
    } /* closed: isCorr */
    
    /* Write  record */
    /* No division by #thinned sublattices:
       This is b/c current density value at a site recieves a contribution only from a single sublattice.
       If we divide by the factor, the current density  will not be properly normalized.
    */
#ifndef TIME_CURR
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_jmu; j++){
	if(j%n_jmu==0) shift=0;
	if(j%n_jmu<numMassPair) shift++;
#ifdef FINDVAR 
	/* Assume: nwrite == nrand */
	int x,y,z,t;
	double msg[NMU];
	for(int k=0;k<nt*NMU;k++) jt_mu[k]=0;
	for(mu=0;mu<NMU;mu++){
	  FORALLFIELDSITES(i){
	    jt_mu[lattice[i].t*NMU+mu]+=jd_mu[j/n_jmu][j%n_jmu][NMU*i+mu];
	  }
	}
	g_vecfloatsum(jt_mu,nt*NMU);
	for(i=0;i<NMU*nt;i++){
	  jt2_mu[j/n_jmu][j%n_jmu][i]-=jt_mu[i]*jt_mu[i]/nwrite/nwrite;
	  node0_printf("jd_t(%d)_%d= %e %e for ns %d mass %d\n",
		       i/NMU,i%NMU,jt_mu[i]/nwrite,sqrt(jt2_mu[j/n_jmu][j%n_jmu][i]/nwrite),j/n_jmu,j%n_jmu);
	}
#endif
        int rep_mass = mass[j%n_jmu+shift];
        if(j%n_jmu<numMassPair){
          rep_mass =-mass[j%n_jmu]; // order switched to make the difference positive
	  node0_printf("For rand %d and mass %g-%g\n", jrand, mass[j%n_jmu+shift-1],mass[j%n_jmu+shift]);fflush(stdout);
	}
        else
	  node0_printf("For rand %d and mass %g\n", jrand,rep_mass);fflush(stdout);
	average_vector_current(nwrite, jd_mu[j/n_jmu][j%n_jmu]);
	// the metadata of mass is read but not used in rcorr
	int status = write_vector_current_record(outfile[j], jrand, nwrite, rep_mass, jd_mu[j/n_jmu][j%n_jmu]);
	clear_r_array_field(jd_mu[j/n_jmu][j%n_jmu], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, outfilenames[j]);
      } /* j: j_nmu */
      wtime += dclock();
    } /* if write */
#endif /* END: ifndef TIME_CURR */
  } /* jrand: nrand */

  /* Compute Statistics */
  if(isCorr){
    int pfactor = parity==EVENANDODD ? 1:2;
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
    /* get the correlation of the current densities at each & over all sites along with val & std of current at each site */
    Real js[n_slp],jp=0,r=0,r_mu[n_slp][NMU],js_mu_sig[n_slp],js_mu_siga[n_slp],jp_mu_sig=0,jp_mu_siga=0,ra[n_slp],rvar[n_slp];
    int x,y,z,t;
    for(ns=0;ns<n_slp;ns++){
      js[ns]=0;
      js_mu_sig[ns]=0;
      js_mu_siga[ns]=0;
      ra[ns]=0;
      rvar[ns]=0;
      for(mu=0;mu<NMU;mu++) r_mu[ns][mu]=0;
    }
    FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
      if(parity==EVENANDODD || ex+ey+ez+et%2==parity%2){
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
	    }
	  }
	  g_vecfloatsum(js,n_slp);
	  g_vecfloatsum(js_mu_sig,n_slp);
	  g_floatsum(&jp);
	  g_floatsum(&jp_mu_sig);
	  for(ns=0;ns<n_slp;ns++){
	    g_vecfloatsum(r_mu[ns],NMU);
	    /* Report */
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
    g_vecfloatsum(ra,n_slp);
    g_vecfloatsum(rvar,n_slp);
    g_vecfloatsum(js_mu_siga,n_slp);
    g_floatsum(&jp_mu_siga);

    /* Report Avg Statistics */
    for(ns=0;ns<n_slp;ns++){
      node0_printf("average of statistics on each sites for each components\n");
      node0_printf("avg_sigma(js,jp): %.16e %.16e\n",js_mu_siga[ns]/NMU/volume*pfactor,jp_mu_siga/NMU/volume*pfactor);
      node0_printf("r_avg: %lf pm %e\n",ra[ns]/NMU/volume*pfactor,sqrt((rvar[ns]-ra[ns]*ra[ns]/NMU/volume*pfactor)/NMU/volume*pfactor));
    }
  }

  node0_printf("\nTime to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);fflush(stdout);
  
  /* Cleanup */  
  for(j = 0; j < n_slp*n_jmu; j++){
#ifndef TIME_CURR
    close_vector_current_file(outfile[j]);
    destroy_r_array_field(jd_mu[j/n_jmu][j%n_jmu], NMU);
#else
    if(this_node==0) fclose(outfile[j]);
    free(outfile);outfile=NULL;
#endif
#if defined(FINDVAR) || defined(TIME_CURR)
    if(j%n_jmu==0)for(ir=0;ir<nr;ir++) destroy_r_array_field(jd_mur[ir][j/n_jmu], NMU);
    free(jt2_mu[j/n_jmu][j%n_jmu]);
  }
  jt2_mu[0][0]=NULL;
#else
  }
#endif
  if(isCorr){
    for(ns=0;ns<n_slp;ns++){
      destroy_r_array_field(cov_mu[ns],NMU);
      for(ir=0;ir<nr;ir++) destroy_r_array_field(js_mur[ir][ns],NMU);
      destroy_r_array_field(jsa_mu[ns],NMU);
      destroy_r_array_field(js_mu_var[ns],NMU);
    }
    for(ir=0;ir<nr;ir++) destroy_r_array_field(jp_mur[ir],NMU);
    destroy_r_array_field(jpa_mu,NMU);
    destroy_r_array_field(jp_mu_var,NMU);
  }

  for(is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr0[is]);
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(gr[is]);
  }
  destroy_v_field(gr_mu); gr_mu = NULL;
  free(gr); gr = NULL;
  free(M_inv_gr); M_inv_gr = NULL;
  free(M_inv_gr0); M_inv_gr0 = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_diff_eig: BLOCKCG version */
#undef FINDVAR

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Return: current densities at EVEN sites
*************************************************************************/
#define TESTEO // tests if the solving on both even and odd sites is working.
#define FINDVAR
void 
f_meas_current( int n_masses, int nrand, int nwrite, int thinning,
		quark_invert_control *qic,
		su3_vector **eigVec, double *eigVal, int Nvecs,
		ks_param *ksp, fermion_links_t *fl,
		Real *precs_sloppy, int n_slp,
		char filenames[][MAXFILENAME],
		int parity, int numMassPair, int isColorThin){
  
  char myname[] = "f_meas_current_multi_eig";

  int i, is, j, jr, n, jrand, mu, ns, p, shift, thin_parity, n_jmu = n_masses-numMassPair;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int ex, ey, ez, et, d = thinning;

  /* Block solver parameters -- temporary */
  int nr = 2;             /* Number of random sources to block */
  int evol = d*d*d*d*(isColorThin*3+(int) pow(0,isColorThin));
  int nsrc = evol*nr;
  su3_vector **gr;        /* Storage for sources */
  su3_vector **M_inv_gr;  /* Storage for solutions */
  su3_vector **M_inv_gr0; /* Storage for solutions */

  Real mass[n_masses];
  Real *j_mu[n_slp*n_jmu];
  Real *jlow_mu[n_jmu];
#if defined(FINDVAR) || defined(TESTEO)
  Real jt_mu[nt*NMU];
#endif
#ifdef FINDVAR
  Real *j_mur;
  Real *jt2_mu[n_slp*n_jmu];
#endif

  imp_ferm_links_t *fn_multi[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

#ifdef TIME_CURR
  FILE **outfile = (FILE**) malloc(sizeof(FILE*)*n_slp*n_jmu);
#else
  QIO_Writer *outfile[n_slp*n_jmu];
#endif
  char outfilenames[n_slp*n_masses][MAXFILENAME];

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
  for(j = 0; j < n_slp*n_jmu; j++){
    j_mu[j] = create_r_array_field(NMU);
    if(j<n_jmu) jlow_mu[j] = create_r_array_field(NMU);
#ifdef FINDVAR
    if(j==0) j_mur = create_r_array_field(NMU);
    jt2_mu[j] = (Real *) malloc(nt*NMU*sizeof(Real));
    for(i=0;i<nt*NMU;i++) jt2_mu[j][i]=0;
#endif
  }

  /* Allocate source vectors */
  gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  M_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(is = 0; is < nsrc; is++){
    gr[is] = create_v_field();
    M_inv_gr[is] = create_v_field();
    if(numMassPair!=0) M_inv_gr0[is] = create_v_field();
  }

  for(j = 0; j < n_masses; j++){
    /* Load masses from ks_param */
    mass[j] = ksp[j].mass;
    /* Load pointers for fermion links, based on Naik epsilon indices */
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  }
  
  /* Open files for writing */
  open_files(myname, filenames, outfile, outfilenames, n_slp, precs_sloppy, n_masses, mass, numMassPair);

  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    shift = 0;
    for(j = 0; j < n_jmu; j++){
      if(j<numMassPair) shift++;
      dslash_fn_field(eigVec[n], gr0, ODD, fn_multi[j+shift]);
      for(int p=1;p<3;p++){
	for(mu = 0; mu < NMU; mu++){
	  if(p == EVEN) spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, gr0);
	  else          spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, eigVec[n]);
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  FORSOMEFIELDPARITY(i,p){
	    complex z;
	    if(p == EVEN) z = su3_dot( eigVec[n] + i, gr_mu + i);
	    else          z = su3_dot( gr0 + i, gr_mu + i);
	    if(j<numMassPair) 
	      jlow_mu[j][NMU*i + mu] += 4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1])*pow(-1,p-1)*z.imag/(eigVal[n]+4.0*mass[j+shift-1]*mass[j+shift-1])/(eigVal[n]+4.0*mass[j+shift]*mass[j+shift]);
	    else
	      jlow_mu[j][NMU*i + mu] += pow(-1,p-1)*z.imag/(eigVal[n]+4.0*mass[j+shift]*mass[j+shift]);
	  } /* i */
	} /* mu */
      } /* p: ODD(1) and EVEN(2) */
    } /* j */
  } /* n */

#ifdef TESTEO
  for(j = 0; j < n_jmu; j++){
    node0_printf("For mass %g\n", mass[j]);
    for(mu = 0; mu < NMU; mu++){
      FORALLFIELDSITES(i){
	node0_printf("j_mu_low  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, jlow_mu[j][NMU*i+mu]);
      }
    }
  }
#endif
#ifdef TESTEO // assume: the current density under consideraton is 2+1+1
  for(int t=0;t<nt*NMU;t++) jt_mu[t]=0;
  int charge[3] = {1,-1,2};
  FORALLFIELDSITES(i){
    for(j=0;j<n_masses;j++) 
      for(mu=0;mu<NMU;mu++) jt_mu[lattice[i].t*NMU+mu]+=charge[j]*jlow_mu[j][NMU*i+mu];
  }
  g_vecfloatsum(jt_mu,nt*NMU);
  for(i=0;i<NMU*nt;i++) node0_printf("j_t_low: %d %d %e\n",i/NMU,i%NMU,jt_mu[i]);
#endif

  dtime += dclock();
  node0_printf("Time for exact low modes %g sec\n", dtime);

  /* HACK to get only result from low modes  */
  if(nrand == 0){
    for(j = 0; j < n_slp*n_jmu; j++){
#ifndef TIME_CURR
      average_vector_current_and_sum(1, j_mu[j], jlow_mu[j%n_jmu]);
      int status = write_vector_current_record(outfile[j], 0, 1, mass[j%n_jmu], j_mu[j]);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, outfilenames[j]);
      } 
#else
      wtime -= dclock();
      write_total_vector_current_record_t(0, jlow_mu[j%n_jmu], outfile, j);
      wtime += dclock();
#endif
      clear_r_array_field(j_mu[j], NMU);
    }
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
      /* Iterate over displacements within a d^4 cube. Use displacements of specified parity */
      /* gr is indexed [(even sites,odd sites) x size of blocked random sources ] 
	 in the order of following for-loop expansion. */
      is = 0;
      for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
	if(parity==EVENANDODD || ex+ey+ez+et%2==parity%2)for(int color=0;color<isColorThin*3+(int) pow(0,isColorThin);color++){
		  
          thin_parity = (ex+ey+et+ez)%2==0?EVEN:ODD;
	
	  /* Apply source thinning */
	  copy_v_field(gr[is+evol/2*(thin_parity%2)+jr*evol], gr0);
	  thin_source_color( gr[is+evol/2*(thin_parity%2)+jr*evol], d, ex, ey, ez, et, color, isColorThin );
	  /* Project out the low mode part, based on the given eigenvectors */
	  project_out(gr[is+evol/2*(thin_parity%2)+jr*evol], eigVec, Nvecs, thin_parity);

	  if((ex+ey+ez+et)%2==0) is++;
	  if(is > nr*evol){
	    node0_printf("collect_evenodd_sources: Internal error: too many sources\n");
	    terminate(1);
	  }
#if 1 	  /* DEBUG */
	  /* Check the norm of the reduced source */
	  double_complex dd;
	  dot_product(gr[is], gr[is], &dd, thin_parity);
	  node0_printf("Deflated source norm %d %d %d %d %g\n", ex, ey, ez, et, dd.real);
#endif
      } /* ex, ey, ez, et, color */
    } /* jr */
    
    /* Start inversion */
    /* M_inv_gr = M^{-1} gr (same random source for each mass) */
    /* Assume: The residual of inversion of mass difference is the same. */
    shift=0;
    for(j = 0; j < n_jmu; j++){
      if(j<numMassPair){
	shift++;      
	node0_printf("Solving for mass %g-%g\n", mass[j+shift-1],mass[j+shift]);
      }else
	node0_printf("Solving for mass %g\n", mass[j+shift]);

      for(is = 0; is < nsrc; is++){
	clear_v_field(M_inv_gr[is]);
	clear_v_field(M_inv_gr0[is]);
      }
      for(ns=0;ns<n_slp;ns++){
	(qic+j+shift)->resid = precs_sloppy[j+shift+ns*n_masses];
	node0_printf("Goal Sloppy Residual: %g\n",(qic+j+shift)->resid);
	if(j<numMassPair){
	  mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr0, qic + j + shift, thin_parity, mass[j+1], fn_multi[j+shift]);
	  mat_invert_block_uml_field_projected( nsrc/2, M_inv_gr0, M_inv_gr, qic + j + shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
	}
	else
	  mat_invert_block_uml_field_projected( nsrc/2, gr, M_inv_gr, qic + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);

	for(is = 0; is < nsrc; is++){
	  /* Apply current in various directions at the sink */
	  for(mu = 0; mu < NMU; mu++){
	    /* Apply the appropriate spin_taste operator for
	       a nearly conserved current. */
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	    /* J_mu = imag[gr.M_inv_gr] */
	    FORALLFIELDSITES(i){
	      complex cc = su3_dot( gr[is]+i, gr_mu+i );
	      if(j<numMassPair) cc.imag*=4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1]);
	      j_mu[n_jmu*ns+j][NMU*i + mu] += cc.imag;
#ifdef TESTEO
	      if(cc.real!=0 || cc.imag!=0) printf("j_mu_%d %d %d %d %d %d %g %g\n",ns, lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
#endif
#ifdef FINDVAR
	      j_mur[NMU*i + mu] += cc.imag;
#endif
	    } /* i: FORALLFIELDSITES*/
	  } /* mu: NMU */
	  if(j<numMassPair) copy_v_field(M_inv_gr[is],M_inv_gr0[is]);
#ifdef FINDVAR
	  /* Compute the variance of the total current density on each time slice */
	  if((is+1)%evol==0){
	    for(i=0;i<NMU*nt;i++) jt_mu[i]=0;
	    for(mu=0;mu<NMU;mu++){
	      FORALLFIELDSITES(i){
		jt_mu[NMU*lattice[i].t+mu]+=jlow_mu[j][NMU*i + mu]+j_mur[NMU*i+mu];
	      }
	    }
	    g_vecfloatsum(jt_mu,NMU*nt);
	    for(i=0;i<NMU*nt;i++) jt2_mu[n_jmu*ns+j][i]+=jt_mu[i]*jt_mu[i]/nwrite;
	    clear_r_array_field(j_mur,NMU);
	  } 
#endif
#ifdef TIME_CURR
	  if((is+1)%evol==0){
	    wtime -= dclock();
	    for(int k=0;k<nt*NMU;k++) jt_mu[k]=0;
	    FORALLFIELDSITES(i){
	      j_mu[ns*n_jmu+j][NMU*i+mu]+=jlow_mu[j][NMU*i + mu];
	    }
	    write_total_vector_current_record_t(jrand+is/evol, j_mu[ns*n_jmu+j], outfile, ns*n_jmu+j);
	    clear_r_array_field(j_mu[ns*n_jmu+j],NMU);
	    wtime += dclock();
	  }
#endif
	} /* is: nsrc */
      } /* ns: n_slp */
    } /* j: n_jmu */

#ifndef TIME_CURR
    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_jmu; j++){
	if(j%n_jmu==0) shift=0;
	if(j%n_jmu<numMassPair) shift++;
#ifdef FINDVAR
	int x,y,z,t;
	double msg[NMU];
	for(int k=0;k<nt*NMU;k++) jt_mu[k]=0;
	for(mu=0;mu<NMU;mu++){
	  FORALLFIELDSITES(i){
	    jt_mu[lattice[i].t*NMU+mu]+=nwrite*jlow_mu[j%n_jmu][NMU*i + mu]+j_mu[j][NMU*i+mu];
	  }
	}
	g_vecfloatsum(jt_mu,nt*NMU);
	for(i=0;i<NMU*nt;i++){
	  jt2_mu[j][i]-=jt_mu[i]*jt_mu[i]/nwrite/nwrite;
	  node0_printf("j_t(%d)_%d= %e %e for ns %d mass %d\n",i/NMU,i%NMU,jt_mu[i]/nwrite,sqrt(jt2_mu[j][i]/nwrite),j/n_jmu,j%n_jmu);
	}
	// print j_mu at all sites
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  for(mu=0;mu<NMU;mu++) msg[mu]=0;
	  if(this_node==node_number(x,y,z,t))for(mu=0;mu<NMU;mu++) msg[mu]=j_mu[j][NMU*node_index(x,y,z,t)+mu];
	  g_vecfloatsum(msg,NMU);
	  for(mu=0;mu<NMU;mu++) 
	    node0_printf("j(%d,%d,%d,%d)_%d= %.16e for mass%d and ns=%d\n",x,y,z,t,mu,msg[mu]/nwrite,j%n_jmu,j/n_jmu);
	}
#endif
#if 0	/* DEBUG */
	FOREVENFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
        int rep_mass = mass[j%n_jmu+shift];
        if(j%n_jmu<numMassPair){
          rep_mass =-mass[j%n_jmu]; // order switched to make the difference positive
          node0_printf("For rand %d and mass %g-%g\n", jrand, mass[j%n_jmu+shift-1],mass[j%n_jmu+shift]);fflush(stdout);
        }
        else
          node0_printf("For rand %d and mass %g\n", jrand,rep_mass);fflush(stdout);
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j%n_jmu]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, rep_mass, j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j: n_slp*n_jmu */
      wtime += dclock();
    } /* if write */
#endif /* END: ifndef TIME_CURR */
  } /* jrand */

  /* Clean up */
  for(j = 0; j < n_slp*n_jmu; j++){
#ifndef TIME_CURR
    close_vector_current_file(outfile[j]);
#else
    if(this_node==0) fclose(outfile[j]);
    free(outfile);outfile=NULL;
#endif
    destroy_r_array_field(j_mu[j], NMU);
    if(j/n_masses==0) destroy_r_array_field(jlow_mu[j], NMU);
#ifdef FINDVAR
    if(j==0) destroy_r_array_field(j_mur,NMU);
    free(jt2_mu[j]);
  }
   jt2_mu[0]=NULL;
#else
  }
#endif

  node0_printf("Time to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);

  for(is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(M_inv_gr0[is]);
    destroy_v_field(gr[is]);
  }

  free(M_inv_gr0); M_inv_gr0 = NULL;
  free(M_inv_gr); M_inv_gr = NULL;
  free(gr); gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig: BLOCKCG version */
#undef FINDVAR
#undef TESTEO

#else //EIGMODE BUT NOT BLOCKCG

/************************************************************************
 * Entry point for multiple masses with deflation and iterated single-mass inverter.
 * This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities with multiple sloppy precisions on sites with the spcified parity.
 * It also computes correlation of sloppy and precise current densities for
   multiple sloppy precisions.  
 * If isCorr != 0, the correlation between sloppy solve at multiple precision and 
   fine solve is evaluated.
 * If numMassPair != 0, the difference of the two current densities for a pair of
   masses (m1,m2) is estimated using the mass-difference trick
   * The difference is J_m2-J_m1
   * In this case, the first numMassPair entries of ksp are assumed to be the 
     pairs for which the trick is applied.
 * If isColorThin == 3, the color thinning is performed.
 * Designed for use with deflation or eigcg. 
 * Requires a set of accurate low-mode eigenpairs
 * Assume: parameter values in ksp for each pair are the same
 * Return: current densities at sites of the specified parity
 ************************************************************************/
#define FINDVAR
void 
f_meas_current_diff( int n_masses, int const nrand, int nwrite, int thinning,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy, 
		     su3_vector **eigVec, double *eigVal, int Nvecs, ks_param *ksp, 
		     fermion_links_t *fl, Real *precs_sloppy, int n_slp,
		     char filenames[][MAXFILENAME],
		     int parity, int isCorr, int numMassPair, int isColorThin){

  char myname[] = "f_meas_current_multi_diff";

  int i, j, jrand, mu, ns, thin_parity, shift, n_jmu=n_masses-numMassPair;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr0 = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  /* Variables for computation of correlation */
  int charge[n_masses]; 
  Real *cov_mu[n_slp];     // covariance btwn slp and fine at each site for each component
  Real *js_mu[n_slp];      // temp container for sloppy current
  Real *jsa_mu[n_slp];     // estimate of sloppy current density
  Real *js_mu_var[n_slp];  // variance of slp current over 4-volume and comp's for each residual
  Real *jp_mu;             // temp container for fine current
  Real *jpa_mu;            // estimate of fine current density
  Real *jp_mu_var;         // variance of fine current over 4-volume and comp's for each residual
#ifdef FINDVAR
  Real jt_mu[nt*NMU];
  Real *jd_mur[n_slp][n_jmu];
  Real *jt2_mu[n_slp][n_jmu];
#endif
  
  /* Variables for computation of current densities at multiple precisions*/
  Real mass[n_masses];
  Real *jd_mu[n_slp][n_jmu];
  imp_ferm_links_t *fn_multi[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

  /* Variables for files */
#ifdef TIME_CURR
  FILE **outfile = (FILE**) malloc(sizeof(FILE*)*n_slp*n_jmu);
#else
  QIO_Writer * outfile[n_slp*n_jmu]; // = malloc( n_slp * n_masses * sizeof(QIO_Writer*) );
#endif
  char outfilenames[n_slp*n_jmu][MAXFILENAME];

  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Initialize fields for current densities, one for each mass */
  for(ns=0;ns<n_slp;ns++){
    for(j = 0; j < n_jmu; j++)
      jd_mu[ns][j] = create_r_array_field(NMU);
    if(isCorr){
      cov_mu[ns]=create_r_array_field(NMU);
      js_mu[ns] = create_r_array_field(NMU);
      jsa_mu[ns]=create_r_array_field(NMU);
      js_mu_var[ns]=create_r_array_field(NMU);
    }
#ifdef FINDVAR
    for(j = 0; j < n_jmu; j++){
      jd_mur[ns][j] = create_r_array_field(NMU);
      jt2_mu[ns][j] = (Real *) malloc(nt*NMU*sizeof(Real));
      for(i=0;i<nt*NMU;i++) jt2_mu[ns][j][i]=0;
    }
#endif  
  }
  if(isCorr){
    jp_mu = create_r_array_field(NMU);
    jpa_mu=create_r_array_field(NMU);
    jp_mu_var=create_r_array_field(NMU);
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

  /* Open files for writing for multiple vcd files */ 
  open_files(myname, filenames, outfile, outfilenames, n_slp, precs_sloppy, n_masses, mass, numMassPair);

  double wtime = 0.0;

  /* Loop over random sources */
  for(jrand = 0; jrand < nrand; jrand++){
     /* Make random source, and do inversion */
#ifndef Z2RSOURCE
    grsource_plain_field( gr0, EVENANDODD );
#else
    z2rsource_plain_field( gr0, EVENANDODD );
#endif
    /* Iterate over displacements within a d^4 cube. Use displacements of specified parity */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if(parity==EVENANDODD || ex+ey+ez+et%2==parity%2)for(int color=0;color<isColorThin*3+(int) pow(0,isColorThin);color++){
      thin_parity = (ex+ey+et+ez)%2==0?EVEN:ODD;
 	      
      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;
      
      /* Apply source thinning */
      copy_v_field(gr, gr0);
      thin_source_color( gr, d, ex, ey, ez, et, color, isColorThin );
      
      /* Project out the low mode part, based on the given eigenvectors */
      project_out(gr, eigVec, Nvecs, thin_parity);
      
      /* Compute the differenvce of sloppy and fine solve */
      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
      /* Assume: The residual of inversion of mass difference is the same. */
      shift=0;
      for(j = 0; j < n_jmu; j++){

	/* Solve sloppily at multiple precisions */
	if(j<numMassPair){
	  shift++;
	  node0_printf("Solving sloppily for %d %d %d %d mass %g-%g\n", ex, ey, ez, et, mass[j+shift-1],mass[j+shift]);
        }else
          node0_printf("Solving sloppily for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j+shift]);

	clear_v_field(M_inv_gr);
	clear_v_field(M_inv_gr0);
	for(ns=0;ns < n_slp; ns++){
	  /* invert */
	  (qic_sloppy+j+shift)->resid = precs_sloppy[j+shift+n_masses*ns];
	  node0_printf("Goal Sloppy Residual: %g\n",(qic_sloppy+j+shift)->resid);
	  if(j<numMassPair){
	    mat_invert_uml_field_projected( gr, M_inv_gr0, qic_sloppy + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);
	    mat_invert_uml_field_projected( M_inv_gr0, M_inv_gr, qic_sloppy + j + shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
	  }
	  else
	    mat_invert_uml_field_projected( gr, M_inv_gr, qic_sloppy + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);

	  /* Apply the appropriate spin_taste operator for a nearly conserved current.  */
	  for(mu = 0; mu < NMU; mu++){
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	    /* J_mu = imag[gr.M_inv_gr] */
	    /* SUBTRACT the sloppy result */
	    FORALLFIELDSITES(i){
	      complex cc = su3_dot( gr+i, gr_mu+i );
	      if(j<numMassPair) cc.imag*=4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1]);
	      jd_mu[ns][j][NMU*i + mu] -= cc.imag;
#ifdef FINDVAR
	      jd_mur[ns][j][NMU*i + mu] -= cc.imag;
#endif	      
	      if(isCorr)
		/* The factor of 2 comes from correcting for MILC normalization convention of M */
		js_mu[ns][NMU*i + mu] += 2*charge[j+shift]*cc.imag/3;
	    }/* i:field sites*/
	  }/* mu */
	  if(j<numMassPair) copy_v_field(M_inv_gr,M_inv_gr0);
	}/* ns:n_slp */
	
	/* Next, continue to a fine solution from the same source */
	node0_printf("Solving precisely for %d %d %d %d\n", ex, ey, ez, et);
	if(j<numMassPair){
	  mat_invert_uml_field_projected( gr, M_inv_gr0, qic_precise + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);
	  mat_invert_uml_field_projected( M_inv_gr0, M_inv_gr, qic_precise + j+shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
	}
	else
	  mat_invert_uml_field_projected( gr, M_inv_gr, qic_precise + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);

	/* Apply the appropriate spin_taste operator for a nearly conserved current. */
	for(mu = 0; mu < NMU; mu++){
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  /* J_mu = imag[gr.M_inv_gr] */
	  /* ADD the precise result, which then gives the difference */
	  FORALLFIELDSITES(i){
	    complex cc = su3_dot( gr+i, gr_mu+i );
	    if(j<numMassPair) cc.imag*=4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1]);
	    for(ns=0;ns<n_slp;ns++) jd_mu[ns][j][NMU*i + mu] += cc.imag;
#ifdef FINDVAR
	    for(ns=0;ns<n_slp;ns++) jd_mur[ns][j][NMU*i + mu] += cc.imag;
#endif	      
	    if(isCorr)
	      /* The factor of 2 comes from correcting for MILC normalization convention of M */
	      jp_mu[NMU*i + mu] += 2*charge[j+shift]*cc.imag/3;
	  }/* i:field sites*/
	} /* mu */
      } /* j: n_jmu */
    } /* ex, ey, ez, et, color */

    /* compute covariance and variance at each site */
    if(isCorr){
      for(mu=0;mu<NMU;mu++){
	//node0_printf("j(jrand=%d)_%d: ",jrand,mu);
	FORALLFIELDSITES(i){
	  for(ns=0;ns<n_slp;ns++){
	    cov_mu[ns][NMU*i+mu] += js_mu[ns][NMU*i + mu]*jp_mu[NMU*i + mu];
	    jsa_mu[ns][NMU*i + mu] += js_mu[ns][NMU*i + mu];
	    js_mu_var[ns][NMU*i + mu] += js_mu[ns][NMU*i + mu]*js_mu[ns][NMU*i + mu];
	    //if(i==node_index(0,0,0,0)) node0_printf("%.16e  ",js_mu[ns][NMU*i + mu]);
	    js_mu[ns][NMU*i + mu]=0;
	  }
	  jpa_mu[NMU*i + mu] += jp_mu[NMU*i + mu];
	  jp_mu_var[NMU*i + mu] += jp_mu[NMU*i + mu]*jp_mu[NMU*i + mu];
	  //if(i==node_index(0,0,0,0)) node0_printf("%.16e at the origin\n",jp_mu[NMU*i + mu]);
	  jp_mu[NMU*i + mu]=0;
	}
      }
    } /* closed: isCorr */
    
#ifdef FINDVAR
    for(j=0;j<n_slp*n_jmu;j++){
      for(i=0;i<NMU*nt;i++) jt_mu[i]=0;
      for(mu=0;mu<NMU;mu++){
	FORALLFIELDSITES(i){
	  jt_mu[NMU*lattice[i].t+mu]+=jd_mur[j/n_jmu][j%n_jmu][NMU*i+mu];
	}
      }
      g_vecfloatsum(jt_mu,NMU*nt);
      for(i=0;i<NMU*nt;i++) jt2_mu[j/n_jmu][j%n_jmu][i]+=jt_mu[i]*jt_mu[i]/nwrite;
      clear_r_array_field(jd_mur[j/n_jmu][j%n_jmu],NMU);
    }
#endif
#ifdef TIME_CURR
    wtime -= dclock();
    for(j=0;j<n_slp*n_jmu;j++)
      write_total_vector_current_record_t(jrand, jd_mur[j/n_jmu][j%n_jmu], outfile, j);
    clear_r_array_field(jd_mu[j/n_jmu][j%n_jmu],NMU);
    wtime += dclock();
#endif

#ifndef TIME_CURR
    /* Write record */
    /* No division by #thinned sublattices:
       This is b/c current density value at a site recieves a contribution only from a single sublattice.
       If we divide by the factor, the current density  will not be properly normalized.
    */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_jmu; j++){
	if(j%n_jmu==0) shift=0;
	if(j%n_jmu<numMassPair) shift++;
#ifdef FINDVAR
	int x,y,z,t;
	double msg[NMU];
	for(int k=0;k<nt*NMU;k++) jt_mu[k]=0;
	for(mu=0;mu<NMU;mu++){
	  FORALLFIELDSITES(i){
	    jt_mu[lattice[i].t*NMU+mu]+=jd_mu[j/n_jmu][j%n_jmu][NMU*i+mu];
	  }
	}
	g_vecfloatsum(jt_mu,nt*NMU);
	for(i=0;i<NMU*nt;i++){
	  jt2_mu[j/n_jmu][j%n_jmu][i]-=jt_mu[i]*jt_mu[i]/nwrite/nwrite;
	  node0_printf("jd_t(%d)_%d= %e %e for ns %d mass %d\n",
		       i/NMU,i%NMU,jt_mu[i]/nwrite,sqrt(jt2_mu[j/n_jmu][j%n_jmu][i]/nwrite),j/n_jmu,j%n_jmu);
	}
#endif
        int rep_mass = mass[j%n_jmu+shift];
        if(j%n_jmu<numMassPair){
          rep_mass =-mass[j%n_jmu]; // order switched to make the difference positive
	  node0_printf("For rand %d and mass %g-%g\n", jrand, mass[j%n_jmu+shift-1],mass[j%n_jmu+shift]);fflush(stdout);
	}
        else
	  node0_printf("For rand %d and mass %g\n", jrand,rep_mass);fflush(stdout);
	average_vector_current(nwrite, jd_mu[j/n_jmu][j%n_jmu]);
	// the metadata of mass is read but not used in rcorr
	int status = write_vector_current_record(outfile[j], jrand, nwrite, rep_mass, jd_mu[j/n_jmu][j%n_jmu]);
	clear_r_array_field(jd_mu[j/n_jmu][j%n_jmu], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, outfilenames[j]);
      } /* j: j_nmu */
      wtime += dclock();
    } /* if write */
#endif /* END: TIME_CURR */
  } /* jrand: nrand */
    
  /* Compute Statistics */
  if(isCorr){
    int pfactor = parity==EVENANDODD ? 1:2;
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
    /* get the correlation of the current densities at each & over all sites along with val & std of current at each site */
    Real js[n_slp],jp=0,r=0,r_mu[n_slp][NMU],js_mu_sig[n_slp],js_mu_siga[n_slp],jp_mu_sig=0,jp_mu_siga=0,ra[n_slp],rvar[n_slp];
    int x,y,z,t;
    for(ns=0;ns<n_slp;ns++){
      js[ns]=0;
      js_mu_sig[ns]=0;
      js_mu_siga[ns]=0;
      ra[ns]=0;
      rvar[ns]=0;
      for(mu=0;mu<NMU;mu++) r_mu[ns][mu]=0;
    }
    FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
      if(parity==EVENANDODD || ex+ey+ez+et%2==parity%2){
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
	    }
	  }
	  g_vecfloatsum(js,n_slp);
	  g_vecfloatsum(js_mu_sig,n_slp);
	  g_floatsum(&jp);
	  g_floatsum(&jp_mu_sig);
	  for(ns=0;ns<n_slp;ns++){
	    g_vecfloatsum(r_mu[ns],NMU);
	    /* Report */
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
    g_vecfloatsum(ra,n_slp);
    g_vecfloatsum(rvar,n_slp);
    g_vecfloatsum(js_mu_siga,n_slp);
    g_floatsum(&jp_mu_siga);

    /* Report Avg Statistics */
    for(ns=0;ns<n_slp;ns++){
      node0_printf("average of statistics on each sites for each components\n");
      node0_printf("avg_sigma(js,jp): %.16e %.16e\n",js_mu_siga[ns]/NMU/volume*pfactor,jp_mu_siga/NMU/volume*pfactor);
      node0_printf("r_avg: %lf pm %e\n",ra[ns]/NMU/volume*pfactor,sqrt((rvar[ns]-ra[ns]*ra[ns]/NMU/volume*pfactor)/NMU/volume*pfactor));
    }
  }

  node0_printf("\nTime to write %d records for %d masses = %e sec\n", nrand/nwrite, n_masses, wtime);fflush(stdout);
  
  /* Cleanup */  
  for(j = 0; j < n_slp*n_jmu; j++){
#ifndef TIME_CURR
    close_vector_current_file(outfile[j]);
#else
    if(this_node==0) fclose(outfile[j]);
    free(outfile);outfile=NULL;
#endif
    destroy_r_array_field(jd_mu[j/n_jmu][j%n_jmu], NMU);
#ifdef FINDVAR
    destroy_r_array_field(jd_mur[j/n_jmu][j%n_jmu], NMU);
    free(jt2_mu[j/n_jmu][j%n_jmu]);
  }
  jt2_mu[0][0]=NULL;
#else
  }
#endif
  if(isCorr){
    for(ns=0;ns<n_slp;ns++){
      destroy_r_array_field(cov_mu[ns],NMU);
      destroy_r_array_field(js_mu[ns],NMU);
      destroy_r_array_field(jsa_mu[ns],NMU);
      destroy_r_array_field(js_mu_var[ns],NMU);
    }
    destroy_r_array_field(jp_mu,NMU);
    destroy_r_array_field(jpa_mu,NMU);
    destroy_r_array_field(jp_mu_var,NMU);
  }

  destroy_v_field(M_inv_gr0); M_inv_gr0 = NULL;
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_diff_eig_corr: iterated single-mass inverter */
#undef FINDVAR

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * It solves for the same source at multiple precisions.
 * If numMassPair != 0, the difference of the two current densities for a pair of
   masses (m1,m2) is estimated using the mass-difference trick
   * The difference is J_m2-J_m1
   * In this case, the first numMassPair entries of ksp are assumed to be the 
     pairs for which the trick is applied.
 * If isColorThin == 3, the color thinning is performed.
 * Designed for use with deflation or eigcg.
 * Does deflation, so requires a set of accurate low-mode eigenpairs
 * Assume: the parameter values for paired masses are identical
 * Return: current density at sites of the specified parity
 ************************************************************************/
//#define TESTEO // tests if the solving on both even and odd sites is working.
//#define PT2PT
#define FINDVAR
//#define FINDCOV
void 
f_meas_current( int n_masses, int nrand, int nwrite, int thinning,
		quark_invert_control *qic,
		su3_vector **eigVec, double *eigVal, int Nvecs,
		ks_param *ksp, fermion_links_t *fl, 
		Real *precs_sloppy, int n_slp,
		char filenames[][MAXFILENAME],
		int parity, int numMassPair, int isColorThin){
  
  char myname[] = "f_meas_current_multi_eig_eo";

  int i, j, n, jrand, mu, ns, p, shift, thin_parity, n_jmu = n_masses-numMassPair;

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *gr0 = create_v_field();
  su3_vector *gr = create_v_field();
  su3_vector *gr_mu = create_v_field();
  su3_vector *M_inv_gr0 = create_v_field();
  su3_vector *M_inv_gr = create_v_field();
  int ex, ey, ez, et, d = thinning;

  Real mass[n_masses];
  Real *j_mu[n_slp*n_jmu];
  Real *jlow_mu[n_jmu];
#if defined(FINDVAR) || defined(FINDCOV) || defined(TESTEO)
  Real jt_mu[nt*NMU];
#endif
#ifdef FINDVAR
  Real *j_mur[n_slp*n_jmu];
  Real *jt2_mu[n_slp*n_jmu];
#endif
#ifdef FINDCOV
  //Real *cov[n_masses][NMU][sites_on_node];
  int svolume=nx*ny*nz;
  Real cov[n_slp*n_jmu][NMU][nt][svolume][svolume];
#endif

  imp_ferm_links_t *fn_multi[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);

#ifdef TIME_CURR
  FILE **outfile = (FILE**) malloc(sizeof(FILE*)*n_slp*n_jmu);
#else
  QIO_Writer *outfile[n_slp*n_jmu];
#endif
  char outfilenames[n_slp*n_masses][MAXFILENAME];

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
  for(j = 0; j < n_slp*n_jmu; j++){
    j_mu[j] = create_r_array_field(NMU);
    if(j<n_jmu) jlow_mu[j] = create_r_array_field(NMU);
#if defined(FINDVAR) || defined(FINDCOV) || defined(TIME_CURR)
    j_mur[j] = create_r_array_field(NMU);
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

  /* Load masses from ks_param */
  for(j = 0; j < n_masses; j++)
    mass[j] = ksp[j].mass;
  
  /* Load pointers for fermion links, based on Naik epsilon indices */
  for(j = 0; j < n_masses; j++)
    fn_multi[j] = fn[ksp[j].naik_term_epsilon_index];
  
  /* Open files for writing */
  open_files(myname, filenames, outfile, outfilenames, n_slp, precs_sloppy, n_masses, mass, numMassPair);

  /* Compute exact low-mode current density */
  double dtime = -dclock();
  for(n = 0; n < Nvecs; n++){
    shift = 0;
    for(j = 0; j < n_jmu; j++){
      if(j<numMassPair) shift++;
      dslash_fn_field(eigVec[n], gr, ODD, fn_multi[j+shift]);
      for(int p=1;p<3;p++){
	for(mu = 0; mu < NMU; mu++){
	  if(p == EVEN) spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, gr);
	  else          spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, eigVec[n]);
	  spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	  FORSOMEFIELDPARITY(i,p){
	    complex z;
	    if(p == EVEN) z = su3_dot( eigVec[n] + i, gr_mu + i);
	    else          z = su3_dot( gr + i, gr_mu + i);
	    if(j<numMassPair) 
	      jlow_mu[j][NMU*i + mu] += 4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1])*pow(-1,p-1)*z.imag/(eigVal[n]+4.0*mass[j+shift-1]*mass[j+shift-1])/(eigVal[n]+4.0*mass[j+shift]*mass[j+shift]);
	    else
	      jlow_mu[j][NMU*i + mu] += pow(-1,p-1)*z.imag/(eigVal[n]+4.0*mass[j+shift]*mass[j+shift]);
	  } /* i */
	} /* mu */
      } /* p: ODD(1) and EVEN(2) */
    } /* j */
  } /* n */

#ifdef TESTEO
  for(j = 0; j < n_jmu; j++){
    node0_printf("For mass %g\n", mass[j]);
    for(mu = 0; mu < NMU; mu++){
      FORALLFIELDSITES(i){
	node0_printf("j_mu_low  %d %d %d %d %d %g\n",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, jlow_mu[j][NMU*i+mu]);
      }
    }
  }
#endif
#ifdef TESTEO // assume: the current density under consideraton is 2+1+1
  Real j_t[nt*NMU];
  for(int t=0;t<nt*NMU;t++) j_t[t]=0;
  int charge[3] = {1,-1,2};
  FORALLFIELDSITES(i){
    for(j=0;j<n_masses;j++) 
      for(mu=0;mu<NMU;mu++) j_t[lattice[i].t*NMU+mu]+=charge[j]*jlow_mu[j][NMU*i+mu];
  }
  g_vecfloatsum(j_t,nt*NMU);
  for(i=0;i<NMU*nt;i++) node0_printf("j_t_low: %d %d %e\n",i/NMU,i%NMU,j_t[i]);
#endif

  dtime += dclock();
  node0_printf("Time for exact low modes %g sec\n", dtime);

  /* HACK to get only result from low modes  */
  if(nrand == 0){
    for(j = 0; j < n_slp*n_jmu; j++){
#ifndef TIME_CURR
      average_vector_current_and_sum(1, j_mu[j], jlow_mu[j%n_jmu]);
      int status = write_vector_current_record(outfile[j], 0, 1, mass[j%n_jmu], j_mu[j]);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, outfilenames[j]);
      } 
#else
      wtime -= dclock();
      write_total_vector_current_record_t(0, jlow_mu[j%n_jmu], outfile, j);
      wtime += dclock();
#endif
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
#ifdef PT2PT /* DEBUG */
    ex=jrand%nx;ey=(jrand/nx)%ny;ez=(jrand/(nx*ny))%nz;et=jrand/(nx*ny*nz);
    thin_parity = (ex+ey+ez+et)%2==0?EVEN:ODD;
    for(int color=0;color<3;color++){
      clear_v_field(gr);
      if(this_node==node_number(ex,ey,ez,et)) gr[node_index(ex,ey,ez,et)].c[color]=(complex){1,0};
#else
    //if(this_node==0) ran+=gr0[0].c[0].real*gr0[0].c[1].real+gr0[0].c[0].imag*gr0[0].c[1].imag;
    //for(int color=0;color<3;color++) node0_printf("gr_%d= (%e,%e)\n",color,gr0[0].c[color].real,gr0[0].c[color].imag);
    /* Iterate over displacements within a d^4 cube. Use even displacements only */
    for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)if(parity==EVENANDODD || (ex+ey+ez+et)%2==parity%2)for(int color=0;color<isColorThin*3+(int) pow(0,isColorThin);color++){
      thin_parity = (ex+ey+et+ez)%2==0?EVEN:ODD;

      // Can't do this now that we are doing deflation.
      // We would need to rephase the eigenvectors
      //	    r_offset[0] = ex; r_offset[1] = ey; r_offset[2] = ez; r_offset[3] = et;

      /* Apply source thinning */
      copy_v_field(gr, gr0);
      thin_source_color(gr, d, ex, ey, ez, et, color, isColorThin);
#endif
      /* Project out the low mode part, based on the given eigenvectors */
      project_out(gr, eigVec, Nvecs, thin_parity);
      
#if 1 /* DEBUG */
      /* Check the norm of the reduced source */
      double_complex dd;
      dot_product(gr, gr, &dd, EVEN);
      node0_printf("Deflated source norm %g\n", dd.real);
#endif

      /* Start inversion */
      /* M_inv_gr = M^{-1} gr (same random source for each mass) */
      /* Assume: The residual of inversion of mass difference is the same. */
      shift=0;
      for(j = 0; j < n_jmu; j++){
	if(j<numMassPair){
	  shift++;      
	  node0_printf("Solving for %d %d %d %d mass %g-%g\n", ex, ey, ez, et, mass[j+shift-1],mass[j+shift]);
	}else
	  node0_printf("Solving for %d %d %d %d mass %g\n", ex, ey, ez, et, mass[j+shift]);
	
	clear_v_field(M_inv_gr);
	clear_v_field(M_inv_gr0);
	for(ns=0;ns<n_slp;ns++){
	  (qic+j+shift)->resid = precs_sloppy[j+shift+ns*n_masses];
	  node0_printf("Goal Sloppy Residual: %g\n",(qic+j+shift)->resid);
	  if(j<numMassPair){
	    mat_invert_uml_field_projected( gr, M_inv_gr0, qic + j + shift, thin_parity, mass[j+1], fn_multi[j+shift]);
	    mat_invert_uml_field_projected( M_inv_gr0, M_inv_gr, qic + j + shift, thin_parity, mass[j+shift-1], fn_multi[j+shift-1]);
	  }
	  else
	    mat_invert_uml_field_projected( gr, M_inv_gr, qic + j+shift, thin_parity, mass[j+shift], fn_multi[j+shift]);
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
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, gr_mu, M_inv_gr);
	    spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
	    /* J_mu = imag[gr.M_inv_gr] */
	    FORALLFIELDSITES(i){
	      complex cc = su3_dot( gr+i, gr_mu+i );
	      if(j<numMassPair) cc.imag*=4*(mass[j+shift]*mass[j+shift]-mass[j+shift-1]*mass[j+shift-1]);
	      j_mu[n_jmu*ns+j][NMU*i + mu] += cc.imag;
#ifdef TESTEO
	      if(cc.real!=0 || cc.imag!=0) printf("j_mu_%d %d %d %d %d %d %g %g\n",ns, lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu, cc.real, cc.imag);
#endif
#if defined(FINDCOV) || defined(FINDVAR) || defined(TIME_CURR)
	      j_mur[n_jmu*ns+j][NMU*i + mu] += cc.imag;
#endif
	    } /* i: FORALLFIELDSITES*/
	  
#if 0		/* DEBUG */
	    su3_vector *gr_mu_test = create_v_field();
	    su3_vector *grp = create_v_field();
	    su3_vector *grpp = create_v_field();
	    
	    for(n = 0; n < Nvecs; n++){
	      dslash_fn_field(eigVec[n], grp, ODD, fn_multi[j]);
	      spin_taste_op_fn(fn_multi[j+shift], spin_taste[mu], r_offset, grpp, grp);
	      spin_taste_op_fn(fn_multi[j+shift], spin_taste_index("pion05"), r_offset, grpp, grpp);
	      
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

	  } /* mu: NMU */
	  /****** ADD DEBUG BEGIN ************************/
#if 0
	  FORSOMEFIELDPARITY(i,p){
	    printf("j_mu %d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	    for(mu = 0; mu < NMU; mu++)
	      printf("%g ",j_mu[n_jmu*ns+j][NMU*i + mu]);
	    printf("\n");
	  }
#endif
	/****** ADD DEBUG END ************************/
	  if(j<numMassPair) copy_v_field(M_inv_gr,M_inv_gr0);
	} /* ns: n_slp */
      } /* j: n_jmu */
    } /* ex, ey, ez, et, color */

    /* Compute the variance of the total current density on each time slice */
#if defined(FINDVAR) || defined(FINDCOV)
    for(j=0;j<n_slp*n_jmu;j++){
      for(i=0;i<NMU*nt;i++) jt_mu[i]=0;
#ifdef FINDVAR
      for(mu=0;mu<NMU;mu++){
	FORALLFIELDSITES(i){
	  jt_mu[NMU*lattice[i].t+mu]+=jlow_mu[j%n_jmu][NMU*i + mu]+j_mur[j][NMU*i+mu];
	}
      }
      g_vecfloatsum(jt_mu,NMU*nt);
      for(i=0;i<NMU*nt;i++) jt2_mu[j][i]+=jt_mu[i]*jt_mu[i]/nwrite;
#endif
#ifdef FINDCOV
#if 0 // ver1
      for(mu=0;mu<NMU;mu++){
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
#if 0 // ver2
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
      clear_r_array_field(j_mur[j],NMU);
    } /* j:n_slp*n_jmu */
#endif
#ifdef  TIME_CURR
    wtime -= dclock();
    for(j=0;j<n_slp*n_jmu;j++)
      write_total_vector_current_record_t(jrand,j_mu[j], outfile, j);
    clear_r_array_field(j_mu[j],NMU);
    wtime += dclock();
#endif

#ifndef TIME_CURR
    /* Write record */
    if((jrand+1) % nwrite == 0){
      wtime -= dclock();
      for(j = 0; j < n_slp*n_jmu; j++){
	if(j%n_jmu==0) shift=0;
	if(j%n_jmu<numMassPair) shift++;
#ifdef PT2PT   /* DEBUG */
	Real msg=0;
	for(mu=0;mu<NMU;mu++){
	  if(this_node==node_number(ex,ey,ez,et)) msg=j_mu[j%n_jmu][NMU*node_index(ex,ey,ez,et)+mu];
	  if(node_number(ex,ey,ez,et)!=0){
	    if(this_node==node_number(ex,ey,ez,et))
	      send_field((char *)&msg,sizeof(Real),0);
	    if(this_node==0)
	      get_field((char *)&msg,sizeof(Real),node_number(ex,ey,ez,et));
	  }
	  g_sync();
	  node0_printf("J( %d %d %d %d )_%d= %g for mass %d\n",ex,ey,ez,et,mu,msg,j);fflush(stdout);
	}
#endif
#if defined(FINDVAR) || defined(FINDCOV)  /* DEBUG */
	int x,y,z,t;
	double msg[NMU];
	for(int k=0;k<nt*NMU;k++) jt_mu[k]=0;
	for(mu=0;mu<NMU;mu++){
	  FORALLFIELDSITES(i){
	    jt_mu[lattice[i].t*NMU+mu]+=nwrite*jlow_mu[j%n_jmu][NMU*i + mu]+j_mu[j][NMU*i+mu];
	  }
	}
	g_vecfloatsum(jt_mu,nt*NMU);
#ifdef FINDVAR
	for(i=0;i<NMU*nt;i++){
	  jt2_mu[j][i]-=jt_mu[i]*jt_mu[i]/nwrite/nwrite;
	  node0_printf("j_t(%d)_%d= %e %e for ns %d mass %d\n",i/NMU,i%NMU,jt_mu[i]/nwrite,sqrt(jt2_mu[j][i]/nwrite),j/n_jmu,j%n_jmu);
	}
#endif
	// print j_mu at all sites
	FORALLFIELDSITESLIN(x,y,z,t,nx,ny,nz,nt){
	  for(mu=0;mu<NMU;mu++) msg[mu]=0;
	  if(this_node==node_number(x,y,z,t))for(mu=0;mu<NMU;mu++) msg[mu]=j_mu[j][NMU*node_index(x,y,z,t)+mu];
	  g_vecfloatsum(msg,NMU);
	  for(mu=0;mu<NMU;mu++) 
	    node0_printf("j(%d,%d,%d,%d)_%d= %.16e for mass%d and ns=%d\n",x,y,z,t,mu,msg[mu]/nwrite,j%n_jmu,j/n_jmu);
	}
#ifdef FINDCOV //computing the covariance matrix
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
        int rep_mass = mass[j%n_jmu+shift];
        if(j%n_jmu<numMassPair){
          rep_mass =-mass[j%n_jmu]; // order switched to make the difference positive
          node0_printf("For rand %d and mass %g-%g\n", jrand, mass[j%n_jmu+shift-1],mass[j%n_jmu+shift]);fflush(stdout);
        }
        else
          node0_printf("For rand %d and mass %g\n", jrand,rep_mass);fflush(stdout);
	average_vector_current_and_sum(nwrite, j_mu[j], jlow_mu[j%n_jmu]);
	int status = write_vector_current_record(outfile[j], jrand, nwrite, rep_mass, j_mu[j]);
	clear_r_array_field(j_mu[j], NMU);
	if(status != QIO_SUCCESS)
	  node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      } /* j: n_slp*n_jmu */
      wtime += dclock();
    } /* if write */
#endif /* END: ifndef TIME_CURR*/
  } /* jrand */
  //node0_printf("grxgr= %e\n",ran/nrand);

  /* Cleanup */
  for(j = 0; j < n_slp*n_jmu; j++){
#ifdef TIME_CURR
    if(this_node==0) fclose(outfile[j]);
    free(outfile);outfile=NULL;
#else
    close_vector_current_file(outfile[j]);
#endif    
    destroy_r_array_field(j_mu[j], NMU);
    if(j/n_masses==0) destroy_r_array_field(jlow_mu[j], NMU);
#ifdef FINDVAR
    destroy_r_array_field(j_mur[j],NMU);
    free(jt2_mu[j]);
  }
   jt2_mu[0]=NULL;
#else
  }
#endif

  destroy_v_field(M_inv_gr0); M_inv_gr0 = NULL;
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_multi_eig_eo: iterated single-mass inverter */

#endif //END: if-else BLOCKCG clause; still within EIGMODE

/************************************************************************
 * Entry point for multiple masses with iterated single-mass inverter.
 * Designed for use with deflation or eigcg
need to be tested
 ************************************************************************/
void 
f_meas_current_pt2pt( int n_masses, quark_invert_control *qic, 
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
f_meas_current_diff( int n_masses, int const nrand, int nwrite, int thinning,
                     quark_invert_control *qic_precise,
                     quark_invert_control *qic_sloppy,
                     su3_vector **eigVec, double *eigVal, int Nvecs, ks_param *ksp,
                     fermion_links_t *fl, Real *precs_sloppy, int n_slp,
                     char filenames[][MAXFILENAME],
                     int parity, int isCorr, int numMassPair, int isColorThin){
  
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
f_meas_current( int n_masses, int nrand, int nwrite, int thinning,
		quark_invert_control *qic,
		su3_vector **eigVec, double *eigVal, int Nvecs,
		ks_param *ksp, fermion_links_t *fl,
		Real *precs_sloppy, int n_slp,
		char filenames[][MAXFILENAME],
		int parity, int numMassPair, int isColorThin){
  
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
    for(j = 0; j < NMU*n_masses; j++) j_mu[j]=0;
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
    g_vecfloatsum(j_mu,NMU*n_masses);

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

