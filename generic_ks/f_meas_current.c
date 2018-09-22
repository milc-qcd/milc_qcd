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

#define MASS_UDLSC /* Calculate densities with ud and ls mass differences */
#define OPT_UDLSC  /* If defined, approximate the u-d difference using a sequential propagator based on ml instead of mu or md */
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

/*Write current record for the accumulated average over random sources */
static void
write_tslice_values(char *tag, int jr, Real mass1, Real mass2, Real *j_mu ){
  double *jtmu = (double *)malloc(sizeof(double)*4*param.nt);

  /* Loop over the set of random sources in this block */
  if(jr < 0){
    if(mass2 == 0){
      node0_printf("JTMU%s LOW MODE MASS %g\n", tag, mass1);
    } else {
      node0_printf("JTMU%s LOW MODE MASS %g minus MASS %g\n", tag, 
		   mass1, mass2);
    }
  } else {
    if(mass2 == 0){
      node0_printf("JTMU%s RANDOM_SOURCE %d MASS %g\n", tag, jr, mass1);
    } else {
      node0_printf("JTMU%s RANDOM_SOURCE %d MASS %g minus MASS %g\n", tag, 
		   jr, mass1, mass2);
    }
  }

  for(int tmu = 0; tmu < 4*param.nt; tmu++){
    jtmu[tmu] = 0.;
  }
  int i;
  FORALLFIELDSITES(i){
    for(int mu = 0; mu < 4; mu++){
      jtmu[4*lattice[i].t + mu] += j_mu[4*i + mu];
    }
  }
  for(int t = 0; t < param.nt; t++){
    if(jr < 0){
      node0_printf("JTMU%s %d ", tag,  t);
    } else {
      node0_printf("JTMU%s %d %d ", tag, jr, t);
    }
    for(int mu = 0; mu < 4; mu++){
      g_doublesum(&jtmu[4*t+mu]);
      node0_printf("%g ", jtmu[4*t+mu]);
    }
    node0_printf("\n");
  }

  free(jtmu);
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
    collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0);
    //    node0_printf("ODD sources\n");
    collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0);
  } /* jr */

  destroy_v_field(gr0);
}

/************************************************************************/
/* Calculate current densities for a given quark mass and source
   parity using a list of thinned stochastic estimators */
/* nsrc sources in gr.  Results in j_mu */

static void
block_current_stochastic( int nr, Real *j_mu_mass[], Real mass, int nsrc, int sign, 
			  int parity, quark_invert_control *qic, 
			  imp_ferm_links_t *fn_mass, su3_vector *gr[]){
  
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
      int evol = nsrc/nr;  /* Number of thinned source sites */
      int ir = is/evol;    /* Random source block index */
      int i;
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_mu+i );
	j_mu_mass[ir][NMU*i + mu] += -sign*cc.imag;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu_mass[ir][NMU*i + mu]);
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
block_current_stochastic_deltam( Real *j_mu01[], Real mass0, Real mass1,
				 Real dm2_014,
				 imp_ferm_links_t *fn01, int nsrc,
				 int sign, int parity, int nr,
				 quark_invert_control *qic, su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic_deltam";

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

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
      int evol = nsrc/nr;  /* Number of thinned source sites */
      int ir = is/evol;    /* Random source block index */
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_mu+i );
	j_mu01[ir][NMU*i + mu] += -sign*cc.imag*dm2_014;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu01[ir][NMU*i + mu]);
#endif
      }
    } /* is, mu */
  destroy_v_field(gr_mu);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M_inv_gr[is]);

} /* block_current_stochastic_deltam */
      
/************************************************************************/
/* Calculate three cobinations of current densities based on four sets of
   masses: u, d, l, s using a list of thinned stochastic
   estimators. Do this for a given parity.
   In particular, calculate J_u - J_d and  J_l - J_s
   nsrc sources in gr.  Results in j_mu[0], j_mu[2] 
*/

static void
block_current_stochastic_delta_udls( Real **j_mu[], Real masses[],
				     imp_ferm_links_t *fn[], int nsrc,
				     int sign, int parity, int nr,
				     quark_invert_control qic[],
				     su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic_delta_udls";

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_l = masses[2];
  Real m_s = masses[3];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_ls = j_mu[2];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_ls = fn[2];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_ud = &qic[0];
  quark_invert_control *qic_ls = &qic[2];

  su3_vector *Mls_inv_gr[nsrc];
  su3_vector *Mud_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    Mls_inv_gr[is] = create_v_field();
    Mud_inv_gr[is] = create_v_field();
  }

  int otherparity;
  if(parity == EVEN)
    otherparity=ODD;
  else
    otherparity=EVEN;

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *Ml_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    Ml_inv_gr[is] = create_v_field();
  }

  qic->parity = parity;
  /* Ml_inv_gr = 1/[D^2 + 4*m_l^2] gr */
  ks_congrad_block_field(nsrc, gr, Ml_inv_gr, qic, m_l, fn_ls);
  /* Mud_inv_gr = 1/[D^2 + 4*m_l^2] Ml_inv_gr */
  /* NOTE: we are approximating ud here */
  ks_congrad_block_field(nsrc, Ml_inv_gr, Mud_inv_gr, qic, m_l, fn_ud);
  /* Mls_inv_gr = 1/[D^2 + 4*m_s^2] Ml_inv_gr */
  ks_congrad_block_field(nsrc, Ml_inv_gr, Mls_inv_gr, qic, m_s, fn_ls);

  for(int is = 0; is < nsrc; is++){
      destroy_v_field(Ml_inv_gr[is]);
  }
  
  /* Might be better to use a block dslash here?? */
  for(int is = 0; is < nsrc; is++){
    dslash_fn_field( Mls_inv_gr[is], Mls_inv_gr[is], otherparity, fn_ls);
    dslash_fn_field( Mud_inv_gr[is], Mud_inv_gr[is], otherparity, fn_ud);
  }
  
  /* For each source, apply current in various directions at the sink */
  su3_vector *gr_ls_mu = create_v_field();
  su3_vector *gr_ud_mu = create_v_field();
  for(int is = 0; is < nsrc; is++)
    for(int mu = 0; mu < NMU; mu++){
      
      /* Apply the appropriate spin_taste operator for
	 a nearly conserved current.  */
      spin_taste_op_fn(fn_ls, spin_taste[mu], r_offset, gr_ls_mu, Mls_inv_gr[is]);
      spin_taste_op_fn(fn_ls, spin_taste_index("pion05"), r_offset, gr_ls_mu, gr_ls_mu);
      spin_taste_op_fn(fn_ls, spin_taste[mu], r_offset, gr_ud_mu, Mud_inv_gr[is]);
      spin_taste_op_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_ud_mu, gr_ud_mu);
      
      /* J_mu = -imag[gr * Mls_inv_gr] */
      /* Multiply result by the mass difference squared times four */
      int i;
      int evol = nsrc/nr;  /* Number of thinned source sites */
      int ir = is/evol;    /* Random source block index */
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_ls_mu+i );
	j_mu_ls[ir][NMU*i + mu] += -sign*cc.imag*dm2_ls4;
	cc = su3_dot( gr[is]+i, gr_ud_mu+i );
	j_mu_ud[ir][NMU*i + mu] += -sign*cc.imag*dm2_ud4;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu_ls[ir][NMU*i + mu], j_mu_ud[ir][NMU*i + mu]);
#endif
      }
    } /* is, mu */
  destroy_v_field(gr_ud_mu);
  destroy_v_field(gr_ls_mu);

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(Mud_inv_gr[is]);
    destroy_v_field(Mls_inv_gr[is]);
  }

} /* block_current_stochastic_delta_udls */
      
/************************************************************************/

/* Calculate the current density difference between fine and sloppy
   solves separately for all masses */

static void
block_current_diff(int n_masses, Real **j_mu[], Real masses[],
		   imp_ferm_links_t *fn_mass[], 
		   quark_invert_control *qic_precise,
		   quark_invert_control *qic_sloppy,
		   int nsrc, int nr, su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "block_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){
    
    /* First, the sloppy high-mode solution */
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", masses[j]);
    
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, -1, EVEN,
			      qic_sloppy + j, fn_mass[j], gr_even);
    //      node0_printf("j_mu[%d][0] = %g\n",j,j_mu[j][0]);
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, -1, ODD,
			      qic_sloppy + j, fn_mass[j], gr_odd);
    //      node0_printf("j_mu[%d][4*node_index(1,0,0,0)] = %g\n",j,j_mu[j][4*node_index(1,0,0,0)]);
    
    /* Next, continue to a "precise" solution from the same sources */
    node0_printf("Solving precisely for all EVEN displacements for mass %g\n", masses[j]);
    
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, +1, EVEN,
			      qic_precise + j, fn_mass[j], gr_even);
    node0_printf("Solving precisely for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, +1, ODD,
			      qic_precise + j, fn_mass[j], gr_odd);
  } /* j */
  
} /* block_current_diff */

/************************************************************************/
/* Calculate the current density difference between fine and sloppy
   solves separately for three masses and take the difference: the
   mass0 result minus the mass1 result. 
   Results are in j_mu[0] and j_mu[2]
*/

static void
block_currents_diff_delta_lsc( int n_masses, Real **j_mu[], Real masses[],
			       imp_ferm_links_t *fn[], int nsrc, int nr,
			       quark_invert_control qic_precise[],
			       quark_invert_control qic_sloppy[],
			       su3_vector **gr_even, su3_vector **gr_odd){

  char myname[] = "block_currents_diff_delta_lsc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  Real m_l = masses[0];
  Real m_s = masses[1];
  Real m_c = masses[2];
  imp_ferm_links_t *fn_ls = fn[0];
  imp_ferm_links_t *fn_c = fn[2];
  Real **j_mu_ls = j_mu[0];
  Real **j_mu_c = j_mu[2];
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_sloppy_ls = &qic_sloppy[0];
  quark_invert_control *qic_sloppy_c = &qic_sloppy[2];
  quark_invert_control *qic_precise_ls = &qic_precise[0];
  quark_invert_control *qic_precise_c = &qic_precise[2];

  /* First, the sloppy high-mode solution */
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, -1, EVEN,
				   nr, qic_sloppy_ls, gr_even);
  if(n_masses == 3){
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, EVEN,
			      qic_sloppy_c, fn_c, gr_even);
  }

  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, -1, ODD,
				   nr, qic_sloppy_ls, gr_odd);
  if(n_masses == 3){
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, ODD,
			      qic_sloppy_c, fn_c, gr_odd);
  }
  
  /* Next, continue to a "precise" solution from the same sources */
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, +1, EVEN,
				   nr, qic_precise_ls, gr_even);
  if(n_masses == 3){
    node0_printf("Solving precisely for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_precise_c, fn_c, gr_even);
  }

  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, +1, ODD,
				   nr, qic_precise_ls, gr_odd);
  if(n_masses == 3){
    node0_printf("Solving precisely for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, ODD,
			      qic_precise_c, fn_c, gr_odd);
  }
  
} /* block_currents_diff_delta_lsc */

/************************************************************************/
/* Calculate the current density difference between fine and sloppy
   solves separately for five masses and take differences: the mass0
   result minus the mass1 result and the mass2 result minus the mass4
   result. 
   Results are in j_mu[0], j_mu[2], and j_mu[4]
*/

static void
block_currents_diff_delta_udlsc(int n_masses, Real **j_mu[], Real masses[],
				imp_ferm_links_t *fn[], int nsrc, int nr,
				quark_invert_control qic_precise[],
				quark_invert_control qic_sloppy[],
				su3_vector **gr_even, su3_vector **gr_odd){

  char myname[] = "block_current_diff_delta_udlsc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_l = masses[2];
  Real m_s = masses[3];
  Real m_c = masses[4];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_ls = j_mu[2];
  Real **j_mu_c  = j_mu[4];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_ls = fn[2];
  imp_ferm_links_t *fn_c  = fn[4];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_sloppy_ud = &qic_sloppy[0];
  quark_invert_control *qic_sloppy_ls = &qic_sloppy[2];
  quark_invert_control *qic_sloppy_c  = &qic_sloppy[4];
  quark_invert_control *qic_precise_ud = &qic_precise[0];
  quark_invert_control *qic_precise_ls = &qic_precise[2];
  quark_invert_control *qic_precise_c  = &qic_precise[4];

  /* The sloppy high-mode solution on even sites */
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, -1, EVEN, nr,
				       qic_sloppy, gr_even);
#else
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, -1, EVEN,
				   nr, qic_sloppy_ud, gr_even);
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, -1, EVEN,
				   nr, qic_sloppy_ls, gr_even);
#endif
  
  if(n_masses == 5){
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, EVEN,
			      qic_sloppy_c, fn_c, gr_even);
  }


  /* The sloppy high-mode solution on odd sites */

#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, -1, ODD, nr, 
				       qic_sloppy, gr_odd);
#else
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, -1, ODD,
				   nr, qic_sloppy_ud, gr_odd);
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, -1, ODD,
				   nr, qic_sloppy_ls, gr_odd);
#endif
  if(n_masses == 5){
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, ODD,
			      qic_sloppy_c, fn_c, gr_odd);
  }
  
  /* The "precise" solution from the same sources on even sites */

#ifdef OPT_UDLSC
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, EVEN, nr, 
				       qic_precise, gr_even);
#else
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, +1, EVEN,
				   nr, qic_precise_ud, gr_even);
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, +1, EVEN,
				   nr, qic_precise_ls, gr_even);
#endif
  if(n_masses == 5){
    node0_printf("Solving precisely for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_precise_c, fn_c, gr_even);
  }

  /* The "precise" solution from the same sources on odd sites */

#ifdef OPT_UDLSC
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, ODD, nr, 
				       qic_precise, gr_odd);
#else
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, +1, ODD,
				   nr, qic_precise_ud, gr_odd);
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc, +1, ODD,
				   nr, qic_precise_ls, gr_odd);
#endif
  if(n_masses == 5){
    node0_printf("Solving precisely for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, ODD,
			      qic_precise_c, fn_c, gr_odd);
  }
  
} /* block_currents_diff_delta_udlsc */

/************************************************************************/
/* Calculate the current density difference between fine and sloppy
   solves separately for three masses and take the difference: the
   mass0 result minus the mass1 result. */

static void
block_currents_diff_deltam(int n_masses, Real **j_mu[], Real masses[],
			   imp_ferm_links_t *fn_mass[], 
			   quark_invert_control *qic_precise,
			   quark_invert_control *qic_sloppy, int nsrc,
			   int nr, su3_vector **gr_even, su3_vector **gr_odd){

  char myname[] = "block_currents_diff_deltam";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  if(n_masses <= 3){
    /* Take the difference between the first two masses, but not the third */
    block_currents_diff_delta_lsc(n_masses, j_mu, masses, fn_mass, nsrc, nr,
				  qic_precise, qic_sloppy, gr_even, gr_odd);
  } else {
    /* Take the difference between the first two masses, and the next two, but not fifth */
    block_currents_diff_delta_udlsc(n_masses, j_mu, masses, fn_mass, nsrc, nr,
				    qic_precise, qic_sloppy, gr_even, gr_odd);
  }
  
} /* block_currents_diff_deltam */

/************************************************************************/
/* Entry point for multiple masses with deflation and iterated single-mass inverter.
   This variant does two solves from the same source -- sloppy and precise --
   and calculates the average of the difference between the resulting current
   densities.  Designed for use with deflation or eigcg. 
   Requires a set of accurate low-mode eigenpairs */

void 
f_meas_current_diff( int n_masses, int nrand, int thinning,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy,
		     ks_param *ksp, fermion_links_t *fl, 
		     char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  int nr = 2;  /* Number of random sources to include in multi-rhs */
  if(nrand < nr)
    nr = nrand;

  /* Create fields for high-mode current densities, one for each
     multi-rhs randm source and mass */
  Real **j_mu[n_masses];
  for(int j = 0; j < n_masses; j++){
    j_mu[j] = (Real **)malloc(n_masses*sizeof(Real *));
    for(int ir = 0; ir < nr; ir++)
      j_mu[j][ir] = create_r_array_field(NMU);
  }

  /* Load arrays with masses and the HISQ link structure for each */
  imp_ferm_links_t *fn_mass[n_masses];
  Real masses[n_masses];
  imp_ferm_links_t **fn = get_fm_links(fl);
  for(int j = 0; j < n_masses; j++){
    masses[j] = ksp[j].mass;
    fn_mass[j] = fn[ksp[j].naik_term_epsilon_index];
  }

  /* Block solver parameters */
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
  for(int jrand = 0; jrand < nrand; jrand += nr){

    /* Create sources in gr_even and gr_odd */
    collect_sources(gr_even, gr_odd, nr, d, evol);
    
#ifdef MASS_UDLSC
    if(n_masses <= 1){
      block_current_diff(n_masses, j_mu, masses, fn_mass, 
			 qic_precise, qic_sloppy, nsrc, nr, gr_even, gr_odd);
    } else {
      block_currents_diff_deltam(n_masses, j_mu, masses, fn_mass, 
				 qic_precise, qic_sloppy, nsrc, nr, gr_even, gr_odd);
    }
#else
    block_current_diff(n_masses, j_mu, masses, fn_mass, 
                       qic_precise, qic_sloppy, nsrc, nr, gr_even, gr_odd);
#endif
    
   
    for(int j = 0; j < n_masses; j++){
      for(int ir = 0; ir < nr; ir++){
#if 0      
	/* DEBUG */
	FORALLFIELDSITES(i){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g ",mu, j_mu[j][ir][NMU*i + mu]);
	  printf("\n");
	}
#endif
#ifdef MASS_UDLSC
	if((j == 0 && n_masses > 1) || (j == 2 && n_masses > 3))
	  {
	    write_tslice_values("DIFF", jrand+ir, masses[j], masses[j+1], j_mu[j][ir]);
	  }
	else if((j == 0 && n_masses == 1) || (j == 2 && n_masses <= 3) || j == 4 ){
	  write_tslice_values("DIFF", jrand+ir, masses[j], 0., j_mu[j][ir]);
	}
#else
	write_tslice_values("DIFF", jrand+ir, masses[j], 0., j_mu[j][ir]);
#endif
	clear_r_array_field(j_mu[j][ir], NMU);
      } /* ir */
    } /* j */

  } /* jrand */

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  free(gr_even); gr_even = NULL;
  free(gr_odd); gr_odd = NULL;

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
exact_current_delta_ls(Real *jlow_mu[], Real masses[],
		       imp_ferm_links_t *fn[]){

  /* Compute the difference in exact low-mode current densities for two masses */
  /* Takes densities for mass0 minus densities for mass1 */

  int r_offset[4] = {0, 0, 0, 0};
  int *spin_taste = get_spin_taste();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  int i;
  Real m_l = masses[0];
  Real m_s = masses[1];
  Real *jlow_mu_ud = jlow_mu[0];
  imp_ferm_links_t *fn_ud = fn[0];
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);

  for(int n = 0; n < Nvecs; n++){
    dslash_fn_field(eigVec[n], gr0, ODD, fn_ud);
    for(int mu = 0; mu < NMU; mu++){
      
      spin_taste_op_fn(fn_ud, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	Real dl = eigVal[n]+4.0*m_l*m_l;
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_ud[NMU*i + mu] += -z.imag*dm2_ls4/(dl*ds);
      } /* i */
      
      spin_taste_op_fn(fn_ud, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FORODDFIELDSITES(i){
	complex z;
	z = su3_dot( gr0 + i, gr_mu + i);
	Real dl = eigVal[n]+4.0*m_l*m_l;
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_ud[NMU*i + mu] += z.imag*dm2_ls4/(dl*ds);
      } /* i */
    } /* mu */
  } /* n */

  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr0); gr0 = NULL;
}

/*********************************************************************/
/* Calculate exact low-mode current density contributions for flavor
   combinations u-d, s-u, and c where the masses are presented as u,
   d, l, s, c 
   Results are as follows:
   jlow_mu[0] = u - d
   jlow_mu[2] = s - u
*/

static void
exact_current_delta_udus(Real *jlow_mu[], Real masses[], imp_ferm_links_t *fn[]){

  int r_offset[4] = {0, 0, 0, 0};
  int *spin_taste = get_spin_taste();
  su3_vector *gr0 = create_v_field();
  su3_vector *gr_mu = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_l = masses[2];
  Real m_s = masses[3];
  Real m_c = masses[4];
  Real *jlow_mu_ud = jlow_mu[0];
  Real *jlow_mu_us = jlow_mu[2];
  imp_ferm_links_t *fn_udls = fn[0];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_us4 = 4*(m_s*m_s - m_u*m_u);
  int i;

  for(int n = 0; n < Nvecs; n++){
    dslash_fn_field(eigVec[n], gr0, ODD, fn_udls);
    for(int mu = 0; mu < NMU; mu++){
      
      spin_taste_op_fn(fn_udls, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_fn(fn_udls, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	Real du = eigVal[n]+4.0*m_u*m_u;
	Real dd = eigVal[n]+4.0*m_d*m_d;
	jlow_mu_ud[NMU*i + mu] += -z.imag*dm2_ud4/(du*dd);
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_us[NMU*i + mu] += -z.imag*dm2_us4/(du*ds);
      } /* i */
      
      spin_taste_op_fn(fn_udls, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_fn(fn_udls, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FORODDFIELDSITES(i){
	complex z;
	z = su3_dot( gr0 + i, gr_mu + i);
	Real du = eigVal[n]+4.0*m_u*m_u;
	Real dd = eigVal[n]+4.0*m_d*m_d;
	jlow_mu_ud[NMU*i + mu] += z.imag*dm2_ud4/(du*dd);
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_us[NMU*i + mu] += z.imag*dm2_us4/(du*ds);
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

  /* Compute exact low-mode current density for all masses 
     No differences taken. 
     Results in jlow_mu
  */

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
exact_currents_deltam(int n_masses, Real *jlow_mu[], Real masses[],
		      imp_ferm_links_t *fn_mass[]){

  /* Compute the difference in exact low-mode current densities 
     If two or three masses, take the difference for 
     mass0 and mass1 and if three, also compute for mass2
     Results in jlow_mu[0], jlow_mu[2]
     If four or five masses, take the difference for
     mass0 and mass1 and the difference for mass0 and mass3
     and if five, also compute for mass 4. 
     Results in jlow_mu[0], jlow_mu[2], and jlow_mu[4]
*/

  char myname[] = "exact_currents_deltam";
  
  double dtime = -dclock();

  switch(n_masses){
  case(1):
    exact_current(jlow_mu[0], masses[0], fn_mass[0]);
    break;
  case(2):
    exact_current_delta_ls(jlow_mu, masses,fn_mass);
    break;
  case(3):
    exact_current_delta_ls(jlow_mu, masses,fn_mass);
    exact_current(jlow_mu[2], masses[2], fn_mass[2]);
    break;
  case(4):
    exact_current_delta_udus(jlow_mu, masses, fn_mass);
    break;
  case(5):
    exact_current_delta_udus(jlow_mu, masses, fn_mass);
    exact_current(jlow_mu[4], masses[4], fn_mass[4]);
    break;
  default:
    node0_printf("%s: wrong number of masses %d\n", myname, n_masses);
    terminate(1);
  }
    
#if 0
  for(int mu = 0; mu < NMU; mu++){
    node0_printf("Exact low modes For mass %g minus mass %g\n",
		 masses[0], masses[1]);
    FORALLFIELDSITES(i){
      node0_printf("j_mu_low  %d %d %d %d %d %g\n",
		   lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		   mu, jlow_mu[0][NMU*i+mu]);
    }
    node0_printf("Exact low modes For mass %g minus mass %g\n", masses[0],
		 masses[2]);
    FORALLFIELDSITES(i){
      node0_printf("j_mu_low  %d %d %d %d %d %g\n",
		   lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		   mu, jlow_mu[2][NMU*i+mu]);
    }
  }
#endif

  dtime += dclock();
#ifdef CGTIME
  node0_printf("Time for exact low modes %g sec\n", dtime);
#endif
}

/*********************************************************************/
/* This version of the block current routine does all masses without
   any differences betweem results of different masses.
   Results in j_mu
 */

static void 
block_current( int n_masses, Real **j_mu[], Real masses[],
	       imp_ferm_links_t *fn_mass[], quark_invert_control *qic,
	       int nsrc, int nr, su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "block_current";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){
    node0_printf("Solving for all EVEN displacements for mass %g\n", masses[j]);
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, +1, EVEN,
			      qic + j, fn_mass[j], gr_even);
    node0_printf("Solving for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic( nr, j_mu[j], masses[j], nsrc, +1, ODD,
			      qic + j, fn_mass[j], gr_odd);
  } /* j */
}

/*********************************************************************/
/* This version takes three masses.
   It takes the differences of current densities for mass0 and mass1 
   and also computes mass2
   Results in j_mu[0] and j_mu[2]
*/
static void 
block_currents_delta_lsc( int n_masses, Real **j_mu[], Real masses[], 
			  imp_ferm_links_t *fn[], int nsrc, int nr,
			  quark_invert_control qic[],
			  su3_vector **gr_even, su3_vector **gr_odd){
  char myname[] = "block_currents_delta_lsc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);
  
  Real m_l = masses[0];
  Real m_s = masses[1];
  Real m_c = masses[2];
  imp_ferm_links_t *fn_ls = fn[0];
  imp_ferm_links_t *fn_c = fn[2];
  Real **j_mu_ls = j_mu[0];
  Real **j_mu_c = j_mu[2];
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_ls = &qic[0];
  quark_invert_control *qic_c = &qic[2];
  
  /* Construct current density from the list of sources */
  node0_printf("Solving for all EVEN displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls,
				   nsrc, +1, EVEN, nr, qic_ls, gr_even);
  if(n_masses == 3){
    node0_printf("Solving for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_c, fn_c, gr_even);
  }

  node0_printf("Solving for all ODD displacements for mass diff %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls,
				   nsrc, +1, ODD, nr, qic_ls, gr_odd);
  if(n_masses == 3){
    node0_printf("Solving for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, ODD,
			      qic_c, fn_c, gr_odd);
  }
}

/*********************************************************************/
/* This version takes five masses
   Takes the differences of mass0 and mass1 and of mass2 and mass3
   and computes mass4
   Results in j_mu[0], j_mu[2], and j_mu[4]
*/

static void 
block_currents_delta_udlsc( int n_masses, Real **j_mu[], Real masses[], 
			    imp_ferm_links_t *fn[], int nsrc, int nr,
			    quark_invert_control qic[],
			    su3_vector **gr_even, su3_vector **gr_odd){

  char myname[] = "block_currents_delta_udlsc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_l = masses[2];
  Real m_s = masses[3];
  Real m_c = masses[4];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_ls = j_mu[2];
  Real **j_mu_c  = j_mu[4];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_ls = fn[2];
  imp_ferm_links_t *fn_c  = fn[4];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_ud = &qic[0];
  quark_invert_control *qic_ls = &qic[2];
  quark_invert_control *qic_c  = &qic[4];

  /* Construct current density from the list of sources */

#ifdef OPT_UDLSC
  node0_printf("Solving for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, EVEN, nr, 
				       qic, gr_even);
#else
  node0_printf("Solving for all EVEN displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, EVEN, nr, qic_ud, gr_even);
  node0_printf("Solving for all EVEN displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc,
				   +1, EVEN, nr, qic_ls, gr_even);
#endif
  if(n_masses == 5){
    node0_printf("Solving for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_c, fn_c, gr_even);
  }
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, ODD, nr, 
				       qic, gr_odd);
#else
  node0_printf("Solving for all ODD displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, ODD, nr, qic_ud, gr_odd);
  node0_printf("Solving for all ODD displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc,
				   +1, ODD, nr, qic_ls, gr_odd);
#endif
  if(n_masses == 5){
    node0_printf("Solving for all ODD displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, ODD,
			      qic_c, fn_c, gr_odd);
  }
}


/*********************************************************************/
/* This version takes three or five masses

   With three, it takes the differences of current densities for
   mass0 and mass1 and computes mass2
   Results in j_mu[0] and j_mu[2]

   With five, takes the differences of current densities for mass0 and
   mass1 and of mass2 and mass3 and computes mass4
   Results in j_mu[0], j_mu[2], and j_mu[4]
   
*/
static void 
block_currents_deltam( int n_masses, Real **j_mu[], Real masses[], 
		       imp_ferm_links_t *fn[], quark_invert_control qic[],
		       int nsrc, int nr, su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "block_currents_deltam";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  if(n_masses <= 3){
    block_currents_delta_lsc(n_masses, j_mu, masses, fn, nsrc, nr,
			     qic, gr_even, gr_odd);
  } else {
    block_currents_delta_udlsc(n_masses, j_mu, masses, fn, nsrc, nr,
			       qic, gr_even, gr_odd);
  }      
}


/*********************************************************************/
/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */

void 
f_meas_current( int n_masses, int nrand, int thinning,
		quark_invert_control *qic, ks_param *ksp,
		fermion_links_t *fl, char filenames[][MAXFILENAME]){

  char myname[] = "f_meas_current";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  int i;
  int Nvecs = param.eigen_param.Nvecs;
  int nr = 2;  /* Number of random sources to block */
  if(nrand < nr)
    nr = nrand;

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

  /* Create fields for low-mode current densities, one for each mass */

  Real *jlow_mu[n_masses];
  for(int j = 0; j < n_masses; j++){
    if(Nvecs > 0)
      jlow_mu[j] = create_r_array_field(NMU);
    else
      jlow_mu[j] = NULL;
  }

  /* Create fields for high-mode current densities, one for each
     multi-rhs randm source and mass */
  Real **j_mu[n_masses];
  for(int j = 0; j < n_masses; j++){
    j_mu[j] = (Real **)malloc(n_masses*sizeof(Real *));
    for(int ir = 0; ir < nr; ir++)
      j_mu[j][ir] = create_r_array_field(NMU);
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
#ifdef MASS_UDLSC
  if(Nvecs > 0){
    if(n_masses <= 1){
      exact_currents(n_masses, jlow_mu, masses, fn_mass);
    } else {
      exact_currents_deltam(n_masses, jlow_mu, masses, fn_mass);
    }
  }
#else
  if(Nvecs > 0){
    exact_currents(n_masses, jlow_mu, masses, fn_mass);
  }
#endif

  /* Print the exact low mode contribution to the current density */
  if(Nvecs > 0)
    for(int j = 0; j < n_masses; j++){
#ifdef MASS_UDLSC
      if((j == 0 && n_masses > 1) || (j == 2 && n_masses > 3))
	{
	  write_tslice_values("LOW", -1, masses[j], masses[j+1], jlow_mu[j]);
	}
      else if((j == 0 && n_masses == 1) || (j == 2 && n_masses <= 3) || j == 4 ){
	write_tslice_values("LOW", -1, masses[j], 0., jlow_mu[j]);
      }
#else
      write_tslice_values("LOW", -1, masses[j], 0., jlow_mu[j]);
#endif
      clear_r_array_field(jlow_mu[j], NMU);
    } /* j */


  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  su3_vector **gr_even = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **gr_odd = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    gr_even[is] = create_v_field();
    gr_odd[is] = create_v_field();
  }

  /* Construct high-mode current density stochastically.  Collect
     random sources for block solves */

  for(int jrand = 0; jrand < nrand; jrand += nr){

    /* Create sources in gr_even and gr_odd */
    collect_sources(gr_even, gr_odd, nr, d, evol);
    
#ifdef MASS_UDLSC
    if(n_masses <= 1){
      block_current( n_masses, j_mu, masses, fn_mass, qic, nsrc, nr, gr_even, gr_odd);
    } else {
      block_currents_deltam( n_masses, j_mu, masses, fn_mass, qic,
			     nsrc, nr, gr_even, gr_odd);
    }
#else
    block_current( n_masses, j_mu, masses, fn_mass, qic, nsrc, nr, gr_even, gr_odd);
#endif
      
    for(int j = 0; j < n_masses; j++){
      for(int ir = 0; ir < nr; ir++){
#if 0
	/* DEBUG */
	node0_printf("For rand %d and mass %g\n", jrand+ir, masses[j]);
	FORSOMEFIELDPARITY(i, parity){
	  printf("%d %d %d %d ",lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t);
	  for(mu = 0; mu < NMU; mu++)
	    printf("%d %g %g ",mu,j_mu[j][ir][NMU*i + mu],jlow_mu[j][NMU*i + mu]);
	  printf("\n");
	}
#endif
	
#ifdef MASS_UDLSC
	if((j == 0 && n_masses > 1) || (j == 2 && n_masses > 3))
	  {
	    write_tslice_values("HI", jrand+ir, masses[j], masses[j+1], j_mu[j][ir]);
	  }
	else if((j == 0 && n_masses == 1) || (j == 2 && n_masses <= 3) || j == 4 ){
	  write_tslice_values("HI", jrand+ir, masses[j], 0., j_mu[j][ir]);
	}
#else
	write_tslice_values("HI", jrand+ir, masses[j], 0., j_mu[j][ir]);
#endif
	clear_r_array_field(j_mu[j][ir], NMU);
      } /* ir */
    } /* j */
  } /* jrand */
  
  for(int j = 0; j < n_masses; j++){
    for(int ir = 0; ir < nr; ir++)
      destroy_r_array_field(j_mu[j][ir], NMU);
    free(j_mu[j]);
  }
  
  if(Nvecs > 0)
    for(int j = 0; j < n_masses; j++)
      destroy_r_array_field(jlow_mu[j], NMU);

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(gr_even[is]);
    destroy_v_field(gr_odd[is]);
  }
  
  free(gr_even); gr_even = NULL;
  free(gr_odd); gr_odd = NULL;

  fflush(stdout);
  
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
f_meas_current_diff( int n_masses, int nrand, int thinning,
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

  //  node0_printf("Entered %s\n", myname); fflush(stdout);
 
  /* Create fields for current densities, one for each mass */
  for(j = 0; j < n_masses; j++)
    j_mu[j] = create_r_array_field(NMU);

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
      write_tslice_values("DIFF", jrand, ksp[j].mass, 0., j_mu[j]);
      clear_r_array_field(j_mu[j], NMU);
    } /* j */

  } /* jrand */
  
  for(j = 0; j < n_masses; j++){
    destroy_r_array_field(j_mu[j], NMU);
  }
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;

} /* f_meas_current_diff -- non-BLOCKCG version */

/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */

void 
f_meas_current( int n_masses, int nrand, int thinning,
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
  
  for(j = 0; j < n_masses; j++){
    write_tslice_values("LOW", 0, mass[j], 0., jlow_mu[j]);
    clear_r_array_field(jlow_mu[j], NMU);
  }

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

  for(j = 0; j < n_masses; j++){
    write_tslice_values("", jrand, mass[j], 0., j_mu[j]);
    clear_r_array_field(j_mu[j], NMU);
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
      clear_r_array_field(j_mu[j], NMU);
      if(status != QIO_SUCCESS){
	node0_printf("%s: Failed to write record to %s\n", myname, filenames[j]);
      }
    } /* j */
  } /* jrand */

  for(j = 0; j < n_masses; j++){
    destroy_r_array_field(j_mu[j], NMU);
    if(Nvecs > 0)
      destroy_r_array_field(jlow_mu[j], NMU);
  }
  
  destroy_v_field(M_inv_gr); M_inv_gr = NULL;
  destroy_v_field(gr_mu); gr_mu = NULL;
  destroy_v_field(gr); gr = NULL;
  destroy_v_field(gr0); gr0 = NULL;
} /* f_meas_current -- non-BLOCKCG version */

#endif  /* End of traditional version */

