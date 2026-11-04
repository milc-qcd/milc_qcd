/**************** f_meas_current.c ***************************************/
/* MIMD version 7 */
/* CD 1/15 */
/* Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.
*/

/* Measure the fermionic observable:

    J_mu and J_sea
    
    Entry points

    f_meas_sea_current_diff
    f_meas_sea_current
    f_meas_valence_current_diff
    f_meas_valence_current

    Quantities measured depend on OPT_UDLSC, MASS_UDLSC  and the n_masses 
   
    With undefined OPT_UDLSC and undefined MASS_UDLSC, compute single results for all given masses

    With OPT_UDLSC and MASS_UDLSC defined, compute as follows...
    n_masses = 1 Result for only one mass
               2 masses l s  Takes the difference j_l - j_s
	       3 masses l s c Takes the difference j_l - j_s and computes j_c
	       4 masses u d s c Takes the difference j_u - j_d and j_u - j_s and computes j_c
	       5 masses u d l s c Does the same as 4 for the low modes (ignoring l), 
	                          but for the high modes, computes j_u - j_l, j_l - j_s, j_s, and j_c
*/

#define MASS_UDLSC /* Calculate densities with ud and ls mass differences */
#define OPT_UDLSC  /* If defined, approximate the u-d difference using a sequential propagator based on ml instead of mu or md */
#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/fn_links.h"
#include "../include/fermion_links.h"
#include "../include/io_scidac.h"
#include "../include/imp_ferm_links.h"
#include "../include/openmp_defs.h"
#define LOOPEND
#include "../include/loopend.h"
#include <qio.h>
#include <string.h>
#ifdef HAVE_QUDA
#include "../include/generic_quda.h"
#else
#include "../include/openmp_defs.h"
#endif
#if defined(U1_VALENCE) || defined(U1_SEA)
#include "../include/generic_u1.h"
#endif

#define NMU 4
#define NRECINFO 128

/* Debugging */
#ifdef FIXEDSOURCE
#define NRMAX 3
#else
#define NRMAX 2
#endif

/*****************************************************************************/
/* Create spin-taste indices for current */
static int *
get_spin_taste(void){
  
  /* Current spin-taste list */
  char *spin_taste_list[NMU] = {"rhoxsfape", "rhoysfape", "rhozsfape", "rhotsfape"};
  static int spin_taste[NMU];
  int mu;
  
  /* Decode spin-taste label */
  for(mu = 0; mu < NMU; mu++){
    char dummy[33];
    strncpy(dummy, spin_taste_list[mu], 32);
    spin_taste[mu] = spin_taste_index(dummy);
  }
  
  return spin_taste;
}

/*****************************************************************************/
/* Thin the random source */
static void
thin_source(su3_vector *src, int thinning, int ex, int ey, int ez, int et){
  site *s;
  int i;

  FORALLSITES_OMP(i,s,) {
    if(s->x % thinning != ex || s->y % thinning != ey ||
       s->z % thinning != ez || s->t % thinning != et){
      clearvec(src+i);
    }
  } END_LOOP_OMP;
}

/*****************************************************************************/
/*Write current record for the accumulated average over random sources */
static void
write_tslice_values_begin(const char *tag){
  node0_printf("BEGIN JTMU%s\n", tag);
}

/*****************************************************************************/
static void
write_tslice_values_end(const char *tag){
  node0_printf("END JTMU%s\n", tag);
}

/*****************************************************************************/
/*Write current record for the accumulated average over random sources */
static void
write_sea_values_begin(const char *tag){
  node0_printf("BEGIN JSEA%s\n", tag);
}

/*****************************************************************************/
/*Write current record for the accumulated average over random sources */
static void
write_sea_values_group(const char *tag, Real mass, Real charge){
  node0_printf("JSEA%s MASS %g CHARGE %g\n", tag, mass, charge);
}

/*****************************************************************************/
static void
write_sea_values_end(const char *tag){
  node0_printf("END JSEA%s\n", tag);
}

/*****************************************************************************/
static void
write_tslice_values(const char *tag, int jr, Real mass1, Real charge1,
		    Real mass2, Real charge2, Real *j_mu ){
  double *jtmu = (double *)malloc(sizeof(double)*4*param.nt);

  /* Write header for this block of values on each time slice */
  /* Negative random index jr signals exact low mode */
  if(jr < 0){
    if(mass2 == 0){
      node0_printf("JTMU%s LOW MODE MASS %g CHARGE %g\n", tag, mass1, charge1);
    } else {
      node0_printf("JTMU%s LOW MODE MASS %g CHARGE %g minus MASS %g CHARGE %g\n",
		   tag, mass1, charge1,  mass2, charge2);
    }
  } else {
    if(mass2 == 0){
      node0_printf("JTMU%s RANDOM_SOURCE %d MASS %g CHARGE %g\n",
		   tag, jr, mass1, charge1);
    } else {
      node0_printf("JTMU%s RANDOM_SOURCE %d MASS %g CHARGE %g minus MASS %g CHARGE %g\n",
		   tag, jr, mass1, charge1, mass2, charge2);
    }
  }

  /* Sum the current density over the time slices */
  for(int tmu = 0; tmu < 4*param.nt; tmu++){
    jtmu[tmu] = 0.;
  }
  int i;
  FORALLFIELDSITES(i){
    for(int mu = 0; mu < 4; mu++){
      jtmu[4*lattice[i].t + mu] += j_mu[4*i + mu];
    }
  }

#if 0
  // DEBUG
  printf("JTMU0 %20.12f %20.12f %20.12f %20.12f\n", 
	 jtmu[0], jtmu[1], jtmu[2], jtmu[3]);
#endif

  /* Write the results for each time slice */
  for(int t = 0; t < param.nt; t++){
    if(jr < 0){
      node0_printf("JTMU%s %d ", tag,  t);
    } else {
      node0_printf("JTMU%s %d %d ", tag, jr, t);
    }
    for(int mu = 0; mu < 4; mu++){
      g_doublesum(&jtmu[4*t+mu]);
      node0_printf("%.10g ", jtmu[4*t+mu]);
    }
    node0_printf("\n");
  }

  free(jtmu);
}

/*****************************************************************************/
static void
write_sea_values(const char *tag, int jr, int jU1_A, Real j_sea ){

  /* Write the values for j_sea summed over the whole lattice */
  Real j_sea_out = j_sea;
  g_doublesum(&j_sea_out);
  /* Negative random index jr signals exact low mode */
  if(jr < 0){
    node0_printf("JSEA%s U1 %d %.10f\n", tag, jU1_A, j_sea_out);
  } else {
    node0_printf("JSEA%s RN %d U1 %d %.10f\n", tag, jr, jU1_A, j_sea_out);
  }
  
#if 0
  // DEBUG
  printf("JSEA0 %20.12f\n", jsea);
#endif
}

/*****************************************************************************/
/* Calculate and write the result for random source jr with the listed
   masses */
static void
write_jdotA_value(const char *tag, int jr, Real mass1, Real charge1,
		  Real mass2, Real charge2, Real *j_mu, Real *u1_A ){

  if(jr < 0){
    if(mass2 == 0){
      node0_printf("JdotA%s LOW MODE MASS %g CHARGE %g ", tag, mass1, charge1);
    } else {
      node0_printf("JdotA%s LOW MODE MASS %g CHARGE %g minus MASS %g CHARGE %g ",
		   tag, mass1, charge1, mass2, charge2);
    }
  } else {
    if(mass2 == 0){
      node0_printf("JdotA%s RANDOM_SOURCE %d MASS %g CHARGE %g ",
		   tag, jr, mass1, charge1);
    } else {
      node0_printf("JdotA%s RANDOM_SOURCE %d MASS %g CHARGE %g minus MASS %g CHARGE %g ",
		   tag, jr, mass1, charge1, mass2, charge2);
    }
  }

  double jdotA = 0.;
  int i;
  FORALLFIELDSITES_OMP(i, ){
    int mu;
    FORALLUPDIR(mu){
      jdotA += j_mu[4*i + mu] * u1_A[4*i + mu];
    }
  } END_LOOP_OMP;
  g_doublesum(&jdotA);
  node0_printf("%.10g\n", jdotA);
} 

#if 1

/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void
dot_product(su3_vector *vec1, su3_vector *vec2, 
	    double_complex *dot, int parity) {
  int i;
  complex cc ;
  double re = 0.0, im = 0.0;

  FORSOMEFIELDPARITY_OMP(i,parity, reduction(+:re,im)){
    complex temp1,temp2;
    CMULJ_(vec1[i].c[0],vec2[i].c[0],temp1) 
    CMULJ_(vec1[i].c[1],vec2[i].c[1],temp2)
    CSUM(temp1,temp2);
    CMULJ_(vec1[i].c[2],vec2[i].c[2],temp2)
    CSUM(temp1,temp2);
    re += temp1.real ;
    im += temp1.imag ;
  } END_LOOP_OMP;

  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

#else

/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void
static dot_product(su3_vector *vec1, su3_vector *vec2, 
		   double_complex *dot, int parity) {
  register double re,im ;
  register  int i;
  complex cc ;
  
  re=im=0.0;
  FORSOMEFIELDPARITY(i,parity){
    cc = su3_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  } END_LOOP;
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
} END_LOOP;

#endif

#if 1

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void
complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
		     su3_vector *vec2, int parity){

  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEFIELDPARITY_OMP(i,parity, ){

    double sr,si,br,bi,cr,ci;

    sr = cc->real; si = cc->imag;

    for(int j=0;j<3;j++){
	br=vec1[i].c[j].real; bi=vec1[i].c[j].imag;

	cr = sr*br - si*bi;
	ci = sr*bi + si*br;

	vec2[i].c[j].real -= cr;
	vec2[i].c[j].imag -= ci;
    }
  }
} END_LOOP_OMP;

#else

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void
complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
		     su3_vector *vec2, int parity){

  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEFIELDPARITY_OMP(i,parity, ){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  } END_LOOP_OMP;
}

#endif

/************************************************************************/
/*  Projects out the set of *vector from the  vec. Num is the number of vectors  *
 * and parity is the parity on which we work on.                           *
 * The vectors are assumed to be orthonormal.                              */
   
#if defined(HAVE_QUDA) && defined(USE_CURRENT_GPU)

/* This QUDA version uses the previously computed eigenvectors kept by QUDA */
/* "vector" is ignored */

static void
project_out(su3_vector *vec, su3_vector *vector[], int Num, int parity){
  register int i ;
  double_complex cc ;
  double ptime = -dclock();

  if(Num == 0)return;

  QudaParity qparity;

  switch(parity){
  case EVEN: qparity = QUDA_EVEN_PARITY; break;
  case ODD:  qparity = QUDA_ODD_PARITY; break;
  default:
    node0_printf("%s: ERROR. Bad parity value\n", __func__);
    terminate(1);
  }

  int nvec = 1;

  su3_vector *invecs[1], *outvecs[1];
  invecs[0] = vec;
  outvecs[0] = create_v_field();
  
  qudaProject(MILC_PRECISION, (void **)invecs, (void **)outvecs, nvec, Num, qparity);

  copy_v_field(vec, outvecs[0]);
  destroy_v_field(outvecs[0]);

  ptime += dclock();
#ifdef CGTIME
  if(parity == EVEN){
    node0_printf("Time to project out low modes from EVEN source %g sec\n", ptime);
  } else {
    node0_printf("Time to project out low modes from ODD source %g sec\n", ptime);
  }
#endif
}

#else

/* This is the original CPU version, which requires vector to be supplied  */

static void
project_out(su3_vector *vec, su3_vector *vector[], int Num, int parity){

  register int i ;
  double_complex cc ;
  double ptime = -dclock();

  if(Num == 0)return;

  int nzero = 0;
  for(i=Num-1;i>-1;i--){
    dot_product(vector[i], vec, &cc, parity) ;
    if(cc.real == 0. && cc.imag == 0.){
      if(nzero < 100){
	node0_printf("project_out got zero dot product\n");
      }
      nzero++;
    }
    complex_vec_mult_sub(&cc, vector[i], vec, parity);
  }

  ptime += dclock();
#ifdef CGTIME
  if(parity == EVEN){
    node0_printf("Time to project out low modes from EVEN source %g sec\n", ptime);
  } else {
    node0_printf("Time to project out low modes from ODD source %g sec\n", ptime);
  }
#endif
}

#endif

#if 0 // unused (debugging) function
/******************************************************************
*									*
*  void dumpvec( su3_vector *vec )					*
*  print out a 3 element complex vector					*
*/

static void
my_dumpvec( su3_vector *v ){
  int j;
  for(j=0;j<3;j++)printf("(%.12e,%.12e)\n",
			 v->c[j].real,v->c[j].imag);
}
#endif

/************************************************************************/
static void
collect_evenodd_sources(su3_vector *gr[], int ns, int parity, int thinning,
			su3_vector *gr0, int project){
  /* Fill even or odd thinned sources from a volume random source gr0 */
  /* Result in gr */

  /* Iterate over displacements within a d^4 cube for this parity. */
  int ex, ey, ez, et;
  int is = 0;
  int d = thinning;

  double dtime = -dclock();

  for(ex=0;ex<d;ex++)for(ey=0;ey<d;ey++)for(ez=0;ez<d;ez++)for(et=0;et<d;et++)
      if ( ((ex+ey+ez+et)%2==0 && parity == EVEN) ||
	   ((ex+ey+ez+et)%2==1 && parity == ODD) ){

#if 0
	node0_printf("Source %d is %d %d %d %d\n", is, ex, ey, ez, et);
#endif
	
	/* Apply source thinning */
	copy_v_field( gr[is], gr0 );
	thin_source( gr[is], d, ex, ey, ez, et );
	
#if 0
	// DEBUG
	int i;
	FORSOMEFIELDPARITY(i, parity){
	  if(lattice[i].t == 0){
	    printf("GR[%d] %d %d %d %d %30.20f %30.20f %30.20f %30.20f %30.20f %30.20f\n", is,
		   lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		   gr[is][i].c[0].real, gr[is][i].c[0].imag,
		   gr[is][i].c[1].real, gr[is][i].c[1].imag,
		   gr[is][i].c[2].real, gr[is][i].c[2].imag);
	  }
	} END_LOOP;
#endif
	
	is++;
	if(is > ns){
	  node0_printf("collect_evenodd_sources: Internal error: too many sources\n");
	  terminate(1);
	}
      } /* ex, ey, ez, et */

  /* Project out the low modes, if requested */

  int Nvecs = param.eigen_param.Nvecs;
  if(project && Nvecs > 0){
    for(int is = 0; is < ns; is++){
      project_out(gr[is], eigVec, Nvecs, parity);
#if 1
    /* DEBUG */
    /* Check the norm of the reduced source */
    double_complex dd;
    dot_product(gr[is], gr[is], &dd, parity);
    node0_printf("Deflated source norm gr[%d] %g\n", is, dd.real);
#endif
    }
  }
    
  dtime += dclock();
  node0_printf("Time to collect sources %g\n", dtime);
}

/************************************************************************/
/* Collect undiluted random sources for mrhs solves */

static su3_vector **
create_source_array(int nr){
  su3_vector **gr0 = (su3_vector **)malloc(nr*sizeof(su3_vector *));
  for(int jr = 0; jr < nr; jr++){
    gr0[jr] = create_v_field();
  }
  return gr0;
}

static void
fill_sources(su3_vector *gr0[], int nr){

  for(int jr = 0; jr < nr; jr++){
    
#if defined(Z2RSOURCE)
    z2rsource_plain_field( gr0[jr], EVENANDODD );
#elif defined(FIXEDSOURCE)
    /* For debugging. Should be run with nr = 3 */
    int i;
    int j = jr % 3;
    FORALLFIELDSITES(i){
      for(int k = 0; k < 3; k++){
	gr0[jr][i].c[k].real = 0.;
	gr0[jr][i].c[k].imag = 0.;
      }
      gr0[jr][i].c[j].real = 1.;
      gr0[jr][i].c[j].imag = 0.;
    }
#else
    grsource_plain_field( gr0[jr], EVENANDODD );
#endif

  } /* jr */
}
 
static void
destroy_source_array(su3_vector *gr0[], int nr){
  if(gr0 == NULL)return;
  for(int jr = 0; jr < nr; jr++)
    if(gr0[jr] != NULL)destroy_v_field(gr0[jr]);
  free(gr0);
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
      spin_taste_op_ape_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
      spin_taste_op_ape_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
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
      } END_LOOP;
    } /* is, mu */
  destroy_v_field(gr_mu);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M_inv_gr[is]);
  free(M_inv_gr);

} /* block_current_stochastic */

#ifdef U1_SEA

/************************************************************************/
static int
check_fn_match(imp_ferm_links_t *fn_U1, imp_ferm_links_t *fn_mass){

  /* The two sets of links have to differ by nothing but the U(1)
     phase, or their difference is not what we want.  scalar_mult_add_fn
     only warns about a phase mismatch, because an existing HISQ caller
     can reach it legitimately, so check it properly here. */

  char myname[] = "check_fn_match";

  if(fn_U1->eps_naik != fn_mass->eps_naik){
    node0_printf("%s: Naik epsilon mismatch: %g with U(1), %g without\n",
		 myname, fn_U1->eps_naik, fn_mass->eps_naik);
    return 1;
  }
  if(fn_U1->phase == NULL || fn_mass->phase == NULL){
    node0_printf("%s: missing link phase info\n", myname);
    return 1;
  }
  if(fn_U1->phase->twist_in != fn_mass->phase->twist_in){
    node0_printf("%s: boundary twist is %s with U(1) and %s without.  "
		 "The U(1) links are built fresh and untwisted, so twist "
		 "them to match before taking the difference.\n", myname,
		 fn_U1->phase->twist_in ? "on" : "off",
		 fn_mass->phase->twist_in ? "on" : "off");
    return 1;
  }
  for(int dir = XUP; dir <= TUP; dir++)
    if(fn_U1->phase->bdry_phase[dir] != fn_mass->phase->bdry_phase[dir] ||
       fn_U1->phase->r0[dir] != fn_mass->phase->r0[dir]){
      node0_printf("%s: link phase conventions differ in direction %d\n",
		   myname, dir);
      return 1;
    }

  return 0;
}

/************************************************************************/
static imp_ferm_links_t *
create_dslash_diff(Real *u1_A_tmp, Real charge, imp_ferm_links_t *fn_mass){

  /* Dummy values.  No support for nonzero Naik epsilon here. */
  int n_naiks = 1;
  double eps_naik[1] = { 0. };
      
  /* Fermion links for the thin links times this U(1) field.  The
     thin links in the site structure are not disturbed. */
  fermion_links_t *fl_U1 =
    create_fermion_links_from_site_u1(MILC_PRECISION, n_naiks, eps_naik, charge, u1_A_tmp);
  imp_ferm_links_t *fn_U1 = get_fm_links(fl_U1, 0);

  /* Apply boundary phase if requested */
  Real bdry_phase[4] = {0.,0.,0.,param.time_bc};
  set_boundary_twist_fn(fn_U1, bdry_phase, param.coord_origin);
  if(twist_status(fn_U1) == OFF)boundary_twist_fn(fn_U1, ON);
  
  if(check_fn_match(fn_U1, fn_mass))
    terminate(1);
  
  /* fn_diff = fn_U1 - fn_mass.  scalar_mult_add_fn checks that the
     two sets of links carry the same phase conventions and copies
     them to fn_diff. */
  
  /* Workspace for the difference links, reused for every U(1) field */
  imp_ferm_links_t *fn_diff = create_fn_links();
  
  scalar_mult_add_fn(fn_U1, fn_mass, -1., fn_diff);
  
  /* Rebuild the back links of the difference if this build uses them */
  /* Might not want back links because we have only a few Dslashes here */
  if(get_fatbacklinks(fn_mass) != NULL)
    load_fn_backlinks(fn_diff);

  destroy_fermion_links_from_site(fl_U1);
  return fn_diff;
}

/************************************************************************/
static void
block_current_stochastic_U1( int nr, Real j_sea[], Real charge, int nsrc,
			     int sign, int parity, 
			     quark_invert_control *qic,
			     su3_vector **M_inv_gr, su3_vector *dM_inv_gr,
			     imp_ferm_links_t *fn_mass, Real *u1_A_tmp, su3_vector *gr[])
{
  /* Dummy values.  No support for nonzero Naik epsilon here. */
  int n_naiks = 1;
  double eps_naik[1] = { 0. };

  imp_ferm_links_t *fn_diff = create_dslash_diff(u1_A_tmp, charge, fn_mass);
      
  /* dM_inv_gr = [ Dslash(fn_U1) - Dslash(fn_mass) ] M_inv_gr */
  /* Using the CPU version avoids wasteful host-device transfers */
  for(int is = 0; is < nsrc; is++){
    dslash_fn_field_cpu( M_inv_gr[is], dM_inv_gr, parity, fn_diff);
    
    /* For each source, close the trace */
    /* J = -imag[gr * M_inv_gr] */
    /* If sign = +1, add the result.  If sign = -1, subtract it */
    int evol = nsrc/nr;  /* Number of thinned source sites */
    int ir = is/evol;    /* Random source block index */
    int i;
    FORSOMEFIELDPARITY(i, parity){
      complex cc = su3_dot( gr[is]+i, dM_inv_gr+i );
      j_sea[ir] += -sign*cc.real;
#if 0
      printf("j_sea src[%d] %d %d %d %d %d %g\n", is,
	     lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
	     j_sea[ir]);
#endif
    } END_LOOP;
  } /* is */

  destroy_fn_links(fn_diff);

} /* block_current_stochastic_U1 */


/*********************************************************************/
static void
block_current_stochastic_U1_T2( int nr, Real j2_sea[], Real charge, int nsrc,
                                int sign, int parity, Real mass,
                                quark_invert_control *qic,
                                su3_vector **p, su3_vector **q, su3_vector **z,
                                imp_ferm_links_t *fn_mass, Real *u1_A, su3_vector *gr[])
{
  imp_ferm_links_t *fn_diff = create_dslash_diff(u1_A, charge, fn_mass);
  int i;

  for(int is = 0; is < nsrc; is++){
    /* q = delta-D p, both parities */
    dslash_fn_field_cpu( p[is], q[is], EVENANDODD, fn_diff );
  }

  /* The sequential solve: z = (M^dag M)^{-1} q, both parities */
  qic->parity = EVENANDODD;
  ks_congrad_block_field( nsrc, q, z, qic, mass, fn_mass );

  su3_vector *y = create_v_field();
  su3_vector *w = create_v_field();
  int evol = nsrc/nr;

  for(int is = 0; is < nsrc; is++){
    /* y = M^{-1} q = 2m z - D z */
    dslash_fn_field_cpu( z[is], y, EVENANDODD, fn_mass );
    FORALLFIELDSITES(i){
      scalar_mult_add_su3_vector( y+i, z[is]+i, -2.*mass, y+i );
      scalar_mult_su3_vector( y+i, -1., y+i );
    }

    /* w = delta-D y on the source parity, then close the trace */
    dslash_fn_field_cpu( y, w, parity, fn_diff );

    int ir = is/evol;
    FORSOMEFIELDPARITY(i, parity){
      complex cc = su3_dot( gr[is]+i, w+i );
      /* NOTE: no leading minus here, unlike T1.  There the stored
         M_inv_gr carries +D rather than the -D of M^dag, and the sign
         is put back at the accumulation.  Here p and y are the true
         M^{-1}, so the dot product is already the estimator. */
      j2_sea[ir] += sign*cc.real;
    } END_LOOP;
  }
  destroy_v_field(w);
  destroy_v_field(y);
  destroy_fn_links(fn_diff);
}

#endif /* U1_SEA */

/*********************************************************************/
/* Create fields for high-mode sea-quark delta D, one for each
   block of nr random sources and mass */
static void
clear_j_sea(int nr, Real j_sea[]){
  for(int ir = 0; ir < nr; ir++) j_sea[ir] = 0.;
}

/*********************************************************************/
/* Create fields for high-mode sea-quark delta D, one for each
   block of nr random sources and mass */
static Real *
create_j_sea(int nr){
  Real *j_sea = (Real *)malloc(nr*sizeof(Real));
  clear_j_sea(nr, j_sea);
  return j_sea;
}

/*********************************************************************/
/* Create fields for high-mode sea-quark delta D, one for each
   block of nr random sources and mass */
static void
combine_j_sea(int nr, Real c,
	      Real j_seaA[], Real j_seaB[], Real j_sea[]){
  for(int ir = 0; ir < nr; ir++)
    j_sea[ir] = j_seaA[ir] + c*j_seaB[ir];
}

/*********************************************************************/
static void
destroy_j_sea(Real j_sea[]){
 free(j_sea);
}

#ifdef U1_SEA

/************************************************************************/
/* Calculate the stochastic difference contribution Tr(delta Dslash
   M^{-1}) caused by introducing a U(1) field, for a given quark mass
   and source parity, using a list of thinned stochastic estimators.
   The calculation is repeated for each of the U(1) fields described
   by u1info.  nsrc sources in gr.  Results in j_mass[iU1][ir][i].

   For each U(1) field we build the fermion links fn_U1 from the thin
   links in the site structure multiplied by the U(1) phase
   exp(i*charge*A_mu(x)), form the difference fn_diff = fn_U1 - fn_mass,
   and apply dslash with the difference links.  Then complete the Tr.

   Dslash is linear in the links, so

       Dslash(fn_U1) - Dslash(fn_mass) = Dslash(fn_U1 - fn_mass) ,

   and differencing the links rather than the two dslash results avoids
   the cancellation that would otherwise dominate the rounding error.
*/

static void
block_current_stochastic_multiU1( int nr, int nU1_A_array, Real *jE_sea[], Real *jO_sea[],
				  Real charge, int nsrc,
				  int sign, int parity, 
				  quark_invert_control *qic, 
				  imp_ferm_links_t *fn_mass,
				  su3_vector *M_inv_gr[],
				  Real *u1_A_array[], su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic_multiU1";

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};
  su3_vector *dM_inv_gr = create_v_field();
  Real *j_sea = create_j_sea(nr);

  for(int iU1_A = 0; iU1_A < nU1_A_array; iU1_A++){

    /* Call with "charge" */
    clear_j_sea(nr, j_sea);
    block_current_stochastic_U1(nr, j_sea, charge, nsrc, sign, parity, qic,
				M_inv_gr, dM_inv_gr, fn_mass, u1_A_array[iU1_A], gr);

    combine_j_sea(nr, +0.5, jE_sea[iU1_A], j_sea, jE_sea[iU1_A]);
    combine_j_sea(nr, +0.5, jO_sea[iU1_A], j_sea, jO_sea[iU1_A]);

    /* Call with opposite charge */
    clear_j_sea(nr, j_sea);
    block_current_stochastic_U1(nr, j_sea, -charge, nsrc, sign, parity, qic,
				M_inv_gr, dM_inv_gr, fn_mass, u1_A_array[iU1_A], gr);

    /* Complete with even charge and odd charge combinations */
    combine_j_sea(nr, +0.5, jE_sea[iU1_A], j_sea, jE_sea[iU1_A]);
    combine_j_sea(nr, -0.5, jO_sea[iU1_A], j_sea, jO_sea[iU1_A]);

  } /* iU1_A */

  destroy_j_sea(j_sea);
  destroy_v_field(dM_inv_gr);

} /* block_current_stochastic_multiU1 */

/************************************************************************/
/* Build the complete M^{-1} gr from the pieces the caller already has.

   M^{-1} = (2m - D)(M^dag M)^{-1}, so on the source parity
   p = 2m (M^dag M)^{-1} gr and on the other parity p = -D (M^dag M)^{-1} gr,
   the latter being exactly the M_inv_gr the caller formed with the
   EVENANDODD dslash.

   NOTE: the 2m term is dropped in the T1 chain, where the single delta-D
   insertion is parity-off-diagonal and the trace is closed on the source
   parity.  In the T2 chain the first insertion moves to the other parity
   and the second comes back, so both parities of M^{-1} gr contribute and
   this term must be kept. */

static void
build_Minv_gr(su3_vector *p, su3_vector *MdagM_inv_gr, su3_vector *M_inv_gr,
	      Real mass, int parity){

  int otherparity = (parity == EVEN) ? ODD : EVEN;
  int i;

  clear_v_field(p);

  FORSOMEFIELDPARITY(i, parity){
    scalar_mult_su3_vector( MdagM_inv_gr+i, 2.*mass, p+i );
  } END_LOOP;

  FORSOMEFIELDPARITY(i, otherparity){
    scalar_mult_su3_vector( M_inv_gr+i, -1., p+i );
  } END_LOOP;
} /* build_Minv_gr */

/************************************************************************/
/* Accumulate the T2 sea-quark loop, Tr[(D^{-1} delta-D)^2], for each of
   the U(1) fields in u1_A_array, at both signs of the charge, forming the
   even and odd combinations of Eq. (16).

   Unlike the T1 version, each (field, charge) here costs a deflated
   EVENANDODD block inversion, not just a dslash.  nU1_A_array is
   therefore expected to be small -- see compute_jhi_sea_T2. */

static void
block_current_stochastic_multiU1_T2( int nr, int nU1_A_array,
				     Real *jE2_sea[], Real *jO2_sea[],
				     Real mass, Real charge, int nsrc,
				     int sign, int parity,
				     quark_invert_control *qic,
				     imp_ferm_links_t *fn_mass,
				     su3_vector *p[], su3_vector *q[], su3_vector *z[],
				     Real *u1_A_array[], su3_vector *gr[]){

  char myname[] = "block_current_stochastic_multiU1_T2";

  Real *j2_sea = create_j_sea(nr);

  for(int iU1_A = 0; iU1_A < nU1_A_array; iU1_A++){

    /* Call with "charge" */
    clear_j_sea(nr, j2_sea);
    block_current_stochastic_U1_T2(nr, j2_sea, charge, nsrc, sign, parity, mass,
				   qic, p, q, z, fn_mass, u1_A_array[iU1_A], gr);

    combine_j_sea(nr, +0.5, jE2_sea[iU1_A], j2_sea, jE2_sea[iU1_A]);
    combine_j_sea(nr, +0.5, jO2_sea[iU1_A], j2_sea, jO2_sea[iU1_A]);

    /* Call with opposite charge.  Both delta-D insertions carry the same
       charge, so the operator is rebuilt, not negated: delta-D(-q) is not
       -delta-D(q). */
    clear_j_sea(nr, j2_sea);
    block_current_stochastic_U1_T2(nr, j2_sea, -charge, nsrc, sign, parity, mass,
				   qic, p, q, z, fn_mass, u1_A_array[iU1_A], gr);

    /* Complete the even and odd charge combinations */
    combine_j_sea(nr, +0.5, jE2_sea[iU1_A], j2_sea, jE2_sea[iU1_A]);
    combine_j_sea(nr, -0.5, jO2_sea[iU1_A], j2_sea, jO2_sea[iU1_A]);

  } /* iU1_A */

  destroy_j_sea(j2_sea);

} /* block_current_stochastic_multiU1_T2 */

#endif
      
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
  /* M_inv_gr = 1/[D^2 + 4*mass1^2] M0_inv_gr */
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
      spin_taste_op_ape_fn(fn01, spin_taste[mu], r_offset, gr_mu, M_inv_gr[is]);
      spin_taste_op_ape_fn(fn01, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
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
      } END_LOOP;
    } /* is, mu */
  destroy_v_field(gr_mu);

  for(int is = 0; is < nsrc; is++)
    destroy_v_field(M_inv_gr[is]);

} /* block_current_stochastic_deltam */
      
/************************************************************************/
/* Calculate the current density difference between fine and sloppy
   solves separately for two combinations of three masses: u, d, s
   using a list of thinned stochastic estimators. Do this for a given
   parity.  In particular, calculate J_u - J_d and J_u - J_s nsrc
   sources in gr.  Results in j_mu[0], j_mu[2]
*/

static void
block_current_stochastic_delta_udus( Real **j_mu[], Real masses[],
				     imp_ferm_links_t *fn[], int nsrc,
				     int sign, int parity, int nr,
				     quark_invert_control qic[],
				     su3_vector *gr[]){
  
  char myname[] = "block_current_stochastic_delta_udus";

  /* Offset for staggered phases in the current definition */
  int r_offset[4] = {0, 0, 0, 0};

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_s = masses[2];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_us = j_mu[2];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_us = fn[2];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_us4 = 4*(m_s*m_s - m_u*m_u);
  quark_invert_control *qic_ud = &qic[0];
  quark_invert_control *qic_us = &qic[2];

  su3_vector *Mud_inv_gr[nsrc];
  su3_vector *Mus_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    Mud_inv_gr[is] = create_v_field();
    Mus_inv_gr[is] = create_v_field();
  }

  int otherparity;
  if(parity == EVEN)
    otherparity=ODD;
  else
    otherparity=EVEN;

  /* Current spin-taste list */
  int *spin_taste = get_spin_taste();

  su3_vector *Mu_inv_gr[nsrc];
  for(int is = 0; is < nsrc; is++){
    Mu_inv_gr[is] = create_v_field();
  }

  qic_ud->parity = parity;
  qic_us->parity = parity;
  /* Mu_inv_gr = 1/[D^2 + 4*m_u^2] gr */
  ks_congrad_block_field(nsrc, gr, Mu_inv_gr, qic_ud, m_u, fn_us);
  /* Mud_inv_gr = 1/[D^2 + 4*m_d^2] Mu_inv_gr */
  ks_congrad_block_field(nsrc, Mu_inv_gr, Mud_inv_gr, qic_ud, m_d, fn_ud);
  /* Mus_inv_gr = 1/[D^2 + 4*m_s^2] Mu_inv_gr */
  ks_congrad_block_field(nsrc, Mu_inv_gr, Mus_inv_gr, qic_us, m_s, fn_us);

  for(int is = 0; is < nsrc; is++){
      destroy_v_field(Mu_inv_gr[is]);
  }
  
  /* Might be better to use a block dslash here?? */
  for(int is = 0; is < nsrc; is++){
    dslash_fn_field( Mus_inv_gr[is], Mus_inv_gr[is], otherparity, fn_us);
    dslash_fn_field( Mud_inv_gr[is], Mud_inv_gr[is], otherparity, fn_ud);
  }
  
  /* For each source, apply current in various directions at the sink */
  su3_vector *gr_us_mu = create_v_field();
  su3_vector *gr_ud_mu = create_v_field();
  for(int is = 0; is < nsrc; is++)
    for(int mu = 0; mu < NMU; mu++){
      
      /* Apply the appropriate spin_taste operator for
	 a nearly conserved current.  */
      spin_taste_op_ape_fn(fn_us, spin_taste[mu], r_offset, gr_us_mu, Mus_inv_gr[is]);
      spin_taste_op_ape_fn(fn_us, spin_taste_index("pion05"), r_offset, gr_us_mu, gr_us_mu);
      spin_taste_op_ape_fn(fn_us, spin_taste[mu], r_offset, gr_ud_mu, Mud_inv_gr[is]);
      spin_taste_op_ape_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_ud_mu, gr_ud_mu);
      
      /* J_mu = -imag[gr * Mus_inv_gr] */
      /* Multiply result by the mass difference squared times four */
      int i;
      int evol = nsrc/nr;  /* Number of thinned source sites */
      int ir = is/evol;    /* Random source block index */
      FORSOMEFIELDPARITY(i, parity){
	complex cc = su3_dot( gr[is]+i, gr_us_mu+i );
	j_mu_us[ir][NMU*i + mu] += -sign*cc.imag*dm2_us4;
	cc = su3_dot( gr[is]+i, gr_ud_mu+i );
	j_mu_ud[ir][NMU*i + mu] += -sign*cc.imag*dm2_ud4;
#if 0
	printf("j_mu src[%d] %d %d %d %d %d %g %g\n", is,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, mu,
	       j_mu_us[ir][NMU*i + mu], j_mu_ud[ir][NMU*i + mu]);
#endif
      } END_LOOP;
    } /* is, mu */
  destroy_v_field(gr_ud_mu);
  destroy_v_field(gr_us_mu);

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(Mud_inv_gr[is]);
    destroy_v_field(Mus_inv_gr[is]);
  }

} /* block_current_stochastic_delta_udus */
      
/************************************************************************/
/* Calculate the current density difference between fine and sloppy
   for two combinations of four masses: u, d, l, s using a list of
   thinned stochastic estimators. Do this for a given parity.  In
   particular, calculate J_u - J_d and J_l - J_s nsrc sources in gr.
   Results in j_mu[0], j_mu[2]
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

  qic_ls->parity = parity;
  int save_deflate = qic_ls->deflate;
  if(parity == ODD)qic_ls->deflate = 0;  // No further deflation for odd parity contribution

  /* Ml_inv_gr = 1/[D^2 + 4*m_l^2] gr */
  ks_congrad_block_field(nsrc, gr, Ml_inv_gr, qic_ls, m_l, fn_ls);
  /* Mud_inv_gr = 1/[D^2 + 4*m_l^2] Ml_inv_gr */
  /* NOTE: we are approximating ud here */
  ks_congrad_block_field(nsrc, Ml_inv_gr, Mud_inv_gr, qic_ls, m_l, fn_ud);
  /* Mls_inv_gr = 1/[D^2 + 4*m_s^2] Ml_inv_gr */
  ks_congrad_block_field(nsrc, Ml_inv_gr, Mls_inv_gr, qic_ls, m_s, fn_ls);
  qic_ls->deflate = save_deflate;

  for(int is = 0; is < nsrc; is++){
      destroy_v_field(Ml_inv_gr[is]);
  }
  
  double onelinktime = -dclock();
  /* Might be better to use a block dslash here?? */
  for(int is = 0; is < nsrc; is++){
    dslash_fn_field( Mls_inv_gr[is], Mls_inv_gr[is], otherparity, fn_ls);
    dslash_fn_field( Mud_inv_gr[is], Mud_inv_gr[is], otherparity, fn_ud);
  }
  onelinktime += dclock();
  node0_printf("Time to do one-link dslash %e\n", onelinktime);
  
  /* For each source, apply current in various directions at the sink */
  su3_vector *gr_ls_mu = create_v_field();
  su3_vector *gr_ud_mu = create_v_field();
  for(int is = 0; is < nsrc; is++)
    for(int mu = 0; mu < NMU; mu++){
      
      /* Apply the appropriate spin_taste operator for
	 a nearly conserved current.  */
      spin_taste_op_ape_fn(fn_ls, spin_taste[mu], r_offset, gr_ls_mu, Mls_inv_gr[is]);
      spin_taste_op_ape_fn(fn_ls, spin_taste_index("pion05"), r_offset, gr_ls_mu, gr_ls_mu);
      spin_taste_op_ape_fn(fn_ls, spin_taste[mu], r_offset, gr_ud_mu, Mud_inv_gr[is]);
      spin_taste_op_ape_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_ud_mu, gr_ud_mu);
      
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
      } END_LOOP;
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
		   int eigstartflags[], char *eigstartfiles[],
		   ks_eigen_param *eigen_param,
		   imp_ferm_links_t *fn_mass[], 
		   quark_invert_control *qic_precise,
		   quark_invert_control *qic_sloppy,
		   int nsrc, int nr, int thinning, su3_vector *gr_even[],
		   su3_vector *gr_odd[], su3_vector *gr0[]){

  char myname[] = "block_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){

    /* Load the appropriate eigenpairs. No provision for saving here. */
    if(eigstartfiles[j] != NULL)
      restore_eigenvectors_from_file(fn_mass[j], eigen_param,
				     eigstartflags[j], eigstartfiles[j], FORGET, NULL);

    int evol = nsrc/nr;  /* Number of thinned source sites */
    for(int jr = 0; jr < nr; jr++){
      collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
      collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
    }

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

#ifdef U1_SEA

/************************************************************************/

/* Calculate the TSM sea quark difference between fine and sloppy
   solves separately for all masses */

static void
compute_jhi_sea_diff(int n_masses, Real **jE_sea[], Real **jO_sea[],
		     Real masses[], Real charges[], imp_ferm_links_t *fn_mass,
		     quark_invert_control *qic_precise,
		     quark_invert_control *qic_sloppy,
		     int nsrc, int nr, int nU1_A_array, Real *u1_A_array[],
		     su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "compute_jhi_sea_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* M_inv_gr is reused for every U(1) field, so the dslash is done out
     of place, into dM_inv_gr */
  su3_vector **MdagM_inv_gr  = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **M_inv_gr  = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    MdagM_inv_gr[is]  = create_v_field();
    M_inv_gr[is]  = create_v_field();
  }

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){
    Real mass = masses[j];
    
    /* First, the sloppy high-mode solution */
    quark_invert_control *qic = qic_sloppy + j;

    /* Calculate (M^\daggerM)^{-1} gr for both parities */
    qic->parity = EVEN;
    ks_congrad_block_field(nsrc, gr_even, MdagM_inv_gr, qic, mass, fn_mass);
    qic->parity = ODD;
    ks_congrad_block_field(nsrc, gr_odd, MdagM_inv_gr, qic, mass, fn_mass);
    /* Calculate M^{-1} gr for both parities */
    for(int is = 0; is < nsrc; is++){
      dslash_fn_field(MdagM_inv_gr[is], M_inv_gr[is] , EVENANDODD, fn_mass);
    }

    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      -1, EVEN, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_even);

    //      node0_printf("jE_sea[%d][0] = %g jO_sea[%d][0] = %g\n",j,jE_sea[j][0],j,jO_sea[j][0]);

    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      -1, ODD, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_odd);
    

    /* Next, continue to a "precise" solution from the same sources */

    /* Calculate (M^\daggerM)^{-1} gr for both parities */
    qic = qic_precise + j;
    qic->parity = EVEN;
    ks_congrad_block_field(nsrc, gr_even, MdagM_inv_gr, qic, mass, fn_mass);
    qic->parity = ODD;
    /* Calculate M^{-1} gr for both parities */
    ks_congrad_block_field(nsrc, gr_odd, MdagM_inv_gr, qic, mass, fn_mass);

    for(int is = 0; is < nsrc; is++){
      dslash_fn_field(MdagM_inv_gr[is], M_inv_gr[is] , EVENANDODD, fn_mass);
    }

    node0_printf("Solving precisely for all EVEN displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      +1, EVEN, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_even);
    
    node0_printf("Solving precisely for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      +1, ODD, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_odd);
  } /* j */
  
  for(int is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(MdagM_inv_gr[is]);
  }
  free(M_inv_gr);
  free(MdagM_inv_gr);

} /* compute_jhi_sea_diff */

#endif

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
   solves separately for four masses, say u d s c  and take differences:
   j_u - j_d, j_u - j_s, and separately, j_c
   Results are in j_mu[0], j_mu[2], and j_mu[3]
*/

static void
block_currents_diff_delta_udusc(int n_masses, Real **j_mu[], Real masses[],
				imp_ferm_links_t *fn[], int nsrc, int nr,
				quark_invert_control qic_precise[],
				quark_invert_control qic_sloppy[],
				su3_vector **gr_even, su3_vector **gr_odd){
  
  char myname[] = "block_current_diff_delta_udusc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_s = masses[2];
  Real m_c = masses[3];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_us = j_mu[2];
  Real **j_mu_c  = j_mu[3];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_us = fn[2];
  imp_ferm_links_t *fn_c  = fn[3];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_us4 = 4*(m_s*m_s - m_u*m_u);
  quark_invert_control *qic_sloppy_ud = &qic_sloppy[0];
  quark_invert_control *qic_sloppy_us = &qic_sloppy[2];
  quark_invert_control *qic_sloppy_c  = &qic_sloppy[3];
  quark_invert_control *qic_precise_ud = &qic_precise[0];
  quark_invert_control *qic_precise_us = &qic_precise[2];
  quark_invert_control *qic_precise_c  = &qic_precise[3];

  /* The sloppy high-mode solution on even sites */
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, -1, EVEN, nr,
				       qic_sloppy, gr_even);
#else
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, -1, EVEN,
				   nr, qic_sloppy_ud, gr_even);
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_s);
  block_current_stochastic_deltam( j_mu_us, m_u, m_s, dm2_us4, fn_us, nsrc, -1, EVEN,
				   nr, qic_sloppy_us, gr_even);
#endif
  
  node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_c);
  block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, EVEN,
			    qic_sloppy_c, fn_c, gr_even);

  /* The sloppy high-mode solution on odd sites */

#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, -1, ODD, nr, 
				       qic_sloppy, gr_odd);
#else
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, -1, ODD,
				   nr, qic_sloppy_ud, gr_odd);
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g\n",
	       m_u, m_s);
  block_current_stochastic_deltam( j_mu_us, m_u, m_s, dm2_us4, fn_us, nsrc, -1, ODD,
				   nr, qic_sloppy_us, gr_odd);
#endif
  node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_c);
  block_current_stochastic( nr, j_mu_c, m_c, nsrc, -1, ODD,
			    qic_sloppy_c, fn_c, gr_odd);
  
  /* The "precise" solution from the same sources on even sites */

#ifdef OPT_UDLSC
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, +1, EVEN, nr, 
				       qic_precise, gr_even);
#else
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, +1, EVEN,
				   nr, qic_precise_ud, gr_even);
  node0_printf("Solving precisely for all EVEN displacements for mass diff %g %g\n",
	       m_u, m_s);
  block_current_stochastic_deltam( j_mu_us, m_u, m_s, dm2_us4, fn_us, nsrc, +1, EVEN,
				   nr, qic_precise_us, gr_even);
#endif
  node0_printf("Solving precisely for all EVEN displacements for mass %g\n", m_c);
  block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			    qic_precise_c, fn_c, gr_even);

  /* The "precise" solution from the same sources on odd sites */

#ifdef OPT_UDLSC
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, +1, ODD, nr, 
				       qic_precise, gr_odd);
#else
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc, +1, ODD,
				   nr, qic_precise_ud, gr_odd);
  node0_printf("Solving precisely for all ODD displacements for mass diff %g %g\n",
	       m_u, m_s);
  block_current_stochastic_deltam( j_mu_us, m_u, m_s, dm2_us4, fn_us, nsrc, +1, ODD,
				   nr, qic_precise_us, gr_odd);
#endif
  node0_printf("Solving precisely for all ODD displacements for mass %g\n", m_c);
  block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, ODD,
			    qic_precise_c, fn_c, gr_odd);
  
} /* block_currents_diff_delta_udusc */

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
  Real **j_mu_s  = j_mu[3];
  Real **j_mu_c  = j_mu[4];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_ls = fn[2];
  imp_ferm_links_t *fn_s  = fn[3];
  imp_ferm_links_t *fn_c  = fn[4];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_sloppy_ud = &qic_sloppy[0];
  quark_invert_control *qic_sloppy_ls = &qic_sloppy[2];
  quark_invert_control *qic_sloppy_s  = &qic_sloppy[3];
  quark_invert_control *qic_sloppy_c  = &qic_sloppy[4];
  quark_invert_control *qic_precise_ud = &qic_precise[0];
  quark_invert_control *qic_precise_ls = &qic_precise[2];
  quark_invert_control *qic_precise_s  = &qic_precise[3];
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
  
  node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, -1, EVEN,
			    qic_sloppy_s, fn_s, gr_even);

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

  node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, -1, ODD,
			    qic_sloppy_s, fn_s, gr_odd);
  
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

  node0_printf("Solving precisely for all EVEN displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, +1, EVEN,
			    qic_precise_s, fn_s, gr_even);

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

  node0_printf("Solving precisely for all ODD displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, +1, ODD,
			    qic_precise_s, fn_s, gr_odd);
  
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
			   int nr, int thinning, su3_vector **gr_even,
			   su3_vector **gr_odd, su3_vector *gr0[]){

  char myname[] = "block_currents_diff_deltam";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  int evol = nsrc/nr;  /* Number of thinned source sites */
  for(int jr = 0; jr < nr; jr++){
    collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
    collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
  }

  if(n_masses <= 3){
    /* Take the difference between the first two masses, but not the third */
    block_currents_diff_delta_lsc(n_masses, j_mu, masses, fn_mass, nsrc, nr,
				  qic_precise, qic_sloppy, gr_even, gr_odd);
  } else if(n_masses == 4) {
    /* Take the difference between the first two masses, and the first and third */
    block_currents_diff_delta_udusc(n_masses, j_mu, masses, fn_mass, nsrc, nr,
				   qic_precise, qic_sloppy, gr_even, gr_odd);
  } else {
    /* Take the difference between the first two masses, and the next two, but not fifth */
    block_currents_diff_delta_udlsc(n_masses, j_mu, masses, fn_mass, nsrc, nr,
				    qic_precise, qic_sloppy, gr_even, gr_odd);
  }
  
} /* block_currents_diff_deltam */

/*********************************************************************/
/* This version of the block current routine does all masses without
   any differences betweem results of different masses.
   Results in j_mu
 */

static void 
block_current( int n_masses, Real **j_mu[], Real masses[],
	       int eigstartflags[], char *eigstartfiles[],
	       ks_eigen_param *eigen_param,
	       imp_ferm_links_t *fn_mass[], quark_invert_control *qic,
	       int nsrc, int nr, int thinning, su3_vector *gr_even[],
	       su3_vector *gr_odd[], su3_vector *gr0[]){

  char myname[] = "block_current";
  // node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){

    /* Load the appropriate eigenpairs. No provision for saving here. */
    if(eigstartfiles[j] != NULL)
      restore_eigenvectors_from_file(fn_mass[j], eigen_param,
				     eigstartflags[j], eigstartfiles[j], FORGET, NULL);

    int evol = nsrc/nr;  /* Number of thinned source sites */
    for(int jr = 0; jr < nr; jr++){
      collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
      collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
    }

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
/* This version takes four masses
   Takes the differences of mass0 and mass1 and of mass0 and mass2
   and computes mass3
   Results in j_mu[0], j_mu[2], and j_mu[3]
*/

static void 
block_currents_delta_udusc( int n_masses, Real **j_mu[], Real masses[], 
			    imp_ferm_links_t *fn[], int nsrc, int nr,
			    quark_invert_control qic[],
			    su3_vector **gr_even, su3_vector **gr_odd){

  char myname[] = "block_currents_delta_udusc";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  Real m_u = masses[0];
  Real m_d = masses[1];
  Real m_s = masses[2];
  Real m_c = masses[3];
  Real **j_mu_ud = j_mu[0];
  Real **j_mu_us = j_mu[2];
  Real **j_mu_c  = j_mu[3];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_us = fn[2];
  imp_ferm_links_t *fn_c  = fn[3];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_us4 = 4*(m_s*m_s - m_u*m_u);
  quark_invert_control *qic_ud = &qic[0];
  quark_invert_control *qic_us = &qic[2];
  quark_invert_control *qic_c  = &qic[3];

  /* Construct current density from the list of sources */

#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, +1, EVEN, nr, 
				       qic, gr_even);
#else
  node0_printf("Solving sloppily for all EVEN displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, EVEN, nr, qic_ud, gr_even);
  node0_printf("Solving sloppily for all EVEN displacements for mass diff  %g %g\n",
	       m_u, m_s);
  block_current_stochastic_deltam( j_mu_us, m_u, m_s, dm2_us4, fn_us, nsrc,
				   +1, EVEN, nr, qic_us, gr_even);
#endif
  if(n_masses == 4){
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_c, fn_c, gr_even);
  }
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_u, m_s);
  block_current_stochastic_delta_udus( j_mu, masses, fn, nsrc, +1, ODD, nr, 
				       qic, gr_odd);
#else
  node0_printf("Solving sloppily for all ODD displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, ODD, nr, qic_ud, gr_odd);
  node0_printf("Solving sloppily for all ODD displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_us, m_l, m_s, dm2_us4, fn_us, nsrc,
				   +1, ODD, nr, qic_us, gr_odd);
#endif
  if(n_masses == 5){
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_c);
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
  Real **j_mu_s  = j_mu[3];
  Real **j_mu_c  = j_mu[4];
  imp_ferm_links_t *fn_ud = fn[0];
  imp_ferm_links_t *fn_ls = fn[2];
  imp_ferm_links_t *fn_s  = fn[3];
  imp_ferm_links_t *fn_c  = fn[4];
  Real dm2_ud4 = 4*(m_d*m_d - m_u*m_u);
  Real dm2_ls4 = 4*(m_s*m_s - m_l*m_l);
  quark_invert_control *qic_ud = &qic[0];
  quark_invert_control *qic_ls = &qic[2];
  quark_invert_control *qic_s  = &qic[3];
  quark_invert_control *qic_c  = &qic[4];

  /* Construct current density from the list of sources */
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all EVEN displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, EVEN, nr, 
				       qic, gr_even);
#else
  node0_printf("Solving sloppily for all EVEN displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, EVEN, nr, qic_ud, gr_even);
  node0_printf("Solving sloppily for all EVEN displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc,
				   +1, EVEN, nr, qic_ls, gr_even);
#endif

  node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, +1, EVEN,
			    qic_s, fn_s, gr_even);
  
  if(n_masses == 5){
    node0_printf("Solving sloppily for all EVEN displacements for mass %g\n", m_c);
    block_current_stochastic( nr, j_mu_c, m_c, nsrc, +1, EVEN,
			      qic_c, fn_c, gr_even);
  }
  
#ifdef OPT_UDLSC
  node0_printf("Solving sloppily for all ODD displacements for mass diff %g %g and %g %g\n",
	       m_u, m_d, m_l, m_s);
  block_current_stochastic_delta_udls( j_mu, masses, fn, nsrc, +1, ODD, nr, 
				       qic, gr_odd);
#else
  node0_printf("Solving sloppily for all ODD displacements for mass diff  %g %g\n",
	       m_u, m_d);
  block_current_stochastic_deltam( j_mu_ud, m_u, m_d, dm2_ud4, fn_ud, nsrc,
				   +1, ODD, nr, qic_ud, gr_odd);
  node0_printf("Solving sloppily for all ODD displacements for mass diff  %g %g\n",
	       m_l, m_s);
  block_current_stochastic_deltam( j_mu_ls, m_l, m_s, dm2_ls4, fn_ls, nsrc,
				   +1, ODD, nr, qic_ls, gr_odd);
#endif

  node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_s);
  block_current_stochastic( nr, j_mu_s, m_s, nsrc, +1, ODD,
			    qic_s, fn_s, gr_odd);

  if(n_masses == 5){
    node0_printf("Solving sloppily for all ODD displacements for mass %g\n", m_c);
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
		       int nsrc, int nr, int thinning, su3_vector *gr_even[],
		       su3_vector *gr_odd[], su3_vector *gr0[]){

  char myname[] = "block_currents_deltam";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  int evol = nsrc/nr;  /* Number of thinned source sites */
  for(int jr = 0; jr < nr; jr++){
    collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
    collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
  }
  
  if(n_masses <= 3){
    block_currents_delta_lsc(n_masses, j_mu, masses, fn, nsrc, nr,
			     qic, gr_even, gr_odd);
  } else if(n_masses == 4) {
    block_currents_delta_udusc(n_masses, j_mu, masses, fn, nsrc, nr,
			       qic, gr_even, gr_odd);
  } else {
    block_currents_delta_udlsc(n_masses, j_mu, masses, fn, nsrc, nr,
			       qic, gr_even, gr_odd);
  }      
}

#ifdef HAVE_QUDA

/*********************************************************************/
/* Calculate exact low-mode current densities using QUDA.
 *
 * For nmass=1, this function emulates MILC's exact_current()
 * For nmass=2, this function emulates MILC's exact_current_delta_ls()
 * For nmass=3, this function emulates MILC's exact_current_delta_udus
 *
 * Both parity eigenvectors must be loaded into QUDA before calling this function.
 */
static void
exact_current_quda(Real *jlow_mu1, Real *jlow_mu2, int nmass, Real masses[], imp_ferm_links_t *fn_mass){

  char myname[] = "exact_current_quda";

  if( nmass<1 || nmass>3 ) {
    node0_printf("%s: wrong number of masses %d\n", myname, nmass);
    terminate(1);
  }

  int refresh  = 0;
  if (fn_mass != get_fn_last() || fresh_fn_links(fn_mass)){
    cancel_quda_notification(fn_mass);
    set_fn_last(fn_mass);
    refresh = 1;
    node0_printf("%s: fn, notify: Signal QUDA to refresh links\n", __func__);
  }

  QudaInvertArgs_t inv_args; // To set invertParams and gaugeParams in QUDA
  inv_args.mixed_precision = 0;
  inv_args.naik_epsilon = fn_mass->eps_naik;
#if (FERM_ACTION==HISQ)
  inv_args.tadpole = 1.0;
#else
  inv_args.tadpole = u0;
#endif

  QudaEigensolverArgs_t eig_args; // To tell QUDA how many eigenvectors to use

#ifdef USE_EIG_GPU
  // Here we use QUDA for the eigensolution or for reading is own eigenvector file
  int quda_does_eigensolve = (param.ks_eigen_startflag == FRESH);
#else
  // Here, eigenvectors were loaded from file(s) by MILC or by a non-QUDA eigensolver
  int quda_does_eigensolve = 0;
#endif

  load_quda_default_eig_args(&eig_args, quda_does_eigensolve);
  eig_args.n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
  eig_args.n_ev = eig_args.n_conv;

  su3_matrix* fatlink = get_fatlinks(fn_mass);
  su3_matrix* longlink = get_lnglinks(fn_mass);

  // Compute exact current via QUDA
  qudaExactCurrent(MILC_PRECISION, MILC_PRECISION, fatlink, longlink, ape_links, nmass, masses, inv_args, eig_args, jlow_mu1, jlow_mu2, refresh);

} // exact_current_quda

#endif

/*********************************************************************/

static void
exact_current(Real *jlow_mu, Real mass, imp_ferm_links_t *fn_mass){

  // node0_printf("Entered exact_current\n"); fflush(stdout);

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
      spin_taste_op_ape_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_ape_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      FOREVENFIELDSITES(i){
	      complex z;
	      z = su3_dot( eigVec[n] + i, gr_mu + i);
	      jlow_mu[NMU*i + mu] += -z.imag/(eigVal[n]+4.0*mass*mass);
      } /* i */
      spin_taste_op_ape_fn(fn_mass, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_ape_fn(fn_mass, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
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

#ifdef HAVE_QUDA

/*********************************************************************/
/* Calculate exact low-mode sea-quark loop using QUDA.
 *
 * For nmass=1, this function emulates MILC's exact_current()
 * For nmass=2, this function emulates MILC's exact_current_delta_ls()
 * For nmass=3, this function emulates MILC's exact_current_delta_udus
 *
 * Both parity eigenvectors must be loaded into QUDA before calling this function.
 */
static void
exact_current_sea_quda(Real *jlow_sea, Real mass, imp_ferm_links_t *fn_diff, imp_ferm_links_t *fn_mass){

  char myname[] = "exact_current_sea_quda";

  QudaEigensolverArgs_t eig_args; // To tell QUDA how many eigenvectors to use

#ifdef USE_EIG_GPU
  // Here we use QUDA for the eigensolution or for reading is own eigenvector file
  int quda_does_eigensolve = (param.ks_eigen_startflag == FRESH);
#else
  // Here, eigenvectors were loaded from file(s) by MILC or by a non-QUDA eigensolver
  int quda_does_eigensolve = 0;
#endif

  load_quda_default_eig_args(&eig_args, quda_does_eigensolve);
  eig_args.n_conv = (param.eigen_param.Nvecs_in > param.eigen_param.Nvecs) ? param.eigen_param.Nvecs_in : param.eigen_param.Nvecs;
  eig_args.n_ev = eig_args.n_conv;

  su3_matrix* dfatlink = get_fatlinks(fn_diff);
  su3_matrix* dlonglink = get_lnglinks(fn_diff);

  su3_matrix* fatlink = get_fatlinks(fn_mass);
  su3_matrix* longlink = get_lnglinks(fn_mass);

  // Compute exact current via QUDA
  qudaExactSeaLoop(MILC_PRECISION, MILC_PRECISION, dfatlink, dlonglink, fatlink, longlink,
		   eig_args, jlow_sea);

} // exact_current_quda

#endif

/*********************************************************************/
static void
exact_current_sea(Real *jlow_sea, Real mass, imp_ferm_links_t *fn_diff, imp_ferm_links_t *fn_mass){

  /* Compute exact sea-quark loop for a single mass */

  su3_vector *dslashV = create_v_field();
  su3_vector *ddslashV = create_v_field();
  int Nvecs = param.eigen_param.Nvecs;
  
  for(int n = 0; n < Nvecs; n++){
     dslash_fn_field(eigVec[n], dslashV, ODD, fn_mass);
     dslash_fn_field(dslashV, ddslashV, EVEN, fn_diff);
     int i;
     FOREVENFIELDSITES(i){
       complex z;
       z = su3_dot( eigVec[n] + i, ddslashV + i);
       *jlow_sea += -2*z.real/(eigVal[n]+4.0*mass*mass);
     }
  } /* n */
  
  destroy_v_field(dslashV);
  destroy_v_field(ddslashV);
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
      
      spin_taste_op_ape_fn(fn_ud, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_ape_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	Real dl = eigVal[n]+4.0*m_l*m_l;
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_ud[NMU*i + mu] += -z.imag*dm2_ls4/(dl*ds);
      } /* i */
      
      spin_taste_op_ape_fn(fn_ud, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_ape_fn(fn_ud, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
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
   jlow_mu[2] = u - s
   The third mass "l" is ignored
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
      
      spin_taste_op_ape_fn(fn_udls, spin_taste[mu], r_offset, gr_mu, gr0);
      spin_taste_op_ape_fn(fn_udls, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
      FOREVENFIELDSITES(i){
	complex z;
	z = su3_dot( eigVec[n] + i, gr_mu + i);
	Real du = eigVal[n]+4.0*m_u*m_u;
	Real dd = eigVal[n]+4.0*m_d*m_d;
	jlow_mu_ud[NMU*i + mu] += -z.imag*dm2_ud4/(du*dd);
	Real ds = eigVal[n]+4.0*m_s*m_s;
	jlow_mu_us[NMU*i + mu] += -z.imag*dm2_us4/(du*ds);
      } /* i */
      
      spin_taste_op_ape_fn(fn_udls, spin_taste[mu], r_offset, gr_mu, eigVec[n]);
      spin_taste_op_ape_fn(fn_udls, spin_taste_index("pion05"), r_offset, gr_mu, gr_mu);
      
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
	       int eigstartflags[], char *eigstartfiles[],
	       ks_eigen_param *eigen_param,
	       imp_ferm_links_t *fn_mass[]){

  // node0_printf("Entered exact_currents\n"); fflush(stdout);

  /* Compute exact low-mode current density for all masses 
     No differences taken. 
     Results in jlow_mu
  */

  double dtime = -dclock();

  for(int j = 0; j < n_masses; j++){

    /* Load the appropriate eigenpairs. No provision for saving here. */
    if(eigstartfiles[j] != NULL)
      restore_eigenvectors_from_file(fn_mass[j], eigen_param,
				     eigstartflags[j], eigstartfiles[j], FORGET, NULL);
    
#if defined(HAVE_QUDA) && defined(USE_CURRENT_GPU)

    Real qmasses[1] = {masses[j]};
    exact_current_quda(jlow_mu[j], NULL, 1, qmasses, fn_mass[j]);

#else

    exact_current(jlow_mu[j], masses[j], fn_mass[j]);

#endif
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
exact_currents_sea(Real mass, Real *jlow_sea, imp_ferm_links_t *fn_diff,
		   imp_ferm_links_t *fn_mass){

  /* Compute exact low-mode sea-quark loop with the insertion of the
     dslash difference operator */
    
#if defined(HAVE_QUDA) && defined(USE_CURRENT_GPU)

  exact_current_sea_quda(jlow_sea, mass, fn_diff, fn_mass);

#else

  exact_current_sea(jlow_sea, mass, fn_diff, fn_mass);

#endif

#if 0
  node0_printf("Exact low mode for mass %g jlow_sea = %g\n", mass, jlow_sea);
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

#if defined(HAVE_QUDA) && defined(USE_CURRENT_GPU)
  Real quda_m[3] = {masses[0], masses[1], masses[3]};
  Real quda_m1[1] = {masses[0]};
  switch(n_masses){
  case(1):
    exact_current_quda(jlow_mu[0], NULL, 1, quda_m1, fn_mass[0]);
    break;
  case(2):
    exact_current_quda(jlow_mu[0], NULL, 2, quda_m, fn_mass[0]);
    break;
  case(3):
    quda_m1[0] = masses[2];
    exact_current_quda(jlow_mu[0], NULL, 2, quda_m, fn_mass[0]);
    exact_current_quda(jlow_mu[2], NULL, 1, quda_m1, fn_mass[2]);
    break;
  case(4):
    exact_current_quda(jlow_mu[0], jlow_mu[2], 3, quda_m, fn_mass[0]);
    break;
  case(5):
    exact_current_quda(jlow_mu[0], jlow_mu[2], 3, quda_m, fn_mass[0]);
    quda_m1[0] = masses[3];
    exact_current_quda(jlow_mu[3], NULL, 1, quda_m1, fn_mass[3]);
    quda_m1[0] = masses[4];
    exact_current_quda(jlow_mu[4], NULL, 1, quda_m1, fn_mass[4]);
    break;
  default:
    node0_printf("%s: wrong number of masses %d\n", myname, n_masses);
    terminate(1);
  }
#else
  Real m[5] = {masses[0], masses[1], 0., masses[2], masses[3]};
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
    exact_current_delta_udus(jlow_mu, m, fn_mass);
    break;
  case(5):
    exact_current_delta_udus(jlow_mu, masses, fn_mass);
    exact_current(jlow_mu[3], masses[3], fn_mass[3]);
    exact_current(jlow_mu[4], masses[4], fn_mass[4]);
    break;
  default:
    node0_printf("%s: wrong number of masses %d\n", myname, n_masses);
    terminate(1);
  }
#endif
  
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

#if 0
/*********************************************************************/
static void
check_eigen(int Nvecs){
  /* Check orthonormality of a few eigenvectors (for debugging) */
  for(int j = 0; j < Nvecs; j += 8)
    for(int i = j; i < Nvecs; i += 8){
      double_complex cc ;
      dot_product(eigVec[i], eigVec[j], &cc, EVEN) ;
      if(((i == j) && (fabs(cc.real - 1) > 1e-8)) || ((i != j && fabs(cc.real) > 1e-8)))
	node0_printf("vec[%d] * vec[%d] = %g %g\n", i, j, cc.real, cc.imag);
    }
}

#endif

/*********************************************************************/
/* Create fields for low- and high-mode current densities */
static void
create_jlow(int n_masses, int Nvecs, Real *jlow_mu[]){
  for(int j = 0; j < n_masses; j++){
    if(Nvecs > 0)
      jlow_mu[j] = create_r_array_field(NMU);
    else
      jlow_mu[j] = NULL;
  }
}

/*********************************************************************/
static void
destroy_jlow(int Nvecs, int n_masses, Real *jlow_mu[]){
  if(Nvecs > 0)
    for(int j = 0; j < n_masses; j++)
      destroy_r_array_field(jlow_mu[j], NMU);
}

/*********************************************************************/
/* Create fields for low-mode current densities */
static void
clear_jlow_sea(int n_masses, int nU1_A, Real *jlow_sea[]){
  for(int j = 0; j < n_masses; j++)
    for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
      jlow_sea[j][iU1_A] = 0.;
}

/*********************************************************************/
/* Create fields for low-mode current densities */
static void
create_jlow_sea(int n_masses, int nU1_A, Real *jlow_sea[]){
  for(int j = 0; j < n_masses; j++){
    jlow_sea[j] = (Real *)malloc(nU1_A*sizeof(Real));
  }
  clear_jlow_sea(n_masses, nU1_A, jlow_sea);
}

/*********************************************************************/
static void
destroy_jlow_sea(int n_masses, Real *jlow_sea[]){
  for(int j = 0; j < n_masses; j++)
    free(jlow_sea[j]);
}

/*********************************************************************/
/* Create fields for high-mode current densities, one for each
   block of nr random sources and mass */
static void
create_jhi(int n_masses, int nr, Real **j_mu[]){
  for(int j = 0; j < n_masses; j++){
    j_mu[j] = (Real **)malloc(nr*sizeof(Real *));
    for(int ir = 0; ir < nr; ir++)
      j_mu[j][ir] = create_r_array_field(NMU);
  }
}

/*********************************************************************/
static void
destroy_jhi(int n_masses, int nr, Real **j_mu[]){
  for(int j = 0; j < n_masses; j++){
    for(int ir = 0; ir < nr; ir++)
      destroy_r_array_field(j_mu[j][ir], NMU);
    free(j_mu[j]);
  }
}

/*********************************************************************/
/* Create fields for high-mode sea-quark delta D, one for each
   block of nr random sources and mass */
static void
clear_jhi_sea(int n_masses, int nU1_A, int nr, Real **j_sea[]){
  for(int j = 0; j < n_masses; j++){
    for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
      for(int ir = 0; ir < nr; ir++)
	j_sea[j][iU1_A][ir] = 0.;
  }
}

/*********************************************************************/
/* Create fields for high-mode sea-quark delta D, one for each
   block of nr random sources and mass */
static void
create_jhi_sea(int n_masses, int nU1_A, int nr, Real **j_sea[]){
  for(int j = 0; j < n_masses; j++){
    j_sea[j] = (Real **)malloc(nU1_A*sizeof(Real *));
    for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
      j_sea[j][iU1_A] = (Real *)malloc(nr*sizeof(Real));
  }
  clear_jhi_sea(n_masses, nU1_A, nr, j_sea);
}

/*********************************************************************/
static void
destroy_jhi_sea(int n_masses, int nU1_A, Real **j_sea[]){
  for(int j = 0; j < n_masses; j++){
    for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
      free(j_sea[j][iU1_A]);
  }
}

/*********************************************************************/
/* Compute exact low-mode current density if we have eigenvectors to do it */
static void
compute_jlow(int Nvecs, int n_masses, Real *jlow_mu[], Real masses[],
	     int eigstartflags[], char *eigstartfiles[],
	     ks_eigen_param *eigen_param, imp_ferm_links_t *fn_mass[]){

#if defined(MASS_UDLSC) && !defined(U1_VALENCE)
  if(Nvecs > 0){
    if(n_masses <= 1){
      exact_currents(n_masses, jlow_mu, masses, eigstartflags, eigstartfiles, eigen_param, fn_mass);
    } else {
      exact_currents_deltam(n_masses, jlow_mu, masses, fn_mass);
    }
  }
#else
  if(Nvecs > 0){
    exact_currents(n_masses, jlow_mu, masses, eigstartflags,
		   eigstartfiles, eigen_param, fn_mass);
  }
#endif
}

#ifdef U1_SEA

/*********************************************************************/
/* Compute exact low-mode sea-quark loop */
static void
compute_jlow_sea(int n_masses, Real *jElow_sea[], Real *jOlow_sea[],
		 Real masses[], Real charges[], int nU1_A_array, Real *u1_A_array[],
		 imp_ferm_links_t *fn_mass){

  double dtime = -dclock();

  for(int j = 0; j < n_masses; j++)
    for(int iU1_A = 0; iU1_A < nU1_A_array; iU1_A++){

      /* Create Dslash difference operator */

      /* Calculate for both signs of the charge and take even and odd combinations */

      imp_ferm_links_t *fn_diff = create_dslash_diff(u1_A_array[iU1_A], charges[j], fn_mass);
      Real jlow_sea = 0.;
      exact_currents_sea(masses[j], &jlow_sea, fn_diff, fn_mass);
      jElow_sea[j][iU1_A] += 0.5*jlow_sea;
      jOlow_sea[j][iU1_A] += 0.5*jlow_sea;
      destroy_fn_links(fn_diff);

      fn_diff = create_dslash_diff(u1_A_array[iU1_A], -charges[j], fn_mass);
      jlow_sea = 0.;
      exact_currents_sea(masses[j], &jlow_sea, fn_diff, fn_mass);
      jElow_sea[j][iU1_A] += 0.5*jlow_sea;
      jOlow_sea[j][iU1_A] -= 0.5*jlow_sea;

      destroy_fn_links(fn_diff);
    }

  dtime += dclock();

#ifdef CGTIME
  node0_printf("Time for exact low modes %g sec\n", dtime);
#endif

}

#endif

/**********************************************************************************/
/* Print the exact low mode contribution to the current density */
static void
write_jlow(int n_masses, Real masses[], Real charges[], Real *jlow_mu[]){
  
  write_tslice_values_begin("LOW");

#ifdef MASS_UDLSC
  switch(n_masses){
    
  case(1): /* u */
    write_tslice_values("LOW", -1, masses[0], charges[0], 0., 0., jlow_mu[0]);
    break;
    
  case(2): /* l - s */
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0]);
    break;
    
  case(3):  /* l - s,  c */
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0]);
    write_tslice_values("LOW", -1, masses[2], charges[2], 0., 0., jlow_mu[2]);
    break;
    
  case(4):  /* u - d,  u - s */
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0]);
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[2], charges[2], jlow_mu[2]);
    break;
    
  case(5):  /* u - d,  u - s,  s, c */
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0]);
    write_tslice_values("LOW", -1, masses[0], charges[0], masses[3], charges[3], jlow_mu[2]);
    write_tslice_values("LOW", -1, masses[3], charges[3], 0., 0., jlow_mu[3]);
    write_tslice_values("LOW", -1, masses[4], charges[4], 0., 0., jlow_mu[4]);
    break;
  }
#else
  for(int j = 0; j < n_masses; j++)
    write_tslice_values("LOW", -1, masses[j], charges[j], 0., 0., jlow_mu[j]);
#endif

  write_tslice_values_end("LOW");
}

/**********************************************************************************/
/* Print the exact low mode contribution to the current density */
static void
write_jlow_sea(char tag[], int n_masses, int nU1_A, 
	       Real masses[], Real charges[], Real *jlow_sea[]){
  
  write_sea_values_begin(tag);

  for(int j = 0; j < n_masses; j++){
    write_sea_values_group(tag, masses[j], charges[j]);
    for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
      write_sea_values(tag, -1, iU1_A, jlow_sea[j][iU1_A]);
  }

  write_sea_values_end(tag);
}

/***********************************************************************************/
/* Reset jlow */
static void
clear_jlow(int n_masses, Real *jlow_mu[]){
  for(int j = 0; j < n_masses; j++){
    /* Reset jlow_mu */
    if(jlow_mu[j] != NULL)clear_r_array_field(jlow_mu[j], NMU);
  }  
}

/***********************************************************************************/
/* Print the exact low mode contribution to the QED current loop */
static void
write_jlowdotA(int n_masses, Real masses[], Real charges[],
	       Real *jlow_mu[], Real *u1_A){

#ifdef MASS_UDLSC
  switch(n_masses){
    
  case(1): /* u */
    write_jdotA_value("LOW", -1, masses[0], charges[0], 0., 0., jlow_mu[0], u1_A);
    break;
    
  case(2): /* l - s */
    write_jdotA_value("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0], u1_A);
    break;
    
  case(3):  /* l - s,  c */
    write_jdotA_value("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0], u1_A);
    write_jdotA_value("LOW", -1, masses[2], charges[2], 0., 0., jlow_mu[2], u1_A);
    break;
    
  case(4):  /* u - d,  u - s */
    write_jdotA_value("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0], u1_A);
    write_jdotA_value("LOW", -1, masses[0], charges[0], masses[2], charges[2], jlow_mu[2], u1_A);
    break;
    
  case(5):  /* u - d,  l - s,  c */
    write_jdotA_value("LOW", -1, masses[0], charges[0], masses[1], charges[1], jlow_mu[0], u1_A);
    write_jdotA_value("LOW", -1, masses[2], charges[2], masses[3], charges[3], jlow_mu[2], u1_A);
    write_jdotA_value("LOW", -1, masses[4], charges[4], 0., 0., jlow_mu[4], u1_A);
    break;
  }

#else
  for(int j = 0; j < n_masses; j++)
    write_jdotA_value("LOW", -1, masses[j], charges[j], 0., 0., jlow_mu[j], u1_A);
#endif
}

/*********************************************************************/
static su3_vector **
create_gr(int nsrc){
  su3_vector **gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++)
    gr[is] = create_v_field();
  return gr;
}

/*********************************************************************/
static void
destroy_gr(int nsrc, su3_vector **gr){
  for(int is = 0; is < nsrc; is++)
    destroy_v_field(gr[is]);
  free(gr);
}

/*********************************************************************/
static Real **
create_U1_A(int nU1_A){
  Real **u1_A = (Real **)malloc(nU1_A*sizeof(Real *));
  for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
    u1_A[iU1_A] = create_r_field();
  return u1_A;
}

/*********************************************************************/
static void
destroy_U1_A(int nU1_A, Real **u1_A){
  for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
    destroy_r_field(u1_A[iU1_A]);
  free(u1_A);
}

/*********************************************************************/
/* Do block solve for high-mode densities */
static void
compute_jhi(int n_masses, Real **j_mu[], Real masses[],
	    int eigstartflags[], char *eigstartfiles[],
	    ks_eigen_param *eigen_param,
	    imp_ferm_links_t *fn_mass[],
	    quark_invert_control *qic, int nsrc, int nr, int thinning,
	    su3_vector **gr_even, su3_vector **gr_odd, su3_vector *gr0[]){    

#if defined(MASS_UDLSC) && !defined(U1_VALENCE)
  if(n_masses <= 1){
    char *null[1] = {NULL};
    block_current( n_masses, j_mu, masses, NULL, null, NULL,
		   fn_mass, qic, nsrc, nr, thinning, gr_even, gr_odd, gr0);
  } else {
    block_currents_deltam( n_masses, j_mu, masses, fn_mass, qic,
			   nsrc, nr, thinning, gr_even, gr_odd, gr0);
  }
#else
  block_current( n_masses, j_mu, masses, eigstartflags, eigstartfiles, eigen_param,
		 fn_mass, qic, nsrc, nr, thinning, gr_even, gr_odd, gr0);
#endif
}

#ifdef U1_SEA

/************************************************************************/

/* Stochastic estimate of T1 = Tr[(D[U]^{-1} delta-D)], Eq. (18). */

static void
compute_jhi_sea(int n_masses, Real **jE_sea[], Real **jO_sea[],
		Real masses[], Real charges[], imp_ferm_links_t *fn_mass,
		quark_invert_control *qic_set,
		int nsrc, int nr, int nU1_A_array, Real *u1_A_array[],
		su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "compute_jhi_sea";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* M_inv_gr is reused for every U(1) field, so the dslash is done out
     of place, into dM_inv_gr */
  su3_vector **MdagM_inv_gr  = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **M_inv_gr  = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    MdagM_inv_gr[is]  = create_v_field();
    M_inv_gr[is]  = create_v_field();
  }

  /* Construct current density from the list of sources */
  for(int j = 0; j < n_masses; j++){
    Real mass = masses[j];
    quark_invert_control *qic = qic_set + j;

    /* Calculate (M^\daggerM)^{-1} gr for both parities */
    qic->parity = EVEN;
    ks_congrad_block_field(nsrc, gr_even, MdagM_inv_gr, qic, mass, fn_mass);
    qic->parity = ODD;
    ks_congrad_block_field(nsrc, gr_odd, MdagM_inv_gr, qic, mass, fn_mass);
    /* Calculate M^{-1} gr for both parities */
    for(int is = 0; is < nsrc; is++){
      dslash_fn_field(MdagM_inv_gr[is], M_inv_gr[is] , EVENANDODD, fn_mass);
    }

    /* The high-mode solution */
    node0_printf("Solving for all EVEN displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      +1, EVEN, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_even);

    //      node0_printf("jE_sea[%d][0] = %g jO_sea[%d][0] = %g\n",j,jE_sea[j][0],j,jO_sea[j][0]);

    node0_printf("Solving for all ODD displacements for mass %g\n", masses[j]);
    block_current_stochastic_multiU1( nr, nU1_A_array, jE_sea[j], jO_sea[j],
				      charges[j], nsrc,
				      +1, ODD, qic, fn_mass,
				      M_inv_gr, u1_A_array, gr_odd);
  } /* j */
  
  for(int is = 0; is < nsrc; is++){
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(MdagM_inv_gr[is]);
  }
  free(M_inv_gr);
  free(MdagM_inv_gr);

} /* compute_jhi_sea */


/************************************************************************/
/* Stochastic estimate of T2 = Tr[(D[U]^{-1} delta-D)^2], Eq. (18).

   The first solve, (M^dag M)^{-1} gr, is shared by every field and charge
   and is done once per mass here.  The sequential solve inside
   block_current_stochastic_U1_T2 is not: it costs one deflated EVENANDODD
   block inversion per field, per charge, per source parity.  With n fields
   that is 4n block solves per mass per source block, against 2 for T1, so
   the caller should pass a much smaller field count than it does for T1.
   That is legitimate: diagram (c) is <C_QCD T2_even>_c and C_QCD carries no
   photon, so the photon average factorises and the fields used here need
   not be the ones used for diagrams (b) and (d).

   The sources are expected to be UNPROJECTED.  Eq. (18) uses the full
   D^{-1} on both sides; with the low modes projected out of gr this
   estimates Tr[X(1-P)] instead, and the missing Tr[XP] would cost one
   solve per eigenvector.  Solver deflation (qic->deflate) is independent
   of source projection and should stay on. */

static void
compute_jhi_sea_T2(int n_masses, Real **jE2_sea[], Real **jO2_sea[],
		   Real masses[], Real charges[], imp_ferm_links_t *fn_mass,
		   quark_invert_control *qic,
		   int nsrc, int nr, int nU1_A_array, Real *u1_A_array[],
		   su3_vector *gr_even[], su3_vector *gr_odd[]){

  char myname[] = "compute_jhi_sea_T2";

  /* Five arrays of nsrc fields.  nsrc = (thinning^4/2)*nr, so this grows
     fast with the source spacing. */
  su3_vector **MdagM_inv_gr = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **M_inv_gr     = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **p            = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **q            = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  su3_vector **z            = (su3_vector **)malloc(nsrc*sizeof(su3_vector *));
  for(int is = 0; is < nsrc; is++){
    MdagM_inv_gr[is] = create_v_field();
    M_inv_gr[is]     = create_v_field();
    p[is]            = create_v_field();
    q[is]            = create_v_field();
    z[is]            = create_v_field();
  }

  for(int j = 0; j < n_masses; j++){

    Real mass = masses[j];
    quark_invert_control *qicj = qic + j;

    /* The shared first solve.  The EVEN and ODD solves write disjoint
       sites of MdagM_inv_gr, and the EVENANDODD dslash then gives
       M_inv_gr|even = D_eo (M^dag M)^{-1} gr_odd and
       M_inv_gr|odd  = D_oe (M^dag M)^{-1} gr_even, so each source parity
       reads back only its own solve. */
    qicj->parity = EVEN;
    ks_congrad_block_field(nsrc, gr_even, MdagM_inv_gr, qicj, mass, fn_mass);
    qicj->parity = ODD;
    ks_congrad_block_field(nsrc, gr_odd,  MdagM_inv_gr, qicj, mass, fn_mass);
    for(int is = 0; is < nsrc; is++)
      dslash_fn_field(MdagM_inv_gr[is], M_inv_gr[is], EVENANDODD, fn_mass);

    /* EVEN sources */
    node0_printf("%s: T2 chain for all EVEN displacements for mass %g\n",
		 myname, mass);
    for(int is = 0; is < nsrc; is++)
      build_Minv_gr(p[is], MdagM_inv_gr[is], M_inv_gr[is], mass, EVEN);

    block_current_stochastic_multiU1_T2( nr, nU1_A_array, jE2_sea[j], jO2_sea[j],
					 mass, charges[j], nsrc,
					 +1, EVEN, qicj, fn_mass,
					 p, q, z, u1_A_array, gr_even);

    /* ODD sources */
    node0_printf("%s: T2 chain for all ODD displacements for mass %g\n",
		 myname, mass);
    for(int is = 0; is < nsrc; is++)
      build_Minv_gr(p[is], MdagM_inv_gr[is], M_inv_gr[is], mass, ODD);

    block_current_stochastic_multiU1_T2( nr, nU1_A_array, jE2_sea[j], jO2_sea[j],
					 mass, charges[j], nsrc,
					 +1, ODD, qicj, fn_mass,
					 p, q, z, u1_A_array, gr_odd);
  } /* j */

  for(int is = 0; is < nsrc; is++){
    destroy_v_field(z[is]);
    destroy_v_field(q[is]);
    destroy_v_field(p[is]);
    destroy_v_field(M_inv_gr[is]);
    destroy_v_field(MdagM_inv_gr[is]);
  }
  free(z); free(q); free(p); free(M_inv_gr); free(MdagM_inv_gr);

} /* compute_jhi_sea_T2 */


#endif

/*********************************************************************/
/* Do block solve for the TSM difference in high-mode densities */
static void
compute_jhi_diff(int n_masses, Real **j_mu[], Real masses[],
		 int eigstartflags[], char *eigstartfiles[], ks_eigen_param *eigen_param,
		 imp_ferm_links_t *fn_mass[],
		 quark_invert_control *qic_precise, quark_invert_control *qic_sloppy,
		 int nsrc, int nr, int thinning, su3_vector **gr_even, su3_vector **gr_odd, su3_vector *gr0[]){
#if defined(MASS_UDLSC) && !defined(U1_VALENCE)
  if(n_masses <= 1){
    block_current_diff(n_masses, j_mu, masses, NULL, NULL, NULL, fn_mass, 
		       qic_precise, qic_sloppy, nsrc, nr, thinning, gr_even, gr_odd, gr0);
  } else {
    block_currents_diff_deltam(n_masses, j_mu, masses, fn_mass, 
			       qic_precise, qic_sloppy, nsrc, nr, thinning, gr_even, gr_odd, gr0);
  }
#else
  block_current_diff(n_masses, j_mu, masses,  eigstartflags, eigstartfiles, eigen_param, fn_mass, 
		     qic_precise, qic_sloppy, nsrc, nr, thinning, gr_even, gr_odd, gr0);
#endif
  
}

/*********************************************************************/
/* Write the high-mode result for this block of nr sources and clear j_mu */
static void
write_jhi(char tag[], int n_masses, int nr, int jrand, Real masses[],
	  Real charges[], Real **j_mu[]){
  
  write_tslice_values_begin(tag);
  
#ifdef MASS_UDLSC
  
  switch(n_masses){
    
  case(1): /* u */
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0], 0., 0., j_mu[0][ir]);
    break;
    
  case(2): /* l - s */
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0],
			  masses[1], charges[1], j_mu[0][ir]);
    break;
    
  case(3):  /* l - s,  c */
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0],
			  masses[1], charges[1], j_mu[0][ir]);
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[2], charges[2], 0., 0., j_mu[2][ir]);
    break;
    
  case(4):  /* u - d,  u - s */
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0],
			  masses[1], charges[1], j_mu[0][ir]);
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0],
			  masses[2], charges[2], j_mu[2][ir]);
    break;
    
  case(5):  /* u - d,  l - s,  s,  c */
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[0], charges[0],
			  masses[1], charges[1], j_mu[0][ir]);
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[2], charges[2],
			  masses[3], charges[3], j_mu[2][ir]);
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[3], charges[3], 0., 0., j_mu[3][ir]);
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[4], charges[4], 0., 0., j_mu[4][ir]);

#if 0
    // DEBUG
    int i;
    FORALLFIELDSITES(i){
      if(lattice[i].t == 0){
	for(int ir = 0; ir < nr; ir++){
	  printf("JMU0[%d] %d %d %d %d %30.20f %30.20f %30.20f %30.20f\n", ir,
		  lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
		  j_mu[0][ir][4*i], j_mu[0][ir][4*i+1],
		  j_mu[0][ir][4*i+2], j_mu[0][ir][4*i+3]);
	}
      }
    }
#endif
    
    break;
    
  }
#else
  for(int j = 0; j < n_masses; j++)
    for(int ir = 0; ir < nr; ir++)
      write_tslice_values(tag, jrand+ir, masses[j], charges[j], 0., 0., j_mu[j][ir]);
#endif

  write_tslice_values_end(tag);
}

#ifdef U1_SEA

/*********************************************************************/
/* Write the high-mode result for this block of nr sources and clear j_mu */
static void
write_jhi_sea(char tag[], int n_masses, int nr, int nU1_A, int jrand,
	      Real masses[], Real charges[], Real **j_sea[]){
  
#if 0
  // DEBUG
  int i;
  FORALLFIELDSITES(i){
    if(lattice[i].t == 0){
      for(int ir = 0; ir < nr; ir++){
	printf("JSEA0[%d] %d %d %d %d %30.20f\n", ir,
	       lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t, j_mu[0][ir][i]);
      }
    }
  }
#endif
  
  write_sea_values_begin(tag);
  
  for(int j = 0; j < n_masses; j++){
    write_sea_values_group(tag, masses[j], charges[j]);
      for(int iU1_A = 0; iU1_A < nU1_A; iU1_A++)
	for(int ir = 0; ir < nr; ir++)
	write_sea_values(tag, jrand+ir, iU1_A, j_sea[j][iU1_A][ir]);
  }
  
  write_sea_values_end(tag);
}

#endif
/*********************************************************************/
/* Print the high-mode contribution to the QED current loop for this
   block of nr sources */
static void
write_jhidotA(char tag[], int n_masses, int nr, int jrand, Real masses[],
	      Real charges[], Real **j_mu[], Real *u1_A){

#ifdef MASS_UDLSC
  switch(n_masses){
    
  case(1): /* u */
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0], 0., 0., j_mu[0][ir], u1_A);
    break;
    
  case(2): /* l - s */
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0],
			masses[1], charges[1], j_mu[0][ir], u1_A);
    break;
    
  case(3):  /* l - s,  c */
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0], 
			masses[1], charges[1], j_mu[0][ir], u1_A);
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[2], charges[2], 0., 0., j_mu[2][ir], u1_A);
    break;
    
  case(4):  /* u - d,  u - s */
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0],
			masses[1], charges[1], j_mu[0][ir], u1_A);
    
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0],
			masses[2], charges[2], j_mu[2][ir], u1_A);
    break;
    
  case(5):  /* u - d,  l - s,  c */
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[0], charges[0],
			masses[1], charges[1], j_mu[0][ir], u1_A);
    
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[2], charges[2],
			masses[3], charges[3], j_mu[2][ir], u1_A);
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[4], charges[4], 0., 0., j_mu[4][ir], u1_A);
    break;
  }
#else
  for(int j = 0; j < n_masses; j++)
    for(int ir = 0; ir < nr; ir++)
      write_jdotA_value(tag, jrand+ir, masses[j], charges[j], 0., 0., j_mu[j][ir], u1_A);
#endif
}

/*********************************************************************/
/* Reset j_mu */
static void
clear_jhi(int n_masses, int nr, Real **j_mu[]){

  for(int j = 0; j < n_masses; j++){
    for(int ir = 0; ir < nr; ir++){
      clear_r_array_field(j_mu[j][ir], NMU);
    } /* ir */
  } /* j */
}

/************************************************************************/
/* Entry point for multiple masses with deflation and iterated
   single-mass inverter.  This variant does two solves from the same
   source -- sloppy and precise -- and calculates the average of the
   difference between the resulting current densities.  Designed for
   use with deflation or eigcg.  The optional set of accurate low-mode
   eigenpairs are in the globals eigVal and eigVec. */
/************************************************************************/

void 
f_meas_valence_current_diff( int n_masses, int nrand, int thinning,
			     quark_invert_control *qic_precise,
			     quark_invert_control *qic_sloppy,
			     Real masses[], Real charges[],
			     int eigstartflags[], char *eigstartfiles[],
			     ks_eigen_param *eigen_param,
			     imp_ferm_links_t *fn_mass[], 
			     Real *u1_A, char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_valence_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Debug */
#ifdef FIXEDSOURCE
  if(nrand != 3){
    node0_printf("%s: FIXEDSOURCE requires npbp_reps = 3\n", myname);
    terminate(1);
  }
#endif

  int nr = NRMAX;  /* Number of random sources to include in multi-rhs */
  if(nrand < nr)
    nr = nrand;

  /* Create fields for high-mode current densities, one for each
     block of nr random sources and mass */
  Real **j_mu[n_masses];
  create_jhi(n_masses, nr, j_mu);

  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  /* Create solver temporaries */
  su3_vector **gr_even = create_gr(nsrc);
  su3_vector **gr_odd = create_gr(nsrc);

  /* Compute high-mode current density stochastically in blocks of size nr */

  /* Create sources in gr0. gr_even and gr_odd are filled in later */
  su3_vector **gr0 = create_source_array(nr);

  for(int jrand = 0; jrand < nrand; jrand += nr){

    double meastime = -dclock();
    
    /* Create sources in gr0. gr_even and gr_odd are filled in later */
    fill_sources(gr0, nr);

    /* Do block solve for the difference in high-mode densities for nr sources */
    compute_jhi_diff(n_masses, j_mu, masses, eigstartflags, eigstartfiles, eigen_param,
		     fn_mass, qic_precise, qic_sloppy,
		     nsrc, nr, thinning, gr_even, gr_odd, gr0);
   
#ifdef CURRENT_DISC
    /* Write the high-mode result for this block of nr sources and clear j_mu */
    write_jhi("DIFF", n_masses, nr, jrand, masses, charges, j_mu);
#endif
    /* Reset j_mu */
    clear_jhi(n_masses, nr, j_mu);

    meastime += dclock();
    node0_printf("Time to do one measurement step %e\n", meastime); fflush(stdout);
    
  } /* jrand */

  /* Clean up */
  destroy_jhi(n_masses, nr, j_mu);
  destroy_source_array(gr0, nr);
  destroy_gr(nsrc, gr_odd);
  destroy_gr(nsrc, gr_even);

  fflush(stdout);
  
} /* f_meas_valence_current_diff */

/*********************************************************************/
/* Entry point for multiple masses with iterated single-mass inverter.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */
/*********************************************************************/

void 
f_meas_valence_current( int n_masses, int nrand, int thinning,
			quark_invert_control *qic, 
			Real masses[], Real charges[],
			int eigstartflags[], char *eigstartfiles[],
			ks_eigen_param *eigen_param,
			imp_ferm_links_t *fn_mass[], 
			Real *u1_A, char filenames[][MAXFILENAME]){

  char myname[] = "f_meas_valence_current";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Debug */
#ifdef FIXEDSOURCE
  if(nrand != 3){
    node0_printf("%s: FIXEDSOURCE requires npbp_reps = 3\n", myname);
    terminate(1);
  }
#endif

  int i;
  int Nvecs = param.eigen_param.Nvecs;
  int nr = NRMAX;  /* Number of random sources to block */
  if(nrand < nr)
    nr = nrand;

  /* Check orthonormality of a few eigenvectors (debugging) */
#if 0
  check_eigen(Nvecs);
#endif

  /* Create fields for low-mode current densities, one for each mass */
  Real *jlow_mu[n_masses];
  create_jlow(n_masses, Nvecs, jlow_mu);

  /* Create fields for high-mode current densities, one for each
     block of nr random sources and mass */
  Real **j_mu[n_masses];
  create_jhi(n_masses, nr, j_mu);

  /* Compute exact low-mode current density if we have eigenvectors to do it */
  compute_jlow(Nvecs, n_masses, jlow_mu, masses, eigstartflags, eigstartfiles, eigen_param, fn_mass);

  /* Print the exact low mode contribution to the current density */
#ifdef CURRENT_DISC
  if(Nvecs > 0) write_jlow(n_masses, masses, charges, jlow_mu);
#endif

  clear_jlow(n_masses, jlow_mu);

  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  /* Create solver temporaries */
  su3_vector **gr_even = create_gr(nsrc);
  su3_vector **gr_odd = create_gr(nsrc);

  /* Compute high-mode current density stochastically in blocks of size nr */

  /* Create sources in gr0. gr_even and gr_odd are filled in later */
  su3_vector **gr0 = create_source_array(nr);

  for(int jrand = 0; jrand < nrand; jrand += nr){

    double meastime = -dclock();

    /* Create sources in gr0. gr_even and gr_odd are filled in later */
    fill_sources(gr0, nr);

    /* Do block solve for high-mode densities for nr sources */
    compute_jhi(n_masses, j_mu, masses, eigstartflags, eigstartfiles, eigen_param,
		fn_mass, qic, nsrc, nr, thinning, gr_even, gr_odd, gr0);
    
#ifdef CURRENT_DISC
    /* Write the high-mode result for this block of nr sources and clear j_mu */
    write_jhi("HI", n_masses, nr, jrand, masses, charges, j_mu);
#endif

    /* Reset j_mu */
    clear_jhi(n_masses, nr, j_mu);

    meastime += dclock();
    node0_printf("Time to do one measurement step %e\n", meastime); fflush(stdout);
  } /* jrand */
  
  /* Clean up */
  destroy_jhi(n_masses, nr, j_mu);
  destroy_jlow(Nvecs, n_masses, jlow_mu);
  destroy_source_array(gr0, nr);
  destroy_gr(nsrc, gr_odd);
  destroy_gr(nsrc, gr_even);

  fflush(stdout);
  
} /* f_meas_valence_current */

#ifdef U1_SEA

/************************************************************************/
/* Entry point for multiple masses with deflation and iterated
   single-mass inverter.  This variant does two solves from the same
   source -- sloppy and precise -- and calculates the average of the
   difference between the resulting current densities.  Designed for
   use with deflation or eigcg.  The optional set of accurate low-mode
   eigenpairs are in the globals eigVal and eigVec. */
/************************************************************************/

void 
f_meas_sea_current_diff( int n_masses, int nrand, int thinning,
			 quark_invert_control *qic_precise,
			 quark_invert_control *qic_sloppy,
			 Real masses[], Real charges[],
			 int nU1_A, Real *u1_A_array[],
			 imp_ferm_links_t *fn_mass,
			 char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_sea_current_diff";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Debug */
#ifdef FIXEDSOURCE
  if(nrand != 3){
    node0_printf("%s: FIXEDSOURCE requires npbp_reps = 3\n", myname);
    terminate(1);
  }
#endif

  int nr = NRMAX;  /* Max number of random sources to include in multi-rhs */
  if(nrand < nr)
    nr = nrand;

  int nU1_A_array = nU1_A;
  if(nU1_A == 0){
    nU1_A_array = 1;
  }

  /* Allocate storage for high-mode sea-quark loops, one for each
     block of nr random sources, mass, and U(1) field */
  Real **jE_sea[n_masses], **jO_sea[n_masses];
  create_jhi_sea(n_masses, nU1_A_array, nr, jE_sea);
  create_jhi_sea(n_masses, nU1_A_array, nr, jO_sea);

  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  /* Allocate source arrays */
  su3_vector **gr_even = create_gr(nsrc);
  su3_vector **gr_odd = create_gr(nsrc);
  su3_vector **gr0 = create_source_array(nr);
  
  /* Compute high-mode sea-loop term stochastically in blocks of size nr */
  
  for(int jrand = 0; jrand < nrand; jrand += nr){
    
    double meastime = -dclock();
    
    /* Create sources in gr_even and gr_odd for block solves */
    fill_sources(gr0, nr);
    for(int jr = 0; jr < nr; jr++){
      collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
      collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
    }
    
    /* Do block solve for the difference in high-mode densities for nr sources */
    compute_jhi_sea_diff(n_masses, jE_sea, jO_sea, masses, charges, fn_mass,
			 qic_precise, qic_sloppy, nsrc, nr, nU1_A_array, u1_A_array,
			 gr_even, gr_odd);
      

#ifdef CURRENT_DISC
    /* Write the high-mode result for this block of nr sources and nU1_A U(1) fields  */
    if(nU1_A == 0){
      write_jhi_sea("DIFF0b", n_masses, nr, 1, jrand, masses, charges, jE_sea);
      write_jhi_sea("DIFF0ad", n_masses, nr, 1, jrand, masses, charges, jO_sea);
    } else {
      write_jhi_sea("DIFFb", n_masses, nr, nU1_A_array, jrand, masses, charges, jE_sea);
      write_jhi_sea("DIFFad", n_masses, nr, nU1_A_array, jrand, masses, charges, jO_sea);
    }
#endif
      
    /* Reset j_sea */
    clear_jhi_sea(n_masses, nU1_A_array, nr, jE_sea);
    clear_jhi_sea(n_masses, nU1_A_array, nr, jO_sea);
    
    meastime += dclock();
    node0_printf("Time to do one measurement step %e\n", meastime); fflush(stdout);
    
  } /* jrand */
  
  /* Clean up */
  destroy_jhi_sea(n_masses, nU1_A_array, jE_sea);
  destroy_jhi_sea(n_masses, nU1_A_array, jO_sea);
  destroy_source_array(gr0, nr);
  destroy_gr(nsrc, gr_odd);
  destroy_gr(nsrc, gr_even);
  
  fflush(stdout);
  
} /* f_meas_sea_current_diff */

/*********************************************************************/
/* Entry point for computing the T1 term. Takes ultiple masses
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */
/*********************************************************************/

void 
f_meas_sea_current( int n_masses, int nrand, int thinning,
		    quark_invert_control *qic, 
		    Real masses[], Real charges[],
		    int nU1_A, Real *u1_A_array[],
		    imp_ferm_links_t *fn_mass, 
		    char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_sea_current";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Debug */
#ifdef FIXEDSOURCE
  if(nrand != 3){
    node0_printf("%s: FIXEDSOURCE requires npbp_reps = 3\n", myname);
    terminate(1);
  }
#endif

  int Nvecs = param.eigen_param.Nvecs;

  int nr = NRMAX;  /* Max number of random sources to include in multi-rhs */
  if(nrand < nr)
    nr = nrand;

  /* Check orthonormality of a few eigenvectors (debugging) */
#if 0
  check_eigen(Nvecs);
#endif

  int nU1_A_array = nU1_A;
  if(nU1_A == 0){
    nU1_A_array = 1;
  }

  /* Allocate storage for high-mode T1 sea-quark loops, one for each
     block of nr random sources, mass, and U(1) field */
  Real **jE_sea[n_masses], **jO_sea[n_masses];
  create_jhi_sea(n_masses, nU1_A_array, nr, jE_sea);
  create_jhi_sea(n_masses, nU1_A_array, nr, jO_sea);

  /* Allocate storage for low-mode sea-quark loops, one for each mass and U(1) field */
  Real *jElow_sea[n_masses], *jOlow_sea[n_masses];
  if(Nvecs > 0){
    create_jlow_sea(n_masses, nU1_A_array, jElow_sea);
    create_jlow_sea(n_masses, nU1_A_array, jOlow_sea);
  }

  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  /* Allocate source arrays */
  su3_vector **gr_even = create_gr(nsrc);
  su3_vector **gr_odd = create_gr(nsrc);
  su3_vector **gr0 = create_source_array(nr);
  
  /* Compute exact low-mode sea-quark loop if we have eigenvectors to do it */
  if(Nvecs > 0){
    
    compute_jlow_sea(n_masses, jElow_sea, jOlow_sea, masses, charges,
		     nU1_A_array, u1_A_array, fn_mass);
      
    /* Print the exact low mode contribution to the current density */
#ifdef CURRENT_DISC
    if(nU1_A == 0){
      write_jlow_sea("LOW0b", n_masses, 1, masses, charges, jElow_sea);
      write_jlow_sea("LOW0ad", n_masses, 1, masses, charges, jOlow_sea);
    } else {
      for(int jU1_A = 0; jU1_A < nU1_A_array; jU1_A++){
	write_jlow_sea("LOWb", n_masses, nU1_A_array, masses, charges, jElow_sea);
	write_jlow_sea("LOWd", n_masses, nU1_A_array, masses, charges, jOlow_sea);
      }
    }
#endif
    clear_jlow_sea(n_masses, nU1_A_array, jElow_sea);
    clear_jlow_sea(n_masses, nU1_A_array, jOlow_sea);
  }

  /* Compute high-mode T1 sea-loop term stochastically in blocks of size nr */
  
  for(int jrand = 0; jrand < nrand; jrand += nr){

    double meastime = -dclock();
    
    /* Create sources in gr_even and gr_odd for block solves */
    fill_sources(gr0, nr);
    for(int jr = 0; jr < nr; jr++){
      collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 1);
      collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 1);
    }
    
    /* Compute the T1 high-mode sea-loops for nr sources */
    compute_jhi_sea(n_masses, jE_sea, jO_sea, masses, charges, fn_mass,
		    qic, nsrc, nr, nU1_A_array, u1_A_array,
		    gr_even, gr_odd);
      
#ifdef CURRENT_DISC
    /* Write the high-mode T1 result for this block of nr sources and nU1_A U(1) fields  */
    if(nU1_A == 0){
      write_jhi_sea("HI0b", n_masses, nr, 1, jrand, masses, charges, jE_sea);
      write_jhi_sea("HI0ad", n_masses, nr, 1, jrand, masses, charges, jO_sea);
    } else {
      write_jhi_sea("HIb", n_masses, nr, nU1_A_array, jrand, masses, charges, jE_sea);
      write_jhi_sea("HIad", n_masses, nr, nU1_A_array, jrand, masses, charges, jO_sea);
    }
#endif

    /* Reset j_sea */
    clear_jhi_sea(n_masses, nU1_A_array, nr, jE_sea);
    clear_jhi_sea(n_masses, nU1_A_array, nr, jO_sea);
      
    meastime += dclock();
    node0_printf("Time to do one T1 measurement step %e\n", meastime); fflush(stdout);

  } /* jrand */

  /* Clean up */
  destroy_jhi_sea(n_masses, nU1_A_array, jE_sea);
  destroy_jhi_sea(n_masses, nU1_A_array, jO_sea);
  destroy_source_array(gr0, nr);
  destroy_gr(nsrc, gr_odd);
  destroy_gr(nsrc, gr_even);
  if(Nvecs > 0){
    destroy_jlow_sea(n_masses, jElow_sea);
    destroy_jlow_sea(n_masses, jOlow_sea);
  }
  
  fflush(stdout);
  
} /* f_meas_sea_current */

/*********************************************************************/
/* Entry point for the T2 term.  Takes multiple masses.
   Designed for use with deflation or eigcg.
   Does deflation, so requires a set of accurate low-mode eigenpairs */
/*********************************************************************/

void 
f_meas_sea_current_T2( int n_masses, int nrand, int thinning,
		       quark_invert_control *qic, 
		       Real masses[], Real charges[],
		       int nU1_A, Real *u1_A_array[],
		       imp_ferm_links_t *fn_mass, 
		       char filenames[][MAXFILENAME]){
  
  char myname[] = "f_meas_sea_current_T2";
  //  node0_printf("Entered %s\n", myname); fflush(stdout);

  /* Debug */
#ifdef FIXEDSOURCE
  if(nrand != 3){
    node0_printf("%s: FIXEDSOURCE requires npbp_reps = 3\n", myname);
    terminate(1);
  }
#endif

  int nr = NRMAX;  /* Max number of random sources to include in multi-rhs */
  if(nrand < nr)
    nr = nrand;

  int nU1_A_array = nU1_A;
  if(nU1_A == 0){
    nU1_A_array = 1;
  }

  /* Allocate storage for high-mode T2 sea-quark loops, one for each
     block of nr random sources, mass, and U(1) field */
  Real **j2E_sea[n_masses], **j2O_sea[n_masses];
  create_jhi_sea(n_masses, nU1_A_array, nr, j2E_sea);
  create_jhi_sea(n_masses, nU1_A_array, nr, j2O_sea);

  /* Block solver parameters */
  int d = thinning;
  int evol = d*d*d*d/2;
  int nsrc = evol*nr;
  
  /* Allocate source arrays */
  su3_vector **gr_even = create_gr(nsrc);
  su3_vector **gr_odd = create_gr(nsrc);
  su3_vector **gr0 = create_source_array(nr);
  
  /* Compute all-mode T2 sea-loop term stochastically in blocks of size nr */

  for(int jrand = 0; jrand < nrand; jrand += nr){

    double meastime = -dclock();

    /* Create sources in gr_even and gr_odd for block solves */
    fill_sources(gr0, nr);
    for(int jr = 0; jr < nr; jr++){
      collect_evenodd_sources(gr_even + jr*evol, evol, EVEN, thinning, gr0[jr], 0);
      collect_evenodd_sources(gr_odd  + jr*evol, evol, ODD,  thinning, gr0[jr], 0);
    }
    
    /* Compute the T2 sea-loops for nr sources (no separation of low and high modes) */
    compute_jhi_sea_T2(n_masses, j2E_sea, j2O_sea, masses, charges, fn_mass,
		       qic, nsrc, nr, nU1_A_array, u1_A_array,
		       gr_even, gr_odd);
    
#ifdef CURRENT_DISC
    /* Write the high-mode T2 result for this block of nr sources and nU1_A U(1) fields  */
    write_jhi_sea("T2c", n_masses, nr, nU1_A_array, jrand, masses, charges, j2E_sea);
    write_jhi_sea("T2odd", n_masses, nr, nU1_A_array, jrand, masses, charges, j2O_sea);
#endif
    
    /* Reset j2_sea */
    clear_jhi_sea(n_masses, nU1_A_array, nr, j2E_sea);
    clear_jhi_sea(n_masses, nU1_A_array, nr, j2O_sea);
    
    meastime += dclock();
    node0_printf("Time to do one T2 measurement step %e\n", meastime); fflush(stdout);
    
  } /* jrand */

  /* Clean up */
  destroy_jhi_sea(n_masses, nU1_A_array, j2E_sea);
  destroy_jhi_sea(n_masses, nU1_A_array, j2O_sea);
  destroy_source_array(gr0, nr);
  destroy_gr(nsrc, gr_odd);
  destroy_gr(nsrc, gr_even);
  fflush(stdout);
  
} /* f_meas_sea_current_T2 */

#endif
