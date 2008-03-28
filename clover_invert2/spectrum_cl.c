/***************** spectrum_cl.c *****************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

/* These procedures compute the zero momentum meson correlators for
   point, rotated, and "smeared" sinks and the zero momentum baryon
   correlators for a point sink.

   Here "point" means both quark and antiquark (or three quarks) are
   tied together at the same sink point.

   "Rotated" means a Fermilab rotation has been applied to the
   quark and antiquark propagators before tying them together.
   
   "Smeared" means the the quark field is convoluted with a sink wave
   function before tying it to the antiquark.  Thus the smearing
   function is a true wave function in the relative coordinate.

*/

#include "cl_inv_includes.h"
#include <string.h>
#include <time.h>

#define MAX_BARYON_PROP 4
static char *bar_kind[MAX_BARYON_PROP] = 
  {"PROTON","PROTON0","DELTA","DELTA0"};
static complex *bar_prop[MAX_BARYON_PROP];
#define MAX_BARYON_PROP_OFFDIAG 6
static char *bar_kind_offdiag[MAX_BARYON_PROP_OFFDIAG] = 
  {"PROTON","PROTON0","DELTA","DELTA0", "LAMBDA","LAMBDA0"};
static complex *bar_prop100[MAX_BARYON_PROP];
static complex *bar_prop011[MAX_BARYON_PROP];

static complex **pmes_prop = NULL;
static complex **rmes_prop = NULL;
static complex **smes_prop = NULL;

/* Table of unique source/sink gamma pairs */
static int *src_gamma, *snk_gamma, num_src_snk_gammas;

/* Table of unique momentum/parity combinations */
static int **momentum, num_mom_parities;
char **parity;

/* Hash table of correlators (indexed as they appear in param file)*/
static int *src_snk_gamma_index, *mom_parity_index;
static int **corr_table, *num_corr_mom_parity;

/* Count of correlator contributions */
static int *num_corr_occur;

/*--------------------------------------------------------------------*/
/* indexing is prop[meson_type][momentum][time] */

static complex ** 
create_mes_prop(int ncor, int ntime){
  complex **prop;
  int m, t;

  prop = (complex **)malloc(ncor*sizeof(complex *));
  if(prop == NULL)return prop;

  for(m = 0; m < ncor; m++){
    prop[m] = (complex *)malloc(ntime*sizeof(complex));
    if(prop[m] == NULL)return NULL;
    
    for(t = 0; t < nt; t++){
      prop[m][t].real = 0.0; prop[m][t].imag = 0.0;
    }
  }
  
  return prop;
}

/*--------------------------------------------------------------------*/

static void 
destroy_mes_prop(complex ***prop, int ncor){
  int m;

  if(*prop == NULL)return;

  for(m = 0; m < ncor; m++){
    if((*prop)[m] != NULL){
      free((*prop)[m]);
    }
  }
  
  free(*prop);
  *prop = NULL;
}

/*--------------------------------------------------------------------*/
/* Construct hash table for the momentum/parity combination.  Look for
   the momentum/parity in the hash table.  If found, return its index.
   Otherwise, add it to the table and assign it the new index */

static int 
hash_meson_mom_parity(int **mom, char **par, int *meson_mom_test,
		      char *par_test, int *nmp){
  int i;
  char myname[] = "hash_meson_mom_parity";

  for(i = 0; i < *nmp; i++){
    if(mom[i][0] == meson_mom_test[0] &&
       mom[i][1] == meson_mom_test[1] &&
       mom[i][2] == meson_mom_test[2] &&
       par[i][0] == par_test[0] &&
       par[i][1] == par_test[1] &&
       par[i][2] == par_test[2])
      return i;
  }

  /* Add to table */

  mom[*nmp] = (int *)malloc(3*sizeof(int));
  if(mom[*nmp] == NULL){
    printf("%s(%d): No room for mom\n",myname, this_node);
    terminate(1);
  }

  par[*nmp] = (char *)malloc(3*sizeof(int));
  if(par[*nmp] == NULL){
    printf("%s(%d): No room for par\n",myname, this_node);
    terminate(1);
  }

  mom[*nmp][0] = meson_mom_test[0];
  mom[*nmp][1] = meson_mom_test[1];
  mom[*nmp][2] = meson_mom_test[2];
  par[*nmp][0] = par_test[0];
  par[*nmp][1] = par_test[1];
  par[*nmp][2] = par_test[2];

  (*nmp)++;
  return *nmp-1;
}

/*--------------------------------------------------------------------*/
/* Create some supporting tables for the correlator calculation */

/* The momentum/parity table list contains all the unique
   momentum/parity combinations. The mom_parity_index maps the
   correlator to its momentum and parity in the momentum/parity table */

static void 
create_mom_parity(int ***mom, char ***par, int **mp_index, int *nmp, int pair){
  int c;
  int num_corr = param.num_corr[pair];
  char myname[] = "create_mom_parity";

  *mom = (int **)malloc(num_corr*sizeof(int *));
  if(*mom == NULL){
    printf("%s(%d): no room for mom table\n",myname,this_node);
    terminate(1);
  }
  
  *par = (char **)malloc(num_corr*sizeof(int *));
  if(*par == NULL){
    printf("%s(%d): no room for parity table\n",myname,this_node);
    terminate(1);
  }

  *mp_index = (int *)malloc(num_corr*sizeof(int));
  if(*mp_index == NULL){
    printf("%s(%d): no room for mp_index\n",myname,this_node);
    terminate(1);
  }

  /* Run through all correlators, tabulating and hashing the unique
     momentum and parity */
  for(c = 0; c < num_corr; c++){
    (*mp_index)[c] = 
      hash_meson_mom_parity(*mom, *par,
			    param.corr_mom[pair][c], 
			    param.corr_parity[pair][c],
			    nmp);
  }
}

/*--------------------------------------------------------------------*/

static void 
destroy_mom_parity(int ***mom, char ***par, int **mp_index, int nmp, int pair){
  
  int p;

  if(*mom != NULL){
    for(p = 0; p < nmp; p++)
      free((*mom)[p]);
    free(*mom); *mom = NULL;
  }

  if(*par != NULL){
    for(p = 0; p < nmp; p++)
      free((*par)[p]);
    free(*par);  *par = NULL;
  }

  if(*mp_index != NULL){
    free(*mp_index);  *mp_index = NULL;
  }
}

/*--------------------------------------------------------------------*/
/* Construct hash table for the momentum/parity combination.  Look for
   the momentum/parity in the hash table.  If found, return its index.
   Otherwise, add it to the table and assign it the new index */

static int 
hash_src_snk_gammas(int *src_g, int *snk_g, int src_g_test,
		    int snk_g_test, int *ng){
  int i;

  for(i = 0; i < *ng; i++){
    if(src_g[i] == src_g_test &&
       snk_g[i] == snk_g_test)
      return i;
  }

  /* Add to table */

  src_g[*ng] = src_g_test;
  snk_g[*ng] = snk_g_test;

  (*ng)++;
  return *ng-1;
}

/*--------------------------------------------------------------------*/
/* Construct hash table of unique source/sink gamma matrices */

static void 
create_src_snk_gamma(int **src_g, int **snk_g, int **g_index, int *ng, int pair)
{
  char myname[] = "create_src_snk_gamma";
  int c;
  int num_corr = param.num_corr[pair];

  *src_g = (int *)malloc(num_corr*sizeof(int));
  if(*src_g == NULL){
    printf("%s(%d): no room for gamma table\n", myname, this_node);
    terminate(1);
  }
  
  *snk_g = (int *)malloc(num_corr*sizeof(int *));
  if(*snk_g == NULL){
    printf("%s(%d): no room for gamma table\n", myname, this_node);
    terminate(1);
  }

  *g_index = (int *)malloc(num_corr*sizeof(int));
  if(*g_index == NULL){
    printf("%s(%d): no room for gamma index table\n", myname, this_node);
    terminate(1);
  }

  /* Run through all correlators, tabulating and hashing the unique
     momentum and parity */
  for(c = 0; c < num_corr; c++){
    (*g_index)[c] = 
      hash_src_snk_gammas(*src_g, *snk_g, param.gam_src[pair][c], 
			  param.gam_snk[pair][c], ng);
  }
}


/*--------------------------------------------------------------------*/

static void 
destroy_src_snk_gamma(int **src_g, int **snk_g, int **g_index)
{
  if(*src_g != NULL){
    free(*src_g); *src_g = NULL;
  }

  if(*snk_g != NULL){
    free(*snk_g);  *snk_g = NULL;
  }

  if(*g_index != NULL){
    free(*g_index);  *g_index = NULL;
  }
}

/*--------------------------------------------------------------------*/

static int *
create_num_corr_mom_parity(int *g_index, int ng, int pair)
{
  char myname[] = "create_num_corr_mom_parity";
  int num_corr = param.num_corr[pair];
  int c, g;
  int *n_cmp;

  n_cmp = (int *)malloc(ng*sizeof(int));
  if(n_cmp == NULL){
    printf("%s(%d): no room for num corr table\n",myname,this_node);
    terminate(1);
  }

  /* Count momentum/parity combinations for each unique source/sink gamma */

  for(g = 0; g < ng; g++)
    n_cmp[g] = 0;

  for(c = 0; c < num_corr; c++){
    g = g_index[c];
    n_cmp[g]++;
  }

  return n_cmp;
}


/*--------------------------------------------------------------------*/

static void
destroy_num_corr_mom_parity(int **n_cmp)
{
  if(*n_cmp != NULL){
    free(*n_cmp); *n_cmp = NULL;
  }
}


/*--------------------------------------------------------------------*/
/* Create a table of correlators.  The table consists of a list of
   unique source/sink gamma pairs, and a list of the correlator
   indices for each pair */

static void
fill_corr_table(int *src_g, int *snk_g, int *g_index, int **c_table, 
		int *n_cmp, int ng, int *mp_index, int pair)
{
  char myname[] = "fill_corr_table";
  int c, g, k;
  int num_corr = param.num_corr[pair];
  int *k_corr_mom_parity;

  /* Set up counters, one for each unique source/sink gamma pair */

  k_corr_mom_parity = (int *)malloc(ng*sizeof(int));
  if(k_corr_mom_parity == NULL){
    printf("%s(%d): no room for counter\n",myname,this_node);
    terminate(1);
  }

  for(g = 0; g < ng; g++)
    k_corr_mom_parity[g] = 0;


  /* Allocate rows of the correlator table, one row for each
     source/sink gamma pair */

  for(g = 0; g < ng; g++){
    c_table[g] = (int *)malloc(sizeof(int)*n_cmp[g]);
    if(c_table[g] == NULL){
      printf("%s(%d): no room for corr table\n",myname, this_node);
       terminate(1);
     }
  }

  /* Populate the rows of the correlator table */
  for(c = 0; c < num_corr; c++){
    g = g_index[c];    /* Hash of gamma matrix combination for this corr */
    k = k_corr_mom_parity[g];
    c_table[g][k] = c;  /* Point table to correlator */
    k_corr_mom_parity[g]++;
  }

  free(k_corr_mom_parity);
}

/*--------------------------------------------------------------------*/

/* Create a table of correlators.  The table consists of a list of
   unique source/sink gamma pairs, and a list of the correlator
   indices for each pair */

static void
create_corr_table(int pair)
{

  char myname[] = "create_corr_table";

  num_src_snk_gammas = 0;
  num_mom_parities = 0;

  create_src_snk_gamma(&src_gamma, &snk_gamma, &src_snk_gamma_index,
		       &num_src_snk_gammas, pair);

  create_mom_parity(&momentum, &parity, &mom_parity_index, 
			  &num_mom_parities, pair);

  num_corr_mom_parity = 
    create_num_corr_mom_parity(src_snk_gamma_index, num_src_snk_gammas,
			       pair);

  /* Create the correlator table */
  
  corr_table = (int **)malloc(num_src_snk_gammas*sizeof(int *));
  if(corr_table == NULL){
    printf("%s(%d): no room for corr table\n",myname, this_node);
    terminate(1);
  }

  /* Fill it */

  fill_corr_table(src_gamma, snk_gamma, src_snk_gamma_index, corr_table,
		  num_corr_mom_parity, num_src_snk_gammas, 
		  mom_parity_index, pair);

}

/*--------------------------------------------------------------------*/
static void
destroy_corr_table(int pair)
{
  int g;

  destroy_src_snk_gamma(&src_gamma, &snk_gamma, &src_snk_gamma_index);

  destroy_mom_parity(&momentum, &parity, &mom_parity_index,
		     num_mom_parities, pair);

  destroy_num_corr_mom_parity(&num_corr_mom_parity);

  if(corr_table != NULL){
    for(g = 0; g < num_src_snk_gammas; g++){
      if(corr_table[g] != NULL)
	free(corr_table[g]);
    }
    free(corr_table);
    corr_table = NULL;
  }
}

/*--------------------------------------------------------------------*/
/* The num_corr_occur array counts the number of contributions to
   a given correlator for normalization purposes    */

static int *
create_num_corr_occur(int pair){
  int c,m;
  int ncr = param.num_corr_report[pair];
  int num_corr = param.num_corr[pair];
  int *nco;
  char myname[] = "create_num_corr_occur";

  nco = (int *)malloc(sizeof(int)*ncr);
  if(nco == NULL){
    printf("%s(%d) no room for num_corr_occur\n", myname, this_node);
    terminate(1);
  }

  for(m = 0; m < ncr; m++)
    nco[m] = 0;
  
  for(c = 0; c < num_corr; c++){
    m = param.corr_index[pair][c];
    nco[m]++;
  }

  return nco;
}
/*--------------------------------------------------------------------*/
static void
destroy_num_corr_occur(int **nco){
  if(*nco != NULL){
    free(*nco); *nco = NULL;
  }
}

/*--------------------------------------------------------------------*/
static void 
spectrum_cl_init(int pair){
  
  int num_corr = 
    param.num_corr_report[pair]; /* number of meson correlators to store */
  int t, num_prop;
  int *qkpair = param.qkpair[pair];
  int key[4];
#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];
  
  if(param.do_point_meson_spect[pair])
    pmes_prop = create_mes_prop(num_corr, nt);

  if(param.do_smear_meson_spect[pair])
    smes_prop = create_mes_prop(num_corr, nt);

  if(param.do_rot_meson_spect[pair])
    rmes_prop = create_mes_prop(num_corr, nt);

  if(qkpair[0] == qkpair[1])
    {
      if(param.do_baryon_spect[pair])
	for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++){
	  bar_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++)
	    bar_prop[num_prop][t] = cmplx(0.0,0.0); 
	}
    } 
  else 
    {

      /* Off-diagonal case */

      if(param.do_baryon_spect[pair])
	for(num_prop=0;num_prop<MAX_BARYON_PROP_OFFDIAG;num_prop++){
	  bar_prop100[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  bar_prop011[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++){
	    bar_prop100[num_prop][t] = cmplx(0.0,0.0); 
	    bar_prop011[num_prop][t] = cmplx(0.0,0.0); 
	  }
	}
    }
  
  /* Set up Fourier transform for smearing */
  if(param.do_smear_meson_spect[pair]){
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 0;
    setup_restrict_fourier(key, restrict);
  }

  /* Set up count of contributions to the correlators */
  num_corr_occur = create_num_corr_occur(pair);
  
  /* Set up table of requested correlators */
  create_corr_table(pair);

}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_baryon(wilson_prop_field qp){
  
  double dtime = start_timing();
  w_baryon(qp, qp, qp, bar_prop);
  print_timing(dtime, "diagonal baryons");
}

/*---------------------------------------------------------------------*/
static void spectrum_cl_offdiag_baryon(wilson_prop_field qp0, 
				       wilson_prop_field qp1)
{
  double dtime = start_timing();
  w_baryon_hl( qp1, qp0, qp0, bar_prop100);
  w_baryon_hl( qp0, qp1, qp1, bar_prop011);
  print_timing(dtime, "offdiag baryons");
}

/*--------------------------------------------------------------------*/
static void rotate_prop(spin_wilson_vector *rp, wilson_prop_field qp, 
			int color, Real d1){
  
  int spin;
  int i;
  site *s;
  wilson_vector *psi, *mp, *tmp;
  
  psi = create_wv_field();
  mp  = create_wv_field();
  tmp = create_wv_field();
  
  /* Construct propagator for "rotated" fields,
     psi_rot = Dslash psi, with Dslash the naive operator. */
  for(spin=0;spin<4;spin++){
    copy_wv_from_wp(psi, qp, color, spin);
    
    /* Do Wilson Dslash on the psi field */
    dslash_w_3D_field(psi, mp,  PLUS, EVENANDODD);
    dslash_w_3D_field(psi, tmp, MINUS, EVENANDODD);

    FORALLSITES(i,s){
      /* From subtraction we get 2*Dslash */
      sub_wilson_vector(mp + i, tmp + i, &rp[i].d[spin]);
      /* Apply rotation */
      scalar_mult_add_wvec(psi + i, &rp[i].d[spin], d1/4., &rp[i].d[spin]);
    }
  }

  cleanup_dslash_w_3D_temps();
  destroy_wv_field(psi); 
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}

/*--------------------------------------------------------------------*/
static void sink_smear_prop(wilson_prop_field qp, wilson_quark_source *wqs){
  
  int color;
  int ci,si,sf,cf;
  int i;
  site *s;
  complex *chi_cs, z;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    /* qps points into qp, so qp is changed here */
    restrict_fourier_field((complex *)qps, sizeof(spin_wilson_vector), 
			   FORWARDS);
  }

  print_timing(dtime,"FFT");

  w_sink_field(chi_cs, wqs);
  
  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function */
  for(ci=0;ci<3;ci++){
    FORALLSITES(i,s)
      for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
	    z = qp[ci][i].d[si].d[sf].c[cf];
	    CMUL(z, chi_cs[i], qp[ci][i].d[si].d[sf].c[cf]);
	  }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    /* qps points into qp, so qp is changed here */
    restrict_fourier_field((complex *)qps, sizeof(spin_wilson_vector), 
			   BACKWARDS);
  }
  print_timing(dtime,"FFT");
  free(chi_cs);
}  


/*--------------------------------------------------------------------*/

static void accum_gen_meson(complex **mp, spin_wilson_vector *qp0, 
			    spin_wilson_vector *qp1, int pair){

  meson_cont_mom(mp, qp0, qp1, 
		 num_mom_parities, momentum, parity,
		 num_src_snk_gammas, num_corr_mom_parity, corr_table, 
		 mom_parity_index, param.gam_snk[pair], param.gam_src[pair], 
		 param.corr_phase[pair], param.corr_factor[pair],
		 param.corr_index[pair]);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_gen_meson(wilson_prop_field qp, complex **mp, 
				       int pair){
  
  int color;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  for(color=0;color<3;color++){
    qps = extract_swv_from_wp(qp, color);
    accum_gen_meson(mp, qps, qps, pair);
  }
  print_timing(dtime, "diagonal mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_meson(wilson_prop_field qp, int pair){
  spectrum_cl_diag_gen_meson(qp, pmes_prop, pair);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_rot_meson(wilson_prop_field qp, int pair){
  spin_wilson_vector *rps;
  int color;
  double dtime = start_timing();
  int iqk = param.qkpair[pair][0];
  /* The MILC sign convention for gamma matrix in Dslash is
     opposite FNAL, so we rotate with -d1 */
  Real d1 = -param.d1[iqk];

  /* Diagonal quark-mass mesons */
  rps = create_swv_field();
  for(color = 0; color < 3; color++){
    rotate_prop(rps, qp, color, d1);
    accum_gen_meson(rmes_prop, rps, rps, pair);
  }
  destroy_swv_field(rps);

  print_timing(dtime,"diagonal mesons with rotns");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_gen_meson(wilson_prop_field qp0, 
					  wilson_prop_field qp1, 
					  complex **mp, int pair){

  spin_wilson_vector *qps0, *qps1;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  for(color = 0; color < 3; color++){
    qps0 = extract_swv_from_wp( qp0, color);
    qps1 = extract_swv_from_wp( qp1, color);

    /* Point meson propagator */
    accum_gen_meson(mp, qps0, qps1, pair);
  }
  print_timing(dtime, "offdiag mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_smeared_meson(wilson_prop_field qp,
					   wilson_quark_source *wqs,
					   int pair){
  wilson_prop_field qpcopy = create_wp_field_copy(qp);

  sink_smear_prop(qpcopy, wqs);
  spectrum_cl_offdiag_gen_meson(qp, qpcopy, smes_prop, pair);
  destroy_wp_field(qpcopy);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_meson(wilson_prop_field qp0, 
				      wilson_prop_field qp1,
				      int pair){

  /* Point meson */
  spectrum_cl_offdiag_gen_meson(qp0, qp1, pmes_prop, pair);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_rot_meson(wilson_prop_field qp0, 
					  wilson_prop_field qp1, int pair)
{
  spin_wilson_vector *rqs0, *rqs1;
  int color;
  double dtime = start_timing();
  int iqk0 = param.qkpair[pair][0];
  int iqk1 = param.qkpair[pair][1];
  Real d10 = -param.d1[iqk0];
  Real d11 = -param.d1[iqk1];

  /* Off diagonal quark mass mesons */
  rqs0 = create_swv_field();
  rqs1 = create_swv_field();
  for(color = 0; color < 3; color++){
    rotate_prop(rqs0, qp0, color, d10);
    rotate_prop(rqs1, qp1, color, d11);
    accum_gen_meson(rmes_prop, rqs0, rqs1, pair);
  }
  destroy_swv_field(rqs0);
  destroy_swv_field(rqs1);

  print_timing(dtime, "offdiag mesons with rotns: two kappa pairs");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_smeared_meson(wilson_prop_field qp0, 
					      wilson_prop_field qp1, 
					      wilson_quark_source *wqs,
					      int pair){
  wilson_prop_field qpcopy = create_wp_field_copy(qp1);

  sink_smear_prop(qpcopy, wqs);
  spectrum_cl_offdiag_gen_meson(qp0, qpcopy, smes_prop, pair);
  destroy_wp_field(qpcopy);
}

/*--------------------------------------------------------------------*/
static void 
print_start_meson_prop(int pair, int m, char sink[]){
  if(param.saveflag_c[pair] != FORGET)return;
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("MOMENTUM: %s\n", param.mom_label[pair][m]);
  printf("MASSES: ");
  if(param.qk_type[iq0] == CLOVER_TYPE)
    printf("kappa[%d]=%g ", iq0, param.dcp[iq0].Kappa);
  else /* KS_TYPE */
    printf("mass[%d]=%g ", iq0, param.ksp[iq0].mass);
  if(param.qk_type[iq1] == CLOVER_TYPE)
    printf("kappa[%d]=%g ", iq1, param.dcp[iq1].Kappa);
  else /* KS_TYPE */
    printf("mass[%d]=%g ", iq1, param.ksp[iq1].mass);
  printf("\n");
  printf("SOURCES: %s %s\n",param.src_wqs[iq0].descrp,
	 param.src_wqs[iq1].descrp);
  printf("SINK: %s %s\n",sink, param.meson_label[pair][m]);
}
		       
/*--------------------------------------------------------------------*/

static char *get_utc_datetime(void)
{
  time_t time_stamp;
  struct tm *gmtime_stamp;
  static char time_string[64];

  time(&time_stamp);
  gmtime_stamp = gmtime(&time_stamp);
  strncpy(time_string,asctime(gmtime_stamp),64);
  
  /* Remove trailing end-of-line character */
  if(time_string[strlen(time_string) - 1] == '\n')
    time_string[strlen(time_string) - 1] = '\0';
  return time_string;
}

/*--------------------------------------------------------------------*/
static FILE* open_fnal_meson_file(int pair, char sink[]){
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  FILE *fp;

  /* Create the FNAL file only for rotated meson.  Only node 0
     writes. */
  if(this_node != 0 || param.saveflag_c[pair] == FORGET )
    return NULL;

  fp = fopen(param.savefile_c[pair],"a");
  if(fp == NULL){
    printf("print_open_fnal_meson_file: ERROR. Can't open %s\n",
	   param.savefile_c[pair]);
    return NULL;
  }
  fprintf(fp,"# Job ID:               %s\n",param.job_id);
  fprintf(fp,"# date:                 %s UTC\n",get_utc_datetime());
  fprintf(fp,"# lattice size:         %d x %d x %d x %d\n", nx, ny, nz, nt);
  fprintf(fp,"# spatial volume:       %f\n",((float)nx)*ny*nz);
  fprintf(fp,"# sources: %s %s # sink: %s\n",param.src_wqs[iq0].descrp,
	  param.src_wqs[iq1].descrp,sink);
  fprintf(fp,"# rotated current: d1[%d]=%f  d1[%d]=%f\n",
	  iq0, param.d1[iq0], iq1, param.d1[iq1]);
  fprintf(fp,"# masses: ");
  if(param.qk_type[iq0] == CLOVER_TYPE)
    fprintf(fp,"kappa[%d]=%g ", iq0, param.dcp[iq0].Kappa);
  else /* KS_TYPE */
    fprintf(fp,"mass[%d]=%g ", iq0, param.ksp[iq0].mass);
  if(param.qk_type[iq1] == CLOVER_TYPE)
    fprintf(fp,"kappa[%d]=%g ", iq1, param.dcp[iq1].Kappa);
  else /* KS_TYPE */
    fprintf(fp,"mass[%d]=%g ", iq1, param.ksp[iq1].mass);
  fprintf(fp, "\n");
  return fp;
}
		       
/*--------------------------------------------------------------------*/
static void print_start_fnal_meson_prop(FILE *fp, int pair, int m)
{
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;

  fprintf(fp,"# operator %d # element: %s  momentum: %s\n",
	  m,param.meson_label[pair][m],param.mom_label[pair][m]);
}
		       
/*--------------------------------------------------------------------*/
static void print_start_baryon_prop(int iq1, int iq2, int iq3, char sink[]){
  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("SOURCES: %s %s %s\n",param.src_wqs[iq1].descrp,
	 param.src_wqs[iq2].descrp, param.src_wqs[iq3].descrp );
  printf("KAPPAS: %g %g %g\n",param.dcp[iq1].Kappa,
	 param.dcp[iq2].Kappa, param.dcp[iq3].Kappa);
  printf("SINK: POINT %s\n",sink);
}
/*--------------------------------------------------------------------*/
static void print_meson_prop(int pair, int t, complex c)
{
  if(param.saveflag_c[pair] != FORGET)return;
  node0_printf("%d %e %e\n",t,(double)c.real,(double)c.imag);
}
/*--------------------------------------------------------------------*/
static void print_fnal_meson_prop(FILE *fp, int pair, int t, complex c)
{
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fprintf(fp, "%d\t%e\t%e\n", t, (double)c.real, (double)c.imag);
}
/*--------------------------------------------------------------------*/
static void print_end_prop(int pair){
  if(param.saveflag_c[pair] != FORGET)return;
  node0_printf("ENDPROP\n");
}
/*--------------------------------------------------------------------*/
static void print_end_fnal_meson_prop(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fprintf(fp, "&\n");
}
/*--------------------------------------------------------------------*/
static void close_fnal_meson_file(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  if(fp != NULL)fclose(fp);
}
/*--------------------------------------------------------------------*/
static void spectrum_cl_print_diag(int pair){

  Real space_vol;
  Real norm_fac;
  FILE *corr_fp;
  int t;
  int m,b;
  int num_report = param.num_corr_report[pair];
  complex prop;
  
  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);
  
  /* Point sink */
  if(param.do_point_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "local sink");
    
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m, "POINT");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = pmes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);

    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
  
  /* Rotated sink */
  if(param.do_rot_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "point rotated");

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m, "ROTATED");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = rmes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
  
  /* Smeared sink */
  if(param.do_smear_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "smeared sink");

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      print_start_meson_prop(pair, m, param.snk_wqs[pair].label);
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = smes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol*norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
  
  /* print baryon propagators */
  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop(pair);
    }
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_print_offdiag(int pair){
  
  Real space_vol;
  Real norm_fac;
  int t;
  int m,b;
  int num_report = param.num_corr_report[pair];
  complex prop;
  FILE *corr_fp;
  
  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);

  /* Point sink */
  if(param.do_point_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "local sink");

    /* print meson propagators */
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      print_start_meson_prop(pair, m, "POINT");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = pmes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
    
    
  /* Rotated sink */
  if(param.do_rot_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "rotated sink");

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m, "ROTATED");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = rmes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    }
    close_fnal_meson_file(corr_fp, pair);
  }
    
  /* Smeared sink */
  if(param.do_smear_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair, "smeared sink");

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m, param.snk_wqs[pair].label);
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	prop = smes_prop[m][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol*norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
  
  /* print baryon propagators */
  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][1], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop100[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop(pair);
    }
    
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][1], 
			      param.qkpair[pair][1], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop011[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop(pair);
    }
  }
}

/*--------------------------------------------------------------------*/
void spectrum_cl_cleanup(int pair){
  int num_prop;
  int num_corr = param.num_corr_report[pair];
  int *qkpair = param.qkpair[pair];

  if(param.do_point_meson_spect[pair])
    destroy_mes_prop(&pmes_prop, num_corr);
  if(param.do_smear_meson_spect[pair])
    destroy_mes_prop(&smes_prop, num_corr);
  if(param.do_rot_meson_spect[pair])
    destroy_mes_prop(&rmes_prop, num_corr);

  if(qkpair[0] == qkpair[1]){
    if(param.do_baryon_spect[pair])
      for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++)
	free(bar_prop[num_prop]);
  } else {
    if(param.do_baryon_spect[pair])
      for(num_prop=0;num_prop<MAX_BARYON_PROP_OFFDIAG;num_prop++){
	free(bar_prop011[num_prop]);
	free(bar_prop100[num_prop]);
      }
  }

  destroy_corr_table(pair);

  destroy_num_corr_occur(&num_corr_occur);
}

/*--------------------------------------------------------------------*/
void spectrum_cl(wilson_prop_field qp0, wilson_prop_field qp1, int pair)
{

  int *qkpair = param.qkpair[pair];

  spectrum_cl_init(pair);

  if(qkpair[0] == qkpair[1]){

    if(param.do_point_meson_spect[pair])
      spectrum_cl_diag_meson(qp0, pair);

    if(param.do_rot_meson_spect[pair])
      spectrum_cl_diag_rot_meson(qp0, pair);

    if(param.do_smear_meson_spect[pair])
      spectrum_cl_diag_smeared_meson(qp0, &param.snk_wqs[pair], pair);

    if(param.do_baryon_spect[pair])
      spectrum_cl_diag_baryon(qp0);

    spectrum_cl_print_diag(pair);

  } else {

    if(param.do_point_meson_spect[pair])
      spectrum_cl_offdiag_meson(qp0, qp1, pair);

    if(param.do_rot_meson_spect[pair])
      spectrum_cl_offdiag_rot_meson(qp0, qp1, pair);

    if(param.do_smear_meson_spect[pair])
      spectrum_cl_offdiag_smeared_meson(qp0, qp1, &param.snk_wqs[pair], pair);

    if(param.do_baryon_spect[pair])
      spectrum_cl_offdiag_baryon(qp0, qp1);
    
    spectrum_cl_print_offdiag(pair);
  }    
	    
  spectrum_cl_cleanup(pair);
}
