/***************** spectrum_cl.c *****************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

/* These procedures compute meson and baryon correlators.
   They also compute "open meson" correlators for later
   recombination.

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
static complex *bar_prop100[MAX_BARYON_PROP_OFFDIAG];
static complex *bar_prop011[MAX_BARYON_PROP_OFFDIAG];

static complex **pmes_prop = NULL;
static wilson_propagator **omes_prop = NULL;

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
void clear_wprop(wilson_propagator *p){
  int c,s;
  for(c = 0; c < 3; c++)
    for(s = 0; s < 4; s++)
      clear_wvec(&p->c[c].d[s]);
}

/*--------------------------------------------------------------------*/
/* indexing is prop[meson_type][momentum][time] */

static wilson_propagator ** 
create_open_mes_prop(int ncor, int ntime){
  wilson_propagator **prop;
  int m, t;

  prop = (wilson_propagator **)malloc(ncor*sizeof(wilson_propagator *));
  if(prop == NULL)return prop;

  for(m = 0; m < ncor; m++){
    prop[m] = (wilson_propagator *)malloc(ntime*sizeof(wilson_propagator));
    if(prop[m] == NULL)return NULL;
    
    for(t = 0; t < nt; t++){
      clear_wprop(&prop[m][t]);
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

static void 
destroy_open_mes_prop(wilson_propagator ***prop, int ncor){
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
/* Construct hash table for the source/sink gamma combinations.  Look for
   the source/sink gamma in the hash table.  If found, return its index.
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

  /* Run through all correlators for this pair, tabulating and hashing
     the unique source-sink gammas */
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

  if(param.do_meson_spect[pair])
    pmes_prop = create_mes_prop(num_corr, nt);

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
  
  if(param.do_open_meson[pair])
    omes_prop = create_open_mes_prop(num_corr, nt);

  /* Set up count of contributions to the correlators */
  num_corr_occur = create_num_corr_occur(pair);
  
  /* Set up table of requested correlators */
  create_corr_table(pair);

}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_baryon(wilson_prop_field *qp){
  
  double dtime = start_timing();
  w_baryon(qp, qp, qp, bar_prop);
  print_timing(dtime, "diagonal baryons");
}

/*---------------------------------------------------------------------*/
static void spectrum_cl_offdiag_baryon(wilson_prop_field *qp0, 
				       wilson_prop_field *qp1)
{
  double dtime = start_timing();
  w_baryon_hl( qp1, qp0, qp0, bar_prop100);
  w_baryon_hl( qp0, qp1, qp1, bar_prop011);
  print_timing(dtime, "offdiag baryons");
}

/*--------------------------------------------------------------------*/

static void accum_gen_meson(complex **mp, spin_wilson_vector *qp0, 
			    spin_wilson_vector *qp1, int pair){

  meson_cont_mom(mp, qp0, qp1, 
		 num_mom_parities, momentum, parity,
		 num_src_snk_gammas, num_corr_mom_parity, corr_table, 
		 mom_parity_index, param.gam_snk[pair], param.gam_src[pair], 
		 param.corr_phase[pair], param.corr_factor[pair],
		 param.corr_index[pair], &param.r_offset[pair][0]);
}

/*--------------------------------------------------------------------*/

static void accum_open_meson(wilson_propagator **mp, spin_wilson_vector *qp0, 
			     spin_wilson_vector *qp1, int pair){

  meson_open_mom(mp, qp0, qp1, 
		 num_mom_parities, momentum, parity,
		 num_src_snk_gammas, num_corr_mom_parity, corr_table, 
		 mom_parity_index, param.gam_snk[pair], param.gam_src[pair], 
		 param.corr_phase[pair], param.corr_factor[pair],
		 param.corr_index[pair], param.r_offset[pair]);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_gen_meson(wilson_prop_field *qp, complex **mp, 
				       int pair){
  
  int color;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  for(color=0;color<qp->nc;color++){
    qps = extract_swv_from_wp(qp, color);
    accum_gen_meson(mp, qps, qps, pair);
  }
  print_timing(dtime, "diagonal mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_meson(wilson_prop_field *qp, int pair){
  spectrum_cl_diag_gen_meson(qp, pmes_prop, pair);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_gen_meson(wilson_prop_field *qp0, 
					  wilson_prop_field *qp1, 
					  complex **mp, int pair){

  spin_wilson_vector *qps0, *qps1;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  for(color = 0; color < qp0->nc; color++){
    qps0 = extract_swv_from_wp( qp0, color);
    qps1 = extract_swv_from_wp( qp1, color);

    /* Point meson propagator */
    accum_gen_meson(mp, qps0, qps1, pair);
  }
  print_timing(dtime, "offdiag mesons");
}

/*---------------------------------------------------------------------*/
static void spectrum_cl_offdiag_open_meson(wilson_prop_field *qp0,
					   wilson_prop_field *qp1,
					   int pair)
{
  spin_wilson_vector *qps0, *qps1;
  int color;
  double dtime = start_timing();

  for(color=0;color<qp0->nc;color++){
    qps0 = extract_swv_from_wp( qp0, color);
    qps1 = extract_swv_from_wp( qp1, color);

    accum_open_meson(omes_prop, qps0, qps1, pair);
  }
  print_timing(dtime, "open meson");
}
/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_meson(wilson_prop_field *qp0, 
				      wilson_prop_field *qp1,
				      int pair){

  /* Point meson */
  spectrum_cl_offdiag_gen_meson(qp0, qp1, pmes_prop, pair);
}

/*--------------------------------------------------------------------*/
/* Quark propagators can be derived from a sequence of sink operators
   acting on a parent propagator.  We want the complete family history
   in the metadata.
   The iq in the call is the index of the quark propagator.
   The history (list of quark propagators and their parents) 
   is returned in in reverse order, as is the number of generations.
   The ancestral propagator is the return value */

#define MAX_HISTORY 16

int get_ancestors(int h[], int *n, int iq){
  int i;

  /* Scan backwards in the linked list until we hit the
     ancestral propagator */
  h[0] = iq;
  for(i = 0; i < MAX_HISTORY-1; i++){
    if(param.parent_type[h[i]] == PROP_TYPE){
      break;
    }
    h[i+1] = param.prop_for_qk[h[i]];
  }
  
  *n = i + 1;

  if(*n == MAX_HISTORY){
    printf("get_ancestors(%d):Out of space for history\n",this_node);
  }
  
  return param.prop_for_qk[h[i]];
}

/*--------------------------------------------------------------------*/
/* Chase linked list to find the original propagator ancestor */

int get_eve(int iq){
  int i, ip;

  for(i = 0; i < MAX_HISTORY-1; i++){
    ip = param.prop_for_qk[iq];
    if(param.parent_type[iq] == PROP_TYPE) break;
    iq = ip;
  }

  return ip;
}
  
/*--------------------------------------------------------------------*/
/* Print to stdout on node 0 in lieu of a separate correlator file */
static void 
print_start_meson_prop(int pair, int m, char sink[]){
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  int ip0 = get_eve(iq0);
  int ip1 = get_eve(iq1);

  if(this_node != 0)return;
  /* No stdout if there is a separate correlator file */
  if(param.saveflag_c[pair] != FORGET)return;
  printf("STARTPROP\n");
  printf("MOMENTUM: %s\n", param.mom_label[pair][m]);
  printf("MASSES: ");
  if(param.prop_type[ip0] == CLOVER_TYPE)
    printf("kappa[%d]=%g ", ip0, param.dcp[ip0].Kappa);
  else /* KS_TYPE */
    printf("mass[%d]=%g ", ip0, param.ksp[ip0].mass);
  if(param.prop_type[ip1] == CLOVER_TYPE)
    printf("kappa[%d]=%g ", ip1, param.dcp[ip1].Kappa);
  else /* KS_TYPE */
    printf("mass[%d]=%g ", ip1, param.ksp[ip1].mass);
  printf("\n");
  printf("SOURCES: %s %s\n",param.src_qs[ip0].descrp,
	 param.src_qs[ip1].descrp);
  printf("SINKS: %s %s %s\n", param.snk_qs_op[iq0].descrp,
         param.snk_qs_op[iq1].descrp,
         param.meson_label[pair][m]);
}
		       
/*--------------------------------------------------------------------*/
static FILE* open_fnal_meson_file(int pair){
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY];
  int nh0, nh1;
  int k;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  FILE *fp;

  /* Create the FNAL file only if requested.  Only node 0 writes. */
  if(this_node != 0 || param.saveflag_c[pair] == FORGET )
    return NULL;

  fp = fopen(param.savefile_c[pair],"a");
  if(fp == NULL){
    printf("open_fnal_meson_file: ERROR. Can't open %s\n",
	   param.savefile_c[pair]);
    return NULL;
  }
  fprintf(fp,"---\n");
  fprintf(fp,"JobID:                        %s\n",param.job_id);
  fprintf(fp,"date:                         \"%s UTC\"\n",utc_date_time);
  fprintf(fp,"lattice_size:                 %d,%d,%d,%d\n", nx, ny, nz, nt);
  //  fprintf(fp,"spatial volume:         %g\n",((float)nx)*ny*nz);

  if(param.prop_type[ip0] == CLOVER_TYPE)
    fprintf(fp,"antiquark_type:               clover\n");
  else
    fprintf(fp,"antiquark_type:               naive\n");

  print_source_info(fp, "antiquark_source", &param.src_qs[ip0]);
  fprintf(fp,"antiquark_source_label:       %s\n",param.src_qs[ip0].label);
  print_field_op_info(fp, "antiquark_source", param.src_qs[ip0].op);

//  fprintf(fp,"antiquark_sink_op:       %s",param.snk_qs_op[ih0[nh0-1]].descrp);
//  for(i = nh0-2; i >=0; i--)
//    fprintf(fp,"/%s",param.snk_qs_op[ih0[i]].descrp);
//  fprintf(fp,"\n");
  fprintf(fp,"antiquark_sink_label:         %s\n",param.snk_qs_op[iq0].label);
  {
    quark_source_sink_op **op_list = (quark_source_sink_op **)
      malloc(sizeof(quark_source_sink_op *)*nh0);
    for(k = 0; k < nh0; k++)
      op_list[k] = &param.snk_qs_op[ih0[nh0-1-k]];
    print_field_op_info_list(fp, "antiquark_sink", op_list, nh0);
    free(op_list);
  }
      
  if(param.prop_type[ip0] == CLOVER_TYPE)
    fprintf(fp,"antiquark_kappa:              \"%s\"\n",param.kappa_label[ip0]);
  else /* KS_TYPE */
    fprintf(fp,"antiquark_mass:               \"%s\"\n",param.mass_label[ip0]);

  if(param.prop_type[ip1] == CLOVER_TYPE)
    fprintf(fp,"quark_type:                   clover\n");
  else
    fprintf(fp,"quark_type:                   naive\n");
    
  print_source_info(fp, "quark_source", &param.src_qs[ip1]);
  fprintf(fp,"quark_source_label:           %s\n",param.src_qs[ip1].label);
  print_field_op_info(fp, "quark_source", param.src_qs[ip1].op);

//  fprintf(fp,"quark_sink_op:           %s",param.snk_qs_op[ih1[nh1-1]].descrp);
//  for(i = nh1-2; i >=0; i--)
//    fprintf(fp,"/%s",param.snk_qs_op[ih1[i]].descrp);
//  fprintf(fp,"\n");
  fprintf(fp,"quark_sink_label:             %s\n",param.snk_qs_op[iq1].label);
  {
    quark_source_sink_op **op_list = (quark_source_sink_op **)
      malloc(sizeof(quark_source_sink_op *)*nh1);
    for(k = 0; k < nh1; k++)
      op_list[k] = &param.snk_qs_op[ih1[nh1-1-k]];
    print_field_op_info_list(fp, "quark_sink", op_list, nh1);
    free(op_list);
  }
  if(param.prop_type[ip1] == CLOVER_TYPE)
    fprintf(fp,"quark_kappa:                  \"%s\"\n",param.kappa_label[ip1]);
  else /* KS_TYPE */
    fprintf(fp,"quark_mass:                   \"%s\"\n",param.mass_label[ip1]);

  fprintf(fp,"...\n");
  return fp;
}
		       
/*--------------------------------------------------------------------*/
/* Open and write a FermiQCD header for the open meson correlator file */

typedef dcomplex mdp_complex;
static FILE* open_open_meson_file(int pair){

  /* Data members of the FermiQCD header class for the open meson correlator */
  struct {
    char  file_id[60];
    char  program_version[60];
    char  creation_date[60];
    u_int32type endianess;
    u_int32type ndim;
    u_int32type box[10];
    u_int32type bytes_per_site;
    u_int64type sites;
  } fermiQCD_header = 
      {
	"File Type: MDP FIELD",
	"milc_qcd version 7",
	"",
	0x87654321,
	1,
	{nt},
	144*sizeof(mdp_complex),
	nt
      };

  FILE *fp;
  char *filename = param.savefile_c[pair];

  if(this_node != 0 || param.saveflag_c[pair] == FORGET )
    return NULL;

  /* Load any additional header values */
  strncpy(fermiQCD_header.creation_date, utc_date_time, 60);

  fp = fopen(filename,"wb");
  if(fp == NULL){
    printf("open_open_meson_file(%d): Can't open %s for writing\n",
	   this_node, filename);
    return NULL;
  }

  /* Write the header */
  if(fwrite(&fermiQCD_header, sizeof(fermiQCD_header), 1, fp) != 1){
    printf("open_open_meson_file(%d): Can't write header to %s\n",
	   this_node, filename);
    return NULL;
  }

  return fp;
}
/*--------------------------------------------------------------------*/
static int lookup_corr_index(int pair, int m){
  int i;

  /* Search corr_index table for first occurrence of correlator index m */
  for(i = 0; i < param.num_corr[pair]; i++){
    if(param.corr_index[pair][i] == m)
      return i;
  }
  printf("lookup_corr_index: Can't locate correlator %d for pair %d\n",
	 m,pair);
  terminate(1);
  return -1;
}
/*--------------------------------------------------------------------*/
static void print_start_fnal_meson_prop(FILE *fp, int pair, int m)
{
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  int ip0 = get_eve(iq0);
  int ip1 = get_eve(iq1);
  int i   = lookup_corr_index(pair,m);

  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;

  fprintf(fp,"---\n");
  fprintf(fp,"correlator:                   %s\n",param.meson_label[pair][m]);
  fprintf(fp,"momentum:                     %s\n",param.mom_label[pair][m]);
  fprintf(fp,"gamma_source:                 %s\n",gamma_label(param.gam_src[pair][i]));
  fprintf(fp,"gamma_sink:                   %s\n",gamma_label(param.gam_snk[pair][i]));

  /* Print correlator key encoding metadata */
  fprintf(fp,"correlator_key:               %s", param.meson_label[pair][m]);

  /* Source labels */
  if(strlen(param.src_qs[ip0].label)>0)
    fprintf(fp,"_%s", param.src_qs[ip0].label);

  if(strlen(param.src_qs[ip1].label)>0)
    fprintf(fp,"_%s", param.src_qs[ip1].label);

  /* Sink labels */
  if(strlen(param.snk_qs_op[iq0].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq0].label);

  if(strlen(param.snk_qs_op[iq1].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq1].label);

  /* Mass or kappa labels */
  if(param.prop_type[ip0] == CLOVER_TYPE)
    fprintf(fp,"_k%s", param.kappa_label[ip0]);
  else /* KS_TYPE */
    fprintf(fp,"_m%s", param.mass_label[ip0]);
  if(param.prop_type[ip1] == CLOVER_TYPE)
    fprintf(fp,"_k%s", param.kappa_label[ip1]);
  else /* KS_TYPE */
    fprintf(fp,"_m%s", param.mass_label[ip1]);
  fprintf(fp, "_%s\n", param.mom_label[pair][m]);

  fprintf(fp,"...\n");
}
		       
/*--------------------------------------------------------------------*/
static void print_start_baryon_prop(int iq0, int iq1, int iq2, char sink[])
{
  int ip0 = param.prop_for_qk[iq0];
  int ip1 = param.prop_for_qk[iq1];
  int ip2 = param.prop_for_qk[iq2];

  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("SOURCES: %s %s %s\n",param.src_qs[ip0].descrp,
	 param.src_qs[ip1].descrp, param.src_qs[ip2].descrp );
  printf("KAPPAS: %g %g %g\n",param.dcp[ip0].Kappa,
	 param.dcp[ip1].Kappa, param.dcp[ip2].Kappa);
  /* Note, the metadata should be updated here to handle the case
     of nontrivial sink operators */
  printf("SINKS: %s %s %s %s\n",param.snk_qs_op[iq0].descrp,
	 param.snk_qs_op[iq1].descrp, param.snk_qs_op[iq2].descrp, sink );
}
/*--------------------------------------------------------------------*/
static void print_meson_prop(int pair, int t, dcomplex c)
{
  if(param.saveflag_c[pair] != FORGET)return;
  node0_printf("%d %e %e\n", t, c.real, c.imag);
}
/*--------------------------------------------------------------------*/
static void print_fnal_meson_prop(FILE *fp, int pair, int t, dcomplex c)
{
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fprintf(fp, "%d\t%e\t%e\n", t, c.real, c.imag);
}
/*--------------------------------------------------------------------*/
static void write_open_meson_prop(FILE *fp, int pair, int t, 
				  dwilson_propagator *wp)
{
  mdp_complex c[4][4][3][3];
  int c0,c1,s0,s1;
  
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;

  /* Remap Wilson propagator elements */
  for(c0=0;c0<3;c0++)
    for(s0=0;s0<4;s0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++){
	  c[s0][s1][c0][c1].real = wp->c[c0].d[s0].d[s1].c[c1].real;
	  c[s0][s1][c0][c1].imag = wp->c[c0].d[s0].d[s1].c[c1].imag;
	}

  /* Write the open propagator element c at time t */
  if(fwrite(&c, sizeof(c), 1, fp) != 1){
    printf("write_open_meson_prop(%d): Error writing open meson correlator\n",
	   this_node);
  }
}
/*--------------------------------------------------------------------*/
static void print_end_prop(int pair){
  if(param.saveflag_c[pair] != FORGET)return;
  node0_printf("ENDPROP\n");
}
/*--------------------------------------------------------------------*/
static void print_end_fnal_meson_prop(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  //  fprintf(fp, "&\n");
}
/*--------------------------------------------------------------------*/
static void close_fnal_meson_file(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  if(fp != NULL)fclose(fp);
}
/*--------------------------------------------------------------------*/
static void close_open_meson_file(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  if(fp != NULL)fclose(fp);
}
/*--------------------------------------------------------------------*/
static void g_wpropsum(wilson_propagator *wp){
  int c,s;

  for(c=0;c<3;c++)
    for(s=0;s<4;s++)
      g_wvectorsumfloat(&wp->c[c].d[s]);
}
/*--------------------------------------------------------------------*/
static void divreal_wprop(wilson_propagator *src, Real x, 
			  wilson_propagator *dest){
  int c0,s0,c1,s1;
  for(c0=0;c0<3;c0++)
    for(s0=0;s0<4;s0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++){
	  CDIVREAL(src->c[c0].d[s0].d[s1].c[c1],x,
		   dest->c[c0].d[s0].d[s1].c[c1]);
	}
}
/*--------------------------------------------------------------------*/
static void mulreal_dwprop(dwilson_propagator *src, Real x, 
			   dwilson_propagator *dest){
  int c0,s0,c1,s1;
  for(c0=0;c0<3;c0++)
    for(s0=0;s0<4;s0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++){
	  CMULREAL(src->c[c0].d[s0].d[s1].c[c1],x,
		   dest->c[c0].d[s0].d[s1].c[c1]);
	}
}
/*--------------------------------------------------------------------*/
static void copy_wprop(dwilson_propagator *dest, wilson_propagator *src){
  int c0,s0,c1,s1;
  for(c0=0;c0<3;c0++)
    for(s0=0;s0<4;s0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++){
	  dest->c[c0].d[s0].d[s1].c[c1].real = src->c[c0].d[s0].d[s1].c[c1].real;
	  dest->c[c0].d[s0].d[s1].c[c1].imag = src->c[c0].d[s0].d[s1].c[c1].imag;
	}

}
/*--------------------------------------------------------------------*/
static double get_meson_scale_factor(int iq0, int iq1){
  int ip0 = get_eve(iq0);
  int ip1 = get_eve(iq1);
  double scale_factor = 1;

  if(param.prop_type[ip0] == CLOVER_TYPE)
    scale_factor *= param.src_qs[ip0].scale_fact;
  if(param.prop_type[ip1] == CLOVER_TYPE)
    scale_factor *= param.src_qs[ip1].scale_fact;

  return 1./scale_factor;
}
/*--------------------------------------------------------------------*/
static double get_baryon_scale_factor(int iq0, int iq1, int iq2){
  int ip0 = get_eve(iq0);
  int ip1 = get_eve(iq1);
  int ip2 = get_eve(iq2);
  double scale_factor = 1;

  if(param.prop_type[ip0] == CLOVER_TYPE)
    scale_factor *= param.src_qs[ip0].scale_fact;
  if(param.prop_type[ip1] == CLOVER_TYPE)
    scale_factor *= param.src_qs[ip1].scale_fact;
  if(param.prop_type[ip2] == CLOVER_TYPE)
    scale_factor *= param.src_qs[ip2].scale_fact;

  return 1./scale_factor;
}
/*--------------------------------------------------------------------*/
static void spectrum_cl_print_diag(int pair){

  Real space_vol;
  Real norm_fac;
  FILE *corr_fp;
  int t, tp;
  int m,b;
  int num_report = param.num_corr_report[pair];
  complex prop;
  dcomplex dprop;
  wilson_propagator wprop;
  dwilson_propagator dwprop;
  double meson_scale, baryon_scale;
  
  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);

  /* Rescaling */
  meson_scale = get_meson_scale_factor(param.qkpair[pair][0],  
				       param.qkpair[pair][1]);
  baryon_scale = get_baryon_scale_factor(param.qkpair[pair][0], 
					 param.qkpair[pair][0], 
					 param.qkpair[pair][0]);
  
  /* Point sink */
  if(param.do_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair);
    if(this_node == 0 && corr_fp == NULL && param.saveflag_c[pair] != FORGET)
      terminate(1);
    
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m, "POINT");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	prop = pmes_prop[m][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	dprop.real = prop.real;	dprop.imag = prop.imag;
	CMULREAL(dprop, meson_scale, dprop);
	print_meson_prop(pair, t, dprop);
	print_fnal_meson_prop(corr_fp, pair, t, dprop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);

    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
  
  /* Open-meson correlator */
  if(param.do_open_meson[pair]){
    corr_fp = open_open_meson_file(pair);

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	wprop = omes_prop[m][tp];
	g_wpropsum( &wprop );
	divreal_wprop(&wprop, norm_fac, &wprop);
	copy_wprop(&dwprop, &wprop);
	mulreal_dwprop(&dwprop, meson_scale, &dwprop);
	write_open_meson_prop(corr_fp, pair, t, &dwprop);
      }
    } /* mesons and momenta */
    close_open_meson_file(corr_fp, pair);
  }
  
  /* print baryon propagators */
  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind[b]);
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	prop = bar_prop[b][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	dprop.real = prop.real;	dprop.imag = prop.imag;
	CMULREAL(dprop, baryon_scale, dprop);
	node0_printf("%d %e %e\n",t, dprop.real,dprop.imag);
      }
      print_end_prop(pair);
    }
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_print_offdiag(int pair){
  
  Real space_vol;
  Real norm_fac;
  FILE *corr_fp;
  int t, tp;
  int m,b;
  int num_report = param.num_corr_report[pair];
  complex prop;
  dcomplex dprop;
  wilson_propagator wprop;
  dwilson_propagator dwprop;
  double meson_scale, baryon_scale1, baryon_scale2;
  
  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);

  /* Rescaling */
  meson_scale = get_meson_scale_factor(param.qkpair[pair][0],  
				       param.qkpair[pair][1]);
  baryon_scale1 = get_baryon_scale_factor(param.qkpair[pair][1], 
					  param.qkpair[pair][0], 
					  param.qkpair[pair][0]);
  baryon_scale2 = get_baryon_scale_factor(param.qkpair[pair][1], 
					  param.qkpair[pair][1], 
					  param.qkpair[pair][0]);
  
  /* Point sink */
  if(param.do_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair);
    if(this_node == 0 && corr_fp == NULL && param.saveflag_c[pair] != FORGET)
      terminate(1);

    /* print meson propagators */
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      print_start_meson_prop(pair, m, "POINT");
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	prop = pmes_prop[m][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	dprop.real = prop.real;	dprop.imag = prop.imag;
	CMULREAL(dprop, meson_scale, dprop);
	print_meson_prop(pair, t, dprop);
	print_fnal_meson_prop(corr_fp, pair, t, dprop);
      }
      print_end_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
    
    
  /* Open-meson correlator */
  if(param.do_open_meson[pair]){
    corr_fp = open_open_meson_file(pair);

    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	wprop = omes_prop[m][tp];
	g_wpropsum( &wprop );
	divreal_wprop(&wprop, norm_fac, &wprop);
	copy_wprop(&dwprop, &wprop);
	mulreal_dwprop(&dwprop, meson_scale, &dwprop);
	write_open_meson_prop(corr_fp, pair, t, &dwprop);
      }
    } /* mesons and momenta */
    close_open_meson_file(corr_fp, pair);
  }
  
  /* print baryon propagators */
  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][1], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	prop = bar_prop100[b][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	dprop.real = prop.real;	dprop.imag = prop.imag;
	CMULREAL(dprop, baryon_scale1, dprop);
	node0_printf("%d %e %e\n",t, dprop.real, dprop.imag);
      }
      print_end_prop(pair);
    }
    
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][1], 
			      param.qkpair[pair][1], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	tp = (t + param.t_offset[pair]) % nt;
	prop = bar_prop011[b][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	dprop.real = prop.real;	dprop.imag = prop.imag;
	CMULREAL(dprop, baryon_scale2, dprop);
	node0_printf("%d %e %e\n",t, dprop.real, dprop.imag);
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

  if(param.do_meson_spect[pair])
    destroy_mes_prop(&pmes_prop, num_corr);

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

  if(param.do_open_meson[pair])
    destroy_open_mes_prop(&omes_prop, num_corr);

  destroy_corr_table(pair);

  destroy_num_corr_occur(&num_corr_occur);
}

/*--------------------------------------------------------------------*/
void spectrum_cl(wilson_prop_field *qp0, wilson_prop_field *qp1, int pair)
{

  int *qkpair = param.qkpair[pair];
  double dtime;

  if(qp0->nc != qp1->nc){
    node0_printf("spectrum_cl: Inconsistent colors %d != %d\n", qp0->nc, qp1->nc);
    terminate(1);
  }
  spectrum_cl_init(pair);

  if(qkpair[0] == qkpair[1]){

    if(param.do_meson_spect[pair])
      spectrum_cl_diag_meson(qp0, pair);

    if(param.do_baryon_spect[pair])
      spectrum_cl_diag_baryon(qp0);

    if(param.do_open_meson[pair]){
      /* We must transpose source and sink color-spin indices */
      transpose_wp_field(qp0);
      spectrum_cl_offdiag_open_meson(qp0, qp0, pair);
      /* Transpose back again */
      transpose_wp_field(qp0);
    }

    dtime = start_timing();
    spectrum_cl_print_diag(pair);
    print_timing(dtime, "printing correlator");

  } else {

    if(param.do_meson_spect[pair])
      spectrum_cl_offdiag_meson(qp0, qp1, pair);

    if(param.do_baryon_spect[pair])
      spectrum_cl_offdiag_baryon(qp0, qp1);
    
    if(param.do_open_meson[pair]){
      /* We must transpose source and sink color-spin indices */
      transpose_wp_field(qp0);
      transpose_wp_field(qp1);
      spectrum_cl_offdiag_open_meson(qp0, qp1, pair);
      /* Transpose back again */
      transpose_wp_field(qp0);
      transpose_wp_field(qp1);
    }

    dtime = start_timing();
    spectrum_cl_print_offdiag(pair);
    print_timing(dtime, "printing correlator");
  }    
	    
  spectrum_cl_cleanup(pair);
}
