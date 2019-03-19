/***************** spectrum_ks.c *****************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

/* These procedures compute meson and baryon correlators.
   They also compute "open meson" correlators for later
   recombination.

*/

#include "ks_spectrum_includes.h"
#include "../include/fermion_links.h"
#include "../db/io_string_stream.h"
#ifdef HAVE_SQLITE
#include "../db/db.h"
#endif
#include <string.h>
#include <time.h>

static complex **baryon_prop = NULL;
static complex **pmes_prop = NULL;

/* Table of unique sink spin-taste indices */
static int *spin_taste_table, num_spin_taste_indices;

/* Table of unique momentum/parity combinations */
static int **momentum, num_mom_parities;
char **parity;

/* Hash table of correlators (indexed as they appear in param file)*/
static int *spin_taste_indices, *mom_parity_index;
static int **corr_table, *num_corr_mom_parity;

/* Count of correlator contributions */
static int *num_corr_occur;

/*--------------------------------------------------------------------*/
/* indexing is prop[meson_type][momentum][time] */

static complex ** 
create_hadron_prop(int ncor, int ntime){
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

static void 
destroy_hadron_prop(complex ***prop, int ncor){
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
  int num_corr = param.num_corr_m[pair];
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
/* Construct hash table for the spin-taste combinations.  Look for
   the spin-taste index in the hash table.  If found, return its index.
   Otherwise, add it to the table and assign it the new index */

static int 
hash_spin_taste_index(int *st, int st_test, int *n){
  int i;

  for(i = 0; i < *n; i++){
    if(st[i] == st_test)
      return i;
  }

  /* Add to table */

  st[*n] = st_test;

  (*n)++;
  return *n-1;
}

/*--------------------------------------------------------------------*/
/* Construct hash table of unique spin/taste indices                  */


static void 
create_spin_taste_index(int **st, int **st_index, int *n, int pair)
{
  char myname[] = "create_spin_tate_index";
  int c;
  int num_corr = param.num_corr_m[pair];

  *st = (int *)malloc(num_corr*sizeof(int));
  if(*st == NULL){
    printf("%s(%d): no room for spin-taste table\n", myname, this_node);
    terminate(1);
  }
  
  *st_index = (int *)malloc(num_corr*sizeof(int));
  if(*st_index == NULL){
    printf("%s(%d): no room for spin-taste index table\n", myname, this_node);
    terminate(1);
  }

  /* Run through all correlators for this pair, tabulating and hashing
     the unique spin-taste indices */
  for(c = 0; c < num_corr; c++){
    (*st_index)[c] = 
      hash_spin_taste_index(*st, param.spin_taste_snk[pair][c], n);
  }
}

/*--------------------------------------------------------------------*/

static void 
destroy_spin_taste_index(int **st, int **st_index)
{
  if(*st != NULL){
    free(*st); *st = NULL;
  }

  if(*st_index != NULL){
    free(*st_index);  *st_index = NULL;
  }
}

/*--------------------------------------------------------------------*/

static int *
create_num_corr_mom_parity(int *st_index, int n, int pair)
{
  char myname[] = "create_num_corr_mom_parity";
  int num_corr = param.num_corr_m[pair];
  int c, g;
  int *n_cmp;

  n_cmp = (int *)malloc(n*sizeof(int));
  if(n_cmp == NULL){
    printf("%s(%d): no room for num corr table\n",myname,this_node);
    terminate(1);
  }

  /* Count momentum/parity combinations for each unique sink spin-taste */

  for(g = 0; g < n; g++)
    n_cmp[g] = 0;

  for(c = 0; c < num_corr; c++){
    g = st_index[c];
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
   unique spin-taste assignments, and a list of the correlator
   indices for each pair */

static void
fill_corr_table(int *st, int *st_index, int **c_table, 
		int *n_cmp, int n, int *mp_index, int pair)
{
  char myname[] = "fill_corr_table";
  int c, g, k;
  int num_corr = param.num_corr_m[pair];
  int *k_corr_mom_parity;

  /* Set up counters, one for each unique sink spin-taste assignment */

  k_corr_mom_parity = (int *)malloc(n*sizeof(int));
  if(k_corr_mom_parity == NULL){
    printf("%s(%d): no room for counter\n",myname,this_node);
    terminate(1);
  }

  for(g = 0; g < n; g++)
    k_corr_mom_parity[g] = 0;


  /* Allocate rows of the correlator table, one row for each
     spin-taste assignment */

  for(g = 0; g < n; g++){
    c_table[g] = (int *)malloc(sizeof(int)*n_cmp[g]);
    if(c_table[g] == NULL){
      printf("%s(%d): no room for corr table\n",myname, this_node);
       terminate(1);
     }
  }

  /* Populate the rows of the correlator table */
  for(c = 0; c < num_corr; c++){
    g = st_index[c];    /* Hash of spin-taste assignment for this corr */
    k = k_corr_mom_parity[g];
    c_table[g][k] = c;  /* Point table to correlator */
    k_corr_mom_parity[g]++;
  }

  free(k_corr_mom_parity);
}

/*--------------------------------------------------------------------*/

/* Create a table of correlators for the current pair.  The table
   consists of a list of unique sink spin-taste indices, and a list of
   the correlator indices for each spin-taste assignment */

static void
create_corr_table(int pair)
{

  char myname[] = "create_corr_table";

  num_spin_taste_indices = 0;
  num_mom_parities = 0;

  create_spin_taste_index(&spin_taste_table, &spin_taste_indices,
		       &num_spin_taste_indices, pair);


  create_mom_parity(&momentum, &parity, &mom_parity_index, 
			  &num_mom_parities, pair);

  num_corr_mom_parity = 
    create_num_corr_mom_parity(spin_taste_indices, num_spin_taste_indices, 
			       pair);

  /* Create the correlator table */
  
  corr_table = (int **)malloc(num_spin_taste_indices*sizeof(int *));
  if(corr_table == NULL){
    printf("%s(%d): no room for corr table\n",myname, this_node);
    terminate(1);
  }

  /* Fill it */

  fill_corr_table(spin_taste_table, spin_taste_indices, corr_table,
		  num_corr_mom_parity, num_spin_taste_indices, 
		  mom_parity_index, pair);

}

/*--------------------------------------------------------------------*/
static void
destroy_corr_table(int pair)
{
  int g;

  destroy_spin_taste_index(&spin_taste_table, &spin_taste_indices);

  destroy_mom_parity(&momentum, &parity, &mom_parity_index,
		     num_mom_parities, pair);

  destroy_num_corr_mom_parity(&num_corr_mom_parity);

  if(corr_table != NULL){
    for(g = 0; g < num_spin_taste_indices; g++){
      if(corr_table[g] != NULL)
	free(corr_table[g]);
    }
    free(corr_table);
    corr_table = NULL;
  }
}

/*--------------------------------------------------------------------*/
/* The num_corr_occur array counts the number of contributions to
   a given meson correlator for normalization purposes    */

static int *
create_num_corr_occur(int pair){
  int c,m;
  int ncr = param.num_corr_report[pair];
  int num_corr = param.num_corr_m[pair];
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
spectrum_ks_init(int pair){
  
  int num_corr = 
    param.num_corr_report[pair]; /* number of meson correlators to store */

  if(param.do_meson_spect[pair])
    pmes_prop = create_hadron_prop(num_corr, nt);

  /* Set up count of contributions to the correlators */
  num_corr_occur = create_num_corr_occur(pair);
  
  /* Set up table of requested correlators */
  create_corr_table(pair);

}

/*--------------------------------------------------------------------*/
static void 
spectrum_ks_baryon_init(int triplet){
  
  int num_corr = param.num_corr_b[triplet];

  if(param.do_baryon_spect[triplet])
    baryon_prop = create_hadron_prop(num_corr, nt);
}

/*---------------------------------------------------------------------*/

static void spectrum_ks_baryon_nd(int nc, int type[], 
				  ks_prop_field *qp0, 
				  ks_prop_field *qp1,
				  ks_prop_field *qp2,
				  int triplet)
{
  double dtime = start_timing();

  ks_baryon_nd( baryon_prop, qp0, qp1, qp2, nc, type, 
		param.baryon_phase[triplet], param.baryon_factor[triplet]  );

  print_timing(dtime, "baryons");
}

/*--------------------------------------------------------------------*/

static void accum_gen_meson(complex **mp, su3_vector *qp0, int naik_index0,
			    su3_vector *qp1, int naik_index1, int pair){

  imp_ferm_links_t **fn = get_fm_links(fn_links);

  ks_meson_cont_mom(mp, qp0, qp1, 
		    num_mom_parities, momentum, parity,
		    num_spin_taste_indices, num_corr_mom_parity, corr_table, 
		    mom_parity_index, fn[naik_index0], fn[naik_index1], 
		    param.spin_taste_snk[pair], 
		    param.meson_phase[pair], param.meson_factor[pair],
		    param.corr_index[pair], &param.r_offset_m[pair][0]);
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_diag_gen_meson(ks_prop_field *qp, int naik_index,
				       complex **mp, int pair){
  
  int color;
  su3_vector *qps;
  double dtime = start_timing();

  for(color=0;color<qp->nc;color++){
    qps = qp->v[color];
    accum_gen_meson(mp, qps, naik_index, qps, naik_index, pair);
  }
  print_timing(dtime, "diagonal mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_diag_meson(ks_prop_field *qp, int naik_index, int pair){
  spectrum_ks_diag_gen_meson(qp, naik_index, pmes_prop, pair);
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_offdiag_gen_meson(ks_prop_field *qp0, 
					  int naik_index0,
					  ks_prop_field *qp1, 
					  int naik_index1,
					  complex **mp, int pair){

  su3_vector *qps0, *qps1;
  int color;
  double dtime = start_timing();

  if(qp0->nc != qp1->nc){
    printf("spectrum_ks_offdiag_gen_meson(%d): mismatched number of colors in pair %d: %d != %d\n", this_node, pair, qp0->nc, qp1->nc);
    terminate(1);
  }

  /* Off diagonal quark mass mesons */
  for(color = 0; color < qp0->nc; color++){
    qps0 = qp0->v[color];
    qps1 = qp1->v[color];

    /* Point meson propagator */
    accum_gen_meson(mp, qps0, naik_index0, qps1, naik_index1, pair);
  }
  print_timing(dtime, "offdiag mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_offdiag_meson(ks_prop_field *qp0, 
				      int naik_index0,
				      ks_prop_field *qp1,
				      int naik_index1,
				      int pair){

  /* Point meson */
  spectrum_ks_offdiag_gen_meson(qp0, naik_index0, qp1, naik_index1, pmes_prop, pair);
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



/*--------------------------------------------------------------------*/

int io_JSON_quark_source_sink_op(io_string_stream *json, quark_source_sink_op *qss_op); // TODO: prototype to header
int io_JSON_source_info(io_string_stream *json, quark_source *qs); // TODO: prototype to header
int io_JSON_field_op_info_list(io_string_stream *json, quark_source_sink_op *qss_op[], int n); // TODO: prototype to header

int extract_tsrc(int qtuple, char hadron)
{
  int tsrc[3];
  int nh0;
  int ih0[MAX_HISTORY];
  int *iq = NULL;
  int nquark = 0;
  switch(hadron) {
  case 'm': case 'M': // meson
    nquark = 2;
    iq = &param.qkpair[qtuple][0];
    break;
  case 'b': case 'B': // baryon
    nquark = 3;
    iq = &param.qktriplet[qtuple][0];
    break;
  default:
    fprintf(stderr,"ERROR: unknown hadron type\n");
    break;
  }

  int qi;
  for(qi=0; qi<nquark; ++qi)
    {
      int iq0 = 0;
      iq0 = iq[qi];
      int ip0 = get_ancestors(ih0, &nh0, iq0);
      int is0 = param.set[ip0];
      tsrc[qi] = param.src_qs[is0].t0;
    }
  if(tsrc[0] != tsrc[1])
    {
      fprintf(stderr, "WARNING: ambiguous mismatched tsrc values quarks! %d %d\n",tsrc[0],tsrc[1]);
    }
  if((nquark > 2) && (tsrc[0] != tsrc[2]))
    {
      fprintf(stderr, "WARNING: ambiguous mismatched tsrc values quarks! %d %d\n",tsrc[0],tsrc[2]);
    }
  return tsrc[0];
}

/* TODO: remove
int extract_tsrc_b(int qtriplet)
{
  int tsrc[3];
  int nh0;
  int ih0[MAX_HISTORY];
  int qi;
  for(qi=0; qi<3; ++qi)
    {
      int iq0 = param.qktriplet[qtriplet][qi];
      int ip0 = get_ancestors(ih0, &nh0, iq0);
      int is0 = param.set[ip0];
      tsrc[qi] = param.src_qs[is0].t0;
    }
  if( (tsrc[0] != tsrc[1]) || (tsrc[0] != tsrc[2]) )
    {
      fprintf(stderr, "WARNING: ambiguous mismatched tsrc values quarks! %d %d %d\n",tsrc[0],tsrc[1],tsrc[2]);
    }
  return tsrc[0];
}
*/

int io_JSON_quark_meta(io_string_stream *json, const char* qtype, char hadron, int qtuple, int qi)
{
  size_t len0 = strlen(json->base);
  int nh0;
  int ih0[MAX_HISTORY];
  int iq0 = 0;
  switch(hadron) {
  case 'm': case 'M': // meson
    iq0 = param.qkpair[qtuple][qi];
    break;
  case 'b': case 'B': // baryon
    iq0 = param.qktriplet[qtuple][qi];
    break;
  default:
    fprintf(stderr,"ERROR: unknown hadron type\n");
    break;
  }
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int is0 = param.set[ip0];
  io_JSON_begin_set(json); io_JSON_key(json,"type"); io_JSON_quoted(json,qtype);
  io_JSON_sep(json); io_JSON_key(json,"source"); io_JSON_source_info(json,&param.src_qs[is0]);
  io_JSON_sep(json); io_JSON_key(json,"sink"); io_JSON_begin_set(json);
  io_JSON_key(json,"label"); io_JSON_quoted(json,param.snk_qs_op[iq0].label);
  io_JSON_sep(json); io_JSON_key(json,"ops");
  {
    int k;
    quark_source_sink_op **op_list = static_cast(quark_source_sink_op**,malloc(sizeof(quark_source_sink_op*)*nh0));
    for(k = 0; k < nh0; k++) {
      op_list[k] = &param.snk_qs_op[ih0[nh0-1-k]];
    }
    io_JSON_field_op_info_list(json,op_list,nh0);
    free(op_list);
  } io_JSON_end_set(json);
  io_JSON_sep(json); io_JSON_key(json,"mass"); io_JSON_quoted(json,param.mass_label[ip0]);
  #if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  io_JSON_sep(json); io_JSON_key(json,"epsilon"); io_JSON_as_text(json,32,"%.6g",param.ksp[ip0].naik_term_epsilon);
  #endif
  #if U1_FIELD
  io_JSON_sep(json); io_JSON_key(json,"charge"); io_JSON_quoted(json,param.charge_label[is0]);
  #endif
  io_JSON_end_set(json);
  return(static_cast(int,strlen(json->base)-len0));
}

/*--------------------------------------------------------------------*/

int io_JSON_meson_unique_key(io_string_stream *json, int qpair, int jmeson)
{
  int nh0, nh1;
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY];
  int iq0 = param.qkpair[qpair][0];
  int iq1 = param.qkpair[qpair][1];
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];

  size_t len0 = strlen(json->base);
  io_JSON_as_text(json,strlen(param.meson_label[qpair][jmeson]),"%s",param.meson_label[qpair][jmeson]);
  int k;
  k = strlen(param.src_qs[is0].label)+1; if(k>0) io_JSON_as_text(json,k,"_%s",param.src_qs[is0].label);
  k = strlen(param.src_qs[is1].label)+1; if(k>0) io_JSON_as_text(json,k,"_%s",param.src_qs[is1].label);
  k = strlen(param.snk_qs_op[iq0].label)+1; if(k>0) io_JSON_as_text(json,k,"_%s",param.snk_qs_op[iq0].label);
  k = strlen(param.snk_qs_op[iq1].label)+1; if(k>0) io_JSON_as_text(json,k,"_%s",param.snk_qs_op[iq1].label);
  io_JSON_as_text(json,strlen(param.mass_label[ip0])+2,"_m%s", param.mass_label[ip0]);
  #if U1_FIELD
  io_JSON_as_text(json,strlen(param.charge_label[is0])+2,"_q%s",param.charge_label[is0]);
  #endif
  io_JSON_as_text(json,strlen(param.mass_label[ip1])+2,"_m%s", param.mass_label[ip1]);
  #if U1_FIELD
  io_JSON_as_text(json,strlen(param.charge_label[is1])+2,"_q%s",param.charge_label[is1]);
  #endif
  io_JSON_as_text(json,strlen(param.mom_label[qpair][jmeson])+2,"_%s",param.mom_label[qpair][jmeson]);
  return(static_cast(int,strlen(json->base)-len0));
}

static int lookup_corr_index(int pair, int m); // foward reference

int prepare_meson_JSON_metadata(io_string_stream *json, int qpair, int jmeson)
{
  size_t len0 = strlen(json->base);
  io_JSON_begin_set(json); io_JSON_key(json,"antiquark"); io_JSON_quark_meta(json,"staggered",'M',qpair,0);
  io_JSON_sep(json); io_JSON_key(json,"quark"); io_JSON_quark_meta(json,"staggered",'M',qpair,1);
  io_JSON_sep(json); io_JSON_key(json,"correlator"); io_JSON_quoted(json,param.meson_label[qpair][jmeson]);
  io_JSON_sep(json); io_JSON_key(json,"momentum"); io_JSON_quoted(json,param.mom_label[qpair][jmeson]);
  io_JSON_sep(json); io_JSON_key(json,"spin_taste_sink"); io_JSON_quoted(json,spin_taste_label(param.spin_taste_snk[qpair][jmeson]));
  io_JSON_sep(json); io_JSON_key(json,"correlator_key");
  io_JSON_quote(json); io_JSON_meson_unique_key(json,qpair,jmeson); io_JSON_quote(json);
  io_JSON_end_set(json);
  return(static_cast(int,strlen(json->base)-len0));
}

/* TODO: remove
int io_JSON_quark_meta_baryon(io_string_stream *json, const char* qtype, int triplet, int qi)
{
  size_t len0 = strlen(json->base);
  int iq0 = param.qktriplet[triplet][qi];
  int ih0[MAX_HISTORY];
  int nh0;
  int i;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int is0 = param.set[ip0];
  io_JSON_begin_set(json); io_JSON_key(json,"type"); io_JSON_quoted(json,qtype);
  io_JSON_sep(json); io_JSON_key(json,"source"); io_JSON_source_info(json,&param.src_qs[is0]);
  io_JSON_sep(json); io_JSON_key(json,"sink"); io_JSON_begin_set(json);
  io_JSON_key(json,"label"); io_JSON_quoted(json,param.snk_qs_op[iq0].label);
  io_JSON_sep(json); io_JSON_key(json,"ops");
  {
    int k;
    quark_source_sink_op **op_list = static_cast(quark_source_sink_op**,malloc(sizeof(quark_source_sink_op*)*nh0));
    for(k = 0; k < nh0; k++) {
      op_list[k] = &param.snk_qs_op[ih0[nh0-1-k]];
    }
    io_JSON_field_op_info_list(json,op_list,nh0);
    free(op_list);
  } io_JSON_end_set(json);
  io_JSON_sep(json); io_JSON_key(json,"mass"); io_JSON_quoted(json,param.mass_label[ip0]);
  #if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  io_JSON_sep(json); io_JSON_key(json,"epsilon"); io_JSON_as_text(json,32,"%.6g",param.ksp[ip0].naik_term_epsilon);
  #endif
  #if U1_FIELD
  io_JSON_sep(json); io_JSON_key(json,"charge"); io_JSON_quoted(json,param.charge_label[is0]);
  #endif
  io_JSON_end_set(json);
  return(static_cast(int,strlen(json->base)-len0));
}
*/

int io_JSON_baryon_unique_key(io_string_stream *json, int triplet, int jbaryon)
{
  size_t len0 = strlen(json->base);
  int iq0 = param.qktriplet[triplet][0];
  int iq1 = param.qktriplet[triplet][1];
  int iq2 = param.qktriplet[triplet][2];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY], ih2[MAX_HISTORY];
  int nh0, nh1, nh2;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int ip2 = get_ancestors(ih2, &nh2, iq2);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  int is2 = param.set[ip2];
  int len;

  io_JSON_as_text(json,strlen(param.baryon_label[triplet][jbaryon]),"%s",param.baryon_label[triplet][jbaryon]);
  // source labels
  len=strlen(param.src_qs[is0].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.src_qs[is0].label);
  len=strlen(param.src_qs[is1].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.src_qs[is1].label);
  len=strlen(param.src_qs[is2].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.src_qs[is2].label);
  // sink labels
  len=strlen(param.snk_qs_op[iq0].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.snk_qs_op[iq0].label);
  len=strlen(param.snk_qs_op[iq1].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.snk_qs_op[iq1].label);
  len=strlen(param.snk_qs_op[iq2].label); if(len>0) io_JSON_as_text(json,len+1,"_%s",param.snk_qs_op[iq2].label);
  // masses and charges
  io_JSON_as_text(json,strlen(param.mass_label[ip0])+1,"_%s",param.mass_label[ip0]);
  #if U1_FIELD
  io_JSON_as_text(json,strlen(param.charge_label[is0])+2,"_q%s",param.charge_label[is0]);
  #endif
  io_JSON_as_text(json,strlen(param.mass_label[ip1])+1,"_%s",param.mass_label[ip1]);
  #if U1_FIELD
  io_JSON_as_text(json,strlen(param.charge_label[is1])+2,"_q%s",param.charge_label[is1]);
  #endif
  io_JSON_as_text(json,strlen(param.mass_label[ip2])+1,"_%s",param.mass_label[ip2]);
  #if U1_FIELD
  io_JSON_as_text(json,strlen(param.charge_label[is2])+2,"_q%s",param.charge_label[is2]);
  #endif
  io_JSON_as_text(json,5,"%s","_p000");
  return(static_cast(int,strlen(json->base)-len0));
}

int prepare_baryon_JSON_metadata(io_string_stream *json, int triplet, int jbaryon)
{
  size_t len0 = strlen(json->base);
  io_JSON_begin_set(json); io_JSON_key(json,"quark0"); io_JSON_quark_meta(json,"staggered",'B',triplet,0);
  io_JSON_sep(json); io_JSON_key(json,"quark1"); io_JSON_quark_meta(json,"staggered",'B',triplet,1);
  io_JSON_sep(json); io_JSON_key(json,"quark2"); io_JSON_quark_meta(json,"staggered",'B',triplet,2);
  io_JSON_sep(json); io_JSON_key(json,"correlator"); io_JSON_quoted(json,param.baryon_label[triplet][jbaryon]);
  io_JSON_sep(json); io_JSON_key(json,"baryon_type"); io_JSON_quoted(json,baryon_type_label(param.baryon_type_snk[triplet][jbaryon]));
  io_JSON_sep(json); io_JSON_key(json,"momentum"); io_JSON_quoted(json,"p000");
  io_JSON_sep(json); io_JSON_key(json,"correlator_key");
  io_JSON_quote(json); io_JSON_baryon_unique_key(json,triplet,jbaryon); io_JSON_quote(json);
  io_JSON_end_set(json);
  return(static_cast(int,strlen(json->base)-len0));
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
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  FILE *fp;

  /* Only node 0 writes, and only if 'SAVE_ASCII' we want the file. */
  if(this_node != 0 || param.saveflag_m[pair] != SAVE_ASCII )
    return NULL;

  /* Always append */
  fp = fopen(param.savefile_m[pair],"a");
  if(fp == NULL){
    printf("open_fnal_meson_file: ERROR. Can't open %s\n",
	   param.savefile_m[pair]);
    return NULL;
  }
  fprintf(fp,"---\n");
  fprintf(fp,"JobID:                        %s\n",param.job_id);
  fprintf(fp,"series:                       %s\n",param.series);
  fprintf(fp,"trajectory:                   %d\n",param.trajectory);
  fprintf(fp,"date:                         \"%s UTC\"\n",utc_date_time);
  fprintf(fp,"lattice_size:                 %d,%d,%d,%d\n", nx, ny, nz, nt);
  //  fprintf(fp,"spatial volume:        %g\n",((float)nx)*ny*nz);

  fprintf(fp,"antiquark_type:               staggered\n");

//  {
//    quark_source_sink_op* qs_op = param.src_qs[is0].op;
  print_source_info(fp, "antiquark_source", &param.src_qs[is0]);
  fprintf(fp,"antiquark_source_label:       %s\n",param.src_qs[is0].label);
  print_field_op_info(fp, "antiquark_source", param.src_qs[is0].op);

//    while(qs_op != NULL){
//      print_field_op_info(fp, "antiquark_source", qs_op);
//      qs_op = qs_op->op;
//    }
//  }

//  fprintf(fp,"antiquark_sink_op:       %s",param.snk_qs_op[ih0[nh0-1]].descrp);
//  for(i = nh0-2; i >=0; i--)
//    fprintf(fp,"/%s",param.snk_qs_op[ih0[i]].descrp);
//  fprintf(fp,"\n");
  fprintf(fp,"antiquark_sink_label:         %s\n",param.snk_qs_op[iq0].label);
//  for(k = 0; k < nh0; k++){
//    int ih = ih0[nh0-1-k];
//    print_field_op_info(fp, "antiquark_sink", &param.snk_qs_op[ih]);
//  }
  
  {
    quark_source_sink_op **op_list = (quark_source_sink_op **)
      malloc(sizeof(quark_source_sink_op *)*nh0);
    for(k = 0; k < nh0; k++)
      op_list[k] = &param.snk_qs_op[ih0[nh0-1-k]];
    print_field_op_info_list(fp, "antiquark_sink", op_list, nh0);
    free(op_list);
  }
      
  fprintf(fp,"antiquark_mass:               \"%s\"\n", param.mass_label[ip0]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fprintf(fp,"antiquark_epsilon:            %g\n", param.ksp[ip0].naik_term_epsilon);
#endif
#if U1_FIELD
  fprintf(fp,"antiquark_charge:             \"%s\"\n", param.charge_label[is0]);
#endif

  fprintf(fp,"quark_type:                   staggered\n");

//  {
//    quark_source_sink_op *wqs_op = param.src_qs[is1].op;
  print_source_info(fp, "quark_source", &param.src_qs[is1]);
  fprintf(fp,"quark_source_label:           %s\n",param.src_qs[is1].label);
  print_field_op_info(fp, "quark_source", param.src_qs[is1].op);
//    while(wqs_op != NULL){
//      print_field_op_info(fp, "quark_source", wqs_op);
//      wqs_op = wqs_op->op;
//    }
//  }

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
//  for(k = 0; k < nh1; k++){
//    int ih = ih1[nh1-1-k];
//    print_field_op_info(fp, "quark_sink", &param.src_qs_op[ih]);
//  }

  fprintf(fp,"quark_mass:                   \"%s\"\n",param.mass_label[ip1]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fprintf(fp,"quark_epsilon:                %g\n",param.ksp[ip1].naik_term_epsilon);
#endif
#if U1_FIELD
  fprintf(fp,"quark_charge:                 \"%s\"\n", param.charge_label[is1]);
#endif
  fprintf(fp,"...\n");
  return fp;
}
		       
/*--------------------------------------------------------------------*/
static FILE* open_fnal_baryon_file(int triplet){
  int iq0 = param.qktriplet[triplet][0];
  int iq1 = param.qktriplet[triplet][1];
  int iq2 = param.qktriplet[triplet][2];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY], ih2[MAX_HISTORY];
  int nh0, nh1, nh2;
  int i;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int ip2 = get_ancestors(ih2, &nh2, iq2);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  int is2 = param.set[ip2];
  FILE *fp;

  /* Only node 0 writes, and only if we want the file. */
  if(this_node != 0 || param.saveflag_b[triplet] != SAVE_ASCII )
    return NULL;

  /* Always append */
  fp = fopen(param.savefile_b[triplet],"a");
  if(fp == NULL){
    printf("open_fnal_baryon_file: ERROR. Can't open %s\n",
	   param.savefile_b[triplet]);
    return NULL;
  }
  fprintf(fp,"---\n");
  fprintf(fp,"JobID:                       %s\n",param.job_id);
  fprintf(fp,"series:                      %s\n",param.series);
  fprintf(fp,"trajectory:                  %d\n",param.trajectory);
  fprintf(fp,"date:                        \"%s UTC\"\n",utc_date_time);
  fprintf(fp,"lattice_size:                %d,%d,%d,%d\n", nx, ny, nz, nt);

  fprintf(fp,"quark0_type:                 staggered\n");
  fprintf(fp,"quark0_source_type:          %s\n",param.src_qs[is0].descrp);
  fprintf(fp,"quark0_source_label:         %s\n",param.src_qs[is0].label);

  fprintf(fp,"quark0_sink_type:            %s",param.snk_qs_op[ih0[nh0-1]].descrp);
  for(i = nh0-2; i >=0; i--)
    fprintf(fp,"/%s",param.snk_qs_op[ih0[i]].descrp);
  fprintf(fp,"\n");
  fprintf(fp,"quark0_sink_label:           %s\n",param.snk_qs_op[iq0].label);
  
  fprintf(fp,"quark0_mass:                 \"%s\"\n",param.mass_label[ip0]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fprintf(fp,"quark0_epsilon:              %g\n",param.ksp[ip0].naik_term_epsilon);
#endif
#if U1_FIELD
  fprintf(fp,"quark0_charge:               \"%s\"\n", param.charge_label[is0]);
#endif

  fprintf(fp,"quark1_type:                 stagered\n");
  fprintf(fp,"quark1_source_type:          %s\n",param.src_qs[is1].descrp);
  fprintf(fp,"quark1_source_label:         %s\n",param.src_qs[is1].label);

  fprintf(fp,"quark1_sink_type:            %s",param.snk_qs_op[ih1[nh1-1]].descrp);
  for(i = nh1-2; i >=0; i--)
    fprintf(fp,"/%s",param.snk_qs_op[ih1[i]].descrp);
  fprintf(fp,"\n");
  fprintf(fp,"quark1_sink_label:           %s\n",param.snk_qs_op[iq1].label);

  fprintf(fp,"quark1_mass:                 \"%s\"\n",param.mass_label[ip1]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fprintf(fp,"quark1_epsilon:              %g\n",param.ksp[ip1].naik_term_epsilon);
#endif
#if U1_FIELD
  fprintf(fp,"quark1_charge:               \"%s\"\n", param.charge_label[is1]);
#endif

  fprintf(fp,"quark2_type:                 staggered\n");
  fprintf(fp,"quark2_source_type:          %s\n",param.src_qs[is2].descrp);
  fprintf(fp,"quark2_source_label:         %s\n",param.src_qs[is2].label);

  fprintf(fp,"quark2_sink_type:            %s",param.snk_qs_op[ih2[nh2-1]].descrp);
  for(i = nh2-2; i >=0; i--)
    fprintf(fp,"/%s",param.snk_qs_op[ih2[i]].descrp);
  fprintf(fp,"\n");
  fprintf(fp,"quark2_sink_label:           %s\n",param.snk_qs_op[iq2].label);
  
  fprintf(fp,"quark2_mass:                 \"%s\"\n",param.mass_label[ip2]);
#if ( FERM_ACTION == HISQ || FERM_ACTION == HYPISQ )
  fprintf(fp,"quark2_epsilon:              %g\n",param.ksp[ip2].naik_term_epsilon);
#endif
#if U1_FIELD
  fprintf(fp,"quark2_charge:               \"%s\"\n", param.charge_label[is2]);
#endif

  fprintf(fp,"...\n");
  return fp;
}
		       
/*--------------------------------------------------------------------*/
static int lookup_corr_index(int pair, int m){
  int i;

  /* Search corr_index table for first occurrence of correlator index m */
  for(i = 0; i < param.num_corr_m[pair]; i++){
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
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY];
  int nh0, nh1;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  int i   = lookup_corr_index(pair,m);

  if(this_node != 0 || param.saveflag_m[pair] != SAVE_ASCII)return;

  fprintf(fp,"---\n");
  fprintf(fp,"correlator:                   %s\n",param.meson_label[pair][m]);
  fprintf(fp,"momentum:                     %s\n",param.mom_label[pair][m]);
  fprintf(fp,"spin_taste_sink:              %s\n",
	  spin_taste_label(param.spin_taste_snk[pair][i]));

  /* Print correlator key encoding metadata */
  fprintf(fp,"correlator_key:               %s", param.meson_label[pair][m]);

  /* Source labels */
  if(strlen(param.src_qs[is0].label)>0)
    fprintf(fp,"_%s", param.src_qs[is0].label);

  if(strlen(param.src_qs[is1].label)>0)
    fprintf(fp,"_%s", param.src_qs[is1].label);

  /* Sink labels */
  if(strlen(param.snk_qs_op[iq0].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq0].label);

  if(strlen(param.snk_qs_op[iq1].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq1].label);

  /* Mass labels */
  fprintf(fp,"_m%s", param.mass_label[ip0]);
#if U1_FIELD
  fprintf(fp,"_q%s", param.charge_label[is0]);
#endif
  fprintf(fp,"_m%s", param.mass_label[ip1]);
#if U1_FIELD
  fprintf(fp,"_q%s", param.charge_label[is1]);
#endif
  fprintf(fp, "_%s\n", param.mom_label[pair][m]);

  fprintf(fp,"...\n");
}
		       
/*--------------------------------------------------------------------*/
static void 
print_start_meson_prop(int pair, int m){
  if(param.saveflag_m[pair] != FORGET)return;
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY];
  int nh0, nh1;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("MOMENTUM: %s\n", param.mom_label[pair][m]);
  printf("MASSES: ");
  printf("%g ", param.ksp[ip0].mass);
  printf("%g ", param.ksp[ip1].mass);
  printf("\n");
#ifdef U1_FIELD
  printf("CHARGES: ");
  printf("%g ", param.charge[is0]);
  printf("%g ", param.charge[is1]);
  printf("\n");
#endif
  printf("SOURCE: %s %s\n",param.src_qs[is0].descrp,
	 param.src_qs[is1].descrp);
  printf("SINKOPS: %s %s\n", param.snk_qs_op[iq0].descrp, 
	 param.snk_qs_op[iq1].descrp);
  printf("SINKS: %s\n", param.meson_label[pair][m]);
}
		       
/*--------------------------------------------------------------------*/
static void print_start_fnal_baryon_prop(FILE *fp, int triplet, int b)
{
  int iq0 = param.qktriplet[triplet][0];
  int iq1 = param.qktriplet[triplet][1];
  int iq2 = param.qktriplet[triplet][2];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY], ih2[MAX_HISTORY];
  int nh0, nh1, nh2;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int ip2 = get_ancestors(ih2, &nh2, iq2);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  int is2 = param.set[ip2];

  if(this_node != 0 || param.saveflag_b[triplet] != SAVE_ASCII)return;

  fprintf(fp,"---\n");
  fprintf(fp,"correlator:                  %s\n",param.baryon_label[triplet][b]);
  fprintf(fp,"baryon_type:                 %s\n",
	  baryon_type_label(param.baryon_type_snk[triplet][b]));

  /* Print correlator key encoding metadata */
  fprintf(fp,"correlator_key:              %s", param.baryon_label[triplet][b]);

  /* Source labels */
  if(strlen(param.src_qs[is0].label)>0)
    fprintf(fp,"_%s", param.src_qs[is0].label);

  if(strlen(param.src_qs[is1].label)>0)
    fprintf(fp,"_%s", param.src_qs[is1].label);

  if(strlen(param.src_qs[is2].label)>0)
    fprintf(fp,"_%s", param.src_qs[is2].label);

  /* Sink labels */
  if(strlen(param.snk_qs_op[iq0].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq0].label);

  if(strlen(param.snk_qs_op[iq1].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq1].label);

  if(strlen(param.snk_qs_op[iq2].label)>0)
    fprintf(fp,"_%s",param.snk_qs_op[iq2].label);

  /* Mass labels */
  fprintf(fp,"_m%s", param.mass_label[ip0]);
#if U1_FIELD
  fprintf(fp,"_q%s", param.charge_label[is0]);
#endif
  fprintf(fp,"_m%s", param.mass_label[ip1]);
#if U1_FIELD
  fprintf(fp,"_q%s", param.charge_label[is1]);
#endif
  fprintf(fp,"_m%s", param.mass_label[ip2]);
#if U1_FIELD
  fprintf(fp,"_q%s", param.charge_label[is2]);
#endif
  fprintf(fp, "\n");

  fprintf(fp,"...\n");
}
		       
/*--------------------------------------------------------------------*/
static void print_start_baryon_prop(int triplet, int b)
{
  if(param.saveflag_b[triplet] != FORGET)return;
  int iq0 = param.qktriplet[triplet][0];
  int iq1 = param.qktriplet[triplet][1];
  int iq2 = param.qktriplet[triplet][2];
  int ih0[MAX_HISTORY], ih1[MAX_HISTORY], ih2[MAX_HISTORY];
  int nh0, nh1, nh2;
  int ip0 = get_ancestors(ih0, &nh0, iq0);
  int ip1 = get_ancestors(ih1, &nh1, iq1);
  int ip2 = get_ancestors(ih2, &nh2, iq2);
  int is0 = param.set[ip0];
  int is1 = param.set[ip1];
  int is2 = param.set[ip2];

  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("SOURCE: %s %s %s\n",param.src_qs[is0].descrp,
	 param.src_qs[is1].descrp, param.src_qs[is2].descrp );
  printf("MASSES: %g %g %g\n",param.ksp[ip0].mass,
	 param.ksp[ip1].mass, param.ksp[ip2].mass);
  /* Note, the metadata should be updated here to handle the case
     of nontrivial sink operators */
  printf("SINKOPS: %s %s %s\n", param.snk_qs_op[iq0].descrp, 
	 param.snk_qs_op[iq1].descrp, param.snk_qs_op[iq2].descrp);
  printf("SINKS: %s\n", param.baryon_label[triplet][b] );
}
/*--------------------------------------------------------------------*/
static void print_meson_prop(int pair, int t, complex c)
{
  if(param.saveflag_m[pair] == FORGET) { node0_printf("%d %e %e\n",t,(double)c.real,(double)c.imag); }
}
/*--------------------------------------------------------------------*/
static void print_baryon_prop(int triplet, int t, complex c)
{
  if(param.saveflag_b[triplet] == FORGET ) { node0_printf("%d %e %e\n",t,(double)c.real,(double)c.imag); }
}
/*--------------------------------------------------------------------*/
static void print_fnal_meson_prop(FILE *fp, int pair, int t, complex c)
{
  if(this_node == 0 && param.saveflag_m[pair] == SAVE_ASCII) {
    fprintf(fp, "%d\t%e\t%e\n", t, (double)c.real, (double)c.imag); }
}
/*--------------------------------------------------------------------*/
static void print_fnal_baryon_prop(FILE *fp, int triplet, int t, complex c)
{
  if(this_node == 0 && param.saveflag_b[triplet] == SAVE_ASCII) {
    fprintf(fp, "%d\t%e\t%e\n", t, (double)c.real, (double)c.imag); }
}
/*--------------------------------------------------------------------*/
static void print_end_meson_prop(int pair){
  if(param.saveflag_m[pair] == FORGET) { node0_printf("ENDPROP\n"); }
}
/*--------------------------------------------------------------------*/
static void print_end_baryon_prop(int pair){
  if(param.saveflag_b[pair] == FORGET) { node0_printf("ENDPROP\n"); }
}
/*--------------------------------------------------------------------*/
static void print_end_fnal_meson_prop(FILE *fp, int pair){
  //if(this_node != 0 || param.saveflag_m[pair] == FORGET)return;
  //  fprintf(fp, "&\n");
}
/*--------------------------------------------------------------------*/
static void close_fnal_meson_file(FILE *fp, int pair){
  if(this_node == 0 && param.saveflag_m[pair] == SAVE_ASCII)
    if(fp != NULL)fclose(fp);
}
/*--------------------------------------------------------------------*/
static void close_fnal_baryon_file(FILE *fp, int triplet){
  if(this_node == 0 && param.saveflag_b[triplet] == SAVE_ASCII)
    if(fp != NULL)fclose(fp);
}

/*--------------------------------------------------------------------*/
#ifdef HAVE_SQLITE

int sql_spectrum_ks(sqlite3 *db, int pair, unsigned long epoch_secs)
{
  if(!param.do_meson_spect[pair]) return(0);

  complex *prop = NULL;
  io_string_stream ckey; ckey.base=NULL;ckey.length=0;
  // JSON arrays of re and im parts of the correlator
  io_string_stream json_re; json_re.base=NULL;json_re.length=0;
  io_string_stream json_im; json_im.base=NULL;json_im.length=0;
  if(this_node == 0){
    prop = static_cast(complex*,calloc(nt,sizeof(complex))); // corr vector
    const int max_keysize = 1000; // starting key size
    io_string_stream_alloc(&ckey,max_keysize); // key space
    const int prec = 6; // default precision for format %e
    const int elem_sz = prec + 7; //  v: +f.ppppppe+dd
    const int sz = nt*(elem_sz + 1)+3;
    io_string_stream_alloc(&json_re, sz); // starting JSON array size for string "[v,v,...,v]"
    io_string_stream_alloc(&json_im, sz);
  }

  int num_report = param.num_corr_report[pair];
  int rc = SQLITE_OK;
  int m;
  for(m=0;m<num_report;m++) {
    // empty correlator object
    db_correlator corr; corr.id = 0; corr.name = NULL; corr.metadata = NULL;
    int free_corr = 0;
    if(this_node == 0)
      {
	// make correlator unique key
	ckey.base[0] = '\0'; // reset string
	io_JSON_meson_unique_key(&ckey,pair,m);

	// lookup key
	rc = db_query_correlator_by_name(db, ckey.base, &corr);
	if(rc != SQLITE_OK) break;
	free_corr = 1; // do free corr pointers
	if(corr.id < 1)
	  {
	    // ckey not found in db, so insert as a new correlator
	    free_corr = 0; // do not free corr pointers
	    corr.name = ckey.base; //copy pointer we manage
	    // generate correlator metadata
	    const int est_meta_size = 4000; // starting size
	    io_string_stream meta; meta.base=NULL;meta.length=0;
	    io_string_stream_alloc(&meta,est_meta_size);
	    prepare_meson_JSON_metadata(&meta,pair,m);
	    corr.metadata = meta.base; // copy pointer we manage
	    // insert new correlator
	    rc = db_insert_correlator(db,&corr);
	    if(rc != SQLITE_OK) break;
	    printf("INSERT name=%s id=%d\n",corr.name,corr.id);
	  }
      }

    // corr complex vector
    Real norm_fac = num_corr_occur[m];
    int t, tp;
    for(t=0; t<nt; t++)
      {
	tp = (t + param.r_offset_m[pair][3]) % nt;
	complex prop_t = pmes_prop[m][tp];
	g_complexsum( &prop_t );
	CDIVREAL(prop_t, norm_fac, prop_t);
	if(this_node == 0) prop[t] = prop_t;
      }

    if(this_node == 0)
      {
	int tsrc = extract_tsrc(pair,'M');
	json_re.base[0] = '\0'; // reset string
	json_im.base[0] = '\0';
	io_JSON_complex_array(&json_re, &json_im, prop, nt);

	db_data data; data.id=0; // unique, will be assigned by db
	data.correlator_id = corr.id;
	data.series=param.series; // copy char pointer
	data.trajectory=param.trajectory;
	data.tsrc=tsrc;
	data.jobid=param.job_id; // copy char pointer
	data.timestamp=epoch_secs;
	data.c_re=json_re.base ;data.c_im=json_im.base; // copy pointers

	// INSERT OR REPLACE numeric data
	rc = db_update_data(db, &data);
	if(free_corr) db_correlator_free(&corr);
	if(rc != SQLITE_OK) break;
      }
  }
  if(this_node == 0)
    {
      // free resources
      io_string_stream_free(&json_im);
      io_string_stream_free(&json_re);
      io_string_stream_free(&ckey);
      free(prop);
    }
  return(rc);
}

int sql_spectrum_ks_baryon(sqlite3 *db, int triplet, unsigned long epoch_secs)
{
  int rc = SQLITE_OK;
  if(!param.do_baryon_spect[triplet]) return(rc);

  complex *prop = NULL;
  io_string_stream ckey; ckey.base=NULL;ckey.length=0;
  // JSON arrays of re and im parts of the correlator
  io_string_stream json_re; json_re.base=NULL;json_re.length=0;
  io_string_stream json_im; json_im.base=NULL;json_im.length=0;
  if(this_node == 0){
    prop = static_cast(complex*,calloc(nt,sizeof(complex))); // corr vector
    const int max_keysize = 1000; // starting key size
    io_string_stream_alloc(&ckey,max_keysize); // key space
    const int prec = 6; // default precision for format %e
    const int elem_sz = prec + 7; //  v: +f.ppppppe+dd
    const int sz = nt*(elem_sz + 1)+3;
    io_string_stream_alloc(&json_re, sz); // starting JSON array size for string "[v,v,...,v]"
    io_string_stream_alloc(&json_im, sz);
  }

  int num_report = param.num_corr_b[triplet];
  int m;
  for(m=0;m<num_report;m++) {
    // empty correlator object
    db_correlator corr; corr.id = 0; corr.name = NULL; corr.metadata = NULL;
    int free_corr = 0;
    if(this_node == 0)
      {
	// make correlator unique key
	ckey.base[0] = '\0'; // reset string
	io_JSON_baryon_unique_key(&ckey,triplet,m);

	// lookup key
	db_query_correlator_by_name(db, ckey.base, &corr);
	free_corr = 1; // do free corr pointers
	if(corr.id < 1)
	  {
	    // ckey not found in db, so insert as a new correlator
	    free_corr = 0; // do not free corr pointers
	    corr.name = ckey.base; //copy pointer we manage
	    // generate correlator metadata
	    const int est_meta_size = 4000; // starting size
	    io_string_stream meta; meta.base=NULL;meta.length=0;
	    io_string_stream_alloc(&meta,est_meta_size);
	    prepare_baryon_JSON_metadata(&meta,triplet,m);
	    corr.metadata = meta.base; // copy pointer we manage
	    // insert new correlator
	    rc = db_insert_correlator(db,&corr);
	    if(rc != SQLITE_OK) break;
	    printf("INSERT name=%s id=%d\n",corr.name,corr.id);
	  }
      }

    // corr complex vector
    int t, tp;
    int r_off_b = param.r_offset_b[triplet][3];
    for(t=0; t<nt; t++)
      {
	tp = (t + param.r_offset_b[triplet][3]) % nt;
	complex prop_t = baryon_prop[m][tp];
	g_complexsum( &prop_t );
	// fix sign for antiperiodic bc
	if( ( ((t+r_off_b)/nt - r_off_b/nt) % 2 ) == 1 ){ CMULREAL(prop_t,-1.,prop_t) };
	if(this_node == 0) prop[t] = prop_t;
      }

    if(this_node == 0)
      {
	int tsrc = extract_tsrc(triplet,'B');
	json_re.base[0] = '\0'; // reset string
	json_im.base[0] = '\0';
	io_JSON_complex_array(&json_re, &json_im, prop, nt);

	db_data data; data.id=0; // unique, will be assigned by db
	data.correlator_id = corr.id;
	data.series=param.series; // copy char pointer
	data.trajectory=param.trajectory;
	data.tsrc=tsrc;
	data.jobid=param.job_id; // copy char pointer
	data.timestamp=epoch_secs;
	data.c_re=json_re.base ;data.c_im=json_im.base; // copy pointers

	// INSERT OR REPLACE numeric data
	rc = db_update_data(db, &data);
	if(free_corr) db_correlator_free(&corr);
	if(rc != SQLITE_OK) break;
      }
  }
  if(this_node == 0)
    {
      // free resources
      io_string_stream_free(&json_im);
      io_string_stream_free(&json_re);
      io_string_stream_free(&ckey);
      free(prop);
    }
  return(0);
}

#endif

/*--------------------------------------------------------------------*/
static void spectrum_ks_print_diag(int pair){

  //  Real space_vol;
  Real norm_fac;
  FILE *corr_fp;
  int t, tp;
  int m;
  int num_report = param.num_corr_report[pair];
  complex prop;
  
  /* Normalization factor */
  //  space_vol = (Real)(nx*ny*nz);
  
  /* Point sink */
  if(param.do_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair);
    
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];

      print_start_meson_prop(pair, m);
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	tp = (t + param.r_offset_m[pair][3]) % nt;
	prop = pmes_prop[m][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_meson_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);

    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_print_offdiag(int pair){
  
  //  Real space_vol;
  Real norm_fac;
  int t, tp;
  int m;
  int num_report = param.num_corr_report[pair];
  complex prop;
  FILE *corr_fp;
  
  /* Normalization factor */
  //  space_vol = (Real)(nx*ny*nz);

  /* Point sink */
  if(param.do_meson_spect[pair]){
    corr_fp = open_fnal_meson_file(pair);

    /* print meson propagators */
    for(m=0;m<num_report;m++) {
      norm_fac = num_corr_occur[m];
      
      print_start_meson_prop(pair, m);
      print_start_fnal_meson_prop(corr_fp, pair, m);
      for(t=0; t<nt; t++){
	tp = (t + param.r_offset_m[pair][3]) % nt;
	prop = pmes_prop[m][tp];
	g_complexsum( &prop );
	CDIVREAL(prop, norm_fac, prop);
	print_meson_prop(pair, t, prop);
	print_fnal_meson_prop(corr_fp, pair, t, prop);
      }
      print_end_meson_prop(pair);
      print_end_fnal_meson_prop(corr_fp, pair);
    } /* mesons and momenta */
    close_fnal_meson_file(corr_fp, pair);
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_ks_print_baryon(int triplet){
  
  //  Real space_vol;
  FILE *corr_fp;
  int t, tp;
  int b;
  complex prop;
  int num_corr = param.num_corr_b[triplet];
  
  /* Normalization factor */
  //  space_vol = (Real)(nx*ny*nz);

  /* print baryon propagator */
  if(param.do_baryon_spect[triplet]){
    corr_fp = open_fnal_baryon_file(triplet);

    for(b=0;b<num_corr;b++){

      print_start_baryon_prop(triplet, b);
      print_start_fnal_baryon_prop(corr_fp, triplet, b);
      for(t=0; t<nt; t++){
	tp = (t + param.r_offset_b[triplet][3]) % nt;
	prop = baryon_prop[b][tp];
	g_complexsum( &prop );
	// CDIVREAL(prop, space_vol, prop);
	/* Fix sign for antiperiodic bc */
	if( (((t+param.r_offset_b[triplet][3])/nt
	      - param.r_offset_b[triplet][3]/nt) %2 ) == 1 ){
	  CMULREAL(prop,-1.,prop);
	}
	print_baryon_prop(triplet, t, prop);
	print_fnal_baryon_prop(corr_fp, triplet, t, prop);
      }
      print_end_baryon_prop(triplet);
    }
    close_fnal_baryon_file(corr_fp, triplet);
  }
}

/*--------------------------------------------------------------------*/
void spectrum_ks_cleanup(int pair){
  int num_corr = param.num_corr_report[pair];

  if(param.do_meson_spect[pair])
    destroy_hadron_prop(&pmes_prop, num_corr);

  destroy_corr_table(pair);

  destroy_num_corr_occur(&num_corr_occur);
}

/*--------------------------------------------------------------------*/
void spectrum_ks_baryon_cleanup(int triplet){
  int num_corr = param.num_corr_b[triplet];

  if(param.do_baryon_spect[triplet])
    destroy_hadron_prop(&baryon_prop, num_corr);
}

/*--------------------------------------------------------------------*/
void spectrum_ks(ks_prop_field *qp0, int naik_index0, 
		 ks_prop_field *qp1, int naik_index1, int pair)
{
  int ret;
  int *qkpair = param.qkpair[pair];
  double dtime;

  spectrum_ks_init(pair);

  if(param.do_meson_spect[pair]){
    if(qkpair[0] == qkpair[1]){
      spectrum_ks_diag_meson(qp0, naik_index0, pair);
    } else {
      spectrum_ks_offdiag_meson(qp0, naik_index0, qp1, naik_index1, pair);
    }      
    dtime = start_timing();
    spectrum_ks_print_diag(pair);
    print_timing(dtime, "printing correlator");

    if(param.saveflag_m[pair] == SAVE_SQLITE)
      #ifdef HAVE_SQLITE
      {
	const char *dbName = &param.savefile_m[pair][0];
	sqlite3 *db = NULL;
	// timestamp: seconds since the epoch
	time_t epoch_secs = time(NULL);
	if(this_node == 0)
	  {
	    // connect to db
	    ret = db_connect(&db,dbName);
	    // init tables if they do not exist
	    ret = db_init_tables(db);
	    // begin transaction
	    ret = db_begin_transaction(db);
	  }
	// insert correlators
	ret = sql_spectrum_ks(db,pair,epoch_secs); // must be called from all nodes to do tie-ups 
	if(this_node == 0)
	  {
	    // end transaction
	    if(ret == SQLITE_OK){
	      ret = db_commit_transaction(db);
	    } else {
	      ret = db_rollback_transaction(db);
	    }
	    // disconnect from db
	    ret = db_disconnect(db);
	  }
      }
      #else
      {
	node0_printf("WARNING: enable sqlite3 in build and recompile to use option 'save_corr_sqlite'\n");
      }
      #endif
  }
	    
  spectrum_ks_cleanup(pair);
}
/*--------------------------------------------------------------------*/
void spectrum_ks_baryon(ks_prop_field *qp0, ks_prop_field *qp1, ks_prop_field *qp2, int triplet)
{
  int ret;
  int do_baryon = param.do_baryon_spect[triplet];

  spectrum_ks_baryon_init(triplet);

  if(do_baryon) {
    spectrum_ks_baryon_nd(param.num_corr_b[triplet], 
			  param.baryon_type_snk[triplet],
			  qp0, qp1, qp2, triplet);

    double dtime = start_timing();
    spectrum_ks_print_baryon(triplet);
    print_timing(dtime, "printing correlator");
  }

  if(do_baryon && param.saveflag_b[triplet]==SAVE_SQLITE)
    #ifdef HAVE_SQLITE
    {
    const char *dbName = &param.savefile_b[triplet][0];
    sqlite3 *db = NULL;
    // timestamp: seconds since the epoch
    time_t epoch_secs = time(NULL);
    if(this_node == 0)
      {
	// connect to db
	ret = db_connect(&db,dbName);
	// init tables if they do not exist
	ret = db_init_tables(db);
	// begin transaction
	ret = db_begin_transaction(db);
      }
    // insert correlators
    ret= sql_spectrum_ks_baryon(db,triplet,epoch_secs); // must be called from all nodes to do tie-ups
    if(this_node == 0)
      {
	// end transaction
	if(ret == SQLITE_OK){
	  ret = db_commit_transaction(db);
	} else {
	  ret = db_rollback_transaction(db);
	}
	// disconnect from db
	ret = db_disconnect(db);
      }
    }
    #else
    {
      node0_printf("WARNING: enable sqlite3 in build and recompile to use option 'save_corr_sqlite'\n");
    }
    #endif

  spectrum_ks_baryon_cleanup(triplet);
}
