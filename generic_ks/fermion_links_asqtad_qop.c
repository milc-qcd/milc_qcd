/****************** fermion_links_asqtad_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */

#include "generic_ks_includes.h"	/* definitions files and prototypes */

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

/* Load QOP_FermionLinksAsqtad object from MILC gauge field */
void create_qop_asqtad_fermion_links_from_gauge_field(
		      QOP_FermionLinksAsqtad** qop_links,
		      Real *act_path_coeff)
{
  su3_matrix **raw_gauge_links;
  QOP_info_t info;
  QOP_asqtad_coeffs_t coeffs;
  QOP_GaugeField *links;

  raw_gauge_links = create_raw_G_from_site_links(EVENANDODD);
  links = QOP_create_G_from_raw((Real **)(raw_gauge_links),QOP_EVENODD);
  destroy_raw_G(raw_gauge_links);   raw_gauge_links = NULL;

  load_qop_asqtad_coeffs(&coeffs, 0.5, act_path_coeff);
  *qop_links = QOP_asqtad_create_L_from_G(&info, &coeffs, links);
}


void load_asqtad_links(int both, su3_matrix **t_fl, su3_matrix **t_ll,
		       Real *act_path_coeff) {

  QOP_FermionLinksAsqtad *qop_links;
  su3_matrix **fatlinks;
  su3_matrix **longlinks;
  char myname[] = "load_asqtad_links";

  if( phases_in != 1){
    node0_printf("load_fermion_links_fn: BOTCH: needs phases in\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("%s(%d): Error initializing QOP\n",myname,this_node);
    terminate(1);
  }

  fatlinks = create_raw_G();
  if(fatlinks == NULL)terminate(1);

  longlinks = create_raw_G();
  if(longlinks == NULL)terminate(1);

  create_qop_asqtad_fermion_links_from_gauge_field(&qop_links,act_path_coeff);
  QOP_asqtad_extract_L_to_raw((Real **)fatlinks, (Real **)longlinks, 
			      qop_links, QOP_EVENODD);

  /* Allocate space for t_fl if NULL */
  if(*t_fl == NULL){
    *t_fl = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_fl==NULL){
      printf("%s(%d): no room for t_fl\n",myname,this_node);
      terminate(1);
    }
  }
  
  /* Allocate space for t_ll if NULL and we are doing both fat and long */
  if(*t_ll == NULL && both){
    *t_ll = (su3_matrix *)special_alloc(sites_on_node*4*sizeof(su3_matrix));
    if(*t_ll==NULL){
      printf("%s(%d): no room for t_ll\n",myname,this_node);
      terminate(1);
    }
  }
  
  unload_raw_G_to_field(*t_fl,  fatlinks,  EVENANDODD);
  if(both)unload_raw_G_to_field(*t_ll, longlinks, EVENANDODD);
  destroy_raw_G(fatlinks);
  destroy_raw_G(longlinks);
  QOP_asqtad_destroy_L(qop_links);

  valid_fn_links = 1;
}

/* Wrappers for MILC call to QOP */
void load_fn_links(){
  load_asqtad_links(1, &t_fatlink, &t_longlink, get_quark_path_coeff());

#ifdef DBLSTORE_FN
  load_fatbacklinks(&t_fatbacklink, t_fatlink);
  load_longbacklinks(&t_longbacklink, t_longlink);
#endif
}

#ifdef DM_DU0
/* Wrappers for MILC call to QOP */
void load_fn_links_dmdu0(){
  su3_matrix *null = NULL;
  load_asqtad_links(0, &t_dfatlink_du0, &null, get_quark_path_coeff_dmdu0());
}
#endif

/* Load QOP_FermionLinksAsqtad object from MILC fat and long links */
void create_qop_asqtad_fermion_links( QOP_FermionLinksAsqtad** qop_links )
{
  su3_matrix **raw_fat_links, **raw_long_links;

  /* Create fat and long links if necessary */

  if( valid_fn_links != 1 ){
    /* Create qop links directly from gauge field */
    create_qop_asqtad_fermion_links_from_gauge_field(qop_links,
						     get_quark_path_coeff());
  }
  else{
    /* Create qop links from existing t_fatlink and t_longlink */
    /* Map fat and long links to raw format */

    raw_fat_links  = create_raw_G_from_field_links(t_fatlink,EVENANDODD);
    if(raw_fat_links == NULL)terminate(1);
    raw_long_links = create_raw_G_from_field_links(t_longlink,EVENANDODD);
    if(raw_long_links == NULL)terminate(1);
    
#if 0
    // Release for memory savings.  Links are recomputed later.
    free_fn_links();
    valid_fn_links = 0;
#endif

    /* Map raw to QOP format */
    
    *qop_links = QOP_asqtad_create_L_from_raw((Real **)raw_fat_links, 
					      (Real **)raw_long_links,
					      QOP_EVENODD);
    destroy_raw_G(raw_fat_links);   raw_fat_links = NULL;
    destroy_raw_G(raw_long_links);  raw_long_links = NULL;
  }
  return;
}
