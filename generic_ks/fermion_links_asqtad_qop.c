/****************** fermion_links_asqtad_qop.c ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for SciDAC Level 3 QOP link smearing */

#include "generic_ks_includes.h"	/* definitions files and prototypes */

/* Load QOP_FermionLinksAsqtad object from MILC gauge field */
void create_qop_asqtad_fermion_links_from_gauge_field(
		      QOP_FermionLinksAsqtad** qop_links)
{
  su3_matrix *raw_gauge_links;
  QOP_info_t info;
  QOP_asqtad_coeffs_t coeffs;
  QOP_GaugeField *links;

  raw_gauge_links = create_raw_G_from_site_links(EVENANDODD);
  links = QOP_create_G_from_raw((Real *)(raw_gauge_links),QOP_EVENODD);
  destroy_raw_G(raw_gauge_links);   raw_gauge_links = NULL;

  load_qop_asqtad_coeffs(&coeffs, 1.);
  *qop_links = QOP_asqtad_create_L_from_G(&info, &coeffs, links);
}

/* load_fermion_links_fn wrapper for QOP */
void load_fermion_links_fn() {

  QOP_FermionLinksAsqtad *qop_links;
  su3_matrix **fatlinks;
  su3_matrix **longlinks;

  if( phases_in != 1){
    node0_printf("load_fermion_links_fn: BOTCH: needs phases in\n");
    terminate(1);
  }

  fatlinks = create_raw_G();
  if(fatlinks == NULL)terminate(1);

  longlinks = create_raw_G();
  if(longlinks == NULL)terminate(1);

  create_qop_asqtad_fermion_links_from_gauge_field(&qop_links);
  QOP_asqtad_extract_L_to_raw((Real *)fatlinks, (Real *)longlinks, qop_links);
  unload_raw_G_to_field(t_longlinks, longlinks, EVENANDODD);
  unload_raw_G_to_field(t_fatlinks,  fatlinks,  EVENANDODD);
  destroy_raw_G(longlinks);
  destroy_raw_G(fatlinks);
  QOP_asqtad_destroy_L(qop_links);
}

/* Load QOP_FermionLinksAsqtad object from MILC fat and long links */
void create_qop_asqtad_fermion_links( QOP_FermionLinksAsqtad** qop_links )
{
  su3_matrix **raw_fat_links, **raw_long_links;

  /* Create fat and long links if necessary */

  if( valid_links_fn  != 1 ){
    /* Create qop links directoyr from gauge field */
    create_qop_asqtad_fermion_links_from_gauge_field(qop_links);
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
    free_links_fn();
    valid_links_fn = 0;
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
