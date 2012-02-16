/******* dslash_fn_qop.c - dslash QOP wrapper for SU3/fermions ****/
/* MIMD version 7 */

/* TEMPORARY ROUTINE UNTIL WE ELIMINATE THE REST OF THE SITE FIELDS */

/* THIS IS NOT MEANT TO SUPPORT congrad.  IT IS ONLY FOR 
   OCCASIONAL USE, SINCE IT INVOLVES REMAPPING MILC FIELDS,
   WHICH IS NOT OPTIMUM. */

/* This is the MILC wrapper for the SciDAC Level 3 QOP dslash */

#include "generic_ks_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"

static int 
opposite_parity(int parity){

  int otherparity;
  if(parity == EVEN)
    otherparity = ODD;
  else if(parity == ODD)
    otherparity = EVEN;
  else
    otherparity = parity;

  return otherparity;
}

static
void qop_to_milc_dslash_normalization(su3_vector *sol, int parity)
{
  int i;

  FORSOMEFIELDPARITY(i,parity){
    scalar_mult_su3_vector(sol+i, 2.0, sol+i);
  }
}


/* "parity" is the destination parity */

static void 
dslash_fn_field_F( su3_vector *src, su3_vector *dest, 
		   int parity, fn_links_qop_t *fn){

  QOP_info_t info = {0., 0., 0, 0, 0};
  int otherparity = opposite_parity(parity);
  QOP_F3_FermionLinksAsqtad *al_F = get_F_asqtad_links(fn);
  QOP_F3_ColorVector *in = create_F_V_from_field( src, otherparity);
  QOP_F3_ColorVector *out = create_F_V_from_field( dest, parity);
  QOP_evenodd_t eo_out, eo_in;

  eo_out = milc2qop_parity(parity);
  eo_in = milc2qop_parity(otherparity);

  /* We want just dslash so we use zero mass */
  QOP_F3_asqtad_dslash(&info, al_F, 0.0, out, in, eo_out, eo_in);

  unload_F_V_to_field( dest, out, parity);

  /* Fix the normalization */
  qop_to_milc_dslash_normalization(dest, parity);

  QOP_F3_destroy_V(out);
  QOP_F3_destroy_V(in);
}  
  
/* "parity" is the destination parity */
  
static void 
dslash_fn_field_D( su3_vector *src, su3_vector *dest, 
		   int parity, fn_links_qop_t *fn){

  QOP_info_t info = {0., 0., 0, 0, 0};
  int otherparity = opposite_parity(parity);
  QOP_D3_FermionLinksAsqtad *al_D = get_D_asqtad_links(fn);
  QOP_D3_ColorVector *in = create_D_V_from_field( src, otherparity);
  QOP_D3_ColorVector *out = create_D_V_from_field( dest, parity);
  QOP_evenodd_t eo_out, eo_in;

  eo_out = milc2qop_parity(parity);
  eo_in = milc2qop_parity(otherparity);

  /* We want just dslash so we use zero mass */
  QOP_D3_asqtad_dslash(&info, al_D, 0.0, out, in, eo_out, eo_in);

  unload_D_V_to_field( dest, out, parity);

  /* Fix the normalization */
  qop_to_milc_dslash_normalization(dest, parity);

  QOP_D3_destroy_V(out);
  QOP_D3_destroy_V(in);
}  
  
  
void 
dslash_fn_field( su3_vector *src, su3_vector *dest, 
		 int parity, fn_links_qop_t *fn){

  if(get_D_asqtad_links(fn) != NULL)
    dslash_fn_field_D( src, dest, parity, fn );
  else
    dslash_fn_field_F( src, dest, parity, fn );
}


void 
dslash_fn_field_special( su3_vector *src, su3_vector *dest, 
			 int parity, msg_tag **tag, int dslash_start, 
			 fn_links_qop_t *fn){
  dslash_fn_field(src, dest, parity, fn);
}

void
cleanup_gathers(msg_tag *tags1[], msg_tag *tags2[]){
}

void
cleanup_dslash_temps(void){
}

void 
dslash_fn_site( field_offset src, field_offset dest, int parity,
		fn_links_qop_t *fn ){

  su3_vector *t_src, *t_dest;

  t_src = create_v_field_from_site_member(src);
  t_dest = create_v_field_from_site_member(dest);

  dslash_fn_field( t_src, t_dest, parity, fn );

  copy_site_member_from_v_field(dest, t_dest);

  destroy_v_field(t_dest);
  destroy_v_field(t_src);
  
}

/* Apply the symmetric shift operator in direction "dir" *
 * Fat and long links are used instead of unfattened links          */

static void 
dslash_fn_dir_F(su3_vector *src, su3_vector *dest, int parity,
		imp_ferm_links_t *fn, int dir, int fb, 
		Real wtfat, Real wtlong)
{
  QOP_info_t info = {0., 0., 0, 0, 0};
  QOP_evenodd_t eo_out = milc2qop_parity(parity);
  int otherparity = opposite_parity(parity);
  QOP_F3_FermionLinksAsqtad *al_F = get_F_asqtad_links(fn);
  QOP_F3_ColorVector *in = create_F_V_from_field( src, otherparity);
  QOP_F3_ColorVector *out = create_F_V_from_field( dest, parity);
  /* Factor of 2 to compensate for QOP normalization */
  double mywtfat = 2.*wtfat;
  double mywtlong = 2.*wtlong;

  QOP_F3_asqtad_dslash_dir(&info, al_F, dir, fb, mywtfat, mywtlong, 
			   out, in, eo_out);
  unload_F_V_to_field( dest, out, parity );

  QOP_F3_destroy_V(out);
  QOP_F3_destroy_V(in);

}  
  
static void 
dslash_fn_dir_D(su3_vector *src, su3_vector *dest, int parity,
		imp_ferm_links_t *fn, int dir, int fb, 
		Real wtfat, Real wtlong)
{
  QOP_info_t info = {0., 0., 0, 0, 0};
  QOP_evenodd_t eo_out = milc2qop_parity(parity);
  int otherparity = opposite_parity(parity);
  QOP_D3_FermionLinksAsqtad *al_D = get_D_asqtad_links(fn);
  QOP_D3_ColorVector *in = create_D_V_from_field( src, otherparity);
  QOP_D3_ColorVector *out = create_D_V_from_field( dest, parity);
  /* Factor of 2 to compensate for QOP normalization */
  double mywtfat = 2.*wtfat;
  double mywtlong = 2.*wtlong;

  QOP_D3_asqtad_dslash_dir(&info, al_D, dir, fb, mywtfat, mywtlong, 
			   out, in, eo_out);
  unload_D_V_to_field( dest, out, parity );

  QOP_D3_destroy_V(out);
  QOP_D3_destroy_V(in);

}  
  
  
void 
dslash_fn_dir(su3_vector *src, su3_vector *dest, int parity,
	      imp_ferm_links_t *fn, int dir, int fb, 
	      Real wtfat, Real wtlong)
{

  if(get_D_asqtad_links(fn) != NULL)
    dslash_fn_dir_D( src, dest, parity, fn, dir, fb, wtfat, wtlong );
  else
    dslash_fn_dir_F( src, dest, parity, fn, dir, fb, wtfat, wtlong );
}

