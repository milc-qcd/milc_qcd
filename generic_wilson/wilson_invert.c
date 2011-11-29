/********************** wilson_invert.c *********************************/
/* MIMD version 7 */

/*  Generic wrapper for quark inverter
 *
 *  This wrapper works for WILSON or CLOVER fermions.
 *  Four variants appear here, depending on whether
 *  (1) data is in the site structure or in flat arrays (_site, _field)
 *  (2) a source should be built or it is preloaded (_wqs or not)
 *  The site-based routines are provided for backward compatibility.
 *  We are phasing out the site structure.
 */

/* Modifications

   4/30/00 Set min = 0 on restart CD
   4/26/98 Wrapper made generic after introducing structures for
           inversion parameters   CD
   11/14/97
   4/27/07 Rearrange to support field-based data.
*/

#include "generic_wilson_includes.h"

int wilson_invert_field( /* Return value is number of iterations taken */
    wilson_vector *src, /* type wilson_vector (where source has been created)*/
    wilson_vector *dest, /* type wilson_vector (answer and initial guess) */
    int (*invert_func_field)(wilson_vector *, wilson_vector *dest,
			     quark_invert_control *qic,
			     void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters;

  /* Do the inversion */
  tot_iters = invert_func_field(src,dest,qic,dmp);

  report_status(qic);
  
  return tot_iters;
} /* wilson_invert_field */

/* This variant builds the source, based on the wqs specification */

int wilson_invert_field_wqs( /* Return value is number of iterations taken */
    quark_source *wqs, /* source parameters */
    int (*source_func_field)(wilson_vector *src, 
			      quark_source *wqs),  /* source function */
    wilson_vector *dest,  /* type wilson_vector (answer and initial guess) */
    int (*invert_func_field)(wilson_vector *src, wilson_vector *dest,
			     quark_invert_control *qic, void *dcp),
    quark_invert_control *qic, /* inverter control */
    void *dmp              /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters;
  wilson_vector *src;

  /* Make the source */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  src = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  if(src == NULL){
    printf("wilson_invert_field_wqs(%d): Can't allocate src\n",this_node);
    terminate(1);
  }
  if(source_func_field(src,wqs)){
    printf("wilson_invert_field_wqs(%d): error getting source\n",this_node);
    terminate(1);
  };

  /* Do the inversion */
  tot_iters = invert_func_field(src,dest,qic,dmp);

  report_status(qic);

  free(src);
  return tot_iters;
} /* wilson_invert_field_wqs */


int wilson_invert_site( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source has been created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int (*invert_func_site)(field_offset src, field_offset dest,
			quark_invert_control *qic,
			void *dmp),
    quark_invert_control *qic, /* inverter control */
    void *dmp                 /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters;

  /* Do the inversion */
  tot_iters = invert_func_site(src,dest,qic,dmp);
  
  report_status(qic);
  
  return tot_iters;
} /* wilson_invert_site */

/* This variant builds the source, based on the wqs specification */

int wilson_invert_site_wqs( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    void (*source_func_site)(field_offset src, 
			quark_source *wqs),  /* source function */
    quark_source *wqs, /* source parameters */
    int (*invert_func_site)(field_offset src, field_offset dest,
			quark_invert_control *qic,void *dcp),
    quark_invert_control *qic, /* inverter control */
    void *dmp              /* Passthrough Dirac matrix parameters */
    )
{
  int tot_iters;

  /* Make the source */
  source_func_site(src,wqs);

  /* Do the inversion */
  tot_iters = invert_func_site(src,dest,qic,dmp);

  report_status(qic);
  
  return tot_iters;
} /* wilson_invert_site_wqs */


