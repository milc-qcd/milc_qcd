/* ************************************************************	*/
/*								*/
/*			    GAUGE_INFO.C			*/
/* For pure SU(3)xU(1) gauge routines, does not contain much.	*/
/* (for details see gauge_info.c file in any of milcv7		*/
/*  application directory)					*/
/*								*/
/* Last updated on 05.03.07					*/
/*								*/
/* ************************************************************	*/

#include "include_u1g.h"

void write_appl_u1gauge_info(FILE *fp, gauge_file *gf)
{

  /* u(1) part of gauge information */
  write_gauge_info_item(fp,"action.description","%s",
			"\"Pure U(1) gauge\"",0,0);
  write_gauge_info_item(fp,"gauge.description","%s",
			"\"Coulomb gauge fixed action\"",0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&echarge,0,0);

} /* end of write_appl_u1gauge_info() */

/* of no use here, only for compilation */
void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{

   write_gauge_info_item(fp,"action.description","%s",
                        "\"Pure U(1), no SU(3) gauge field\"",0,0);

} /* end of write_appl_gauge_info() */

/* ************************************************************	*/

