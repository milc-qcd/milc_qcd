/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For ks_eigen */

/* Application-dependent routine for writing gauge info file */
/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_lat.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

#include "ks_eig_includes.h"
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  printf("write_appl_gauge_info: HEY, I'M NOT SUPPOSED TO BE USED HERE!\n");
}

/* Print string formatting with digits based on intended precision */
void print_prec(char string[], size_t n, Real value, int prec){
  if(prec == 1){
    snprintf(string,n,"%.6e",value);  /* single precision */
  }
  else
    snprintf(string,n,"%.15e",value); /* double precision */
}

static char MILC_info[] = "";
char *create_MILC_info(){
  return MILC_info;
}

void destroy_MILC_info(char *info){
}

#ifdef HAVE_QIO
static QIO_String *xml_record;

/* Follow USQCD style for record XML */

char *create_QCDML(){
  QIO_USQCDLatticeInfo *record_info;
  char *milc_info;
  char plaqstring[32];
  char linktrstring[32];
  
  /* Build the components of the USQCD info string */
  print_prec(plaqstring, 32, (g_ssplaq+g_stplaq)/6., 1);
  print_prec(linktrstring, 32, linktrsum.real/3., 1);
  milc_info = create_MILC_info();  /* The MILC info data as a string */

  /* Stuff the data structure */
  record_info = QIO_create_usqcd_lattice_info(plaqstring, linktrstring, 
					      milc_info);

  destroy_MILC_info(milc_info);

  /* Convert the data structure to a QIO character string */
  xml_record = QIO_string_create();
  QIO_encode_usqcd_lattice_info(xml_record, record_info);
  QIO_destroy_usqcd_lattice_info(record_info);

  /* Return a pointer to the actual character string */
  return QIO_string_ptr(xml_record);
}

void free_QCDML(char *info){
  QIO_string_destroy(xml_record);
}

#endif
