/*********************** io_scidac_ks_eigen.c *************************/
/* MIMD version 7 */

/* 1/16  C. DeTar created */

/* For QIO-formatted eigenvalue and eigenvector files */

#include "generic_ks_includes.h"
#ifndef HAVE_QIO
# error REQUIRES QIO
#else
#include <qio.h>
#endif
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_ks_eigen.h"
#include <string.h>
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"

#define FILEINFOSTRING_MAX 512
#define RECINFOSTRING_MAX  256

static char *
create_file_xml(int Nvecs, int packed){
  char begin_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info><title>KS eigenvalues and vectors</title>";
  char end_xml[] = "</info>";
  char *xml = (char *)malloc(FILEINFOSTRING_MAX);
  size_t bytes = 0;
  size_t max = FILEINFOSTRING_MAX;
  
  /* Create the file XML */
  snprintf(xml+bytes, max-bytes, "%s", begin_xml);
  bytes = strlen(xml);

  snprintf(xml+bytes, max-bytes, "<Nvecs>%d</Nvecs>", Nvecs);
  bytes = strlen(xml);
  
  if(packed){
    snprintf(xml+bytes, max-bytes, "<format>Packed</format>");
    bytes = strlen(xml);
  }

  snprintf(xml+bytes, max-bytes, "%s", end_xml);
  bytes = strlen(xml);

  return xml;
}

/* Extract the number of vectors from the file XML */

/* A real XML parser would be nice! */

static void
parse_file_xml_Nvec(int *Nvecs, char *xml){

  char begtag[] = "<Nvecs>";
  char endtag[] = "</Nvecs>";
  char *pb, *pe;
  int status;

  /* Find the end tag */
  pe = strstr(xml, endtag);
  if(pe == NULL){
    *Nvecs = 0;
    return;
  }

  /* Truncate the string at the end tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtag);
  
  if(pb == NULL){
    *Nvecs = 0;
    return;
  }

  /* Read the tagged integer */
  status = sscanf(pb + strlen(begtag), "%d", Nvecs);

  /* Restore the xml */
  *pe = '<';
  
  if(status != 1)
    *Nvecs = 0;
}

/* Check whether the eigenvectors are packed */

/* A real XML parser would be nice! */

static void
parse_file_xml_packed(int *packed, char *xml){

  char begtag[] = "<format>";
  char endtag[] = "</format>";
  char *pb, *pe;
  int status;

  /* Find the end tag */
  pe = strstr(xml, endtag);
  if(pe == NULL){
    *packed = 0;
    return;
  }

  /* Truncate the string at the end tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtag);
  
  if(pb == NULL){
    *packed = 0;
    return;
  }

  /* Read the tagged integer */
  status = strcmp(pb + strlen(begtag), "Packed");

  /* Restore the xml */
  *pe = '<';

  if(status==0) *packed = 1;
  else *packed = 0;
}

/* Create the record xml, encoding the eigenvalue */

static char *
create_record_xml(double eigVal, double resid){
  char begin_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  char end_xml[] = "</info>";
  char *xml = (char *)malloc(FILEINFOSTRING_MAX);
  size_t bytes = 0;
  size_t max = FILEINFOSTRING_MAX;
  
  /* Create the file XML */
  snprintf(xml+bytes, max-bytes, "%s", begin_xml);
  bytes = strlen(xml);

  snprintf(xml+bytes, max-bytes, "<eigVal>%.18e</eigVal>", eigVal);
  bytes = strlen(xml);

  snprintf(xml+bytes, max-bytes, "<resid>%e</resid>", resid);
  bytes = strlen(xml);

  snprintf(xml+bytes, max-bytes, "%s", end_xml);
  bytes = strlen(xml);

  return xml;

}

/* Extract the eigenvalue from the file XML */

static void
parse_record_xml(double *eigVal, char *xml){

  char begtag[] = "<eigVal>";
  char endtag[] = "</eigVal>";
  char *pb, *pe;
  int status;

  /* Find the end tag */
  pe = strstr(xml, endtag);
  if(pe == NULL){
    *eigVal = 0;
    return;
  }

  /* Truncate the string at the ending tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtag);
  
  if(pb == NULL){
    *eigVal = 0;
    return;
  }

  status = sscanf(pb + strlen(begtag), "%lf", eigVal);
  
  if(status != 1)
    *eigVal = 0;

}


/* Open a file for writing eigenvectors */

QIO_Writer *
open_ks_eigen_outfile(const char *filename, int Nvecs, int volfmt, int serpar, int packed){
  char *xml;

  QIO_String *filexml = QIO_string_create();
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Create the file XML */
  xml = create_file_xml(Nvecs, packed);
  QIO_string_set(filexml, xml);

  /* Open the file for output */
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO, NULL,
			       &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  free(xml);

  return outfile;
}

/* Write an eigenvector and its eigenvalue */

int
write_ks_eigenvector(QIO_Writer *outfile, int packed, su3_vector *eigVec, double eigVal, 
		     double resid){
  int status;
  char *xml;
  QIO_String *recxml = QIO_string_create();

  xml = create_record_xml(eigVal, resid);
  QIO_string_set(recxml, xml);

  if(packed){
    pack_field(eigVec, sizeof(su3_vector));
    if(MILC_PRECISION == 1)
      status = write_F3_V_from_half_field(outfile, recxml, eigVec, 1);
    else
      status = write_D3_V_from_half_field(outfile, recxml, eigVec, 1);
  } else {
    if(MILC_PRECISION == 1)
      status = write_F3_V_from_field(outfile, recxml, eigVec, 1);
    else
      status = write_D3_V_from_field(outfile, recxml, eigVec, 1);
  }

  QIO_string_destroy(recxml);
  free(xml);

  return status;
}

/* Close the output eigenvector file */

void
close_ks_eigen_outfile(QIO_Writer *outfile){
  QIO_close_write(outfile);
}

/* Open the eigenvector file for reading */

QIO_Reader *
open_ks_eigen_infile(const char *filename, int *Nvecs, int *packed, int *file_type, int serpar){
  char myname[] = "open_ks_eigen_infile";
  char *xml;

  QIO_String *filexml = QIO_string_create();
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  int Nvecs_test;

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O system */
  build_qio_filesystem(&fs);

  /* Open the file for input */
  node0_printf("%s: Opening %s for reading\n", myname, filename);
  infile = open_scidac_input_xml(filename, &layout, &fs, serpar, filexml);
  if(infile == NULL) return infile;

  /* Interpret the file XML */

  xml = QIO_string_ptr(filexml);
  parse_file_xml_Nvec(&Nvecs_test, xml);
  parse_file_xml_packed(packed, xml);

  if(Nvecs_test > 0){
    /* Warn if the number of eigenvectors doesn't match expectations,
       and if the number of eigenvectors is less than expected, reduce
       the expected number */
    *file_type = 0;
    if(Nvecs_test != *Nvecs){
      node0_printf("WARNING: Called for %d vectors, but found %d\n", 
		   *Nvecs, Nvecs_test);
      if(*Nvecs > Nvecs_test){
	node0_printf("WARNING: Resetting Nvecs = %d\n", Nvecs_test);
	*Nvecs = Nvecs_test;
      }
    }
  } else {
    node0_printf("%s: WARNING: Nvecs tag is either missing or 0\n", myname);
    node0_printf("%s: WARNING: Treating it as a QUDA file\n", myname);
    *file_type = 1;
  }

  QIO_string_destroy(filexml);

  return infile;
}

/* Read an eigenvector and its eigenvalue */

int
read_ks_eigenvector(QIO_Reader *infile, int packed, su3_vector *eigVec, double *eigVal){
  int status;
  char *xml;
  QIO_String *recxml = QIO_string_create();
  QIO_RecordInfo recinfo;

  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS){
    QIO_string_destroy(recxml);
    return status;
  }

  int typesize = QIO_get_typesize(&recinfo);
  if(typesize == 24)
    status = read_F3_V_to_field(infile, recxml, eigVec, 1);
  else if (typesize == 48 )
    status = read_D3_V_to_field(infile, recxml, eigVec, 1);
  else
    {
      node0_printf("read_ks_eigenvector: Bad typesize %d\n",typesize);
      terminate(1);
    }

  if(packed){
    double dt = -dclock();
    unpack_field(eigVec, sizeof(su3_vector));
    dt += dclock();
    node0_printf("%s unpack time %0.2f\n",__func__,dt);
  }

  if(status != QIO_EOF){
    xml = QIO_string_ptr(recxml);
    parse_record_xml(eigVal, xml);
  }

  QIO_string_destroy(recxml);

  return status;
}

/* Read a set of eigenvectors from a QUDA-formated color-spin-field file */

int
read_quda_ks_eigenvectors(QIO_Reader *infile, su3_vector *eigVec[], double *eigVal, int *Nvecs,
			  int parity){
  char myname[] = "read_quda_ks_eigenvectors";
  int status;
  char *xml;
  QIO_String *recxml = QIO_string_create();
  QIO_RecordInfo recinfo;

  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS){
    QIO_string_destroy(recxml);
    return status;
  }

  const char *datatype = QIO_get_datatype(&recinfo);
  if(strcmp("QUDA_DNs1Nc3_ColorSpinorField", datatype) != 0){
    node0_printf("%s: WARNING: Unexpected datatype.  Found %s\n",
		 myname, datatype);
  }
  
  int typesize = QIO_get_typesize(&recinfo);
  int datacount = QIO_get_datacount(&recinfo);

  if(datacount < *Nvecs){
    node0_printf("%s WARNING: Requested %d eigenvectors but the file has %d.\n",
		 myname, *Nvecs, datacount);
    node0_printf("%s WARNING: Reducing the request.\n", myname);
    *Nvecs = datacount;
  }

  //  su3_vector *eigVecs = (su3_vector *)malloc((*Nvecs)*sizeof(su3_vector)*sites_on_node);
  su3_vector *eigVecs = (su3_vector *)malloc(datacount*sizeof(su3_vector)*sites_on_node);
  if(eigVecs == NULL){
    node0_printf("%s FATAL: No room for a temporary array for %d eigenvectors\n",
		 myname, datacount);
    terminate(1);
  }
  
  node0_printf("Reading %d eigenvectors\n",datacount);
  if(typesize == 24)
    status = read_F3_V_to_field(infile, recxml, eigVecs, datacount);
  else if (typesize == 48 )
    status = read_D3_V_to_field(infile, recxml, eigVecs, datacount);
  else
    {
      node0_printf("read_ks_eigenvector: Bad typesize %d\n",typesize);
      terminate(1);
    }

  if(status != QIO_EOF){
    xml = QIO_string_ptr(recxml);
    parse_record_xml(eigVal, xml);
  }

  QIO_string_destroy(recxml);

  /* Map eigenvectors to our dynamic array */

  int i;
  FORALLFIELDSITES_OMP(i,){
    for(int j = 0; j < *Nvecs; j++){
      eigVec[j][i] = eigVecs[*Nvecs*i+j];
    }
  } END_LOOP_OMP;

  free(eigVecs);
  
  return status;
}

/* Close the input eigenvector file */

void
close_ks_eigen_infile(QIO_Reader *infile){
  QIO_close_read(infile);
}

