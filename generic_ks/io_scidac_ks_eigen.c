/*********************** io_scidac_ks_eigen.c *************************/
/* MIMD version 7 */

/* 1/16  C. DeTar created */

/* For QIO-formatted eigenvalue and eigenvector files */

#include "generic_ks_includes.h"
#ifndef HAVE_QIO
# error REQUIRES QIO
#else
#include <lime.h>
#endif
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_ks_eigen.h"
#include <string.h>
#define LOOPEND
#include "../include/loopend.h"
#include "../include/openmp_defs.h"
#include <zlib.h>

off_t ftello(FILE *stream);
int fseeko(FILE *stream, off_t offset, int whence);
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
parse_xml_eigVal(Real *eigVal, char *xml){

  char begtagMILC[] = "<eigVal>";
  char endtagMILC[] = "</eigVal>";
  char begtagGRID[] = "<eval>";
  char endtagGRID[] = "</eval>";
  char *pb, *pe, *p;
  int status;

  /* Find the end tag */
  pe = strstr(xml, endtagMILC);
  if(pe == NULL){
    pe = strstr(xml, endtagGRID);
    if(pe == NULL){
      *eigVal = 0;
      return;
    }
  }

  /* Truncate the string at the ending tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtagMILC);
  if(pb != NULL){
    p = pb + strlen(begtagMILC);
  } else {
    pb = strstr(xml, begtagGRID);
    if(pb != NULL){
      p = pb + strlen(begtagGRID);
    }
  }

  if(pb == NULL){
    *eigVal = 0;
    return;
  }

#if MILC_PRECISION == 1
  status = sscanf(p, "%f", eigVal);
#else
  status = sscanf(p, "%lf", eigVal);
#endif
  
  if(status != 1) *eigVal = 0;

}

/* Extract the typesize  from the file XML */

static void
parse_xml_typesize(int *typesize, char *xml){

  char begtagMILC[] = "<typesize>";
  char endtagMILC[] = "</typesize>";
  char *pb, *pe, *p;
  int status;

  /* Find the end tag */
  pe = strstr(xml, endtagMILC);
  if(pe == NULL){
    *typesize = 0;
    return;
  }

  /* Truncate the string at the ending tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtagMILC);
  if(pb != NULL){
    p = pb + strlen(begtagMILC);
  }

  if(pb == NULL){
    *typesize = 0;
    return;
  }

  status = sscanf(p, "%d", typesize);
  if(status != 1) *typesize = 0;
}


/* Extract the checksums from the file XML */

static int
parse_xml_checksums(uint32_t *csuma, uint32_t *csumb, char *xml){

  char begtaga[] = "<suma>";
  char endtaga[] = "</suma>";
  char begtagb[] = "<sumb>";
  char endtagb[] = "</sumb>";
  char *pb, *pe, *p, *pnew;
  int status;

  *csuma = 0;
  *csumb = 0;

  /* Find the end tag */
  pe = strstr(xml, endtaga);
  if(pe == NULL)  return 1;

  /* Truncate the string at the ending tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(xml, begtaga);
  if(pb == NULL) return 1;

  p = pb + strlen(begtaga);

  status = sscanf(p, "%x", csuma);
  if(status != 1) *csuma = 0;

  /* Find the end tag */
  pnew = pe+1;
  pe = strstr(pnew, endtagb);
  if(pe == NULL) return 1;

  /* Truncate the string at the ending tag */
  *pe = '\0';

  /* Find the beginning tag */
  pb = strstr(pnew, begtagb);
  if(pb == NULL) return 1;

  p = pb + strlen(begtagb);
  
  status = sscanf(p, "%x", csumb);
  if(status != 1) *csumb = 0;

  return 0;
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
write_ks_eigenvector(QIO_Writer *outfile, int packed, su3_vector *eigVec, Real eigVal, 
		     Real resid){
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
QIO_create_reader(const char *filename, 
		  QIO_Layout *layout, QIO_Iflag *iflag,
		  int (*io_node)(const int), int (*master_io_node)());


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
  QIO_verbose(QIO_VERB_OFF);
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

/*------------------------------------------------------------------*/
/* Convert rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], 
		       const size_t rank)
{
  int d;
  size_t r = rank;
  
  for(d = 0; d < dim; d++){
    coords[d] = r % size[d];
    r /= size[d];
  }
}

/*------------------------------------------------------------------*/
/* Convert coordinate to linear lexicographic PE rank (inverse of
   lex_coords) */

static size_t lex_rank(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
  }
  return rank;
}

/*------------------------------------------------------------------*/
/* Read an eigenvector and its eigenvalue */

int
read_ks_eigenvector(QIO_Reader *infile, int packed, su3_vector *eigVec, Real *eigVal, int parity){
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
    parse_xml_eigVal(eigVal, xml);
  }

  QIO_string_destroy(recxml);

  return status;
}

/*---------------------------------------------------------------*/

/* The Grid eigenvector file format is as follows...

    Message:        1
    Record:         1
    Type:           scidac-file-xml
    
    Message:        2
    Record:         1
    Type:           grid-format
    
    Message:        2
    Record:         2
    Type:           scidac-record-xml
    Data Length:    104
    Padding Length: 0
    MB flag:        0
    ME flag:        0
    Data:           "<?xml version="1.0"?>
    <VecRecord>
      <index>0</index>
      <eval>2.9215128862583114e-10</eval>
    </VecRecord>
    " 
    
    Message:        2
    Record:         3
    Type:           scidac-private-record-xml
    Data Length:    313
    Padding Length: 7
    MB flag:        0
    ME flag:        0
    Data:           "<?xml version="1.0"?>
    <scidacRecord>
      <version>1</version>
      <date>Mon Oct 24 01:11:28 2022 EDT</date>
      <recordtype>0</recordtype>
      <datatype>GRID_D_ColourVector3</datatype>
      <precision>D</precision>
      <colors>3</colors>
      <spins>1</spins>
      <typesize>48</typesize>
      <datacount>1</datacount>
    </scidacRecord>
    " 
    Message:        2
    Record:         4
    Type:           ildg-binary-data
    Data Length:    4076863488
    Padding Length: 0
    MB flag:        0
    ME flag:        0
    Data:           [Long record skipped]
    
    Message:        2
    Record:         5
    Type:           FieldNormMetaData
    Data Length:    82
    Padding Length: 6
    MB flag:        0
    ME flag:        0
    Data:           "<?xml version="1.0"?>
    <FieldNormMetaData>
      <norm2>1</norm2>
    </FieldNormMetaData>
    " 
    Message:        2
    Record:         6
    Type:           scidac-checksum
    Data Length:    128
    Padding Length: 0
    MB flag:        0
    ME flag:        1
    Data:           "<?xml version="1.0"?>
    <scidacChecksum>
      <version>1</version>
      <suma>47ed713a</suma>
      <sumb>fb019e7f</sumb>
    </scidacChecksum>
    " 
*/

/* We need the user and private record XMLs, the binary data, and the checksum.
   Node 0 reads the XML strings via LIME.
   The payload data are read using MPI parallel I/O
*/

/* Seek an XML LIME record and read the XML string */
char *
readLimeXml(LimeReader *LimeR, char *lime_type){
  char *xmlstring = NULL;
  
  while ( limeReaderNextRecord(LimeR) == LIME_SUCCESS ) {

    if ( strncmp( limeReaderType(LimeR), lime_type, strlen(lime_type) ) == 0 ) {
      uint64_t nbytes = limeReaderBytes(LimeR);
      xmlstring = (char *)malloc(sizeof(char)*(nbytes+1));
      memset(xmlstring, '\0', nbytes+1);
      limeReaderReadData((void *)xmlstring, &nbytes, LimeR);
      break;
    }

  }  

  return xmlstring;
}

/* Seek a LIME record, but don't read the payload.  Just return the
   byte count */
uint64_t
skipToLimeRecord(LimeReader *LimeR, char *lime_type){
  uint64_t nbytes = 0;
  
  while ( limeReaderNextRecord(LimeR) == LIME_SUCCESS ) {

    if ( strncmp( limeReaderType(LimeR), lime_type, strlen(lime_type) ) == 0 ) {
      nbytes = limeReaderBytes(LimeR);
    }
    break;
  }  

  return nbytes;
}

/* Read the Grid eigenector file */
int
read_grid_ks_eigenvector(char *eigfile, int *Nvecs, su3_vector *eigVec, Real *eigVal){
  int packed = 0, file_type = 0;
  char myname[] = "read_grid_ks_eigenvector";

  /* Open the file and read the private file XML, but ignore it */

  FILE *fpt = fopen(eigfile, "r");
  LimeReader *LimeR = limeCreateReader(fpt);
  if(LimeR == NULL)return 1;

  /* Skip to the scidac-record-xml and read the eigenvalue */
  if(this_node == 0){

    char lime_type[] = "scidac-record-xml";
    char *record_xml = readLimeXml(LimeR, lime_type);
    printf("record_xml = %s\n", record_xml);
    if(record_xml == NULL){
      printf("Couldn't find %s in %s\n", lime_type, eigfile);
      terminate(1);
    }
    parse_xml_eigVal(eigVal, record_xml);
    free(record_xml);
  }

  /* Broadcast the eigenvalue to the other nodes */
  broadcast_bytes((char *)eigVal, sizeof(*eigVal));
  node0_printf("%s(%d): eigVal = %g\n", myname, this_node, *eigVal);
  
  /* Read the scidac-private-record-xml and get the typesize */
  int typesize;
  if(this_node == 0){

    char lime_type[] = "scidac-private-record-xml";
    char *record_xml = readLimeXml(LimeR, lime_type);
    printf("record_xml = %s\n", record_xml);
    if(record_xml == NULL){
      printf("Couldn't find %s in %s\n", lime_type, eigfile);
      terminate(1);
    }
    parse_xml_typesize(&typesize, record_xml);
    free(record_xml);
  }

  /* Broadcast the typesize to the other nodes */
  broadcast_bytes((char *)&typesize, sizeof(typesize));
  // printf("%s(%d): typesize = %d\n", myname, this_node, typesize);
  
  /* Skip to the binary record data and read it */
  off_t offset = 0;
  if(this_node == 0){

    char lime_type[] = "ildg-binary-data";
    uint64_t nbytes = skipToLimeRecord(LimeR, lime_type);
    if(nbytes == 0){
      printf("Couldn't find %s in %s\n", lime_type, eigfile);
      terminate(1);
    }

    uint64_t bytes_wanted = typesize*volume/2;
    if(nbytes != bytes_wanted){
      printf("ERROR: LIME record byte count %lu is not equal to the wanted number %lu\n",
	     nbytes, bytes_wanted);
      terminate(1);
    }

    offset = ftello(fpt);
  }

  /* Broadcast the file offset to the other nodes */
  broadcast_bytes((char *)&offset, sizeof(offset));
  // printf("%s(%d): offset = %lu\n", myname, this_node, offset);

  /* Read the eigenvector. All nodes participate. */
  uint32_t data_suma, data_sumb;
  read_grid_eigenvector_data(eigVec, typesize, eigfile, offset, &data_suma, &data_sumb);


  /* Skip ahead and read the checksum record */
  uint32_t csuma, csumb;
  if(this_node == 0){

    /* Reset file pointer */
    /* Is this necessary? */
    fseeko(fpt, offset, SEEK_SET);

    char lime_type[] = "scidac-checksum";
    char *record_xml = readLimeXml(LimeR, lime_type);
    // printf("record_xml = %s\n", record_xml);
    if(record_xml == NULL){
      printf("Couldn't find %s in %s\n", lime_type, eigfile);
      terminate(1);
    }
    parse_xml_checksums(&csuma, &csumb, record_xml);
    free(record_xml);
  }

  int status = 0;
  if(this_node == 0){
    if(data_suma != csuma || data_sumb != csumb){
      printf("Checksum mismatch.  Expected %x %x. Found %x %x\n",
	     csuma, csumb, data_suma, data_sumb); fflush(stdout);
      status = 1;
    } else {
      printf("Checksums %x %x OK\n", data_suma, data_sumb); fflush(stdout);
    }
  }
  /* Broadcast the status to the other nodes */
  broadcast_bytes((char *)&status, sizeof(status));
  // printf("%s(%d): status = %d\n", myname, this_node, status);

  fclose(fpt);

  return status;

} /* read_grid_ks_eigenvector */

/* Read a set of eigenvectors from a QUDA-formated color-spin-field file */

int
read_quda_ks_eigenvectors(QIO_Reader *infile, su3_vector *eigVec[], Real *eigVal, int *Nvecs,
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
    parse_xml_eigVal(eigVal, xml);
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

