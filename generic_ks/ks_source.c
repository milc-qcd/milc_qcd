/************************** ks_source.c *****************************/
/* MIMD version 7 */

/* Initialize a source for the inverter */
/* Initialize a sink */
/* Interpret source and sink parameters */

#include "generic_ks_includes.h"
#include "../include/io_ksprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#endif

/* Utilities */

/* Initialize by calling this function before using ksqs */
void init_ksqs(ks_quark_source *ksqs){
  ksqs->type             = UNKNOWN;
  ksqs->c_src            = NULL;
  ksqs->cv_src           = NULL;
  ksqs->file_initialized = 0;
#ifdef HAVE_QIO
  ksqs->infile           = NULL;
  ksqs->outfile          = NULL;
#endif
  ksqs->x0               = 0;
  ksqs->y0               = 0;
  ksqs->z0               = 0;
  ksqs->t0               = 0;
  ksqs->descrp[0]        = '\0';
}

/* Must be called to free */
static void reset_ksqs(ks_quark_source *ksqs){
  if(ksqs->c_src != NULL){free(ksqs->c_src); ksqs->c_src = NULL;}
  if(ksqs->cv_src != NULL){free(ksqs->cv_src); ksqs->cv_src = NULL;}
}

void clear_ksqs(ks_quark_source *ksqs){
  reset_ksqs(ksqs);
#ifdef HAVE_QIO
  if(ksqs->infile != NULL){
    if(ksqs->type == COMPLEX_FIELD_FILE)
      r_close_complex_scidac_file(ksqs->infile);
    else if(ksqs->type == VECTOR_FIELD_FILE)
      r_close_ks_vector_scidac_file(ksqs->infile);
  }
#endif
}

void alloc_ksqs_cv_src(ks_quark_source *ksqs)
{
  if(ksqs->cv_src == NULL)
    ksqs->cv_src = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(ksqs->cv_src == NULL){
    printf("alloc_ksqs_cv_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(ksqs->cv_src,'\0', sizeof(su3_vector)*sites_on_node);
}

void alloc_ksqs_c_src(ks_quark_source *ksqs)
{
  if(ksqs->c_src == NULL)
    ksqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(ksqs->c_src == NULL){
    printf("alloc_ksqs_c_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(ksqs->c_src, '\0', sizeof(complex)*sites_on_node);
}

/* remove "if 0" when we need it */
#if 0
/* Copy c-source from vector */
static void extract_C_from_V(complex *dest, su3_vector *src, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i].c[color];
  }
}
#endif

#ifdef HAVE_QIO
/* Copy vector from c-source */
static void insert_V_from_C(su3_vector *dest, complex *src, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i].c[color] = src[i];
  }
}

/* Copy color vector source */
static void copy_V(su3_vector *dest, su3_vector *src)
{
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i];
  }
}
#endif

/* Choose the specific USQCD format appropriate to the source type */
int choose_usqcd_ks_file_type(int source_type){
  int file_type;

  switch(source_type){
  case UNKNOWN:
  case POINT:
  case CORNER_WALL:
  case EVEN_WALL:
  case EVENANDODD_WALL:
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_STORE:
    file_type = FILE_TYPE_KS_USQCD_C1V3;
    break;
  case RANDOM_VECTOR_WALL:
  case VECTOR_FIELD_FILE:
  case VECTOR_FIELD_STORE:
    file_type = FILE_TYPE_KS_USQCD_VV_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

static void insert_mom_V(su3_vector *v, int mom[3], int t0){
  int i; site *s;
  int c;
  Real pi = M_PI;
  Real th;
  Real px = 2*pi*mom[0]/nx;
  Real py = 2*pi*mom[1]/ny;
  Real pz = 2*pi*mom[2]/nz;
  complex z,y;
  
  if(mom[0] == 0 && mom[1] == 0 && mom[2] == 0)return;

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*s->x + py*s->y + pz*s->z;
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++){
      z = v[i].c[c];
      CMUL(z,y,v[i].c[c]);
    }
  }
}

static void insert_mom_C(complex *c, int mom[3], int t0){
  int i; site *s;
  Real pi = M_PI;
  Real th;
  Real px = 2*pi*mom[0]/nx;
  Real py = 2*pi*mom[1]/ny;
  Real pz = 2*pi*mom[2]/nz;
  complex z,y;
  
  if(mom[0] == 0 && mom[1] == 0 && mom[2] == 0)return;

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*s->x + py*s->y + pz*s->z;
    y.real = cos(th); y.imag = sin(th);
    z = c[i];
    CMUL(z,y,c[i]);
  }
}

#ifdef HAVE_QIO
/********************************************************************/
/* Parse the record XML to get the color and check it               */
/********************************************************************/
/* For a KS source file we use the same XML encoding as with the KS
   propagator file  */
static int check_color(QIO_String *recxml, int color){
  int status;
  int input_color;
  QIO_USQCDKSPropRecordInfo recinfo;
  char myname[] = "check_color";

  status = QIO_decode_usqcd_ksproprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS) 
    return qio_status(status);
  input_color = QIO_get_usqcd_ksproprecord_color(&recinfo);
  if(color != input_color){
    node0_printf("%s(%d): Error: expected color %d got %d\n",
		 myname, this_node, color, input_color);
    return 1;
  }
  return 0;
}

#endif

/* Translate QIO status to our conventions */
/* 0 = success
  -1 = end of file
   1 = other failure
*/

/********************************************************************/
/* Construct the source field */
/********************************************************************/

void r_open_ks_source(ks_quark_source *ksqs){

#ifdef HAVE_QIO
  char *source_file         = ksqs->source_file;
  QIO_String *xml_file;
#endif
  char myname[] = "r_open_ks_source";

  if(ksqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(ksqs->type == COMPLEX_FIELD_FILE){
    xml_file = QIO_string_create();
    ksqs->infile = r_open_complex_scidac_file_xml(source_file, QIO_SERIAL,
						  xml_file);
  }
  else if(ksqs->type == VECTOR_FIELD_FILE){
    xml_file = QIO_string_create();
    ksqs->infile = r_open_ks_vector_scidac_file_xml(source_file, QIO_SERIAL,
						    xml_file);
  }
  else {
    node0_printf("%s: a file-type source is required\n",myname);
    return;
  }
  if(ksqs->infile == NULL)terminate(1);
  ksqs->file_initialized = 1;
  node0_printf("Source file info\n%s\n",QIO_string_ptr(xml_file));
  QIO_string_destroy(xml_file);
  
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

int ks_source_field(su3_vector *src, ks_quark_source *ksqs)
{
  register int i;
  register site *s; 
  char myname[] = "ks_source_field";
  
  int status = 0;
  
  /* Unpack structure */
  int color                 = ksqs->color;
  int source_type           = ksqs->type;
  int x0                    = ksqs->x0; 
  int y0                    = ksqs->y0; 
  int z0                    = ksqs->z0; 
  int t0                    = ksqs->t0;
#ifdef HAVE_QIO
  int *mom                  = ksqs->mom;
  int file_initialized      = ksqs->file_initialized;
#endif
  
  /* zero src to be safe */
  if(src != NULL){
    FORALLSITES(i,s) { clearvec( src+i ); }
  }
  
  /* Unless we are taking the source from storage, clear it before
     reconstructing it */
  if(source_type != COMPLEX_FIELD_STORE &&
     source_type != VECTOR_FIELD_STORE)
    reset_ksqs(ksqs);

  if(source_type == UNKNOWN) {
    /* Default source */
    alloc_ksqs_c_src(ksqs);
  }
  else if(source_type == POINT) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    /* Save a copy of the source in ksqs->c_src */

    alloc_ksqs_c_src(ksqs);

    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      ksqs->c_src[i].real = 1.0;
      if(src != NULL)src[i].c[color].real = 1.0;
    }
  }
  else if(source_type == COMPLEX_FIELD_FILE){
#ifdef HAVE_QIO
    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    alloc_ksqs_c_src(ksqs);
    status = qio_status(read_complex_scidac(ksqs->infile, ksqs->c_src, 1 ));
    
    /* Do momentum insertion */
    if(status == 0)insert_mom_C(ksqs->c_src, mom, t0);
    if(src != NULL)insert_V_from_C(src, ksqs->c_src, color);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }
  else if(source_type == VECTOR_FIELD_FILE){
#ifdef HAVE_QIO
    QIO_String *recxml;

    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    recxml = QIO_string_create();
    /* Load the SU(3) vector from the file */
    alloc_ksqs_cv_src(ksqs);
    status = qio_status(
       read_ks_vector_scidac_xml(ksqs->infile, ksqs->cv_src, 1, recxml) );

    /* Do momentum insertion */
    if(status == 0)insert_mom_V(ksqs->cv_src, mom, t0);

    /* Verify that the requested color matches the color encoded in
       the record XML */
    if(status == 0)
      status = check_color(recxml, color);

    QIO_string_destroy(recxml);
    /* Copy to requested location */
    if(src != NULL && status == 0)
      copy_V(src, ksqs->cv_src);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }
  else if(source_type == COMPLEX_FIELD_STORE){
    if(ksqs->c_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    /* Load to the specified color and timeslice(s) */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i].c[color] = ksqs->c_src[i];
      }
    }
  }
  else if(source_type == VECTOR_FIELD_STORE){
    if(ksqs->cv_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    /* Ignore color and load the whole SU(3) vector from storage */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i] = ksqs->cv_src[i];
      }
    }
  }
  else {
    node0_printf("%s: Unrecognized source type %d\n",myname, source_type);
    terminate(1);
  }
  return status;

} /* ks_source_field */

void r_close_ks_source(ks_quark_source *ksqs){

  char myname[] = "r_close_ks_source";

  if(ksqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(ksqs->type == COMPLEX_FIELD_FILE){
    r_close_complex_scidac_file(ksqs->infile);
    ksqs->file_initialized = 0;
    ksqs->infile = NULL;
  }
  else if(ksqs->type == VECTOR_FIELD_FILE){
    r_close_ks_vector_scidac_file(ksqs->infile);
    ksqs->file_initialized = 0;
    ksqs->infile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}


void w_open_ks_source(ks_quark_source *ksqs, char *fileinfo){

  char myname[] = "w_open_ks_source";

#ifdef HAVE_QIO
  int volfmt, serpar;
  char *source_file = ksqs->source_file;

  interpret_usqcd_ks_save_flag(&volfmt, &serpar, ksqs->flag);

  if(ksqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }

  if(ksqs->type == COMPLEX_FIELD_FILE){
    ksqs->outfile = w_open_complex_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    ksqs->file_initialized = 1;
  }
  if(ksqs->type == VECTOR_FIELD_FILE){
    ksqs->outfile = w_open_ks_vector_scidac_file(source_file, fileinfo,
						 volfmt, serpar);
    ksqs->file_initialized = 1;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

/* For now we support writing only a color vector source field */

int ks_source_write(su3_vector *src, ks_quark_source *ksqs)
{
  char myname[] = "ks_source_write";
  
  int status = 0;
  double dtime = 0;
  
  /* Unpack structure */
  int source_type           = ksqs->type;
  int file_initialized      = ksqs->file_initialized;
#ifdef HAVE_QIO
  int t0                    = ksqs->t0;
  int color                 = ksqs->color;
  QIO_String *recxml;
  QIO_USQCDKSPropRecordInfo *recinfo;
#endif
  
  dtime = -dclock();

  if(source_type != VECTOR_FIELD_FILE){
    node0_printf("%s: Unrecognized source type for writing\n", myname);
    return 1;
  }

  if(file_initialized == 0){
    node0_printf("%s: File has not been initialized\n", myname);
    return 1;
  }
#ifdef HAVE_QIO
  recxml = QIO_string_create();

  /* Construct the record XML */
  recinfo = QIO_create_usqcd_ksproprecord_c_info(color, "");
  QIO_encode_usqcd_ksproprecord_info(recxml, recinfo);
  QIO_destroy_usqcd_ksproprecord_info(recinfo);
  /* Write the SU(3) source vector to the file */
  status = qio_status(
      write_kspropsource_V_usqcd_xml(ksqs->outfile, recxml, src, t0) );
  node0_printf("Wrote source for color %d time slice %d\n", color, t0);
  QIO_string_destroy(recxml);
  dtime += dclock();
  node0_printf("Time to save source color %d = %e\n", color,dtime);
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
}

void w_close_ks_source(ks_quark_source *ksqs){

  char myname[] = "w_close_ks_source";

  if(ksqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(ksqs->type == COMPLEX_FIELD_FILE){
    w_close_complex_scidac_file(ksqs->outfile);
    ksqs->file_initialized = 0;
    ksqs->outfile = NULL;
  }
  else if(ksqs->type == VECTOR_FIELD_FILE){
    w_close_ks_vector_scidac_file(ksqs->outfile);
    ksqs->file_initialized = 0;
    ksqs->outfile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

int ks_source_site(field_offset src, ks_quark_source *ksqs)
{
  int i;
  site *s;
  su3_vector *t_src;
  int status;

#define PAD 0
  t_src  = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  
  if(t_src == NULL){
    printf("ks_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  status = ks_source_field(t_src, ksqs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(su3_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);

  return status;

} /* ks_source_site */


/* 
   Initialize a sink using an SU(3) vector as the storage container .
 */

void ks_sink_field(complex *snk, ks_quark_source *ksqs)
{
  int i;
  
  /* Unpack structure.  We don't use member t0 here. */
  int sink_type     = ksqs->type;
  int x0            = ksqs->x0; 
  int y0            = ksqs->y0; 
  int z0            = ksqs->z0; 
  int t0            = ALL_T_SLICES;
  char *sink_file   = ksqs->source_file;
  int *mom          = ksqs->mom;

  if(sink_type == POINT) {
    /* load 1.0 into sink at cooordinates given by sink_coord */
    /* initialize snk to be a delta function at point x0,y0,z0 */
    
    for(t0=0;t0<nt;t0++){
      if(node_number(x0,y0,z0,t0) == mynode()){
	i = node_index(x0,y0,z0,t0);
	snk[i].real = 1.0;
      }
    }
  }
  else if(sink_type == COMPLEX_FIELD_FILE){
  /* CAREFUL: The SciDAC file must have an identical sink wave
     function on each time slice */
#ifdef HAVE_QIO
    restore_complex_scidac_to_field(sink_file, QIO_SERIAL, snk, 1);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, t0);
#else
    node0_printf("ks_sink_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(sink_type == COMPLEX_FIELD_FM_FILE){
    r_source_cmplx_fm_to_field(sink_file, snk, t0, sink_type);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, t0);
  }
  else {
    node0_printf("ks_sink_field: Unrecognized sink type %d\n",sink_type);
    terminate(1);
  }
  
} /* ks_sink_field */


/* snk must be type complex */
void ks_sink_site(field_offset snk, ks_quark_source *ksqs)
{
  int i;
  site *s;
  complex *t_snk;

#define PAD 0
  t_snk  = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));
  
  if(t_snk == NULL){
    printf("ks_sink_site(%d): Can't allocate t_snk\n",this_node);
    terminate(1);
  }
  
  ks_sink_field(t_snk, ksqs);

  /* copy snk back */
  FORALLSITES(i,s) {
    *(complex *)F_PT(s,snk) = t_snk[i];
  }

  free(t_snk);

} /* ks_sink_site */

int ask_ks_quark_source( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_ks_quark_source";

  if (prompt!=0){
    printf("enter 'point', 'corner_wall', 'even_wall',");
    printf("'evenandodd_wall', 'random_vector_wall', 'complex_field',");
    printf("'vector_field'");
    printf(", for source type\n");
  }

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("corner_wall",savebuf) == 0 ) {
    *source_type = CORNER_WALL;
    strcpy(descrp,"corner wall");
  }
  else if(strcmp("even_wall",savebuf) == 0 ) {
    *source_type = EVEN_WALL;
    strcpy(descrp,"even wall");
  }
  else if(strcmp("evenandodd_wall",savebuf) == 0 ) {
    *source_type = EVENANDODD_WALL;
    strcpy(descrp,"even and odd wall");
  }
  /* Source file is a series of random color vector fields over
     the entire lattice: SciDAC format */
  else if(strcmp("random_vector_wall",savebuf) == 0 ) {
    *source_type = RANDOM_VECTOR_WALL;
    strcpy(descrp,"random vector wall");
  }
  /* Source file is a complex field over the entire lattice: SciDAC
     format */
  else if(strcmp("complex_field",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FILE;
    strcpy(descrp,"complex_field");
  }
  else if(strcmp("complex_field_fm",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FM_FILE;
    strcpy(descrp,"complex_field_fm");
  }
  /* Source file is an SU(3) vector field over the entire lattice: SciDAC
     format */
  else if(strcmp("vector_field",savebuf) == 0 ) {
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
  }
  else{
    printf("%s: ERROR IN INPUT: KS source command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_ks_quark_source */

#define IF_OK if(status==0)

int ask_output_ks_quark_source_file( FILE *fp, int prompt, 
				     int *flag, int *source_type,
				     int *t0, char *descrp, char *filename)
{
  char *savebuf;
  int status = 0;
  char myname[] = "ask_output_ks_quark_source_file";

  if (prompt!=0){
    printf("'save_serial_scidac_ks_source' or ");
    printf("'save_multifile_scidac_ks_source' or ");
    printf("'save_partfile_scidac_ks_source' or ");
    printf("'save_serial_scidac_w_source' or ");
    printf("'save_multifile_scidac_w_source' or ");
    printf("'save_partfile_scidac_w_source'? ");
  }

  savebuf = get_next_tag(fp, "output quark source command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("save_serial_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_serial_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else{
    printf("\n%s: ERROR IN INPUT: KS source file command %s is invalid\n",
	   myname, savebuf); 
    return 1;
  }
  
  /* Get file name */
  if( *flag != FORGET ){
    if(prompt!=0)printf("enter filename\n");
    if(scanf("%s",filename) != 1){
      printf("\n%s(%d): ERROR IN INPUT: Can't read filename\n",
	     myname, this_node); 
      status++;
    }
    else
      printf("%s\n",filename);
  }

  /* Get time slice*/
  IF_OK status += get_i(stdin, prompt, "t0", t0);

  return status;
} /* ask_output_ks_quark_source_file */

/* Get the additional input parameters needed to specify the source */
int get_ks_quark_source(FILE *fp, int prompt, ks_quark_source *ksqs){
  
  int  source_type;
  int  source_loc[4] = { 0,0,0,0 };
  char source_file[MAXFILENAME] = "";
  char source_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_ks_quark_source(fp, prompt,&source_type,
				      ksqs->descrp);
  IF_OK ksqs->type  = source_type;
  
  IF_OK {
    if ( source_type == POINT ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
    }
    else if ( source_type == CORNER_WALL ||
	      source_type == EVEN_WALL ||
	      source_type == EVENANDODD_WALL ||
	      source_type == RANDOM_VECTOR_WALL ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
    }
    else if ( source_type == COMPLEX_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == COMPLEX_FIELD_FM_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == VECTOR_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else {
      printf("Source type not supported in this application\n");
      status++;
    }
  }
  
  IF_OK status += get_s(stdin, prompt, "source_label", source_label);

  ksqs->x0    = source_loc[0];
  ksqs->y0    = source_loc[1];
  ksqs->z0    = source_loc[2];
  ksqs->t0    = source_loc[3];
  strcpy(ksqs->source_file,source_file);
  strncpy(ksqs->label,source_label,MAXSRCLABEL);
  
  return status;
} /* get_ks_quark_source */

int get_ks_quark_sink(FILE *fp, int prompt, ks_quark_source *ksqs){

  int  sink_type;
  int  sink_loc[3] = { 0,0,0 };
  char sink_file[MAXFILENAME] = "";
  char sink_label[MAXSRCLABEL];
  int  status = 0;

  /* Get antiquark source type */
  IF_OK status += ask_ks_quark_source(fp, prompt,&sink_type,
				      ksqs->descrp);
  IF_OK ksqs->type  = sink_type;
  
  IF_OK {
    if ( sink_type == POINT ){
      IF_OK status += get_vi(stdin, prompt, "origin", sink_loc, 3);
    }
    else if ( sink_type == COMPLEX_FIELD_FILE ){
      IF_OK status += get_s(stdin, prompt, "load_source", sink_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( sink_type == COMPLEX_FIELD_FM_FILE ){
      IF_OK status += get_s(stdin, prompt, "load_source", sink_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else {
      printf("Sink type not supported in this application\n");
      status++;
    }
  }

  IF_OK status += get_s(stdin, prompt, "sink_label", sink_label);
  
  ksqs->x0    = sink_loc[0];
  ksqs->y0    = sink_loc[1];
  ksqs->z0    = sink_loc[2];
  ksqs->t0    = ALL_T_SLICES;
  ksqs->color = 0;
  /* For a sink "source_file" is a misnomer */
  strcpy(ksqs->source_file,sink_file);
  strncpy(ksqs->label,sink_label,MAXSRCLABEL);
  
  return status;
} /* get_ks_quark_sink */


