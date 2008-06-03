/************************** w_source.c *****************************/
/* MIMD version 7 */

/* Initialize a source for the inverter */
/* Initialize a sink */
/* Interpret source and sink parameters */

#include "generic_wilson_includes.h"
#include "../include/io_wprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#endif

/* Utilities */


void init_wqs(wilson_quark_source *wqs){
  wqs->type             = UNKNOWN;
  wqs->c_src            = NULL;
  wqs->wv_src           = NULL;
  wqs->file_initialized = 0;
#ifdef HAVEQIO
  wqs->infile           = NULL;
  wqs->outfile          = NULL;
#endif
  wqs->x0               = 0;
  wqs->y0               = 0;
  wqs->z0               = 0;
  wqs->t0               = 0;
  wqs->descrp[0]        = '\0';
}


/* Must be called to free */
static void reset_wqs(wilson_quark_source *wqs){
  if(wqs->c_src != NULL){free(wqs->c_src); wqs->c_src = NULL;}
  if(wqs->wv_src != NULL){free(wqs->wv_src); wqs->wv_src = NULL;}
}

void clear_wqs(wilson_quark_source *wqs){
  reset_wqs(wqs);
#ifdef HAVE_QIO
  if(wqs->infile != NULL){
    if(wqs->type == COMPLEX_FIELD_FILE)
      r_close_complex_scidac_file(wqs->infile);
    else if(wqs->type == VECTOR_FIELD_FILE)
      r_close_w_vector_scidac_file(wqs->infile);
  }
#endif
}

void alloc_wqs_wv_src(wilson_quark_source *wqs)
{
  if(wqs->wv_src == NULL)
    wqs->wv_src = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(wqs->wv_src == NULL){
    printf("alloc_wqs_wv_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->wv_src, 0, sizeof(wilson_vector)*sites_on_node);
}

void alloc_wqs_c_src(wilson_quark_source *wqs)
{
  if(wqs->c_src == NULL)
    wqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(wqs->c_src == NULL){
    printf("alloc_wqs_c_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->c_src, 0, sizeof(complex)*sites_on_node);
}

/* Copy c-source from wv */
static void extract_C_from_D(complex *dest, wilson_vector *src, 
			     int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i].d[spin].c[color];
  }
}

/* Copy wv from c-source */
static void insert_D_from_C(wilson_vector *dest, complex *src, 
			    int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i].d[spin].c[color] = src[i];
  }
}

/* Copy wv-source */
static void copy_D(wilson_vector *dest, wilson_vector *src)
{
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i];
  }
}

/* Choose the specific USQCD format appropriate to the source type */
int choose_usqcd_w_file_type(int source_type){
  int file_type;

  switch(source_type){
  case UNKNOWN:
  case POINT:
  case GAUSSIAN:
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_FM_FILE:
  case COMPLEX_FIELD_STORE:
    file_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case ROTATE_3D:
  case COVARIANT_GAUSSIAN:
  case DIRAC_FIELD_FILE:
  case DIRAC_FIELD_FM_FILE:
  case DIRAC_FIELD_STORE:
    file_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

static void insert_mom_D(wilson_vector *wv, int mom[3], int t0){
  int i; site *s;
  int c,d;
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
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	z = wv[i].d[d].c[c];
	CMUL(z,y,wv[i].d[d].c[c]);
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
/* For the source file we borrow the XML encoding from the USQCD
   Wilson propagator file  */
static int check_color_spin(QIO_String *recxml, int color, int spin){
  int status;
  int input_color, input_spin;
  QIO_USQCDPropRecordInfo recinfo;
  char myname[] = "check_color_spin";

  status = QIO_decode_usqcd_proprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS) 
    return qio_status(status);
  input_color = QIO_get_usqcd_proprecord_color(&recinfo);
  input_spin = QIO_get_usqcd_proprecord_spin(&recinfo);
  if(color != input_color || spin  != input_spin ){
    node0_printf("%s(%d): Error: expected color %d and spin %d got %d and %d\n",
		 myname, this_node, color, spin, 
		 input_color, input_spin);
    return 1;
  }
  return 0;
}

#endif

/*--------------------------------------------------------------------*/
/* Do the 3D Fermilab rotation on a Wilson vector field.  That is,
   compute src <- (1 + a d1 * \gamma * D/2) src where D is a 3D
   covariant difference normalized to give twice the covariant
   derivative in the continuum limit. */
static void rotate_3D_wvec(wilson_vector *src, Real d1)
{
  int i;
  site *s;
  wilson_vector *mp, *tmp;
  
  mp  = create_wv_field();
  tmp = create_wv_field();
  
  /* Do Wilson Dslash on the source field */
  dslash_w_3D_field(src, mp,  PLUS, EVENANDODD);
  dslash_w_3D_field(src, tmp, MINUS, EVENANDODD);
  
  FORALLSITES(i,s){
    /* tmp <- mp - tmp = 2*Dslash*src */
    sub_wilson_vector(mp + i, tmp + i, tmp + i);
    /* src <- d1/4 * tmp + src */
    scalar_mult_add_wvec(src + i, tmp + i, d1/4., src + i);
  }

  cleanup_dslash_w_3D_temps();
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}

/********************************************************************/
/* Construct the source field */
/********************************************************************/

void r_open_w_source(wilson_quark_source *wqs){
  
#ifdef HAVE_QIO
  char *source_file         = wqs->source_file;
  QIO_String *xml_file;
#endif
  char myname[] = "r_open_w_source";

  if(wqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    xml_file = QIO_string_create();
    wqs->infile = r_open_complex_scidac_file_xml(source_file, QIO_SERIAL,
						 xml_file);
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    xml_file = QIO_string_create();
    wqs->infile = r_open_w_vector_scidac_file_xml(source_file, QIO_SERIAL,
						  xml_file);
  }
  else{
    node0_printf("%s: a file-type source is required\n",myname);
    return;
  }
  if(wqs->infile == NULL)terminate(1);
  wqs->file_initialized = 1;
  node0_printf("Source file info\n%s\n",QIO_string_ptr(xml_file));
  QIO_string_destroy(xml_file);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}


int w_source_field(wilson_vector *src, wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  
  short my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  int status = 0;
  
  /* Unpack structure */
  int color                 = wqs->color;
  int spin                  = wqs->spin;
  int source_type           = wqs->type;
  int x0                    = wqs->x0; 
  int y0                    = wqs->y0; 
  int z0                    = wqs->z0; 
  int t0                    = wqs->t0;
  Real r0                   = wqs->r0;
  int iters                 = wqs->iters;
  char *source_file         = wqs->source_file;
  int *mom                  = wqs->mom;
#ifdef HAVE_QIO
  int file_initialized      = wqs->file_initialized;
#endif
  
  /* zero src to be safe */
  if(src != NULL){
    FORALLSITES(i,s) { clear_wvec( src+i ); }
  }
  
  /* Unless we are taking the source from storage, clear it before
     reconstructing it */
  if(source_type != COMPLEX_FIELD_STORE &&
     source_type != DIRAC_FIELD_STORE)
    reset_wqs(wqs);
  
  if(source_type == UNKNOWN) {
    /* Dummy source */
    alloc_wqs_c_src(wqs);
  }
  else if(source_type == POINT) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    /* Save a copy of the source in wqs->c_src */

    alloc_wqs_c_src(wqs);

    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      wqs->c_src[i].real = 1.0;
      if(src != NULL)src[i].d[spin].c[color].real = 1.0;
    }
  }
  else if(source_type == GAUSSIAN) {
    /* Gaussian trial source centered on  x0,y0,z0,t0 */
    /* Save a copy of the source in wqs->c_src */

    alloc_wqs_c_src(wqs);

    FORALLSITES(i,s) {
      if(s->t != t0)continue;	/* only do this if t==t0 */
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);

      wqs->c_src[i].real = (Real)exp((double)(- radius2));
      if(src != NULL)
	src[i].d[spin].c[color].real = wqs->c_src[i].real;
    }
  }
  else if(source_type == ROTATE_3D){
    alloc_wqs_wv_src(wqs);

    /* Set delta function source */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      wqs->wv_src[i].d[spin].c[color].real = 1.;
    }
    /* Then do 3D rotation */
    rotate_3D_wvec(wqs->wv_src, wqs->d1);

    /* Copy to requested location */
    if(src != NULL)copy_D(wqs->wv_src, src);
  }
  else if(source_type == COVARIANT_GAUSSIAN){
    alloc_wqs_wv_src(wqs);

    /* Set delta function source */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      wqs->wv_src[i].d[spin].c[color].real = 1.;
    }
    /* Then smear */
    gauss_smear_field(wqs->wv_src, r0, iters, t0);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == COMPLEX_FIELD_FM_FILE){
    /* Load to the specified spin, color and timeslice */
    /* zero src to be safe */
    alloc_wqs_wv_src(wqs);
    r_source_w_fm_to_field(source_file, wqs->wv_src, spin, color, t0, 
			   source_type);
    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, t0);

    /* Save a copy as a complex field */
    alloc_wqs_c_src(wqs);
    extract_C_from_D(wqs->c_src, wqs->wv_src, spin, color);

    /* Copy to requested location */
    if(src != NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == COMPLEX_FIELD_FILE){
#ifdef HAVE_QIO
    if(file_initialized == 0)
      r_open_w_source(wqs);

    alloc_wqs_c_src(wqs);
    status = qio_status(read_complex_scidac(wqs->infile, wqs->c_src, 1 ));

    /* Do momentum insertion */
    if(status == 0)insert_mom_C(wqs->c_src, mom, t0);

    if(src!=NULL)insert_D_from_C(src, wqs->c_src, spin, color);
#else
    node0_printf("w_source_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(source_type == DIRAC_FIELD_FM_FILE){
    alloc_wqs_wv_src(wqs);
    r_source_w_fm_to_field(source_file, wqs->wv_src, spin, color, 
			   t0, source_type);
    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, t0);

    /* Copy to requested location */
    if(src!=NULL)copy_D(src, wqs->wv_src);
  }
  else if(source_type == DIRAC_FIELD_FILE){
#ifdef HAVE_QIO
    QIO_String *recxml;

    if(file_initialized == 0)
      r_open_w_source(wqs);

    recxml = QIO_string_create();
    /* Load the Dirac vector from the file */
    alloc_wqs_wv_src(wqs);
    status = qio_status(
       read_w_vector_scidac_xml(wqs->infile, wqs->wv_src, 1, recxml));

    /* Do momentum insertion */
    if(status == 0)insert_mom_D(wqs->wv_src, mom, t0);

    /* Verify that the requested color matches the color encoded in
       the record XML */
    if(status == 0)
      status = check_color_spin(recxml, color, spin);

    QIO_string_destroy(recxml);

    /* Copy to requested location */
    if(src != NULL && status == 0)
      copy_D(src, wqs->wv_src);
#else
    node0_printf("w_source_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(source_type == COMPLEX_FIELD_STORE){
    if(wqs->c_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Load to the specified spin, color and timeslice */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i].d[spin].c[color] = wqs->c_src[i];
      }
    }
  }
  else if(source_type == DIRAC_FIELD_STORE){
    if(wqs->wv_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Ignore spin, color and load the whole wilson vector from storage */
    if(src != NULL){
      FORALLSITES(i,s){
	if(t0 == ALL_T_SLICES || s->t == t0)
	  src[i] = wqs->wv_src[i];
      }
    }
  }
  else {
    node0_printf("w_source: Unrecognized source type %d\n",source_type);
    terminate(1);
  }

  return status;

} /* w_source_field */

void r_close_w_source(wilson_quark_source *wqs){

  char myname[] = "r_close_w_source";

  if(wqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    r_close_complex_scidac_file(wqs->infile);
    wqs->file_initialized = 0;
    wqs->infile = NULL;
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    r_close_w_vector_scidac_file(wqs->infile);
    wqs->file_initialized = 0;
    wqs->infile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}


void w_open_w_source(wilson_quark_source *wqs, char *fileinfo){

  char myname[] = "w_open_w_source";
#ifdef HAVE_QIO
  char *source_file = wqs->source_file;
  int volfmt, serpar;
#endif

  if(wqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  interpret_usqcd_w_save_flag(&volfmt, &serpar, wqs->flag);

  if(wqs->type == COMPLEX_FIELD_FILE){
    wqs->outfile = w_open_complex_scidac_file(source_file, fileinfo,
					      volfmt, serpar);
    wqs->file_initialized = 1;
  }
  if(wqs->type == DIRAC_FIELD_FILE){
    wqs->outfile = w_open_w_vector_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    wqs->file_initialized = 1;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

/* For now we support writing only a Dirac vector source field */

int w_source_write(wilson_vector *src, wilson_quark_source *wqs)
{
  char myname[] = "w_source_write";
  
  int status = 0;
  
  /* Unpack structure */
  int source_type           = wqs->type;
  int file_initialized      = wqs->file_initialized;
#ifdef HAVE_QIO
  int t0                    = wqs->t0;
  int color                 = wqs->color;
  int spin                  = wqs->spin;
  QIO_String *recxml;
  QIO_USQCDPropRecordInfo *recinfo;
#endif
  
  if(source_type != DIRAC_FIELD_FILE){
    node0_printf("%s: Unrecognized source type for writing\n", myname);
    return 1;
  }

  if(file_initialized == 0){
    node0_printf("%s: File has not been initialized\n", myname);
    return 1;
  }
#ifdef HAVE_QIO
  recxml = QIO_string_create();

  /* Construct the record XML - we borrow the USQCD prop record XML */
  recinfo = QIO_create_usqcd_proprecord_sc_info(spin, color, "");
  QIO_encode_usqcd_proprecord_info(recxml, recinfo);
  QIO_destroy_usqcd_proprecord_info(recinfo);
  /* Write the Dirac source vector to the file */
  status = qio_status(
      write_wpropsource_D_usqcd_xml(wqs->outfile, recxml, src, t0) );
  node0_printf("Wrote source for color %d spin %d time slice %d\n", 
	       color, spin, t0);
  QIO_string_destroy(recxml);
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
}

void w_close_w_source(wilson_quark_source *wqs){

  char myname[] = "w_close_w_source";

  if(wqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(wqs->type == COMPLEX_FIELD_FILE){
    w_close_complex_scidac_file(wqs->outfile);
    wqs->file_initialized = 0;
    wqs->outfile = NULL;
  }
  else if(wqs->type == DIRAC_FIELD_FILE){
    w_close_w_vector_scidac_file(wqs->outfile);
    wqs->file_initialized = 0;
    wqs->outfile = NULL;
  }
  else
    node0_printf("%s: a file-type source is required\n",myname);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

int w_source_site(field_offset src, wilson_quark_source *wqs)
{
  int i;
  site *s;
  wilson_vector *t_src;
  int status;

#define PAD 0
  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  
  if(t_src == NULL){
    printf("w_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  status = w_source_field(t_src, wqs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);
  return status;

} /* w_source_site */


/* 
   Initialize a sink.  The complex sink wave function is replicated on
   each time slice.
 */

void w_sink_field(complex *snk, wilson_quark_source *wqs)
{
  int i;
  site *s; 
  
  int my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  
  /* Unpack structure. */
  int sink_type     = wqs->type;
  int x0            = wqs->x0; 
  int y0            = wqs->y0; 
  int z0            = wqs->z0; 
  int t0            = ALL_T_SLICES;
  Real r0           = wqs->r0;
  char *sink_file   = wqs->source_file;
  int *mom          = wqs->mom;

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
  else if(sink_type == GAUSSIAN) {
    /* Gaussian trial sink centered on (x0,y0,z0) of each timeslice */
    
    FORALLSITES(i,s) {
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);
      
      snk[i].real = (Real)exp((double)(- radius2));
      snk[i].imag = 0;
    }
  }
  else if(sink_type == COMPLEX_FIELD_FM_FILE){
    r_source_cmplx_fm_to_field(sink_file, snk, t0, sink_type);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, t0);
  }

  else if(sink_type == COMPLEX_FIELD_FILE){
  /* CAREFUL: This SciDAC file must have an identical sink wave
     function on each time slice */
#ifdef HAVE_QIO
    restore_complex_scidac_to_field(sink_file, QIO_SERIAL, snk, 1);

    /* Do momentum insertion */
    insert_mom_C(snk, mom, t0);
#else
    node0_printf("w_sink_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else{
    node0_printf("w_sink_field: Unknown sink type %d\n",sink_type);
    terminate(1);
  }
} /* w_sink_field */


/* 
   Initialize a sink for the inverter, using a wilson number as the
   storage container.
 */

/* snk must be type complex */
void w_sink_site(field_offset snk, wilson_quark_source *wqs)
{
  int i;
  site *s;
  complex *t_snk;

#define PAD 0
  t_snk  = (complex *) malloc((sites_on_node+PAD)*sizeof(complex));
  
  if(t_snk == NULL){
    printf("w_sink_site(%d): Can't allocate t_snk\n",this_node);
    terminate(1);
  }
  
  w_sink_field(t_snk, wqs);

  /* copy snk back */
  FORALLSITES(i,s) {
    *(complex *)F_PT(s,snk) = t_snk[i];
  }

  free(t_snk);

} /* w_sink_site */


int ask_w_quark_source( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_w_quark_source";

  if (prompt!=0)
    printf("enter 'point', 'rotate_3D', 'gaussian', 'covariant_gaussian', 'complex_field', 'complex_field_fm', 'dirac_field', 'dirac_field_fm' for source type\n");

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("rotate_3D",savebuf) == 0 ){
    *source_type = ROTATE_3D;
    strcpy(descrp,"rotate_3D");
  }
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else if(strcmp("covariant_gaussian",savebuf) == 0 ) {
    *source_type = COVARIANT_GAUSSIAN;
    strcpy(descrp,"covariant_gaussian");
  }
  else if(strcmp("complex_field",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FILE;
    strcpy(descrp,"complex_field");
  }
  else if(strcmp("complex_field_fm",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FM_FILE;
    strcpy(descrp,"complex_field_fm");
  }
  else if(strcmp("dirac_field",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
  }
  else if(strcmp("dirac_field_fm",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FM_FILE;
    strcpy(descrp,"dirac_field_fm");
  }
  else{
    printf("%s: ERROR IN INPUT: Wilson source command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_w_quark_source */

#define IF_OK if(status==0)

int ask_output_w_quark_source_file( FILE *fp, int prompt, 
				    int *flag, int *source_type,
				    int *t0, char *descrp, char *filename)
{
  char *savebuf;
  int status = 0;
  char myname[] = "ask_output_w_quark_source_file";

  if (prompt!=0){
    printf("enter 'save_serial_scidac_w_source'");
    printf(", for source type\n");
  }

  savebuf = get_next_tag(fp, "output quark source command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("save_serial_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else{
    printf("\n%s: ERROR IN INPUT: Wilson source file command %s is invalid\n",
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

  return 0;

} /* ask_output_w_quark_source_file */

/* Get the additional input parameters needed to specify the source */
int get_w_quark_source(FILE *fp, int prompt, wilson_quark_source *wqs){
  
  Real source_r0 = 0;
  Real d1 = 0;
  int  source_type;
  int  source_loc[4] = { 0,0,0,0 };
  int  source_iters = 0;
  char source_file[MAXFILENAME] = "";
  char source_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_w_quark_source(fp, prompt,&source_type,
				   wqs->descrp);
  IF_OK wqs->type  = source_type;
  
  IF_OK {
    if ( source_type == GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      /* width: psi=exp(-(r/r0)^2) */
      IF_OK status += get_f(stdin, prompt,"r0", &source_r0 );
    }
    else if ( source_type == POINT ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
    }
    else if ( source_type == ROTATE_3D ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "d1", &d1);
    }
    else if ( source_type == COVARIANT_GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);
    }
    else if ( source_type == COMPLEX_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else if ( source_type == COMPLEX_FIELD_FM_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
    }
    else if ( source_type == DIRAC_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else if ( source_type == DIRAC_FIELD_FM_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else {
      printf("Source type not supported in this application\n");
      status++;
    }
  }

  IF_OK status += get_s(stdin, prompt, "source_label", source_label);
  
  wqs->r0    = source_r0;
  wqs->x0    = source_loc[0];
  wqs->y0    = source_loc[1];
  wqs->z0    = source_loc[2];
  wqs->t0    = source_loc[3];
  wqs->iters = source_iters;
  wqs->d1    = d1;
  strcpy(wqs->source_file,source_file);
  strncpy(wqs->label,source_label,MAXSRCLABEL);
  
  return status;
} /* get_w_quark_source */

int get_w_quark_sink(FILE *fp, int prompt, wilson_quark_source *wqs){

  Real sink_r0 = 0;
  int  sink_type;
  int  sink_loc[3] = { 0,0,0 };
  char sink_file[MAXFILENAME] = "";
  char sink_label[MAXSRCLABEL];
  int  status = 0;

  /* Get antiquark source type */
  IF_OK status += ask_w_quark_source(fp, prompt,&sink_type,
				     wqs->descrp);
  IF_OK wqs->type  = sink_type;
  
  IF_OK {
    if ( sink_type == GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", sink_loc, 3);
      /* width: psi=exp(-(r/r0)^2) */
	IF_OK status += get_f(stdin, prompt,"r0", &sink_r0 );
      }
      else if ( sink_type == POINT ){
	IF_OK status += get_vi(stdin, prompt, "origin", sink_loc, 3);
      }
      else if ( sink_type == COMPLEX_FIELD_FILE ||
		sink_type == COMPLEX_FIELD_FM_FILE){
	IF_OK status += get_s(stdin, prompt, "load_sink", sink_file);
	IF_OK status += get_vi(stdin, prompt, "momentum", wqs->mom, 3);
      }
      else {
	printf("Sink type not supported in this application\n");
	status++;
      }
    }
    
    IF_OK status += get_s(stdin, prompt, "sink_label", sink_label);

    wqs->r0    = sink_r0;
    wqs->x0    = sink_loc[0];
    wqs->y0    = sink_loc[1];
    wqs->z0    = sink_loc[2];
    wqs->t0    = ALL_T_SLICES;
    wqs->spin  = 0;
    wqs->color = 0;
    /* For a sink "source_file" is a misnomer */
    strcpy(wqs->source_file,sink_file);
    strncpy(wqs->label,sink_label,MAXSRCLABEL);

    return status;
} /* get_w_quark_sink */


