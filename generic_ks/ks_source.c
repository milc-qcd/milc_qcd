OBSOLETE!!
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
void init_qs(quark_source *ksqs){
  ksqs->type             = UNKNOWN;
  ksqs->c_src            = NULL;
  ksqs->v_src            = NULL;
  ksqs->file_initialized = 0;
#ifdef HAVE_QIO
  ksqs->infile           = NULL;
  ksqs->outfile          = NULL;
#endif
  ksqs->kssf             = NULL;
  /* default values for some */
  ksqs->x0               = 0;
  ksqs->y0               = 0;
  ksqs->z0               = 0;
  ksqs->t0               = 0;
  ksqs->ncolor           = 3;
  ksqs->descrp[0]        = '\0';
}

static void reset_ksqs_color(quark_source *ksqs, int color){
  if(ksqs->v_src != NULL){
    if(ksqs->v_src[color] != NULL)
      free(ksqs->v_src[color]); 
    ksqs->v_src[color] = NULL;
  }
}

/* Must be called to free */
static void reset_ksqs(quark_source *ksqs){
  int color;
  if(ksqs->c_src != NULL){free(ksqs->c_src); ksqs->c_src = NULL;}
  if(ksqs->v_src != NULL){
    for(color = 0; color < ksqs->ncolor; color++){
      reset_ksqs_color(ksqs, color);
    }
    free(ksqs->v_src);
    ksqs->v_src = NULL;
  }
}

void clear_ksqs(quark_source *ksqs){
  reset_ksqs(ksqs);
  /* For VECTOR_FIELD_FM_FILE */
  if(ksqs->kssf != NULL){
    r_source_ks_fm_f(ksqs->kssf);
    ksqs->kssf = NULL;
  }
#ifdef HAVE_QIO
  if(ksqs->infile != NULL){
    if(ksqs->type == COMPLEX_FIELD_FILE)
      r_close_complex_scidac_file(ksqs->infile);
    else if(ksqs->type == VECTOR_FIELD_FILE)
      r_close_ks_vector_scidac_file(ksqs->infile);
  }
#endif
  ksqs->file_initialized = 0;
}

static void alloc_ksqs_color(quark_source *ksqs, int color)
{
  int ic;
  if(ksqs->v_src == NULL){
    ksqs->v_src = (su3_vector **)malloc(sizeof(su3_vector *)*ksqs->ncolor);
    if(ksqs->v_src == NULL){
      printf("alloc_ksqs_color(%d): No room for source field\n",this_node);
      terminate(1);
    }
    for(ic=0;ic<ksqs->ncolor;ic++)
      ksqs->v_src[ic] = NULL;
  }
  ksqs->v_src[color] = create_v_field();
  memset(ksqs->v_src[color],'\0', sizeof(su3_vector)*sites_on_node);
}

void alloc_ksqs_v_src(quark_source *ksqs)
{
  int color;
  if(ksqs->v_src == NULL){
    for(color = 0; color < ksqs->ncolor; color++){
      alloc_ksqs_color(ksqs, color);
    }
  }
}

void alloc_ksqs_c_src(quark_source *ksqs)
{
  if(ksqs->c_src == NULL)
    ksqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(ksqs->c_src == NULL){
    printf("alloc_ksqs_c_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(ksqs->c_src, '\0', sizeof(complex)*sites_on_node);
}

/* Choose the specific USQCD format appropriate to the source type */
int choose_usqcd_ks_file_type(int source_type){
  int file_type;

  switch(source_type){
    /* complex field sources */
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_FM_FILE:
  case COMPLEX_FIELD_STORE:
  case CORNER_WALL:
  case EVEN_WALL:
  case EVENANDODD_WALL:
  case EVENMINUSODD_WALL:
  case GAUSSIAN:
  case POINT:
  case UNKNOWN:
    file_type = FILE_TYPE_KS_USQCD_C1V3;
    break;
    /* vector field sources */
  case FAT_COVARIANT_GAUSSIAN:
  case RANDOM_CORNER_COLOR_WALL:
  case RANDOM_COLOR_WALL:
  case VECTOR_FIELD_FILE:
  case VECTOR_FIELD_FM_FILE:
  case VECTOR_FIELD_STORE:
  case WAVEFUNCTION_FILE:
    file_type = FILE_TYPE_KS_USQCD_VV_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

static void insert_mom_V(su3_vector *v, int mom[3], 
			 int x0, int y0, int z0, int t0){
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
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++){
      z = v[i].c[c];
      CMUL(z,y,v[i].c[c]);
    }
  }
}

static void insert_mom_C(complex *c, int mom[3], 
			 int x0, int y0, int z0, int t0){
  int i; site *s;
  Real pi = M_PI;
  Real th;
  Real px = 2*pi*mom[0]/nx;
  Real py = 2*pi*mom[1]/ny;
  Real pz = 2*pi*mom[2]/nz;
  complex z,y;
  
  if(mom[0] == 0 && mom[1] == 0 && mom[2] == 0)return;

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
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

static void r_open_ks_source(quark_source *ksqs){

  char *source_file         = ksqs->source_file;
#ifdef HAVE_QIO
  QIO_String *xml_file;
#endif
  char myname[] = "r_open_ks_source";

  if(ksqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }
  if(ksqs->type == VECTOR_FIELD_FM_FILE){
    ksqs->kssf = r_source_ks_fm_i(source_file);
    if(ksqs->kssf == NULL){
      node0_printf("%s: Failed to open source %s\n", myname, source_file);
    }
  }
#ifdef HAVE_QIO
  else if(ksqs->type == COMPLEX_FIELD_FILE){
    xml_file = QIO_string_create();
    ksqs->infile = r_open_complex_scidac_file_xml(source_file, QIO_SERIAL,
						  xml_file);
    node0_printf("Source file info\n%s\n",QIO_string_ptr(xml_file));
    QIO_string_destroy(xml_file);
  }
  else if(ksqs->type == VECTOR_FIELD_FILE){
    xml_file = QIO_string_create();
    ksqs->infile = r_open_ks_vector_scidac_file_xml(source_file, QIO_SERIAL,
						    xml_file);
    QIO_string_destroy(xml_file);
  }
#endif
  else {
    node0_printf("%s: can't open a file for source type %d\n", myname, ksqs->type);
    return;
  }
#ifdef HAVE_QIO
  if(ksqs->infile == NULL && ksqs->kssf == NULL)terminate(1);
#else
  if(ksqs->kssf == NULL)terminate(1);
#endif
  ksqs->file_initialized = 1;
}

static void random_corner_wall(su3_vector *v, int t0){
  int i,jc;
  site *s;
  Real x;

  FORALLSITES(i,s){
    if( s->t==t0 && s->x%2==0 && s->y%2==0 && s->z%2==0 ){
      for(jc=0;jc<3;jc++){
	v[i].c[jc].real = gaussian_rand_no(&(s->site_prn));
	v[i].c[jc].imag = gaussian_rand_no(&(s->site_prn));
      }
      x = 2.0/sqrt( magsq_su3vec( v+i ) );
      /* "2" because we invert 2m + 2 Dslash */
      scalar_mult_su3_vector( v+i, x, v+i );
    }
  }
}
    

static void random_wall(su3_vector *v, int t0){
  int i,jc;
  site *s;
  Real x;

  FORALLSITES(i,s){
    if( s->t==t0 ){
      for(jc=0;jc<3;jc++){
	v[i].c[jc].real = gaussian_rand_no(&(s->site_prn));
	v[i].c[jc].imag = gaussian_rand_no(&(s->site_prn));
      }
      x = 2.0/sqrt( magsq_su3vec( v+i ) );
      /* "2" because we invert 2m + 2 Dslash */
      scalar_mult_su3_vector( v+i, x, v+i );
    }
  }
}
    

/* Build a single source vector */

int ks_base_source(su3_vector *src, quark_source *ksqs)
{
  register int i;
  register site *s; 
  char myname[] = "ks_base_source";
  
  int status = 0;
#ifdef HAVE_QIO
  int rshift[4];
#endif
  
  /* Unpack structure */
  int color                 = ksqs->color;
  int source_type           = ksqs->type;
  int x0                    = ksqs->x0; 
  int y0                    = ksqs->y0; 
  int z0                    = ksqs->z0; 
  int t0                    = ksqs->t0;
  int *mom                  = ksqs->mom;
  int file_initialized      = ksqs->file_initialized;
  
  /* zero src to be safe */
  if(src != NULL){
    FORALLSITES(i,s) { clearvec( src+i ); }
  }
  

  if(source_type == UNKNOWN) {
    /* Default source */
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
  }

  else if(source_type == COMPLEX_FIELD_FILE){
    reset_ksqs(ksqs);
#ifdef HAVE_QIO
    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    alloc_ksqs_c_src(ksqs);
    status = qio_status(read_complex_scidac(ksqs->infile, ksqs->c_src, 1 ));

    /* Translate the origin as requested */
    rshift[0] = x0; rshift[1] = y0; rshift[2] = z0; rshift[3] = t0;
    shift_complex(ksqs->c_src, rshift);
    
    /* Do momentum insertion */
    if(status == 0)insert_mom_C(ksqs->c_src, mom, x0, y0, z0, t0);
    if(src != NULL)insert_v_from_c(src, ksqs->c_src, color);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }
  
  else if(source_type == COMPLEX_FIELD_FM_FILE){
    alloc_ksqs_v_src(ksqs);
    r_source_cmplx_fm_to_field(ksqs->source_file, &(ksqs->v_src[color][0].c[color]), 
			       sizeof(su3_vector)/sizeof(complex), x0, y0, z0, t0);
    /* Do momentum insertion */
    if(status == 0)insert_mom_V(ksqs->v_src[color], mom, x0, y0, z0, t0);

    /* Save a copy as a complex field */
    alloc_ksqs_c_src(ksqs);
    extract_c_from_v(ksqs->c_src, ksqs->v_src[color], color);

    /* Change source type to storage */
    ksqs->type = COMPLEX_FIELD_STORE;

    /* Copy to requested location */
    if(src != NULL)copy_v_field(src, ksqs->v_src[color]);
  }

  else if(source_type == COMPLEX_FIELD_STORE){
    if(ksqs->c_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    
    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }
  
  else if(source_type == CORNER_WALL){
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
    /* Build a corner wall source on time slice t0 */
    FORALLSITES(i,s){
      if( s->t==t0 && s->x % 2 == 0 && s->y % 2 == 0 && s->z % 2 == 0 ){
	(ksqs->c_src+i)->real  = 1.0;
	(ksqs->c_src+i)->imag  = 0.0;
	/* "2" because we invert 2m + 2 Dslash */
      }
    }

    //insert_mom_C(ksqs->c_src, mom, x0, y0, z0, t0);  /* DEBUG ! */

    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }
    
  else if(source_type == EVEN_WALL){
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
    /* Build a corner wall source on time slice t0 */
    FORALLSITES(i,s){
      if( s->t==t0 && (s->x + s->y + s->z + s->t) % 2 == 0 ){
	(ksqs->c_src+i)->real  = 1.0;
	(ksqs->c_src+i)->imag  = 0.0;
	/* "2" because we invert 2m + 2 Dslash */
      }
    }
    
    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }

  else if(source_type == EVENANDODD_WALL){
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
    /* Build an even and odd wall source on time slice t0 */
    FORALLSITES(i,s){
      if( s->t==t0 ){
	(ksqs->c_src+i)->real  = 1.0;
	(ksqs->c_src+i)->imag  = 0.0;
      }
    }

    //insert_mom_C(ksqs->c_src, mom, x0, y0, z0, t0);  /* DEBUG ! */

    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }
  
  else if(source_type == EVENMINUSODD_WALL){
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
    /* Build an even and odd wall source on time slice t0 */
    FORALLSITES(i,s){
      if( s->t==t0 ){
	if( (s->x + s->y + s->z + s->t) % 2 == 0 ){
	  (ksqs->c_src+i)->real  = 1.0;
	  (ksqs->c_src+i)->imag  = 0.0;
	  /* "2" because we invert 2m + 2 Dslash */
	} else {
	  (ksqs->c_src+i)->real  = -1.0;
	  (ksqs->c_src+i)->imag  = 0.0;
	}
      }
    }

    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }
  
  else if(source_type == GAUSSIAN) {
    /* only even sites on time slice t0 */

    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);

    FORALLSITES(i,s) {
      short my_x,my_y,my_z;
      Real rx,ry,rz,radius2;
      Real r0 = ksqs->r0;
      if(s->t != t0 || (s->x + s->y + s->z + s->t) % 2 != 0)
	continue;
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);

      ksqs->c_src[i].real = (Real)exp((double)(- radius2));
    }

    if(src != NULL)
      insert_v_from_c(src, ksqs->c_src, color);
  }

  else if(source_type == POINT) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    /* Save a copy of the source in ksqs->c_src */

    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);

    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      ksqs->c_src[i].real = 1.0;
      if(src != NULL)src[i].c[color].real = 1.0;
    }
  }

  else if(source_type == RANDOM_CORNER_COLOR_WALL){
    reset_ksqs_color(ksqs,color);
    alloc_ksqs_color(ksqs,color);

    /* Build a single random wall source on time slice t0 */

    random_corner_wall(ksqs->v_src[color], t0);
    insert_mom_V(ksqs->v_src[color], mom, x0, y0, z0, t0);

    /* Copy to requested location */
    if(src != NULL){
      copy_v_field(src, ksqs->v_src[color]);
    }
  }

  else if(source_type == RANDOM_COLOR_WALL){
    reset_ksqs_color(ksqs,color);
    alloc_ksqs_color(ksqs,color);

    /* Build a single random wall source on time slice t0 */

    random_wall(ksqs->v_src[color], t0);
    insert_mom_V(ksqs->v_src[color], mom, x0, y0, z0, t0);

    /* Copy to requested location */
    if(src != NULL){
      copy_v_field(src, ksqs->v_src[color]);
    }
  }

  else if(source_type == VECTOR_FIELD_FILE){
    /* Source file contains a list of color vectors in USQCD format,
       i.e., a succession of records with one color vector field per
       record.  Each color vector field is considered to have support
       on the entire lattice, although the record might specify values
       for just one time slice. */
    /* We read and deliver the next color vector in the file. */
#ifdef HAVE_QIO
    QIO_String *recxml;
    
    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    recxml = QIO_string_create();
    /* Load the SU(3) vector from the file.  Zero all components before reading. */
    reset_ksqs_color(ksqs,color);
    alloc_ksqs_color(ksqs,color);
    status = qio_status(
       read_ks_vector_scidac_xml(ksqs->infile, ksqs->v_src[color], 1, recxml) );

    /* Translate the origin as requested */
    rshift[0] = x0; rshift[1] = y0; rshift[2] = z0; rshift[3] = t0;
    shift_su3_vector(ksqs->v_src[color], rshift);

    /* Do momentum insertion */
    if(status == 0)
      insert_mom_V(ksqs->v_src[color], mom, x0, y0, z0, t0);

    /* Verify that the requested color matches the color encoded in
       the record XML */
    if(status == 0)
      status = check_color(recxml, color);

    QIO_string_destroy(recxml);

    /* Copy to requested location */
    if(src != NULL && status == 0)
      copy_v_field(src, ksqs->v_src[color]);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }

  else if(source_type == VECTOR_FIELD_FM_FILE){
    /* Source file is in FNAL format. 
       We read the next record, store it in v_src and copy to src */

    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    reset_ksqs_color(ksqs,color);
    alloc_ksqs_color(ksqs,color);

    /* Read the next source to v_src[color] */
    r_source_ks_fm(ksqs->kssf, ksqs->v_src[color], x0, y0, z0, t0);

    /* Do momentum insertion on all colors */
    insert_mom_V(ksqs->v_src[color], mom, x0, y0, z0, t0);

    /* Deliver the requested color */
    if(src != NULL)copy_v_field(src, ksqs->v_src[color]);
  }
  
  else if(source_type == VECTOR_PROPAGATOR_FILE){
    /* We don't build the source now.  We wait to read it from the file */
  }
  else if(source_type == WAVEFUNCTION_FILE){

    /* Read the file, convert to lattice units, and copy to c_src */
    reset_ksqs(ksqs);
    alloc_ksqs_c_src(ksqs);
    fnal_wavefunction(ksqs->c_src, x0, y0, z0, t0, ksqs->a, ksqs->source_file);

    /* Do momentum insertion */
    if(status == 0)insert_mom_C(ksqs->c_src, mom, x0, y0, z0, t0);

    /* Load to the specified color and timeslice */
    /* zero src to be safe */
    alloc_ksqs_v_src(ksqs);

    insert_v_from_c(ksqs->v_src[color], ksqs->c_src, color);

    /* Copy to requested location */
    if(src != NULL)
      copy_v_field(src, ksqs->v_src[color]);
  }

  else if(source_type == VECTOR_FIELD_STORE){
    if(ksqs->v_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    
    if(src != NULL)
      copy_v_field(src, ksqs->v_src[color]);
    
  }
  
  else {
    node0_printf("%s: Unrecognized source type %d\n",myname, source_type);
    terminate(1);
  }
  return status;

} /* ks_base_source */


/* Build multiple source vectors at once */

int ks_multi_source_field(su3_vector **src, quark_source *ksqs)
{
  char myname[] = "ks_multi_source_field";
  
  int status = 0;
#ifdef HAVE_QIO
  int rshift[4];
#endif
  
  /* Unpack structure */
  int nsource               = ksqs->ncolor;
  int source_type           = ksqs->type;
  int t0                    = ksqs->t0;
#ifdef HAVE_QIO
  int x0                    = ksqs->x0; 
  int y0                    = ksqs->y0; 
  int z0                    = ksqs->z0; 
  int *mom                  = ksqs->mom;
  int file_initialized      = ksqs->file_initialized;
#endif
  int ic;

  /* zero src to be safe */
  if(src != NULL){
    for(ic = 0; ic < nsource; ic++){
      if(src[ic] != NULL)
	clear_v_field( src[ic] );
    }
  }
  
  /* Unless we are taking the source from storage, clear it before
     reconstructing it */
  if(source_type != COMPLEX_FIELD_STORE &&
     source_type != VECTOR_FIELD_STORE)
    reset_ksqs(ksqs);

  if(source_type == VECTOR_FIELD_FILE){
#ifdef HAVE_QIO
    QIO_String *recxml;

    if(file_initialized == 0)
      r_open_ks_source(ksqs);

    recxml = QIO_string_create();
    /* Load the SU(3) vector from the file */
    alloc_ksqs_v_src(ksqs);
    rshift[0] = x0; rshift[1] = y0; rshift[2] = z0; rshift[3] = t0;
    for(ic = 0; ic < nsource; ic++){
      status = qio_status(
        read_ks_vector_scidac_xml(ksqs->infile, ksqs->v_src[ic], 1, recxml) );
      /* Apply the requested translation */
      shift_su3_vector(ksqs->v_src[ic], rshift);
      /* Do momentum insertion */
      if(status == 0)
	insert_mom_V(ksqs->v_src[ic], mom, x0, y0, z0, t0);
      /* Verify that the requested color matches the color encoded in
	 the record XML */
      if(status == 0)status = check_color(recxml, ic);

      /* Copy to requested location */
      if(src[ic] != NULL && status == 0)
	copy_v_field(src[ic], ksqs->v_src[ic]);
    }

    /* Change source type to storage */
    ksqs->type = VECTOR_FIELD_STORE;

    QIO_string_destroy(recxml);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }

  else if(source_type == VECTOR_FIELD_STORE){
    if(ksqs->v_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }

    if(src != NULL){
      for(ic = 0; ic < nsource; ic++){
	if(src[ic] != NULL)
	  copy_v_field(src[ic], ksqs->v_src[ic]);
      }
    }
  }
  else if(source_type == RANDOM_CORNER_COLOR_WALL){

    /* Build random wall sources on time slice t0 */
    int ic;

    for(ic=0;ic<nsource;ic++){
      random_corner_wall(ksqs->v_src[ic], t0);
      insert_mom_V(ksqs->v_src[ic], mom, x0, y0, z0, t0);
    }

    /* Copy to requested location */
    if(src != NULL){
      for(ic = 0; ic < nsource; ic++){
	if(src[ic] != NULL)
	  copy_v_field(src[ic], ksqs->v_src[ic]);
      }
    }

    /* Change source type to storage */
    ksqs->type = VECTOR_FIELD_STORE;

  }
  else if(source_type == RANDOM_COLOR_WALL){

    /* Build random wall sources on time slice t0 */
    int ic;

    for(ic=0;ic<nsource;ic++){
      random_wall(ksqs->v_src[ic], t0);
      insert_mom_V(ksqs->v_src[ic], mom, x0, y0, z0, t0);
    }

    /* Copy to requested location */
    if(src != NULL){
      for(ic = 0; ic < nsource; ic++){
	if(src[ic] != NULL)
	  copy_v_field(src[ic], ksqs->v_src[ic]);
      }
    }

    /* Change source type to storage */
    ksqs->type = VECTOR_FIELD_STORE;

  }
  else {
    node0_printf("%s: Unrecognized multi source type %d\n",myname, source_type);
    terminate(1);
  }
  return status;

} /* ks_multi_source_field */

void r_close_ks_source(quark_source *ksqs){

  char myname[] = "r_close_ks_source";

  if(ksqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
  if(ksqs->type == VECTOR_FIELD_FM_FILE){
    r_source_ks_fm_f(ksqs->kssf);
    ksqs->file_initialized = 0;
    ksqs->kssf = NULL;
  }
#ifdef HAVE_QIO
  else if(ksqs->type == COMPLEX_FIELD_FILE){
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


void w_open_ks_source(quark_source *ksqs, char *fileinfo){

  char myname[] = "w_open_ks_source";

#ifdef HAVE_QIO
  int volfmt, serpar;
  char *source_file = ksqs->source_file;

  interpret_usqcd_ks_save_flag(&volfmt, &serpar, ksqs->flag);

  if(ksqs->file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }

  if(ksqs->type == COMPLEX_FIELD_FILE || ksqs->type == COMPLEX_FIELD_STORE ||
     ksqs->type == COMPLEX_FIELD_FM_FILE){
    ksqs->outfile = w_open_complex_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    ksqs->file_initialized = 1;
  }
  if(ksqs->type == VECTOR_FIELD_FILE || ksqs->type == VECTOR_FIELD_STORE ||
     ksqs->type == VECTOR_FIELD_FM_FILE){
    ksqs->outfile = w_open_ks_vector_scidac_file(source_file, fileinfo,
						 volfmt, serpar);
    ksqs->file_initialized = 1;
  }
  else
    node0_printf("%s: Can't open a file for source type %d\n",myname, ksqs->type);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
}

/* For now we support writing only a color vector source field */

int ks_source_write(su3_vector *src, quark_source *ksqs)
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

  if(source_type != VECTOR_FIELD_FILE && source_type != VECTOR_FIELD_STORE){
    node0_printf("%s: Don't know how to write source type %d\n", myname, source_type);
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
#ifdef IOTIME
  node0_printf("Time to save source color %d = %e\n", color,dtime);
#endif
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
}

void w_close_ks_source(quark_source *ksqs){

  char myname[] = "w_close_ks_source";

  if(ksqs->file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }
#ifdef HAVE_QIO
  if(ksqs->type == COMPLEX_FIELD_FILE || ksqs->type == COMPLEX_FIELD_STORE ||
     ksqs->type == COMPLEX_FIELD_FM_FILE){
    w_close_complex_scidac_file(ksqs->outfile);
    ksqs->file_initialized = 0;
    ksqs->outfile = NULL;
  }
  else if(ksqs->type == VECTOR_FIELD_FILE || ksqs->type == VECTOR_FIELD_STORE){
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

int ks_source_site(field_offset src, quark_source *ksqs)
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
  
  status = ks_base_source(t_src, ksqs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(su3_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);

  return status;

} /* ks_source_site */


/* 
   Initialize a sink.  Result is a complex field.
 */

void ks_sink_field(complex *snk, quark_source *ksqs)
{
  int i;
  site *s;
#ifdef HAVE_QIO
  int rshift[4];
#endif
  
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
    
    /* Translate the origin as requested (but no translation in t for sink ) */
    rshift[0] = x0; rshift[1] = y0; rshift[2] = z0; rshift[3] = 0;
    shift_complex(snk, rshift);
    
    /* Do momentum insertion */
    insert_mom_C(snk, mom, x0, y0, z0, t0);
#else
    node0_printf("ks_sink_field: QIO compilation required for this source\n");
    terminate(1);
#endif
  }
  else if(sink_type == COMPLEX_FIELD_FM_FILE){
    r_source_cmplx_fm_to_field(sink_file, snk, 1, x0, y0, z0, t0);
    
    /* Do momentum insertion */
    insert_mom_C(snk, mom, x0, y0, z0, t0);
  }
  else if(sink_type == EVENANDODD_WALL){
    /* Build an even and odd wall source on all time slices */
    FORALLSITES(i,s){
      snk[i].real  = 1.0;
      snk[i].imag  = 0.0;
    }
  }
  
  else if(sink_type == GAUSSIAN) {
    /* only even sites on time slice t0 */

    FORALLSITES(i,s) {
      short my_x,my_y,my_z;
      Real rx,ry,rz,radius2;
      Real r0 = ksqs->r0;
      if(s->t != t0 || (s->x + s->y + s->z + s->t) % 2 != 0)
	continue;
      
      my_x = ((s->x)-x0+nx) % nx;
      rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
      my_y = ((s->y)-y0+ny) % ny;
      ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
      my_z = ((s->z)-z0+nz) % nz;
      rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
      
      radius2 = rx*rx + ry*ry + rz*rz;
      radius2 /= (r0*r0);

      snk[i].real = (Real)exp((double)(- radius2));
    }
  }

  else {
    node0_printf("ks_sink_field: Unrecognized sink type %d\n",sink_type);
    terminate(1);
  }
  
} /* ks_sink_field */


/* snk must be type complex */
void ks_sink_site(field_offset snk, quark_source *ksqs)
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

/*--------------------------------------------------------------------*/
static void sink_smear_prop(ks_prop_field *ksprop, quark_source *ksqs){
  
  int color;
  int ci,cf;
  int i;
  site *s;
  complex *chi_cs, z;
  Real x;
  su3_vector *v;
  double dtime = start_timing();
  int key[4] = {1,1,1,0};  /* 3D Fourier transform */
  
  /* Set up Fourier transform for smearing */
  setup_restrict_fourier(key, NULL);

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark_propagator (in place) */
  for(color = 0; color < ksprop->nc; color++){
    v = ksprop->v[color];
    /* v points into ksprop, so ksprop is changed here */
    restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			   FORWARDS);
  }

  print_timing(dtime,"FFT");

  dtime = start_timing();

  /* Build sink smearing wave function as a complex field repeated on
     each time slice */
  ks_sink_field(chi_cs, ksqs);

  /* Normalize (for FFT) */
  x = 1./(((Real)nx)*((Real)ny)*(Real)nz);
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i],x,chi_cs[i]);
  }
  
  print_timing(dtime,"build sink field");

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function */
  for(ci=0;ci<3;ci++){
    FORALLSITES(i,s)
      for(cf=0;cf<ksprop->nc;cf++){
	z = ksprop->v[ci][i].c[cf];
	CMUL(z, chi_cs[i], ksprop->v[ci][i].c[cf]);
      }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark_propagator (in place) */
  for(color = 0; color < ksprop->nc; color++){
    v = ksprop->v[color];
    /* v points into ksprop, so ksprop is changed here */
    restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			   BACKWARDS);
  }
  print_timing(dtime,"FFT");
  cleanup_restrict_fourier();
  free(chi_cs);
}  


/* 
   Apply a source operation
 */

void ks_source_op(quark_source *mod_ksqs, quark_source_op *ksqs_op,
		  quark_source *ksqs )
{
  /* The source operator is specified by mod_ksqs.  The parent source
     is in ksqs. The result goes in mod_ksqs. */
  
  char myname[] = "ks_source_op";
  int i;

  /* Unpack structure. */
  int op_type     = ksqs_op->type;
  int *mom        = ksqs_op->mom;
  int nsource     = ksqs->ncolor;
  int x0          = ksqs->x0; 
  int y0          = ksqs->y0; 
  int z0          = ksqs->z0; 
  int t0          = ksqs->t0;

  /* Copy some attributes */
  mod_ksqs->ncolor           = ksqs->ncolor;
  mod_ksqs->x0               = ksqs->x0; 
  mod_ksqs->y0               = ksqs->y0; 
  mod_ksqs->z0               = ksqs->z0; 
  mod_ksqs->t0               = ksqs->t0;
  mod_ksqs->file_initialized = 0;
  strcpy(mod_ksqs->descrp, ksqs->descrp);
  strcpy(mod_ksqs->label,  ksqs->label);
  /* Append the operator info to the description and label if the operator is nontrivial */
  if(op_type != POINT){
    strncat(mod_ksqs->descrp, "/", MAXDESCRP-strlen(mod_ksqs->descrp)-1);
    strncat(mod_ksqs->descrp, ksqs_op->descrp, MAXDESCRP-strlen(mod_ksqs->descrp)-1);
    strncat(mod_ksqs->label,  ksqs_op->label, MAXSRCLABEL-strlen(mod_ksqs->label)-1);
  }

  /* Identity source operator */
  if(op_type == POINT) {
    if(ksqs->c_src != NULL){
      alloc_ksqs_c_src(mod_ksqs);
      copy_c_field(mod_ksqs->c_src, ksqs->c_src);
      mod_ksqs->type = COMPLEX_FIELD_STORE;
    }
    else if(ksqs->v_src != NULL){
      alloc_ksqs_v_src(mod_ksqs);
      for(i = 0; i < nsource; i++)
	copy_v_field(mod_ksqs->v_src[i], ksqs->v_src[i]);
      mod_ksqs->type = VECTOR_FIELD_STORE;
    }
    else {
      printf("%s(%d): Base source should not be NULL\n",myname, this_node);
      terminate(1);
    }
  }
  
  /* (couples to pion5, pioni5, pioni, pions, rhoi, rhos) */
  else if(op_type == FUNNYWALL1) {
    su3_vector *tvec = create_v_field();

    /* Convert complex field source to color vector field if necessary */
    if(ksqs->v_src == NULL){
      alloc_ksqs_v_src(ksqs);
      if(nsource != 3){
	node0_printf("%s: Requires three colors.\n",myname);
	terminate(1);
      }
      for(i = 0; i < nsource; i++)
	insert_v_from_c(ksqs->v_src[i], ksqs->c_src, i);
      free(ksqs->c_src); ksqs->c_src = NULL;
      ksqs->type = VECTOR_FIELD_STORE;
    }
    
    alloc_ksqs_v_src(mod_ksqs);
    mod_ksqs->type = VECTOR_FIELD_STORE;

    for(i = 0; i < nsource; i++){
      mult_pion5_field( ksqs->v_src[i], tvec );
      copy_v_field( mod_ksqs->v_src[i], tvec );
      mult_pioni_field( ZUP, ksqs->v_src[i], tvec );
      add_v_fields( mod_ksqs->v_src[i], tvec, mod_ksqs->v_src[i] );
      mult_rhoi_field( ZUP, ksqs->v_src[i], tvec );
      add_v_fields( mod_ksqs->v_src[i], tvec, mod_ksqs->v_src[i] );
    }
    
    destroy_v_field(tvec);
  }

  /* (couples to pion05, pionij, pioni0, pion0, rhoi0, rho0) */
  else if(op_type == FUNNYWALL2) {
    su3_vector *tvec = create_v_field();

    /* Convert complex field source to color vector field if necessary */
    if(ksqs->v_src == NULL){
      alloc_ksqs_v_src(ksqs);
      if(nsource != 3){
	node0_printf("%s: Requires three colors.\n",myname);
	terminate(1);
      }
      for(i = 0; i < nsource; i++)
	insert_v_from_c(ksqs->v_src[i], ksqs->c_src, i);
      free(ksqs->c_src); ksqs->c_src = NULL;
      ksqs->type = VECTOR_FIELD_STORE;
    }
    
    alloc_ksqs_v_src(mod_ksqs);
    mod_ksqs->type = VECTOR_FIELD_STORE;

    for(i = 0; i < nsource; i++){
      mult_pion05_field( ksqs->v_src[i], tvec );
      copy_v_field( mod_ksqs->v_src[i], tvec );
      mult_pioni0_field( ZUP, ksqs->v_src[i], tvec );
      add_v_fields( mod_ksqs->v_src[i], tvec, mod_ksqs->v_src[i] );
      mult_rhoi0_field( ZUP, ksqs->v_src[i], tvec );
      add_v_fields( mod_ksqs->v_src[i], tvec, mod_ksqs->v_src[i] );
    }
    
    destroy_v_field(tvec);
  }

  else{
    node0_printf("%s: Unknown source_op type %d\n",myname,op_type);
    terminate(1);
  }

  /* Do momentum insertion */
  if(mod_ksqs->c_src != NULL)
    insert_mom_C(mod_ksqs->c_src, mom, x0, y0, z0, t0);
  else if(mod_ksqs->v_src != NULL){
    for(i = 0; i < nsource; i++)
      insert_mom_V(mod_ksqs->v_src[i], mom, x0, y0, z0, t0);
  }
  else{
    printf("%s(%d): Modified source should not be NULL\n",myname, this_node);
    terminate(1);
  }
    
  
} /* ks_source_op */


/* 
   Apply a sink operation on the given propagator.
 */

ks_prop_field *ks_sink_op(quark_source *ksqs, ks_prop_field *src )
{
  
  /* Unpack structure. */
  int sink_type     = ksqs->type;
  ks_prop_field *dst = NULL;
  char myname[] = "ks_sink_op";

  /* Identity sink operator */
  if(sink_type == POINT) {
    dst = create_ksp_field_copy(src);
  }
  /* Various wave-function-based smearings requiring convolution */
  else if(sink_type == GAUSSIAN || 
	  sink_type == COMPLEX_FIELD_FM_FILE ||
	  sink_type == COMPLEX_FIELD_FILE ||
	  sink_type == EVENANDODD_WALL ||
	  sink_type == WAVEFUNCTION_FILE) {
    dst = create_ksp_field_copy(src);
    sink_smear_prop(dst, ksqs);
  }
  
  else if(sink_type == FAT_COVARIANT_GAUSSIAN){

    /* The following operations are performed as in the source,
       but in reverse order */
    
    /* dst <- src */
    dst = create_ksp_field_copy(src);
    
    gauss_smear_ks_prop_field(dst, ape_links, ksqs->r0, ksqs->iters, ALL_T_SLICES);
  }
  else{
    node0_printf("%s: Unknown sink type %d\n",myname,sink_type);
    terminate(1);
  }

  return dst;
} /* ks_sink_op */


int ask_quark_source( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_quark_source";

  if (prompt==1){
    printf("enter ");
    printf("'complex_field', 'complex_field_fm', 'corner_wall',");
    printf("'even_wall', 'evenandodd_wall', 'evenminusodd_wall', ");
    printf("'fat_covariant_gaussian', 'gaussian'");
    printf("'point', 'random_corner_vector_wall', 'random_vector_wall',");
    printf("'vector_field','vector_field_fm', 'vector_propagator_file',");
    printf("'wavefunction', for source type\n");
  }

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  /* Source file is a complex field over the entire lattice: SciDAC
     format */
  if(strcmp("complex_field",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FILE;
    strcpy(descrp,"complex_field");
  }
  else if(strcmp("complex_field_fm",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FM_FILE;
    strcpy(descrp,"complex_field_fm");
  }
  else if(strcmp("corner_wall",savebuf) == 0 ) {
    *source_type = CORNER_WALL;
    strcpy(descrp,"CORNER");
  }
  else if(strcmp("even_wall",savebuf) == 0 ) {
    *source_type = EVEN_WALL;
    strcpy(descrp,"even_wall");
  }
  else if(strcmp("evenandodd_wall",savebuf) == 0 ) {
    *source_type = EVENANDODD_WALL;
    strcpy(descrp,"even_and_odd_wall");
  }
  else if(strcmp("evenminusodd_wall",savebuf) == 0 ) {
    *source_type = EVENMINUSODD_WALL;
    strcpy(descrp,"even_minus_odd_wall");
  }
  else if(strcmp("fat_covariant_gaussian",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN;
    strcpy(descrp,"fat_covariant_gaussian");
  }
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  /* Source file is a series of random color vector fields on
     corners of the hypercubes: SciDAC format */
  else if(strcmp("random_corner_vector_wall",savebuf) == 0 ) {
    *source_type = RANDOM_CORNER_COLOR_WALL;
    strcpy(descrp,"RANDOM_CORNER_WALL");
  }
  /* Source file is a series of random color vector fields over
     the entire lattice: SciDAC format */
  else if(strcmp("random_vector_wall",savebuf) == 0 ) {
    *source_type = RANDOM_COLOR_WALL;
    strcpy(descrp,"RANDOM_WALL");
  }
  /* Source file is an SU(3) vector field over the entire lattice: SciDAC
     format */
  else if(strcmp("vector_field",savebuf) == 0 ) {
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
  }
  /* Source file is an SU(3) vector field over a single time slice: FNAL
     format */
  else if(strcmp("vector_field_fm",savebuf) == 0 ) {
    *source_type = VECTOR_FIELD_FM_FILE;
    strcpy(descrp,"vector_field_fm");
  }
  else if(strcmp("vector_propagator_file",savebuf) == 0 ) {
    *source_type = VECTOR_PROPAGATOR_FILE;
    strcpy(descrp,"vector_propagator_file");
  }
  else if(strcmp("wavefunction",savebuf) == 0 ){
    *source_type = WAVEFUNCTION_FILE;
    strcpy(descrp,"wavefunction");
  }
  else{
    printf("%s: ERROR IN INPUT: KS source command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_source */

int ask_quark_source_op( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_quark_source_op";

  if (prompt==1){
    printf("enter ");
    printf("'identity', ");
    printf("'funnywall1', ");
    printf("'funnywall2', ");
    printf(", for source type\n");
  }

  savebuf = get_next_tag(fp, "quark source op command", myname);
  if (savebuf == NULL)return 1;

  /* Identity operator */
  if(strcmp("identity",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  /* Funnywall1 operator (couples to pion5, pioni5, pioni, pions, rhoi, rhos) */
  else if(strcmp("funnywall1",savebuf) == 0 ){
    *source_type = FUNNYWALL1;
    strcpy(descrp,"FUNNYWALL1");
  }
  /* Funnywall2 operator (couples to pion05, pionij, pioni0, pion0, rhoi0, rho0) */
  else if(strcmp("funnywall2",savebuf) == 0 ){
    *source_type = FUNNYWALL2;
    strcpy(descrp,"FUNNYWALL2");
  }
  else{
    printf("%s: ERROR IN INPUT: KS source_op command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;

} /* ask_quark_source_op */

#define IF_OK if(status==0)

int ask_output_quark_source_file( FILE *fp, int prompt, 
				     int *flag, int *source_type,
				     int *t0, char *descrp, char *filename)
{
  char *savebuf;
  int status = 0;
  char myname[] = "ask_output_quark_source_file";

  if (prompt==1){
    printf("'forget_ks_source' or ");
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

  if(strcmp("forget_ks_source",savebuf) == 0 ) {
    *flag=FORGET;
    *source_type = UNKNOWN;
    strcpy(descrp,"");
  }
  else if(strcmp("save_serial_scidac_ks_source",savebuf) == 0 ) {
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
    if(prompt==1)printf("enter filename\n");
    if(scanf("%s",filename) != 1){
      printf("\n%s(%d): ERROR IN INPUT: Can't read filename\n",
	     myname, this_node); 
      status++;
    }
    else
      printf("%s\n",filename);
  } else {
    printf("\n");
  }

  /* Get time slice*/
  if(t0 != NULL && *flag != FORGET)
    IF_OK status += get_i(stdin, prompt, "t0", t0);

  return status;
} /* ask_output_quark_source_file */

/* Get the additional input parameters needed to specify the source */
int get_quark_source(FILE *fp, int prompt, quark_source *ksqs){
  
  int  source_type;
  int  source_loc[4] = { 0,0,0,0 };
  char source_file[MAXFILENAME] = "";
  char source_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_quark_source(fp, prompt,&source_type,
				      ksqs->descrp);
  IF_OK ksqs->type  = source_type;
  
  IF_OK {
    if ( source_type == POINT ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
    }
    else if ( source_type == CORNER_WALL ||
	      source_type == EVEN_WALL ||
	      source_type == EVENANDODD_WALL ||
	      source_type == EVENMINUSODD_WALL ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      //IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3); /* DEBUG! */
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
    else if ( source_type == RANDOM_CORNER_COLOR_WALL ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_i(stdin, prompt, "nsource", &(ksqs->ncolor));
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == RANDOM_COLOR_WALL ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_i(stdin, prompt, "nsource", &(ksqs->ncolor));
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == VECTOR_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_i(stdin, prompt, "nsource", &(ksqs->ncolor));
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == VECTOR_FIELD_FM_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
      IF_OK status += get_i(stdin, prompt, "nsource", &(ksqs->ncolor));
      IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
    }
    else if ( source_type == VECTOR_PROPAGATOR_FILE ){
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
} /* get_quark_source */

/* Get the additional input parameters needed to specify the source */
int get_quark_source_op(FILE *fp, int prompt, quark_source_op *ksqs){
  
  int  source_type;
  char source_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_quark_source_op(fp, prompt, &source_type,
					 ksqs->descrp);
  IF_OK ksqs->type  = source_type;
  
  IF_OK {
    if ( source_type == POINT ){
    }
    else if (source_type == FUNNYWALL1 ){
    }
    else if (source_type == FUNNYWALL2 ){
    }
    else {
      printf("Source_op type not supported in this application\n");
      status++;
    }
  }

  IF_OK status += get_vi(stdin, prompt, "momentum", ksqs->mom, 3);
  IF_OK status += get_s(stdin, prompt, "source_label", source_label);

  strncpy(ksqs->label,source_label,MAXSRCLABEL);
  
  return status;
} /* get_quark_source_op */

int get_ks_quark_sink(FILE *fp, int prompt, quark_source *ksqs){

  int  sink_type;
  int  sink_loc[3] = { 0,0,0 };
  char sink_file[MAXFILENAME] = "";
  char sink_label[MAXSRCLABEL];
  int  status = 0;

  /* Get antiquark source type */
  IF_OK status += ask_quark_source(fp, prompt,&sink_type,
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
    else if ( sink_type == EVENANDODD_WALL ){
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


