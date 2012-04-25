/********************** quark_source.c ******************************/
/* MIMD version 7 */

/* Initialize a source for the inverter */
/* Modify a source by applying source operators */
/* Apply sink operators */
/* Interpret source and sink parameters */

/* Choices of base sources

   Complex field options (diagonal in source colors and spins):

   point                      Kronecker delta at a specified source point 
                              (same for all colors and spins.)
   gaussian                   Gaussian distribution with center at a 
                              specified point for all colors and spins.
   wavefunction_file          File containing an ASCII wave function, 
                              to be converted via interpolation into a 
			      "complex field" based on a specified 
			      lattice spacing.
   complex_field_file         File is a single complex field in QIO format.
                              The same source is used for each source color 
			      and spin.
   complex_field_fm_file      File is a single complex field in FNAL format.
                              The same source is used for each source 
			      color and spin.
   Color vector field options:
   
   random_color_wall          Same as above, but on all sites on a time slice.
   vector_field_file          List of color fields, replicated for each 
                              of four source spins.
   even_wall                  +1 on even sites.  zero elsewhere.
   evenandodd_wall            +1 on all sites.
   evenminusodd_wall          +1 on even sites, -1 on odd sites.

   Dirac field options:

   dirac_field_file           File contains a series of Dirac fields in 
                              QIO format. For a color-spin interpretation, 
			      the color varies most rapidly. For a source 
			      list interpretation, the sources are read 
			      sequentially.
   dirac_field_fm_file        File contains a series of Dirac fields in 
                              FNAL format.
                              Color-spin interpretation only, 
			      the color varies most rapidly.

   See quark_source_sink_ops.c for operations on sources and sinks.

   Immutable source options:

   propagator_file            Take source directly from QIO propagator file
                              (i.e. one that contains both sources and 
			      solutions.  No modifications permitted.
*/

/*
  Entry points

     init_qs
     alloc_cached_c_source
     alloc_cached_v_source
     alloc_cached_wv_source
     ask_starting_quark_source
     get_cached_c_source
     get_cached_v_source
     get_cached_wv_source
     clear_qs
     convert_ksource_to_color
     convert_ksource_to_spin
     corner_wall
     even_wall
     even_and_odd_wall
     even_minus_odd_wall
     gaussian_source
     point_source
     random_color_wall
     subset_mask_v
     subset_mask_wv
     get_complex_source
     v_base_source
     wv_base_source
     v_source_field
     v_source_site
     wv_source_field
     get_v_quark_source
     get_wv_quark_source
     print_source_info

 */

#include "generic_includes.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846;
#endif

/*-------------------------------------------------------------*/
/* Quark source utilities */
/*-------------------------------------------------------------*/

void init_qs(quark_source *qs){
  qs->type             = UNKNOWN;
  qs->subset           = FULL;
  qs->scale_fact       = 1.;
  qs->color            = -1;  /* Counter will be preincremented */
  qs->ncolor           =  3;
  qs->ksource          = -1;  /* Counter will be preincremented */
  qs->nsource          = 12;
  qs->c_src            = NULL;
  qs->v_src            = NULL;
  qs->wv_src           = NULL;
  qs->source_file[0]   = '\0';
  qs->sourceflag       = RELOAD_SERIAL;
  qs->source_file_initialized = 0;
  qs->save_file_initialized = 0;
#ifdef HAVE_QIO
  qs->infile           = NULL;
  qs->outfile          = NULL;
#endif
  qs->kssf             = NULL;
  qs->x0               = 0;
  qs->y0               = 0;
  qs->z0               = 0;
  qs->t0               = 0;
  qs->r0               = 0.;
  qs->descrp[0]        = '\0';
  qs->mom[0]           = 0;
  qs->mom[1]           = 0;
  qs->mom[2]           = 0;
  qs->op               = NULL;
}

/* Allocate cached source members */

void alloc_cached_c_source(quark_source *qs)
{
  if(qs->c_src == NULL)
    qs->c_src = create_c_field();
  clear_c_field(qs->c_src);
}

void alloc_cached_v_source(quark_source *qs)
{
  if(qs->v_src == NULL)
    qs->v_src = create_v_field();
  clear_v_field(qs->v_src);
}

void alloc_cached_wv_source(quark_source *qs)
{
  if(qs->wv_src == NULL)
    qs->wv_src = create_wv_field();
  clear_wv_field(qs->wv_src);
}

/* Access cached source members */

complex *get_cached_c_source(quark_source *qs){
  return qs->c_src;
}

su3_vector *get_cached_v_source(quark_source *qs){
  return qs->v_src;
}

wilson_vector *get_cached_wv_source(quark_source *qs){
  return qs->wv_src;
}

/* Free all allocated space */

static void reset_qs(quark_source *qs){
  if(qs->c_src != NULL){
    free(qs->c_src); 
    qs->c_src = NULL;
  }
  if(qs->v_src != NULL){
    free(qs->v_src);
    qs->v_src = NULL;
  }
  if(qs->wv_src != NULL){
    free(qs->wv_src);
    qs->wv_src = NULL;
  }
}

/* Close files and free allocated space */

void clear_qs(quark_source *qs){
  reset_qs(qs);
#ifdef HAVE_QIO
  if(qs->infile != NULL){
    if(qs->type == COMPLEX_FIELD_FILE)
      r_close_complex_scidac_file(qs->infile);
#ifdef HAVE_KS
    else if(qs->type == VECTOR_FIELD_FILE)
      r_close_ks_vector_scidac_file(qs->infile);
#endif
#ifdef HAVE_DIRAC
    else if(qs->type == DIRAC_FIELD_FILE)
      r_close_w_vector_scidac_file(qs->infile);
#endif
  }
#endif

#ifdef HAVE_KS
  /* For VECTOR_FIELD_FM_FILE */
  if(qs->kssf != NULL){
    r_source_ks_fm_f(qs->kssf);
    qs->kssf = NULL;
  }
#endif
  qs->source_file_initialized = 0;
}

/*--------------------------------------------------------------------*/
/* Color and spin translation */
/*--------------------------------------------------------------------*/

/* The sources are listed in a flat array indexed by ksource.  When color and
   spin are relevant, we let the spin index vary most rapidly. */

int convert_ksource_to_color(int ksource){
  return ksource/4;
}

int convert_ksource_to_spin(int ksource){
  return ksource % 4;
}

static int convert_ncolor_to_nsource(int ncolor){
  return ncolor * 4;
}


/*--------------------------------------------------------------------*/
/* Momentum insertion */
/*--------------------------------------------------------------------*/

static void insert_c_mom(complex *c, int mom[3], 
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

#ifdef HAVE_KS

static void insert_v_mom(su3_vector *v, int mom[3], 
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

#endif

#ifdef HAVE_DIRAC

static void insert_wv_mom(wilson_vector *wv, int mom[3], 
			 int x0, int y0, int z0, int t0){
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
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	z = wv[i].d[d].c[c];
	CMUL(z,y,wv[i].d[d].c[c]);
      }
  }
}

#endif

/*--------------------------------------------------------------------*/
/* Build complex-field sources */
/*--------------------------------------------------------------------*/

static void corner_wall(complex *c, int t0){
  int i;
  site *s;

  /* Build a corner wall source on time slice t0 */
  FORALLSITES(i,s){
    if( (s->t==t0 || t0 == ALL_T_SLICES) &&
	s->x % 2 == 0 && s->y % 2 == 0 && s->z % 2 == 0 ){
      c[i].real  = 1.0;     c[i].imag  = 0.0;
    }
  }
}
    
static void even_wall(complex *c, int t0){
  int i;
  site *s;

  FORALLSITES(i,s){
    if( (s->t==t0 || t0 == ALL_T_SLICES) && 
	(s->x + s->y + s->z + s->t) % 2 == 0 ){
      c[i].real  = 1.0;     c[i].imag  = 0.0;
    }
  }
}
    
void even_and_odd_wall(complex *c, int t0){
  int i;
  site *s;

  FORALLSITES(i,s){
    if( s->t==t0 || t0 == ALL_T_SLICES){
      c[i].real  = 1.0;     c[i].imag  = 0.0;
    }
  }
}
    
static void even_minus_odd_wall(complex *c, int t0){
  int i;
  site *s;

  FORALLSITES(i,s){
    if( s->t==t0 || t0 == ALL_T_SLICES){
      if( (s->x + s->y + s->z + s->t) % 2 == 0 ){
	c[i].real  = 1.0;   c[i].imag  = 0.0;
      } else {
	c[i].real  = -1.0;  c[i].imag  = 0.0;
      }
    }
  }
}
    
void gaussian_source(complex *src, Real r0, 
		     int x0, int y0, int z0, int t0)
{
  short my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  int i;
  site *s;

  /* Gaussian trial source centered on  x0,y0,z0,t0 */
  
  FORALLSITES(i,s) {
    if(t0 != ALL_T_SLICES && s->t != t0)continue;
    my_x = ((s->x)-x0+nx) % nx;
    rx = (my_x < (nx-my_x)) ? (Real) my_x : (Real) (my_x-nx);
    my_y = ((s->y)-y0+ny) % ny;
    ry = (my_y < (ny-my_y)) ? (Real) my_y : (Real) (my_y-ny);
    my_z = ((s->z)-z0+nz) % nz;
    rz = (my_z < (nz-my_z)) ? (Real) my_z : (Real) (my_z-nz);
    
    radius2 = rx*rx + ry*ry + rz*rz;
    radius2 /= (r0*r0);
    
    src[i].real = (Real)exp((double)(- radius2));
  }
}

static void point_source(complex *src, int x0, int y0, int z0, int t0){
  int i;
  
  /* load 1.0 into source at cooordinates given by source_coord */
  /* initialize src to be a delta function at point x0,y0,z0,t0 */
  /* Save a copy of the source in qs->c_src */
  
  if(node_number(x0,y0,z0,t0) == mynode()){
    i = node_index(x0,y0,z0,t0);
    src[i].real = 1.0;
  }
}

#ifdef HAVE_KS

/*--------------------------------------------------------------------*/
/* Build color-vector-field sources */
/*--------------------------------------------------------------------*/
/* Generate a random color vector, each component of which is orthogonal
   to the other and has variance 1/3 */

static void random_color_wall(su3_vector *v, int t0){
  int i,jc;
  site *s;
  Real x;

  FORALLSITES(i,s){
    if( s->t==t0 || t0 == ALL_T_SLICES){
      for(jc=0;jc<3;jc++){
	v[i].c[jc].real = gaussian_rand_no(&(s->site_prn));
	v[i].c[jc].imag = gaussian_rand_no(&(s->site_prn));
      }
      x = 1.0/sqrt( magsq_su3vec( v+i ) );
      scalar_mult_su3_vector( v+i, x, v+i );
    }
  }
}
    
#if 0
static void random_corner_wall(su3_vector *v, int t0){
  int i,jc;
  site *s;
  Real x;

  FORALLSITES(i,s){
    if( (s->t==t0 || t0 == ALL_T_SLICES) && 
	s->x%2==0 && s->y%2==0 && s->z%2==0 ){
      for(jc=0;jc<3;jc++){
	v[i].c[jc].real = gaussian_rand_no(&(s->site_prn));
	v[i].c[jc].imag = gaussian_rand_no(&(s->site_prn));
      }
      x = 1.0/sqrt( magsq_su3vec( v+i ) );
      scalar_mult_su3_vector( v+i, x, v+i );
    }
  }
}
#endif

#endif
    
/*--------------------------------------------------------------------*/
/* Subset mask operation                                              */
/*--------------------------------------------------------------------*/

/* Apply subset mask to color vector field */

void subset_mask_v(su3_vector *src, int subset, int t0)
{
  int i;
  site *s;

  if(subset == FULL)
    return;

  else if(subset == HYPERCUBE){
    FORALLSITES(i,s) {
      if(t0 != ALL_T_SLICES && s->t != t0)continue;
      if(s->x % 2 != 0 || s->y %2 != 0 || s->z %2 != 0){
	clearvec(src+i);
      }
    }
  }

  else{
    node0_printf("subset_mask_wv: Unrecognized subset type %d\n", subset);
    terminate(1);
  }
}

/* Apply subset mask to Dirac field */

void subset_mask_wv(wilson_vector *src, int subset, int t0)
{
  int i;
  site *s;

  if(subset == FULL)
    return;

  else if(subset == HYPERCUBE){
    FORALLSITES(i,s) {
      if(t0 != ALL_T_SLICES && s->t != t0)continue;
      if(s->x % 2 != 0 || s->y %2 != 0 || s->z %2 != 0){
	clear_wvec(src+i);
      }
    }
  }

  else{
    node0_printf("subset_mask_wv: Unrecognized subset type %d\n", subset);
    terminate(1);
  }
}

int is_complex_source(int source_type){

  return
    source_type == COMPLEX_FIELD_FILE ||
    source_type == COMPLEX_FIELD_FM_FILE ||
    source_type == COMPLEX_FIELD_STORE ||
    source_type == CORNER_WALL ||
    source_type == EVEN_WALL ||
    source_type == EVENANDODD_WALL ||
    source_type == EVENMINUSODD_WALL ||
    source_type == GAUSSIAN ||
    source_type == POINT ||
    source_type == WAVEFUNCTION_FILE;
}

/* Construct or read a complex source */
/* Return 1 if this is actually a complex field source and 0 if not */
int get_complex_source(quark_source *qs){

  int source_type           = qs->type;
  Real a                    = qs->a;
  int x0                    = qs->x0; 
  int y0                    = qs->y0; 
  int z0                    = qs->z0; 
  int t0                    = qs->t0;
  Real r0                   = qs->r0;
  char *source_file         = qs->source_file;
  int *mom                  = qs->mom;
  int status = 0;
#ifndef HAVE_QIO
  char myname[] = "get_complex_source";
#endif
  

  /* See if source is already cached */
  if(source_type == COMPLEX_FIELD_STORE)
    return 1;

  /* Create in the c_src member */
  reset_qs(qs);
  alloc_cached_c_source(qs);
  
  if(source_type == COMPLEX_FIELD_FILE){
#ifdef HAVE_QIO
    if(qs->source_file_initialized)
      qs->infile = r_source_cmplx_scidac_open(qs->source_file);
    r_source_cmplx_scidac(qs->infile, qs->c_src, 
			  qs->x0, qs->y0, qs->z0, qs->t0);
#else
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
#endif
  }
  else if(source_type == COMPLEX_FIELD_FM_FILE)
    r_source_cmplx_fm_to_field(source_file, qs->c_src, 1, x0, y0, z0, t0);
  
  else if(source_type == CORNER_WALL)
    corner_wall(qs->c_src, t0);
  
  else if(source_type == EVEN_WALL)
    even_wall(qs->c_src, t0);
  
  else if(source_type == EVENANDODD_WALL)
    even_and_odd_wall(qs->c_src, t0);
  
  else if(source_type == EVENMINUSODD_WALL)
    even_minus_odd_wall(qs->c_src, t0);
  
  else if(source_type == GAUSSIAN) {
    gaussian_source(qs->c_src, r0, x0, y0, z0, t0);
  }      
  else if(source_type == POINT) {
    point_source(qs->c_src, x0, y0, z0, t0);
  }
  else if(source_type == WAVEFUNCTION_FILE){
    fnal_wavefunction(qs->c_src, x0, y0, z0, t0, a, source_file);
  }
  else {
    return 0;
  }
  
  /* Do momentum insertion */
  if(status == 0)insert_c_mom(qs->c_src, mom, x0, y0, z0, t0);
  
  /* Change source type to storage so we repeat the same field for
     other colors */
  
  qs->type = COMPLEX_FIELD_STORE;
  
  return 1;
}

int is_vector_source(int source_type){

  return
    source_type == RANDOM_COLOR_WALL ||
    source_type == VECTOR_FIELD_FILE ||
    source_type == VECTOR_FIELD_FM_FILE ||
    source_type == VECTOR_FIELD_STORE;
}

#ifdef HAVE_KS

static int get_vector_source(quark_source *qs){

  int source_type           = qs->type;
  int x0                    = qs->x0; 
  int y0                    = qs->y0; 
  int z0                    = qs->z0; 
  int t0                    = qs->t0;
  int *mom                  = qs->mom;
  char myname[] = "get_vector_source";

  /* See if source is already cached */
  if(source_type == VECTOR_FIELD_STORE)
    return 1;

  reset_qs(qs);
  alloc_cached_v_source(qs);
  
  /* Get a new color vector source on time slice t0 whenever the spin
     index is zero */
  
  if(source_type == RANDOM_COLOR_WALL)
    random_color_wall(qs->v_src, t0);

#ifdef HAVE_KS  
  else if(source_type == VECTOR_FIELD_FILE){
#ifdef HAVE_QIO
    if(r_source_vector(qs) != 0){
      node0_printf("%s: Failed to read source\n", myname);
      terminate(1);
    }
#else
    node0_printf("%s: QIO compilation required for this source\n",myname);
    terminate(1);
#endif
  }
  else if(source_type == VECTOR_FIELD_FM_FILE){
    /* Source file is in FNAL format. 
       We read the next record, store it in v_src and copy to src */
    
    if(qs->source_file_initialized == 0)
      r_source_open(qs);
    r_source_ks_fm(qs->kssf, qs->v_src, x0, y0, z0, t0);
  }
#endif
  else {
    return 0;
  }
  
  /* Insert momentum */
  insert_v_mom(qs->v_src, mom, x0, y0, z0, t0);

  return 1;
  
}

#endif

#ifdef HAVE_DIRAC

int is_dirac_source(int source_type){
  return
    source_type == DIRAC_FIELD_FILE ||
    source_type == DIRAC_FIELD_FM_FILE ||
    source_type == DIRAC_FIELD_STORE;

}

static int get_dirac_source(quark_source *qs, int spin, int color){

  int source_type           = qs->type;
  int x0                    = qs->x0; 
  int y0                    = qs->y0; 
  int z0                    = qs->z0; 
  int t0                    = qs->t0;
  char *source_file         = qs->source_file;
  int *mom                  = qs->mom;
  char myname[] = "get_dirac_source";
  
  /* See if source is already cached */
  if(source_type == DIRAC_FIELD_STORE)
    return 1;
  
  reset_qs(qs);
  alloc_cached_wv_source(qs);
  
  if(source_type == DIRAC_FIELD_FILE){
#ifdef HAVE_QIO
    if(r_source_dirac(qs) != 0){
      node0_printf("%s: Failed to read source\n", myname);
      terminate(1);
    }
#else
    node0_printf("%s: QIO compilation required for this source\n",myname);
    terminate(1);
#endif
  }
  else if(source_type == DIRAC_FIELD_FM_FILE){
    r_source_w_fm_to_field(source_file, qs->wv_src, spin, color, 
			   x0, y0, z0, t0);
  }
  else {
    return 0;
  }
  
  /* Do momentum insertion */
  insert_wv_mom(qs->wv_src, mom, x0, y0, z0, t0);

  return 1;
}  
#endif /* HAVE_DIRAC */

#ifdef HAVE_KS
/********************************************************************/
/* Construct the next basic color-vector field */
/********************************************************************/

/* Build a single source vector */

static int v_base_source(su3_vector *src, quark_source *qs)
{
  char myname[] = "v_base_source";
  
  /* Unpack structure */
  int source_type           = qs->type;
  int color                 = qs->color;
  int t0                    = qs->t0;
  int status = 0;

  /* zero src to be safe */
  clear_v_field(src);

  /* Sources from storage */

  if(source_type == COMPLEX_FIELD_STORE){
    if(qs->c_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    /* Load to the specified color and timeslice */
    if(src != NULL){
      insert_v_from_c(src, qs->c_src, color);
    }
  }

  else if(source_type == VECTOR_FIELD_STORE){
    if(qs->v_src == NULL){
      printf("%s: Can't copy from null field\n", myname);
      terminate(1);
    }
    copy_v_field(src, qs->v_src);
  }

  /* Sources built from a complex field */

  else if(is_complex_source(source_type)){
    get_complex_source(qs);
    insert_v_from_c(src, qs->c_src, color);
  }

  /* Sources built from a color vector field */

  else if(is_vector_source(source_type)){
      get_vector_source(qs);
    copy_v_field(src, qs->v_src);
  }

  else if(source_type == VECTOR_PROPAGATOR_FILE){
    /* In this case we aren't using the quark_source utilities to
       create the source.  It will be taken from the file we are going
       to read. This file is specified on the reload_serial_ksprop
       line.  The specified file must have source records interleaved
       with solution records. */
  }

  /* Dirac vector sources */
  /* In this case we extract one color vector for each spin.
     Thus one Dirac field source yields four color sources.
     This routine must be called for four spins for each source */

#ifdef HAVE_DIRAC
  else if(is_dirac_source(source_type)){
    int spin_snk = qs->spin_snk;

    /* Fetch a new source only when spin_snk = 0 */
    if(spin_snk == 0){
      int spin_src, color_src;
      spin_src = convert_ksource_to_spin(qs->ksource);
      color_src = convert_ksource_to_color(qs->ksource);
      get_dirac_source(qs, spin_src, color_src);
    }

    extract_v_from_wv(src, qs->wv_src, spin_snk);
  }
#endif

  else {
    node0_printf("%s: Unrecognized source type %d for a color-vector field\n",
		 myname, source_type);
    terminate(1);
  }

  /* Apply subset mask */

  subset_mask_v(src, qs->subset, t0);

  return status;

} /* v_base_source */

#endif /* HAVE_KS */

/********************************************************************/
/* Construct the next basic Dirac source field                            */
/********************************************************************/

#ifdef HAVE_DIRAC

static int wv_base_source(wilson_vector *src, quark_source *qs)
{
  char myname[] = "wv_base_source";
  
  int status = 0;

  /* Unpack structure */
  int source_type             = qs->type;
  int ksource                 = qs->ksource;
  
  int t0                      = qs->t0;

  /* Get the current "color" and spin  from ksource */
  int color = convert_ksource_to_color(ksource);
  int spin  = convert_ksource_to_spin(ksource);
  
  /* zero src to be safe */
  clear_wv_field(src);
  
  /* Complex field sources */

  if(is_complex_source(source_type)){

    get_complex_source(qs);
    insert_wv_from_c(src, qs->c_src, spin, color);
  }
  
#ifdef HAVE_KS
  /* Color vector sources */
  else if(is_vector_source(source_type)){
    /* Get a new color vector source on time slice t0 whenever the spin
       index is zero */
    if(spin == 0)
      get_vector_source(qs);
    insert_wv_from_v(src, qs->v_src, spin);
  }
#endif

  /* Dirac vector sources */

  else if(is_dirac_source(source_type)){

    get_dirac_source(qs, spin, color);
    copy_wv_field(src, qs->wv_src);
  }

  else if(source_type == DIRAC_PROPAGATOR_FILE){
    /* In this case we aren't using the quark_source utilities to
       create the source.  It will be taken from the file
       when we read the propagator itself. */
  }

  else {
    node0_printf("%s: Unrecognized source type %d\n", myname, source_type);
    terminate(1);
  }

  /* Rescale the source */
  if(qs->scale_fact != 1.0){
    int i; site *s;
    FORALLSITES(i,s){
      if(t0 == ALL_T_SLICES || s->t == t0)
	scalar_mult_wvec(src+i, qs->scale_fact, src+i);
    }
  }

  /* Apply subset mask */

  subset_mask_wv(src, qs->subset, t0);

  return status;

} /* wv_base_source */

#endif

#ifdef HAVE_KS

/********************************************************************/
/* Generate the next color source field in the list */
/********************************************************************/

int v_source_field(su3_vector *src, quark_source *qs)
{
  char myname[] = "v_source_field";
  int status = 0;
  quark_source_sink_op* op;

  /* Increment the source counter, depending on source type */

#ifdef HAVE_DIRAC

  if(is_dirac_source(qs->type)){
    if(qs->spin_snk == 0){
      qs->ksource++;
      if(qs->ksource >= qs->nsource){
	node0_printf("%s: Exceeded source count %d\n", myname, qs->nsource);
	status++;
      }
    }
  }
  else{
    qs->color++;
    if(qs->color >= qs->ncolor){
      node0_printf("%s: Exceeded source color %d\n", myname, qs->ncolor);
      status++;
    }
  }

#else

  qs->color++;
  if(qs->color >= qs->ncolor){
    node0_printf("%s: Exceeded source color %d\n", myname, qs->ncolor);
    status++;
  }

#endif

#ifdef HAVE_KS

  /* Get the next base source */

  status += v_base_source(src, qs);

  /* Apply the chain of operators to the base source */

  op = qs->op;
  while(op != NULL){
    v_field_op(src, op, qs->subset, qs->t0);
    op = op->op;
  }

#endif

  return status;
}

/* Deprecated entry point */

int v_source_site(field_offset src, quark_source *qs)
{
  int i;
  site *s;
  su3_vector *t_src;
  int status = 0;

#define PAD 0
  t_src  = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  
  if(t_src == NULL){
    printf("v_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  status = v_source_field(t_src, qs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(su3_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);

  return status;

} /* v_source_site */

#endif /* HAVE_KS */

/********************************************************************/
/* Generate the next source field in the list */
/********************************************************************/

#ifdef HAVE_DIRAC

int wv_source_field(wilson_vector *src, quark_source *qs)
{
  char myname[] = "wv_source_field";
  int status = 0;
  quark_source_sink_op* op;

  /* Increment the source counter */

  qs->ksource++;

  if(qs->ksource >= qs->nsource){
    node0_printf("%s: Exceeded source number %d\n", myname, qs->nsource);
    status++;
  }

  qs->color = convert_ksource_to_color(qs->ksource);

  /* Get the next base source */

  status += wv_base_source(src, qs);

  /* Apply the chain of operators to the base source */

  op = qs->op;
  while(op != NULL){
    wv_field_op(src, op, qs->subset, qs->t0);
    op = op->op;
  }

  return status;
}

/* Deprecated entry point */

int wv_source_site(field_offset src, quark_source *qs)
{
  int i;
  site *s;
  wilson_vector *t_src;
  int status = 0;

#define PAD 0
  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  
  if(t_src == NULL){
    printf("wv_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  status = wv_source_field(t_src, qs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);

  return status;

} /* wv_source_site */

#endif /* HAVE_DIRAC */

static int ask_quark_source( FILE *fp, int prompt, int *source_type, 
			     char *descrp )
{
  char *savebuf;
  char myname[] = "ask_quark_source";

  if (prompt==1){
    printf("enter ");
    printf("'complex_field', ");
    printf("'complex_field_fm', ");
    printf("'corner_wall', ");
    printf("\n     ");
    printf("'even_wall', ");
    printf("'evenandodd_wall', ");
    printf("'evenminusodd_wall', ");
    printf("\n     ");
    printf("'gaussian', ");
    printf("'point', ");
    printf("'wavefunction', ");
    printf("\n     ");
    printf("'random_color_wall', ");
    printf("'vector_field', ");
    printf("'vector_field_fm' ");
    printf("\n     ");
    printf("'vector_propagator_file', ");
    printf("'dirac_field', ");
    printf("'dirac_field_fm', ");
    printf("'dirac_propagator_file' ");
    printf("\nfor source type\n");
  }

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  /* Complex field sources */
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
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("wavefunction",savebuf) == 0 ){
    *source_type = WAVEFUNCTION_FILE;
    strcpy(descrp,"wavefunction");
  }
  /* Vector field sources */
  else if(strcmp("random_color_wall",savebuf) == 0 ){
    *source_type = RANDOM_COLOR_WALL;
    strcpy(descrp,"random_color_wall");
  }
  else if(strcmp("vector_field",savebuf) == 0 ) {
    *source_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
  }
  else if(strcmp("vector_field_fm",savebuf) == 0 ) {
    *source_type = VECTOR_FIELD_FM_FILE;
    strcpy(descrp,"vector_field_fm");
  }
  else if(strcmp("vector_propagator_file",savebuf) == 0 ) {
    *source_type = VECTOR_PROPAGATOR_FILE;
    strcpy(descrp,"vector_propagator_file");
  }

  /* Dirac field sources */
  else if(strcmp("dirac_field",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
  }
  else if(strcmp("dirac_field_fm",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FM_FILE;
    strcpy(descrp,"dirac_field_fm");
  }
  else if(strcmp("dirac_propagator_file",savebuf) == 0 ) {
    *source_type = DIRAC_PROPAGATOR_FILE;
    strcpy(descrp,"dirac_propagator_file");
  }
  else{
    printf("%s: ERROR IN INPUT: source command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_source */

/* For parsing the subset mask */
static int encode_mask(int *mask, char c_mask[]){
  int status = 0;
  if(strcmp(c_mask,"full")==0)
    *mask =  FULL;
  else if(strcmp(c_mask,"corner")==0)
    *mask = HYPERCUBE;
  else{
    node0_printf("Expected 'full' or 'corner' but got %s\n",c_mask);
    status++;
  }
  return status;
}

static char *decode_mask(int mask){
  if(mask == FULL)
    return "full";
  else if (mask == HYPERCUBE)
    return "corner";
  else
    return "??";
}

int ask_starting_source( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_starting_source";

  if (prompt==1) printf("enter 'load_source_serial', or 'load_source_parallel'\n");
  savebuf = get_next_tag(fp, "load source command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("load_source_serial",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("load_source_parallel",savebuf) == 0 ) {
    *flag = RELOAD_PARALLEL;
  }
  /* Backward compatibility */
  else if(strcmp("load_source",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }

  if(prompt==1)printf("enter name of file containing the source\n");
  status=fscanf(fp," %s",filename);
  if(status !=1) {
    printf("\n%s(%d): ERROR IN INPUT: error reading file name\n",
	   myname, this_node); 
    return 1;
  }
  printf("%s\n",filename);
  return 0;
}

/* Get the additional input parameters needed to specify the source */
#define IF_OK if(status==0)
static int get_quark_source(int *status_p, FILE *fp, int prompt, 
			    quark_source *qs){
  
  int  source_loc[4] = { 0,0,0,0 };
  char source_file[MAXFILENAME] = "";
  Real source_r0 = 0.;
  int  source_type = qs->type;
  int  status = *status_p;
  
  /* Complex field sources */
  if ( source_type == POINT ){
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
  }
  else if ( source_type == CORNER_WALL ||
	    source_type == EVEN_WALL ||
	    source_type == EVENANDODD_WALL ||
	    source_type == EVENMINUSODD_WALL ){
    IF_OK status += get_i(fp, prompt, "t0", &source_loc[3]);
  }
  else if ( source_type == COMPLEX_FIELD_FILE ||
	    source_type == COMPLEX_FIELD_FM_FILE){
    //    IF_OK status += get_i(fp, prompt, "t0", &source_loc[3]);
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
    IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
  }
  else if ( source_type == GAUSSIAN ){
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
    /* width: psi=exp(-(r/r0)^2) */
    IF_OK status += get_f(fp, prompt,"r0", &source_r0 );
  }
  else if ( source_type == WAVEFUNCTION_FILE ){
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
    IF_OK status += get_f(fp, prompt, "a", &qs->a);
    IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
  }
  
  /* Vector field sources */
  else if ( source_type == RANDOM_COLOR_WALL ){
    IF_OK status += get_i(fp, prompt, "t0", &source_loc[3]);
    IF_OK status += get_i(fp, prompt, "ncolor", &(qs->ncolor));
    IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
  }
  else if ( source_type == VECTOR_FIELD_FILE ){
    //    IF_OK status += get_i(fp, prompt, "t0", &source_loc[3]);
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
    IF_OK status += ask_starting_source(fp, prompt, &qs->sourceflag, source_file);
    IF_OK status += get_i(fp, prompt, "ncolor", &(qs->ncolor));
    IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
  }
  else if ( source_type == VECTOR_FIELD_FM_FILE ){
    //    IF_OK status += get_i(fp, prompt, "t0", &source_loc[3]);
    IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
    IF_OK status += ask_starting_source(fp, prompt, &qs->sourceflag, source_file);
    IF_OK status += get_i(fp, prompt, "ncolor", &(qs->ncolor));
    IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
  }
  else{
    return 0;
  }
  
  qs->x0    = source_loc[0];
  qs->y0    = source_loc[1];
  qs->z0    = source_loc[2];
  qs->t0    = source_loc[3];
  qs->r0    = source_r0;
  strcpy(qs->source_file,source_file);
  
  *status_p = status;
  return 1;

} /* get_quark_source */

/* Get the additional input parameters needed to specify the source */
int get_v_quark_source(FILE *fp, int prompt, quark_source *qs){
  
  char c_mask[16];
  char source_label[MAXSRCLABEL];
  int  source_type;
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_quark_source(fp, prompt, &source_type,
				     qs->descrp);
  IF_OK qs->type  = source_type;
  IF_OK qs->orig_type  = source_type;  /* In case we change the type */
  
  /* Specify the subset */
  IF_OK status += get_s(fp, prompt, "subset", c_mask);
  IF_OK status += encode_mask(&qs->subset, c_mask);
  
  /* Complex field sources */
  IF_OK {

    if ( get_quark_source(&status, fp, prompt, qs) );
    else if ( source_type == VECTOR_PROPAGATOR_FILE ){
      IF_OK status += get_i(fp, prompt, "ncolor", &(qs->ncolor));
    }
    else {
      printf("Source type not supported in this application\n");
      status++;
    }
  }
  
  IF_OK status += get_s(stdin, prompt, "source_label", source_label);
  /* Source label '(null)' suppresses the label */
  if(strcmp(source_label,"(null)")==0)source_label[0]='\0';
  strncpy(qs->label,source_label,MAXSRCLABEL);

  return status;

} /* get_v_quark_source */

/* Get the input parameters needed to specify the base source */
int get_wv_quark_source(FILE *fp, int prompt, quark_source *qs){
  
  char myname[] = "get_wv_quark_source";
  char c_mask[16];
  int  source_type;
  char source_label[MAXSRCLABEL];
  int  source_loc[4] = { 0,0,0,0 };
  char source_file[MAXFILENAME] = "";
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_quark_source(fp, prompt,&source_type,
				     qs->descrp);
  IF_OK qs->type  = source_type;
  IF_OK qs->orig_type  = source_type;  /* In case we change the type */
  
  IF_OK {

    /* Specify the subset */
    IF_OK status += get_s(fp, prompt, "subset", c_mask);
    IF_OK status += encode_mask(&qs->subset, c_mask);

    /* Complex field sources */
    if ( get_quark_source( &status, fp, prompt, qs) );

    /* Dirac field sources */
    else if ( source_type == DIRAC_FIELD_FILE ){
      //      IF_OK status += get_i(fp, prompt, "t0", &qs->t0);
      IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
      IF_OK status += ask_starting_source(fp, prompt, &qs->sourceflag, source_file);
      IF_OK status += get_i(fp, prompt, "nsource", &(qs->nsource));
      IF_OK {
	int ncolor = convert_ksource_to_color(qs->nsource);
	if(qs->nsource != convert_ncolor_to_nsource(ncolor)){
	  printf("%s: ERROR nsource must be a multiple of the number of spins\n",
		 myname);
	  status++;
	}
	qs->ncolor = ncolor;
      }
      IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
      strcpy(qs->source_file,source_file);
      qs->x0    = source_loc[0];
      qs->y0    = source_loc[1];
      qs->z0    = source_loc[2];
      qs->t0    = source_loc[3];
    }
    else if ( source_type == DIRAC_FIELD_FM_FILE ){
      //      IF_OK status += get_i(fp, prompt, "t0", &qs->t0);
      IF_OK status += get_vi(fp, prompt, "origin", source_loc, 4);
      IF_OK status += ask_starting_source(fp, prompt, &qs->sourceflag, source_file);
      IF_OK status += get_vi(fp, prompt, "momentum", qs->mom, 3);
      strcpy(qs->source_file,source_file);
      qs->nsource = 12;  /* Required */
      qs->ncolor  = 3;
      qs->x0    = source_loc[0];
      qs->y0    = source_loc[1];
      qs->z0    = source_loc[2];
      qs->t0    = source_loc[3];
    }
    else if ( source_type == VECTOR_PROPAGATOR_FILE ){
      IF_OK status += get_i(fp, prompt, "ncolor", &(qs->ncolor));
    }
    else if ( source_type == DIRAC_PROPAGATOR_FILE ){
      IF_OK status += get_i(fp, prompt, "nsource", &(qs->nsource));
      IF_OK {
	int ncolor = convert_ksource_to_color(qs->nsource);
	if(qs->nsource != convert_ncolor_to_nsource(ncolor)){
	  printf("%s: ERROR nsource must be a multiple of the number of spins\n",
		 myname);
	  status++;
	}
	qs->ncolor = ncolor;
      }
    }
    else {
      printf("%s: Source type not supported in this application\n",myname);
      status++;
    }
  }

  if ( qs->type == RANDOM_COLOR_WALL ||
       qs->type == VECTOR_FIELD_FILE ||
       qs->type == VECTOR_FIELD_FM_FILE )
    qs->nsource = convert_ncolor_to_nsource(qs->ncolor);

  /* Provision for an arbitrary scale factor.  The scale factor is
     applied to the source and then removed in the hadron correlator. */

  /* When SCALE_PROP is defined, all Dirac sources must specify a scale_factor */
#ifdef SCALE_PROP
  IF_OK status += get_f(fp, prompt, "scale_factor", &qs->scale_fact);
#endif

  IF_OK status += get_s(stdin, prompt, "source_label", source_label);

  /* Source label '(null)' suppresses the label */
  if(strcmp(source_label,"(null)")==0)source_label[0]='\0';
  strncpy(qs->label,source_label,MAXSRCLABEL);
  
  return status;

} /* get_wv_quark_source */

/*--------------------------------------------------------------------*/
/* Print the parameters of the source used.
   This is intended to be part of the metadata
   in the FNAL-style correlator file*/
/*--------------------------------------------------------------------*/

#define NTAG 30
/* Create a fixed-width tag for tidy output */
static char *make_tag(char prefix[], char tag[]){
  static char full_tag[NTAG];
  full_tag[0] = '\0';
  strncat(full_tag, prefix, NTAG-1);
  strncat(full_tag, "_", NTAG-1-strlen(full_tag));
  strncat(full_tag, tag, NTAG-1-strlen(full_tag));
  strncat(full_tag, ":", NTAG-1-strlen(full_tag));
  memset(full_tag+strlen(full_tag), ' ', NTAG-1-strlen(full_tag));
  full_tag[NTAG-1] = '\0';
  return full_tag;
}

/*--------------------------------------------------------------------*/
void print_source_info(FILE *fp, char prefix[], quark_source *qs){
  int source_type = qs->orig_type;

  fprintf(fp,"%s %s\n", make_tag(prefix, "type"), qs->descrp);
  fprintf(fp,"%s %s\n", make_tag(prefix, "subset"), decode_mask(qs->subset));

  if ( source_type == POINT ){
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
  }
  else if ( source_type == CORNER_WALL ||
	    source_type == EVEN_WALL ||
	    source_type == EVENANDODD_WALL ||
	    source_type == EVENMINUSODD_WALL ){
    fprintf(fp,"%s %d\n", make_tag(prefix, "t0"), qs->t0);
  }
  else if ( source_type == COMPLEX_FIELD_FILE ||
	    source_type == COMPLEX_FIELD_FM_FILE){
    //    fprintf(fp,"%s %d\n", make_tag(prefix, "t0"), qs->t0);
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
    fprintf(fp,"%s %s\n", make_tag(prefix, "file"), qs->source_file);
    fprintf(fp,"%s %d %d %d\n", make_tag(prefix, "mom"), qs->mom[0],
	    qs->mom[1], qs->mom[2]);
  }
  else if ( source_type == GAUSSIAN ){
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
    fprintf(fp,"%s %g\n", make_tag(prefix, "r0"), qs->r0);
  }
  else if ( source_type == WAVEFUNCTION_FILE ){
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
    fprintf(fp,"%s %s\n", make_tag(prefix, "file"), qs->source_file);
    fprintf(fp,"%s %g\n", make_tag(prefix, "a"), qs->a);
    fprintf(fp,"%s %d %d %d\n", make_tag(prefix, "mom"), qs->mom[0],
	    qs->mom[1], qs->mom[2]);
  }
  else if ( source_type == RANDOM_COLOR_WALL ){
    fprintf(fp,"%s %d\n", make_tag(prefix, "t0"), qs->t0);
    fprintf(fp,"%s %d\n", make_tag(prefix, "ncolor"), qs->ncolor);
    fprintf(fp,"%s %d %d %d\n", make_tag(prefix, "mom"), qs->mom[0],
	    qs->mom[1], qs->mom[2]);
  }
  else if ( source_type == VECTOR_FIELD_FILE ||
	    source_type == VECTOR_FIELD_FM_FILE ){
    //    fprintf(fp,"%s %d\n", make_tag(prefix, "t0"), qs->t0);
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
    fprintf(fp,"%s %s\n", make_tag(prefix, "file"), qs->source_file);
    fprintf(fp,"%s %d\n", make_tag(prefix, "ncolor"), qs->ncolor);
    fprintf(fp,"%s %d %d %d\n", make_tag(prefix, "mom"), qs->mom[0],
	    qs->mom[1], qs->mom[2]);
  }
  else if ( source_type == VECTOR_PROPAGATOR_FILE ){
    fprintf(fp,"%s %d\n", make_tag(prefix, "ncolor"), qs->ncolor);
  }
  else if ( source_type == DIRAC_FIELD_FILE ||
	    source_type == DIRAC_FIELD_FM_FILE ){
    //    fprintf(fp,"%s %d\n", make_tag(prefix, "t0"), qs->t0);
    fprintf(fp,"%s %d %d %d %d\n", make_tag(prefix, "origin"), 
	    qs->x0, qs->y0, qs->z0, qs->t0);
    fprintf(fp,"%s %s\n", make_tag(prefix, "file"), qs->source_file);
    fprintf(fp,"%s %d\n", make_tag(prefix, "nsource"), qs->nsource);
    fprintf(fp,"%s %d %d %d\n", make_tag(prefix, "mom"), qs->mom[0],
	    qs->mom[1], qs->mom[2]);
  }
  else if ( source_type == DIRAC_PROPAGATOR_FILE ){
    fprintf(fp,"%s %d\n", make_tag(prefix, "nsource"), qs->nsource);
  }
} /* print_source_info */
