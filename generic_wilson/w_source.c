/************************** w_source.c *****************************/
/* MIMD version 7 */

/*  2/15/98 to fail on illegal source type CD */

/* Initialize a source for the inverter */

#include "generic_wilson_includes.h"
#include "../include/io_wprop.h"
#include <string.h>

/* Utilities */

/* Initialize by calling this function before using wqs */
void init_wqs(wilson_quark_source *wqs){
  wqs->type = 0;
  wqs->c_src = NULL;
  wqs->wv_src = NULL;
}

void alloc_wqs_wv_src(wilson_quark_source *wqs)
{
  if(wqs->wv_src != NULL){
    printf("alloc_wqs_wv_src(%d): source is already allocated? \n",this_node);
    terminate(1);
  }
  wqs->wv_src = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(wqs->wv_src == NULL){
    printf("alloc_wqs_wv_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->wv_src, 0, sizeof(wilson_vector)*sites_on_node);
}

void alloc_wqs_c_src(wilson_quark_source *wqs)
{
  if(wqs->c_src != NULL){
    printf("alloc_wqs_c_src(%d): source is already allocated? \n",this_node);
    terminate(1);
  }
  wqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(wqs->c_src == NULL){
    printf("alloc_wqs_c_src(%d): No room for source field\n",this_node);
    terminate(1);
  }
  memset(wqs->c_src, 0, sizeof(complex)*sites_on_node);
}

/* Must be called to free */
void reset_wqs(wilson_quark_source *wqs){
  if(wqs->c_src != NULL){free(wqs->c_src); wqs->c_src = NULL;}
  if(wqs->wv_src != NULL){free(wqs->wv_src); wqs->wv_src = NULL;}
}

/* Copy c-source */
void copy_c_src(complex *dest, wilson_vector *src, int spin, int color){
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i].d[spin].c[color];
  }
}

/* Copy wv-source */
void copy_wv_src(wilson_vector *dest, wilson_vector *src)
{
  int i; site *s;
  FORALLSITES(i,s){
    dest[i] = src[i];
  }
}

/* Choose the specific USQCD format appropriate to the source type */
int choose_usqcd_file_type(int source_type){
  int file_type;

  switch(source_type){
  case POINT:
  case GAUSSIAN:
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_STORE:
    file_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case COVARIANT_GAUSSIAN:
  case DIRAC_FIELD_FILE:
  case DIRAC_FIELD_STORE:
    file_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

/********************************************************************/
/* Construct the source field */
/********************************************************************/

void w_source_field(wilson_vector *src, wilson_quark_source *wqs)
{
  /* src has size wilson_vector */
  register int i;
  register site *s; 
  
  short my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  
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
  complex *c_src            = wqs->c_src;
  wilson_vector *wv_src     = wqs->wv_src;
  
  
  /* zero src to be safe */
  FORALLSITES(i,s) {
    clear_wvec( src+i ); 
  }
  
  /* If we are taking the source from a propagator file, set the
     appropriate action */
  if(source_type == PROP_FILE){
    if(wqs->c_src != NULL && wqs->wv_src == NULL)
      source_type = COMPLEX_FIELD_STORE;
    if(wqs->wv_src != NULL && wqs->c_src == NULL)
      source_type = DIRAC_FIELD_STORE;
    else{
      node0_printf("w_source_field: propagator file source ambiguous or missing.\n");
      terminate(1);
    }
  }

  /* Unless we are taking the source from storage, clear it before
     reconstructing it */
  if(source_type != COMPLEX_FIELD_STORE &&
     source_type != DIRAC_FIELD_STORE)
    reset_wqs(wqs);
  
  if(source_type == POINT) {
    /* load 1.0 into source at cooordinates given by source_coord */
    /* initialize src to be a delta function at point x0,y0,z0,t0 */
    /* Save a copy of the source in wqs->c_src */

    alloc_wqs_c_src(wqs);

    if(node_number(x0,y0,z0,t0) == mynode()){
      i = node_index(x0,y0,z0,t0);
      wqs->c_src[i].real = 1.0;
      src[i].d[spin].c[color].real = 1.0;
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

      wqs->c_src[i].real = 
	src[i].d[spin].c[color].real = (Real)exp((double)(- radius2));
    }
  }
  else if(source_type == COVARIANT_GAUSSIAN){
    /* Set delta function source */
    if(node_number(x0,y0,z0,t0) == this_node){
      i = node_index(x0,y0,z0,t0);
      src[i].d[spin].c[color].real = 1.;
    }
    /* Then smear */
    gauss_smear_field(src, r0, iters, t0);
    /* Save a copy */
    alloc_wqs_wv_src(wqs);
    copy_wv_src(wqs->wv_src, src);
  }
  else if(source_type == COMPLEX_FIELD_FILE){
    /* Load to the specified spin, color and timeslice */
    r_source_w_fm_to_field(source_file, src, spin, color, t0, source_type);
    /* Save a copy */
    alloc_wqs_c_src(wqs);
    copy_c_src(wqs->c_src, src, spin, color);
  }
  else if(source_type == DIRAC_FIELD_FILE){
    /* Ignore spin, color and load the whole wilson vector from the file */
    r_source_w_fm_to_field(source_file, src, spin, color, t0, source_type);
    /* Save a copy */
    alloc_wqs_wv_src(wqs);
    copy_wv_src(wqs->wv_src, src);
  }
  else if(source_type == COMPLEX_FIELD_STORE){
    if(c_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Load to the specified spin, color and timeslice */
    FORALLSITES(i,s){
      if(t0 == ALL_T_SLICES || s->t == t0)
	src[i].d[spin].c[color] = c_src[i];
    }
  }
  else if(source_type == DIRAC_FIELD_STORE){
    if(wv_src == NULL){
      printf("w_source: Can't copy from null field\n");
      terminate(1);
    }
    /* Ignore spin, color and load the whole wilson vector from storage */
    FORALLSITES(i,s){
      if(t0 == ALL_T_SLICES || s->t == t0)
	src[i] = wv_src[i];
    }
  }
  else {
    node0_printf("w_source: Unrecognized source type %d\n",source_type);
    terminate(1);
  }
} /* w_source_field */

void w_source_site(field_offset src, wilson_quark_source *wqs)
{
  int i;
  site *s;
  wilson_vector *t_src;

#define PAD 0
  t_src  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  
  if(t_src == NULL){
    printf("w_source_site(%d): Can't allocate t_src\n",this_node);
    terminate(1);
  }
  
  w_source_field(t_src, wqs);

  /* copy src back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,src) = t_src[i];
  }

  free(t_src);

} /* w_source_site */


/* 
   Initialize a sink for the inverter, using a wilson vector as the
   storage container .
 */

void w_sink_field(wilson_vector *snk,wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  
  int my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  
  /* Unpack structure.  We don't use member t0 here. */
  int color         = wqs->color;
  int spin          = wqs->spin;
  int sink_type     = wqs->type;
  int x0            = wqs->x0; 
  int y0            = wqs->y0; 
  int z0            = wqs->z0; 
  int t0            = wqs->t0;
  Real r0           = wqs->r0;
  char *source_file = wqs->source_file;

    if(sink_type == POINT) {
	/* load 1.0 into sink at cooordinates given by sink_coord */
	/* initialize snk to be a delta function at point x0,y0,z0 */

	for(t0=0;t0<nt;t0++){
	    if(node_number(x0,y0,z0,t0) == mynode()){
		i = node_index(x0,y0,z0,t0);
		snk[i].d[spin].c[color].real = 1.0;
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

	    snk[i].d[spin].c[color].real = 
	      (Real)exp((double)(- radius2));
	}
    }
    else if(sink_type == COMPLEX_FIELD_FILE){
      r_source_w_fm_to_field(source_file, snk, spin, color, t0, sink_type);
    }
    else if(sink_type == DIRAC_FIELD_FILE){
      r_source_w_fm_to_field(source_file, snk, spin, color, t0, sink_type);
    }
    else {
      node0_printf("w_sink_field: Unrecognized sink type %d\n",sink_type);
      terminate(1);
    }
    
} /* w_sink_field */


/* snk must be type wilson_vector */
void w_sink_site(field_offset snk, wilson_quark_source *wqs)
{
  int i;
  site *s;
  wilson_vector *t_snk;

#define PAD 0
  t_snk  = (wilson_vector *) malloc((sites_on_node+PAD)*sizeof(wilson_vector));
  
  if(t_snk == NULL){
    printf("w_sink_site(%d): Can't allocate t_snk\n",this_node);
    terminate(1);
  }
  
  w_sink_field(t_snk, wqs);

  /* copy snk back */
  FORALLSITES(i,s) {
    *(wilson_vector *)F_PT(s,snk) = t_snk[i];
  }

  free(t_snk);

} /* w_sink_site */

/* 
   Initialize a sink for the inverter, using a wilson number as the
   storage container.
 */

void w_sink_scalar(field_offset snk,wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  
  int my_x,my_y,my_z,t0;
  Real rx,ry,rz,radius2;
  int sink_type;
  int x0, y0, z0;
  Real r0;
  complex zero ;
  /****************************************/

  zero.real = zero.imag = 0.0 ; 

  /** zero the smearing function first ****/
  FORALLSITES(i,s)
  {
    *((complex *)F_PT(&(lattice[i]),snk)) = zero ;
  }


  /* Unpack structure.  We don't use member t0 here. */
  x0 = wqs->x0; y0 = wqs->y0; z0 = wqs->z0;
  sink_type = wqs->type;
  r0 = wqs->r0;

    if(sink_type == POINT) {
	/* load 1.0 into sink at cooordinates given by sink_coord */
	/* initialize snk to be a delta function at point x0,y0,z0 */

	for(t0=0;t0<nt;t0++){
	    if(node_number(x0,y0,z0,t0) == mynode()){
		i = node_index(x0,y0,z0,t0);
		((complex *)F_PT(&(lattice[i]),snk))->real = 1.0;
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

	    ((complex *)F_PT(s,snk))->real = (Real)exp((double)(- radius2));
	}
    }



} /* w_sink_scalar */


int ask_quark_source( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_quark_source";

  if (prompt!=0)
    printf("enter 'point', 'gaussian', 'covariant_gaussian', 'complex_field', 'dirac_field', 'prop_file' for source type\n");

  savebuf = get_next_tag(fp, "quark source command", myname);
  if (savebuf == NULL)return 1;

  if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
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
  else if(strcmp("dirac_field",savebuf) == 0 ) {
    *source_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
  }
  else if(strcmp("prop_file",savebuf) == 0 ) {
    *source_type = PROP_FILE;
    strcpy(descrp,"propagator_file");
  }
  else{
    printf("ask_source: ERROR IN INPUT: source command %s is invalid\n",
	   savebuf); 
    return 1;
  }
  
  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_source */

#define IF_OK if(status==0)

/* Get the additional input parameters needed to specify the source */
int get_quark_source(FILE *fp, int prompt, wilson_quark_source *wqs){
  
  Real source_r0 = 0;
  int  source_type;
  int  source_loc[4] = { 0,0,0,0 };
  int  source_iters = 0;
  char source_file[MAXFILENAME] = "";
  int  status = 0;
  
  /* Get antiquark source type */
  IF_OK status += ask_quark_source(fp, prompt,&source_type,
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
    else if ( source_type == COVARIANT_GAUSSIAN ){
      IF_OK status += get_vi(stdin, prompt, "origin", source_loc, 4);
      IF_OK status += get_f(stdin, prompt, "r0", &source_r0);
      IF_OK status += get_i(stdin, prompt, "source_iters", &source_iters);
    }
    else if ( source_type == COMPLEX_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else if ( source_type == DIRAC_FIELD_FILE ){
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
      IF_OK status += get_s(stdin, prompt, "load_source", source_file);
    }
    else if (source_type == PROP_FILE ){
      /* This information is not taken from the file */
      IF_OK status += get_i(stdin, prompt, "t0", &source_loc[3]);
    }
    else {
      printf("Source type not supported in this application\n");
      status++;
    }
  }
  
  wqs->r0    = source_r0;
  wqs->x0    = source_loc[0];
  wqs->y0    = source_loc[1];
  wqs->z0    = source_loc[2];
  wqs->t0    = source_loc[3];
  wqs->iters = source_iters;
  strcpy(wqs->source_file,source_file);
  
  return status;
} /* get_quark_source */

int get_quark_sink(FILE *fp, int prompt, wilson_quark_source *wqs){

  Real sink_r0 = 0;
  int  sink_type;
  int  sink_loc[3] = { 0,0,0 };
  char sink_file[MAXFILENAME] = "";
  int  status = 0;

  /* Get antiquark source type */
  IF_OK status += ask_quark_source(fp, prompt,&sink_type,
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
      else if ( sink_type == COMPLEX_FIELD_FILE ){
	IF_OK status += get_s(stdin, prompt, "load_source", sink_file);
      }
      else {
	printf("Sink type not supported in this application\n");
	status++;
      }
    }
    
    wqs->r0    = sink_r0;
    wqs->x0    = sink_loc[0];
    wqs->y0    = sink_loc[1];
    wqs->z0    = sink_loc[2];
    wqs->t0    = ALL_T_SLICES;
    strcpy(wqs->source_file,sink_file);

    return status;
} /* get_quark_sink */


