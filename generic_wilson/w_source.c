/************************** w_source.c *****************************/
/* MIMD version 7 */

/*  2/15/98 to fail on illegal source type CD */

/* Initialize a source for the inverter */

#include "generic_wilson_includes.h"
#include <string.h>

void w_source(field_offset src,wilson_quark_source *wqs)
{
  /* src has size wilson_vector */
  register int i;
  register site *s; 
  
  short my_x,my_y,my_z;
  Real rx,ry,rz,radius2;
  int color, spin, source_type;
  int x0, y0, z0, t0;
  Real r0;
  
  /* Unpack structure */
  x0 = wqs->x0; y0 = wqs->y0; z0 = wqs->z0; t0 = wqs->t0;
  color = wqs->color; spin = wqs->spin;
  source_type = wqs->type;
  r0 = wqs->r0;
  
/*printf("WSOURCE: source = %d\n",source_type); */
	
    /* zero src to be safe */
    FORALLSITES(i,s) {
	clear_wvec((wilson_vector *)F_PT(s,src)); 
    }

    if(source_type == POINT) {
	/* load 1.0 into source at cooordinates given by source_coord */
	/* initialize src to be a delta function at point x0,y0,z0,t0 */

	if(node_number(x0,y0,z0,t0) == mynode()){
	    i = node_index(x0,y0,z0,t0);
	    ((wilson_vector *)F_PT(&(lattice[i]),src))->
		d[spin].c[color].real = 1.0;
	}
    }
    else if(source_type == GAUSSIAN) {
	/* Gaussian trial source centered on  x0,y0,z0,t0 */

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

	    ((wilson_vector *)F_PT(s,src))->d[spin].c[color].real = 
		(Real)exp((double)(- radius2));
	}
      }
  else {
    node0_printf("w_source: Unrecognized source type %d\n",source_type);
    terminate(1);
  }
    

} /* w_source */

/* 
   Initialize a sink for the inverter, using a wilson vector as the
   storage container .
 */

void w_sink(field_offset snk,wilson_quark_source *wqs)
{
  register int i;
  register site *s; 
  
  int my_x,my_y,my_z,t0;
  Real rx,ry,rz,radius2;
  int color, spin, sink_type;
  int x0, y0, z0;
  Real r0;
  
/*printf("WSINK: sink = %d\n",sink); */
	

  /* Unpack structure.  We don't use member t0 here. */
  x0 = wqs->x0; y0 = wqs->y0; z0 = wqs->z0;
  color = wqs->color; spin = wqs->spin;
  sink_type = wqs->type;
  r0 = wqs->r0;

    if(sink_type == POINT) {
	/* load 1.0 into sink at cooordinates given by sink_coord */
	/* initialize snk to be a delta function at point x0,y0,z0 */

	for(t0=0;t0<nt;t0++){
	    if(node_number(x0,y0,z0,t0) == mynode()){
		i = node_index(x0,y0,z0,t0);
		((wilson_vector *)F_PT(&(lattice[i]),snk))->
		    d[spin].c[color].real = 1.0;
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

	    ((wilson_vector *)F_PT(s,snk))->d[spin].c[color].real = 
		(Real)exp((double)(- radius2));
	}
    }

} /* w_sink */






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


int ask_quark_source( int prompt, int *source_type, char *descrp)
{
  char savebuf[256];
  int status;

  if (prompt!=0)
    printf("enter 'point', or 'gaussian' for source type\n");
  status = scanf("%s",savebuf);
  if(status !=1) {
    printf("ask_quark_source: ERROR IN INPUT: source type command\n");
    return 1;
  }
  if(strcmp("point",savebuf) == 0 ){
    *source_type = POINT;
    strcpy(descrp,"point");
  }
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else{
    printf("ask_source: ERROR IN INPUT: source command is invalid\n"); 
    return 1;
  }

  printf("%s\n",savebuf);
  return 0;
} /* ask_quark_source */
