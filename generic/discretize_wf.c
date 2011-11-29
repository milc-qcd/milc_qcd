/*********************  discretize_wf.c *********************************/
/* MILC Version 7 

   Read a real ASCII wave function and interpolate.
   The wave function file has r, U(r) pairs, where r is in fm
   and U(r) is the reduced wave function.

   Rescale to lattice spacing a and store result in a MILC complex
   field over the entire specified time slice.

   C. DeTar 11/07/2008: Stolen from Canopy code from a file of the
   same name.

*/

#include "generic_includes.h"
#define ALL_T_SLICES -1

typedef struct {
  double x;
  double y;
} point;

typedef struct {
  long   size;
  point* e;
} point_vector;

static point_vector* new_point_vector ( long size )
{
  point_vector* self = (point_vector *) malloc ( sizeof (point_vector) );
  self -> size = size;
  self -> e = (point *) malloc ( size * sizeof (point) );

  return self;
}

static void delete_point_vector ( point_vector* self )
{
  if ( self )
    {
      self -> size = 0;

      if ( self -> e )
        free ( self -> e );

      self -> e = 0;
      free ( self );
    }
}

void resize_point_vector ( point_vector* self, long size )
{
  self -> e = (point *)realloc ( self -> e, size * sizeof (point) );
  self -> size = size;
}

static point_vector* read_ascii_points ( const char* filename )
{
  point p;
  long cnt;
  long size = 1000; // starting size
  point_vector* v = new_point_vector ( size );
  FILE* fp;
  int status;

  if(this_node == 0){
    fp = fopen ( filename, "r" );
    if ( fp == NULL )
      {
	printf ( "Open failed: %s\n", filename );
	terminate ( 1 );
      }
    
    cnt = 0;
    while ( (status = fscanf ( fp, "%lf%lf", &p.x, &p.y )) != EOF && status != 0 )
      {
	if ( cnt >= size )
	  {
	    size *= 2; // double the size
	    printf ( "resizing vector %ld\n", size );
	    resize_point_vector ( v, size );
	  }
	
	v -> e [cnt] = p;
	
	++cnt;
      }

    if(status != EOF){
      printf("read_ascii_points: format error reading %s\n", filename);
      terminate(1);
    }
    
    resize_point_vector ( v, cnt );
    
    fclose ( fp );
  }

  /* Broadcast to all the nodes */

  broadcast_bytes((char *)&(v -> size), sizeof(long));
  broadcast_bytes((char *)( v -> e), v->size * sizeof(point) );

  return v;
}

// compute psi (r) = U (r) / r
static void compute_wf ( point_vector* wf )
{
  long cnt;

  for ( cnt = 0; cnt < wf -> size; ++cnt )
    {
      wf -> e [cnt].y /= wf -> e [cnt].x;
    }
}



typedef struct {
  point lower;
  point upper;
} bracket;

// return bracket: bracket.lower.x <= x <= bracket.upper.x
// best for equally spaced x
static bracket search ( const double x, const point_vector* f )
{
  bracket b;
  double dx = f -> e [1].x - f -> e [0].x; // spacing

  double r = ( x - f -> e [0].x ) / dx;

  long low = (long) floor ( r );

  if ( low < 0 )
    low = 0;

  if ( low > -2 + f -> size )
    low = -2 + f -> size;

  b.lower = f -> e [low];
  b.upper = f -> e [low+1];

  return b;
}

static double interpolate ( const double r, const point_vector* f )
{
  bracket b = search ( r, f );

  // lagrange linear
  double y = ( r - b.lower.x ) * b.upper.y + ( b.upper.x - r ) * b.lower.y;

  return y / ( b.upper.x - b.lower.x );
}

static void setwf(complex *wf, int x, int y, int z, int t, Real f)
{
  /* Adjust to the periodic home cube */
  x = (x + nx) % nx; y = (y + ny) % ny; z = (z + nz) % nz;
  /* Store the value if this node has it */
  if( this_node == node_number(x,y,z,t) )
    {
      wf[ node_index(x,y,z,t) ] .real = f;
      wf[ node_index(x,y,z,t) ] .imag = 0.;
    }
}

static void discretize_wf ( double const a, int x0, int y0, int z0, int t0, 
			    const point_vector* R, complex* wf )
{
  double norm = 0;
  double r;
  Real f;
  long x, y, z;
  long mx, my, mz;
  int t, tmin, tmax;

  if(t0 == ALL_T_SLICES){ tmin = 0; tmax = nt-1; }
  else                  { tmin = t0; tmax = t0; }
  
  for ( z = 0; z <= nz / 2; ++z )
    {
      mz = z ? nz - z : 0;
      
      for ( y = 0; y <= ny / 2; ++y )
        {
          my = y ? ny - y : 0;
	  
          for ( x = 0; x <= nx / 2; ++x )
            {
              mx = x ? nx - x : 0;

              r =  x * x + y * y + z * z;
              r = a * sqrt ( r );
              f = interpolate ( r, R );
	      for ( t = tmin; t <= tmax; t++ )
		{
		  setwf ( wf,  x+x0,  y+y0,  z+z0, t, f );
		  setwf ( wf, mx+x0,  y+y0,  z+z0, t, f );
		  setwf ( wf,  x+x0, my+y0,  z+z0, t, f );
		  setwf ( wf,  x+x0,  y+y0, mz+z0, t, f );
		  setwf ( wf, mx+x0, my+y0,  z+z0, t, f );
		  setwf ( wf, mx+x0,  y+y0, mz+z0, t, f );
		  setwf ( wf,  x+x0, my+y0, mz+z0, t, f );
		  setwf ( wf, mx+x0, my+y0, mz+z0, t, f );
		}
	    }
        }
    }

  // integrate
  t = tmin;
  for ( z = 0; z < nz; ++z )
    for ( y = 0; y < ny; ++y )
      for ( x = 0; x < nx; ++x )
	{
	  if(node_number(x,y,z,t) == this_node){
	    complex tmp = wf[ node_index(x,y,z,t) ];
	    norm += tmp.real * tmp.real + tmp.imag * tmp.imag;
	  }
	}
  g_doublesum( &norm );

  norm = 1. / sqrt ( norm );

  // normalize
  for ( z = 0; z < nz; ++z )
    for ( y = 0; y < ny; ++y )
      for ( x = 0; x < nx; ++x )
	for ( t = tmin; t <= tmax; t++ )
	  {
	    if(node_number(x,y,z,t) == this_node){
	      complex* tmp = wf + node_index(x,y,z,t);
	      tmp -> real *= norm;
	      tmp -> imag *= norm;
	    }
	  }
}

/*---------------------------------------------------------------*/

/* t0 specifies the time slice.  (t0 = ALL_T_SLICES is for all t)

   a is the lattice spacing in fm

   wf_file is the file name

   wf is returned as a MILC complex field.  It must be allocated and
   zeroed first.
*/


void fnal_wavefunction(complex *wf, int x0, int y0, int z0, int t0, 
		       Real a, char wf_file[]){
  point_vector *radial;

  // read reduced wavefunction U (r)
  radial = read_ascii_points ( wf_file );

  // compute R ( r ) = U (r) / r
  compute_wf ( radial );

  discretize_wf ( a, x0, y0, z0, t0, radial, wf );

  delete_point_vector( radial );

}
