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

static void discretize_wf ( double const a, int stride, 
			    int x0, int y0, int z0, int t0, 
			    const point_vector* R, complex* wf )
{
  double norm = 0;
  double r;
  Real f;
  long x, y, z;
  long mx, my, mz;
  int t, tmin, tmax;

  if(stride <= 0){
    node0_printf("discretize_wf: illegal stride parameter %d\n", stride);
    return;
  }

  if(t0 == ALL_T_SLICES){ tmin = 0; tmax = nt-1; }
  else                  { tmin = t0; tmax = t0; }
  
  /* Assumes nx, ny, nz are even */
  for ( z = 0; z <= nz / 2; ++z )
    {
      if ( z % stride > 0 ) continue;
      mz = z ? nz - z : 0;
      
      for ( y = 0; y <= ny / 2; ++y )
        {
	  if ( y % stride > 0 ) continue;
          my = y ? ny - y : 0;
	  
          for ( x = 0; x <= nx / 2; ++x )
            {
	      if ( x % stride > 0 ) continue;
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

static complex *wf_cache = NULL;
static int x0_cache = 0;
static int y0_cache = 0;
static int z0_cache = 0;
static int t0_cache = 0;
static int stride_cache = 0;
static Real a_cache = 0.;
static char *wf_file_cache = NULL;

void copy_cached(complex *wf){
  if(wf_cache == NULL){
    node0_printf("discretize_wf: FATAL: Attempted to copy from a NULL cached field\n");
    terminate(1);
  }
  node0_printf("Using cached wf\n");
  copy_c_field(wf, wf_cache);
}

/* If all the parameters are the same, copy the wf from cache and return 0
   If there is no cache or at least one parameter differs, return nonzero */

static int copy_from_cache(complex *wf, int x0, int y0, int z0, int t0, 
			   int stride, Real a, const char wf_file[]){
  int status;

  if ( wf_file_cache == NULL )
    status = 1;
  else
    status = strcmp(wf_file, wf_file_cache);

  if ( status == 0 )
    status = x0 != x0_cache || y0 != y0_cache || z0 != z0_cache ||
      t0 != t0_cache || stride != stride_cache || a != a_cache;

  /* Copy wf from cache if parameters are the same */
  if ( status == 0 ) copy_cached(wf);

  return status;
}


static void save_to_cache(complex *wf, int x0, int y0, int z0, int t0,
			  int stride, Real a, const char wf_file[]){

  if ( wf_cache == NULL ) wf_cache = create_c_field();
  copy_c_field(wf_cache, wf);
  x0_cache = x0;
  y0_cache = y0;
  z0_cache = z0;
  t0_cache = t0;
  stride_cache = stride;
  a_cache = a;
  wf_file_cache = (char *)realloc(wf_file_cache, (strlen(wf_file)+1)*sizeof(char));
  strcpy(wf_file_cache, wf_file);
}

/*---------------------------------------------------------------*/

/* t0 specifies the time slice.  (t0 = ALL_T_SLICES is for all t)

   a is the lattice spacing in fm

   wf_file is the file name

   wf is returned as a MILC complex field.  It must be allocated and
   zeroed first.
*/

void fnal_wavefunction(complex *wf, int x0, int y0, int z0, int t0, 
		       int stride, Real a, const char wf_file[]){
  point_vector *radial;

  int status = copy_from_cache(wf, x0, y0, z0, t0, stride, a, wf_file);
  if ( status == 0 ) return;

  // read reduced wavefunction U (r)
  radial = read_ascii_points ( wf_file );

  // compute R ( r ) = U (r) / r
  compute_wf ( radial );

  discretize_wf ( a, stride, x0, y0, z0, t0, radial, wf );

  delete_point_vector( radial );

  if(status != 0)
    save_to_cache(wf, x0, y0, z0, t0, stride, a, wf_file);
}
