/************************* remap_fftw_fourier.c **************************/
/* MIMD version 7 */

/* Fourier transform using remapping and scalar 1-D transforms */

/* At present this method supports only layouts for which each
   transformed direction is a divisor of the local volume. */

/* The data is remapped for each of the transformed directions, one
   dimension at a time, so that the entire dimension resides
   temporarily on a single processor while it is subjected to a 1D
   scalar FFT.  The resulting transform is then remapped to the
   original layout.
*/

#include "generic_includes.h"
#include <fftw3.h>

#if PRECISION==1
#define FFTWP(x) fftwf_##x
#else
#define FFTWP(x) fftw_##x
#endif

#ifdef CHECK_MALLOC

#define FFTWP(malloc)(_size) \
  (( (_malloc_ptr = FFTWP(malloc)(_size)), \
   (this_node == 0 ? \
   printf("%x = FFTWP(malloc)(%d) %s:%d\n",_malloc_ptr,_size,__func__,__LINE__) \
   && fflush(stdout) : 0 )), _malloc_ptr)

#endif

/* Data structure for the layout */

typedef struct {
  int *squaresize;
  int *nsquares;
  int *dirp;
  int (*node_number)(int *, int *, int *, int *);
  int (*node_index)(int *, int *, int *);
  void (*get_coords)(int *, int, int, int *, int *, int *);
  int nxfm;           /* The dimension of the transformed coordinate */
} ft_layout;

/* Data structure for the FT routine */

typedef struct {
  FFTWP(complex) *data;
  FFTWP(complex) *tmp;
  int size;           /* The size in bytes of the site data */
  int dir;
} ft_data;

/* List of possible prime factors for dividing the lattice, as in
 * layout_hyper_prime.c */

static int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 
		      43, 47, 53};
# define MAXPRIMES ( sizeof(prime) / sizeof(int) )

/* Table of layout map dirs from make_gather.  Filled in as they are created */
/* remap_dir[i][j] gives the map from layout i to layout j */
/* Layout index order is x,y,z,t,MILC=layout_*.c */

/* Even though it looks like you could just change dimensions here, the code
 * must still be changed below when you change this. */
#define NDIM 4
#define MILC_DIR NDIM

static int remap_dir[5][5] = { {NODIR, NODIR, NODIR, NODIR, NODIR},
			       {NODIR, NODIR, NODIR, NODIR, NODIR},
			       {NODIR, NODIR, NODIR, NODIR, NODIR},
			       {NODIR, NODIR, NODIR, NODIR, NODIR},
			       {NODIR, NODIR, NODIR, NODIR, NODIR}  };

/* layout[NDIM] is the MILC layout */
static ft_layout *layout[NDIM+1] = { NULL, NULL, NULL, NULL, NULL};

static int fwd_map[NDIM+1], bck_map[NDIM+1];
static FFTWP(plan) fwd_plan[NDIM], bck_plan[NDIM];

/*------------------------------------------------------------------*/
/* Layout functions for FT */
/*------------------------------------------------------------------*/

/* Convert rank to coordinates */
/* Here the lexicographic order of directions is specified by dirp */

static void lex_coords(const int dirp[], int coords[], const int dim, 
		       const int dims[], const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    coords[dirp[d]] = r % dims[dirp[d]];
    r /= dims[dirp[d]];
  }
}

/*------------------------------------------------------------------*/
/* The ft_node_number, ft_node_index, and ft_get_coords have the same
   purpose as the conventional layout routines with similar names.
   The FT layouts do not use EVENFIRST */

static int ft_node_number(int coords[], int dirp[], int squaresize[], 
			  int nsquares[]) {
  int i, dir;

  i = 0;
  for(dir = NDIM-1; dir >= 0; dir--)
    i = i*nsquares[dirp[dir]] + coords[dirp[dir]]/squaresize[dirp[dir]];

  return i;
}

/*------------------------------------------------------------------*/
static int ft_node_index(int coords[], int dirp[], int squaresize[]) {
  int i,dir;

  i = 0;
  for(dir = NDIM-1; dir >=0; dir--)
    i = i*squaresize[dirp[dir]] + coords[dirp[dir]]%squaresize[dirp[dir]];

  return i;
}

/*------------------------------------------------------------------*/
/* Map node number and index to coordinates  */
/* (The inverse of node_number and node_index) */

static void ft_get_coords(int coords[], int node, int index, 
			  int dirp[], int squaresize[], int nsquares[]){
  int mc[NDIM];
  int dir;

  /* coords = the sublattice coordinate */
  lex_coords(dirp, coords, NDIM, squaresize, index);

  /* mc = the machine coordinates for node k */
  lex_coords(dirp, mc, NDIM, nsquares, node);

  /* Add offset to get full lattice coordinate */
  for(dir = 0; dir < NDIM; dir++)
    coords[dir] += mc[dir]*squaresize[dir];
}

/*------------------------------------------------------------------*/
/* These are wrappers for the MILC layout routines.  We need versions
   that pretend to take the dirp, squaresize and nsquares arguments, but
   they are ignored */

static int milc_node_number(int coords[], int dirp[], int squaresize[], 
			    int nsquares[]) {
  return node_number(coords[0], coords[1], coords[2], coords[3]);
}

/*------------------------------------------------------------------*/
static int milc_node_index(int coords[], int dirp[], int squaresize[]) {
  return node_index(coords[0], coords[1], coords[2], coords[3]);
}

/*------------------------------------------------------------------*/
static void milc_get_coords(int coords[], int node, int index, 
			    int dirp[], int squaresize[], int nsquares[]){
  get_coords(coords, node, index); 
}

/*----------------------------------------------------------------------*/
/* Create the layout for an FFT in the direction ft_dir

   The layout is defined by squaresize, nsquares, and dirp.
   The first two have the same meaning as in standard MILC layout_*.c.
   The last is used to permute axes in the lexicographic ordering.

   We want this direction to be completely contained on a single node.
   Also, we want all data to be contiguous for a given value of the
   transform coordinate -- i.e. the transform coordinate should have
   slowest variation.  Since our layout functions make the fourth
   (time) coordinate vary slowest, we permute the lexicographic
   ordering of coordinates so the transform coordinate is last.

   Note that for this map to exist, the local volume must be divisible
   by the transform dimension.

   This algorithm is the same as that of
   layout_hyper_prime:setup_hyper_prime, except that we fix the
   squaresize for the transform dimension and do the divisions on the
   remaining dimensions.
*/

static int ft_setup_layout(int nsquares[], int squaresize[], int dirp[],
			   int ndim, int dims[], int ft_dir){
  int i,j,k,dir;
  char myname[] = "ft_setup_layout";
  int nodes = numnodes();

  /* Starting layout */
  for(dir = 0; dir < ndim; dir++){
    squaresize[dir] = dims[dir];
    nsquares[dir] = 1;
    dirp[dir] = dir;
  }
  
  /* dirp is used to redefine the lexicographic ordering of
     coordinates so the slowest varying (last one) is the FT
     coordinate */
  if(ndim-1 != ft_dir){
    k = dirp[ft_dir];
    dirp[ft_dir] = dirp[ndim-1];
    dirp[ndim-1] = k;
  }
  
  /* We don't divide the coordinate in the ft_dir */
  i = 1;	/* current number of hypercubes */
  while(i<nodes){
    /* figure out which prime to divide by starting with largest */
    k = MAXPRIMES-1;
    while( (nodes/i)%prime[k] != 0 && k>0 ) --k;
    /* figure out which direction to divide */
    
    /* find largest dimension of h-cubes divisible by prime[k] */
    for(j=0,dir=0;dir<=ndim-1;dir++)if(dir!=ft_dir)
      if( squaresize[dir]>j && squaresize[dir]%prime[k]==0 )
	j=squaresize[dir];
    
    /* if one direction with largest dimension has already been
       divided, divide it again.  Otherwise divide first direction
       with largest dimension. */
    for(dir=0;dir<=ndim-1;dir++)if(dir!=ft_dir)
      if( squaresize[dir]==j && nsquares[dir]>1 )break;
    if( dir >= ndim)for(dir=0;dir<=ndim-1;dir++)if(dir!=ft_dir)
      if( squaresize[dir]==j )break;
    /* This can fail if I run out of prime factors in the dimensions */
    if(dir >= ndim){
      if(mynode()==0)
	printf("%s: Can't lay out this lattice, not enough factors of %d\n",
	       myname,prime[k]);
      return 1;
    }
    
    /* do the surgery */
    i*=prime[k]; squaresize[dir] /= prime[k]; nsquares[dir] *= prime[k];
  }

  return 0;
}

/*----------------------------------------------------------------------*/
/* Create FT layout structure */

static ft_layout *ft_create_ft_layout(int ndim){
  ft_layout *ftl;
  
  ftl = (ft_layout *)malloc(sizeof(ft_layout));
  if(ftl == NULL){
    node0_printf("ft_create_ft_layout: No room\n");
    return NULL;
  }
  ftl->nsquares = (int *)malloc(ndim*sizeof(int));
  ftl->squaresize = (int *)malloc(ndim*sizeof(int));
  ftl->dirp = (int *)malloc(ndim*sizeof(int));

  if(ftl->nsquares == NULL || ftl->squaresize == NULL ||
     ftl->dirp == NULL){
    node0_printf("ft_create_ft_layout: No room\n");
    free(ftl);
    return NULL;
  }
  return ftl;
}

/*----------------------------------------------------------------------*/
/* Fill FT layout structure for MILC layout */

void ft_fill_milc_layout(int nsquares[], int squaresize[], 
			 int dirp[], int ndim, int dims[]){
  const int *milc_nsquares;
  int dir;

  /* Copy in the current MILC hypercubic layout dimensions */
  milc_nsquares = get_logical_dimensions();
  for(dir = 0; dir < ndim; dir++){
    nsquares[dir] = milc_nsquares[dir];
    squaresize[dir] = dims[dir]/milc_nsquares[dir];
    dirp[dir] = dir;
  }
}

/*----------------------------------------------------------------------*/
/* Make FFTW plans */

void make_fftw_plans(int size, ft_data *ftd){
  int ncmp;
  int rank, howmany, istride, idist, ostride, odist;
  int n[1], inembed[1], onembed[1];
  int nxfm;
  unsigned flags;
  int dir;
  //  double dtime = start_timing();

  flags = FFTW_ESTIMATE;  /* Could try FFTW_MEASURE */
  rank = 1;
  /* Number of complex values in a 4D site datum */
  ncmp = size/sizeof(complex);
  idist = odist = 1;

  for(dir = 0; dir < NDIM; dir++)
    if(layout[dir] != NULL){
      
      nxfm = layout[dir]->nxfm;
      
      /* The FT dimension */
      n[0] = inembed[0] = onembed[0] = nxfm;
      
      /* Number of contiguous complex values per 1D coordinate being
	 transformed */
      howmany = (sites_on_node*ncmp)/nxfm;
      ostride = istride = howmany;
      
      fwd_plan[dir] = 
	FFTWP(plan_many_dft)(rank, n, howmany, 
			    ftd->data, inembed, istride, idist, 
			    ftd->tmp, onembed, ostride, odist, 
			    FFTW_FORWARD, flags);
      bck_plan[dir] = 
	FFTWP(plan_many_dft)(rank, n, howmany, 
			    ftd->data, inembed, istride, idist, 
			    ftd->tmp, onembed, ostride, odist, 
			    FFTW_BACKWARD, flags);
    }

  // print_timing(dtime, "make FFTW plans");
}

/*----------------------------------------------------------------------*/
void destroy_fftw_plans(){
  int dir;
  for(dir = 0; dir < NDIM; dir++){
    FFTWP(destroy_plan)(fwd_plan[dir]);
    FFTWP(destroy_plan)(bck_plan[dir]);
  }
}


/*----------------------------------------------------------------------*/
/* Create all FT layouts needed */

void ft_create_layouts(ft_layout *ftl[], ft_layout **ft_milc, 
		       int ndim, int dims[], int key[]){
  int dir;
  //  int dtime = start_timing();

  /* Set up the FT layout structure for each dir needed */
  /* We don't remake the FT layout if it already exists, i.e.  the
     pointer is nonnull.  This provision allows a user to call
     setup_restrict_fourier multiple times without risking a memory
     leak */
  for(dir = 0; dir < ndim; dir++)
    if(ftl[dir] == NULL && key[dir] != 0){
      ftl[dir] = ft_create_ft_layout(ndim);
      ft_setup_layout(ftl[dir]->nsquares, ftl[dir]->squaresize,
		      ftl[dir]->dirp, ndim, dims, dir);
      ftl[dir]->node_number = ft_node_number;
      ftl[dir]->node_index = ft_node_index;
      ftl[dir]->get_coords = ft_get_coords;
      ftl[dir]->nxfm = dims[dir];
    }

  /* Set up the MILC layout structure */
  if(*ft_milc == NULL){
    *ft_milc = ft_create_ft_layout(ndim);
    /* Copy in the current MILC hypercubic layout dimensions and
       function pointers */
    ft_fill_milc_layout((*ft_milc)->nsquares, (*ft_milc)->squaresize, 
			(*ft_milc)->dirp, ndim, dims);
    (*ft_milc)->node_number = milc_node_number;
    (*ft_milc)->node_index = milc_node_index;
    (*ft_milc)->get_coords = milc_get_coords;
    (*ft_milc)->nxfm = 0;  /* We don't run transforms with the MILC layout */
  }

  //  print_timing(dtime, "create FFT layouts");
}

/*----------------------------------------------------------------------*/
/* Destroy FT layout */
void ft_destroy_ft_layout(ft_layout *ftl){
  if(ftl == NULL)return;
  if(ftl->nsquares != NULL){
    free(ftl->nsquares);
    ftl->nsquares = NULL;
  }
  if(ftl->squaresize != NULL){
    free(ftl->squaresize);
    ftl->squaresize = NULL;
  }
}

/*----------------------------------------------------------------------*/
/* Destroy all FT layouts */
void ft_destroy_ft_layouts(ft_layout **ft_milc, ft_layout *ftl[], int ndim){
  int dir;
  for(dir = 0; dir < ndim; dir++){
    ft_destroy_ft_layout(ftl[dir]);
    ftl[dir] = NULL;
  }
  ft_destroy_ft_layout(*ft_milc);
  *ft_milc = NULL;
}

/*----------------------------------------------------------------------*/
/* Remapping function for 4D lattices */
/* The mapping is defined so the initial and final coordinates
   have the same node number and index */

void ft_map_layouts(int x, int y, int z, int t, int *args, int fb, 
		    int *xp, int *yp, int *zp, int *tp){
  /* args[0] gives the final layout and args[1] the initial layout for
     forward mapping */
  ft_layout **ftl = (ft_layout **)args;
  ft_layout *ftl_src, *ftl_dst;
  int node, index;
  int coords[NDIM] = {x, y, z, t};

  if(fb == FORWARDS){
    ftl_dst = ftl[0];
    ftl_src = ftl[1];
  } else {  /* BACKWARDS */
    ftl_dst = ftl[1];
    ftl_src = ftl[0];
  }

  /* x, y, z, t are fake MILC coordinates that specify the MILC node
     and index to which we want to map */

  node = milc_node_number(coords, NULL, NULL, NULL);
  index = milc_node_index(coords, NULL, NULL);

  /* The actual coordinate that we want at this node and index based
     on the destination layout */

  ftl_dst->get_coords(coords, node, index, ftl_dst->dirp,
		      ftl_dst->squaresize, ftl_dst->nsquares);

  /* The node and index where the needed coordinate currently resides */
  
  node = ftl_src->node_number(coords,ftl_src->dirp, ftl_src->squaresize, 
			      ftl_src->nsquares);
  index = ftl_src->node_index(coords, ftl_src->dirp, ftl_src->squaresize);

  /* The fake MILC coordinate for the needed source node and index */

  milc_get_coords(coords, node, index, NULL, NULL, NULL);
  
  *xp = coords[0];
  *yp = coords[1];
  *zp = coords[2];
  *tp = coords[3];
}

/*----------------------------------------------------------------------*/
/* Make gather from layout to layout */

static int ft_make_gather(ft_layout *ftl_dst, ft_layout *ftl_src){
  ft_layout* ftl[2] = {ftl_dst, ftl_src};
  int *args = (int *)ftl;
  int dir;

  dir =  make_gather(ft_map_layouts, args, WANT_INVERSE,
		     ALLOW_EVEN_ODD, SCRAMBLE_PARITY);

  return dir;
}

/*----------------------------------------------------------------------*/
/* Make one FT map */

static void ft_make_map(int dirold, int dir){
  int index;

  fwd_map[dirold] = dir;
  bck_map[dir] = dirold;

  /* Make the gathers unless we already have them */
  if(remap_dir[dirold][dir] == NODIR ||
     remap_dir[dir][dirold] == NODIR){
    index = ft_make_gather(layout[dir],layout[dirold]);
    remap_dir[dirold][dir] = index;    /* forward map */
    remap_dir[dir][dirold] = index+1;  /* backward map */
  }
}
/*----------------------------------------------------------------------*/
/* Make all gather layouts needed for FT */

/* We do FT's in the order x, y, z, t for each key[dir] = 1 */
/* So we need to map from MILC to the first FT and from there
   to the second, etc, and finally back to MILC. */
/* We record the gather indexes in the matrix remap_dir */

static void ft_make_maps(ft_layout *ftl[], int key[], int ndim){

  int dir, dirold;
  //  double dtime = start_timing();

  dirold = MILC_DIR;  /* Start from MILC layout */
  for(dir = 0; dir < ndim; dir++){
    if(key[dir] != 0){
      ft_make_map(dirold, dir);
      dirold = dir;
    }
  }
  /* End with the MILC layout */
  dir = MILC_DIR;
  ft_make_map(dirold, dir);

  //  print_timing(dtime, "make FFTW gathers");
}

/*----------------------------------------------------------------------*/
/* This is the setup API.  It must be called before doing any FT */
/*----------------------------------------------------------------------*/
void setup_restrict_fourier( int *key, int *slice){
  /* "key" is a four component array.  If a component is 1, the Fourier
     transform is done in that direction, if it is 0 that direction is
     left alone. 
     If it is 2, then do the FT on a subset of the lattice with a
     fixed value slice[dir] of the coordinate in that direction. 
     "slice" is a four component array.  Not used unless key[dir]=2. */
  
  int dims[NDIM] = {nx, ny, nz, nt};
  int ndim = NDIM;
  int dir;
  
  /* No support for key[dir] = 2 */
  for(dir = 0; dir < ndim; dir++){
    if(key[dir] == 2){
      node0_printf("setup_restrict_fourier: No support for remapped slice FT's\n");
      terminate(1);
    }
  }

  /* Create the layouts for the 1D FT's */

  ft_create_layouts(layout, &layout[MILC_DIR], ndim, dims, key);

  /* Create the maps for switching layouts */

  ft_make_maps(layout, key, ndim);
}

/*----------------------------------------------------------------------*/

void cleanup_restrict_fourier(void){
  /* We have nothing to clean up */
}

/*----------------------------------------------------------------------*/
/* Copy from MILC field to FFTW field */

static void ft_copy_from_milc(FFTWP(complex) *data, complex *src, int size){
  int i,j;
  int ncmp = size/sizeof(complex);

  /* Copy data from src */
  for(i = 0; i < sites_on_node; i++)
    for(j = 0; j < ncmp; j++){
      data[i*ncmp+j][0] = src[i*ncmp+j].real;
      data[i*ncmp+j][1] = src[i*ncmp+j].imag;
    }
}

/*----------------------------------------------------------------------*/
/* Copy from FFTW field to MILC field */

static void ft_copy_to_milc(complex *src, FFTWP(complex) *data, int size){
  int i,j;
  int ncmp = size/sizeof(complex);

  /* Copy data from src */
  for(i = 0; i < sites_on_node; i++)
    for(j = 0; j < ncmp; j++){
      src[i*ncmp+j].real = data[i*ncmp+j][0];
      src[i*ncmp+j].imag = data[i*ncmp+j][1];
    }
}

/*----------------------------------------------------------------------*/
/* Create ft_data structure */

ft_data *create_ft_data(complex *src, int size){
  char myname[] = "create_ft_data";
  ft_data *ftd;
  int ncmp;
  //  double dtime = start_timing();

  ftd = (ft_data *)malloc(sizeof(ft_data));
  if(ftd == NULL){
    printf("%s: No room\n", myname);
    terminate(1);
  }

  ncmp = size/sizeof(complex);
  ftd->size = ncmp*sizeof(FFTWP(complex));
  ftd->dir = MILC_DIR;
  ftd->data 
    = (FFTWP(complex)*) FFTWP(malloc)(sizeof(FFTWP(complex))*sites_on_node*ncmp);
  ftd->tmp 
    = (FFTWP(complex)*) FFTWP(malloc)(sizeof(FFTWP(complex))*sites_on_node*ncmp);

  if(ftd->data == NULL || ftd->tmp == NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }
  /* Copy data in */
  ft_copy_from_milc(ftd->data, src, size);

  //  print_timing(dtime, "REMAP FFTW copy MILC");
  return ftd;
}

/*----------------------------------------------------------------------*/

void destroy_ft_data(ft_data *ftd){
  if(ftd != NULL){
    if(ftd->data != NULL)
      free(ftd->data);
    if(ftd->tmp != NULL)
      free(ftd->tmp);
    free(ftd);
  }
}

/*----------------------------------------------------------------------*/
/* Remap data according to map specified by "index" */

static void remap_data(int index, ft_data *ftd){
  msg_tag *mtag;
  char *temp;
  int i;
  //  double dtime = start_timing();

  temp = (char *)malloc(sites_on_node*ftd->size);
  if(temp==NULL){
    printf("remap_data: No room\n");
    terminate(1);
  }

  mtag = start_gather_field(ftd->data, ftd->size, index, EVENANDODD,
			    gen_pt[0]);
  wait_gather(mtag);

  /* First copy gathered data to temporary */
  for(i = 0; i < sites_on_node; i++)
    memcpy(temp + ftd->size*i, gen_pt[0][i], ftd->size);

  cleanup_gather(mtag);

  /* Then copy temp back to field */
  memcpy((char *)ftd->data, temp, sites_on_node*ftd->size);

  free(temp);

  //  print_timing(dtime, "REMAP FFTW remap");
}

/*----------------------------------------------------------------------*/
/* Return the next direction for a FT */

static int next_dir(int dirold, int isign){
  int dir;

  if(isign == 1)
    dir = fwd_map[dirold];
  else
    dir = bck_map[dirold];
  return dir;
}

/*----------------------------------------------------------------------*/
/* Return the direction of the last FT for the direction isign */

static int last_dir(int isign){

  /* The forward map is the reverse of the backward map
     so the last forward direction is the first one in the
     backward direction and vice versa.  We always start from
     the MILC_DIR. */

  if(isign == 1)
    return bck_map[MILC_DIR];
  else
    return fwd_map[MILC_DIR];
}

/*----------------------------------------------------------------------*/
/* Remap data for next step in chain */

static int remap_data_next(ft_data *ftd, int isign){
  int dirold, dirnew, index;

  dirold = ftd->dir;
  dirnew = next_dir(dirold, isign);

  index = remap_dir[dirold][dirnew];
  if(index == NODIR){
    printf("fourier_ftdata_alldir: Bad map %d to %d\n",dirold, dirnew);
    terminate(1);
  }
  remap_data(index, ftd);
  ftd->dir = dirnew;

  return dirnew;
}

/*----------------------------------------------------------------------*/
/* Do FT in one direction using FFTW */

void fourier_ftdata( ft_data *ftd, int isign ){

  //  double dtime = start_timing();

  if(isign == 1)
    FFTWP(execute)(fwd_plan[ftd->dir]);
  else
    FFTWP(execute)(bck_plan[ftd->dir]);

  // print_timing(dtime, "FFTW transform");
  //  dtime = start_timing();

  /* Copy the result from "tmp" back to "data" */

  memcpy((char *)ftd->data, (char *)ftd->tmp, ftd->size*sites_on_node);
  //  print_timing(dtime, "REMAP FFTW copy back");
}

/*----------------------------------------------------------------------*/
/* Do the whole FT, but stop before remapping to MILC  */

/* We do this so for convolutions we can avoid remapping to MILC by
   using this call instead of the standard FT call */

/* Pattern is

   MILC -> first -> next -> ... -> last -> MILC

   On entry we assume the state is "first" and untransformed. 
   On exit the state is "last" and transformed
*/

void fourier_ftdata_alldir(ft_data *ftd, int isign)
{

  int last = last_dir(isign);  /* Last direction before MILC */

  while( ftd->dir != last ){
    fourier_ftdata(ftd, isign);
    remap_data_next(ftd, isign);
  }

  fourier_ftdata(ftd, isign);
}

/*----------------------------------------------------------------------*/
void restrict_fourier_field(
  complex *src,	 /* src is field to be transformed */
  int size,	 /* Size of field in bytes.  The field must consist of
		    size/sizeof(complex) consecutive complex numbers.
		    For example, an su3_vector is 3 complex
		    numbers. */
   int isign)	 /* 1 for x -> k, -1 for k -> x */
{
  ft_data *ftd;

  /* Set up ft_data structure */

  ftd = create_ft_data(src, size);

  /* Create plans */
  make_fftw_plans(size, ftd);

  /* Map MILC to first FT dir */

  remap_data_next(ftd, isign);

  if(ftd->dir != MILC_DIR){
    /* Do FT, leaving data in FT layout for last direction */
    fourier_ftdata_alldir( ftd, isign );

    /* Map data back to MILC */
    remap_data_next(ftd, isign);
  }

  /* Copy result back to src */

  ft_copy_to_milc(src, ftd->data, size);

  destroy_fftw_plans();
  destroy_ft_data(ftd);
}

/*----------------------------------------------------------------------*/
void restrict_fourier_site(
     field_offset src,	 /* src is field to be transformed */
     int size,		 /* Size of field in bytes.  The field must
			    consist of size/sizeof(complex) consecutive
			    complex numbers.  For example, an su3_vector
			    is 3 complex numbers. */
     int isign)		 /* 1 for x -> k, -1 for k -> x */
{
  int i;
  site *s;
  complex *t_src;
  int ncomp = size/sizeof(complex);

  t_src  = (complex *) malloc(sites_on_node*size);

  if(t_src == NULL){
    printf("restrict_fourier_site(%d): Can't allocate src\n",this_node);
    terminate(1);
  }

  /* copy src to temporary */
  FORALLSITES(i,s) {
    memcpy( t_src + i*ncomp, F_PT(s,src), size );
  }

  restrict_fourier_field( t_src, size, isign);

  /* copy src back */
  FORALLSITES(i,s) {
    memcpy( F_PT(s,src), t_src + i*ncomp, size );
  }

  free(t_src);
}

/*----------------------------------------------------------------------*/
/* Multiply two FT's */

void mult_ft_by_ft(ft_data *dest, ft_data *src1, ft_data *src2){
  printf("mult_ft_by_ft not implemented\n");
}

