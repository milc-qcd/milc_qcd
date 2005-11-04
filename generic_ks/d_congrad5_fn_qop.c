/******* d_congrad5_fn_qop.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */
/* 2/2005 D. Renner and C. Jung */

#include "generic_ks_includes.h"
#include <qop.h>

static int *machine_dimensions;

// These values are set in initialize_congrad
// Dimension of sub lattice on this node
int sub_lattice_nx;
int sub_lattice_ny;
int sub_lattice_nz;
int sub_lattice_nt;

static int sub_lattice_volume;

// These values are set in initialize_congrad
// Logical mesh coordinates of this node
int machine_x;
int machine_y;
int machine_z;
int machine_t;

// Used by initialize_congrad and finalize_congrad:
static int is_congrad_initialized = 0;

/* Generic MILC interface for the Asqtad inverter */
int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, Real rsqmin, int milc_parity, Real* final_rsq_ptr )
{
  Real* qop_fat_links = NULL;
  Real* qop_long_links = NULL;
  Real* qop_src = NULL;
  Real* qop_sol = NULL;

  //printf( "MILC: ks_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  //printf( "MILC: level 3 conjugate gradient being used\n" );

  #ifdef CGTIME

    double dtimec;
    double nflop = 1187;
    if( milc_parity == EVENANDODD ) nflop *= 2;

  #endif

  ///////////////////////////////////////////////////////
  // load fat and long links                           //
  ///////////////////////////////////////////////////////

  if( valid_fatlinks  != 1 ) load_fatlinks();
  if( valid_longlinks != 1 ) load_longlinks();

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  // Initialize geometry variables
  initialize_congrad();

  ///////////////////////////////////////////////////////
  // allocate qop fields                               //
  ///////////////////////////////////////////////////////

  congrad_fn_allocate_qop_fields( & qop_fat_links, & qop_long_links, 
				 & qop_src, & qop_sol );

  ///////////////////////////////////////////////////////
  // map milc fields to qop fields                     //
  ///////////////////////////////////////////////////////

  // The milc_fat_links and milc_long_links are passed implicitly.
  congrad_fn_map_milc_to_qop_raw( milc_src, milc_sol, qop_fat_links, 
			  qop_long_links, qop_src, qop_sol, milc_parity );

  // For memory savings.  Links may need to be recomputed later.
  free_fatlinks();
  free_longlinks();
  valid_fatlinks = 0;
  valid_longlinks = 0;

  ///////////////////////////////////////////////////////
  // set qop_invert_arg                                //
  ///////////////////////////////////////////////////////

  QOP_invert_arg qop_invert_arg;

  congrad_fn_set_qop_invert_arg( & qop_invert_arg, mass, niter, 
				 rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  int iterations_used = ks_congrad_qop( qop_src, qop_sol, qop_fat_links, 
		qop_long_links, & qop_invert_arg, final_rsq_ptr );

  ///////////////////////////////////////////////////////
  // map qop fields to milc fields                     //
  ///////////////////////////////////////////////////////

  congrad_fn_map_qop_raw_to_milc( qop_sol, milc_sol, milc_parity );

  ///////////////////////////////////////////////////////
  // free qop fields                                   //
  ///////////////////////////////////////////////////////

  free( qop_fat_links );   qop_fat_links  = NULL;
  free( qop_long_links );  qop_long_links = NULL;
  free( qop_src );         qop_src        = NULL;
  free( qop_sol );         qop_sol        = NULL;

  #ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
    {
      printf("CONGRAD5(total): time = %e iters = %d mflops = %e\n", 
	     dtimec, iterations_used,
        (double)( nflop * volume * iterations_used / 
		  ( 1.0e6 * dtimec * numnodes() ) ) );
      fflush(stdout);
    }
  }
  #endif

  total_iters += iterations_used;

  //printf( "MILC: ks_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return( iterations_used );
}

/* Local QOP interface */
int ks_congrad_qop( Real* qop_src, Real* qop_sol,
		    Real* qop_fat_links, Real* qop_long_links,
		    QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr )
{

  //printf( "MILC: ks_congrad_qop called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // This will only initialize congrad if it is not already done.
  initialize_congrad();

  QOP_asqtad_invert_load_links_raw( qop_fat_links, qop_long_links );

  // QOP remaps the links again, so this copy isn't needed
  // free( qop_fat_links );   qop_fat_links  = NULL;
  // free( qop_long_links );  qop_long_links = NULL;

  #ifdef CGTIME

    double dtimec;
    double nflop = 1187;
    if( qop_invert_arg->evenodd == QOP_EVENODD ) nflop *= 2;
    dtimec = -dclock(); 

  #endif


  // This is hard wired to use at most 5 restarts.
  int number_restarts = 5;
  int iterations_used = 0;
  int restart;

  for( restart = 0; restart < number_restarts; restart ++ )
  {
    int iterations_used_this_restart = 
      QOP_asqtad_inv_raw( qop_invert_arg, qop_sol, qop_src );

    iterations_used += iterations_used_this_restart;
    
    if( iterations_used_this_restart < qop_invert_arg->max_iter )
      {
	break;
      }
    else
      {
	if( restart == ( number_restarts - 1 ) )
	  {
	    printf( "MILC: The number of restarts, %i, was saturated.\n", 
		    number_restarts );
	  }
      }
  }
  
#ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
      {
	printf("CONGRAD5(level 3 only): time = %e iters = %d mflops = %e\n", 
	       dtimec, iterations_used,
	       (double)( nflop * volume * iterations_used / 
			 ( 1.0e6 * dtimec * numnodes() ) ) );
	fflush(stdout);
      }
  }
#endif

  QOP_asqtad_invert_unload_links();

  // This will only close congrad if it is not already done.
  // Eventually it should be moved somewhere else.
  finalize_congrad();

  //printf( "MILC: ks_congrad_qop finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return( iterations_used );
}

void initialize_congrad( void )
{
  //printf( "MILC: initialize_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  if( is_congrad_initialized == 1 )
    {
      //printf( "MILC: congrad is initialized already\n" );
      //printf( "MILC: initialize_congrad finished: FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
      return;
    }

  machine_dimensions = get_logical_machine_dimensions();
  sub_lattice_nx = nx/machine_dimensions[XUP];
  sub_lattice_ny = ny/machine_dimensions[YUP];
  sub_lattice_nz = nz/machine_dimensions[ZUP];
  sub_lattice_nt = nt/machine_dimensions[TUP];
  sub_lattice_volume = sub_lattice_nx*sub_lattice_ny*
    sub_lattice_nz*sub_lattice_nt;

  machine_coordinates = get_logical_machine_coordinates();
  machine_x = machine_coordinates[XUP];
  machine_y = machine_coordinates[YUP];
  machine_z = machine_coordinates[ZUP];
  machine_t = machine_coordinates[TUP];

  if( sub_lattice_nt % 2 != 0 )
    {
      printf( "MILC: sub_lattice_nt must be even : FILE %s at LINE %i\n", 
	      __FILE__, __LINE__ );
      terminate( 0 );
    }
  
  ///////////////////////////////////////////////////////
  // set qop layout                                    //
  ///////////////////////////////////////////////////////
  
  QOP_layout qop_layout;
  
  // The number of dimensions is hardwired to 4.
  qop_layout.ndims = 4;
  
  // Set the sub-lattice length in each direction.
  qop_layout.sites[ 0 ] = sub_lattice_nx;
  qop_layout.sites[ 1 ] = sub_lattice_ny;
  qop_layout.sites[ 2 ] = sub_lattice_nz;
  qop_layout.sites[ 3 ] = sub_lattice_nt;
  
  // Set the boundary condition in each direction.
  qop_layout.bc[ 0 ] = 0;  // x direction
  qop_layout.bc[ 1 ] = 0;  // y direction
  qop_layout.bc[ 2 ] = 0;  // z direction
  qop_layout.bc[ 3 ] = 1;  // t direction
  
  QOP_asqtad_invert_init( & qop_layout );

  is_congrad_initialized = 1;

  //printf( "MILC: initialize_congrad finished: FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void finalize_congrad( void )
{
  //printf( "MILC: finalize_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  if( is_congrad_initialized == 0 )
  {
    return;
  }

  QOP_asqtad_invert_finalize();

  is_congrad_initialized = 0;

  //printf( "MILC: finalize_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void congrad_fn_allocate_qop_fields( Real** qop_fat_links, 
     Real** qop_long_links, Real** qop_src, Real** qop_sol )
{
  //printf( "MILC: congrad_fn_allocate_qop_fields called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  // This is a little paranoid.
  if( *qop_fat_links != NULL )
  {
    printf( "MILC: qop_fat_links should not be allocated twice in FILE %s at LINE %i\n", __FILE__, __LINE__ );
    fflush( NULL );
    terminate( 0 );
  }
  if( *qop_long_links != NULL )
    {
      printf( "MILC: qop_long_links should not be allocated twice in FILE %s at LINE %i\n", __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  if( *qop_src != NULL )
    {
      printf( "MILC: qop_src should not be allocated twice in FILE %s at LINE %i\n", __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  if( *qop_sol != NULL )
    {
      printf( "MILC: qop_sol should not be allocated twice in FILE %s at LINE %i\n", __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  
  // 6 = 3 * 2
  *qop_src = (Real*) malloc( sub_lattice_volume * 6 * sizeof(Real) );
  *qop_sol = (Real*) malloc( sub_lattice_volume * 6 * sizeof(Real) );
  // 72 = 4 * 3 * 3 * 2
  *qop_fat_links  = (Real*) malloc( sub_lattice_volume * 72 * sizeof(Real) );
  *qop_long_links = (Real*) malloc( sub_lattice_volume * 72 * sizeof(Real) );
  
  if( *qop_src == NULL )
    {
      printf( "MILC: failed to allocate qop_src in FILE %s at LINE %i\n", 
	      __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  if( *qop_sol == NULL )
    {
      printf( "MILC: failed to allocate qop_sol in FILE %s at LINE %i\n", 
	      __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  if( *qop_fat_links == NULL )
    {
      printf("MILC: failed to allocate qop_fat_links in FILE %s at LINE %i\n", 
	      __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  if( *qop_long_links == NULL )
    {
      printf("MILC: failed to allocate qop_long_links in FILE %s at LINE %i\n", 
	     __FILE__, __LINE__ );
      fflush( NULL );
      terminate( 0 );
    }
  
  // This is a little paranoid.
  //memset( (void*) *qop_fat_links,  0, sub_lattice_volume * 72 * sizeof(Real) );
  //memset( (void*) *qop_long_links, 0, sub_lattice_volume * 72 * sizeof(Real) );
  //memset( (void*) *qop_src, 0, sub_lattice_volume * 6 * sizeof(Real) );
  //memset( (void*) *qop_sol, 0, sub_lattice_volume * 6 * sizeof(Real) );
  
  //printf( "MILC: congrad_fn_allocate_qop_fields finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void congrad_fn_set_qop_invert_arg( QOP_invert_arg* qop_invert_arg, Real mass, 
			 int max_iters, Real min_resid_sq, int milc_parity )
{
  qop_invert_arg->mass     = mass;
  qop_invert_arg->max_iter = max_iters;
  qop_invert_arg->rsqmin   = min_resid_sq;
  
  switch( milc_parity )
    {
    case( EVEN ):  qop_invert_arg->evenodd = QOP_EVEN;     break;
    case( ODD ):   qop_invert_arg->evenodd = QOP_ODD;      break;
    default:       qop_invert_arg->evenodd = QOP_EVENODD;  break;
    }
  
  if( milc_parity != EVEN )
    {
      printf( "MILC: only even parity is allowed in congrad for now\n" );
      normal_exit( 0 );
    }
  
  return;
}

void congrad_fn_map_milc_to_qop_raw( field_offset milc_src, 
				     field_offset milc_sol, 
				     Real* qop_fat_links, 
				     Real* qop_long_links, 
				     Real* qop_src, Real* qop_sol, 
				     int milc_parity )
{
  //printf( "MILC: congrad_fn_map_milc_to_qop_raw called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  
  // These are the lattice coordinates of the point with sub-lattice coordinates (0,0,0,0).
  int corner_x = machine_x * sub_lattice_nx;
  int corner_y = machine_y * sub_lattice_ny;
  int corner_z = machine_z * sub_lattice_nz;
  int corner_t = machine_t * sub_lattice_nt;
  
  // The following are qop fields.
  
  Real* qop_even_src = qop_src;
  Real* qop_odd_src  = qop_src + 6*even_sites_on_node;
  
  Real* qop_even_sol = qop_sol;
  Real* qop_odd_sol  = qop_sol + 6*even_sites_on_node;
  
  Real* qop_fat_t_link = qop_fat_links;
  Real* qop_fat_x_link = qop_fat_t_link + 18 * sub_lattice_volume;
  Real* qop_fat_y_link = qop_fat_x_link + 18 * sub_lattice_volume;
  Real* qop_fat_z_link = qop_fat_y_link + 18 * sub_lattice_volume;
  
  Real* qop_long_t_link = qop_long_links;
  Real* qop_long_x_link = qop_long_t_link + 18 * sub_lattice_volume;
  Real* qop_long_y_link = qop_long_x_link + 18 * sub_lattice_volume;
  Real* qop_long_z_link = qop_long_y_link + 18 * sub_lattice_volume;
  
  int sub_lattice_t;
  int sub_lattice_z;
  int sub_lattice_y;
  int sub_lattice_x;
  int c;

  // This loops over all the sub-lattice coordinates.
  for( sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; sub_lattice_t ++ )
  {
    int lattice_t = corner_t + sub_lattice_t;

  for( sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; sub_lattice_z ++ )
  {
    int lattice_z = corner_z + sub_lattice_z;

  for( sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; sub_lattice_y ++ )
  {
    int lattice_y = corner_y + sub_lattice_y;

  for( sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; sub_lattice_x ++ )
  {
    int lattice_x = corner_x + sub_lattice_x;

    int site_index = node_index( lattice_x, lattice_y, lattice_z, lattice_t );
    site* site_variable = & lattice[ site_index ];

    ///////////////////////////////////////////////////////
    // remap even source and solution vectors            //
    ///////////////////////////////////////////////////////

    if( ( ( milc_parity == EVENANDODD )                                     ) ||
        ( ( milc_parity == EVEN       ) & ( site_variable->parity == EVEN ) )    )
    {
      Real* milc_even_src = (Real*) F_PT( site_variable, milc_src );
      Real* milc_even_sol = (Real*) F_PT( site_variable, milc_sol );

      for( c = 0; c < 6; c++ )
      {
        qop_even_src[ c ] = (Real) milc_even_src[ c ];
        qop_even_sol[ c ] = (Real) milc_even_sol[ c ];
      }

      qop_even_src += 6;
      qop_even_sol += 6;
    }

    ///////////////////////////////////////////////////////
    // remap odd source and solution vectors             //
    ///////////////////////////////////////////////////////

    if( ( milc_parity == EVENANDODD ) ||
	( ( milc_parity == ODD  ) && ( site_variable->parity == ODD  ) ) )
    {
      Real* milc_odd_src = (Real*) F_PT( site_variable, milc_src );
      Real* milc_odd_sol = (Real*) F_PT( site_variable, milc_sol );

      for( c = 0; c < 6; c++ )
      {
        qop_odd_src[ c ] = (Real) milc_odd_src[ c ];
        qop_odd_sol[ c ] = (Real) milc_odd_sol[ c ];
      }

      qop_odd_src += 6;
      qop_odd_sol += 6;
    }

    ///////////////////////////////////////////////////////
    // remap fat and long links                          //
    ///////////////////////////////////////////////////////

    {
      Real* milc_fat_x_link = NULL;
      Real* milc_fat_y_link = NULL;
      Real* milc_fat_z_link = NULL;
      Real* milc_fat_t_link = NULL;

      Real* milc_long_x_link = NULL;
      Real* milc_long_y_link = NULL;
      Real* milc_long_z_link = NULL;
      Real* milc_long_t_link = NULL;

      {
        int site_index_4 = site_index * 4;

        milc_fat_x_link = (Real*) &( t_fatlink[ 0 + site_index_4 ] );
        milc_fat_y_link = (Real*) &( t_fatlink[ 1 + site_index_4 ] );
        milc_fat_z_link = (Real*) &( t_fatlink[ 2 + site_index_4 ] );
        milc_fat_t_link = (Real*) &( t_fatlink[ 3 + site_index_4 ] );

        milc_long_x_link = (Real*) &( t_longlink[ 0 + site_index_4 ] );
        milc_long_y_link = (Real*) &( t_longlink[ 1 + site_index_4 ] );
        milc_long_z_link = (Real*) &( t_longlink[ 2 + site_index_4 ] );
        milc_long_t_link = (Real*) &( t_longlink[ 3 + site_index_4 ] );
      }

      for( c = 0; c < 18; c ++ )
      {
        qop_fat_x_link[ c ] = milc_fat_x_link[ c ];
        qop_fat_y_link[ c ] = milc_fat_y_link[ c ];
        qop_fat_z_link[ c ] = milc_fat_z_link[ c ];
        qop_fat_t_link[ c ] = milc_fat_t_link[ c ];

        qop_long_x_link[ c ] = milc_long_x_link[ c ];
        qop_long_y_link[ c ] = milc_long_y_link[ c ];
        qop_long_z_link[ c ] = milc_long_z_link[ c ];
        qop_long_t_link[ c ] = milc_long_t_link[ c ];
      }

      qop_fat_x_link += 18;
      qop_fat_y_link += 18;
      qop_fat_z_link += 18;
      qop_fat_t_link += 18;

      qop_long_x_link += 18;
      qop_long_y_link += 18;
      qop_long_z_link += 18;
      qop_long_t_link += 18;
    }
  }
  }
  }
  }

  //printf( "MILC: congrad_fn_map_milc_to_qop_raw finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void congrad_fn_map_qop_raw_to_milc( Real* qop_sol, 
				     field_offset milc_sol, int milc_parity )
{
  //printf( "MILC: congrad_fn_map_qop_raw_to_milc called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );


  // These are the lattice coordinates of the point with sub-lattice coordinates (0,0,0,0).
  int corner_x = machine_x * sub_lattice_nx;
  int corner_y = machine_y * sub_lattice_ny;
  int corner_z = machine_z * sub_lattice_nz;
  int corner_t = machine_t * sub_lattice_nt;
  int sub_lattice_t;
  int sub_lattice_z;
  int sub_lattice_y;
  int sub_lattice_x;
  int c;

  // The following are qop fields.
  Real* qop_even_sol = qop_sol;
  Real* qop_odd_sol  = qop_sol + 6*even_sites_on_node;

  // This loops over all the sub-lattice coordinates.
  for( sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; sub_lattice_t ++ )
  {
    int lattice_t = corner_t + sub_lattice_t;

  for( sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; sub_lattice_z ++ )
  {
    int lattice_z = corner_z + sub_lattice_z;

  for( sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; sub_lattice_y ++ )
  {
    int lattice_y = corner_y + sub_lattice_y;

  for( sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; sub_lattice_x ++ )
  {
    int lattice_x = corner_x + sub_lattice_x;

    site* site_variable = & lattice[ node_index( lattice_x, lattice_y, lattice_z, lattice_t ) ];

    ///////////////////////////////////////////////////////
    // remap even source and solution vectors            //
    ///////////////////////////////////////////////////////

    if( ( milc_parity == EVENANDODD ) ||
        ( ( milc_parity == EVEN     ) && ( site_variable->parity == EVEN ) ) )
    {
      Real* milc_even_sol = (Real*) F_PT( site_variable, milc_sol );

      for( c = 0; c < 6; c++ )
      {
        milc_even_sol[ c ] = (Real) qop_even_sol[ c ];
      }

      qop_even_sol += 6;
    }

    ///////////////////////////////////////////////////////
    // remap odd source and solution vectors             //
    ///////////////////////////////////////////////////////

    if( ( milc_parity == EVENANDODD ) ||
	( ( milc_parity == ODD      ) & ( site_variable->parity == ODD  ) ) )
    {
      Real* milc_odd_sol = (Real*) F_PT( site_variable, milc_sol );

      for( c = 0; c < 6; c++ )
      {
        milc_odd_sol[ c ] = (Real) qop_odd_sol[ c ];
      }

      qop_odd_sol += 6;
    }
  }
  }
  }
  }

  //printf( "MILC: congrad_fn_map_qop_raw_to_milc finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}


/***************************/
/* original MILC functions */
/***************************/

#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

/* clear an su3_vector in the lattice */
void clear_latvec(field_offset v,int parity){
register int i,j;
register site *s;
register su3_vector *vv;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case ODD: FORODDSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
    } 
}

/* copy an su3_vector in the lattice */
void copy_latvec(field_offset src,field_offset sol,int parity){
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		*dpt = *spt;
	    } break;
	case ODD: FORODDSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		*dpt = *spt;
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		*dpt = *spt;
	    } break;
    } 
}

/* scalar multiply and add an SU3 vector in the lattice */
void scalar_mult_add_latvec(field_offset src1,field_offset src2,
			    Real scalar,field_offset sol,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
               spt1 = (su3_vector *)F_PT(s,src1);
                spt2 = (su3_vector *)F_PT(s,src2);
                dpt = (su3_vector *)F_PT(s,sol);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV( (su3_vector *)F_PT((s+FETCH_UP),src1),
				(su3_vector *)F_PT((s+FETCH_UP),src2),
				(su3_vector *)F_PT((s+FETCH_UP),sol) );
		}
                scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt);
       } END_LOOP
}

void scalar2_mult_add_su3_vector(su3_vector *a, Real s1, su3_vector *b, 
				 Real s2, su3_vector *c){
register int i;
    for(i=0;i<3;i++){
        c->c[i].real = s1*a->c[i].real + s2*b->c[i].real;
        c->c[i].imag = s1*a->c[i].imag + s2*b->c[i].imag;
    }
}

/* scalar multiply two SU3 vector and add through the lattice */
void scalar2_mult_add_latvec(field_offset src1,Real scalar1,
			     field_offset src2,Real scalar2,
			     field_offset sol,int parity)
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
        FORSOMEPARITY(i,s,parity){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt  = (su3_vector *)F_PT(s,sol);
		if( i < loopend-FETCH_UP ){
		  prefetch_VVV((su3_vector *)F_PT((s+FETCH_UP),src1),
			       (su3_vector *)F_PT((s+FETCH_UP),src2),
			       (su3_vector *)F_PT((s+FETCH_UP),sol) );
		}
		scalar2_mult_add_su3_vector( spt1, scalar1, spt2, scalar2, dpt);
       } END_LOOP
}

/* scalar multiply an SU3 vector in the lattice */
void scalar_mult_latvec(field_offset src,Real scalar,
			field_offset sol,int parity)
{
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case ODD: FORODDSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,sol);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
    } 
}
