
#ifndef __cplusplus
#error d_congrad5_fn_qop.c requires a c++ compiler
#endif

// The file qop.h is a c++ header file.
#include "qop.h"

// The file generic_ks_includes.h is a c header.  This
// syntax may be dangerous.
extern "C" {
#include "generic_ks_includes.h"
}

// This requires that the choosen layout supplies the following global
// variables.  Current acceptable options for the layout are layout_qdp
// and layout_qcdoc.  The following list of globals can be made smaller
// if needed, but is convenient for now.

extern "C" int lattice_nx;  // This is another name for nx.
extern "C" int lattice_ny;  // This is another name for ny.
extern "C" int lattice_nz;  // This is another name for nz.
extern "C" int lattice_nt;  // This is another name for nt.

extern "C" int sub_lattice_nx;
extern "C" int sub_lattice_ny;
extern "C" int sub_lattice_nz;
extern "C" int sub_lattice_nt;

extern "C" int sub_lattice_volume;

extern "C" int machine_nx;
extern "C" int machine_ny;
extern "C" int machine_nz;
extern "C" int machine_nt;

extern "C" int machine_x;
extern "C" int machine_y;
extern "C" int machine_z;
extern "C" int machine_t;

// The function ks_congrad is a wrapper over the level 3, "qop",
// conjugate gradient calls.
int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, Real rsqmin, int parity, Real *final_rsq_ptr );

static int is_congrad_initialized = 0;  // This is used by initialize_congrad and finalize_congrad.

extern "C" {
void initialize_congrad( void );
void finalize_congrad( void );

void allocate_qop_fields( Float** qop_fat_links, Float** qop_long_links, Float** qop_src, Float** qop_sol );

void map_milc_to_qop( field_offset milc_src, field_offset milc_sol,
		      Float* qop_fat_links, Float* qop_long_links,
		      Float* qop_src, Float* qop_sol, int milc_parity );

void set_qop_invert_arg( QOP_invert_arg* qop_invert_arg, Real mass, int max_iterations, Real min_resid_sq, int milc_parity );

void map_qop_to_milc( Float* qop_sol, field_offset milc_sol, int milc_parity );

int ks_congrad_qop( Float* qop_source, Float* qop_solution,
                           Float* qop_fat_links, Float* qop_long_links,
                           QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr );
}

int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, Real rsqmin, int milc_parity, Real* final_rsq_ptr )
{
  //printf( "MILC: ks_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  //printf( "MILC: level 3 conjugate gradient being used\n" );

  #ifdef CGTIME

    double dtimec;
    double nflop = 1187;
    if( milc_parity == EVENANDODD ) nflop *= 2;
    dtimec = -dclock(); 

  #endif

  ///////////////////////////////////////////////////////
  // load fat and long links                           //
  ///////////////////////////////////////////////////////

  if( valid_fatlinks  != 1 ) load_fatlinks();
  if( valid_longlinks != 1 ) load_longlinks();

  ///////////////////////////////////////////////////////
  // allocate qop fields                               //
  ///////////////////////////////////////////////////////

  // Float is the QOP precision.  It reduces to float if SINGLE
  // is defined and defaults to double otherwise.
  Float* qop_fat_links = NULL;
  Float* qop_long_links = NULL;
  Float* qop_src = NULL;
  Float* qop_sol = NULL;

  allocate_qop_fields( & qop_fat_links, & qop_long_links, & qop_src, & qop_sol );

  ///////////////////////////////////////////////////////
  // map milc fields to qop fields                     //
  ///////////////////////////////////////////////////////

  // The milc_fat_links and milc_long_links are passed implicitly.
  map_milc_to_qop( milc_src, milc_sol, qop_fat_links, qop_long_links, qop_src, qop_sol, milc_parity );

  // For memory savings.  Links may need to be recomputed later.
  free_fatlinks();
  free_longlinks();
  valid_fatlinks = 0;
  valid_longlinks = 0;

  ///////////////////////////////////////////////////////
  // set qop_invert_arg                                //
  ///////////////////////////////////////////////////////

  QOP_invert_arg qop_invert_arg;

  set_qop_invert_arg( & qop_invert_arg, mass, niter, rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  int iterations_used = ks_congrad_qop( qop_src, qop_sol, qop_fat_links, qop_long_links, & qop_invert_arg, final_rsq_ptr );

  ///////////////////////////////////////////////////////
  // map qop fields to milc fields                     //
  ///////////////////////////////////////////////////////

  map_qop_to_milc( qop_sol, milc_sol, milc_parity );

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
      printf("CONGRAD5(total): time = %e iters = %d mflops = %e\n", dtimec, iterations_used,
        (double)( nflop * volume * iterations_used / ( 1.0e6 * dtimec * numnodes() ) ) );
      fflush(stdout);
    }
  }
  #endif

  total_iters += iterations_used;

  //printf( "MILC: ks_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return( iterations_used );
}

int ks_congrad_qop( Float* qop_src, Float* qop_sol,
                           Float* qop_fat_links, Float* qop_long_links,
                           QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr )
{
  //printf( "MILC: ks_congrad_qop called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // This will only initialize congrad if it is not already done.
  initialize_congrad();

  //printf( "MILC: QOP_asqtad_invert_load_links_raw called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  QOP_asqtad_invert_load_links_raw( qop_fat_links, qop_long_links );
  //printf( "MILC: QOP_asqtad_invert_load_links_raw finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

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
  for( int restart = 0; restart < number_restarts; restart ++ )
  {
    //printf( "MILC: QOP_asqtad_inv_raw called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
    int iterations_used_this_restart = QOP_asqtad_inv_raw( qop_invert_arg, qop_sol, qop_src );
    //printf( "MILC: QOP_asqtad_inv_raw finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

    iterations_used += iterations_used_this_restart;

    if( iterations_used_this_restart < qop_invert_arg->max_iter )
    {
      break;
    }
    else
    {
      if( restart == ( number_restarts - 1 ) )
      {
        printf( "MILC: The number of restarts, %i, was saturated.\n", number_restarts );
      }
    }
  }

  #ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
    {
      printf("CONGRAD5(level 3 only): time = %e iters = %d mflops = %e\n", dtimec, iterations_used,
        (double)( nflop * volume * iterations_used / ( 1.0e6 * dtimec * numnodes() ) ) );
      fflush(stdout);
    }
  }
  #endif

  //printf( "MILC: QOP_asqtad_invert_unload_links called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  QOP_asqtad_invert_unload_links();
  //printf( "MILC: QOP_asqtad_invert_unload_links finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // This will only close congrad if it is not already done.  Eventually it should be
  // moved somewhere else.
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

  if( sub_lattice_nt % 2 != 0 )
  {
    printf( "MILC: sub_lattice_nt must be even : FILE %s at LINE %i\n", __FILE__, __LINE__ );
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

  //printf( "MILC: QOP_asqtad_invert_init( & qop_layout ) called\n" );  fflush( NULL );
  QOP_asqtad_invert_init( & qop_layout );
  //printf( "MILC: QOP_asqtad_invert_init( & qop_layout ) finished\n" );  fflush( NULL );

  is_congrad_initialized = 1;

  //printf( "MILC: initialize_congrad finished: FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void finalize_congrad( void )
{
  //printf( "MILC: finalize_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  if( is_congrad_initialized == 0 )
  {
    //printf( "MILC: congrad is finalized already\n" );
    //printf( "MILC: finalize_congrad finished: FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
    return;
  }

  //printf( "MILC: QOP_asqtad_invert_finalize() called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  QOP_asqtad_invert_finalize();
  //printf( "MILC: QOP_asqtad_invert_finalize() finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  is_congrad_initialized = 0;

  //printf( "MILC: finalize_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void allocate_qop_fields( Float** qop_fat_links, Float** qop_long_links, Float** qop_src, Float** qop_sol )
{
  //printf( "MILC: allocate_qop_fields called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

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
  *qop_src = (Float*) malloc( sub_lattice_volume * 6 * sizeof(Float) );
  *qop_sol = (Float*) malloc( sub_lattice_volume * 6 * sizeof(Float) );
  // 72 = 4 * 3 * 3 * 2
  *qop_fat_links  = (Float*) malloc( sub_lattice_volume * 72 * sizeof(Float) );
  *qop_long_links = (Float*) malloc( sub_lattice_volume * 72 * sizeof(Float) );

  if( *qop_src == NULL )
  {
    printf( "MILC: failed to allocate qop_src in FILE %s at LINE %i\n", __FILE__, __LINE__ );
    fflush( NULL );
    terminate( 0 );
  }
  if( *qop_sol == NULL )
  {
    printf( "MILC: failed to allocate qop_sol in FILE %s at LINE %i\n", __FILE__, __LINE__ );
    fflush( NULL );
    terminate( 0 );
  }
  if( *qop_fat_links == NULL )
  {
    printf( "MILC: failed to allocate qop_fat_links in FILE %s at LINE %i\n", __FILE__, __LINE__ );
    fflush( NULL );
    terminate( 0 );
  }
  if( *qop_long_links == NULL )
  {
    printf( "MILC: failed to allocate qop_long_links in FILE %s at LINE %i\n", __FILE__, __LINE__ );
    fflush( NULL );
    terminate( 0 );
  }

  // This is a little paranoid.
  //memset( (void*) *qop_fat_links,  0, sub_lattice_volume * 72 * sizeof(Float) );
  //memset( (void*) *qop_long_links, 0, sub_lattice_volume * 72 * sizeof(Float) );
  //memset( (void*) *qop_src, 0, sub_lattice_volume * 6 * sizeof(Float) );
  //memset( (void*) *qop_sol, 0, sub_lattice_volume * 6 * sizeof(Float) );

  //printf( "MILC: allocate_qop_fields finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void set_qop_invert_arg( QOP_invert_arg* qop_invert_arg, Real mass, int max_iters, Real min_resid_sq, int milc_parity )
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

void map_milc_to_qop( field_offset milc_src, field_offset milc_sol, Float* qop_fat_links, Float* qop_long_links, Float* qop_src, Float* qop_sol, int milc_parity )
{
  //printf( "MILC: map_milc_to_qop called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // These are the lattice coordinates of the point with sub-lattice coordinates (0,0,0,0).
  int corner_x = machine_x * sub_lattice_nx;
  int corner_y = machine_y * sub_lattice_ny;
  int corner_z = machine_z * sub_lattice_nz;
  int corner_t = machine_t * sub_lattice_nt;

  // The following are qop fields.

  Float* qop_even_src = qop_src;
  Float* qop_odd_src  = qop_src + 6*even_sites_on_node;

  Float* qop_even_sol = qop_sol;
  Float* qop_odd_sol  = qop_sol + 6*even_sites_on_node;

  Float* qop_fat_t_link = qop_fat_links;
  Float* qop_fat_x_link = qop_fat_t_link + 18 * sub_lattice_volume;
  Float* qop_fat_y_link = qop_fat_x_link + 18 * sub_lattice_volume;
  Float* qop_fat_z_link = qop_fat_y_link + 18 * sub_lattice_volume;

  Float* qop_long_t_link = qop_long_links;
  Float* qop_long_x_link = qop_long_t_link + 18 * sub_lattice_volume;
  Float* qop_long_y_link = qop_long_x_link + 18 * sub_lattice_volume;
  Float* qop_long_z_link = qop_long_y_link + 18 * sub_lattice_volume;

  // This loops over all the sub-lattice coordinates.
  for( int sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; sub_lattice_t ++ )
  {
    int lattice_t = corner_t + sub_lattice_t;

  for( int sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; sub_lattice_z ++ )
  {
    int lattice_z = corner_z + sub_lattice_z;

  for( int sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; sub_lattice_y ++ )
  {
    int lattice_y = corner_y + sub_lattice_y;

  for( int sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; sub_lattice_x ++ )
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

      for( int c = 0; c < 6; c++ )
      {
        qop_even_src[ c ] = (Float) milc_even_src[ c ];
        qop_even_sol[ c ] = (Float) milc_even_sol[ c ];
      }

      qop_even_src += 6;
      qop_even_sol += 6;
    }

    ///////////////////////////////////////////////////////
    // remap odd source and solution vectors             //
    ///////////////////////////////////////////////////////

    if( ( ( milc_parity == EVENANDODD )                                     ) ||
	( ( milc_parity == ODD        ) & ( site_variable->parity == ODD  ) )    )
    {
      Real* milc_odd_src = (Real*) F_PT( site_variable, milc_src );
      Real* milc_odd_sol = (Real*) F_PT( site_variable, milc_sol );

      for( int c = 0; c < 6; c++ )
      {
        qop_odd_src[ c ] = (Float) milc_odd_src[ c ];
        qop_odd_sol[ c ] = (Float) milc_odd_sol[ c ];
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

      for( int c = 0; c < 18; c ++ )
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

  //printf( "MILC: map_milc_to_qop finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

void map_qop_to_milc( Float* qop_sol, field_offset milc_sol, int milc_parity )
{
  //printf( "MILC: map_qop_to_milc called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // These are the lattice coordinates of the point with sub-lattice coordinates (0,0,0,0).
  int corner_x = machine_x * sub_lattice_nx;
  int corner_y = machine_y * sub_lattice_ny;
  int corner_z = machine_z * sub_lattice_nz;
  int corner_t = machine_t * sub_lattice_nt;

  // The following are qop fields.
  Float* qop_even_sol = qop_sol;
  Float* qop_odd_sol  = qop_sol + 6*even_sites_on_node;

  // This loops over all the sub-lattice coordinates.
  for( int sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; sub_lattice_t ++ )
  {
    int lattice_t = corner_t + sub_lattice_t;

  for( int sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; sub_lattice_z ++ )
  {
    int lattice_z = corner_z + sub_lattice_z;

  for( int sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; sub_lattice_y ++ )
  {
    int lattice_y = corner_y + sub_lattice_y;

  for( int sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; sub_lattice_x ++ )
  {
    int lattice_x = corner_x + sub_lattice_x;

    site* site_variable = & lattice[ node_index( lattice_x, lattice_y, lattice_z, lattice_t ) ];

    ///////////////////////////////////////////////////////
    // remap even source and solution vectors            //
    ///////////////////////////////////////////////////////

    if( ( ( milc_parity == EVENANDODD )                                     ) ||
        ( ( milc_parity == EVEN       ) & ( site_variable->parity == EVEN ) )    )
    {
      Real* milc_even_sol = (Real*) F_PT( site_variable, milc_sol );

      for( int c = 0; c < 6; c++ )
      {
        milc_even_sol[ c ] = (Real) qop_even_sol[ c ];
      }

      qop_even_sol += 6;
    }

    ///////////////////////////////////////////////////////
    // remap odd source and solution vectors             //
    ///////////////////////////////////////////////////////

    if( ( ( milc_parity == EVENANDODD )                                     ) ||
	( ( milc_parity == ODD        ) & ( site_variable->parity == ODD  ) )    )
    {
      Real* milc_odd_sol = (Real*) F_PT( site_variable, milc_sol );

      for( int c = 0; c < 6; c++ )
      {
        milc_odd_sol[ c ] = (Real) qop_odd_sol[ c ];
      }

      qop_odd_sol += 6;
    }
  }
  }
  }
  }

  //printf( "MILC: map_qop_to_milc finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

