/******* d_congrad5_fn_qop_two_src.c - conjugate gradient for SU3/fermions **/
/* MIMD version 7 */

/* This is the two-source MILC wrapper for the SciDAC Level 3 QOP inverter 
   using the raw interface */
/* 2/2005 D. Renner and C. Jung */
/* 5/2005 C. DeTar two source version eliminates one remapping */
/* 9/2005 C. DeTar converted to C code */

#include "generic_ks_includes.h"

#include "/host/cdetar/qop/asqtad-2.6.0-CJ-8-16-05/include/qop.h"

// These values are set in initialize_congrad
// Dimension of sub lattice on this node
extern int sub_lattice_nx;
extern int sub_lattice_ny;
extern int sub_lattice_nz;
extern int sub_lattice_nt;

// These values are set in initialize_congrad
// Logical mesh coordinates of this node
extern int machine_x;
extern int machine_y;
extern int machine_z;
extern int machine_t;


static void map_milc_src_sol_to_qop( field_offset milc_src, 
				     field_offset milc_sol,
				     Real* qop_src, Real* qop_sol, 
				     int milc_parity );


int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset milc_src1,     /* source vector (type su3_vector) */
    field_offset milc_src2,
    field_offset milc_sol1,	/* solution vectors */
    field_offset milc_sol2,
    Real mass1,
    Real mass2,
    int niter,		        /* maximal number of CG interations */
    Real rsqmin,	        /* desired residue squared */
    int milc_parity,		/* parity to be worked on */
    Real  *final_rsq_ptr 	/* final residue squared */
    )
{
#ifdef CGTIME
  
  // Real is the QOP precision.  It reduces to float if SINGLE
  // is defined and defaults to double otherwise.
  Real* qop_fat_links = NULL;
  Real* qop_long_links = NULL;
  Real* qop_src = NULL;
  Real* qop_sol = NULL;
  
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
  // map milc fields for src1 sol1 to qop fields       //
  ///////////////////////////////////////////////////////
  
  // The milc_fat_links and milc_long_links are passed implicitly.
  congrad_fn_map_milc_to_qop_raw( milc_src1, milc_sol1, 
		   qop_fat_links, qop_long_links, 
		   qop_src, qop_sol, milc_parity );
  
  // For memory savings.  Links may need to be recomputed later.
  free_fatlinks();
  free_longlinks();
  valid_fatlinks = 0;
  valid_longlinks = 0;
  
  ///////////////////////////////////////////////////////
  // set qop_invert_arg                                //
  ///////////////////////////////////////////////////////
  
  QOP_invert_arg qop_invert_arg;
  
  congrad_fn_set_qop_invert_arg( & qop_invert_arg, mass1, niter, 
				 rsqmin, milc_parity );
  
  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////
  
  int iterations_used = ks_congrad_qop( qop_src, qop_sol, qop_fat_links, 
					qop_long_links, & qop_invert_arg, 
					final_rsq_ptr );
  
  total_iters += iterations_used;
  
  ///////////////////////////////////////////////////////
  // map qop field to milc sol1 field                  //
  ///////////////////////////////////////////////////////
  
  congrad_fn_map_qop_raw_to_milc( qop_sol, milc_sol1, milc_parity );
  
#ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
      {
	printf("CONGRAD5(total): time = %e iters = %d mflops = %e\n", 
	       dtimec, iterations_used,
	       (double)( nflop * volume * iterations_used / ( 1.0e6 * dtimec * numnodes() ) ) );
	fflush(stdout);
      }
    
    /* Restart the timing */
    dtimec = -dclock(); 
  }
  
#endif

  ///////////////////////////////////////////////////////
  // map milc fields for src2 sol2 to qop fields       //
  ///////////////////////////////////////////////////////

  map_milc_src_sol_to_qop( milc_src2, milc_sol2, qop_src, qop_sol, milc_parity );

  ///////////////////////////////////////////////////////
  // reset qop_invert_arg                                //
  ///////////////////////////////////////////////////////

  congrad_fn_set_qop_invert_arg( & qop_invert_arg, mass2, niter, 
				 rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  iterations_used = ks_congrad_qop( qop_src, qop_sol, qop_fat_links, 
			    qop_long_links, & qop_invert_arg, final_rsq_ptr );

  total_iters += iterations_used;

  ///////////////////////////////////////////////////////
  // map qop fields to milc sol2 fields                //
  ///////////////////////////////////////////////////////

  congrad_fn_map_qop_raw_to_milc( qop_sol, milc_sol2, milc_parity );

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

  //printf( "MILC: ks_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return( iterations_used );
}

static void map_milc_src_sol_to_qop( 
   field_offset milc_src, 
   field_offset milc_sol, 
   Real* qop_src, 
   Real* qop_sol, 
   int milc_parity )
{
  //printf( "MILC: map_milc_src_sol_to_qop called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  
  // These are the lattice coordinates of the point with sub-lattice
  // coordinates (0,0,0,0).
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
  
  Real* qop_even_src = qop_src;
  Real* qop_odd_src  = qop_src + 6*even_sites_on_node;
  
  Real* qop_even_sol = qop_sol;
  Real* qop_odd_sol  = qop_sol + 6*even_sites_on_node;
  
  // This loops over all the sub-lattice coordinates.
  for( sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; 
       sub_lattice_t ++ )
    {
      int lattice_t = corner_t + sub_lattice_t;
      
      for( sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; 
	   sub_lattice_z ++ )
	{
	  int lattice_z = corner_z + sub_lattice_z;
	  
	  for( sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; 
	       sub_lattice_y ++ )
	    {
	      int lattice_y = corner_y + sub_lattice_y;
	      
	      for( sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; 
		   sub_lattice_x ++ )
		{
		  int lattice_x = corner_x + sub_lattice_x;
		  
		  int site_index = node_index( lattice_x, lattice_y, lattice_z, lattice_t );
		  site* site_variable = & lattice[ site_index ];
		  
		  ///////////////////////////////////////////////////////
		  // remap even source and solution vectors            //
		  ///////////////////////////////////////////////////////
		  
		  if( ( ( milc_parity == EVENANDODD ) ) ||
		      ( ( milc_parity == EVEN       ) && 
			( site_variable->parity == EVEN ) )    )
		    {
		      Real* milc_even_src = (Real*) F_PT( site_variable, 
							  milc_src );
		      Real* milc_even_sol = (Real*) F_PT( site_variable, 
							  milc_sol );
		      
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
		  
		  if( ( ( milc_parity == EVENANDODD ) ) ||
		      ( ( milc_parity == ODD        ) && 
			( site_variable->parity == ODD  ) )    )
		    {
		      Real* milc_odd_src = (Real*) F_PT( site_variable, 
							 milc_src );
		      Real* milc_odd_sol = (Real*) F_PT( site_variable, 
							 milc_sol );
		      
		      for( c = 0; c < 6; c++ )
			{
			  qop_odd_src[ c ] = (Real) milc_odd_src[ c ];
			  qop_odd_sol[ c ] = (Real) milc_odd_sol[ c ];
			}
		      
		      qop_odd_src += 6;
		      qop_odd_sol += 6;
		    }
		}
	    }
	}
    }
  
  //printf( "MILC: map_milc_src_sol_to_qop finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return;
}

