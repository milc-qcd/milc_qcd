NOT MAINTAINED !
#ifndef __cplusplus
#error d_congrad5_fn_qop_two_src.c requires a c++ compiler
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

extern "C" {
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


static void map_milc_src_sol_to_qop( field_offset milc_src, field_offset milc_sol,
                             Float* qop_src, Float* qop_sol, int milc_parity );


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

  ///////////////////////////////////////////////////////
  // allocate qop fields                               //
  ///////////////////////////////////////////////////////
  
  // Float is the QOP precision.  It reduces to float if SINGLE
  // is defined and defaults to double otherwise.
  Float* qop_fat_links = NULL;
  Float* qop_long_links = NULL;
  Float* qop_src = NULL;
  Float* qop_sol = NULL;
  
  allocate_qop_fields( & qop_fat_links, & qop_long_links, 
		       & qop_src, & qop_sol );
  
  ///////////////////////////////////////////////////////
  // map milc fields for src1 sol1 to qop fields       //
  ///////////////////////////////////////////////////////
  
  // The milc_fat_links and milc_long_links are passed implicitly.
  map_milc_to_qop( milc_src1, milc_sol1, 
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
  
  set_qop_invert_arg( & qop_invert_arg, mass1, niter, rsqmin, milc_parity );
  
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
  
  map_qop_to_milc( qop_sol, milc_sol1, milc_parity );
  
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

  set_qop_invert_arg( & qop_invert_arg, mass2, niter, rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  iterations_used = ks_congrad_qop( qop_src, qop_sol, qop_fat_links, qop_long_links, & qop_invert_arg, final_rsq_ptr );

  total_iters += iterations_used;

  ///////////////////////////////////////////////////////
  // map qop fields to milc sol2 fields                //
  ///////////////////////////////////////////////////////

  map_qop_to_milc( qop_sol, milc_sol2, milc_parity );

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

  //printf( "MILC: ks_congrad finished in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
  return( iterations_used );
}

static void map_milc_src_sol_to_qop( 
   field_offset milc_src, 
   field_offset milc_sol, 
   Float* qop_src, 
   Float* qop_sol, 
   int milc_parity )
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
  
  // This loops over all the sub-lattice coordinates.
  for( int sub_lattice_t = 0 ; sub_lattice_t < sub_lattice_nt ; 
       sub_lattice_t ++ )
    {
      int lattice_t = corner_t + sub_lattice_t;
      
      for( int sub_lattice_z = 0 ; sub_lattice_z < sub_lattice_nz ; 
	   sub_lattice_z ++ )
	{
	  int lattice_z = corner_z + sub_lattice_z;
	  
	  for( int sub_lattice_y = 0 ; sub_lattice_y < sub_lattice_ny ; 
	       sub_lattice_y ++ )
	    {
	      int lattice_y = corner_y + sub_lattice_y;
	      
	      for( int sub_lattice_x = 0 ; sub_lattice_x < sub_lattice_nx ; 
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
		  
		  if( ( ( milc_parity == EVENANDODD ) ) ||
		      ( ( milc_parity == ODD        ) && 
			( site_variable->parity == ODD  ) )    )
		    {
		      Real* milc_odd_src = (Real*) F_PT( site_variable, 
							 milc_src );
		      Real* milc_odd_sol = (Real*) F_PT( site_variable, 
							 milc_sol );
		      
		      for( int c = 0; c < 6; c++ )
			{
			  qop_odd_src[ c ] = (Float) milc_odd_src[ c ];
			  qop_odd_sol[ c ] = (Float) milc_odd_sol[ c ];
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

