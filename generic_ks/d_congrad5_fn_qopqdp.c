NOT MAINTAINED!
/******* d_congrad5_fn_qopqdp.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter 
   for QDP format fields */
/* 2/2005 D. Renner and C. Jung */
/* 9/2005 C. DeTar modified for QDP/C */
/* 6/2006 NEEDS UPDATING */

#include "generic_ks_includes.h"
#include <qop_qdp.h>
#include <qop.h>

// Used by initialize_congrad and finalize_congrad:
static int is_congrad_initialized = 0;

void set_M_from_strided_parity_field(QDP_ColorMatrix *dest, su3_matrix *src,
				    int stride, int parity)
{
  int i;
  site *s;
  QLA_ColorMatrix *temp;
  temp = QDP_expose_M (dest);
  FORSOMEPARITY(i,s,parity) {
    memcpy((void *)&temp[i], (void *)&src[stride*i], sizeof(QLA_ColorMatrix));
  }
  QDP_reset_M (dest);
}

void set_V_from_parity_site(QDP_ColorVector *dest, field_offset src,
			     int parity)
{
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_expose_V (dest);
  FORSOMEPARITY(i,s,parity) {
    memcpy((void *)&temp[i], F_PT(s,src), sizeof(QLA_ColorVector));
  }
  QDP_reset_V (dest);
}

void set_parity_site_from_V(field_offset dest, QDP_ColorVector *src,
			     int parity)
{
  int i;
  site *s;
  QLA_ColorVector *temp;
  temp = QDP_expose_V (src);
  FORSOMEPARITY(i,s,parity) {
    memcpy(F_PT(s,dest), (void *)&temp[i], sizeof(QLA_##TYPE));
  }
  QDP_reset_V (src);
}

/* Generic MILC interface for the Asqtad inverter -- single mass */
/* (prec argument is ignored) */
int ks_congrad( field_offset milc_src, field_offset milc_sol, Real mass,
	        int niter, int nrestart, Real rsqmin, int prec, 
		int milc_parity, Real* final_rsq_ptr,
		ferm_links_t *fn)
{
  QDP_ColorMatrix *qop_fat_links[4];
  QDP_ColorMatrix *qop_long_links[4];
  QDP_ColorVector *qop_src;
  QDP_ColorVector *qop_sol;
  QOP_invert_arg qop_invert_arg;
  int iterations_used;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;

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

  t_fatlink = fn->fat;
  t_longlink = fn->lng;

#ifdef CGTIME
  dtimec = -dclock(); 
#endif

  // Initialize geometry variables
  initialize_congrad();

  ///////////////////////////////////////////////////////
  // map milc fields for src1 sol1 to qopqdp fields       //
  ///////////////////////////////////////////////////////

  FORALLUPDIR(dir){
    qop_fat_links[dir]  = QDP_create_M();
    set_M_from_strided_parity_field(qop_fat_links[dir], t_fatlink + dir, 
				   4, milc_parity);
    qop_long_links[dir] = QDP_create_M();
    set_M_from_strided_parity_field(qop_long_links[dir], t_longlink + dir, 
				   4, milc_parity);
  }
  qop_src = QDP_create_V();
  set_V_from_parity_site(qop_src, milc_src1, milc_parity);
  qop_sol = QDP_create_V();
  set_V_from_parity_site(qop_sol, milc_sol1, milc_parity);

  // For memory savings.  Links may need to be recomputed later.
  free_fn_links(fn);
  invalidate_all_ferm_links(fn);

  ///////////////////////////////////////////////////////
  // set qop_invert_arg                                //
  ///////////////////////////////////////////////////////

  congrad_fn_set_qop_invert_arg( & qop_invert_arg, niter, nrestart,
				 rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  iterations_used = ks_congrad_qopqdp( qop_src, qop_sol, qop_fat_links, 
		qop_long_links, mass, & qop_invert_arg, final_rsq_ptr );

  ///////////////////////////////////////////////////////
  // map qop fields to milc fields                     //
  ///////////////////////////////////////////////////////

  set_parity_site_from_V(milc_sol, qop_sol, milc_parity);

  ///////////////////////////////////////////////////////
  // free qop fields                                   //
  ///////////////////////////////////////////////////////

  FORALLUPDIR(dir){
    QDP_destroy_M( qop_fat_links[dir] );   qop_fat_links[dir] = NULL;
    QDP_destroy_M( qop_long_links[dir] );  qop_long_links[dir] = NULL;
  }
  QDP_destroy_V( qop_src );  qop_src = NULL;
  QDP_destroy_V( qop_sol );  qop_sol = NULL;

  #ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
    {
      printf("CONGRAD5(total): time = %e (fn_qopqdp) masses = 1 iters = %d mflops = %e\n", 
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
int ks_congrad_qopqdp( QDP_ColorVector *qop_src, QDP_ColorVector *qop_sol,
		       QDP_ColorMatrix *qop_fat_links[4], 
		       QDP_ColorMatrix *qop_long_links[4], Real mass,
		       QOP_invert_arg* qop_invert_arg, Real* final_rsq_ptr )
{

  //printf( "MILC: ks_congrad_qop called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  // This will only initialize congrad if it is not already done.
  initialize_congrad();

  QOP_asqtad_invert_load_links_qdp( qop_fat_links, qop_long_links );


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
      QOP_asqtad_inv_qdp( qop_invert_arg, qop_sol, qop_src );

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
	printf("CONGRAD5(level 3 only): time = %e (fn_qopqdp) masses = 1 iters = %d mflops = %e\n", 
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
  QOP_layout qop_layout;
  
  //printf( "MILC: initialize_congrad called in FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );

  if( is_congrad_initialized == 1 )
    {
      //printf( "MILC: congrad is initialized already\n" );
      //printf( "MILC: initialize_congrad finished: FILE %s at LINE %i\n", __FILE__, __LINE__ );  fflush( NULL );
      return;
    }

  ///////////////////////////////////////////////////////
  // set qop layout                                    //
  ///////////////////////////////////////////////////////
  
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

void congrad_fn_set_qop_invert_arg( QOP_invert_arg* qop_invert_arg, 
	    int max_iters, int restarts, Real min_resid_sq, int milc_parity )
{

  /* Set the inversion argument structure */
  qop_invert_arg->max_iter = max_iters;
  qop_invert_arg->restart  = restarts;
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
