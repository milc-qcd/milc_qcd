NOT MAINTAINED!
/******* d_congrad5_fn_qop_two_src.c - conjugate gradient for SU3/fermions **/
/* MIMD version 7 */

/* This is the two-source MILC wrapper for the SciDAC Level 3 QOP inverter 
   using the QDP interface */
/* 2/2005 D. Renner and C. Jung */
/* 5/2005 C. DeTar two source version eliminates one remapping */
/* 9/2005 C. DeTar converted to C code */
/* 9/2005 C. DeTar modified for QDP/C */
/* 6/2006 NEEDS UPDATING */

#include "generic_ks_includes.h"

#include <qopqdp.h>

/* Standard MILC interface for two-source inversion */
int ks_congrad_two_src(	/* Return value is number of iterations taken */
    field_offset milc_src1,     /* source vector (type su3_vector) */
    field_offset milc_src2,
    field_offset milc_sol1,	/* solution vectors */
    field_offset milc_sol2,
    Real mass1,
    Real mass2,
    int niter,		        /* maximal number of CG interations */
    int nrestart,	        /* maximal number of CG restarts */
    Real rsqmin,	        /* desired residue squared */
    int prec,                   /* internal precision for the inversion 
				   (ignored) */
    int milc_parity,		/* parity to be worked on */
    Real  *final_rsq_ptr, 	/* final residue squared */
    ferm_links_t *fn             /* Storage for fermion links */
    )
{
#ifdef CGTIME
  
  QDP_ColorMatrix *qop_fat_links[4];
  QDP_ColorMatrix *qop_long_links[4];
  QDP_ColorVector *qop_src;
  QDP_ColorVector *qop_sol;
  int dir;
  double dtimec;
  double nflop = 1187;
  QOP_invert_arg qop_invert_arg;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;

  if( milc_parity == EVENANDODD ) nflop *= 2;
  
#endif
  
  ///////////////////////////////////////////////////////
  // load fat and long links                           //
  ///////////////////////////////////////////////////////
  
  t_fatlink = fn->fl.fat;
  t_longlink = fn->fl.lng;
  
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

  // Freeing memory is not necessary, but done for memory savings.
  // Links must be recomputed later.
  free_fn_links(fn);
  invalidate_all_ferm_links(fn);
  
  ///////////////////////////////////////////////////////
  // set qop_invert_arg                                //
  ///////////////////////////////////////////////////////
  
  congrad_fn_set_qop_invert_arg( & qop_invert_arg, mass1, niter, 
				 nrestart, rsqmin, milc_parity );
  
  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////
  
  int iterations_used = ks_congrad_qopqdp( qop_src, qop_sol, qop_fat_links, 
					   qop_long_links, & qop_invert_arg, 
					   final_rsq_ptr );
  
  total_iters += iterations_used;
  
  ///////////////////////////////////////////////////////
  // map qop field to milc sol1 field                  //
  ///////////////////////////////////////////////////////
  
  set_parity_site_from_V(milc_sol1, qop_sol, milc_parity);
  
#ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
      {
	printf("CONGRAD5(total): time = %e (qopqdp_two_src) masses = 1 iters = %d mflops = %e\n", 
	       dtimec, iterations_used,
	       (double)( nflop * volume * iterations_used / ( 1.0e6 * dtimec * numnodes() ) ) );
	fflush(stdout);
      }
    
    /* Restart the timing */
    dtimec = -dclock(); 
  }
  
#endif

  ////////////////////////////////////////////////////////
  // map second milc fields for src2 sol2 to qop fields //
  ////////////////////////////////////////////////////////

  set_V_from_parity_site( qop_src, milc_src2,  milc_parity );
  set_V_from_parity_site( qop_sol, milc_sol2,  milc_parity );

  ///////////////////////////////////////////////////////
  // reset qop_invert_arg                                //
  ///////////////////////////////////////////////////////

  congrad_fn_set_qop_invert_arg( & qop_invert_arg, mass2, niter, nrestart,
				 rsqmin, milc_parity );

  ///////////////////////////////////////////////////////
  // qop conjugate gradient                            //
  ///////////////////////////////////////////////////////

  iterations_used = ks_congrad_qopqdp( qop_src, qop_sol, qop_fat_links, 
			    qop_long_links, & qop_invert_arg, final_rsq_ptr );

  total_iters += iterations_used;

  ///////////////////////////////////////////////////////
  // map second qop solution to milc                   //
  ///////////////////////////////////////////////////////

  set_parity_site_from_V(milc_sol2, qop_sol, milc_parity);

  ///////////////////////////////////////////////////////
  // free qop fields                                   //
  ///////////////////////////////////////////////////////

  FORALLUPDIR(dir){
    QDP_destroy_M( qop_fat_links[dir] );
    QDP_destroy_M( qop_long_links );    
    qop_long_links = NULL;
  }
  QDP_destroy_V( qop_src ); qop_src = NULL;
  QDP_destroy_V( qop_sol ); qop_sol = NULL;

  #ifdef CGTIME
  {
    dtimec += dclock();
    if( this_node == 0 )
    {
      printf("CONGRAD5(total): time = %e (qopqdp_two_src) masses = 1 iters = %d mflops = %e\n", 
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


