/******************* d_congrad5_fn_qphix.c ************************/
/* For the QPhiX interface */
/* MIMD version 7 */

#include "../include/generic_qphix.h"
#include "../include/generic_ks_qphix.h"
#include "../include/generic.h"
#include <lattice.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

/*! \brief call to the qphix_ks_congrad_parity.
 *
 * Choose the inversion precision
 */
int
ks_congrad_parity_qphix ( su3_vector *src
			, su3_vector *sol
			, quark_invert_control *qic
			, Real mass
			, fn_links_t *fn)			 
{
  int iterations_used;
  
  if(qic->prec == 1)
    iterations_used = 
      ks_congrad_parity_qphix_F( src, sol, qic, mass, fn );
  else
    iterations_used = 
      ks_congrad_parity_qphix_D( src, sol, qic, mass, fn );
  
  total_iters += iterations_used;
  return iterations_used;
}

#if 0

    int niter        = qic->max;      /* maximum number of iters per restart */
    int nrestart     = qic->nrestart; /* maximum restarts */
    Real rsqmin      = qic->resid * qic->resid;    /* desired residual - 
						      normalized as sqrt(r*r)/sqrt(src_e*src_e) */
    Real relrsqmin   = qic->relresid * qic->relresid; /* desired relative residual (FNAL)*/
    int prec         = qic->prec;     /* precision */
    int parity       = qic->parity;   /* EVEN, ODD */

    char *prec_label[2] = {"F", "D"};
    double ttime, dctime, tot_cg_time;
    double dtime;
    double nflop = 1187;
    int iters;
    int otherparity;
    struct QuarkInvertControl qphix_qic;

    assert(parity != EVENANDODD && "EVENANDODD not yet implemented");

#ifdef CGTIME
    tot_cg_time = -dclock();
#endif   
    /* Pack structure */
    qphix_qic.prec      = prec;       /* Currently ignored */
    qphix_qic.parity    = parity;
    qphix_qic.max       = niter;
    qphix_qic.nrestart  = nrestart;
    qphix_qic.resid     = rsqmin;
    qphix_qic.relresid  = 0;          /* Suppresses this test */
    
    /* Check if the mbench object has been created */
    initialize_qphix();
    
    /* \fixme - how to handle EVENANDODD? Get the opposite parity */
    if (parity == EVEN) 
        otherparity = ODD;
    else if (parity == ODD)
        otherparity = EVEN;
    else
        otherparity = parity; // BOTH should be EVENANDODD    
    
#ifdef CG_DEBUG
    dctime = -dclock();
#endif
   
 /* Data layout conversions */
#if UNOPTIMIZED_PACK_UNPACK
#warning using unoptimized gathers
#if CG_DEBUG
    double t_sp1, t_sp2, t_ll1, t_ll2, t_fl1, t_fl2;
    t_sp1 = -dclock();
#endif    
    get_ks_spinors_from_lattice(src, g_qphix_env_obj->ks_src1, parity);
#if CG_DEBUG
    t_sp1 += dclock();
    t_sp2 = -dclock();
#endif    
    get_ks_spinors_from_lattice(sol, g_qphix_env_obj->ks_dest1, parity);
#if CG_DEBUG
    t_sp2 += dclock();
    t_ll1 = -dclock(); 
#endif     
    get_longlinks_from_lattice(g_qphix_env_obj->gll[0], EVEN, fn);
#if CG_DEBUG
    t_ll1 += dclock(); 
    t_ll2 = -dclock();
#endif       
    get_longlinks_from_lattice(g_qphix_env_obj->gll[1], ODD, fn);
#if CG_DEBUG
    t_ll2 += dclock();
    t_fl1 = -dclock(); 
#endif       
    get_fatlinks_from_lattice(g_qphix_env_obj->gfl[0], EVEN, fn);
#if CG_DEBUG
    t_fl1 += dclock();
    t_fl2 = -dclock(); 
#endif       
    get_fatlinks_from_lattice(g_qphix_env_obj->gfl[1], ODD, fn);
#if CG_DEBUG
    t_fl2 += dclock();
    node0_printf("MILC-->QPhiX data layout conversion timings"
                 " (Unoptimized Gathers).\n"
                 "\t src-spinor  = %e\n"
                 "\t dest-spinor = %e\n"
                 "\t long-link1  = %e\n"
                 "\t long-link2  = %e\n"
                 "\t fat-link1   = %e\n"
                 "\t fat-link2   = %e\n"
                 "\t total       = %e\n"
                 , t_sp1, t_sp2, t_ll1, t_ll2, t_fl1, t_fl2
                 , t_sp1 + t_sp2 + t_ll1 + t_ll2 + t_fl1 + t_fl2
                );
    fflush(stdout);
#endif       
#else
#warning using optimized gather
#if CG_DEBUG
    double t_sp1, t_sp2, t_ll;
    t_sp1 = -dclock();
#endif    
    //gather_su3vectors_from_lattice(src_arg, g_qphix_env_obj->ks_src1, parity);
    get_ks_spinors_from_lattice(src, g_qphix_env_obj->ks_src1, parity);
#if CG_DEBUG
    t_sp1 += dclock();
    t_sp2 = -dclock();
#endif    
    get_ks_spinors_from_lattice(sol, g_qphix_env_obj->ks_dest1, parity);
    //gather_su3vectors_from_lattice(dest_arg, g_qphix_env_obj->ks_dest1, parity);
#if CG_DEBUG
    t_sp2 += dclock();
    t_ll = -dclock(); 
#endif         
    get_links_from_lattice(g_qphix_env_obj->gfl, g_qphix_env_obj->gll, fn);
#if CG_DEBUG
    t_ll += dclock();
    node0_printf("MILC-->QPhiX data layout conversion timings"
                 " (Optimized Gathers).\n"
                 "\t src-spinor  = %e\n"
                 "\t dest-spinor = %e\n"
                 "\t all-links   = %e\n"
                 "\t total   = %e\n"
                 , t_sp1, t_sp2, t_ll, t_sp1 + t_sp2 + t_ll
                );
#endif 
    
#endif

#ifdef CG_DEBUG
    dctime +=dclock();
    dtime = -dclock();
#endif    
    iters = qphix_ks_congrad_parity ( g_qphix_env_obj->ks_src1
                                    , g_qphix_env_obj->ks_dest1
                                    , &qphix_qic
                                    , mass
                                    , g_qphix_env_obj->gll
                                    , g_qphix_env_obj->gfl
                                    );
#ifdef CG_DEBUG    
    dtime += dclock();
#endif

    /* Unpack the results */
    qic->final_rsq    = qphix_qic.final_rsq;
    
#ifdef CG_DEBUG
    ttime = -dclock();
#endif
    /* Copy results back to su3_vector */
    set_ks_spinors_into_lattice (sol, g_qphix_env_obj->ks_dest1, parity);
    
#ifdef CG_DEBUG
    ttime +=dclock();
    dctime +=ttime;
#endif
#ifdef CGTIME
    tot_cg_time +=dclock();
    if(this_node==0) {
        node0_printf("CONGRAD5: total cg-time = %e "
#if CG_DEBUG
               "solve-time = %e "
               "layout-conversion-time = %e "
#endif            
               "(QPHIX %s) masses = 1 iters = "
               "%d mflops = %e "
#ifdef CG_DEBUG
               "mflops(ignore data-conv.) = %e"
#endif
               "\n"
               , tot_cg_time
#ifdef CG_DEBUG
               , dtime, dctime
#endif
               , prec_label[PRECISION-1], iters
               , (double)(nflop*volume*iters/(1.0e6*tot_cg_time*numnodes()))
#ifdef CG_DEBUG
               , (double)(nflop*volume*iters/(1.0e6*dtime*numnodes()))
#endif
               );
        fflush(stdout);
    }
#endif

    return iters;
}


#endif
