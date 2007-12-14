/******* d_bicgilu_cl_qop_P.c - BiCG ILU ***********************/
/* MIMD version 7 */

/* This is the MILC wrapper for the SciDAC Level 3 QOP inverter */

/* NOTE: This code is actually an include file for d_bicgilu_cl_F.c
   and d_bicgilu_cl_D.c, so any edits should be consistent with this
   purpose. 

   NOTE: Several local procedures are defined "static".  That makes it
   possible to use the same name without confusion in the "D" and "F"
   variants of the compilation.  If you need to make them global
   externals, they must be given precision-specific names as we did
   for BICGILU_MILC2QOP.  Otherwise, the linker will find that they
   are multiply defined.

*/

/* Precision-dependent types

   NYREAL
   MYSU3_MATRIX
*/

#include "generic_clover_includes.h"
#include "../include/generic_qop.h"
#include "../include/generic_clover_qop.h"

/*#define CGDEBUG*/

/* Redefinitions according to requested precision */

#if ( QOP_Precision == 1 )

#define BGCGILU_MILC2QOP   bicgilu_cl_milc2qop_F
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix

#else

#define BGCGILU_MILC2QOP   bicgilu_cl_milc2qop_D
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix

#endif

/* 4/29/07 C. DeTar */

/*
 * $Log: d_bicgilu_cl_qop_P.c,v $
 * Revision 1.2  2007/12/14 04:37:48  detar
 * Add final restart member to qic.
 *
 * Revision 1.1  2007/05/21 04:43:58  detar
 * Support Level 3 inversion
 *
 */

#ifdef CGTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_clover/d_bicgilu_cl_qop_P.c,v 1.2 2007/12/14 04:37:48 detar Exp $";

/********************************************************************/
/* Load Wilson clover parameters                                    */
/********************************************************************/

static void 
load_qop_wilson_coeffs(QOP_wilson_coeffs_t *c, Real clov_c)
{
  c->clov_c       = clov_c;
}

/********************************************************************/
/* Create fermion links                                             */
/********************************************************************/

static QOP_FermionLinksWilson *
create_qop_wilson_fermion_links( Real clov )
{
  QOP_FermionLinksWilson *qop_links = NULL;
  QOP_info_t info;
  QOP_GaugeField *links;
  QOP_wilson_coeffs_t coeffs;
  double remaptime;

  /* Load coeffs structure */
  load_qop_wilson_coeffs(&coeffs, clov);

  /* Map SU(3) gauge field to G type */
  remaptime = -dclock(); 
  links = create_G_from_site4(F_OFFSET(link),EVENANDODD);
  remaptime += dclock();

  /* Create links */
  qop_links = QOP_wilson_create_L_from_G(&info, &coeffs, links);

  QOP_destroy_G(links);

#ifdef FFTIME
#ifdef REMAP
    node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
  node0_printf("FFTIME:  time = %e (cl_qop) terms = 1 mflops = %e\n",
	       info.final_sec, (Real)info.final_flop/(1e6*info.final_sec) );
#endif
  return qop_links;
}

/********************************************************************/
/* Destroy fermion links                                            */
/********************************************************************/

static void
destroy_qop_wilson_fermion_links( QOP_FermionLinksWilson *qop_links )
{
  QOP_wilson_destroy_L(qop_links);
  qop_links = NULL;
}

/********************************************************************/
/* Load inversion args for Level 3 inverter                         */
/********************************************************************/

static void 
set_qop_invert_arg( QOP_invert_arg_t* qop_invert_arg, 
		    quark_invert_control *qic, int milc_parity )
{
  qop_invert_arg->max_iter     = qic->nrestart*qic->max;
  qop_invert_arg->restart      = qic->max;
  qop_invert_arg->max_restarts = qic->nrestart;
  qop_invert_arg->evenodd      = milc2qop_parity(milc_parity);
}

/********************************************************************/
/* Load residual values for stopping                                */
/********************************************************************/

static QOP_resid_arg_t ***
create_qop_resid_arg( int nsrc, int nmass[], Real min_resid_sq )
{
  QOP_resid_arg_t ***res_arg;
  char myname[] = "create_qop_resid_arg";
  int isrc,imass;

  /* Pointers for residual errors */
  res_arg = (QOP_resid_arg_t ***)malloc(sizeof(QOP_resid_arg_t **)*nsrc);
  if(res_arg == NULL){
    printf("%s(%d): Can't allocate res_arg*\n",myname,this_node);
    terminate(1);
  }
  for(isrc = 0; isrc < nsrc; isrc++){
    res_arg[isrc] = 
      (QOP_resid_arg_t **)malloc(sizeof(QOP_resid_arg_t *)*nmass[isrc]);
    if(res_arg[isrc] == NULL){
      printf("%s(%d): Can't allocate res_arg*\n",myname,this_node);
      terminate(1);
    }
    for(imass = 0; imass < nmass[isrc]; imass++){
      res_arg[isrc][imass] = 
	(QOP_resid_arg_t *)malloc(sizeof(QOP_resid_arg_t ));
      if(res_arg[isrc][imass] == NULL){
	printf("%s(%d): Can't allocate res_arg\n",myname,this_node);
	terminate(1);
      }
      /* For now the residuals are the same for all sources and masses */
      res_arg[isrc][imass]->rsqmin = min_resid_sq;
    }
  }
  return res_arg;
}

static void
destroy_qop_resid_arg(QOP_resid_arg_t ***res_arg, int nsrc, int nmass[])
{
  int isrc, imass;

  for(isrc = 0; isrc < nsrc; isrc++){
    for(imass = 0; imass < nmass[isrc]; imass++){
      free(res_arg[isrc][imass]);
    }
    free(res_arg[isrc]);
  }
  free(res_arg);
}

/********************************************************************/
/* General MILC wrapper for Level 3 inverter                        */
/********************************************************************/

static int 
bicgilu_cl_qop_generic( QOP_FermionLinksWilson *qop_links, 
			QOP_invert_arg_t *qop_invert_arg,
			QOP_resid_arg_t  ***qop_resid_arg,
			MYREAL *kappas[], int nkappa[], 
			QOP_DiracFermion **qop_sol[], 
			QOP_DiracFermion *qop_src[], 
			int nsrc,		    
			int *final_nrestart,
			Real *final_rsq_ptr )
{
  int isrc, ikappa;
  int iters;
  QOP_info_t info;
  
  if(nsrc == 1 && nkappa[0] == 1)
    QOP_wilson_invert( &info, qop_links, qop_invert_arg, qop_resid_arg[0][0],
		       kappas[0][0], qop_sol[0][0], qop_src[0] );
  else
    QOP_wilson_invert_multi( &info, qop_links, qop_invert_arg, qop_resid_arg,
			     kappas, nkappa, qop_sol, qop_src, nsrc );

  /* For now we return the largest value and total iterations */
  *final_rsq_ptr = 0;
  *final_restart = 0;
  iters = 0;
  for(isrc = 0; isrc < nsrc; isrc++)
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      if(*final_rsq_ptr < qop_resid_arg[isrc][ikappa]->final_rsq)
	*final_rsq_ptr = qop_resid_arg[isrc][ikappa]->final_rsq;
      if(*final_restart < qop_resid_arg[isrc][ikappa]->final_restart)
	*final_restart = qop_resid_arg[isrc][ikappa]->final_restart;
      iters += qop_resid_arg[isrc][ikappa]->final_iter;
#ifdef CGDEBUG
      if(nsrc > 1 || nkappa[isrc] > 1)
	node0_printf("CONGRAD5(src %d,kappa %d): iters = %d resid = %e\n",
	       isrc, ikappa,
	       qop_resid_arg[isrc][ikappa]->final_iter,
	       qop_resid_arg[isrc][ikappa]->final_rsq);
#endif
    }

#ifdef CGTIME
  node0_printf("CGTIME: time = %e (wilson_qop %s) ",
	       info.final_sec,qop_prec[QOP_Precision-1]);
  for(isrc = 0; isrc < nsrc; isrc++)
    node0_printf("nkappa[%d] = %d iters = %d ",
		 isrc,nkappa[isrc],qop_resid_arg[isrc][0]->final_iter);
  node0_printf("mflops = %e\n", info.final_flop/(1.0e6*info.final_sec) );
  fflush(stdout);
#endif

  return iters;
}

#define MAXSRC 20

/********************************************************************/
/* Map MILC fields to QOP format and call generic QOP driver        */
/********************************************************************/

static int
bicgilu_cl_qop(quark_invert_control *qic, Real clov,
	       MYREAL *kappas[], int nkappa[], 
	       wilson_vector *milc_srcs[], 
	       wilson_vector **milc_sols[],
	       int nsrc, int *final_restart,
               Real* final_rsq_ptr, int milc_parity )
{
  int isrc, ikappa;
  QOP_FermionLinksWilson *qop_links;
  QOP_DiracFermion **qop_sol[MAXSRC], *qop_src[MAXSRC];
  int iterations_used = 0;
  QOP_invert_arg_t qop_invert_arg;
  QOP_resid_arg_t  ***qop_resid_arg;
  double remaptime;

  if(nsrc > MAXSRC){
    printf("bicgilu_cl_qop: too many sources\n");
    terminate(1);
  }

  /* Initialize QOP */
  if(initialize_qop() != QOP_SUCCESS){
    printf("bicbilu_cl_qop: Error initializing QOP\n");
    terminate(1);
  }

  /* Create QOP links object */

  qop_links = create_qop_wilson_fermion_links( clov );

  /* Set qop_invert_arg */
  set_qop_invert_arg( & qop_invert_arg, qic, milc_parity );
  
  /* Pointers for residual errors */
  qop_resid_arg = create_qop_resid_arg( nsrc, nkappa, qic->resid );

  remaptime = -dclock(); 

  /* Pointers for solution vectors */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_sol[isrc] = 
      (QOP_DiracFermion **)malloc(sizeof(QOP_DiracFermion *)*nkappa[isrc]);
    if(qop_sol[isrc] == NULL){
      printf("bicgilu_cl_qop: Can't allocate qop_sol\n");
      terminate(1);
    }
  }

  /* Map MILC source and sink to QOP fields */
  for(isrc = 0; isrc < nsrc; isrc++){
    qop_src[isrc] = create_D_from_field( milc_srcs[isrc], milc_parity);
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      qop_sol[isrc][ikappa] = 
	create_D_from_field( milc_sols[isrc][ikappa], milc_parity);
    }
  }
  
  /* Call QOP inverter */

  remaptime += dclock();
  iterations_used = bicgilu_cl_qop_generic( qop_links, &qop_invert_arg,
    qop_resid_arg, kappas, nkappa, qop_sol, qop_src, nsrc, 
    final_restart, final_rsq_ptr );
  remaptime -= dclock();
  
  /* Map qop solutions to MILC fields   */

  for(isrc = 0; isrc < nsrc; isrc++)
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++)
      unload_D_to_field( milc_sols[isrc][ikappa], 
			 qop_sol[isrc][ikappa], milc_parity );

  /* Free QOP fields  */

  destroy_qop_wilson_fermion_links( qop_links );

  for(isrc = 0; isrc < nsrc; isrc++){
    QOP_destroy_D(qop_src[isrc]);    
    qop_src[isrc] = NULL;
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      QOP_destroy_D(qop_sol[isrc][ikappa]);     
      free(qop_sol[isrc]);
      qop_sol[isrc] = NULL;
    }
  }

  remaptime += dclock();
#ifdef CGTIME
#ifdef REMAP
    node0_printf("CGREMAP:  time = %e\n",remaptime);
#endif
#endif
  destroy_qop_resid_arg(qop_resid_arg, nsrc, nkappa);
  qop_resid_arg = NULL;

  return iterations_used;
}

/********************************************************************/
/* Inverter interface for specific precision                        */
/********************************************************************/
int 
BGCGILU_MILC2QOP( wilson_vector *milc_src, wilson_vector *milc_sol, 
		  quark_invert_control *qic, void *dmp)
{
  int iterations_used;
  static MYREAL t_kappa;
  MYREAL *kappas[1];
  int nkappa[1], nsrc;
  wilson_vector *milc_srcs[1], *milc_sols0[1], **milc_sols[1];

  dirac_clover_param *dcp 
    = (dirac_clover_param *)dmp; /* Cast pass-through pointer */

  Real Kappa = dcp->Kappa;     /* hopping */
  Real Clov_c = dcp->Clov_c;   /* Perturbative clover coeff */
  Real U0 = dcp->U0;           /* Tadpole correction to Clov_c */
  Real clov = Clov_c/(U0*U0*U0); /* Full clover coefficient */
  Real final_rsq_ptr;           
  int final_restart;

  /* Set up general source and solution pointers for one mass, one source */
  nsrc = 1;
  milc_srcs[0] = milc_src;

  nkappa[0] = 1;
  t_kappa = Kappa;
  kappas[0] = &t_kappa;

  milc_sols0[0] = milc_sol;
  milc_sols[0]  = milc_sols0;

  iterations_used = 
    bicgilu_cl_qop( qic, clov, kappas, nkappa, milc_srcs,
		    milc_sols, nsrc, &final_restart, &final_rsq_ptr,
		    EVENANDODD );

  qic->size_r = final_rsq_ptr;
  qic->size_relr = 0;
  return  iterations_used;
}
