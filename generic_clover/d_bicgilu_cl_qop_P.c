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
#include "../include/generic_qopqdp.h"
#include <string.h>

/*#define CG_DEBUG*/

/* Redefinitions according to requested precision */

#if ( QOP_Precision == 1 )

#define BGCGILU_MILC2QOP   bicgilu_cl_milc2qop_F
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix
#define CREATE_RAW4_G_FROM_SITE create_raw4_F_G_from_site
#define DESTROY_RAW4_G destroy_raw4_F_G

#else

#define BGCGILU_MILC2QOP   bicgilu_cl_milc2qop_D
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix
#define CREATE_RAW4_G_FROM_SITE create_raw4_D_G_from_site
#define DESTROY_RAW4_G destroy_raw4_D_G

#endif

/* 4/29/07 C. DeTar */

/*
 * $Log: d_bicgilu_cl_qop_P.c,v $
 * Revision 1.7  2009/06/04 16:37:09  detar
 * Make clover term persistent. Accommodate changes to generic_clover/make_clov2.c
 *
 * Revision 1.6  2009/04/05 18:12:39  detar
 * Move #if 0
 *
 * Revision 1.5  2008/04/18 23:11:54  detar
 * Remove QOP_verbose setting
 *
 * Revision 1.4  2008/04/18 15:37:03  detar
 * Support qopqdp-0.11.2
 *
 * Revision 1.3  2008/03/28 15:48:09  detar
 * Report absolute and relative residuals
 *
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

static char* cvsHeader = "$Header: /lqcdproj/detar/cvsroot/milc_qcd/generic_clover/d_bicgilu_cl_qop_P.c,v 1.7 2009/06/04 16:37:09 detar Exp $";

#if 0

/********************************************************************/
/* Load Wilson clover parameters                                    */
/********************************************************************/

static void 
load_qop_wilson_coeffs(QOP_wilson_coeffs_t *c, Real clov_c)
{
  c->clov_c       = clov_c;
  c->aniso        = 0.;
}


/* The QOP/QDP routine is not yet available */

/********************************************************************/
/* Create fermion links using QOP routine                           */
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

#else

#define DI(x) (-(1 - (x))/2.)
#define TR(x) ((x)/2.)

static void map_milc_clov_to_qop_raw(MYREAL *raw_clov, clover *milc_clov ){
  int i,c;
  site *s;
  MYREAL *r;

  r = raw_clov;
  FORALLSITES(i,s){

    for(c = 0; c < 3; c++){
      r[2*c]   = DI(milc_clov->clov_diag[i].di[0][c]);      /* c0 c0 */
      r[2*c+1] = DI(milc_clov->clov_diag[i].di[0][c+3]);    /* c1 c1 */
    }
      
    r += 6;

    r[ 0] = TR( milc_clov->clov[i].tr[0][ 3].real);         /* 01 00 */
    r[ 1] = TR( milc_clov->clov[i].tr[0][ 3].imag);         /* 01 00 */
    r[ 2] = TR( milc_clov->clov[i].tr[0][ 0].real);         /* 10 00 */
    r[ 3] = TR( milc_clov->clov[i].tr[0][ 0].imag);         /* 10 00 */
    r[ 4] = TR( milc_clov->clov[i].tr[0][ 6].real);         /* 11 00 */
    r[ 5] = TR( milc_clov->clov[i].tr[0][ 6].imag);         /* 11 00 */
    r[ 6] = TR( milc_clov->clov[i].tr[0][ 1].real);         /* 20 00 */
    r[ 7] = TR( milc_clov->clov[i].tr[0][ 1].imag);         /* 20 00 */
    r[ 8] = TR( milc_clov->clov[i].tr[0][10].real);         /* 21 00 */
    r[ 9] = TR( milc_clov->clov[i].tr[0][10].imag);         /* 21 00 */

    r[10] = TR( milc_clov->clov[i].tr[0][ 4].real);         /* 10 01 */
    r[11] = TR(-milc_clov->clov[i].tr[0][ 4].imag);         /* 10 01 */
    r[12] = TR( milc_clov->clov[i].tr[0][ 9].real);         /* 11 01 */
    r[13] = TR( milc_clov->clov[i].tr[0][ 9].imag);         /* 11 01 */
    r[14] = TR( milc_clov->clov[i].tr[0][ 5].real);         /* 20 01 */
    r[15] = TR(-milc_clov->clov[i].tr[0][ 5].imag);         /* 20 01 */
    r[16] = TR( milc_clov->clov[i].tr[0][13].real);         /* 21 01 */
    r[17] = TR( milc_clov->clov[i].tr[0][13].imag);         /* 21 01 */

    r[18] = TR( milc_clov->clov[i].tr[0][ 7].real);         /* 11 10 */
    r[19] = TR( milc_clov->clov[i].tr[0][ 7].imag);         /* 11 10 */
    r[20] = TR( milc_clov->clov[i].tr[0][ 2].real);         /* 20 10 */
    r[21] = TR( milc_clov->clov[i].tr[0][ 2].imag);         /* 20 10 */
    r[22] = TR( milc_clov->clov[i].tr[0][11].real);         /* 21 10 */
    r[23] = TR( milc_clov->clov[i].tr[0][11].imag);         /* 21 10 */

    r[24] = TR( milc_clov->clov[i].tr[0][ 8].real);         /* 20 11 */
    r[25] = TR(-milc_clov->clov[i].tr[0][ 8].imag);         /* 20 11 */
    r[26] = TR( milc_clov->clov[i].tr[0][14].real);         /* 21 11 */
    r[27] = TR( milc_clov->clov[i].tr[0][14].imag);         /* 21 11 */

    r[28] = TR( milc_clov->clov[i].tr[0][12].real);         /* 21 20 */
    r[29] = TR( milc_clov->clov[i].tr[0][12].imag);         /* 21 20 */

    r += 30;

    for(c = 0; c < 3; c++){
      r[2*c]   = DI(milc_clov->clov_diag[i].di[1][c]);      /* c2 c2 */
      r[2*c+1] = DI(milc_clov->clov_diag[i].di[1][c+3]);    /* c3 c3 */
    }

    r += 6;

    r[ 0] = TR( milc_clov->clov[i].tr[1][ 3].real);         /* 03 02 */
    r[ 1] = TR( milc_clov->clov[i].tr[1][ 3].imag);         /* 03 02 */
    r[ 2] = TR( milc_clov->clov[i].tr[1][ 0].real);         /* 12 02 */
    r[ 3] = TR( milc_clov->clov[i].tr[1][ 0].imag);         /* 12 02 */
    r[ 4] = TR( milc_clov->clov[i].tr[1][ 6].real);         /* 13 02 */
    r[ 5] = TR( milc_clov->clov[i].tr[1][ 6].imag);         /* 13 02 */
    r[ 6] = TR( milc_clov->clov[i].tr[1][ 1].real);         /* 22 02 */
    r[ 7] = TR( milc_clov->clov[i].tr[1][ 1].imag);         /* 22 02 */
    r[ 8] = TR( milc_clov->clov[i].tr[1][10].real);         /* 23 02 */
    r[ 9] = TR( milc_clov->clov[i].tr[1][10].imag);         /* 23 02 */

    r[10] = TR( milc_clov->clov[i].tr[1][ 4].real);         /* 12 03 */
    r[11] = TR(-milc_clov->clov[i].tr[1][ 4].imag);         /* 12 03 */
    r[12] = TR( milc_clov->clov[i].tr[1][ 9].real);         /* 13 03 */
    r[13] = TR( milc_clov->clov[i].tr[1][ 9].imag);         /* 13 03 */
    r[14] = TR( milc_clov->clov[i].tr[1][ 5].real);         /* 22 03 */
    r[15] = TR(-milc_clov->clov[i].tr[1][ 5].imag);         /* 22 03 */
    r[16] = TR( milc_clov->clov[i].tr[1][13].real);         /* 23 03 */
    r[17] = TR( milc_clov->clov[i].tr[1][13].imag);         /* 23 03 */

    r[18] = TR( milc_clov->clov[i].tr[1][ 7].real);         /* 13 12 */
    r[19] = TR( milc_clov->clov[i].tr[1][ 7].imag);         /* 13 12 */
    r[20] = TR( milc_clov->clov[i].tr[1][ 2].real);         /* 22 12 */
    r[21] = TR( milc_clov->clov[i].tr[1][ 2].imag);         /* 22 12 */
    r[22] = TR( milc_clov->clov[i].tr[1][11].real);         /* 23 12 */
    r[23] = TR( milc_clov->clov[i].tr[1][11].imag);         /* 23 12 */

    r[24] = TR( milc_clov->clov[i].tr[1][ 8].real);         /* 22 13 */
    r[25] = TR(-milc_clov->clov[i].tr[1][ 8].imag);         /* 22 13 */
    r[26] = TR( milc_clov->clov[i].tr[1][14].real);         /* 23 13 */
    r[27] = TR( milc_clov->clov[i].tr[1][14].imag);         /* 23 13 */

    r[28] = TR( milc_clov->clov[i].tr[1][12].real);         /* 23 22 */
    r[29] = TR( milc_clov->clov[i].tr[1][12].imag);         /* 23 22 */

    r += 30;
  }
}

/* Use MILC code to construct clover term for now */

/* Extract pure clover from MILC */

static QOP_FermionLinksWilson *
create_qop_wilson_fermion_links( Real clov ){

  clover *milc_clov = gen_clov;
  MYSU3_MATRIX **raw_links;
  MYREAL *raw_clov;
  QOP_FermionLinksWilson *qop_links;

  /* Construct raw QOP clover term from MILC clover term */
  if(clov == 0){
    raw_clov = NULL;
  }
  else{
    raw_clov = (MYREAL *)malloc(72*sites_on_node*sizeof(MYREAL));
    if(raw_clov == NULL){
      printf("create_qop_wilson_fermion_links(%d): no room for raw_clov\n",
	     this_node);
      terminate(1);
    }
    // milc_clov = create_clov(); /* Note Real clov has no kappa factor! */
    if(milc_clov == NULL) terminate(1);
    compute_clov(milc_clov,clov);
    map_milc_clov_to_qop_raw(raw_clov, milc_clov);
    //free_this_clov(milc_clov);
  }

  raw_links = CREATE_RAW4_G_FROM_SITE(F_OFFSET(link), EVENANDODD);
  if(raw_links == NULL)terminate(1);

  /* Map QOP/QDP raw to QOP/QDP structure */

  qop_links = QOP_wilson_create_L_from_raw((MYREAL **)raw_links, raw_clov, 
					   QOP_EVENODD);
  
  DESTROY_RAW4_G(raw_links); raw_links = NULL;
  free(raw_clov);

  return qop_links;
}

#endif

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
			int *final_restart,
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
#ifdef CG_DEBUG
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

/* temporary hack until we get a more flexible QOP inverter that can
   handle our gamma matrix conventions */

static void
gamma5_flip(wilson_vector *milc, int parity){
  int i,c;
  site *s;
  FORSOMEPARITY(i,s,parity){
    for(c = 0; c < 3; c++){
      CNEGATE(milc[i].d[2].c[c],milc[i].d[2].c[c]);
      CNEGATE(milc[i].d[3].c[c],milc[i].d[3].c[c]);
    }
  }
}

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
  int i;
  site *s;

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
  qop_resid_arg = create_qop_resid_arg( nsrc, nkappa, (qic->resid)*(qic->resid));

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
    gamma5_flip(milc_srcs[isrc], milc_parity);  /* compensate for QOP gamma */
    qop_src[isrc] = create_D_from_field( milc_srcs[isrc], milc_parity);
    gamma5_flip(milc_srcs[isrc], milc_parity);  /* restore the source */
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      /* Adjust normalization for MILC conventions */
      gamma5_flip(milc_sols[isrc][ikappa], milc_parity);  /* compensate for QOP gamma */
      FORALLSITES(i,s){
	scalar_mult_wvec( milc_sols[isrc][ikappa]+i, 2.*kappas[isrc][ikappa],
			  milc_sols[isrc][ikappa]+i);
      }
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
    for(ikappa = 0; ikappa < nkappa[isrc]; ikappa++){
      unload_D_to_field( milc_sols[isrc][ikappa], 
			 qop_sol[isrc][ikappa], milc_parity );
      /* Adjust normalization for MILC conventions */
      gamma5_flip(milc_sols[isrc][ikappa], milc_parity);  /* compensate for QOP gamma */
      FORALLSITES(i,s){
	scalar_mult_wvec( milc_sols[isrc][ikappa]+i, 1/(2.*kappas[isrc][ikappa]),
			  milc_sols[isrc][ikappa]+i);
      }
    }

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
  Real final_rsq_val;           
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
		    milc_sols, nsrc, &final_restart, &final_rsq_val,
		    EVENANDODD );

  qic->size_r = 0;  /* We don't see the cumulative resid with QOP*/
  qic->size_relr = 0;
  qic->final_rsq = final_rsq_val;
  qic->final_relrsq = 0;
  qic->final_iters = iterations_used;
  qic->final_restart = final_restart;
  return  iterations_used;
}
