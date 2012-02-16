/****** fermion_force_asqtad_qop_milc_P.c  -- ******************/
/* MIMD version 7 */
/* fermion force optimized for the Asqtad action

 * This version implements the QOP API with the standard MILC algorithm
 * It is intended for testing other QOP routines.
 * It is based on fermion_force_asqtad3.c

 * Uses restart-gathers and a bit more memory for better performance
 * The algorithm is the same as fermion_force_asqtad2.c, except that
 * the shifts are consistent with Dslash.
 * D.T. 1/28/98, starting from gauge_stuff.c
 * K.O. 3/99 Added optimized fattening for Asq actions
 * D.T. 4/99 Combine force calculations for both mass quarks
 * K.O. 4/99 Optimized force for Asq action
 * S.G. 7/01, modified to use t_longlink and t_fatlink
 * C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
 *
 * J.O. 3/04 Rearranged loops for optimization
 * J.O. C.D. 3/04 Copied forward links for optimization and 
 *                kept mtags open for restart_gather
 *                Worked with pointers where possible to avoid copying.
 * C.D. 6/04 Corrected memory leak
 * C.D. 3/05 Separated from quark_stuff4.c
 * C.D. 11/05 Converted to the QOP API 
 *
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
 */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED - C. DeTar
 * Fermion force: 253935 for QOP_asqtad_force()
 * Fermion force: 433968 for QOP_asqtad_force_two()
 */

#if ( QOP_Precision == 1 )
#define MYREAL float
#define MYSU3_MATRIX fsu3_matrix
#else
#define MYREAL double
#define MYSU3_MATRIX dsu3_matrix
#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/qop_milc.h"
#include <string.h>

/* This routine is valid only for Asqtad, so requires the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

static void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) ;
static void add_force_to_mom(su3_vector *back, su3_vector *forw, int dir, Real coef);
static void side_link_force(int mu, int nu, Real coeff, su3_vector *Path,
		     su3_vector *Path_nu, su3_vector *Path_mu, 
		     su3_vector *Path_numu) ;

static void u_shift_hw_fermion(half_wilson_vector *src, half_wilson_vector *dest, 
			int dir, msg_tag** mtag, half_wilson_vector *tmpvec ) ;
static void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, int dir, Real coeff[2]) ;
static void side_link_3f_force(int mu, int nu, Real coeff[2], 
			half_wilson_vector *Path   , 
			half_wilson_vector *Path_nu, 
			half_wilson_vector *Path_mu, 
			half_wilson_vector *Path_numu) ;

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

static su3_matrix *forwardlink[4];
static su3_matrix *tmpmom[4];
static su3_matrix *backwardlink[4];
static su3_vector *temp_x;

/*---------------------------------------------------------------------------*/
/* Convert (or copy) four single precision su3_matrices to generic precision */

#if (PRECISION==QOP_Precision)


static void 
f2p_mat(su3_matrix *dest, fsu3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
p2f_mat(fsu3_matrix *dest, su3_matrix *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_matrix));
}

static void 
d2p_mat(su3_matrix *dest, dsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
p2d_mat(dsu3_matrix *dest, su3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

static void 
f2p_vec(su3_vector *dest, fsu3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
p2f_vec(fsu3_vector *dest, su3_vector *src){
  memcpy((void *)dest, (void *)src, sizeof(fsu3_vector));
}

static void 
d2p_vec(su3_vector *dest, dsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

static void 
p2d_vec(dsu3_vector *dest, su3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Just use pointer */

static su3_matrix *get_matrix_from_qop(su3_matrix *g){
  return g;
}

static su3_vector *get_vector_from_qop(su3_vector *v){
  return v;
}

/* Do nothing */

static void copy_qop_from_matrix(su3_matrix *g, su3_matrix *mat){
}

static void destroy_matrix(su3_matrix *g){
}

static void destroy_vector(su3_vector *v){
}

#else

static void 
f2d_mat(dsu3_matrix *dest, fsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}

/* Convert (or copy) su3_matrix from prevailing to single precision */
static void 
d2f_mat(fsu3_matrix *dest, dsu3_matrix *src){
  int i,j;
  
  for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
    dest->e[i][j].real = src->e[i][j].real;
    dest->e[i][j].imag = src->e[i][j].imag;
  }
}


static void 
f2d_vec(dsu3_vector *dest, fsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

/* Convert (or copy) su3_vector from prevailing to single precision */
static void 
d2f_vec(fsu3_vector *dest, dsu3_vector *src){
  int i;
  
  for(i = 0; i < 3; i++){
    dest->c[i].real = src->c[i].real;
    dest->c[i].imag = src->c[i].imag;
  }
}

#if (QOP_Precision == 1)

static dsu3_matrix *get_matrix_from_qop(fsu3_matrix *g){
  int i;
  site *s;
  dsu3_matrix *mat;
 
  mat = (dsu3_matrix *)malloc(sizeof(dsu3_matrix)*sites_on_node);
  if(mat == NULL)terminate(1);

  FORALLSITES(i,s){
    f2d_mat(mat+i, g+i);
  }
}

static dsu3_vector *get_vector_from_qop(fsu3_vector *v){
  int i;
  site *s;
  dsu3_vector *vec;
 
  vec = (dsu3_vector *)malloc(sizeof(dsu3_vector)*sites_on_node);
  if(vec == NULL)terminate(1);

  FORALLSITES(i,s){
    f2d_vec(vec+i, v+i);
  }
}

static void copy_qop_from_matrix(fsu3_matrix *g, dsu3_matrix *mat){
  int i;
  site *s;

  FORALLSITES(i,s){
    d2f_mat(g+i, mat+i);
  }
}

static void destroy_matrix(dsu3_matrix *mat){
  free(mat);
}

static void destroy_vector(dsu3_vector *vec){
  free(vec);
}

#else

static fsu3_matrix *get_matrix_from_qop(dsu3_matrix *g){
  int i;
  site *s;
  fsu3_matrix *mat;
 
  mat = (fsu3_matrix *)malloc(sizeof(fsu3_matrix)*sites_on_node);
  if(mat == NULL)terminate(1);

  FORALLSITES(i,s){
    d2f_mat(mat+i, g+i);
  }
  return mat;
}

static fsu3_vector *get_vector_from_qop(dsu3_vector *v){
  int i;
  site *s;
  fsu3_vector *vec;
 
  vec = (fsu3_vector *)malloc(sizeof(fsu3_vector)*sites_on_node);
  if(vec == NULL)terminate(1);

  FORALLSITES(i,s){
    d2f_vec(vec+i, v+i);
  }
  return vec;
}

static void copy_qop_from_matrix(dsu3_matrix *g, fsu3_matrix *mat){
  int i;
  site *s;

  FORALLSITES(i,s){
    f2d_mat(g+i, mat+i);
  }
}

static void destroy_matrix(fsu3_matrix *mat){
  free(mat);
}

static void destroy_vector(fsu3_vector *vec){
  free(vec);
}
#endif

#endif

/**********************************************************************/
/*   Version for a single set of degenerate flavors                   */
/**********************************************************************/

/* Optimized force code for the Asq and Asqtad actions                 *
 * I assume that path 0 is the one link path 2 the 3-staple            *
 * path 3 the 5-staple path 4 the 7-staple and path 5 the Lapage term. *
 * Path 1 is the Naik term.                                            */
#define Pmu          tempvec[0] 
#define Pnumu        tempvec[1]
#define Prhonumu     tempvec[2]
#define P7           tempvec[3]
#define P7rho        tempvec[4]              
#define P7rhonu      tempvec[5]
#define P5           tempvec[6]
#define P3           tempvec[7]
#define P5nu         tempvec[3]
#define P3mu         tempvec[3]
#define Popmu        tempvec[4]
#define Pmumumu      tempvec[4]
void QOP_asqtad_force(QOP_info_t *info,
		      QOP_GaugeField *gauge,
		      QOP_Force *force,
		      QOP_asqtad_coeffs_t *coeffs, MYREAL eps,
		      QOP_ColorVector *in_pt){
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  Real coeff;
  Real OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt ;
  su3_vector *tempvec[8] ;
  int dir;

  /* Timing */

  Real final_flop;
  Real nflop = 253935;
  double dtime = -dclock();

  info->status = QOP_FAIL; 
  
  /* Parity requirements */
  if(gauge->evenodd != QOP_EVENODD ||
     force->evenodd != QOP_EVENODD ||
     in_pt->evenodd != QOP_EVENODD
     )
    {
      printf("QOP_asqtad_force: Bad parity src %d gauge %d force %d\n",
	     in_pt->evenodd, gauge->evenodd, force->evenodd);
      return;
    }

  /* Map field pointers to local static pointers */

  FORALLUPDIR(dir){
    forwardlink[dir] = get_matrix_from_qop(gauge->g + dir*sites_on_node);
    tmpmom[dir]      = get_matrix_from_qop(force->f + dir*sites_on_node);
  }
  temp_x = get_vector_from_qop(in_pt->v);

  /* Load path coefficients from table */

  /* Path coefficients times fermion epsilon */
  OneLink = coeffs->one_link*eps ; 
  Naik    = coeffs->naik*eps ;
  ThreeSt = coeffs->three_staple*eps ;
  FiveSt  = coeffs->five_staple*eps ;
  SevenSt = coeffs->seven_staple*eps ;
  Lepage  = coeffs->lepage*eps ;       
  /* *************************************** */

  /* Initialize the DirectLink flags */
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

  /* Allocate temporary vectors */
  for(mu=0;mu<8;mu++)
    tempvec[mu] = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

  for(sig=0;sig<8;sig++)
    {
      for(mu=0;mu<8;mu++)if((mu!=sig)&&(mu!=OPP_DIR(sig)))
	{
	  u_shift_fermion(temp_x, Pmu, OPP_DIR(mu));
	  u_shift_fermion(Pmu, P3, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
	      add_force_to_mom(P3, Pmu, sig, -ThreeSt) ;
	    }
	  for(nu=0;nu<8;nu++)if((nu!=mu )&&(nu!=OPP_DIR(mu ))&&
				(nu!=sig)&&(nu!=OPP_DIR(sig)))
	    {
	      u_shift_fermion(Pmu, Pnumu, OPP_DIR(nu));
	      u_shift_fermion(Pnumu, P5, sig);
	      if(GOES_FORWARDS(sig))
		{
		  /* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
		  add_force_to_mom(P5, Pnumu, sig, FiveSt);
		}
	      for(rho=0;rho<8;rho++)if((rho!=mu )&&(rho!=OPP_DIR(mu ))&&
				       (rho!=nu )&&(rho!=OPP_DIR(nu ))&&
				       (rho!=sig)&&(rho!=OPP_DIR(sig)))
		{
		  u_shift_fermion(Pnumu, Prhonumu, OPP_DIR(rho));
		  /* Length 7 paths */
		  u_shift_fermion(Prhonumu, P7,sig);
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      add_force_to_mom(P7, Prhonumu, sig, -SevenSt ) ;
		    }
		  /*Add the force F_rho the 2(4) link in the path: +     */
		  u_shift_fermion(P7, P7rho, rho);
		  side_link_force(rho,sig,SevenSt, Pnumu, P7, Prhonumu, P7rho);
		  /* Add the P7rho vector to P5 */
		  if(FiveSt != 0)coeff = SevenSt/FiveSt ; else coeff = 0;
		  FORALLSITES(i,s)
		    scalar_mult_add_su3_vector(&(P5[i]),&(P7rho[i]),coeff,
					       &(P5[i]));
		}/* rho */
	      /* Length 5 paths */
	      /*Add the force F_nu the 1(3) link in the path: -     */
	      u_shift_fermion(P5,P5nu, nu);
	      side_link_force(nu,sig,-FiveSt,Pmu,P5, 
			      Pnumu,P5nu) ;
	      /* Add the P5nu vector to P3 */
	      if(ThreeSt != 0)coeff = FiveSt/ThreeSt ; else coeff = 0;
	      FORALLSITES(i,s)
		scalar_mult_add_su3_vector(&(P3[i]),&(P5nu[i]),coeff,&(P3[i]));
	    }/* nu */

	  /* Now the Lepage term... It is the same with 5-link paths with
             nu=mu and FiveSt=Lepage. So Pnumu is really Pmumu */
	  u_shift_fermion(Pmu, Pnumu, OPP_DIR(mu));
	  u_shift_fermion(Pnumu, P5, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      add_force_to_mom(P5, Pnumu, sig, Lepage) ;
	    }
	  /*Add the force F_nu the 1(3) link in the path: -     */
	  u_shift_fermion(P5,P5nu, mu);
	  side_link_force(mu, sig, -Lepage, Pmu, P5, Pnumu, P5nu) ;
	  /* Add the P5nu vector to P3 */
	  if(ThreeSt != 0) coeff = Lepage/ThreeSt ; else coeff = 0;
	  FORALLSITES(i,s)
	    scalar_mult_add_su3_vector(&(P3[i]),&(P5nu[i]),coeff,&(P3[i]));

	  /* Length 3 paths (Not the Naik term) */
	  /*Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    u_shift_fermion(P3,P3mu, mu );
	  /* The above shift is not needed if mu is backwards */
	  side_link_force(mu, sig, ThreeSt, temp_x, P3, Pmu, P3mu);

	  /* Finally the OneLink and the Naik term */
	  /* Check if this direction is not already done */
	  if( (!DirectLinks[mu]) ){
	    DirectLinks[mu]=1 ;
	    if(GOES_BACKWARDS(mu))/* Do only the forward terms in the Dslash */
	      {
		/* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
		 * shift.                                                   */
		/* The one link */
		add_force_to_mom(Pmu, temp_x, OPP_DIR(mu), OneLink) ;
		/* For the same reason Pnumu is the forward double link */

		/* Popmu is a backward shift */
		u_shift_fermion(temp_x, Popmu, mu);
		/* The Naik */
		/* link no 1: - */
		add_force_to_mom(Pnumu, Popmu, OPP_DIR(mu), -Naik) ;
		/*Pmumumu can overwrite Popmu which is no longer needed */
		u_shift_fermion(Pnumu, Pmumumu, OPP_DIR(mu));
		/* link no 0: + */
		add_force_to_mom(Pmumumu, temp_x, OPP_DIR(mu), Naik);
	      }
	    else /* The rest of the Naik terms */
	      {
		u_shift_fermion(temp_x, Popmu, mu);
		/* link no 2: + */
		/* Pnumu is double backward shift */
		add_force_to_mom(Popmu, Pnumu, mu, Naik) ;
	      }
	  }
	}/* mu */
      /* Here we have to do together the Naik term and the one link term */
    }/*sig */

  /* Copy force back and clean up */

  FORALLUPDIR(dir){
    copy_qop_from_matrix(force->f + dir*sites_on_node, tmpmom[dir]);
    destroy_matrix(tmpmom[dir]);
    destroy_matrix(forwardlink[dir]);
  }
  destroy_vector(temp_x);

  /* Free temporary vectors */
  for(mu=0;mu<8;mu++)
    free(tempvec[mu]) ;

  final_flop = (Real)nflop*volume/numnodes();
  dtime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (qop_milc) terms = 1 mflops = %e\n",dtime,
	       final_flop/(1.0e6*dtime) );
#endif

 info->final_sec = dtime; 
 info->final_flop = final_flop;
 info->status = QOP_SUCCESS; 
} /* QOP_asqtad_force(version 7) */
#undef Pmu          
#undef Pnumu        
#undef Prhonumu     
#undef P7           
#undef P7rho        
#undef P7rhonu      
#undef P5           
#undef P3           
#undef P5nu         
#undef P3mu         
#undef Popmu        
#undef Pmumumu      

/**********************************************************************/
/*   Version for multiple sources - Asqtad only                       */
/**********************************************************************/

static su3_matrix *backwardlink[4];
static su3_matrix *tmpmom[4];
static su3_vector *temp_x1;
static su3_vector *temp_x2;

void QOP_asqtad_force_multi(QOP_info_t *info,
			    QOP_GaugeField *gauge,
			    QOP_Force *force,
			    QOP_asqtad_coeffs_t *coef, MYREAL eps[],
			    QOP_ColorVector *in_pt[], int nsrc){

  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* 4/15/99 combine force from two different mass quarks, (eg 2+1flavors) */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  int i;
  site *s;
  int mu,nu,rho,sig;
  int dir;
  Real coeff[2];
  Real OneLink[2], Lepage[2], Naik[2], FiveSt[2], ThreeSt[2], SevenSt[2];
  Real mNaik[2], mLepage[2], mFiveSt[2], mThreeSt[2], mSevenSt[2];
  half_wilson_vector *Pnumu, *Prhonumu, *P7, *P7rho, *P5nu, 
    *P3mu=NULL, *P5sig, *Popmu, *Pmumumu;
  half_wilson_vector *P3[8];
  half_wilson_vector *P5[8];
  half_wilson_vector *temp_x;
  half_wilson_vector *Pmu;
  half_wilson_vector *Pmumu;
  half_wilson_vector *temp_hw[8];
  msg_tag *mt[8];
  msg_tag *mtag[4];
  half_wilson_vector *hw[8];
  int isrc;

  /* Timing */

  Real final_flop;
  double dtime = -dclock();
  Real nflop = 433968;
  QOP_info_t info_one = {0., 0., 0, 0, 0};

  /* Parity requirements */
  for(isrc = 0; isrc < nsrc; isrc++){
    if(in_pt[isrc]->evenodd != QOP_EVENODD){
      printf("QOP_asqtad_force_multi: Bad parity in_pt[%d] %d\n",
	     isrc, in_pt[isrc]->evenodd);
      return;
    }
  }
  if(gauge->evenodd != QOP_EVENODD ||
     force->evenodd != QOP_EVENODD 
     )
    {
      printf("QOP_asqtad_force_multi: Bad parity gauge %d force %d\n",
	     gauge->evenodd, force->evenodd);
      return;
    }

  /* If not two sources, just do single source repeatedly */

  if(nsrc != 2){
    /* Initialize status */
    info->final_sec = 0;
    info->final_flop = 0;

    for(isrc = 0; isrc < nsrc; isrc++){
      QOP_asqtad_force(&info_one, gauge, force, coef, eps[isrc], in_pt[isrc]);

      /* Accumulate status */
      info->final_sec += info_one.final_sec;
      info->final_flop += info_one.final_flop;
      info->status = info_one.status;
      if(info->status != QOP_SUCCESS)break;
    }
    return;
  }

  /* Two-source version from here on */

  FORALLUPDIR(dir){
    forwardlink[dir] = get_matrix_from_qop(gauge->g + dir*sites_on_node);
    tmpmom[dir]      = get_matrix_from_qop(force->f + dir*sites_on_node);
  }
  temp_x1 = get_vector_from_qop(in_pt[0]->v);
  temp_x2 = get_vector_from_qop(in_pt[1]->v);

  /* Allocate temporary hw vector */
  for(mu = 0; mu < 8; mu++){
    half_wilson_vector *pt;
    pt = (half_wilson_vector *)
      special_alloc(sites_on_node*sizeof(half_wilson_vector));
    if(pt == NULL){
      printf("QOP_asqtad_force_multi: No room for hw\n");
      return;
    }
    hw[mu] = pt;
  }

  Pmu = 
    (half_wilson_vector *)special_alloc(sites_on_node*sizeof(half_wilson_vector));
  if(Pmu == NULL){
    printf("QOP_asqtad_force_multi: No room for Pmu\n");
    return;
  }
  
  Pmumu = 
    (half_wilson_vector *)special_alloc(sites_on_node*sizeof(half_wilson_vector));
  if(Pmumu == NULL){
    printf("QOP_asqtad_force_multi: No room for Pmumu\n");
    return;
  }
  
  /* Allocate temporary vectors */
  for(mu=0; mu<8; mu++) {
    P3[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(P3[mu] == NULL){
      printf("QOP_asqtad_fermion_force: No room for P3\n");
      return;
    }
    P5[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(P5[mu] == NULL){
      printf("QOP_asqtad_fermion_force: No room for P5\n");
      return;
    }
  }

  /* Initialize message pointers */
  for(mu = 0; mu < 8; mu++)mt[mu] = NULL;

  /* Double store backward gauge links */
  FORALLUPDIR(dir){
    su3_matrix *pt;
    pt = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(pt == NULL){
      printf("QOP_asqtad_fermion_force: No room for backwardlink\n");
      return;
    }
    backwardlink[dir] = pt;
  }
  
  /* Gather backward links */
  FORALLUPDIR(dir){
    mtag[dir] = start_gather_field( forwardlink[dir], sizeof(su3_matrix), 
			      OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
  }
  FORALLUPDIR(dir){
    wait_gather(mtag[dir]);
    FORALLSITES(i,s){
      backwardlink[dir][i] = *((su3_matrix *)gen_pt[dir][i]);
    }
    cleanup_gather(mtag[dir]);
  }
  
  /* Path coefficients times fermion epsilon */

  OneLink[0] = coef->one_link*eps[0];
  Naik[0]    = coef->naik*eps[0];         mNaik[0]    = -Naik[0];
  ThreeSt[0] = coef->three_staple*eps[0]; mThreeSt[0] = -ThreeSt[0];
  FiveSt[0]  = coef->five_staple*eps[0];  mFiveSt[0]  = -FiveSt[0];
  SevenSt[0] = coef->seven_staple*eps[0]; mSevenSt[0] = -SevenSt[0];
  Lepage[0]  = coef->lepage*eps[0];       mLepage[0]  = -Lepage[0];

  OneLink[1] = coef->one_link*eps[1];
  Naik[1]    = coef->naik*eps[1];         mNaik[1]    = -Naik[1];
  ThreeSt[1] = coef->three_staple*eps[1]; mThreeSt[1] = -ThreeSt[1];
  FiveSt[1]  = coef->five_staple*eps[1];  mFiveSt[1]  = -FiveSt[1];
  SevenSt[1] = coef->seven_staple*eps[1]; mSevenSt[1] = -SevenSt[1];
  Lepage[1]  = coef->lepage*eps[1];       mLepage[1]  = -Lepage[1];

  /* *************************************** */

  /* Allocate temporary vectors */

  for(mu = 0; mu < 8; mu++){
    temp_hw[mu] = 
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(temp_hw[mu] == NULL){
      printf("QOP_asqtad_fermion_force: No room for temp_hw\n");
      return;
    }
  }

  /* copy x_off to a temporary vector */
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  FORALLSITES(i,s)
    {
      temp_x[i].h[0] = temp_x1[i];
      temp_x[i].h[1] = temp_x2[i];
    }

  for(mu=0; mu<8; mu++)
    {
      u_shift_hw_fermion(temp_x, Pmu, OPP_DIR(mu), 
			 &mt[OPP_DIR(mu)], temp_hw[OPP_DIR(mu)]);

      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  u_shift_hw_fermion(Pmu, P3[sig], sig, &mt[sig], temp_hw[sig]);

	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
	      add_3f_force_to_mom(P3[sig], Pmu, sig, mThreeSt);
	    }
	}
      for(nu=0; nu<8; nu++) if( (nu!=mu)&&(nu!=OPP_DIR(mu)) )
	{
	  Pnumu = hw[OPP_DIR(nu)];
	  u_shift_hw_fermion(Pmu, Pnumu, OPP_DIR(nu), &mt[OPP_DIR(nu)], 
			     temp_hw[OPP_DIR(nu)]);
	  for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
				       (sig!=nu)&&(sig!=OPP_DIR(nu)) )
	    {
	      u_shift_hw_fermion(Pnumu, P5[sig], sig, &mt[sig], 
				 temp_hw[sig]);

	      if(GOES_FORWARDS(sig))
		{
		  /* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
		  add_3f_force_to_mom(P5[sig], Pnumu, sig, FiveSt);
		}
	    }
	  for(rho=0; rho<8; rho++) if( (rho!=mu)&&(rho!=OPP_DIR(mu)) &&
				       (rho!=nu)&&(rho!=OPP_DIR(nu)) )
	    {
	      Prhonumu = hw[OPP_DIR(rho)];
	      u_shift_hw_fermion(Pnumu, Prhonumu, OPP_DIR(rho), 
				 &mt[OPP_DIR(rho)], temp_hw[OPP_DIR(rho)] );
	      for(sig=0; sig<8; sig++) if( (sig!=mu )&&(sig!=OPP_DIR(mu )) &&
					   (sig!=nu )&&(sig!=OPP_DIR(nu )) &&
					   (sig!=rho)&&(sig!=OPP_DIR(rho)) )
		{
		  /* Length 7 paths */
		  P7 = hw[sig];
		  u_shift_hw_fermion(Prhonumu, P7, sig, &mt[sig], 
				     temp_hw[sig] );
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      add_3f_force_to_mom(P7, Prhonumu, sig, mSevenSt);
		    }
		  /* Add the force F_rho the 2(4) link in the path: +     */
		  P7rho = hw[rho];
		  u_shift_hw_fermion(P7, P7rho, rho, &mt[rho], 
				     temp_hw[rho]);
		  side_link_3f_force(rho,sig,SevenSt,Pnumu,P7,Prhonumu,P7rho);
		  /* Add the P7rho vector to P5 */
		  if(FiveSt[0] != 0)coeff[0] = SevenSt[0]/FiveSt[0];
		  else coeff[0] = 0;
		  if(FiveSt[1] != 0)coeff[1] = SevenSt[1]/FiveSt[1];
		  else coeff[1] = 0;
#ifdef FFSTIME
	time = -dclock();
#endif
		  FORALLSITES(i,s)
		    {
		      scalar_mult_add_su3_vector(&(P5[sig][i].h[0]),
						 &(P7rho[i].h[0]),coeff[0],
						 &(P5[sig][i].h[0]));
		      scalar_mult_add_su3_vector(&(P5[sig][i].h[1]),
						 &(P7rho[i].h[1]),coeff[1],
						 &(P5[sig][i].h[1]));
		    }
#ifdef FFSTIME
      time += dclock();
      node0_printf("FFSHIFT time4 = %e\n",time);
#endif
		} /* sig */
	    } /* rho */
	  for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
				       (sig!=nu)&&(sig!=OPP_DIR(nu)) )
	    {
	      /* Length 5 paths */
	      /* Add the force F_nu the 1(3) link in the path: -     */
	      P5nu = hw[nu];
	      u_shift_hw_fermion(P5[sig], P5nu, nu, &mt[nu], temp_hw[nu]);
	      side_link_3f_force(nu, sig, mFiveSt, Pmu, P5[sig], Pnumu, P5nu);
	      /* Add the P5nu vector to P3 */
	      if(ThreeSt[0] != 0)coeff[0] = FiveSt[0]/ThreeSt[0]; 
	      else coeff[0] = 0;
	      if(ThreeSt[1] != 0)coeff[1] = FiveSt[1]/ThreeSt[1]; 
	      else coeff[1] = 0;
#ifdef FFSTIME
	time = -dclock();
#endif
#if 0
	      FORALLSITES(i,s)
		{
		  scalar_mult_add_su3_vector(&(P3[sig][i].h[0]),
					     &(P5nu[i].h[0]), coeff[0],
					     &(P3[sig][i].h[0]));
		  scalar_mult_add_su3_vector(&(P3[sig][i].h[1]),
					     &(P5nu[i].h[1]), coeff[1],
					     &(P3[sig][i].h[1]));
		}
#else
	      scalar_mult_add_lathwvec(P3[sig], P5nu, coeff);
#endif
#ifdef FFSTIME
	      time += dclock();
	      node0_printf("FFSHIFT time4 = %e\n",time);
#endif
	    } /* sig */
	} /* nu */

      /* Now the Lepage term... It is the same as 5-link paths with
	 nu=mu and FiveSt=Lepage. */
      u_shift_hw_fermion(Pmu, Pmumu, OPP_DIR(mu), 
			 &mt[OPP_DIR(mu)], temp_hw[OPP_DIR(mu)] );

      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  P5sig = hw[sig];
	  u_shift_hw_fermion(Pmumu, P5sig, sig, &mt[sig], 
			     temp_hw[sig]);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      add_3f_force_to_mom(P5sig, Pmumu, sig, Lepage);
	    }
	  /* Add the force F_nu the 1(3) link in the path: -     */
	  P5nu = hw[mu];
	  u_shift_hw_fermion(P5sig, P5nu, mu, &mt[mu], temp_hw[mu]);
	  side_link_3f_force(mu, sig, mLepage, Pmu, P5sig, Pmumu, P5nu);
	  /* Add the P5nu vector to P3 */
	  if(ThreeSt[0] != 0)coeff[0] = Lepage[0]/ThreeSt[0];
	  else coeff[0] = 0;
	  if(ThreeSt[1] != 0)coeff[1] = Lepage[1]/ThreeSt[1];
	  else coeff[1] = 0;
#ifdef FFSTIME
	  time = -dclock();
#endif
#if 0
	  FORALLSITES(i,s)
	    {
	      scalar_mult_add_su3_vector(&(P3[sig][i].h[0]),
					 &(P5nu[i].h[0]),coeff[0],
					 &(P3[sig][i].h[0]));
	      scalar_mult_add_su3_vector(&(P3[sig][i].h[1]),
					 &(P5nu[i].h[1]),coeff[1],
					 &(P3[sig][i].h[1]));
	    }
#else
	  scalar_mult_add_lathwvec(P3[sig], P5nu, coeff);
#endif
#ifdef FFSTIME
	      time += dclock();
	      node0_printf("FFSHIFT time4 = %e\n",time);
#endif

	  /* Length 3 paths (Not the Naik term) */
	  /* Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    P3mu = hw[mu];  /* OK to clobber P5nu */
	    u_shift_hw_fermion(P3[sig], P3mu, mu, &mt[mu], temp_hw[mu]);
	    /* The above shift is not needed if mu is backwards */
	  side_link_3f_force(mu, sig, ThreeSt, temp_x, P3[sig], Pmu, P3mu);
	}

      /* Finally the OneLink and the Naik term */
      if(GOES_BACKWARDS(mu))/* Do only the forward terms in the Dslash */
	{
	  /* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
	   * shift.                                                   */
	  /* The one link */
	  add_3f_force_to_mom(Pmu, temp_x, OPP_DIR(mu), OneLink);
	  /* For the same reason Pmumu is the forward double link */

	  /* Popmu is a backward shift */
	  Popmu = hw[mu]; /* OK to clobber P3mu */
	  u_shift_hw_fermion(temp_x, Popmu, mu, &mt[mu], temp_hw[mu]);
	  /* The Naik */
	  /* link no 1: - */
	  add_3f_force_to_mom(Pmumu, Popmu, OPP_DIR(mu), mNaik);
	  /* Pmumumu can overwrite Popmu which is no longer needed */
	  Pmumumu = hw[OPP_DIR(mu)];
	  u_shift_hw_fermion(Pmumu, Pmumumu, OPP_DIR(mu), 
			     &mt[OPP_DIR(mu)], temp_hw[OPP_DIR(mu)] );
	  /* link no 0: + */
	  add_3f_force_to_mom(Pmumumu, temp_x, OPP_DIR(mu), Naik);
	}
      else /* The rest of the Naik terms */
	{
	  Popmu = hw[mu]; /* OK to clobber P3mu */
	  u_shift_hw_fermion(temp_x, Popmu, mu, &mt[mu], temp_hw[mu]);
	  /* link no 2: + */
	  /* Pmumu is double backward shift */
	  add_3f_force_to_mom(Popmu, Pmumu, mu, Naik);
	}
      /* Here we have to do together the Naik term and the one link term */
    }/* mu */

  /* Copy force back and clean up */

  FORALLUPDIR(dir){
    copy_qop_from_matrix(force->f + dir*sites_on_node, tmpmom[dir]);
    destroy_matrix(tmpmom[dir]);
    destroy_matrix(forwardlink[dir]);
  }
  destroy_vector(temp_x1);
  destroy_vector(temp_x2);

  /* Free temporary vectors */
  free(temp_x) ;
  special_free(Pmu) ;
  special_free(Pmumu) ;

  for(mu=0; mu<8; mu++) {
    free(P3[mu]);
    free(P5[mu]);
    special_free(hw[mu]);
    free(temp_hw[mu]);
  }

  FORALLUPDIR(dir){
    free(backwardlink[dir]);
  }

  /* Cleanup gathers */
  for(mu = 0; mu < 8; mu++)
    if(mt[mu] != NULL)cleanup_gather(mt[mu]);
  
  final_flop = (Real)nflop*volume/numnodes();
  dtime += dclock();
#ifdef FFTIME
node0_printf("FFTIME:  time = %e (qop_milc) terms = %d mflops = %e\n",dtime,
	     nsrc, (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif

 info->final_sec = dtime; 
 info->final_flop = final_flop;
 info->status = QOP_SUCCESS; 
} /* QOP_asqtad_force_multi */

/*   Covariant shift of the src fermion field in the direction dir  *
 *  by one unit. The result is stored in dest.                       */ 
static void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) {
  su3_vector *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(GOES_FORWARDS(dir)) /* forward shift */
    {
      mtag = start_gather_field(src, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i,s)
	mult_su3_mat_vec(forwardlink[dir]+i, (su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      cleanup_gather(mtag);
    }
  else /* backward shift */
    {
      tmpvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
      FORALLSITES(i,s)
	mult_adj_su3_mat_vec(forwardlink[OPP_DIR(dir)]+i, &(src[i]), &tmpvec[i]);
      mtag = start_gather_field(tmpvec, sizeof(su3_vector), dir, 
				    EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      /* copy the gen_pt to the dest */
      FORALLSITES(i,s)
	dest[i] = *(su3_vector *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

/*#define FN_DEBUG */

/*  Covariant shift of the src half wilson fermion field in the *
 *  direction dir by one unit.  The result is stored in dest.  */
static void u_shift_hw_fermion(half_wilson_vector *src, 
			half_wilson_vector *dest, 
			int dir, msg_tag **mtag, 
			half_wilson_vector *tmpvec) {

#ifdef FN_DEBUG
  site *s ;
  int i ;
#endif

#ifdef FFSTIME
  double time0, time1;
  time1 = -dclock();
#endif
  /* Copy source to temporary vector */
  memcpy( (char *)tmpvec, (char *)src, 
	 sites_on_node*sizeof(half_wilson_vector) );

  if(*mtag == NULL)
    *mtag = start_gather_field(tmpvec, sizeof(half_wilson_vector), 
			       dir, EVENANDODD, gen_pt[dir]);
  else
    restart_gather_field(tmpvec, sizeof(half_wilson_vector), 
			 dir, EVENANDODD, gen_pt[dir], *mtag);
  wait_gather(*mtag);

#ifdef FFSTIME
  time0 = -dclock();
  time1 -= time0;
#endif

  if(GOES_FORWARDS(dir)) /* forward shift */
    {
#ifdef FN_DEBUG
      FORALLSITES(i,s)
	mult_su3_mat_hwvec( forwardlink[dir]+i, 
			    (half_wilson_vector *)gen_pt[dir][i], 
			    dest + i );
#else
      mult_su3_fieldlink_lathwvec(forwardlink[dir], 
				  (half_wilson_vector **)gen_pt[dir],
				  dest);
#endif
    }
  else /* backward shift */
    {
#ifdef FN_DEBUG
      FORALLSITES(i,s)
	mult_adj_su3_mat_hwvec( backwardlink[OPP_DIR(dir)] + i, 
				(half_wilson_vector *)gen_pt[dir][i], 
				dest + i );
#else
      mult_adj_su3_fieldlink_lathwvec( backwardlink[OPP_DIR(dir)],
				       (half_wilson_vector **)gen_pt[dir],
				       dest);
#endif
    }

#ifdef FFSTIME
  time0 += dclock();
  node0_printf("FFSHIFT time0 = %e\nFFSHIFT time1 = %e\n",time0,time1);
#endif
}

/* Add in contribution to the force */
/* Put antihermitian traceless part into momentum */
static void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,Real coeff) {
  register site *s ;
  register int i ;  
  register Real tmp_coeff ;

  su3_matrix tmat;
  su3_matrix *tmat2;

  if(GOES_BACKWARDS(dir))
    {
      dir = OPP_DIR(dir) ; 
      coeff = -coeff ;
    }
  FORALLSITES(i,s){
    if(s->parity==ODD) 
      tmp_coeff = -coeff ;
    else
      tmp_coeff = coeff ;

    tmat2 = tmpmom[dir] + i;
    su3_projector(&(back[i]), &(forw[i]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff, tmat2 );
  }
}


/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
static void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, 
			 int dir, Real coeff[2]) {
#ifdef FN_DEBUG
  register site *s ;
  register int i ;  
  Real tmp_coeff[2] ;
  su3_matrix tmat, *tmat2;
#endif
  Real my_coeff[2];
  int mydir;
#ifdef FFSTIME
  double time;
  time = -dclock();
#endif

  if(GOES_BACKWARDS(dir))
    {mydir = OPP_DIR(dir); my_coeff[0] = -coeff[0]; my_coeff[1] = -coeff[1];}
  else
    {mydir = dir; my_coeff[0] = coeff[0]; my_coeff[1] = coeff[1]; }
#ifdef FN_DEBUG
  FORALLSITES(i,s){
    if(s->parity==ODD)
      {	tmp_coeff[0] = -my_coeff[0] ; tmp_coeff[1] = -my_coeff[1] ; }
    else
      { tmp_coeff[0] = my_coeff[0] ; tmp_coeff[1] = my_coeff[1] ; }
    tmat2 = tmpmom[mydir] + i;
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[0], tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[1], tmat2 );
  }
#else
  scalar_mult_add_lathwvec_proj_su3mat(tmpmom[mydir], back, forw, my_coeff);
#endif

#ifdef FFSTIME
  time += dclock();
  node0_printf("FFSHIFT time2 = %e\n", time);
#endif
}

/*  This routine is needed in order to add the force on the side link *
 * of the paths in the Asq and Asqtad actions. It gets as inputs the  *
 * direction mu of the side link and the direction nu of the Dslash   *
 * term we are dealing with. Then it takes also 4 fermion fields:     *
 * Path: the piece of the path with no hop in the nu or mu direction  *
 * Path_nu: the piece of the path with a hop in nu  but not in mu     *
 * Path_mu: is Path times the link mu                                 *
 * Path_numu: is Path_nu times the link mu                            */
static void side_link_force(int mu, int nu, Real coeff, 
		     su3_vector *Path   , su3_vector *Path_nu, 
		     su3_vector *Path_mu, su3_vector *Path_numu) {
  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu))
	add_force_to_mom(Path_numu, Path, mu, coeff ) ;
      else
	add_force_to_mom(Path, Path_numu, OPP_DIR(mu), -coeff );/* ? extra - */
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      if(GOES_FORWARDS(nu))
	add_force_to_mom(Path_nu, Path_mu, mu, -coeff) ; /* ? extra - */
      else
	add_force_to_mom(Path_mu, Path_nu, OPP_DIR(mu), coeff) ;
    }
}
 
/*  The 3 flavor version of side_link_force used *
 * to optimize fermion transports                */
static void side_link_3f_force(int mu, int nu, Real coeff[2], 
			half_wilson_vector *Path   , 
			half_wilson_vector *Path_nu, 
			half_wilson_vector *Path_mu, 
			half_wilson_vector *Path_numu) {
  Real m_coeff[2] ;

  m_coeff[0] = -coeff[0] ;
  m_coeff[1] = -coeff[1] ;

  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom(Path_numu, Path, mu, coeff ) ;
      else
	add_3f_force_to_mom(Path,Path_numu,OPP_DIR(mu),m_coeff);/* ? extra - */
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom(Path_nu, Path_mu, mu, m_coeff) ; /* ? extra - */
      else
	add_3f_force_to_mom(Path_mu, Path_nu, OPP_DIR(mu), coeff) ;
    }
}

/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
