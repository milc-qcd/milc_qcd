/****** fermion_force_asqtad.c  -- ******************/
/* MIMD version 6 */
/* fermion force optimized for the Asqtad action
 * (algorithm previously found in quark_stuff.c)
 * D.T. 1/28/98, starting from gauge_stuff.c
 * K.O. 3/99 Added optimized fattening for Asq actions
 * D.T. 4/99 Combine force calculations for both mass quarks
 * K.O. 4/99 Optimized force for Asq action
 * S.G. 7/01, modified to use t_longlink and t_fatlink
 * C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
 * T.B. 11/01, Added d(M)/d(u0) dependencies for equation of state calc's with
 *             Asqtad action - ATTN: site structure needs 'dfatlink_du0' in
 *             addition to 'fatlink': #define DM_DU0
 *
 * J.O. 3/04 Rearranged loops for optimization
 * J.O. C.D. 3/04 Copied forward links for optimization and 
 *                kept mtags open for restart_gather
 *                Worked with pointers where possible to avoid copying.
 * C.D. 6/04 Corrected memory leak
 * C.D. 3/05 Separated from quark_stuff?.c
 
 * This code combines quark_stuff.c and quark_stuff_tmp.c
 * with DSLASH_TMP_LINKS defined, puts links in field-major variables
 * t_longlink and t_fatlink. Otherwise, puts them in the site structure.
 *
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
 */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED - C. DeTar
 * Fermion force: 253935 for eo_fermion_force()
 * Fermion force: 433968 for eo_fermion_force_3f()
 */

/**#define FFTIME**/
/**#define FFSTIME**/

#include "generic_ks_includes.h"	/* definitions files and prototypes */

/* This routine is valid only for Asqtad, so requires the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) ;
void add_force_to_mom(su3_vector *back, su3_vector *forw, int dir, Real coef);
void side_link_force(int mu, int nu, Real coeff, su3_vector *Path,
		     su3_vector *Path_nu, su3_vector *Path_mu, 
		     su3_vector *Path_numu) ;

void u_shift_hw_fermion(half_wilson_vector *src, 
			half_wilson_vector *dest, int dir ) ;
void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, int dir, Real coeff[2]) ;
void side_link_3f_force(int mu, int nu, Real coeff[2], 
			half_wilson_vector *Path   , 
			half_wilson_vector *Path_nu, 
			half_wilson_vector *Path_mu, 
			half_wilson_vector *Path_numu) ;

#ifdef QSINLINE
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) {\
\
  Real _a0r,_a0i,_a1r,_a1i,_a2r,_a2i;\
  Real _b0r,_b0i,_b1r,_b1i,_b2r,_b2i;\
  \
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[0].c[0].real;  _b0i=(src)->h[0].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[0].c[1].real;  _b1i=(src)->h[0].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[0].c[2].real;  _b2i=(src)->h[0].c[2].imag;\
\
  (dest)->h[0].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[0].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
\
  _a0r=(mat)->e[0][0].real;    _a0i=(mat)->e[0][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[0][1].real;    _a1i=(mat)->e[0][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[0][2].real;    _a2i=(mat)->e[0][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[0].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[0].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
  \
  _a0r=(mat)->e[1][0].real;    _a0i=(mat)->e[1][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[1][1].real;    _a1i=(mat)->e[1][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[1][2].real;    _a2i=(mat)->e[1][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[1].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[1].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
  _a0r=(mat)->e[2][0].real;    _a0i=(mat)->e[2][0].imag;\
  _b0r=(src)->h[1].c[0].real;  _b0i=(src)->h[1].c[0].imag;\
  _a1r=(mat)->e[2][1].real;    _a1i=(mat)->e[2][1].imag;\
  _b1r=(src)->h[1].c[1].real;  _b1i=(src)->h[1].c[1].imag;\
  _a2r=(mat)->e[2][2].real;    _a2i=(mat)->e[2][2].imag;\
  _b2r=(src)->h[1].c[2].real;  _b2i=(src)->h[1].c[2].imag;\
\
  (dest)->h[1].c[2].real = _a0r*_b0r - _a0i*_b0i + _a1r*_b1r - _a1i*_b1i + _a2r*_b2r - _a2i*_b2i;\
  (dest)->h[1].c[2].imag = _a0r*_b0i + _a0i*_b0r + _a1r*_b1i + _a1i*_b1r + _a2r*_b2i + _a2i*_b2r;\
\
}

#define su3_projector_for_inline( a, b, dest ){\
  int _i,_j;\
  Real _tmp,_tmp2;\
    for(_i=0;_i<3;_i++)for(_j=0;_j<3;_j++){\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].real;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].imag;\
	(dest)->e[_i][_j].real = _tmp + _tmp2;\
	_tmp2 = (a)->c[_i].real * (b)->c[_j].imag;\
	_tmp = (a)->c[_i].imag * (b)->c[_j].real;\
	(dest)->e[_i][_j].imag = _tmp - _tmp2;\
    }\
}

#else /* External versions */
#define mult_su3_mat_hwvec_for_inline( mat, src, dest ) mult_su3_mat_hwvec( mat, src, dest ) 
#define su3_projector_for_inline( a, b, dest ) su3_projector( a, b, dest )
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
void eo_fermion_force( Real eps, int nflavors, field_offset x_off ){
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  Real ferm_epsilon, coeff;
  Real OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt ;
  su3_vector *tempvec[8] ;
  su3_vector *temp_x ;
  Real *act_path_coeff;

#ifdef FFTIME
  int nflop = 253935;
  double dtime;

  dtime=-dclock();
#endif
  ferm_epsilon = 2.0*(nflavors/4.0)*eps;
  
  /* Load path coefficients from table */
  act_path_coeff = get_quark_path_coeff();

  /* Path coefficients times fermion epsilon */
  OneLink = act_path_coeff[0]*ferm_epsilon ; 
  Naik    = act_path_coeff[1]*ferm_epsilon ;
  ThreeSt = act_path_coeff[2]*ferm_epsilon ;
  FiveSt  = act_path_coeff[3]*ferm_epsilon ;
  SevenSt = act_path_coeff[4]*ferm_epsilon ;
  Lepage  = act_path_coeff[5]*ferm_epsilon ;
  /* *************************************** */

  /* Initialize the DirectLink flags */
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

  /* Allocate temporary vectors */
  for(mu=0;mu<8;mu++)
    tempvec[mu] = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

  /*copy x_off to a temporary vector */
  temp_x = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
  FORALLSITES(i,s) temp_x[i] = *(su3_vector *)F_PT(s,x_off) ;

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
		  coeff = SevenSt/FiveSt ; 
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
	      coeff = FiveSt/ThreeSt ; 
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
	  coeff = Lepage/ThreeSt ; 
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

  /* Free temporary vectors */
  free(temp_x) ;
  for(mu=0;mu<8;mu++)
    free(tempvec[mu]) ;

#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force(version 6) */
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
/*   Version for two sets of flavors with distinct masses             */
/**********************************************************************/

#define Pmu          hwvec[0] 
#define Pnumu        hwvec[1]
#define Prhonumu     hwvec[2]
#define P7           hwvec[3]
#define P7rho        hwvec[4]              
#define P7rhonu      hwvec[5]
#define P5           hwvec[6]
#define P3           hwvec[7]
#define P5nu         hwvec[3]
#define P3mu         hwvec[3]
#define Popmu        hwvec[4]
#define Pmumumu      hwvec[4]
void eo_fermion_force_3f( Real eps, int nflav1, field_offset x1_off, 
			  int nflav2, field_offset x2_off){
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* 4/15/99 combine force from two different mass quarks, (eg 2+1flavors) */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  Real coeff[2],ferm_epsilon ;
  Real *act_path_coeff;
  Real OneLink[2], Lepage[2], Naik[2], FiveSt[2], ThreeSt[2], SevenSt[2] ;
  Real mNaik[2], mLepage[2], mFiveSt[2], mThreeSt[2], mSevenSt[2] ;
  half_wilson_vector *hwvec[8] ;
  half_wilson_vector *temp_x ;

#ifdef FFTIME
  int nflop = 433968;
  double dtime;

  dtime=-dclock();
#endif
  
  /* Path coefficients times fermion epsilon */

  /* Load path coefficients from table */
  act_path_coeff = get_quark_path_coeff();

  ferm_epsilon = 2.0*(nflav1/4.0)*eps;
  OneLink[0] = act_path_coeff[0]*ferm_epsilon ;
  Naik[0]    = act_path_coeff[1]*ferm_epsilon ; mNaik[0]    = -Naik[0] ;
  ThreeSt[0] = act_path_coeff[2]*ferm_epsilon ; mThreeSt[0] = -ThreeSt[0] ;
  FiveSt[0]  = act_path_coeff[3]*ferm_epsilon ; mFiveSt[0]  = -FiveSt[0] ;
  SevenSt[0] = act_path_coeff[4]*ferm_epsilon ; mSevenSt[0] = -SevenSt[0] ;
  Lepage[0]  = act_path_coeff[5]*ferm_epsilon ; mLepage[0]  = -Lepage[0] ;

  ferm_epsilon = 2.0*(nflav2/4.0)*eps;
  OneLink[1] = act_path_coeff[0]*ferm_epsilon ;
  Naik[1]    = act_path_coeff[1]*ferm_epsilon ; mNaik[1]    = -Naik[1] ;
  ThreeSt[1] = act_path_coeff[2]*ferm_epsilon ; mThreeSt[1] = -ThreeSt[1] ;
  FiveSt[1]  = act_path_coeff[3]*ferm_epsilon ; mFiveSt[1]  = -FiveSt[1] ;
  SevenSt[1] = act_path_coeff[4]*ferm_epsilon ; mSevenSt[1] = -SevenSt[1] ;
  Lepage[1]  = act_path_coeff[5]*ferm_epsilon ; mLepage[1]  = -Lepage[1] ;
  /* *************************************** */

  /* Initialize the DirectLink flags */
  for(mu=0;mu<8;mu++)
    DirectLinks[mu] = 0 ;

  /* Allocate temporary vectors */
  for(mu=0;mu<8;mu++)
    hwvec[mu]= 
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));

  /*copy x_off to a temporary vector */
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  FORALLSITES(i,s)
    {
      temp_x[i].h[0] = *(su3_vector *)F_PT(s,x1_off) ;
      temp_x[i].h[1] = *(su3_vector *)F_PT(s,x2_off) ;
    }

  for(sig=0;sig<8;sig++)
    {
      for(mu=0;mu<8;mu++)if((mu!=sig)&&(mu!=OPP_DIR(sig)))
	{
	  u_shift_hw_fermion(temp_x, Pmu, OPP_DIR(mu));
	  u_shift_hw_fermion(Pmu, P3, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
	      add_3f_force_to_mom(P3, Pmu, sig, mThreeSt) ;
	    }
	  for(nu=0;nu<8;nu++)if((nu!=mu )&&(nu!=OPP_DIR(mu ))&&
				(nu!=sig)&&(nu!=OPP_DIR(sig)))
	    {
	      u_shift_hw_fermion(Pmu, Pnumu, OPP_DIR(nu));
	      u_shift_hw_fermion(Pnumu, P5, sig);
	      if(GOES_FORWARDS(sig))
		{
		  /* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
		  add_3f_force_to_mom(P5, Pnumu, sig, FiveSt);
		}
	      for(rho=0;rho<8;rho++)if((rho!=mu )&&(rho!=OPP_DIR(mu ))&&
				       (rho!=nu )&&(rho!=OPP_DIR(nu ))&&
				       (rho!=sig)&&(rho!=OPP_DIR(sig)))
		{
		  u_shift_hw_fermion(Pnumu, Prhonumu, OPP_DIR(rho));
		  /* Length 7 paths */
		  u_shift_hw_fermion(Prhonumu, P7,sig);
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      add_3f_force_to_mom(P7, Prhonumu, sig, mSevenSt ) ;
		    }
		  /*Add the force F_rho the 2(4) link in the path: +     */
		  u_shift_hw_fermion(P7, P7rho, rho);
		  side_link_3f_force(rho,sig,SevenSt,Pnumu,P7,Prhonumu,P7rho);
		  /* Add the P7rho vector to P5 */
		  coeff[0] = SevenSt[0]/FiveSt[0] ; 
		  coeff[1] = SevenSt[1]/FiveSt[1] ; 
		  FORALLSITES(i,s)
		    {
		      scalar_mult_add_su3_vector(&(P5[i].h[0]),
						 &(P7rho[i].h[0]),coeff[0],
						 &(P5[i].h[0]));
		      scalar_mult_add_su3_vector(&(P5[i].h[1]),
						 &(P7rho[i].h[1]),coeff[1],
						 &(P5[i].h[1]));
		    }
		}/* rho */
	      /* Length 5 paths */
	      /*Add the force F_nu the 1(3) link in the path: -     */
	      u_shift_hw_fermion(P5,P5nu, nu);
	      side_link_3f_force(nu,sig,mFiveSt,Pmu,P5, 
			      Pnumu,P5nu) ;
	      /* Add the P5nu vector to P3 */
	      coeff[0] = FiveSt[0]/ThreeSt[0] ; 
	      coeff[1] = FiveSt[1]/ThreeSt[1] ; 
	      FORALLSITES(i,s)
		{
		  scalar_mult_add_su3_vector(&(P3[i].h[0]),
					     &(P5nu[i].h[0]), coeff[0],
					     &(P3[i].h[0]));
		  scalar_mult_add_su3_vector(&(P3[i].h[1]),
					     &(P5nu[i].h[1]), coeff[1],
					     &(P3[i].h[1]));
		}
	    }/* nu */

	  /* Now the Lepage term... It is the same with 5-link paths with
             nu=mu and FiveSt=Lepage. So Pnumu is really Pmumu */
	  u_shift_hw_fermion(Pmu, Pnumu, OPP_DIR(mu));
	  u_shift_hw_fermion(Pnumu, P5, sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      add_3f_force_to_mom(P5, Pnumu, sig, Lepage) ;
	    }
	  /*Add the force F_nu the 1(3) link in the path: -     */
	  u_shift_hw_fermion(P5,P5nu, mu);
	  side_link_3f_force(mu, sig, mLepage, Pmu, P5, Pnumu, P5nu) ;
	  /* Add the P5nu vector to P3 */
	  coeff[0] = Lepage[0]/ThreeSt[0] ; 
	  coeff[1] = Lepage[1]/ThreeSt[1] ; 
	  FORALLSITES(i,s)
	    {
	      scalar_mult_add_su3_vector(&(P3[i].h[0]),
					 &(P5nu[i].h[0]),coeff[0],
					 &(P3[i].h[0]));
	      scalar_mult_add_su3_vector(&(P3[i].h[1]),
					 &(P5nu[i].h[1]),coeff[1],
					 &(P3[i].h[1]));
	    }

	  /* Length 3 paths (Not the Naik term) */
	  /*Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    u_shift_hw_fermion(P3,P3mu, mu );
	  /* The above shift is not needed if mu is backwards */
	  side_link_3f_force(mu, sig, ThreeSt, temp_x, P3, Pmu, P3mu);

	  /* Finally the OneLink and the Naik term */
	  /* Check if this direction is not already done */
	  if( (!DirectLinks[mu]) ){
	    DirectLinks[mu]=1 ;
	    if(GOES_BACKWARDS(mu))/* Do only the forward terms in the Dslash */
	      {
		/* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
		 * shift.                                                   */
		/* The one link */
		add_3f_force_to_mom(Pmu, temp_x, OPP_DIR(mu), OneLink) ;
		/* For the same reason Pnumu is the forward double link */

		/* Popmu is a backward shift */
		u_shift_hw_fermion(temp_x, Popmu, mu);
		/* The Naik */
		/* link no 1: - */
		add_3f_force_to_mom(Pnumu, Popmu, OPP_DIR(mu), mNaik) ;
		/*Pmumumu can overwrite Popmu which is no longer needed */
		u_shift_hw_fermion(Pnumu, Pmumumu, OPP_DIR(mu));
		/* link no 0: + */
		add_3f_force_to_mom(Pmumumu, temp_x, OPP_DIR(mu), Naik);
	      }
	    else /* The rest of the Naik terms */
	      {
		u_shift_hw_fermion(temp_x, Popmu, mu);
		/* link no 2: + */
		/* Pnumu is double backward shift */
		add_3f_force_to_mom(Popmu, Pnumu, mu, Naik) ;
	      }
	  }
	}/* mu */
      /* Here we have to do together the Naik term and the one link term */
    }/*sig */

  /* Free temporary vectors */
  free(temp_x) ;
  for(mu=0;mu<8;mu++)
    free(hwvec[mu]) ;

#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_3f(version 6) */
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
 
/*   Covariant shift of the src fermion field in the direction dir  *
 *  by one unit. The result is stored in dest.                       */ 
void u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) {
  su3_vector *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(GOES_FORWARDS(dir)) /* forward shift */
    {
      mtag = start_gather_from_temp(src, sizeof(su3_vector), 
				    dir, EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i,s)
	mult_su3_mat_vec(&(s->link[dir]),(su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      cleanup_gather(mtag);
    }
  else /* backward shift */
    {
      tmpvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
      FORALLSITES(i,s)
	mult_adj_su3_mat_vec(&(s->link[OPP_DIR(dir)]),&(src[i]), &tmpvec[i]);
      mtag = start_gather_from_temp(tmpvec, sizeof(su3_vector), dir, 
				    EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      /* copy the gen_pt to the dest */
      FORALLSITES(i,s)
	dest[i] = *(su3_vector *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

/*  Covariant shift of the src half wilson fermion field in the  *
 * direction dir   by one unit. The result is stored in dest.    */
void u_shift_hw_fermion(half_wilson_vector *src, 
			half_wilson_vector *dest, int dir ) {
  half_wilson_vector *tmpvec ; 
  msg_tag *mtag ;
  register site *s ;
  register int i ;
  
  if(GOES_FORWARDS(dir)) /* forward shift */
    {
      mtag = start_gather_from_temp(src, sizeof(half_wilson_vector), 
				    dir, EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      FORALLSITES(i,s)
	mult_su3_mat_hwvec(&(s->link[dir]),
			   (half_wilson_vector *)(gen_pt[0][i]), &(dest[i]));
      cleanup_gather(mtag);
    }
  else /* backward shift */
    {
      tmpvec = 
	(half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
      FORALLSITES(i,s)
	mult_adj_su3_mat_hwvec(&(s->link[OPP_DIR(dir)]),&(src[i]), &tmpvec[i]);
      mtag = start_gather_from_temp(tmpvec, sizeof(half_wilson_vector), dir, 
				    EVENANDODD, gen_pt[0]);
      wait_gather(mtag);
      /* copy the gen_pt to the dest */
      FORALLSITES(i,s)
	dest[i] = *(half_wilson_vector *)gen_pt[0][i];
      cleanup_gather(mtag);
      free(tmpvec) ;
    }
}

/* Add in contribution to the force */
/* Put antihermitian traceless part into momentum */
void add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,Real coeff) {
  register site *s ;
  register int i ;  
  register Real tmp_coeff ;

  su3_matrix tmat,tmat2;

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
    uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
    su3_projector(&(back[i]), &(forw[i]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff, &tmat2 );
    make_anti_hermitian( &tmat2, &(s->mom[dir]) ); 
  }
}

/* Add in contribution to the force ( 3flavor case ) */
/* Put antihermitian traceless part into momentum */
void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, int dir, Real coeff[2]) {
  register site *s ;
  register int i ;  
  Real tmp_coeff[2] ;

  su3_matrix tmat,tmat2;

  if(GOES_BACKWARDS(dir))
    {
      dir = OPP_DIR(dir) ; 
      coeff[0] = -coeff[0] ;
      coeff[1] = -coeff[1] ;
    }
  FORALLSITES(i,s){
    if(s->parity==ODD)
      {
	tmp_coeff[0] = -coeff[0] ;
	tmp_coeff[1] = -coeff[1] ;
      }
    else
      {
	tmp_coeff[0] = coeff[0] ;
	tmp_coeff[1] = coeff[1] ;
      }
    uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff[0], &tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat,  tmp_coeff[1], &tmat2 );
    make_anti_hermitian( &tmat2, &(s->mom[dir]) ); 
  }
}

/*  This routine is needed in order to add the force on the side link *
 * of the paths in the Asq and Asqtad actions. It gets as inputs the  *
 * direction mu of the side link and the direction nu of the Dslash   *
 * term we are dealing with. Then it takes also 4 fermion fields:     *
 * Path: the piece of the path with no hop in the nu or mu direction  *
 * Path_nu: the piece of the path with a hop in nu  but not in mu     *
 * Path_mu: is Path times the link mu                                 *
 * Path_numu: is Path_nu times the link mu                            */
void side_link_force(int mu, int nu, Real coeff, 
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
void side_link_3f_force(int mu, int nu, Real coeff[2], 
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
