/****** fermion_force_asqtad.c  -- ******************/
/* MIMD version 7 */
/* fermion force optimized for the Asqtad action 
 * formerly fermion_force_asqtad3.c
 * Uses restart-gathers and a bit more memory for better performance
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
 * D.T. 12/05 First try at RHMC version
 * C.D. 5/07  Renamed from fermion_force_asqtad3.c
 *            Collected all asqtad: single, double, multi in this file.
 *
 * In this directory, assume all paths connect even to odd sites, etc.
 * Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
 * as separate parity transforms of the fundamental paths.  They will
 * generally need a negative sign in Dslash.  See bottom for a long
 * comment on sign conventions.
 */

/* External entry points in this file
   
    eo_fermion_force_oneterm
    eo_fermion_force_twoterms
    fermion_force_asqtad_block
    fermion_force_asqtad_multi
   
 */

/* Compile with fermion_force_fn_multi.c */

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED - C. DeTar
 * Fermion force: 253935 for eo_fermion_force_oneterm()
 * Fermion force: 433968 for eo_fermion_force_twoterms()
 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include <string.h>

/* This routine is valid only for Asqtad, so requires the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

/* Forward declarations */

static void 
u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) ;

static void 
add_force_to_mom(su3_vector *back, su3_vector *forw, int dir,  Real coef);

static void 
side_link_force(int mu, int nu, Real coeff, su3_vector *Path,
		su3_vector *Path_nu, su3_vector *Path_mu, 
		su3_vector *Path_numu) ;

static void 
u_shift_hw_fermion(half_wilson_vector *src, half_wilson_vector *dest, 
		   int dir, msg_tag** mtag, half_wilson_vector *tmpvec ) ;

static void 
add_3f_force_to_mom(half_wilson_vector *back, half_wilson_vector *forw, 
		    int dir, Real coeff[2]) ;

static void 
side_link_3f_force(int mu, int nu, Real coeff[2], half_wilson_vector *Path   , 
		   half_wilson_vector *Path_nu, half_wilson_vector *Path_mu, 
		   half_wilson_vector *Path_numu) ;

void 
fermion_force_asqtad_multi( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, int prec,
			    fermion_links_t *fl );

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
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

void 
eo_fermion_force_oneterm( Real eps, Real weight, su3_vector *temp_x,
			  int prec, fermion_links_t *fl )
{
  /* prec is ignored for now */
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  /* For example weight = nflavors/4 */
  ks_action_paths *ap = get_action_paths(fl);
  register int i ;
  register site *s;
  int mu,nu,rho,sig ;
  int DirectLinks[8] ;
  Real ferm_epsilon, coeff;
  Real OneLink, Lepage, Naik, FiveSt, ThreeSt, SevenSt ;
  su3_vector *tempvec[8] ;

#ifdef FFTIME
  int nflop = 253935;
  double dtime;

  dtime=-dclock();
#endif
  ferm_epsilon = 2.0*weight*eps;
  
  /*  node0_printf("STARTING fn_fermion_force_oneterm() nterms = 1\n");*/

  /* Load path coefficients from table */
  //act_path_coeff = ap->p.act_path_coeff;

  /* Path coefficients times fermion epsilon */
  OneLink = ap->p.act_path_coeff.one_link*ferm_epsilon ; 
  Naik    = ap->p.act_path_coeff.naik*ferm_epsilon ;
  ThreeSt = ap->p.act_path_coeff.three_staple*ferm_epsilon ;
  FiveSt  = ap->p.act_path_coeff.five_staple*ferm_epsilon ;
  SevenSt = ap->p.act_path_coeff.seven_staple*ferm_epsilon ;
  Lepage  = ap->p.act_path_coeff.lepage*ferm_epsilon ;
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

  /* Free temporary vectors */
  for(mu=0;mu<8;mu++)
    free(tempvec[mu]) ;

#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e (asqtad3) terms = 1 mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
/**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_oneterm (version 7) */
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

void 
eo_fermion_force_oneterm_site( Real eps, Real weight, field_offset x_off,
			       int prec, fermion_links_t *fl )
{
  su3_vector *temp_x  = create_v_field_from_site_member(x_off);
  eo_fermion_force_oneterm(eps, weight, temp_x, prec, fl );
  
  destroy_v_field(temp_x);
}

/**********************************************************************/
/*   Version for two sets of flavors with distinct masses             */
/**********************************************************************/

static su3_matrix *backwardlink[4];
static anti_hermitmat *tempmom[4];

static void 
eo_fermion_force_twoterms_hwv( Real eps, Real weight1, Real weight2, 
			       half_wilson_vector *temp_x, int prec,
			       fermion_links_t *fl )
{
   /* prec is ignored for now */
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* 4/15/99 combine force from two different mass quarks, (eg 2+1flavors) */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  /* For example weight1 = nflavor1/4; weight2 = nflavor2/4 */
  ks_action_paths *ap = get_action_paths(fl);
  int i;
  site *s;
  int mu,nu,rho,sig;
  int dir;
  Real coeff[2],ferm_epsilon;
  Real OneLink[2], Lepage[2], Naik[2], FiveSt[2], ThreeSt[2], SevenSt[2];
  Real mNaik[2], mLepage[2], mFiveSt[2], mThreeSt[2], mSevenSt[2];
  half_wilson_vector *Pnumu, *Prhonumu, *P7, *P7rho, *P5nu, 
    *P3mu = NULL, *P5sig, *Popmu, *Pmumumu;
  half_wilson_vector *P3[8];
  half_wilson_vector *P5[8];
  half_wilson_vector *Pmu;
  half_wilson_vector *Pmumu;
  half_wilson_vector *temp_hw[8];
  msg_tag *mt[8];
  msg_tag *mtag[4];
  half_wilson_vector *hw[8];
  char myname[] = "eo_fermion_force_twoterms_hwv";

#ifdef FFSTIME
  double time;
#endif

#ifdef FFTIME
  int nflop = 433968;
  double dtime;

  dtime=-dclock();
#endif

  /*  node0_printf("STARTING fn_fermion_force_twoterms() nterms = 2\n");*/

  /* Allocate temporary hw vector */
  for(mu = 0; mu < 8; mu++){
    half_wilson_vector *pt;
    pt = (half_wilson_vector *)
      special_alloc(sites_on_node*sizeof(half_wilson_vector));
    if(pt == NULL){
      printf("%s: No room for hw\n", myname);
      terminate(1);
    }
    hw[mu] = pt;
  }

  Pmu = 
    (half_wilson_vector *)special_alloc(sites_on_node*sizeof(half_wilson_vector));
  if(Pmu == NULL){
    printf("%s: No room for Pmu\n", myname);
    terminate(1);
  }
  
  Pmumu = 
    (half_wilson_vector *)special_alloc(sites_on_node*sizeof(half_wilson_vector));
  if(Pmumu == NULL){
    printf("%s: No room for Pmumu\n", myname);
    terminate(1);
  }
  
  /* Allocate temporary vectors */
  for(mu=0; mu<8; mu++) {
    P3[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(P3[mu] == NULL){
      printf("%s: No room for P3\n", myname);
      terminate(1);
    }
    P5[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(P5[mu] == NULL){
      printf("%s: No room for P5\n", myname);
      terminate(1);
    }
  }

  /* Initialize message pointers */
  for(mu = 0; mu < 8; mu++)mt[mu] = NULL;

  /* Double store backward gauge links */
  FORALLUPDIR(dir){
    su3_matrix *pt;
    pt = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(pt == NULL){
      printf("%s: No room for backwardlink\n", myname);
      terminate(1);
    }
    backwardlink[dir] = pt;
  }
  
  /* Gather backward links */
  FORALLUPDIR(dir){
    mtag[dir] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix), 
			      OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
  }
  FORALLUPDIR(dir){
    wait_gather(mtag[dir]);
    FORALLSITES(i,s){
      backwardlink[dir][i] = *((su3_matrix *)gen_pt[dir][i]);
    }
    cleanup_gather(mtag[dir]);
  }
  
  /* Copy gauge momenta */
  FORALLUPDIR(dir){
    anti_hermitmat *pt;
    pt = (anti_hermitmat *)malloc(sites_on_node*sizeof(anti_hermitmat));
    if(pt == NULL){
      printf("%s: No room for tempmom\n", myname);
      terminate(1);
    }
    tempmom[dir] = pt;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      tempmom[dir][i] = s->mom[dir];
    }
  }
  
  /* Load path coefficients from table */
  // act_path_coeff = ap->p.act_path_coeff;

  /* Path coefficients times fermion epsilon */
  ferm_epsilon = 2.0*weight1*eps;
  OneLink[0] = ap->p.act_path_coeff.one_link*ferm_epsilon;
  Naik[0]    = ap->p.act_path_coeff.naik*ferm_epsilon; mNaik[0]    = -Naik[0];
  ThreeSt[0] = ap->p.act_path_coeff.three_staple*ferm_epsilon; mThreeSt[0] = -ThreeSt[0];
  FiveSt[0]  = ap->p.act_path_coeff.five_staple*ferm_epsilon; mFiveSt[0]  = -FiveSt[0];
  SevenSt[0] = ap->p.act_path_coeff.seven_staple*ferm_epsilon; mSevenSt[0] = -SevenSt[0];
  Lepage[0]  = ap->p.act_path_coeff.lepage*ferm_epsilon; mLepage[0]  = -Lepage[0];

  ferm_epsilon = 2.0*weight2*eps;
  OneLink[1] = ap->p.act_path_coeff.one_link*ferm_epsilon;
  Naik[1]    = ap->p.act_path_coeff.naik*ferm_epsilon; mNaik[1]    = -Naik[1];
  ThreeSt[1] = ap->p.act_path_coeff.three_staple*ferm_epsilon; mThreeSt[1] = -ThreeSt[1];
  FiveSt[1]  = ap->p.act_path_coeff.five_staple*ferm_epsilon; mFiveSt[1]  = -FiveSt[1];
  SevenSt[1] = ap->p.act_path_coeff.seven_staple*ferm_epsilon; mSevenSt[1] = -SevenSt[1];
  Lepage[1]  = ap->p.act_path_coeff.lepage*ferm_epsilon; mLepage[1]  = -Lepage[1];
  /* *************************************** */

  /* Allocate temporary vectors */

  for(mu = 0; mu < 8; mu++){
    temp_hw[mu] = 
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(temp_hw[mu] == NULL){
      printf("%s: No room for temp_hw\n", myname);
      terminate(1);
    }
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

  /* Repack momenta */

  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      s->mom[dir] = tempmom[dir][i];
    }
  }

  /* Free temporary vectors */
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
    free(tempmom[dir]);
  }

  /* Cleanup gathers */
  for(mu = 0; mu < 8; mu++)
    if(mt[mu] != NULL)cleanup_gather(mt[mu]);
  
#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e (asqtad3) terms = 2 mflops = %e\n",dtime,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
  /**printf("TLENGTH: %d\n",tlength);**/
#endif

} /* eo_fermion_force_twoterms */

void 
eo_fermion_force_twoterms( Real eps, Real weight1, Real weight2, 
			   su3_vector *x1_off, su3_vector *x2_off,
			   int prec, fermion_links_t *fl ){

  half_wilson_vector *temp_x;
  int i;
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  if(temp_x == NULL){
    printf("fermion_force_twoterms: No room for temporary\n");
    terminate(1);
  }
  
  FORALLFIELDSITES(i)
    {
      temp_x[i].h[0] = x1_off[i] ;
      temp_x[i].h[1] = x2_off[i] ;
    }
  
  eo_fermion_force_twoterms_hwv(eps, weight1, weight2, temp_x, prec, fl);

  free(temp_x);
}  



void 
eo_fermion_force_twoterms_site( Real eps, Real weight1, Real weight2, 
				field_offset x1_off, field_offset x2_off,
				int prec, fermion_links_t *fl )
{
  int i;
  site *s;
  half_wilson_vector *temp_x;

  /* copy x_off to a temporary vector */
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  if(temp_x == NULL){
    printf("eo_fermion_force_twoterms_site: No room for temporary\n");
    terminate(1);
  }

  FORALLSITES(i,s)
    {
      temp_x[i].h[0] = *(su3_vector *)F_PT(s,x1_off);
      temp_x[i].h[1] = *(su3_vector *)F_PT(s,x2_off);
    }
  
  eo_fermion_force_twoterms_hwv(eps, weight1, weight2, temp_x, prec, fl );
  
  free(temp_x) ;
}  


/**********************************************************************/
/*   Parallel transport vectors in blocks of veclength. Asqtad only   */
/**********************************************************************/
void 
fermion_force_asqtad_block( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, int veclength, 
			    int prec, fermion_links_t *fl ) 
{

  int j;

  /* First do blocks of size veclength */
  for( j = 0;  j <= nterms-veclength; j += veclength )
    fermion_force_asqtad_multi( eps, &(residues[j]), xxx+j, veclength, prec, fl );
  
  /* Continue with pairs if needed */
  if(j <= nterms-2){
//    half_wilson_vector *temp_x;
//    temp_x= 
//      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
//    if(temp_x == NULL){
//      printf("fermion_force_asqtad_block: No room for temporary\n");
//      terminate(1);
//    }
    
    for( ; j <= nterms-2 ; j+=2 ){
      
//      FORALLSITES(i,s)
//	{
//	  temp_x[i].h[0] = xxx[j  ][i] ;
//	  temp_x[i].h[1] = xxx[j+1][i] ;
//	}
      
      eo_fermion_force_twoterms( eps, residues[j], residues[j+1],
				 xxx[j], xxx[j + 1], prec, fl );
    }
    //    free(temp_x);
  }

  /* Finish with a single if needed */
  for( ; j <= nterms-1; j++ ){
    eo_fermion_force_oneterm( eps, residues[j], xxx[j], prec, fl );
  }
}

/**********************************************************************/
/*   Version for asqtad.  Parallel transport nterms source vectors    */
/**********************************************************************/

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED
 * Fermion force: 21698*nterms
 */

static void 
u_shift_veclist_fermion(veclist *src, veclist *dest, 
			int dir, msg_tag **mtag, veclist *tmpvec, 
			int listlength );

static void 
add_3f_force_to_mom_list(veclist *back,	 veclist *forw, int dir, 
			 Real coeff[2], int listlength ) ;

static void 
side_link_3f_force_list(int mu, int nu, Real *coeff, veclist *Path, 
			veclist *Path_nu, veclist *Path_mu, 
			veclist *Path_numu, int listlength ) ;

static su3_matrix *backwardlink[4];
static anti_hermitmat *tempmom[4];

void 
fermion_force_asqtad_multi( Real eps, Real *residues, 
			    su3_vector **xxx, int nterms, int prec,
			    fermion_links_t *fl ) 
{
  /* prec is ignored for now */
  // note CG_solution and Dslash * solution are combined in "*xxx"
  // New version 1/21/99.  Use forward part of Dslash to get force
  // see long comment at end
  // For each link we need x_off transported from both ends of path.
  // xxx[n][site] is the n'th inverse solution at location index "site"
  //
  // combine an arbitrary number of force terms, with weights "residues[j]"
  ks_action_paths *ap = get_action_paths(fl);
  int i,j;
  site *s;
  int mu,nu,rho,sig;
  int dir;
  Real *coeff,ferm_epsilon;
  Real *OneLink, *Lepage, *Naik, *FiveSt, *ThreeSt, *SevenSt;
  Real *mNaik, *mLepage, *mFiveSt, *mThreeSt, *mSevenSt;
  veclist *Pnumu, *Prhonumu, *P7, *P7rho, *P5nu, *P3mu = NULL, *P5sig, *Popmu;
  veclist *Pmu, *Pmumu, *Pmumumu;
  veclist *P3[8], *P5[8];
  veclist *temp_x;
  veclist *vl[8], *temp_vl[8];
  msg_tag *mt[8];
  msg_tag *mtag[4];
  char myname[] = "fermion_force_asqtad_multi";

#ifdef FFSTIME
  double time;
#endif

#ifdef FFTIME
  int nflop;
  double dtime;

  nflop = 216984*nterms;
  dtime=-dclock();
#endif

  /*node0_printf("STARTING %s nterms = %d\n",myname,nterms);*/
  if( nterms==0 )return;

  if(nterms > VECLENGTH){
    printf("%s(%d): too many terms %d > %d\n",myname,this_node,
	   nterms,VECLENGTH);
    terminate(1);
  }

  /* Allocate temporary hw vector */
  for(mu = 0; mu < 8; mu++){
    veclist *pt;
    pt = (veclist *) special_alloc(sites_on_node*sizeof(veclist));
    if(pt == NULL){
      printf("%s No room for pt\n",myname);
      terminate(1);
    }
    vl[mu] = pt;
  }

  Pmu = (veclist *)special_alloc(sites_on_node*sizeof(veclist));
  if(Pmu == NULL){
    printf("%s: No room for Pmu\n",myname);
    terminate(1);
  }
  
  Pmumu = (veclist *)special_alloc(sites_on_node*sizeof(veclist));
  if(Pmumu == NULL){
    printf("%s: No room for Pmumu\n",myname);
    terminate(1);
  }
  
  /* Allocate temporary vectors */
  for(mu=0; mu<8; mu++) {
    P3[mu]= (veclist *)malloc(sites_on_node*sizeof(veclist));
    if(P3[mu] == NULL){
      printf("%s: No room for P3\n",myname);
      terminate(1);
    }
    P5[mu]= (veclist *)malloc(sites_on_node*sizeof(veclist));
    if(P5[mu] == NULL){
      printf("%s: No room for P5\n",myname);
      terminate(1);
    }
  }

  /* Initialize message pointers */
  for(mu = 0; mu < 8; mu++)mt[mu] = NULL;

  /* Double store backward gauge links */
  FORALLUPDIR(dir){
    su3_matrix *pt;
    pt = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(pt == NULL){
      printf("%s: No room for backwardlink\n",myname);
      terminate(1);
    }
    backwardlink[dir] = pt;
  }
  
  /* Gather backward links */
  FORALLUPDIR(dir){
    mtag[dir] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix), 
			      OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
  }
  FORALLUPDIR(dir){
    wait_gather(mtag[dir]);
    FORALLSITES(i,s){
      backwardlink[dir][i] = *((su3_matrix *)gen_pt[dir][i]);
    }
    cleanup_gather(mtag[dir]);
  }
  
  /* Copy gauge momenta */
  FORALLUPDIR(dir){
    anti_hermitmat *pt;
    pt = (anti_hermitmat *)malloc(sites_on_node*sizeof(anti_hermitmat));
    if(pt == NULL){
      printf("%s: No room for tempmom\n",myname);
      terminate(1);
    }
    tempmom[dir] = pt;
  }
  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      tempmom[dir][i] = s->mom[dir];
    }
  }
  
  /* Load path coefficients from table */
  // act_path_coeff = ap->p.act_path_coeff;

  /* Path coefficients times fermion epsilon */
  coeff = (Real *)malloc( nterms*sizeof(Real) );
  OneLink = (Real *)malloc( nterms*sizeof(Real) );
  Naik = (Real *)malloc( nterms*sizeof(Real) ); mNaik = (Real *)malloc( nterms*sizeof(Real) );
  ThreeSt = (Real *)malloc( nterms*sizeof(Real) ); mThreeSt = (Real *)malloc( nterms*sizeof(Real) );
  FiveSt = (Real *)malloc( nterms*sizeof(Real) ); mFiveSt = (Real *)malloc( nterms*sizeof(Real) );
  SevenSt = (Real *)malloc( nterms*sizeof(Real) ); mSevenSt = (Real *)malloc( nterms*sizeof(Real) );
  Lepage = (Real *)malloc( nterms*sizeof(Real) ); mLepage = (Real *)malloc( nterms*sizeof(Real) );
  for( j=0; j<nterms;j++ ){
    ferm_epsilon = 2.0*residues[j]*eps;
    OneLink[j] = ap->p.act_path_coeff.one_link*ferm_epsilon;
    Naik[j]    = ap->p.act_path_coeff.naik*ferm_epsilon; mNaik[j]    = -Naik[j];
    ThreeSt[j] = ap->p.act_path_coeff.three_staple*ferm_epsilon; mThreeSt[j] = -ThreeSt[j];
    FiveSt[j]  = ap->p.act_path_coeff.five_staple*ferm_epsilon; mFiveSt[j]  = -FiveSt[j];
    SevenSt[j] = ap->p.act_path_coeff.seven_staple*ferm_epsilon; mSevenSt[j] = -SevenSt[j];
    Lepage[j]  = ap->p.act_path_coeff.lepage*ferm_epsilon; mLepage[j]  = -Lepage[j];
  }

  /* *************************************** */

  /* Allocate temporary vectors */

  for(mu = 0; mu < 8; mu++){
    temp_vl[mu] = 
      (veclist *)malloc(sites_on_node*sizeof(veclist));
    if(temp_vl[mu] == NULL){
      printf("%s: No room for temp_vl\n",myname);
      terminate(1);
    }
  }

  /* copy x_off to a temporary vector */
  temp_x= 
    (veclist *)malloc(sites_on_node*sizeof(veclist));
  FORALLSITES(i,s) {
      for( j=0; j<nterms; j++ ){
         temp_x[i].v[j] = xxx[j][i];
      }
    }

  for(mu=0; mu<8; mu++)
    {
      u_shift_veclist_fermion( temp_x, Pmu, OPP_DIR(mu), 
	 &mt[OPP_DIR(mu)],  temp_vl[OPP_DIR(mu)], nterms );

      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  u_shift_veclist_fermion( Pmu, P3[sig], sig, &mt[sig],  temp_vl[sig], nterms );

	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu]:         x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 1 link in the path: - (numbering starts form 0) */
	      add_3f_force_to_mom_list( P3[sig], Pmu, sig, mThreeSt, nterms );
	    }
	}
      for(nu=0; nu<8; nu++) if( (nu!=mu)&&(nu!=OPP_DIR(mu)) )
	{
	  Pnumu = vl[OPP_DIR(nu)];
	  u_shift_veclist_fermion( Pmu, Pnumu, OPP_DIR(nu), &mt[OPP_DIR(nu)], 
	     temp_vl[OPP_DIR(nu)], nterms );
	  for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
				       (sig!=nu)&&(sig!=OPP_DIR(nu)) )
	    {
	      u_shift_veclist_fermion( Pnumu, P5[sig], sig, &mt[sig], 
			  temp_vl[sig], nterms );

	      if(GOES_FORWARDS(sig))
		{
		  /* Add the force F_sig[x+mu+nu]:      x--+             *
		   *                                   |   |             *
		   *                                   o   o             *
		   * the 2 link in the path: + (numbering starts form 0) */
		  add_3f_force_to_mom_list( P5[sig], Pnumu, sig, FiveSt, nterms );
		}
	    }
	  for(rho=0; rho<8; rho++) if( (rho!=mu)&&(rho!=OPP_DIR(mu)) &&
				       (rho!=nu)&&(rho!=OPP_DIR(nu)) )
	    {
	      Prhonumu = vl[OPP_DIR(rho)];
	      u_shift_veclist_fermion( Pnumu, Prhonumu, OPP_DIR(rho), 
		 &mt[OPP_DIR(rho)],  temp_vl[OPP_DIR(rho)], nterms );
	      for(sig=0; sig<8; sig++) if( (sig!=mu )&&(sig!=OPP_DIR(mu )) &&
					   (sig!=nu )&&(sig!=OPP_DIR(nu )) &&
					   (sig!=rho)&&(sig!=OPP_DIR(rho)) )
		{
		  /* Length 7 paths */
		  P7 = vl[sig];
		  u_shift_veclist_fermion( Prhonumu, P7, sig, &mt[sig], 
		      temp_vl[sig], nterms );
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      add_3f_force_to_mom_list( P7, Prhonumu, sig, mSevenSt, nterms );
		    }
		  /* Add the force F_rho the 2(4) link in the path: +     */
		  P7rho = vl[rho];
		  u_shift_veclist_fermion( P7, P7rho, rho, &mt[rho], 
			     temp_vl[rho], nterms );
		  side_link_3f_force_list( rho, sig, SevenSt, Pnumu, P7, Prhonumu, P7rho, nterms );

		  /* Add the P7rho vector to P5 */
		  for( j=0; j<nterms; j++ ){
		     if(FiveSt[j] != 0.0)coeff[j] = SevenSt[j]/FiveSt[j];
		     else coeff[j] = 0.0;
		  }
#ifdef FFSTIME
	time = -dclock();
#endif
		  FORALLSITES(i,s) {
		     for( j=0; j< nterms; j++ ){
		      scalar_mult_add_su3_vector(&(P5[sig][i].v[j]),
			   &(P7rho[i].v[j]),coeff[j], &(P5[sig][i].v[j]));
		     }
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
	      P5nu = vl[nu];
	      u_shift_veclist_fermion( P5[sig], P5nu, nu, &mt[nu], temp_vl[nu], nterms );
	      side_link_3f_force_list(nu, sig, mFiveSt, Pmu, P5[sig], Pnumu, P5nu, nterms );
	      /* Add the P5nu vector to P3 */
	      for( j=0; j<nterms; j++ ){
	         if(ThreeSt[j] != 0.0)coeff[j] = FiveSt[j]/ThreeSt[j]; 
	         else coeff[j] = 0.0;
	      }
#ifdef FFSTIME
	time = -dclock();
#endif
	      scalar_mult_add_latveclist( P3[sig], P5nu, coeff, nterms);
#ifdef FFSTIME
	      time += dclock();
	      node0_printf("FFSHIFT time4 = %e\n",time);
#endif
	    } /* sig */
	} /* nu */

      /* Now the Lepage term... It is the same as 5-link paths with
	 nu=mu and FiveSt=Lepage. */
      u_shift_veclist_fermion( Pmu, Pmumu, OPP_DIR(mu), 
		 &mt[OPP_DIR(mu)], temp_vl[OPP_DIR(mu)], nterms );

      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  P5sig = vl[sig];
	  u_shift_veclist_fermion( Pmumu, P5sig, sig, &mt[sig], 
		 temp_vl[sig], nterms );
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      add_3f_force_to_mom_list( P5sig, Pmumu, sig, Lepage, nterms );
	    }
	  /* Add the force F_nu the 1(3) link in the path: -     */
	  P5nu = vl[mu];
	  u_shift_veclist_fermion( P5sig, P5nu, mu, &mt[mu], temp_vl[mu], nterms );
	  side_link_3f_force_list(mu, sig, mLepage, Pmu, P5sig, Pmumu, P5nu, nterms );
	  /* Add the P5nu vector to P3 */
	  for( j=0; j<nterms; j++ ){
	     if(ThreeSt[j] != 0.0)coeff[j] = Lepage[j]/ThreeSt[j];
	     else coeff[j] = 0.0;
	  }
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
	  scalar_mult_add_latveclist( P3[sig], P5nu, coeff, nterms );
#endif
#ifdef FFSTIME
	      time += dclock();
	      node0_printf("FFSHIFT time4 = %e\n",time);
#endif

	  /* Length 3 paths (Not the Naik term) */
	  /* Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    P3mu = vl[mu];  /* OK to clobber P5nu */
	    u_shift_veclist_fermion( P3[sig], P3mu, mu, &mt[mu], temp_vl[mu], nterms );
	    /* The above shift is not needed if mu is backwards */
	  side_link_3f_force_list(mu, sig, ThreeSt, temp_x, P3[sig], Pmu, P3mu, nterms );
	}

      /* Finally the OneLink and the Naik term */
      if(GOES_BACKWARDS(mu))/* Do only the forward terms in the Dslash */
	{
	  /* Because I have shifted with OPP_DIR(mu) Pmu is a forward *
	   * shift.                                                   */
	  /* The one link */
	  add_3f_force_to_mom_list( Pmu, temp_x, OPP_DIR(mu), OneLink, nterms );
	  /* For the same reason Pmumu is the forward double link */

	  /* Popmu is a backward shift */
	  Popmu = vl[mu]; /* OK to clobber P3mu */
	  u_shift_veclist_fermion( temp_x, Popmu, mu, &mt[mu], temp_vl[mu], nterms );
	  /* The Naik */
	  /* link no 1: - */
	  add_3f_force_to_mom_list( Pmumu, Popmu, OPP_DIR(mu), mNaik, nterms );
	  /* Pmumumu can overwrite Popmu which is no longer needed */
	  Pmumumu = vl[OPP_DIR(mu)];
	  u_shift_veclist_fermion( Pmumu, Pmumumu, OPP_DIR(mu), 
		&mt[OPP_DIR(mu)], temp_vl[OPP_DIR(mu)], nterms );
	  /* link no 0: + */
	  add_3f_force_to_mom_list( Pmumumu, temp_x, OPP_DIR(mu), Naik, nterms );
	}
      else /* The rest of the Naik terms */
	{
	  Popmu = vl[mu]; /* OK to clobber P3mu */
	  u_shift_veclist_fermion( temp_x, Popmu, mu, &mt[mu], temp_vl[mu], nterms );
	  /* link no 2: + */
	  /* Pmumu is double backward shift */
	  add_3f_force_to_mom_list( Popmu, Pmumu, mu, Naik, nterms );
	}
      /* Here we have to do together the Naik term and the one link term */
    }/* mu */

  /* Repack momenta */

  FORALLUPDIR(dir){
    FORALLSITES(i,s){
      s->mom[dir] = tempmom[dir][i];
    }
  }

  /* Free temporary vectors */
  free(temp_x) ;
  special_free(Pmu) ;
  special_free(Pmumu) ;

  for(mu=0; mu<8; mu++) {
    free(P3[mu]);
    free(P5[mu]);
    special_free(vl[mu]);
    free(temp_vl[mu]);
  }

  FORALLUPDIR(dir){
    free(backwardlink[dir]);
    free(tempmom[dir]);
  }

  free(coeff);
  free(OneLink);
  free(ThreeSt); free(mThreeSt);
  free(FiveSt); free(mFiveSt);
  free(SevenSt); free(mSevenSt);
  free(Lepage); free(mLepage);

  /* Cleanup gathers */
  for(mu = 0; mu < 8; mu++)
    if(mt[mu] != NULL)cleanup_gather(mt[mu]);
  
#ifdef FFTIME
  dtime += dclock();
node0_printf("FFTIME:  time = %e (ASVEC) terms = %d mflops = %e\n",dtime,
	     nterms,(Real)nflop*volume/(1e6*dtime*numnodes()) );
  /**printf("TLENGTH: %d\n",tlength);**/
#endif


} /* fermion_force_asqtad_multi */

/*   Covariant shift of the src fermion field in the direction dir  *
 *  by one unit. The result is stored in dest.                       */ 
static void 
u_shift_fermion(su3_vector *src, su3_vector *dest, int dir ) {
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
	mult_su3_mat_vec(&(s->link[dir]),(su3_vector *)(gen_pt[0][i]),
			 &(dest[i]));
      cleanup_gather(mtag);
    }
  else /* backward shift */
    {
      tmpvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );
      FORALLSITES(i,s)
	mult_adj_su3_mat_vec(&(s->link[OPP_DIR(dir)]),&(src[i]), &tmpvec[i]);
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

/*  Covariant shift of the src half wilson fermion field in the *
 *  direction dir by one unit.  The result is stored in dest.  */
static void 
u_shift_hw_fermion(half_wilson_vector *src, half_wilson_vector *dest, 
		   int dir, msg_tag **mtag, half_wilson_vector *tmpvec) {
#if 0
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
#if 0
      FORALLSITES(i,s)
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector *)gen_pt[dir][i], 
			    dest + i );
#else
      mult_su3_sitelink_lathwvec(dir, (half_wilson_vector **)gen_pt[dir],
				 dest);
#endif
    }
  else /* backward shift */
    {
#if 0
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

static void 
add_force_to_mom(su3_vector *back,su3_vector *forw, int dir,Real coeff) {

  register site *s ;
  register int i ;  
  register Real tmp_coeff ;

  su3_matrix tmat, tmat2;

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
static void add_3f_force_to_mom(half_wilson_vector *back,
			 half_wilson_vector *forw, 
			 int dir, Real coeff[2]) {
#if 0
  register site *s ;
  register int i ;  
  Real tmp_coeff[2] ;
  su3_matrix tmat, *tmat2;
#endif
  Real my_coeff[2] ;
  int mydir;
#ifdef FFSTIME
  double time;
  time = -dclock();
#endif

  if(GOES_BACKWARDS(dir))
    {mydir = OPP_DIR(dir); my_coeff[0] = -coeff[0]; my_coeff[1] = -coeff[1];}
  else
    {mydir = dir; my_coeff[0] = coeff[0]; my_coeff[1] = coeff[1]; }
#if 0
  FORALLSITES(i,s){
    if(s->parity==ODD)
      {	tmp_coeff[0] = -my_coeff[0] ; tmp_coeff[1] = -my_coeff[1] ; }
    else
      { tmp_coeff[0] = my_coeff[0] ; tmp_coeff[1] = my_coeff[1] ; }
    tmat2 = tempmom[mydir] + i;
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[0], tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[1], tmat2 );
  }
#else
  scalar_mult_add_lathwvec_proj(tempmom[mydir], back, forw, my_coeff);
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

static void 
side_link_force(int mu, int nu, Real coeff, su3_vector *Path, 
		su3_vector *Path_nu, su3_vector *Path_mu, 
		su3_vector *Path_numu) {

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

static void 
side_link_3f_force(int mu, int nu, Real coeff[2], half_wilson_vector *Path   , 
		   half_wilson_vector *Path_nu, half_wilson_vector *Path_mu, 
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

static void 
u_shift_veclist_fermion(veclist *src, veclist *dest, int dir, 
			msg_tag **mtag, veclist *tmpvec, int listlength ) {
#if 0
  site *s ;
  int i ;
#endif

#ifdef FFSTIME
  double time0, time1;
  time1 = -dclock();
#endif
  /* Copy source to temporary vector */
  memcpy( (char *)tmpvec, (char *)src, 
	 sites_on_node*sizeof(veclist) );

  if(*mtag == NULL)
    *mtag = start_gather_field(tmpvec, sizeof(veclist), 
			       dir, EVENANDODD, gen_pt[dir]);
  else
    restart_gather_field(tmpvec, sizeof(veclist), 
			 dir, EVENANDODD, gen_pt[dir], *mtag);
  wait_gather(*mtag);

#ifdef FFSTIME
  time0 = -dclock();
  time1 -= time0;
#endif

  if(GOES_FORWARDS(dir)) /* forward shift */
    {
#if 0
      FORALLSITES(i,s)
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector *)gen_pt[dir][i], 
			    dest + i );
#else
      mult_su3_sitelink_latveclist(dir, (veclist **)gen_pt[dir], dest, listlength );
#endif
    }
  else /* backward shift */
    {
#if 0
      FORALLSITES(i,s)
	mult_adj_su3_mat_hwvec( backwardlink[OPP_DIR(dir)] + i, 
				(half_wilson_vector *)gen_pt[dir][i], 
				dest + i );
#else
     // mult_adj_su3_fieldlink_lathwvec( backwardlink[OPP_DIR(dir)],
				       //(half_wilson_vector **)gen_pt[dir],
				       //dest);
      mult_adj_su3_fieldlink_latveclist( backwardlink[OPP_DIR(dir)],
	(veclist **)gen_pt[dir], dest, listlength );
#endif
    }

#ifdef FFSTIME
  time0 += dclock();
  node0_printf("FFSHIFT time0 = %e\nFFSHIFT time1 = %e\n",time0,time1);
#endif
}



static void 
add_3f_force_to_mom_list(veclist *back, veclist *forw, int dir, 
			 Real *coeff, int listlength ) {

#if 0
  register site *s ;
  register int i;
  su3_matrix tmat, *tmat2;
  Real *tmp_coeff = (Real *)malloc(listlength*sizeof(Real) );
#endif
  int j ;  
  Real *my_coeff = (Real *)malloc(listlength*sizeof(Real) );
  int mydir;
#ifdef FFSTIME
  double time;
  time = -dclock();
#endif

  if(GOES_BACKWARDS(dir)) {
    mydir = OPP_DIR(dir); for( j=0; j<listlength; j++ ){my_coeff[j] = -coeff[j];}
  }
  else {
    mydir = dir; for( j=0; j<listlength; j++ ){my_coeff[j] = coeff[j];}
  }
#if 0
  FORALLSITES(i,s){
    if(s->parity==ODD)
      {	tmp_coeff[0] = -my_coeff[0] ; tmp_coeff[1] = -my_coeff[1] ; }
    else
      { tmp_coeff[0] = my_coeff[0] ; tmp_coeff[1] = my_coeff[1] ; }
    tmat2 = tempmom[mydir] + i;
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[0], tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[1], tmat2 );
  }
#else
  scalar_mult_add_latveclist_proj( tempmom[mydir], back,
	forw, my_coeff, listlength );
#endif

#ifdef FFSTIME
  time += dclock();
  node0_printf("FFSHIFT time2 = %e\n", time);
#endif
  free(my_coeff); 
#if 0
  free(tmp_coeff);
#endif
}

static void 
side_link_3f_force_list(int mu, int nu, Real *coeff, veclist *Path, 
			veclist *Path_nu, veclist *Path_mu, 
			veclist *Path_numu, int listlength ) {

  Real *m_coeff;
  int j;

  m_coeff = (Real *)malloc( listlength*sizeof(Real) );
  for( j=0; j<listlength; j++){
     m_coeff[j] = -coeff[j] ;
  }

  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom_list( Path_numu, Path, mu, coeff, listlength ) ;
      else
	add_3f_force_to_mom_list( Path, Path_numu,OPP_DIR(mu),m_coeff, listlength);/* ? extra - */
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom_list( Path_nu, Path_mu, mu, m_coeff, listlength) ; /* ? extra - */
      else
	add_3f_force_to_mom_list( Path_mu, Path_nu, OPP_DIR(mu), coeff, listlength) ;
    }
    free(m_coeff);
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
