#define Pmu          hwvec[0] 
#define Pnumu        hwvec[1]
#define Prhonumu     hwvec[2]
#define P7           hwvec[3]
#define P7rho        hwvec[4]              
#define P7rhonu      hwvec[5]
//#define P5           hwvec[6]
//#define P3           hwvec[7]
#define P5nu         hwvec[3]
#define P3mu         hwvec[3]
#define Popmu        hwvec[4]
#define Pmumumu      hwvec[4]
void eo_fermion_force_two( Real eps, Real res1, Real res2,
			   field_offset x1_off, field_offset x2_off,
			   ferm_links_t *fn, ks_action_paths *ap ) {
  /* note CG_solution and Dslash * solution are combined in "x_off" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* 4/15/99 combine force from two different mass quarks, (eg 2+1flavors) */
  /* see long comment at end */
  /* For each link we need x_off transported from both ends of path. */
  int i;
  site *s;
  int mu,nu,rho,sig;
  Real coeff[2],ferm_epsilon;
  Real OneLink[2], Lepage[2], Naik[2], FiveSt[2], ThreeSt[2], SevenSt[2];
  Real mNaik[2], mLepage[2], mFiveSt[2], mThreeSt[2], mSevenSt[2];
  half_wilson_vector *hwvec[6];
  half_wilson_vector *P3[8];
  half_wilson_vector *P5[8];
  half_wilson_vector *temp_x;

#ifdef FFTIME
  double dtime;

  dtime=-dclock();
#endif

  /* Path coefficients times fermion epsilon */
  ferm_epsilon = 2.0*res1*eps;
  OneLink[0] = act_path_coeff[0]*ferm_epsilon;
  Naik[0]    = act_path_coeff[1]*ferm_epsilon; mNaik[0]    = -Naik[0];
  ThreeSt[0] = act_path_coeff[2]*ferm_epsilon; mThreeSt[0] = -ThreeSt[0];
  FiveSt[0]  = act_path_coeff[3]*ferm_epsilon; mFiveSt[0]  = -FiveSt[0];
  SevenSt[0] = act_path_coeff[4]*ferm_epsilon; mSevenSt[0] = -SevenSt[0];
  Lepage[0]  = act_path_coeff[5]*ferm_epsilon; mLepage[0]  = -Lepage[0];

  ferm_epsilon = 2.0*res2*eps;
  OneLink[1] = act_path_coeff[0]*ferm_epsilon;
  Naik[1]    = act_path_coeff[1]*ferm_epsilon; mNaik[1]    = -Naik[1];
  ThreeSt[1] = act_path_coeff[2]*ferm_epsilon; mThreeSt[1] = -ThreeSt[1];
  FiveSt[1]  = act_path_coeff[3]*ferm_epsilon; mFiveSt[1]  = -FiveSt[1];
  SevenSt[1] = act_path_coeff[4]*ferm_epsilon; mSevenSt[1] = -SevenSt[1];
  Lepage[1]  = act_path_coeff[5]*ferm_epsilon; mLepage[1]  = -Lepage[1];
  /* *************************************** */

  /* Allocate temporary vectors */
  for(mu=0; mu<6; mu++) {
    hwvec[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  }
  for(mu=0; mu<8; mu++) {
    P3[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    P5[mu]=
      (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  }

  /* copy x_off to a temporary vector */
  temp_x= 
    (half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
  FORALLSITES(i,s)
    {
      temp_x[i].h[0] = *(su3_vector *)F_PT(s,x1_off);
      temp_x[i].h[1] = *(su3_vector *)F_PT(s,x2_off);
    }

  for(mu=0; mu<8; mu++)
    {
      u_shift_hw_fermion(temp_x, Pmu, OPP_DIR(mu));
      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  u_shift_hw_fermion(Pmu, P3[sig], sig);
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
	  u_shift_hw_fermion(Pmu, Pnumu, OPP_DIR(nu));
	  for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
				       (sig!=nu)&&(sig!=OPP_DIR(nu)) )
	    {
	      u_shift_hw_fermion(Pnumu, P5[sig], sig);
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
	      u_shift_hw_fermion(Pnumu, Prhonumu, OPP_DIR(rho));
	      for(sig=0; sig<8; sig++) if( (sig!=mu )&&(sig!=OPP_DIR(mu )) &&
					   (sig!=nu )&&(sig!=OPP_DIR(nu )) &&
					   (sig!=rho)&&(sig!=OPP_DIR(rho)) )
		{
		  /* Length 7 paths */
		  u_shift_hw_fermion(Prhonumu, P7, sig);
		  if(GOES_FORWARDS(sig))
		    {
		      /* Add the force F_sig[x+mu+nu+rho]:  x--+             *
		       *                                   |   |             *
		       *                                   o   o             *
		       * the 3 link in the path: - (numbering starts form 0) */
		      add_3f_force_to_mom(P7, Prhonumu, sig, mSevenSt);
		    }
		  /* Add the force F_rho the 2(4) link in the path: +     */
		  u_shift_hw_fermion(P7, P7rho, rho);
		  side_link_3f_force(rho,sig,SevenSt,Pnumu,P7,Prhonumu,P7rho);
		  /* Add the P7rho vector to P5 */
		  coeff[0] = SevenSt[0]/FiveSt[0];
		  coeff[1] = SevenSt[1]/FiveSt[1];
		  FORALLSITES(i,s)
		    {
		      scalar_mult_add_su3_vector(&(P5[sig][i].h[0]),
						 &(P7rho[i].h[0]),coeff[0],
						 &(P5[sig][i].h[0]));
		      scalar_mult_add_su3_vector(&(P5[sig][i].h[1]),
						 &(P7rho[i].h[1]),coeff[1],
						 &(P5[sig][i].h[1]));
		    }
		} /* sig */
	    } /* rho */
	  for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) &&
				       (sig!=nu)&&(sig!=OPP_DIR(nu)) )
	    {
	      /* Length 5 paths */
	      /* Add the force F_nu the 1(3) link in the path: -     */
	      u_shift_hw_fermion(P5[sig], P5nu, nu);
	      side_link_3f_force(nu, sig, mFiveSt, Pmu, P5[sig], Pnumu, P5nu);
	      /* Add the P5nu vector to P3 */
	      coeff[0] = FiveSt[0]/ThreeSt[0]; 
	      coeff[1] = FiveSt[1]/ThreeSt[1]; 
	      FORALLSITES(i,s)
		{
		  scalar_mult_add_su3_vector(&(P3[sig][i].h[0]),
					     &(P5nu[i].h[0]), coeff[0],
					     &(P3[sig][i].h[0]));
		  scalar_mult_add_su3_vector(&(P3[sig][i].h[1]),
					     &(P5nu[i].h[1]), coeff[1],
					     &(P3[sig][i].h[1]));
		}
	    } /* sig */
	} /* nu */

      /* Now the Lepage term... It is the same with 5-link paths with
	 nu=mu and FiveSt=Lepage. So Pnumu is really Pmumu */
      u_shift_hw_fermion(Pmu, Pnumu, OPP_DIR(mu));
      for(sig=0; sig<8; sig++) if( (sig!=mu)&&(sig!=OPP_DIR(mu)) )
	{
	  u_shift_hw_fermion(Pnumu, P5[sig], sig);
	  if(GOES_FORWARDS(sig))
	    {
	      /* Add the force F_sig[x+mu+nu]:      x--+             *
	       *                                   |   |             *
	       *                                   o   o             *
	       * the 2 link in the path: + (numbering starts form 0) */
	      add_3f_force_to_mom(P5[sig], Pnumu, sig, Lepage);
	    }
	  /* Add the force F_nu the 1(3) link in the path: -     */
	  u_shift_hw_fermion(P5[sig], P5nu, mu);
	  side_link_3f_force(mu, sig, mLepage, Pmu, P5[sig], Pnumu, P5nu);
	  /* Add the P5nu vector to P3 */
	  coeff[0] = Lepage[0]/ThreeSt[0];
	  coeff[1] = Lepage[1]/ThreeSt[1];
	  FORALLSITES(i,s)
	    {
	      scalar_mult_add_su3_vector(&(P3[sig][i].h[0]),
					 &(P5nu[i].h[0]),coeff[0],
					 &(P3[sig][i].h[0]));
	      scalar_mult_add_su3_vector(&(P3[sig][i].h[1]),
					 &(P5nu[i].h[1]),coeff[1],
					 &(P3[sig][i].h[1]));
	    }

	  /* Length 3 paths (Not the Naik term) */
	  /* Add the force F_mu the 0(2) link in the path: +     */
	  if(GOES_FORWARDS(mu)) 
	    u_shift_hw_fermion(P3[sig], P3mu, mu);
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
	  /* For the same reason Pnumu is the forward double link */

	  /* Popmu is a backward shift */
	  u_shift_hw_fermion(temp_x, Popmu, mu);
	  /* The Naik */
	  /* link no 1: - */
	  add_3f_force_to_mom(Pnumu, Popmu, OPP_DIR(mu), mNaik);
	  /* Pmumumu can overwrite Popmu which is no longer needed */
	  u_shift_hw_fermion(Pnumu, Pmumumu, OPP_DIR(mu));
	  /* link no 0: + */
	  add_3f_force_to_mom(Pmumumu, temp_x, OPP_DIR(mu), Naik);
	}
      else /* The rest of the Naik terms */
	{
	  u_shift_hw_fermion(temp_x, Popmu, mu);
	  /* link no 2: + */
	  /* Pnumu is double backward shift */
	  add_3f_force_to_mom(Popmu, Pnumu, mu, Naik);
	}
      /* Here we have to do together the Naik term and the one link term */
    }/* mu */

  /* Free temporary vectors */
  free(temp_x) ;
  for(mu=0; mu<6; mu++) {
    free(hwvec[mu]);
  }
  for(mu=0; mu<8; mu++) {
    free(P3[mu]);
    free(P5[mu]);
  }

#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  %e\n",dtime);
  /**printf("TLENGTH: %d\n",tlength);**/
#endif
} /* eo_fermion_force_two(version 6) */
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
