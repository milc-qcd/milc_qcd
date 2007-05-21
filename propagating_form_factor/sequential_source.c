/******************************** sequential_source.c *****************/  
/* MIMD version 6 */

/* The src is overwritten
 

  Calculate the sequential source 

   src(x,t) =  \gamma_5 S(x, t = tf) exp( -i px)     t = tf
            = 0   otherwise                          t != tf

  This function assumes that the quark propagator is alraedy smeared.
 

  The source is then used as input to an inverter.

 */

/* MIMD version 6 */

/* Initialize a source for the inverter */

#include "prop_form_includes.h"


void sequential_source(
  field_offset ans,          /* size wilson_vector
				on input: 1st propagator in sequence
				on output: result of inversion */
  field_offset work,         /* for building seq source - size wilson_vector */
  int px, int py, int pz,    /* momentum injected at sequential source */
  int t_final,               /* time slice of sequential source */
  int color,int spin,        /* color and spin of 1st source (unused) */
  Real Kappa,               /* for seq source inversion */
  int inverter_type,         /* HOPILU or CGILU */
  int MaxMR,                 /* MR or CG max iterations per restart */
  int nrestart,              /* MR or CG max restarts */
  Real RsdMR ,              /* MR or CG residual tolerance */
  int p_insert               /* (unused) */
  )
{
  register int i;
  int MinMR;
  register site *s; 
  complex phase_fact ; 
  double theta ; 
  double fact = 2.0*PI/(1.0*nx) ; 
  double tt ; 

  tt =-1*dclock() ; 

	
  /* Zero the src over the lattice *****/
  FORALLSITES(i,s) 
  {
    if( s->t == t_final)
    {   
      mult_by_gamma((wilson_vector *)F_PT(s,ans), 
		    (wilson_vector *)F_PT(s,work), GAMMAFIVE );

      /*** multiply in the momentum phase factor ****/
      theta = fact*(s->x * px + s->y * py + s->z * pz  ) ; 
      phase_fact = cmplx((Real) cos(theta)  , (Real) -sin(theta)) ; 

      c_scale_wilson_vector2( (wilson_vector *)F_PT(s,work), &phase_fact);

    } 
    else
    {
      clear_wvec((wilson_vector *)F_PT(s,work)); 
    } 
  }



  /*
   *   do the inversion
   */

  IF_MASTER 
  {
    printf("Sequential source inversion Kappa = %f \n",Kappa); 
  }

  /* To be sure we "hop" far enough in the inversion */
  /* Note: one iteration advances two time slices 
     so nt/2 covers the lattice twice */
  MinMR = nt;

  /* Load inversion control structure */
  qic_sequential.prec = PRECISION;
  qic_sequential.min = MinMR;
  qic_sequential.max = MaxMR;
  qic_sequential.nrestart = nrestart;
  qic_sequential.resid = RsdMR;
  qic_sequential.start_flag = START_NONZERO_GUESS;
  
#ifdef CLOVER
  /* Load Dirac matrix parameters */
  dcp.Kappa = Kappa;
  dcp.Clov_c = clov_c;
  dcp.U0 = u0;
  
  /* Use specified inverter */
  if(inverter_type == HOPILU)
    {
      qic_sequential.nrestart = 1;   /* No restarts with hopping inverter */
      wilson_invert_site(work, ans,
			 hopilu_cl_site,&qic_sequential,(void *)&dcp);
    }
  else if(inverter_type == BICGILU)
    wilson_invert_site(work, ans
		       bicgilu_cl_site,&qic_sequential,(void *)&dcp);
  else
    {
      printf("sequential_source: ERROR inverter type %d unknown\n",
	     inverter_type);
      terminate(1);
    }

#else
  /* Load Dirac matrix parameters */
  dwp.Kappa = Kappa;

  wilson_invert_site(work, ans,
		     mrilu_w_site,&qic_sequential,(void *)&dwp);
#endif
			    
  /*
   *
   */
  /**  IF_MASTER 
  {
    printf("------------------------------------------\n"); 
    printf("END of the sequential source inversion \n"); 
    printf("------------------------------------------\n"); 
    } **/


  tt += dclock() ; 

  IF_VERBOSE_ON(1)
    printf("Time in sequential_source = %e sec\n",tt) ;




}  /*** end of the sequential source function *****/





