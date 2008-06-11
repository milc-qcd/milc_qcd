/******** spectrum_hybrids5.c *************/
/* MIMD version 7 */
/* DT 11/21/95 started */
/* DT 2/21/96 version 2.  Make it easy to add operators.  Keep collection
   of "mult_by_XXX" routines around.  Select only operators that look
   interesting in previous tests. 
   DT 7/10/96 multiple sink operators for each source.
   New output format:  "SINK:" line lists sink operators, 
    "operator_SRC:" lines give propagators from each source.
   DT 10/4/96 add 0+- "baryon number" operator
   DT 1/21/97 use defines to select desired source operators, add
	1-+ 4 quark operator QQQQ

   cmn 11/05/97 I have added calls to the full meson and baryon
                spectrum. I have also modified the code so that clover
                inverter can be slotted in.
*/



/* Spectrum for Wilson hybrid mesons. 
   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie boundary_flip(MINUS) )

   Also does spectrum for a few conventional mesons, as reference.

    Gauge fixing should be done before calling this function
    (Coulomb gauge, probably).
*/


/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* operators:
   pion:	conventional 0-+:  qbar gamma_5 q
   pion2:	conventional 0-+:  qbar gamma_0 gamma_5 q
   0mp:		0-+ hyrid
   rho:		conventional 1--: qbar gamma_i q
   1mm:		1-- hybrid
   a1:		conventional 1++: qbar gamma_5 gamma_i q
   a1P:		P wave a_1: qbar eps_ijk gamma_j partial_k q
   1pp:		1++ hybrid
   0pm:		0+- hybrid
   0pmP		0+- P-wave hybrid
   0pmB		0+- bilinear, qbar gamma_0 q
   0mm		0-- hybrid
   0mmP:	0-- P-wave hybrid
   1mp:		1-+ hybrid, qbar gamma_i q eps_ijk B_j
   1mp2:	1-+ hybrid, qbar gamma_0 q E_i
   qqqq	1-+, a_1 plus pion
*/
enum prop_name { 
    prop_pion_pion,
    prop_pion_pion2,
    prop_pion_0mp,
    prop_pion2_pion,
    prop_pion2_pion2,
    prop_pion2_0mp,
    prop_0mp_pion,
    prop_0mp_pion2,
    prop_0mp_0mp,

    prop_rho_rho,
    prop_rho_rho2,
    prop_rho_1mm,
    prop_rho2_rho,
    prop_rho2_rho2,
    prop_rho2_1mm,
    prop_1mm_rho,
    prop_1mm_rho2,
    prop_1mm_1mm,

    prop_a1_a1,
    prop_a1_a1P,
    prop_a1_1pp,
    prop_a1P_a1,
    prop_a1P_a1P,
    prop_a1P_1pp,
    prop_1pp_a1,
    prop_1pp_a1P,
    prop_1pp_1pp,

    prop_0pm_0pm,
    prop_0pm_0pmP,
    prop_0pm_0pmB,
    prop_0pmP_0pm,
    prop_0pmP_0pmP,
    prop_0pmP_0pmB,
    prop_0pmB_0pm,
    prop_0pmB_0pmP,
    prop_0pmB_0pmB,

    prop_0mm_0mm,
    prop_0mmP_0mmP,
    prop_0mm_0mmP,
    prop_0mmP_0mm,

    prop_1mp_1mp,
    prop_1mp2_1mp2,
    prop_1mp2_1mp,
    prop_1mp_1mp2,
    prop_qqqq_1mp,
    prop_qqqq_1mp2,
prop_qqqq_junk,
prop_qqqq_junk2,

    prop_2pm_2pm,

    nprops		/* nprops = number of propagators */
};

#include "cl_hyb_includes.h"

/* Various meson operators.  Source code later in this file */
/* All of these operators contain an extra gamma_5, entering when
    antiquark propagators are computed from the same source as quark
    propagators */
   /* quark-antiquark operators*/
void mult_pion( field_offset src, field_offset dest );
void mult_pion2( field_offset src, field_offset dest );
void mult_rho( int pdir, field_offset src, field_offset dest );
void mult_rho2( int pdir, field_offset src, field_offset dest );
void mult_a1( int pdir, field_offset src, field_offset dest );
void mult_a1_P( int pdir, field_offset src, field_offset dest );
void mult_b1( int pdir, field_offset src, field_offset dest );
  /* exotic hybrids, and one exotic bilinear */
void mult_zero_pm( field_offset src, field_offset dest );
void mult_zero_pm_P( field_offset src, field_offset dest );
void mult_zero_pm_B( field_offset src, field_offset dest );
void mult_zero_mm( field_offset src, field_offset dest );
void mult_zero_mm_P( field_offset src, field_offset dest );
void mult_one_mp( int pdir, field_offset src, field_offset dest );
void mult_one_mp2( int pdir, field_offset src, field_offset dest );
void mult_qqqq1mp( int src_t, int pdir, field_offset src, field_offset dest,
   field_offset work );
   /* non-exotic hybrids*/
void mult_zero_mp( field_offset src, field_offset dest );
void mult_one_mm( int pdir, field_offset src, field_offset dest );
void mult_one_pp( int pdir, field_offset src, field_offset dest );

void mult_by_field_strength( int dir1, int dir2,
    field_offset src, field_offset dest );
void smear_links( field_offset src, field_offset dest);

void parallel_transport( int dir, field_offset src, field_offset dest,
    field_offset work, int gdir );
void check_invert_cl( field_offset src, field_offset dest );
int test_converge(int t_source);
void copy_site_wilson_vector(field_offset src, field_offset dest) ;


void light_meson_spectrum(int t_source) ;
void light_baryon_spectrum(field_offset quark, int t_source) ;


int spectrum_hybrids(){ /* return the C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc;

  register int t_source;
  int dir;	/* direction in lattice */
  int spin;	/* spin for source */
  int color;	/* color for source */
  int src_count; /* number of source time slices used */
  complex **props;	/* arrays of propagators */

  cgn=0; /* number of CG iterations */

  /* allocate arrays to accumulate propagators */
  props = (complex **)malloc(nprops*sizeof(complex *));
  props[0] = (complex *)malloc(nprops*nt*sizeof(complex));
  for(i=1;i<nprops;i++)props[i]=props[i-1]+nt;

  /* compute the field strength tensor, uses smeared links if SMEAR
     is defined */

if( smearing_level != 0 ) 
{
  if(this_node==0)printf("SMEARING IS ON, level %d\n",smearing_level);
  if(smearing_level >=1)smear_links( F_OFFSET(link[0]), F_OFFSET(smearlink[0]) );
  else FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
    s->smearlink[dir] = s->link[dir];
  }
  for(j=2;j<= smearing_level   ;j++){
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      s->templink[dir] = s->smearlink[dir];
    }
    smear_links( F_OFFSET(templink[0]), F_OFFSET(smearlink[0]) );
  }
}  /** end of smearing level <> 0 **/
else
{
  if(this_node==0)printf("NO SMEARING OF f_mu_mu in the hybrid operators\n") ;
}



#ifdef SMEAR
  make_field_strength(F_OFFSET(smearlink[0]), F_OFFSET(field_strength[0]));
#else
  make_field_strength(F_OFFSET(link[0]), F_OFFSET(field_strength[0]));
#endif

  /**temp_test();**/

  /* loop over "source" time slice */
  for(src_count=0,t_source=source_start; t_source<nt && src_count<n_sources;
    t_source += source_inc,src_count++){
    /* Wall source */
    /* Use quark_source for quark source */
    if(this_node==0)printf("spectrum_hybrids(): source time = %d\n",t_source);

    /* set propagators to zero */
    for(cc.real=cc.imag=0.0,i=0;i<nprops;i++)for(j=0;j<nt;j++){
      props[i][j]=cc;
    }

    for(spin=0;spin<4;spin++)
      for(color=0;color<3;color++)
      {

	switch (wot_src) 
	{
	case WALL_QUARK_SRC : 

	  FORALLSITES(i,s)
	  {
	    clear_wvec( &(s->quark_source) );
            if(s->t==t_source) s->quark_source.d[spin].c[color].real=1.0;
	  }
	  break ; 

	case LOCAL_QUARK_SRC :

	  FORALLSITES(i,s)
	  {
	    clear_wvec( &(s->quark_source) );
            if(s->t==t_source && s->x == 0 && s->y == 0 && s->z == 0 ) s->quark_source.d[spin].c[color].real=1.0;
	  }
	  break ; 



	default :
	  if( this_node == 0 ) printf("ERROR: wot_src flag = %d is out of range\n",wot_src) ; 
	  terminate(1) ; 
	  break ; 
	}



        /* compute M^-1 * quark_source */
        cgn += mat_invert( F_OFFSET(quark_source), F_OFFSET(quark_prop) );
	/*if(t_source==0)test_converge(t_source);*/ /*TEMP*/
	/* TEMP: test inversion, */
	/**check_invert_cl( F_OFFSET(quark_prop), F_OFFSET(quark_source) );**/

        /* Begin with the pion source and sink operator */
        /* Use mat_invert to construct anti_prop */
	/* ACTUALLY, nothing to do for pion - anti_prop = quark_prop */
        /* tie propagators together at sink end to project out desired mesons */
        /* add into propagator arrays at distance from source */


	if( oper_PION_SOURCE == CALCULATE ) 
	{
	/* 0-+ (pion) PROPAGATORS ***************************************/
        /**
        mult_pion( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_pion( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion_pion][(s->t+nt-t_source)%nt], cc );
        }
        **/
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->quark_prop), &(s->quark_prop) );
	    CSUM( props[prop_pion_pion][(s->t+nt-t_source)%nt], cc );
        }
        mult_pion2( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->quark_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion_pion2][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_mp( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->quark_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion_0mp][(s->t+nt-t_source)%nt], cc );
        }
} /*PION_SOURCE*/

	if( oper_PION2_SOURCE == CALCULATE ) 
	{
	/* now the pion2 source */
        mult_pion2( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_pion( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion2_pion][(s->t+nt-t_source)%nt], cc );
        }
        mult_pion2( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion2_pion2][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_mp( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_pion2_0mp][(s->t+nt-t_source)%nt], cc );
        }
      } /*PION2_SOURCE*/

if ( oper_ZEROMP_SOURCE == CALCULATE ) 
{
        /* The 0-+ hybrid source.  */
        mult_zero_mp( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_pion( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mp_pion][(s->t+nt-t_source)%nt], cc );
        }
        mult_pion2( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mp_pion2][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_mp( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mp_0mp][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROMP_SOURCE*/

if(  oper_RHO_SOURCE == CALCULATE ) 
{
	/* 1-- (rho) PROPAGATORS ***************************************/
        /* Now the rho source */
        mult_rho( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_rho( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho_rho][(s->t+nt-t_source)%nt], cc );
        }
        mult_rho2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho_rho2][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mm( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho_1mm][(s->t+nt-t_source)%nt], cc );
        }
} /*RHO_SOURCE*/

if( oper_RHO2_SOURCE == CALCULATE ) 
{
        /* Now the rho2 source */
        mult_rho2( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_rho( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho2_rho][(s->t+nt-t_source)%nt], cc );
        }
        mult_rho2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho2_rho2][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mm( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_rho2_1mm][(s->t+nt-t_source)%nt], cc );
        }
} /*RHO2_SOURCE*/

if(  oper_ONEMM_SOURCE == CALCULATE ) 
{
        /* Now do the 1-- _1-- hybrid source. 
		For the moment, Z component only */
        mult_one_mm( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_rho( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mm_rho][(s->t+nt-t_source)%nt], cc );
        }
        mult_rho2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mm_rho2][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mm( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mm_1mm][(s->t+nt-t_source)%nt], cc );
        }
} /*ONEMM_SOURCE*/

if(  oper_A1_SOURCE == CALCULATE ) 
{
	/* 1++ (a1) PROPAGATORS ***************************************/
        /* Now the a1 source (1++, gamma_5 gamma_z) */
        mult_a1( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_a1( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1_a1][(s->t+nt-t_source)%nt], cc );
        }
        mult_a1_P( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1_a1P][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_pp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1_1pp][(s->t+nt-t_source)%nt], cc );
        }
} /*A1_SOURCE*/

if(  oper_A1P_SOURCE == CALCULATE ) 
{
        /* The a1P source (1++, epsilon_ijk gamma_j deriv_k ) */
        mult_a1_P( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_a1( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1P_a1][(s->t+nt-t_source)%nt], cc );
        }
        mult_a1_P( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1P_a1P][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_pp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_a1P_1pp][(s->t+nt-t_source)%nt], cc );
        }
} /*A1P_SOURCE*/

if(  oper_ONEPP_SOURCE == CALCULATE ) 
{
        /* the 1++  hybrid source.  For the moment, Z component only */
        mult_one_pp( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_a1( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1pp_a1][(s->t+nt-t_source)%nt], cc );
        }
        mult_a1_P( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1pp_a1P][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_pp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1pp_1pp][(s->t+nt-t_source)%nt], cc );
        }
} /*ONEPP_SOURCE*/

if(  oper_ZEROPM_SOURCE == CALCULATE ) 
{
	/* 0+- (exotic) PROPAGATORS ***************************************/
        /* Now the 0+-_0+- propagator */
        mult_zero_pm( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_pm( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pm_0pm][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pm_0pmP][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_B( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pm_0pmB][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROPM_SOURCE*/

if( oper_ZEROPMP_SOURCE == CALCULATE ) 
{
        /* The 0+-P exotic hybrid P wave source */
        mult_zero_pm_P( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_pm( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmP_0pm][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmP_0pmP][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_B( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmP_0pmB][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROPMP_SOURCE*/

if( oper_ZEROPMB_SOURCE == CALCULATE ) 
{
        /* The 0+- exotic quark bilinear source */
        mult_zero_pm_B( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_pm( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmB_0pm][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmB_0pmP][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_pm_B( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0pmB_0pmB][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROPMB_SOURCE*/

if(  oper_ZEROMM_SOURCE == CALCULATE ) 
{
	/* 0-- (exotic) PROPAGATORS ***************************************/
        /* Now do the 0--_0-- propagator */
        mult_zero_mm( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_mm( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mm_0mm][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_mm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mm_0mmP][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROMM_SOURCE*/

if(  oper_ZEROMMP_SOURCE == CALCULATE ) 
{
        /* The 0--P_0--P exotic hybrid P wave propagator */
        mult_zero_mm_P( F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_zero_mm( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mmP_0mm][(s->t+nt-t_source)%nt], cc );
        }
        mult_zero_mm_P( F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_0mmP_0mmP][(s->t+nt-t_source)%nt], cc );
        }
} /*ZEROMMP_SOURCE*/

if ( oper_ONEMP_SOURCE == CALCULATE ) 
{
	/* 1-+ (exotic) PROPAGATORS ***************************************/
        /* Now do the 1-+_1-+ source.  For the moment, Z component only */
        mult_one_mp( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_one_mp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
	/* Take both 1mp and 1mp2 sinks */
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mp_1mp][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mp2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mp_1mp2][(s->t+nt-t_source)%nt], cc );
        }
} /*ONEMP_SOURCE*/

if( oper_ONEMP2_SOURCE == CALCULATE ) 
{
        /* Now the 1-+2_1+-2 source.  For the moment, Z component only */
        mult_one_mp2( ZUP, F_OFFSET(quark_source), F_OFFSET(G_RAND) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        mult_one_mp2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mp2_1mp2][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_1mp2_1mp][(s->t+nt-t_source)%nt], cc );
        }
} /*ONEMP2_SOURCE*/

if(  oper_QQQQ_SOURCE == CALCULATE ) 
{
        /* Now the 4-quark (a1-pi) source.  For the moment, Z component only */
        mult_qqqq1mp( t_source,  ZUP,
	    F_OFFSET(quark_source), F_OFFSET(G_RAND), F_OFFSET(anti_prop) );
        cgn += mat_invert( F_OFFSET(G_RAND), F_OFFSET(anti_prop) );

        mult_one_mp2( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_qqqq_1mp2][(s->t+nt-t_source)%nt], cc );
        }
        mult_one_mp( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_qqqq_1mp][(s->t+nt-t_source)%nt], cc );
        }
#ifdef NONSENSE
        mult_a1( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_qqqq_junk][(s->t+nt-t_source)%nt], cc );
        }
        mult_b1( ZUP, F_OFFSET(quark_prop), F_OFFSET(G_RAND) );
        FORALLSITES(i,s){
	    cc = wvec_dot( &(s->anti_prop), &(s->G_RAND) );
	    CSUM( props[prop_qqqq_junk2][(s->t+nt-t_source)%nt], cc );
        }
#endif /*NONSENSE*/
} /*QQQQ_SOURCE*/

	
	/* 2+- (exotic) PROPAGATORS ***************************************/

	/** store the source colour/spin part of the propagator for future use ***/
	copy_site_wilson_vector(F_OFFSET(quark_prop),F_OFFSET(quark_store.c[color] .d[spin]) ) ;


      }  /*** end the loop over colour and spin **/

    /* Sum propagator arrays over nodes */
    /* print out propagators */
    g_veccomplexsum( props[0] , nprops*nt );
    for(i=0;i<nprops;i++)for(j=0;j<nt;j++){
      CDIVREAL(props[i][j],nx*ny*nz,props[i][j]);
    }
    if(this_node==0){
  

if(  oper_PION_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: PION\n");
      printf("SINKS: PION PION2 0MP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_pion_pion][j].real, props[prop_pion_pion][j].imag,
            props[prop_pion_pion2][j].real, props[prop_pion_pion2][j].imag,
            props[prop_pion_0mp][j].real, props[prop_pion_0mp][j].imag );
      }
      printf("ENDPROP\n");
} /*PION_SOURCE*/
  
if(  oper_PION2_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: PION2\n");
      printf("SINKS: PION PION2 0MP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_pion2_pion][j].real, props[prop_pion2_pion][j].imag,
            props[prop_pion2_pion2][j].real, props[prop_pion2_pion2][j].imag,
            props[prop_pion2_0mp][j].real, props[prop_pion2_0mp][j].imag );
      }
      printf("ENDPROP\n");
} /*PION2_SOURCE*/
  
if(  oper_ZEROMP_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0MP\n");
      printf("SINKS: PION PION2 0MP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_0mp_pion][j].real, props[prop_0mp_pion][j].imag,
            props[prop_0mp_pion2][j].real, props[prop_0mp_pion2][j].imag,
            props[prop_0mp_0mp][j].real, props[prop_0mp_0mp][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROMP_SOURCE*/
  
if(  oper_RHO_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: RHO\n");
      printf("SINKS: RHO RHO2 1MM\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_rho_rho][j].real, props[prop_rho_rho][j].imag,
            props[prop_rho_rho2][j].real, props[prop_rho_rho2][j].imag,
            props[prop_rho_1mm][j].real, props[prop_rho_1mm][j].imag );
      }
      printf("ENDPROP\n");
} /*RHO_SOURCE*/
  
if(  oper_RHO2_SOURCE  == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: RHO2\n");
      printf("SINKS: RHO RHO2 1MM\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_rho2_rho][j].real, props[prop_rho2_rho][j].imag,
            props[prop_rho2_rho2][j].real, props[prop_rho2_rho2][j].imag,
            props[prop_rho2_1mm][j].real, props[prop_rho2_1mm][j].imag );
      }
      printf("ENDPROP\n");
} /*RHO2_SOURCE*/
  
if(  oper_ONEMM_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 1MM\n");
      printf("SINKS: RHO RHO2 1MM\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_1mm_rho][j].real, props[prop_1mm_rho][j].imag,
            props[prop_1mm_rho2][j].real, props[prop_1mm_rho2][j].imag,
            props[prop_1mm_1mm][j].real, props[prop_1mm_1mm][j].imag );
      }
      printf("ENDPROP\n");
} /*ONEMM_SOURCE*/

if(  oper_A1_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: A1\n");
      printf("SINKS: A1 A1P 1PP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_a1_a1][j].real, props[prop_a1_a1][j].imag,
            props[prop_a1_a1P][j].real, props[prop_a1_a1P][j].imag,
            props[prop_a1_1pp][j].real, props[prop_a1_1pp][j].imag );
      }
      printf("ENDPROP\n");
} /*A1_SOURCE*/
  
if(  oper_A1P_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: A1P\n");
      printf("SINKS: A1 A1P 1PP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_a1P_a1][j].real, props[prop_a1P_a1][j].imag,
            props[prop_a1P_a1P][j].real, props[prop_a1P_a1P][j].imag,
            props[prop_a1P_1pp][j].real, props[prop_a1P_1pp][j].imag );
      }
      printf("ENDPROP\n");
} /*A1P_SOURCE*/

if(  oper_ONEPP_SOURCE  == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 1PP\n");
      printf("SINKS: A1 A1P 1PP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_1pp_a1][j].real, props[prop_1pp_a1][j].imag,
            props[prop_1pp_a1P][j].real, props[prop_1pp_a1P][j].imag,
            props[prop_1pp_1pp][j].real, props[prop_1pp_1pp][j].imag );
      }
      printf("ENDPROP\n");
} /*ONEPP_SOURCE*/
  
if(  oper_ZEROPM_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0PM\n");
      printf("SINKS: 0PM 0PMP 0PMB\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_0pm_0pm][j].real, props[prop_0pm_0pm][j].imag,
            props[prop_0pm_0pmP][j].real, props[prop_0pm_0pmP][j].imag,
            props[prop_0pm_0pmB][j].real, props[prop_0pm_0pmB][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROPM_SOURCE*/
  
if (  oper_ZEROPMP_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0PMP\n");
      printf("SINKS: 0PM 0PMP 0PMB\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_0pmP_0pm][j].real, props[prop_0pmP_0pm][j].imag,
            props[prop_0pmP_0pmP][j].real, props[prop_0pmP_0pmP][j].imag,
            props[prop_0pmP_0pmB][j].real, props[prop_0pmP_0pmB][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROPMP_SOURCE*/
  
if(  oper_ZEROPMB_SOURCE  == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0PMB\n");
      printf("SINKS: 0PM 0PMP 0PMB\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e %e %e\n",j,
            props[prop_0pmB_0pm][j].real, props[prop_0pmB_0pm][j].imag,
            props[prop_0pmB_0pmP][j].real, props[prop_0pmB_0pmP][j].imag,
            props[prop_0pmB_0pmB][j].real, props[prop_0pmB_0pmB][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROPMB_SOURCE*/

if( oper_ZEROMM_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0MM\n");
      printf("SINKS: 0MM 0MMP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_0mm_0mm][j].real, props[prop_0mm_0mm][j].imag,
            props[prop_0mm_0mmP][j].real, props[prop_0mm_0mmP][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROMM_SOURCE*/
  
if(  oper_ZEROMMP_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: 0MMP\n");
      printf("SINKS: 0MM 0MMP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_0mmP_0mm][j].real, props[prop_0mmP_0mm][j].imag,
            props[prop_0mmP_0mmP][j].real, props[prop_0mmP_0mmP][j].imag );
      }
      printf("ENDPROP\n");
} /*ZEROMMP_SOURCE*/
  
if(  oper_ONEMP_SOURCE  == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: ONEMP\n");
      printf("SINKS: ONEMP ONEMP2\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_1mp_1mp][j].real, props[prop_1mp_1mp][j].imag,
            props[prop_1mp_1mp2][j].real, props[prop_1mp_1mp2][j].imag );
      }
      printf("ENDPROP\n");
} /*ONEMP_SOURCE*/
  
if(  oper_ONEMP2_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: ONEMP2\n");
      printf("SINKS: ONEMP2 ONEMP\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_1mp2_1mp2][j].real, props[prop_1mp2_1mp2][j].imag,
            props[prop_1mp2_1mp][j].real, props[prop_1mp2_1mp][j].imag );
      }
      printf("ENDPROP\n");
} /*ONEMP2_SOURCE*/

if( oper_QQQQ_SOURCE == CALCULATE ) 
{
      printf("STARTPROP\n");
      printf("SOURCE: QQQQ\n");
      printf("SINKS: ONEMP ONEMP2\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_qqqq_1mp][j].real, props[prop_qqqq_1mp][j].imag,
            props[prop_qqqq_1mp2][j].real, props[prop_qqqq_1mp2][j].imag );
      }
      printf("ENDPROP\n");

#ifdef NONSENSE
      printf("STARTPROP\n");
      printf("SOURCE: QQQQ\n");
      printf("SINKS: A1 B1\n");
      for(j=0;j<nt;j++){
          printf("%d %e %e %e %e\n",j,
            props[prop_qqqq_junk][j].real, props[prop_qqqq_junk][j].imag,
            props[prop_qqqq_junk2][j].real, props[prop_qqqq_junk2][j].imag );
      }
      printf("ENDPROP\n");
#endif /*NONSENSE*/
} /*QQQQ_SOURCE*/
  
      fflush(stdout);
    } /* end if(this_node==0) */


    /*** calculate the full meson/baryon spectrum ***/

    light_meson_spectrum( t_source) ;
    light_baryon_spectrum(F_OFFSET(quark_store),t_source) ;

  } /* end loop on t_source */

  /* free arrays */
  free(props[0]); free(props);
  
  return(cgn);
} /* spectrum_hybrids */


/* "Multiply by" the quark-antiquark local pion operator */
void mult_pion( field_offset src, field_offset dest ){
   /* operator is gamma_5, another gamma_5 for antiquark, so nothing to do */
    register int i;
    register site *s;

    FORALLSITES(i,s){
	*(wilson_vector *)F_PT(s,dest) = *(wilson_vector *)F_PT(s,src);
    }
}

/* "Multiply by" the second quark-antiquark local pion operator */
void mult_pion2( field_offset src, field_offset dest ){
   /* operator is gamma_0 gamma_5, another gamma_5 for antiquark, 
	so just multiply by gamma_0 */
    register int i;
    register site *s;
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src),
	    (wilson_vector *)F_PT(s,dest), TUP );
    }
}

/* "Multiply by" the quark-antiquark local rho operator */
void mult_rho( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, another gamma_5 for antiquark */
    register int i;
    register site *s;
    wilson_vector tvec;
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec, pdir );
	mult_by_gamma( &tvec, (wilson_vector *)F_PT(s,dest), GAMMAFIVE );
    }
}

/* "Multiply by" the second quark-antiquark local rho operator */
void mult_rho2( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_0 gamma_pdir, another gamma_5 for antiquark */
    register int i;
    register site *s;
    wilson_vector tvec1,tvec2;
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, TUP );
	mult_by_gamma( &tvec1, &tvec2, pdir );
	mult_by_gamma( &tvec2, (wilson_vector *)F_PT(s,dest), GAMMAFIVE );
    }
}

/* "Multiply by" the quark-antiquark local a1 operator */
void mult_a1( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, gamma_5, another gamma_5 for antiquark */
    register int i;
    register site *s;

    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src),
	    (wilson_vector *)F_PT(s,dest), pdir );
    }
}

/* "Multiply by" the quark-antiquark local b1 operator */
void mult_b1( int pdir,  field_offset src, field_offset dest ){
   /* operator is gamma_pdir, gamma_0, gamma_5, another gamma_5 for antiquark */
    register int i;
    register site *s;
    wilson_vector tvec;
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec, pdir );
	mult_by_gamma( &tvec, (wilson_vector *)F_PT(s,dest), TUP );
    }
}

/* "Multiply by" the exotic 0-+ quark bilinear operator */
void mult_zero_pm_B( field_offset src, field_offset dest ){
   /* operator is gamma_0, another gamma_5 for antiquark */
    register int i;
    register site *s;
    wilson_vector tvec;
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec, TUP );
	mult_by_gamma( &tvec, (wilson_vector *)F_PT(s,dest), GAMMAFIVE );
    }
}

/* "Multiply by" the "S-wave" zero-plus-minus operator 
	quark operator is a1 = gamma_5 gamma_i = 1++,
	gluon operator is B_i = 1+-,  take dot product */
void mult_zero_pm( field_offset src, field_offset dest ){
    register int dir,i;
    register site *s;
    wilson_vector tvec;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir=XUP;dir<=ZUP;dir++){
        /* Multiply by dir+1,dir+2 component of magnetic field,
           multiply by gamma_dir gamma_5(for a1) and gamma_five (for
           antiquark propagator ) */
        mult_by_field_strength( (dir+1)%3, (dir+2)%3, src, F_OFFSET(MP) );
        FORALLSITES(i,s){
            mult_by_gamma( &(s->MP), &tvec, dir );
            add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec,
                (wilson_vector *)F_PT(s,dest) );
        }
    } /* end loop on dir */
}


/* "Multiply by" the zero-minus-minus operator 
	quark operator is a1 = gamma_5 gamma_i = 1++,
	gluon operator is E_i = 1--,  take dot product */
void mult_zero_mm( field_offset src, field_offset dest ){
    register int dir,i;
    register site *s;
    wilson_vector tvec;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir=XUP;dir<=ZUP;dir++){
        /* Multiply by dir component of electric field,
           multiply by gamma_dir gamma_5(for a1) and gamma_five (for
           antiquark propagator ) */
        mult_by_field_strength( dir, TUP, src, F_OFFSET(MP) );
        FORALLSITES(i,s){
            mult_by_gamma( &(s->MP), &tvec, dir );
            add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec,
                (wilson_vector *)F_PT(s,dest) );
        }
    } /* end loop on dir */
}



/* "Multiply by" the one-minus-plus operator 
    quark operator is "rho", gluon operator is magnetic field B.
    Cross product is spin one.
    Z component = "rho_X*F_XZ + rho_Y*F_YZ.
   "pdir" is the polarization direction of the meson */
void mult_one_mp( int pdir, field_offset src, field_offset dest ){
    /* use mp as temporary storage */
    register int dir,i;
    register site *s;
    wilson_vector tvec1,tvec2;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions orthogonal to pdir */
    for(dir=XUP;dir<=ZUP;dir++)if(dir != pdir){
	/* Multiply by dir,pdir component of magnetic field, 
	   multiply by gamma_dir (for rho) and gamma_five (for
	   antiquark propagator ) */
	mult_by_field_strength( dir, pdir, src, F_OFFSET(MP) );
	FORALLSITES(i,s){
	    mult_by_gamma( &(s->MP), &tvec1, dir );
	    mult_by_gamma( &tvec1, &tvec2, GAMMAFIVE );
	    
	    add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec2,
		(wilson_vector *)F_PT(s,dest) );
	}
    } /* end loop on dir */
} /* end mult_one_mp */

/* "Multiply by" another one-minus-plus operator 
    quark operator is zero-plus-minus bilinear
    gluon operator is electric field E.
    Z component = "psibar gamma_0 psi * F_0Z
   "pdir" is the polarization direction of the meson */
void mult_one_mp2( int pdir, field_offset src, field_offset dest ){
    /* use mp as temporary storage */
    register int i;
    register site *s;
    wilson_vector tvec1;

    /* Multiply by pdir component of electric field, 
    multiply by gamma_0 and gamma_five (for antiquark propagator ) */
    mult_by_field_strength( TUP, pdir, src, F_OFFSET(MP) );
    FORALLSITES(i,s){
	mult_by_gamma( &(s->MP), &tvec1, TUP );
	mult_by_gamma( &tvec1, (wilson_vector *)F_PT(s,dest), GAMMAFIVE );
    }
} /* end mult_one_mp2 */

/* "Multiply by" another one-minus-plus operator 
   a 4-quark operator, a_1 plus pion.  All quark lines different
   flavors, so don't worry about multiple contractions.
   "src_t" is the operator time slice
   "pdir" is the polarization direction of the meson */
void mult_qqqq1mp( int src_t, int pdir, field_offset src, field_offset dest,
    field_offset work ){
    register int i;
    register site *s;
    /* use "work" as temporary storage */

    /* Multiply by a_1 operator, extra gamma_5 for direction of
	quark line */
    /* compute propagator, back to same time slice */
    /* multiply by pion operator */
    FORALLSITES(i,s){
	mult_by_gamma( (wilson_vector *)F_PT(s,src),
	    (wilson_vector *)F_PT(s,dest), pdir );
	/* two gamma_5's = 1 */
    }
    mat_invert( dest, work );
    FORALLSITES(i,s){
	if(s->t== src_t ){
	    mult_by_gamma( (wilson_vector *)F_PT(s,work),
		(wilson_vector *)F_PT(s,dest), GAMMAFIVE );
	    /* really three gamma_5's - both propagators have one end here */
	}
	else{
	    sub_wilson_vector( (wilson_vector *)F_PT(s,dest),
		(wilson_vector *)F_PT(s,dest), (wilson_vector *)F_PT(s,dest) );
	}
    }
} /* end mult_qqqq1mp */

/* "Multiply by" a zero-plus-minus P wave operator 
	quark operator is a1 = epsilon_ijk gamma_j deriv_k = 1++,
	deriv is forward_deriv - backward_deriv
	gluon operator is B_i = 1+-,  take dot product 
	see "sources.tex" for details. */
    /* Uses mp, vtmp and sss for temporary storage */
void mult_zero_pm_P( field_offset src, field_offset dest ){
    register int dir1,dir2,i;
    register site *s;
    wilson_vector tvec1,tvec2;
    msg_tag *tag0,*tag1;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }

    /* Loop over spatial directions */
    for(dir1=XUP;dir1<=ZUP;dir1++)for(dir2=XUP;dir2<=ZUP;dir2++)if(dir1!=dir2){

        /* multiply by gamma_dir1 (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(i,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, dir1 );
            mult_by_gamma( &tvec1, &(s->sss), GAMMAFIVE );
        }

	/* parallel transport sss from positive dir2 direction */
	/* parallel transport sss, from negative dir2 direction */
        tag0=start_gather_site( F_OFFSET(sss), sizeof(wilson_vector),
            dir2, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){
	    mult_adj_mat_wilson_vec( &(s->link[dir2]), &(s->sss),  &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[dir2]),
		(wilson_vector *)gen_pt[0][i], &tvec2 );
            sub_wilson_vector( &tvec2, (wilson_vector *)gen_pt[1][i], &(s->MP));
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
        mult_by_field_strength( dir1, dir2, F_OFFSET(MP), F_OFFSET(vtmp) );
	FORALLSITES(i,s){
	    add_wilson_vector( (wilson_vector *)F_PT(s,dest), &(s->vtmp),
		(wilson_vector *)F_PT(s,dest) );
	}

        /* Multiply sss by dir1,dir2 component of magnetic field, */
	/* keep in temp. vector mp */
        mult_by_field_strength( dir1, dir2, F_OFFSET(sss), F_OFFSET(MP) );

	/* parallel transport mp from positive dir2 direction */
	/* parallel transport mp from negative dir2 direction */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            dir2, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){
	    mult_adj_mat_wilson_vec( &(s->link[dir2]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[dir2]),
		(wilson_vector *)gen_pt[0][i], &tvec1 );
            sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][i], &tvec1 );
            add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec1,
                (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);

    } /* end loops on dir1 and dir2 */
}


/* "Multiply by" a zero-minus-minus P wave operator 
	quark operator is a1 = epsilon_ijk gamma_j  deriv_k = 1++,
	gluon operator is E_i = 1--,  take dot product */
void mult_zero_mm_P( field_offset src, field_offset dest ){
    register int in,i,j,k;
    register site *s;
    wilson_vector tvec1,tvec2;
    msg_tag *tag0,*tag1;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(i=XUP;i<=ZUP;i++)for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==i || k== i )continue;
	else if ( j == ((i+1)%3) )sign = 1.0;
	else sign = -1.0;

        /* multiply by gamma_j (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(in,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, j );
            mult_by_gamma( &tvec1, &(s->sss), GAMMAFIVE );
        }

	/* parallel transport sss from positive k direction */
	/* parallel transport sss, from negative k direction */
	/* subtract, and multiply by field strength at site */
        tag0=start_gather_site( F_OFFSET(sss), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(in,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->sss),  &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(in,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][in], &tvec2 );
            sub_wilson_vector( &tvec2, (wilson_vector *)gen_pt[1][in],&(s->MP));
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
        mult_by_field_strength( TUP, i, F_OFFSET(MP), F_OFFSET(vtmp) );
	FORALLSITES(in,s){
	    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &(s->vtmp),
		sign, (wilson_vector *)F_PT(s,dest) );
	}

        /* Multiply sss by 0,i component of magnetic field, */
	/* keep in temp. vector mp */
        mult_by_field_strength( TUP, i, F_OFFSET(sss), F_OFFSET(MP) );

	/* parallel transport mp from positive k direction */
	/* parallel transport mp from negative k direction */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(in,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(in,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][in], &tvec1 );
            sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][in], &tvec1 );
            scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &tvec1, sign,
                (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);
    } /* end loops on directions i,j,k */

}


/* "Multiply by" a a1 P wave operator */
void mult_a1_P( int pdir,  field_offset src, field_offset dest ){
    register int i,j,k;
    register site *s;
    wilson_vector tvec1;
    Real sign;
    msg_tag *tag0,*tag1;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==pdir || k== pdir )continue;
	else if ( j == ((pdir+1)%3) )sign = 1.0;
	else sign = -1.0;

        /* multiply by gamma_j (for a1) and gamma_five (for
           antiquark propagator ) */
        FORALLSITES(i,s){
            mult_by_gamma( (wilson_vector *)F_PT(s,src), &tvec1, j );
            mult_by_gamma( &tvec1, &(s->MP), GAMMAFIVE );
        }
	/* parallel transport mp from forwards and backwards. */
        tag0=start_gather_site( F_OFFSET(MP), sizeof(wilson_vector),
            k, EVENANDODD, gen_pt[0] );
	FORALLSITES(i,s){
	    mult_adj_mat_wilson_vec( &(s->link[k]), &(s->MP), &(s->vtmp) );
	}
        tag1=start_gather_site( F_OFFSET(vtmp), sizeof(wilson_vector),
            OPP_DIR(k), EVENANDODD, gen_pt[1] );
	wait_gather(tag0);
	wait_gather(tag1);
	FORALLSITES(i,s){
	    mult_mat_wilson_vec( &(s->link[k]),
		(wilson_vector *)gen_pt[0][i], &tvec1 );
	    sub_wilson_vector( &tvec1, (wilson_vector *)gen_pt[1][i], &tvec1);
            scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest),
		&tvec1, sign, (wilson_vector *)F_PT(s,dest) );
	}
	cleanup_gather(tag0);
	cleanup_gather(tag1);

    } /* end loops on j and k (directions) */
}

/* "Multiply by" a one-minus-minus hybrid operator 
    quark operator is "pion", gluon operator is magnetic field.
   "pdir" is the polarization direction of the meson */
void mult_one_mm( int pdir, field_offset src, field_offset dest ){
    /* use mp as temporary storage */

    /* multiply by magnetic field.  gamma_5 for antiquark cancels
	gamma_5 in pion operator. */
    mult_by_field_strength( (pdir+1)%3, (pdir+2)%3, src, dest );
} /* end mult_one_mm */


/* "Multiply by" a one-plus-plus hybrid operator 
    quark operator is "rho", gluon operator is electric field.
    1++ = epsilon_ijk \psibar gamma_j \psi F_{0,k}
   "pdir" is the polarization direction of the meson */
void mult_one_pp( int pdir, field_offset src, field_offset dest ){
    /* use mp as temporary storage */
    register int i,j,k;
    register site *s;
    wilson_vector tvec;
    register Real sign;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(j=XUP;j<=ZUP;j++)for(k=XUP;k<=ZUP;k++){
	if( j==k || j==pdir || k== pdir )continue;
	else if ( j == ((pdir+1)%3) )sign = 1.0;
	else sign = -1.0;

        mult_by_field_strength( TUP, k, src, F_OFFSET(MP) );
	FORALLSITES(i,s){
            mult_by_gamma( &(s->MP), &tvec, j );
	    scalar_mult_add_wvec( (wilson_vector *)F_PT(s,dest), &tvec, sign,
		(wilson_vector *)F_PT(s,dest) );
	}

    }/* end loop over j,k (directions) */

    /* multiply everything by gamma_5 for antiquarks */
    FORALLSITES(i,s){
        mult_by_gamma( (wilson_vector *)F_PT(s,dest), &tvec, GAMMAFIVE );
	*(wilson_vector *)F_PT(s,dest) = tvec;
    }
} /* end mult_one_pp */

/* "Multiply by" a zero-minus-plus hybrid operator 
    quark operator is "rho", gluon operator is magnetic field.
    0-+ = epsilon_ijk \psibar gamma_i \psi F_{j,k} */
void mult_zero_mp( field_offset src, field_offset dest ){
    /* use mp as temporary storage */
    register int i,j,k,in;
    register site *s;
    wilson_vector tvec;

    /* set destination to zero */
    FORALLSITES(i,s){
	clear_wvec( (wilson_vector *)F_PT(s,dest) );
    }
    /* loop over directions, sign is sign of epsilon tensor */
    for(i=XUP;i<=ZUP;i++){
	/* antisymmetry of epsilon and F_{jk} means no need to sum all terms */
	j=(i+1)%3; k = (i+2)%3;

        mult_by_field_strength( j, k, src, F_OFFSET(MP) );
	FORALLSITES(in,s){
            mult_by_gamma( &(s->MP), &tvec, i );
	    add_wilson_vector( (wilson_vector *)F_PT(s,dest), &tvec,
		(wilson_vector *)F_PT(s,dest) );
	}

    }/* end loop over i,j,k (directions) */

    /* multiply everything by gamma_5 for antiquarks */
    FORALLSITES(in,s){
        mult_by_gamma( (wilson_vector *)F_PT(s,dest), &tvec, GAMMAFIVE );
	*(wilson_vector *)F_PT(s,dest) = tvec;
    }
} /* end mult_zero_mp */





/* Multiply by the field strength.  Arguments are two indices on
  field strength tensor, source and dest vectors */
void mult_by_field_strength( int dir1, int dir2,
    field_offset src, field_offset dest ){
    wilson_vector tvec;

    register site *s;
    register int i;
    int fs_dir = -1;	/* index in field strength tensor, FS_XY to FS_ZT */
    Real sign;	/* minus one if indices in "wrong" order */

    /* figure out what fs_dir is */
    if(dir1==dir2){
	if(this_node==0)printf("Bad call to mult_by_field_strength\n");
	terminate(0);
    }
    if(dir1>dir2){ 
	sign= -1.0;
	i=dir1; dir1=dir2; dir2=i;	/* switch dirs */
    }
    else sign= 1.0;
    switch( dir1+10*dir2 ){
	case XUP+10*YUP: fs_dir = FS_XY; break;
	case XUP+10*ZUP: fs_dir = FS_XZ; break;
	case XUP+10*TUP: fs_dir = FS_XT; break;
	case YUP+10*ZUP: fs_dir = FS_YZ; break;
	case YUP+10*TUP: fs_dir = FS_YT; break;
	case ZUP+10*TUP: fs_dir = FS_ZT; break;
	default: 
            if(this_node==0)printf("Bad call to mult_by_field_strength\n");
            terminate(0); break;
    }

    /* multiply by matrix */
    FORALLSITES(i,s){
	mult_mat_wilson_vec( &(s->field_strength[fs_dir]), 
	    (wilson_vector *)F_PT(s,src), &tvec );
	scalar_mult_wvec( &tvec, sign, (wilson_vector *)F_PT(s,dest) );
    }
}




/* FOR TESTING: multiply src by matrix and check against dest */
void check_invert_cl( field_offset src, field_offset dest ){
    register int i,j,k,flag;
    register site *s;
    dslash_w_site( src, F_OFFSET(MP), PLUS, EVENANDODD);
    cleanup_dslash_wtemps();
    FORALLSITES(i,s){
	scalar_mult_add_wvec( (wilson_vector *)F_PT(s,src),
	    &(s->MP), -kappa, &(s->MP) );
    }
    FORALLSITES(i,s){
        for(flag=0,j=0;j<4;j++)for(k=0;k<3;k++){
	    if( fabs( ((wilson_vector *)F_PT(s,dest))->d[j].c[k].real
		- s->MP.d[j].c[k].real) > 2e-5 )flag=1;
	    if( fabs( ((wilson_vector *)F_PT(s,dest))->d[j].c[k].imag
		- s->MP.d[j].c[k].imag) > 2e-5 )flag=1;
	    if(flag)printf("%d %d %d  ( %.4e , %.4e )  ( %.4e , %.4e )\n",
	        i,j,k,
		((wilson_vector *)F_PT(s,dest))->d[j].c[k].real,
		((wilson_vector *)F_PT(s,dest))->d[j].c[k].imag,
	        s->MP.d[j].c[k].real,s->MP.d[j].c[k].imag);
	    if(flag)terminate(0);
        }
    }
    g_sync(); if(this_node==0){printf("Inversion checked\n");fflush(stdout);}
}



/*TEMP TEST change gauge, see if things change appropriately*
temp_test(){
int dir,fs_dir;
complex cc;
site *s;
s = &(lattice[1]); 
    for(dir=XUP;dir<=TUP;dir++){
	cc = det_su3( &(s->smearlink[dir]) );
	printf("BEFORE: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,dir,cc.real,cc.imag);
    }
    for(fs_dir=FS_XY;fs_dir<=FS_ZT;fs_dir++){
	cc = det_su3( &(s->field_strength[fs_dir]) );
	printf("BEFORE: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,fs_dir,cc.real,cc.imag);
    }
    gaugefix(XUP,(Real)1.8,500,GAUGE_FIX_TOL);
    smear_links(F_OFFSET(link[0]), F_OFFSET(smearlink[0]));
#ifdef SMEAR
    make_field_strength(F_OFFSET(smearlink[0]), F_OFFSET(field_strength[0]));
#else
    make_field_strength(F_OFFSET(link[0]), F_OFFSET(field_strength[0]));
#endif
    for(dir=XUP;dir<=TUP;dir++){
	cc = det_su3( &(s->smearlink[dir]) );
	printf("AFTER: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,dir,cc.real,cc.imag);
    }
    for(fs_dir=FS_XY;fs_dir<=FS_ZT;fs_dir++){
	cc = det_su3( &(s->field_strength[fs_dir]) );
	printf("AFTER: %d %d %d %d dir= %d  Det = %f  %f\n",
	s->x,s->y,s->z,s->t,fs_dir,cc.real,cc.imag);
    }
}
*END TEMP TEST */



