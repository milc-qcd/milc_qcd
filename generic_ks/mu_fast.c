/******* mu_fast.c - conjugate gradient for conventional fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions */

/* 2007 Ludmila Levkova */
/* 10/11 upgrade to version 7.7.2 */


/* Calculates the derivatives of the Asqtad quark matrix up to 6th order*/
 
#include "generic_ks_includes.h"    /* definitions files and prototypes */

#ifndef FN
BOMB the compilation.  Works only for FN actions
#endif

static su3_vector *tfat, *tfat0, *tlong, *tlong0;
static su3_vector *deriv[6];

/* take a derivative of M of n order, result is a vector field in out_off */
static void dnM_dmun_R(int n, su3_vector *xxx_off, su3_vector *out_off,
		imp_ferm_links_t *fn)
{

  site *st;
  int i=0;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  su3_matrix *t_fatlink = get_fatlinks(fn);
  su3_matrix *t_longlink = get_lnglinks(fn);


  m = pow(-1.0, n%2);
  m3 = pow(3.0, n);

  /* Start gathers from positive t-direction */
  tag0 = start_gather_field( xxx_off, sizeof(su3_vector), TUP,
		       EVENANDODD, gen_pt[0] );
  tag1 = start_gather_field( xxx_off, sizeof(su3_vector), T3UP,
		       EVENANDODD, gen_pt[1] );
  FORALLSITES(i,st){
    mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
			  xxx_off+i, tfat+i );
    mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
			  xxx_off+i, tlong+i );
  }
  
  /* Start gathers from negative t-direction */
  tag2 = start_gather_field( tfat, sizeof(su3_vector),
		       OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
  tag3 = start_gather_field( tlong, sizeof(su3_vector),
		       OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );
  
  /* Wait gathers from positive t-direction and multiply by matrix */
  wait_gather(tag0);
  wait_gather(tag1);
  
  /* Do the Uo^F(x)M^-1R(x+0) and Uo^L(x)M^-1R(x+0) multiplication*/
  FORALLSITES(i,st){
    mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
		      (su3_vector *)gen_pt[0][i], tfat0+i );
    mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
		      (su3_vector *)gen_pt[1][i], tlong0+i );
  }
  
  /* Wait gathers from negative t-direction */
  wait_gather(tag2);
  wait_gather(tag3);
  
  /* Do the dM_M_inv = Uo^F(x)M^-1R(x+0) - (-1)^m*Uo^F(x-0)M^-1R(x-0) 
     + (3)^m*[ Uo^L(x)M^-1R(x+0) - (-1)^m*Uo^L(x-0)M^-1R(x-0) ] */
  
  FORALLSITES(i,st){
    scalar_mult_sub_su3_vector( tfat0+i, (su3_vector *)gen_pt[2][i],
				m, tfat0+i );
    scalar_mult_sub_su3_vector( tlong0+i, (su3_vector *)gen_pt[3][i],
				m, tlong0+i );
    scalar_mult_add_su3_vector( tfat0+i, tlong0+i,
				m3, out_off+i );
  }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
    cleanup_gather(tag3);

}

// the pointer fn is just a dummy, it is not actually used in this function
#ifdef DM_DU0
static void dn_dMdu_dmun (int n, su3_vector *xxx_off, su3_vector *xxx1_off,
		   imp_ferm_links_t *fn, imp_ferm_links_t *fn_dmdu0)
{
  site *st;
  int i;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  su3_matrix *t_dfatlink_du0 = get_fatlinks(fn_dmdu0);
  su3_matrix *t_longlink = get_lnglinks(fn_dmdu0);


  // node0_printf("\nINITIAL values inside tr_dnMdmun_term xxx1=%e, grand=%e\n", lattice[1].xxx1.c[0].imag,lattice[1].g_rand.c[0].imag );

  //Loop taking up to 6 derivatives
    m = pow(-1.0, n%2);
    m3 = pow(3.0, n);


    //node0_printf("j=%d, n[j] = %d, m = %e, m3 = %e\n", j , n[j], m, m3);
    if(n>0){
    /* Start gathers from positive t-direction */
        tag0 = start_gather_field( xxx_off, sizeof(su3_vector), TUP,
                             EVENANDODD, gen_pt[0] );
        tag1 = start_gather_field( xxx_off, sizeof(su3_vector), T3UP,
                             EVENANDODD, gen_pt[1] );

        FORALLSITES(i,st){
          mult_adj_su3_mat_vec( &(t_dfatlink_du0[4*i+TUP]),
                                xxx_off+i, tfat+i );
          mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
                                xxx_off+i, tlong+i );
        }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_field( tfat, sizeof(su3_vector),
                           OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_field( tlong, sizeof(su3_vector),
                           OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);
      /* Do the Uo^F(x)M^-1R(x+0) and Uo^L(x)M^-1R(x+0) multiplication*/
      FORALLSITES(i,st){
        mult_su3_mat_vec( &(t_dfatlink_du0[4*i+TUP]),
                          (su3_vector *)gen_pt[0][i], tfat0+i );
        mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
                          (su3_vector *)gen_pt[1][i], tlong0+i );
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);
     /* Do the dM_M_inv = Uo^F(x)M^-1R(x+0) - (-1)^m*Uo^F(x-0)M^-1R(x-0)
         + (3)^m*[ Uo^L(x)M^-1R(x+0) - (-1)^m*Uo^L(x-3 0)M^-1R(x-0) ] */

      FORALLSITES(i,st){
        scalar_mult_sub_su3_vector( tfat0+i, (su3_vector *)gen_pt[2][i],
                                    m, tfat0+i );
        scalar_mult_sub_su3_vector( tlong0+i, (su3_vector *)gen_pt[3][i],
                                    m, tlong0+i );
        scalar_mult_add_su3_vector( tfat0+i, tlong0+i,
                                    m3, xxx1_off+i );
      }

      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
   } //closes if n>0
   else {
          dslash_fn_field(xxx_off, xxx1_off, EVENANDODD,fn_dmdu0);
      }
 }
#endif

#ifdef HISQ
 
static void dn_dMde_dmun (int n, su3_vector *xxx_off, su3_vector *xxx1_off,
		   imp_ferm_links_t *fn_deps)
{
  site *st;
  int i;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  su3_matrix *t_longlink = get_lnglinks(fn_deps);
  su3_matrix *t_fatlink = get_fatlinks(fn_deps);

  if(n>0){
  //Loop taking up to 6 derivatives
    m = pow(-1.0, n%2);
    m3 = pow(3.0, n);

        /* Start gathers from positive t-direction */
        tag0 = start_gather_field( xxx_off, sizeof(su3_vector), TUP,
                             EVENANDODD, gen_pt[0] );

        tag1 = start_gather_field( xxx_off, sizeof(su3_vector), T3UP,
                             EVENANDODD, gen_pt[1] );

        FORALLSITES(i,st){
          mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
                                xxx_off+i, tfat+i );
          mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
                                xxx_off+i, tlong+i );
        }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_field( tfat, sizeof(su3_vector),
                           OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_field( tlong, sizeof(su3_vector),
                           OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );
      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      /* Do the  Uo^L(x)M^-1R(x+0) multiplication*/
      FORALLSITES(i,st){
        mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
                          (su3_vector *)gen_pt[0][i], tfat0+i );
        mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
                          (su3_vector *)gen_pt[1][i], tlong0+i );
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);
     /* Do the dM_M_inv = 
         + (3)^n*[ Uo^L(x)M^-1R(x+0) - (-1)^m*Uo^L(x-3 0)M^-1R(x-3 0) ] */

      FORALLSITES(i,st){
        scalar_mult_sub_su3_vector( tfat0+i, (su3_vector *)gen_pt[2][i],
                                    m, tfat0+i );
        scalar_mult_sub_su3_vector( tlong0+i, (su3_vector *)gen_pt[3][i],
                                    m, tlong0+i );
        scalar_mult_add_su3_vector( tfat0+i, tlong0+i,
                                    m3, xxx1_off+i );

      }

      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
    }
   else { //n==0
    dslash_fn_field( xxx_off, xxx1_off, EVENANDODD, fn_deps );
   }

}


#endif






static double_complex trace( su3_vector *g_rand, su3_vector *xxx_off) {
  site *st;
  int i;
  complex cc;
  double_complex trace;

  trace.real=0.0; trace.imag = 0.0;

  FORALLSITES(i,st){
    cc = su3_dot( g_rand+i, xxx_off+i);
    trace.real += cc.real;
    trace.imag += cc.imag;
  }
  g_doublesum(&trace.real);
  g_doublesum(&trace.imag);
 
  trace.real *= (1.0/(double)volume);
  trace.imag *= (1.0/(double)volume);
 
  return(trace);

}



static void derivatives (su3_vector *phi_off, su3_vector *xxx_off, 
			 su3_vector *xxx1_off, quark_invert_control *qic,
			 Real mass,
			 int jpbp_reps, int npbp_reps,
			 imp_ferm_links_t *fn, imp_ferm_links_t *fn_dmdu0, imp_ferm_links_t *fn_deps, Real eps)
{

  double derivatives[6][2]; 
  double_complex tmp, temp[7][7];
  int i,j;
  site *st;
  su3_vector *g_rand = create_v_field();
  su3_vector *dMdu_x = create_v_field();
  su3_vector *dM_M_inv = create_v_field();
    
  /* Make random source, and do inversion */
  grsource_plain_field( g_rand, EVENANDODD );
    clear_v_field(xxx_off);
    //common starting vector M^-1* R in xxx_off
    mat_invert_uml_field( g_rand, xxx_off, qic, mass, 
		    fn ); 
    tmp = trace(g_rand, xxx_off);
    node0_printf ("trM_inv: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
#ifdef DM_DU0  
    for(i=0;i<7;i++){
      dn_dMdu_dmun (i, xxx_off, dMdu_x, fn, fn_dmdu0);
      temp[0][i] = trace(g_rand, dMdu_x);
      //node0_printf ("TR_d%dMdu%d_M_inv: mass %e,  R: %e  Im: %e ( %d of %d )\n", i,i,mass,
      //              tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
    }
#endif
 
#ifdef HISQ
   if(eps!=0){  
    for(i=0;i<7;i++){
      dn_dMde_dmun (i, xxx_off, dMdu_x, fn_deps);
      temp[0][i] = trace(g_rand, dMdu_x);
      //node0_printf ("TR_d%dMdu%d_M_inv: mass %e,  R: %e  Im: %e ( %d of %d )\n", i,i,mass,
      //              tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
    }
   }
#endif

  
    dnM_dmun_R(1, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    tmp = trace(g_rand, dM_M_inv);
    derivatives[0][0] = tmp.real;
    derivatives[0][1] = tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass, 
		     fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, -1, deriv[0]+i);    

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn); // result in dM_M_inv
    tmp = trace(g_rand, dM_M_inv);
    derivatives[1][0] = -tmp.real;
    derivatives[1][1] = -tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field( dM_M_inv, xxx1_off, qic, mass,
		    fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, 2.0, deriv[1]+i);
 
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[2][0] = 2*tmp.real;
    derivatives[2][1] = 2*tmp.imag;    

   // node0_printf("TR_d3trlnM_dmu3 term 3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
     //          derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, -6, deriv[2]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);   
    tmp = trace(g_rand, dM_M_inv);
    derivatives[3][0] = -6*tmp.real;
    derivatives[3][1] = -6*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, 24, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] = 24*tmp.real;
    derivatives[4][1] = 24*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass, 
		     fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, -120, deriv[4]+i);
    
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] = -120*tmp.real;
    derivatives[5][1] = -120*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_su3_vector( xxx1_off+i, 720, deriv[5]+i);

      //derivatives starting with d2M/dmu2

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    tmp = trace(g_rand, dM_M_inv);
    derivatives[1][0] += tmp.real;
    derivatives[1][1] += tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[1]+i,  xxx1_off+i,
                                                  -1.0, deriv[1]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[2][0] += -3*tmp.real;
    derivatives[2][1] += -3*tmp.imag;
    //node0_printf("d3trlnM_dmu3 term 2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
      //         derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[2]+i,  xxx1_off+i,
                                                  3, deriv[2]+i);
 
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[3][0] += 12*tmp.real;
    derivatives[3][1] += 12*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass, 
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                                  -12, deriv[3]+i);

 

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += -60*tmp.real;
    derivatives[4][1] += -60*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                                  60, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += 360*tmp.real;
    derivatives[5][1] += 360*tmp.imag;          

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  -360, deriv[5]+i);
    //derivatives starting with d2M/dmu2 M_inv d2M/dmu2

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[3][0] += -3*tmp.real;
    derivatives[3][1] += -3*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                                  6, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += 30*tmp.real;
    derivatives[4][1] += 30*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv 
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                                  -30, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -180*tmp.real;
    derivatives[5][1] += -180*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  180, deriv[5]+i);

    //derivatives starting with d2M/dmu2 M_inv d2M/dmu2 M_inv d2M/dmu2

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += 30*tmp.real;
    derivatives[5][1] += 30*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  -90, deriv[5]+i);


    //derivatives starting with d2M/dmu2 M_inv dM/dmu M_inv d2M/dmu2

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                                  -30, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -90*tmp.real;
    derivatives[5][1] += -90*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn );
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  180, deriv[5]+i);


    //derivatives starting with d3M/dmu3
     
    dnM_dmun_R(3, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    tmp = trace(g_rand, dM_M_inv);
    derivatives[2][0] += tmp.real;
    derivatives[2][1] += tmp.imag;
       //node0_printf("TR_d3trlnM_dmu3 term 1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
         //      derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[2]+i,  xxx1_off+i,
                                                  -1, deriv[2]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[3][0] += -4*tmp.real;
    derivatives[3][1] += -4*tmp.imag;
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                                  4, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += 20*tmp.real;
    derivatives[4][1] += 20*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                                  -20, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -120*tmp.real;
    derivatives[5][1] += -120*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  120, deriv[5]+i);

    //derivatives starting with d3M/dmu3 M-inv d2M/dmu2
     
    dnM_dmun_R(3, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += -10*tmp.real;
    derivatives[4][1] += -10*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                                  10, deriv[4]+i);
  
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += 60*tmp.real;
    derivatives[5][1] += 60*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  -60, deriv[5]+i);


    //derivative starting with d3M/dmu3 M-inv d3M/dmu3

    dnM_dmun_R(3, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -10*tmp.real;
    derivatives[5][1] += -10*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  20, deriv[5]+i);



    //derivative starting with d3M/dmu3 M-inv d1M/dmu1 M-inv d2M/dmu2

    dnM_dmun_R(3, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    
    clear_v_field(xxx1_off);
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += 60*tmp.real;
    derivatives[5][1] += 60*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                                  -60, deriv[5]+i);


    //derivative starting with d4M/dmu4 M-inv 

    dnM_dmun_R(4, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    tmp = trace(g_rand, dM_M_inv);
    derivatives[3][0] += tmp.real;
    derivatives[3][1] += tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                             -1, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += -5*tmp.real;
    derivatives[4][1] += -5*tmp.imag;


    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             5, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += 30*tmp.real;
    derivatives[5][1] += 30*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -30, deriv[5]+i);

    //derivative starting with d4M/dmu4 M-inv d2M/dmu2 

    dnM_dmun_R(4, xxx_off, dM_M_inv, fn); // dM/dmu * starting vector
    
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -15*tmp.real;
    derivatives[5][1] += -15*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             15, deriv[5]+i);

    //derivative starting with d5M/dmu5 

    dnM_dmun_R(5, xxx_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[4][0] += tmp.real;
    derivatives[4][1] += tmp.imag;

    clear_v_field(xxx1_off );
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             -1, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += -6*tmp.real;
    derivatives[5][1] += -6*tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             6, deriv[5]+i);


    //derivative starting with d6M/dmu6 
    dnM_dmun_R(6, xxx_off, dM_M_inv, fn);
    tmp = trace(g_rand, dM_M_inv);
    derivatives[5][0] += tmp.real;
    derivatives[5][1] += tmp.imag;

    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -1, deriv[5]+i);

// branch (2,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);  
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv    
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[2]+i,  xxx1_off+i,
                                             3, deriv[2]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                             -12, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             60, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -360, deriv[5]+i);

// branch (2,1,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
   FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                             -12, deriv[3]+i);


   dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             60, deriv[4]+i);


    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -360, deriv[5]+i);

// branch (3,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[3]+i,  xxx1_off+i,
                                             4, deriv[3]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             -20, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             120, deriv[5]+i);

//branch (4,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(4, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             5, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -30, deriv[5]+i);

// branch (3,2)

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             10, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -60, deriv[5]+i);

// branch (2,2,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             -30, deriv[4]+i);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             180, deriv[5]+i);

// branch (3,1,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             -20, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             120, deriv[5]+i);

//branch (2,1,1,1)

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[4]+i,  xxx1_off+i,
                                             60, deriv[4]+i);

    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -360, deriv[5]+i);

//other terms

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx_off, qic, mass,
		     fn);

    dnM_dmun_R(5, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             6 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(4, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -30 , deriv[5]+i);

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -60 , deriv[5]+i);

    dnM_dmun_R(3, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -60 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             120 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             180 , deriv[5]+i);

    dnM_dmun_R(2, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             180 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -360 , deriv[5]+i);


    clear_v_field(xxx_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  g_rand, xxx_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field( dM_M_inv, xxx_off, qic, mass,
		    fn);

    dnM_dmun_R(4, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             15 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(3, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             -60 , deriv[5]+i);

    dnM_dmun_R(1, xxx_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    dnM_dmun_R(2, xxx1_off, dM_M_inv, fn);
    clear_v_field(xxx1_off);
    // M^-1 * dM_M_inv   dnM_dmun_R(1, xxx1_off, dM_M_inv, fn);
    mat_invert_uml_field(  dM_M_inv, xxx1_off, qic, mass,
		     fn);
    FORALLSITES(i,st) scalar_mult_add_su3_vector(deriv[5]+i,  xxx1_off+i,
                                             180 , deriv[5]+i);

// print the results

    for(i=0;i<6;i++)
      node0_printf ("d%d_trlnM_dmu%d: mass %e,  R: %e  Im: %e ( %d of %d )\n", i+1,i+1,mass,
		    derivatives[i][0], derivatives[i][1], jpbp_reps+1, npbp_reps);

   for(i=0;i<6;i++){
      tmp = trace(g_rand, deriv[i]);
      node0_printf ("d%d_trMi_dmu%d: mass %e,  R: %e  Im: %e ( %d of %d )\n", i+1,i+1,mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
   } 

#ifdef DM_DU0
   for(i=0;i<6;i++){
    for(j=0;j<6-i;j++){
      dn_dMdu_dmun (j, deriv[i], dMdu_x, fn,
		    fn_dmdu0);
      temp[i+1][j] = trace(g_rand, dMdu_x);
      //node0_printf ("d%dMdu%d_d%dM_inv_dmu%d: mass %e,  R: %e  Im: %e ( %d of %d )\n",j,j, i+1,i+1,mass,
                    //tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
    }
   }


  node0_printf ("d0_trdMduMi_dmu0: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    temp[0][0].real, temp[0][0].imag, jpbp_reps+1, npbp_reps); 
  tmp.real = temp[1][0].real + temp[0][1].real; 
  tmp.imag = temp[1][0].imag + temp[0][1].imag;
  node0_printf ("d1_trdMduMi_dmu1: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[2][0].real + 2*temp[1][1].real + temp[0][2].real; 
 tmp.imag = temp[2][0].imag + 2*temp[1][1].imag + temp[0][2].imag;
 node0_printf ("d2_trdMduMi_dmu2: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[3][0].real + 3*temp[2][1].real + 3*temp[1][2].real +temp[0][3].real;
 tmp.imag = temp[3][0].imag + 3*temp[2][1].imag + 3*temp[1][2].imag +temp[0][3].imag;
 node0_printf ("d3_trdMduMi_dmu3: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);


 tmp.real = temp[4][0].real + 4*temp[3][1].real + 6*temp[2][2].real + 4*temp[1][3].real + temp[0][4].real;
 tmp.imag = temp[4][0].imag + 4*temp[3][1].imag + 6*temp[2][2].imag + 4*temp[1][3].imag + temp[0][4].imag;
 node0_printf ("d4_trdMduMi_dmu4: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[5][0].real + 5*temp[4][1].real + 10*temp[3][2].real + 5*temp[1][4].real + 10*temp[2][3].real 
           +temp[0][5].real;
 tmp.imag = temp[5][0].imag + 5*temp[4][1].imag + 10*temp[3][2].imag + 5*temp[1][4].imag + 10*temp[2][3].imag 
           +temp[0][5].imag;
 node0_printf ("d5_trdMduMi_dmu5: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[6][0].real + 6*temp[5][1].real + 6*temp[1][5].real + 15*temp[2][4].real + 15*temp[4][2].real 
           +20*temp[3][3].real + temp[0][6].real;
 tmp.imag = temp[6][0].imag + 6*temp[5][1].imag + 6*temp[1][5].imag + 15*temp[2][4].imag + 15*temp[4][2].imag 
           +20*temp[3][3].imag + temp[0][6].imag;
 node0_printf ("d6_trdMduMi_dmu6: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

#endif

#ifdef HISQ
 if(eps != 0){

   for(i=0;i<6;i++){
    for(j=0;j<6-i;j++){
      dn_dMde_dmun (j, deriv[i], dMdu_x, fn_deps);
      temp[i+1][j] = trace(g_rand, dMdu_x);
    }
   }


  node0_printf ("d0_trdMdeMi_dmu0: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    temp[0][0].real, temp[0][0].imag, jpbp_reps+1, npbp_reps); 
  tmp.real = temp[1][0].real + temp[0][1].real; 
  tmp.imag = temp[1][0].imag + temp[0][1].imag;
  node0_printf ("d1_trdMdeMi_dmu1: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[2][0].real + 2*temp[1][1].real + temp[0][2].real; 
 tmp.imag = temp[2][0].imag + 2*temp[1][1].imag + temp[0][2].imag;
 node0_printf ("d2_trdMdeMi_dmu2: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[3][0].real + 3*temp[2][1].real + 3*temp[1][2].real +temp[0][3].real;
 tmp.imag = temp[3][0].imag + 3*temp[2][1].imag + 3*temp[1][2].imag +temp[0][3].imag;
 node0_printf ("d3_trdMdeMi_dmu3: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);


 tmp.real = temp[4][0].real + 4*temp[3][1].real + 6*temp[2][2].real + 4*temp[1][3].real + temp[0][4].real;
 tmp.imag = temp[4][0].imag + 4*temp[3][1].imag + 6*temp[2][2].imag + 4*temp[1][3].imag + temp[0][4].imag;
 node0_printf ("d4_trdMdeMi_dmu4: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[5][0].real + 5*temp[4][1].real + 10*temp[3][2].real + 5*temp[1][4].real + 10*temp[2][3].real 
           +temp[0][5].real;
 tmp.imag = temp[5][0].imag + 5*temp[4][1].imag + 10*temp[3][2].imag + 5*temp[1][4].imag + 10*temp[2][3].imag 
           +temp[0][5].imag;
 node0_printf ("d5_trdMdeMi_dmu5: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[6][0].real + 6*temp[5][1].real + 6*temp[1][5].real + 15*temp[2][4].real + 15*temp[4][2].real 
           +20*temp[3][3].real + temp[0][6].real;
 tmp.imag = temp[6][0].imag + 6*temp[5][1].imag + 6*temp[1][5].imag + 15*temp[2][4].imag + 15*temp[4][2].imag 
           +20*temp[3][3].imag + temp[0][6].imag;
 node0_printf ("d6_trdMdeMi_dmu6: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
}//close if eps!=0
#endif



  destroy_v_field(dM_M_inv);
  destroy_v_field(dMdu_x);
  destroy_v_field(g_rand);
}


void 
Deriv_O6_field( int npbp_reps, quark_invert_control *qic, 
		Real mass, fermion_links_t *fl, int naik_term_epsilon_index, Real eps){

   imp_ferm_links_t *fn = get_fm_links(fl)[naik_term_epsilon_index];

#ifdef DM_DU0
   imp_ferm_links_t *fn_dmdu0 = get_fm_du0_links(fl)[naik_term_epsilon_index];
#else
   imp_ferm_links_t *fn_dmdu0 = NULL;
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
   imp_ferm_links_t *fn_deps = get_fn_deps_links(fl);
#else 
   imp_ferm_links_t *fn_deps = NULL;
#endif


  su3_vector *phi_off, *xxx_off, *xxx1_off;
  phi_off = create_v_field();
  xxx_off = create_v_field();
  xxx1_off = create_v_field();
  tfat = create_v_field();
  tlong = create_v_field();
  tfat0 = create_v_field();
  tlong0 = create_v_field();
  int i;

#ifndef FN              /* FN is assumed for quark number susc. */
  node0_printf("Problem with FN definition\n");
  terminate(1);
#endif
 
  int jpbp_reps;

  for(i=0;i<6;i++){
    deriv[i] = create_v_field();
  }

  for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++)
    derivatives (phi_off, xxx_off, xxx1_off, qic, mass,
		 jpbp_reps, npbp_reps, fn, fn_dmdu0, fn_deps, eps);
  for(i=0;i<6;i++){
    destroy_v_field(deriv[i]);
  }
  destroy_v_field(tfat);
  destroy_v_field(tlong);
  destroy_v_field(tfat0);
  destroy_v_field(tlong0);
  destroy_v_field(xxx1_off);
  destroy_v_field(xxx_off);
  destroy_v_field(phi_off);

/*fflush(stdout);
for(x=0;x<nx;x++)
for(y=0;y<ny;y++)
for(z=0;z<nz;z++)
for(t=0;t<nt;t++)
{
  i=node_index(x,y,z,t);
  for(j=0;j<3;j++){
    node0_printf("(%d %d %d %d) %e %e\n", lattice[i].x, lattice[i].y, lattice[i].z, lattice[i].t,
    lattice[i].g_rand.c[j].real, lattice[i].g_rand.c[j].imag);
    fflush(stdout);
  }

}*/
}

/* DEPRECATED */
void Deriv_O6(int npbp_reps, int prec, field_offset phi_off, field_offset xxx_off, 
	      field_offset xxx1_off, Real mass, fermion_links_t *fl){

  quark_invert_control qic;
  qic.prec       = prec;
  qic.min        = 0;
  qic.max        = niter;
  qic.nrestart   = nrestart;
  qic.parity     = EVENANDODD;
  qic.start_flag = 0;
  qic.nsrc       = 1;
  qic.resid      = sqrt(rsqprop);
  qic.relresid   = 0;

  Deriv_O6_field( npbp_reps, &qic, mass, fl, 0, 0. );
}

void Deriv_O6_multi(int num_masses, int npbp_reps, quark_invert_control *qic,
		            ks_param *ksp, fermion_links_t *fl){

   int m;
   imp_ferm_links_t **fn = get_fm_links(fl);

#ifdef DM_DU0
   imp_ferm_links_t **Fn_dmdu0 = get_fm_du0_links(fl);
#else
  imp_ferm_links_t **Fn_dmdu0 = NULL; 
#endif

#if FERM_ACTION == HISQ & defined(DM_DEPS)
   imp_ferm_links_t *fn_deps = get_fn_deps_links(fl);
#else 
   imp_ferm_links_t *fn_deps = NULL;
#endif
   imp_ferm_links_t *fn_dmdu0;
  
  su3_vector *phi_off, *xxx_off, *xxx1_off;
  phi_off = create_v_field();
  xxx_off = create_v_field();
  xxx1_off = create_v_field();
  tfat = create_v_field();
  tlong = create_v_field();
  tfat0 = create_v_field();
  tlong0 = create_v_field();
  int i;

#ifndef FN              /* FN is assumed for quark number susc. */
  node0_printf("Problem with FN definition\n");
  terminate(1);
#endif
 
  int jpbp_reps;

  for(i=0;i<6;i++){
    deriv[i] = create_v_field();
  }
  for(m=0;m<num_masses;m++){
    if(Fn_dmdu0==NULL) fn_dmdu0=NULL;
    else fn_dmdu0 = Fn_dmdu0[(ksp+m)->naik_term_epsilon_index]; 
 
    for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){
     
      derivatives (phi_off, xxx_off, xxx1_off, qic+m, (ksp+m)->mass,
	            	 jpbp_reps, npbp_reps, fn[(ksp+m)->naik_term_epsilon_index], 
                     fn_dmdu0, fn_deps, (ksp+m)->naik_term_epsilon);
    }
  }


  for(i=0;i<6;i++){
    destroy_v_field(deriv[i]);
  }
  destroy_v_field(tfat);
  destroy_v_field(tlong);
  destroy_v_field(tfat0);
  destroy_v_field(tlong0);
  destroy_v_field(xxx1_off);
  destroy_v_field(xxx_off);
  destroy_v_field(phi_off);

}


