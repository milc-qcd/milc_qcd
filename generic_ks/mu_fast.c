
/* Calculates the derivatives of the Asqtad quark matrix up to 6th order*/
 
#include "generic_ks_includes.h"    /* definitions files and prototypes */

void initialize (field_offset xxx_off){
  site *st;
  int i,j;
  su3_vector * v;
  FORALLSITES(i,st) {
   for(j=0;j<3;j++){
     v = (su3_vector *)F_PT(st,xxx_off);
     v->c[j].real =0;
     v = (su3_vector *)F_PT(st,xxx_off);
     v->c[j].imag =0;
   }
  }
}

/* take a derivative of n order, result is a vector field */
void dnM_dmun_R(int n, field_offset xxx_off, field_offset out_off){

  site *st;
  int i=0;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;


  m = pow(-1.0, n%2);
  m3 = pow(3.0, n);

  /* Start gathers from positive t-direction */
  tag0 = start_gather_site( xxx_off, sizeof(su3_vector), TUP,
		       EVENANDODD, gen_pt[0] );
  tag1 = start_gather_site( xxx_off, sizeof(su3_vector), T3UP,
		       EVENANDODD, gen_pt[1] );
  FORALLSITES(i,st){
    mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
			  (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) );
    mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
			  (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) );
  }
  
  /* Start gathers from negative t-direction */
  tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(su3_vector),
		       OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
  tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(su3_vector),
		       OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );
  
  /* Wait gathers from positive t-direction and multiply by matrix */
  wait_gather(tag0);
  wait_gather(tag1);
  
  /* Do the Uo^F(x)M^-1R(x+0) and Uo^L(x)M^-1R(x+0) multiplication*/
  FORALLSITES(i,st){
    mult_su3_mat_vec( &(t_fatlink[4*i+TUP]),
		      (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) );
    mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
		      (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) );
  }
  
  /* Wait gathers from negative t-direction */
  wait_gather(tag2);
  wait_gather(tag3);
  
  /* Do the dM_M_inv = Uo^F(x)M^-1R(x+0) - (-1)^m*Uo^F(x-0)M^-1R(x-0) 
     + (3)^m*[ Uo^L(x)M^-1R(x+0) - (-1)^m*Uo^L(x-0)M^-1R(x-0) ] */
  
  FORALLSITES(i,st){
    scalar_mult_sub_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
				m, &(st->tempvec[0]) );
    scalar_mult_sub_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
				m, &(st->templongvec[0]) );
    scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
				m3, (su3_vector *)F_PT(st, out_off) );
  }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
    cleanup_gather(tag3);

}

void dn_dMdu_dmun (int n, field_offset xxx_off, field_offset xxx1_off)
{
  site *st;
  int i;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  double_complex trace;


  trace.real = 0.0;
  trace.imag = 0.0;
  // node0_printf("\nINITIAL values inside tr_dnMdmun_term xxx1=%e, grand=%e\n", lattice[1].xxx1.c[0].imag,lattice[1].g_rand.c[0].imag );

  //Loop taking up to 6 derivatives
    m = pow(-1.0, n%2);
    m3 = (-2.0/u0)*pow(3.0, n);


    //node0_printf("j=%d, n[j] = %d, m = %e, m3 = %e\n", j , n[j], m, m3);
    if(n>0){
    /* Start gathers from positive t-direction */
        tag0 = start_gather_site( xxx_off, sizeof(su3_vector), TUP,
                             EVENANDODD, gen_pt[0] );
        tag1 = start_gather_site( xxx_off, sizeof(su3_vector), T3UP,
                             EVENANDODD, gen_pt[1] );

        FORALLSITES(i,st){
          mult_adj_su3_mat_vec( &(t_dfatlink_du0[4*i+TUP]),
                                (su3_vector *)F_PT(st,xxx_off), &(st->tempvec[TUP]) );
          mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
                                (su3_vector *)F_PT(st,xxx_off), &(st->templongvec[TUP]) );
        }

      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(su3_vector),
                           OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(su3_vector),
                           OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );

      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);
      /* Do the Uo^F(x)M^-1R(x+0) and Uo^L(x)M^-1R(x+0) multiplication*/
      FORALLSITES(i,st){
        mult_su3_mat_vec( &(t_dfatlink_du0[4*i+TUP]),
                          (su3_vector *)gen_pt[0][i], &(st->tempvec[0]) );
        mult_su3_mat_vec( &(t_longlink[4*i+TUP]),
                          (su3_vector *)gen_pt[1][i], &(st->templongvec[0]) );
      }

      /* Wait gathers from negative t-direction */
      wait_gather(tag2);
      wait_gather(tag3);
     /* Do the dM_M_inv = Uo^F(x)M^-1R(x+0) - (-1)^m*Uo^F(x-0)M^-1R(x-0)
         + (3)^m*[ Uo^L(x)M^-1R(x+0) - (-1)^m*Uo^L(x-3 0)M^-1R(x-0) ] */

      FORALLSITES(i,st){
        scalar_mult_sub_su3_vector( &(st->tempvec[0]), (su3_vector *)gen_pt[2][i],
                                    m, &(st->tempvec[0]) );
        scalar_mult_sub_su3_vector( &(st->templongvec[0]), (su3_vector *)gen_pt[3][i],
                                    m, &(st->templongvec[0]) );
        scalar_mult_add_su3_vector( &(st->tempvec[0]), &(st->templongvec[0]),
                                    m3, (su3_vector *)F_PT(st,xxx1_off) );
      }

      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
   } //closes if n>0
   else {
         ddslash_fn_du0_site( xxx_off, xxx1_off, EVENANDODD );
       }
 }













double_complex trace( field_offset g_rand, field_offset xxx_off) {
  site *st;
  int i;
  complex cc;
  double_complex trace;

  trace.real=0.0; trace.imag = 0.0;

  FORALLSITES(i,st){
    cc = su3_dot( (su3_vector *)F_PT(st, g_rand), (su3_vector *)F_PT(st, xxx_off));
    trace.real += cc.real;
    trace.imag += cc.imag;
  }
  g_doublesum(&trace.real);
  g_doublesum(&trace.imag);
 
  trace.real *= (1.0/(double)volume);
  trace.imag *= (1.0/(double)volume);
 
  return(trace);

}



void derivatives (field_offset phi_off, field_offset xxx_off, field_offset xxx1_off, Real mass,
		    int jpbp_reps, int npbp_reps){

  double derivatives[6][2]; 
  double_complex tmp, temp[7][7];
  int i,j;
  site *st;
    
    /* Make random source, and do inversion */
    grsource_imp( phi_off, mass, EVENANDODD );
    initialize(xxx_off);
    mat_invert_uml( F_OFFSET(g_rand), xxx_off, phi_off, mass ); //common starting vector M^-1* R in xxx_off
    tmp = trace(F_OFFSET(g_rand), xxx_off);
    node0_printf ("trM_inv: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
   
    for(i=0;i<7;i++){
      dn_dMdu_dmun (i, xxx_off, F_OFFSET(dMdu_x));
      temp[0][i] = trace(F_OFFSET(g_rand), F_OFFSET(dMdu_x));
      //node0_printf ("TR_d%dMdu%d_M_inv: mass %e,  R: %e  Im: %e ( %d of %d )\n", i,i,mass,
      //              tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
    }

  
   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[0][0] = tmp.real;
    derivatives[0][1] = tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), -1, &(st->deriv[0]));    

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv)); // result in dM_M_inv
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[1][0] = -tmp.real;
    derivatives[1][1] = -tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml( F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), 2.0, &(st->deriv[1]));
 
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[2][0] = 2*tmp.real;
    derivatives[2][1] = 2*tmp.imag;    

   // node0_printf("TR_d3trlnM_dmu3 term 3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
     //          derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), -6, &(st->deriv[2]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));   
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[3][0] = -6*tmp.real;
    derivatives[3][1] = -6*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), 24, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] = 24*tmp.real;
    derivatives[4][1] = 24*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), -120, &(st->deriv[4]));
    
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] = -120*tmp.real;
    derivatives[5][1] = -120*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_su3_vector((su3_vector *)F_PT(st, xxx1_off), 720, &(st->deriv[5]));

      //derivatives starting with d2M/dmu2

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[1][0] += tmp.real;
    derivatives[1][1] += tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[1]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -1.0, &(st->deriv[1]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[2][0] += -3*tmp.real;
    derivatives[2][1] += -3*tmp.imag;
    //node0_printf("d3trlnM_dmu3 term 2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
      //         derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[2]), (su3_vector *)F_PT(st, xxx1_off),
                                                  3, &(st->deriv[2]));
 
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[3][0] += 12*tmp.real;
    derivatives[3][1] += 12*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -12, &(st->deriv[3]));

 

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += -60*tmp.real;
    derivatives[4][1] += -60*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                                  60, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += 360*tmp.real;
    derivatives[5][1] += 360*tmp.imag;          

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -360, &(st->deriv[5]));
    //derivatives starting with d2M/dmu2 M_inv d2M/dmu2

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[3][0] += -3*tmp.real;
    derivatives[3][1] += -3*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                                  6, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += 30*tmp.real;
    derivatives[4][1] += 30*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv 
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -30, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -180*tmp.real;
    derivatives[5][1] += -180*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  180, &(st->deriv[5]));

    //derivatives starting with d2M/dmu2 M_inv d2M/dmu2 M_inv d2M/dmu2

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += 30*tmp.real;
    derivatives[5][1] += 30*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -90, &(st->deriv[5]));


    //derivatives starting with d2M/dmu2 M_inv dM/dmu M_inv d2M/dmu2

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -30, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -90*tmp.real;
    derivatives[5][1] += -90*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  180, &(st->deriv[5]));


    //derivatives starting with d3M/dmu3
     
   dnM_dmun_R(3, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[2][0] += tmp.real;
    derivatives[2][1] += tmp.imag;
       //node0_printf("TR_d3trlnM_dmu3 term 1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
         //      derivatives[2][0], derivatives[2][1], jpbp_reps+1, npbp_reps);
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[2]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -1, &(st->deriv[2]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[3][0] += -4*tmp.real;
    derivatives[3][1] += -4*tmp.imag;
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                                  4, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += 20*tmp.real;
    derivatives[4][1] += 20*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -20, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -120*tmp.real;
    derivatives[5][1] += -120*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  120, &(st->deriv[5]));

    //derivatives starting with d3M/dmu3 M-inv d2M/dmu2
     
   dnM_dmun_R(3, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += -10*tmp.real;
    derivatives[4][1] += -10*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                                  10, &(st->deriv[4]));
  
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += 60*tmp.real;
    derivatives[5][1] += 60*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -60, &(st->deriv[5]));


    //derivative starting with d3M/dmu3 M-inv d3M/dmu3

   dnM_dmun_R(3, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -10*tmp.real;
    derivatives[5][1] += -10*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  20, &(st->deriv[5]));



    //derivative starting with d3M/dmu3 M-inv d1M/dmu1 M-inv d2M/dmu2

   dnM_dmun_R(3, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass );
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += 60*tmp.real;
    derivatives[5][1] += 60*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                                  -60, &(st->deriv[5]));


    //derivative starting with d4M/dmu4 M-inv 

   dnM_dmun_R(4, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[3][0] += tmp.real;
    derivatives[3][1] += tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                             -1, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += -5*tmp.real;
    derivatives[4][1] += -5*tmp.imag;


    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             5, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += 30*tmp.real;
    derivatives[5][1] += 30*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -30, &(st->deriv[5]));

    //derivative starting with d4M/dmu4 M-inv d2M/dmu2 

   dnM_dmun_R(4, xxx_off, F_OFFSET(dM_M_inv)); // dM/dmu * starting vector
    
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -15*tmp.real;
    derivatives[5][1] += -15*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             15, &(st->deriv[5]));

    //derivative starting with d5M/dmu5 

   dnM_dmun_R(5, xxx_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[4][0] += tmp.real;
    derivatives[4][1] += tmp.imag;

    initialize(xxx1_off );
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             -1, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += -6*tmp.real;
    derivatives[5][1] += -6*tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             6, &(st->deriv[5]));


    //derivative starting with d6M/dmu6 
   dnM_dmun_R(6, xxx_off, F_OFFSET(dM_M_inv));
    tmp = trace(F_OFFSET(g_rand), F_OFFSET(dM_M_inv));
    derivatives[5][0] += tmp.real;
    derivatives[5][1] += tmp.imag;

    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -1, &(st->deriv[5]));

// branch (2,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));  
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv    
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[2]), (su3_vector *)F_PT(st, xxx1_off),
                                             3, &(st->deriv[2]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                             -12, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             60, &(st->deriv[4]));

  dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -360, &(st->deriv[5]));

// branch (2,1,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                             -12, &(st->deriv[3]));


   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             60, &(st->deriv[4]));


  dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -360, &(st->deriv[5]));

// branch (3,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[3]), (su3_vector *)F_PT(st, xxx1_off),
                                             4, &(st->deriv[3]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             -20, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             120, &(st->deriv[5]));

//branch (4,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(4, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             5, &(st->deriv[4]));

  dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -30, &(st->deriv[5]));

// branch (3,2)

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             10, &(st->deriv[4]));

  dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -60, &(st->deriv[5]));

// branch (2,2,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             -30, &(st->deriv[4]));
 dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             180, &(st->deriv[5]));

// branch (3,1,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             -20, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             120, &(st->deriv[5]));

//branch (2,1,1,1)

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[4]), (su3_vector *)F_PT(st, xxx1_off),
                                             60, &(st->deriv[4]));

   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -360, &(st->deriv[5]));

//other terms

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx_off, phi_off, mass ); // M^-1 * dM_M_inv

   dnM_dmun_R(5, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             6 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(4, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -30 , &(st->deriv[5]));

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -60 , &(st->deriv[5]));

   dnM_dmun_R(3, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -60 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             120 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             180 , &(st->deriv[5]));

   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             180 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -360 , &(st->deriv[5]));


    initialize(xxx_off);
    mat_invert_uml(  F_OFFSET(g_rand), xxx_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx_off);
    mat_invert_uml( F_OFFSET(dM_M_inv), xxx_off, phi_off, mass ); // M^-1 * dM_M_inv

  dnM_dmun_R(4, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             15 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(3, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             -60 , &(st->deriv[5]));

   dnM_dmun_R(1, xxx_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv
   dnM_dmun_R(2, xxx1_off, F_OFFSET(dM_M_inv));
    initialize(xxx1_off);
    mat_invert_uml(  F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass ); // M^-1 * dM_M_inv   dnM_dmun_R(1, xxx1_off, F_OFFSET(dM_M_inv));
    FORALLSITES(i,st) scalar_mult_add_su3_vector(&(st->deriv[5]), (su3_vector *)F_PT(st, xxx1_off),
                                             180 , &(st->deriv[5]));



    for(i=0;i<6;i++)
      node0_printf ("d%d_trlnM_dmu%d: mass %e,  R: %e  Im: %e ( %d of %d )\n", i+1,i+1,mass,
		    derivatives[i][0], derivatives[i][1], jpbp_reps+1, npbp_reps);

   for(i=0;i<6;i++){
      tmp = trace(F_OFFSET(g_rand), F_OFFSET(deriv[i]));
      node0_printf ("d%d_trMi_dmu%d: mass %e,  R: %e  Im: %e ( %d of %d )\n", i+1,i+1,mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);
   } 
   for(i=0;i<6;i++){
    for(j=0;j<6-i;j++){
      dn_dMdu_dmun (j, F_OFFSET(deriv[i]), F_OFFSET(dMdu_x));
      temp[i+1][j] = trace(F_OFFSET(g_rand), F_OFFSET(dMdu_x));
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

}


void Deriv_O6(field_offset phi_off, field_offset xxx_off, field_offset xxx1_off, Real mass ){

#ifndef FN              /* FN is assumed for quark number susc. */
  node0_printf("Problem with FN definition\n");
  terminate(1);
#endif
#ifndef NPBP_REPS       /* Need multiple repetitions for susceptibilities! */
  node0_printf("Problem with NPBP_REP definition\n");
  terminate(1);
#endif

 
  int npbp_reps = npbp_reps_in;  /* Number of repetitions of stochastic estimate */

  int jpbp_reps;

#ifdef FN
  if( !(valid_fn_links==1))  load_fn_links();
#endif
  
  for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++)
    derivatives (phi_off, xxx_off, xxx1_off, mass,
		   jpbp_reps, npbp_reps);
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


