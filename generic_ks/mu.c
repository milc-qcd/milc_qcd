/******************* mu.c *****************************/
/* MIMD version 7 */
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

void pbp(field_offset xxx_off, Real mass){
  site *st;
  int i;
  double_complex cc;
  complex cc1;

  /*FORALLSITES(i,st){
    for(j=0;j<3;j++)
      node0_printf("grand %e %e; xxx1 %e %e\n", lattice[i].g_rand.c[j].real,
		   lattice[i].g_rand.c[j].imag, lattice[i].xxx1.c[j].real,
		   lattice[i].xxx1.c[j].imag);
    
		   }*/
  cc.real =0.0; cc.imag=0.0;

  FORALLSITES(i,st){
    cc1 = su3_dot( &(st->g_rand), (su3_vector *)F_PT(st,xxx_off) );
    cc.real += cc1.real; cc.imag += cc1.imag;
    //node0_printf("cc %e %e\n", cc.real, cc.imag);
  }

  g_doublesum(&cc.real);
  g_doublesum(&cc.imag);
  node0_printf(" PBP: mass %e, Re %e, IM %e\n", mass, cc.real*(1.0/(double)volume), cc.imag*(1.0/(double)volume));
} 

// result ends up in dM_M_inv
void dn_dMdu_dmun (int n, field_offset xxx_off,
		   ferm_links_t *fn, ferm_links_t *fn_dmdu0 )
{
  site *st;
  int i,j;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  double_complex trace;
  complex cc;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;
  su3_matrix *t_dfatlink_du0;

  t_fatlink = fn->fl.fat;
  t_longlink = fn->fl.lng;
  t_dfatlink_du0 = fn_dmdu0->fat;

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
                                    m3, &(st->dMdu_x) );
      }

      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);
   } //closes if n>0
   else {
         ddslash_fn_du0_site( xxx_off, F_OFFSET(dMdu_x), EVENANDODD,
			      fn_dmdu0);
       }
 }





// source vector in xxx_off, flag =0 no last inversion, flag =1 with last inversion
double_complex tr_dnMdmun_term (int n[6], field_offset g_rand, field_offset phi_off, 
				field_offset xxx_off, field_offset xxx1_off, Real mass, int flag, int u){
  site *st;
  int i,j;
  msg_tag *tag0, *tag1, *tag2, *tag3;
  double m, m3;
  double_complex trace;
  complex cc;
  

  trace.real = 0.0;
  trace.imag = 0.0;  
  // node0_printf("\nINITIAL values inside tr_dnMdmun_term xxx1=%e, grand=%e\n", lattice[1].xxx1.c[0].imag,lattice[1].g_rand.c[0].imag );

  //Loop taking up to 6 derivatives

  for(j=0;j<6;j++){
    m = pow(-1.0, n[j]%2);
    m3 = pow(3.0, n[j]);
 
    
    //node0_printf("j=%d, n[j] = %d, m = %e, m3 = %e\n", j , n[j], m, m3);
    if(n[j]>0){
    /* Start gathers from positive t-direction */
      if(j==0){
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

	//node0_printf("fat link %e %e\n", lattice[0].fatlink[3].e[1][0].real, lattice[0].fatlink[3].e[1][0].imag);

      }

      else {
	tag0 = start_gather_site( xxx1_off, sizeof(su3_vector), TUP,
			     EVENANDODD, gen_pt[0] );
	tag1 = start_gather_site( xxx1_off, sizeof(su3_vector), T3UP,
			     EVENANDODD, gen_pt[1] );
	/* do the Uo^F(x)M^-1R(x) and  Uo^L(x)M^-1R(x) multiplication */ 
	
	FORALLSITES(i,st){
	  mult_adj_su3_mat_vec( &(t_fatlink[4*i+TUP]),
				(su3_vector *)F_PT(st,xxx1_off), &(st->tempvec[TUP]) );
	  mult_adj_su3_mat_vec( &(t_longlink[4*i+TUP]),
				(su3_vector *)F_PT(st,xxx1_off), &(st->templongvec[TUP]) );
	}
      }
      
      /* Start gathers from negative t-direction */
      tag2 = start_gather_site( F_OFFSET(tempvec[TUP]), sizeof(su3_vector),
			   OPP_DIR(TUP), EVENANDODD, gen_pt[2] );
      tag3 = start_gather_site( F_OFFSET(templongvec[TUP]), sizeof(su3_vector),
			   OPP_3_DIR(T3UP), EVENANDODD, gen_pt[3] );
      
      /* Wait gathers from positive t-direction and multiply by matrix */
      wait_gather(tag0);
      wait_gather(tag1);

      //tmp =*(su3_vector*)gen_pt[0][0];
      //node0_printf("gen_pt[0] content %e %e\n",tmp.c[0].real, tmp.c[0].imag);
      //node0_printf("xxx1 content %e %e\n", lattice[288].xxx1.c[0].real, lattice[288].xxx1.c[0].imag);      

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
				    m3, &(st->dM_M_inv) );
      }
      //  node0_printf("FINISHED MULTIPLICATION val =%e %e\n", lattice[0].dM_M_inv.c[0].real, lattice[0].dM_M_inv.c[0].imag);      
      if( (j+1 < 6 ) && (n[j+1]>0)){
        initialize(xxx1_off);
	mat_invert_uml( F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass, PRECISION,
			fn);
      }
      else
      {
        if (flag==1) {
         initialize(xxx1_off);
         mat_invert_uml( F_OFFSET(dM_M_inv), xxx1_off, phi_off, mass, PRECISION,
			 fn);
         if (u>-1) dn_dMdu_dmun (u, xxx1_off, fn, fn_dmdu0 );     
        }
      }
      //node0_printf("FINISHED INVERSION val =%e %e \n", lattice[0].xxx2.c[0].real, lattice[0].xxx2.c[0].imag);
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      cleanup_gather(tag2);
      cleanup_gather(tag3);   
      
    } //finish the if(n[j]>0 statement which determines whether the procedure should be done
  } // finish loop over derivatives j in the term

  // in dM_M_inv now we should have dM^i* M_inv *dM^j ....*M_inv R
  FORALLSITES(i,st){
   if(flag==0)
    cc = su3_dot( &(st->g_rand), &(st->dM_M_inv));
   else 
    if (u<0)
     cc = su3_dot( &(st->g_rand), (su3_vector*)F_PT(st, xxx1_off));
    else
     cc = su3_dot( &(st->g_rand), &(st->dMdu_x));
    trace.real += cc.real;
    trace.imag += cc.imag;
  }
  trace.real *= (1.0/(double)volume);
  trace.imag *= (1.0/(double)volume);
  g_doublesum(&trace.real);
  g_doublesum(&trace.imag);
  
  //node0_printf("TR_dtrlnM_dmu: mass %e,  R: %e  Im: %e\n", mass,trace.real, trace.imag);  
  return(trace);
}





void d1trlnM_dmu1   (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{
  int n[6] = {1,0,0,0,0,0};
  double_complex  trace;
  //  node0_printf("INSIDE d1trlnM_dmu1\n\n");
  trace = tr_dnMdmun_term(n, g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);

  node0_printf("d1trlnM_dmu1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace.real, trace.imag, jpbp_reps+1, npbp_reps);  

}

void d2trlnM_dmu2 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{
  int n[2][6] ={{2,0,0,0,0,0},
                {1,1,0,0,0,0}};
  double_complex trace[2]; // first index is the term index, the second is re/im
  double trace_re, trace_im;  
  
  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass,0, -1);
  node0_printf("TR_d2M_dmu2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
         trace[0].real, trace[0].imag, jpbp_reps+1, npbp_reps); 
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  node0_printf("TR_MidM_MidM: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
         trace[1].real, trace[1].imag, jpbp_reps+1, npbp_reps);
  trace_re = trace[0].real - trace[1].real;
  trace_im = trace[0].imag - trace[1].imag;

  node0_printf("d2trlnM_dmu2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace_re, trace_im, jpbp_reps+1, npbp_reps); 

}

void d3trlnM_dmu3 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{

  int n[3][6] ={{3,0,0,0,0,0},
		{2,1,0,0,0,0},
		{1,1,1,0,0,0}};
  double_complex trace[3]; // first index is the term index, the second is re/im
  double trace_re, trace_im;  

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[2] = tr_dnMdmun_term(n[2], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);

  trace_re = trace[0].real - 3*trace[1].real + 2*trace[2].real;
  trace_im = trace[0].imag - 3*trace[1].imag + 2*trace[2].imag;

  node0_printf("d3trlnM_dmu3 term 1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               trace[0].real, trace[0].imag, jpbp_reps+1, npbp_reps);
  node0_printf("d3trlnM_dmu3 term 2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               -3.0*trace[1].real, -3.0*trace[1].imag, jpbp_reps+1, npbp_reps);
  node0_printf("d3trlnM_dmu3 term 3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               2.0*trace[2].real, 2.0*trace[2].imag, jpbp_reps+1, npbp_reps);

  node0_printf("d3trlnM_dmu3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace_re, trace_im, jpbp_reps+1, npbp_reps);

}


void d4trlnM_dmu4 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{
  int n[5][6] ={{4,0,0,0,0,0},
		{3,1,0,0,0,0},
		{2,2,0,0,0,0},
		{2,1,1,0,0,0},
		{1,1,1,1,0,0}};

  double_complex trace[5]; // first index is the term index, the second is re/im
  double trace_re, trace_im;  

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[2] = tr_dnMdmun_term(n[2], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[3] = tr_dnMdmun_term(n[3], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[4] = tr_dnMdmun_term(n[4], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);

  trace_re = trace[0].real - 4*trace[1].real - 3*trace[2].real + 12*trace[3].real - 6*trace[4].real;
  trace_im = trace[0].imag - 4*trace[1].imag - 3*trace[2].imag + 12*trace[3].imag - 6*trace[4].imag;

  node0_printf("d4trlnM_dmu4: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace_re, trace_im, jpbp_reps+1, npbp_reps);


}

void d5trlnM_dmu5 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{
  int n[7][6] ={{5,0,0,0,0,0},
		{4,1,0,0,0,0},
		{3,2,0,0,0,0},
		{3,1,1,0,0,0},
		{2,2,1,0,0,0},
		{2,1,1,1,0,0},
		{1,1,1,1,1,0}};

  double_complex trace[7]; // first index is the term index, the second is re/im
  double trace_re, trace_im;  

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[2] = tr_dnMdmun_term(n[2], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[3] = tr_dnMdmun_term(n[3], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[4] = tr_dnMdmun_term(n[4], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[5] = tr_dnMdmun_term(n[5], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[6] = tr_dnMdmun_term(n[6], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);

  trace_re = trace[0].real - 5*trace[1].real - 10*trace[2].real + 20*trace[3].real + 30*trace[4].real
    - 60*trace[5].real +24*trace[6].real;
  trace_im = trace[0].imag - 5*trace[1].imag - 10*trace[2].imag + 20*trace[3].imag + 30*trace[4].imag
    - 60*trace[5].imag +24*trace[6].imag;

  node0_printf("d5trlnM_dmu5: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace_re, trace_im, jpbp_reps+1, npbp_reps);

}


void d6trlnM_dmu6 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
			     field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{


  int n[13][6] ={{6,0,0,0,0,0},
		 {5,1,0,0,0,0},
		 {4,2,0,0,0,0},
		 {3,3,0,0,0,0},
		 {4,1,1,0,0,0},
		 {3,2,1,0,0,0},
		 {3,1,2,0,0,0},
		 {2,2,2,0,0,0},
		 {3,1,1,1,0,0},
		 {2,2,1,1,0,0},
		 {2,1,2,1,0,0},
		 {2,1,1,1,1,0},
		 {1,1,1,1,1,1}};

  double_complex trace[13]; // first index is the term index, the second is re/im
  double trace_re, trace_im;  

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[2] = tr_dnMdmun_term(n[2], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[3] = tr_dnMdmun_term(n[3], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[4] = tr_dnMdmun_term(n[4], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[5] = tr_dnMdmun_term(n[5], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[6] = tr_dnMdmun_term(n[6], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[7] = tr_dnMdmun_term(n[7], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[8] = tr_dnMdmun_term(n[8], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[9] = tr_dnMdmun_term(n[9], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[10] = tr_dnMdmun_term(n[10], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[11] = tr_dnMdmun_term(n[11], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);
  trace[12] = tr_dnMdmun_term(n[12], g_rand, phi_off, xxx_off, xxx1_off, mass, 0, -1);



  trace_re = trace[0].real - 6*trace[1].real - 15*trace[2].real - 10*trace[3].real + 30*trace[4].real
        + 60*trace[5].real + 60*trace[6].real + 30*trace[7].real - 120*trace[8].real - 180*trace[9].real
        - 90*trace[10].real + 360*trace[11].real - 120*trace[12].real;

  trace_im = trace[0].imag - 6*trace[1].imag - 15*trace[2].imag - 10*trace[3].imag + 30*trace[4].imag
        + 60*trace[5].imag + 60*trace[6].imag + 30*trace[7].imag - 120*trace[8].imag - 180*trace[9].imag
        - 90*trace[10].imag + 360*trace[11].imag - 120*trace[12].imag;

  node0_printf("d6trlnM_dmu6: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
	       trace_re, trace_im, jpbp_reps+1, npbp_reps);


}

//the derivatives of M_inv require M^-1*M^-1 R as initial value of xxx_off
 
double_complex d1trM_inv_dmu1 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{
  int n[6] = {1,0,0,0,0,0};
  double_complex  trace;
  //  node0_printf("INSIDE d1trlnM_dmu1\n\n");
  trace = tr_dnMdmun_term(n, g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  trace.real *=-1;
  trace.imag *=-1;
  node0_printf("d1trM_inv_dmu1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               trace.real, trace.imag, jpbp_reps+1, npbp_reps);
  return(trace);
}

double_complex d2trM_inv_dmu2 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{
  int n[2][6] ={{2,0,0,0,0,0},
                {1,1,0,0,0,0}};
  double_complex trace[2]; // first index is the term index, the second is re/im
  double_complex Trace;

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  node0_printf("TR_d2M_dmu2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
         trace[0].real, trace[0].imag, jpbp_reps+1, npbp_reps);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  node0_printf("TR_MidM_MidM: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
       trace[1].real, trace[1].imag, jpbp_reps+1, npbp_reps);
  Trace.real = -trace[0].real +2*trace[1].real;
  Trace.imag = -trace[0].imag +2*trace[1].imag;

  node0_printf("d2trM_inv_dmu2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               Trace.real, Trace.imag, jpbp_reps+1, npbp_reps);
 return(Trace);
}

double_complex d3trM_inv_dmu3 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{

  int n[4][6] ={{3,0,0,0,0,0},
                {2,1,0,0,0,0},
                {1,2,0,0,0,0},
                {1,1,1,0,0,0}};
  double_complex trace[4]; // first index is the term index, the second is re/im
  double_complex Trace;

  trace[0] = tr_dnMdmun_term(n[0], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  trace[1] = tr_dnMdmun_term(n[1], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  trace[2] = tr_dnMdmun_term(n[2], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  trace[3] = tr_dnMdmun_term(n[3], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);

  Trace.real = -trace[0].real + 3*trace[1].real + 3*trace[2].real -6*trace[3].real;
  Trace.imag = -trace[0].imag + 3*trace[1].imag + 3*trace[2].imag -6*trace[3].imag;

  /*node0_printf("TR_d3trlnM_dmu3 term 1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               trace[0].real, trace[0].imag, jpbp_reps+1, npbp_reps);
  node0_printf("TR_d3trlnM_dmu3 term 2: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               -3.0*trace[1].real, -3.0*trace[1].imag, jpbp_reps+1, npbp_reps);
  node0_printf("TR_d3trlnM_dmu3 term 3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               2.0*trace[2].real, 2.0*trace[2].imag, jpbp_reps+1, npbp_reps);
*/
  node0_printf("d3trM_inv_dmu3: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               Trace.real, Trace.imag, jpbp_reps+1, npbp_reps);
return(Trace);
}

double_complex d4trM_inv_dmu4 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{ int i;
  int n[8][6] ={{4,0,0,0,0,0},
                {3,1,0,0,0,0},
                {1,3,0,0,0,0},
                {2,2,0,0,0,0},
                {2,1,1,0,0,0},
                {1,2,1,0,0,0},
                {1,1,2,0,0,0},
                {1,1,1,1,0,0}};

  double_complex trace[8]; // first index is the term index, the second is re/im
  double_complex Trace;
  for(i=0;i<8;i++)
   trace[i] = tr_dnMdmun_term(n[i], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);

  Trace.real = -trace[0].real + 4*trace[1].real +4*trace[2].real + 6*trace[3].real - 12*trace[4].real
             - 12*trace[5].real - 12*trace[6].real + 24*trace[7].real;
  Trace.imag = -trace[0].imag + 4*trace[1].imag +4*trace[2].imag + 6*trace[3].imag - 12*trace[4].imag
             - 12*trace[5].imag - 12*trace[6].imag + 24*trace[7].imag;

  node0_printf("d4trM_inv_dmu4: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               Trace.real, Trace.imag, jpbp_reps+1, npbp_reps);
 return(Trace);
}

double_complex d5trM_inv_dmu5 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{ int i;
  int n[16][6] ={{5,0,0,0,0,0},
                 {4,1,0,0,0,0},
                 {1,4,0,0,0,0},
                 {3,2,0,0,0,0},
                 {2,3,0,0,0,0},
                 {3,1,1,0,0,0},
                 {1,3,1,0,0,0},
                 {1,1,3,0,0,0},
                 {2,2,1,0,0,0},
                 {2,1,2,0,0,0},
                 {1,2,2,0,0,0},
                 {2,1,1,1,0,0},
                 {1,2,1,1,0,0},
                 {1,1,2,1,0,0},
                 {1,1,1,2,0,0},
                 {1,1,1,1,1,0}};

  double_complex trace[16]; // first index is the term index, the second is re/im
  double_complex Trace;
  for(i=0;i<16;i++)
   trace[i] = tr_dnMdmun_term(n[i], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
 /*for(i=0;i<16;i++)
    node0_printf("TERM1: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               trace[i].real, trace[i].imag, jpbp_reps+1, npbp_reps);
*/
  Trace.real = -trace[0].real +5*trace[1].real + 5*trace[2].real + 10*trace[3].real + 10*trace[4].real
    - 20*trace[5].real -20*trace[6].real -20*trace[7].real
    - 30*trace[8].real -30*trace[9].real -30*trace[10].real 
    + 60*trace[11].real +60*trace[12].real +60*trace[13].real +60*trace[14].real
    -120 * trace[15].real ;
  Trace.imag = -trace[0].imag +5*trace[1].imag + 5*trace[2].imag + 10*trace[3].imag + 10*trace[4].imag
    - 20*trace[5].imag -20*trace[6].imag -20*trace[7].imag
    - 30*trace[8].imag -30*trace[9].imag -30*trace[10].imag
    + 60*trace[11].imag +60*trace[12].imag +60*trace[13].imag +60*trace[14].imag
    -120 * trace[15].imag;

  node0_printf("d5trM_inv_dmu5: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               Trace.real, Trace.imag, jpbp_reps+1, npbp_reps);

return(Trace);
}

 
double_complex d6trM_inv_dmu6 (field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                             field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps, int u)
{
  int i;
  int n[32][6] ={{6,0,0,0,0,0},
                 {5,1,0,0,0,0},{1,5,0,0,0,0},
                 {2,4,0,0,0,0},{4,2,0,0,0,0},
                 {3,3,0,0,0,0},
                 {4,1,1,0,0,0},{1,4,1,0,0,0},{1,1,4,0,0,0},
                 {3,2,1,0,0,0},{2,3,1,0,0,0},{1,2,3,0,0,0},{1,3,2,0,0,0},{2,1,3,0,0,0},{3,1,2,0,0,0},
                 {2,2,2,0,0,0},
                 {3,1,1,1,0,0},{1,3,1,1,0,0},{1,1,3,1,0,0},{1,1,1,3,0,0},
                 {2,2,1,1,0,0},{2,1,2,1,0,0},{1,2,1,2,0,0},{1,1,2,2,0,0},{2,1,1,2,0,0},{1,2,2,1,0,0},
                 {2,1,1,1,1,0},{1,2,1,1,1,0},{1,1,2,1,1,0},{1,1,1,2,1,0},{1,1,1,1,2,0},
                 {1,1,1,1,1,1}};

  double_complex trace[32]; // first index is the term index, the second is re/im
  double_complex Trace;

  for(i=0;i<32;i++) 
     trace[i] = tr_dnMdmun_term(n[i], g_rand, phi_off, xxx_off, xxx1_off, mass, 1, u);
  

  Trace.real = -trace[0].real + 6*trace[1].real +6*trace[2].real +15*trace[3].real + 15*trace[4].real
        +20*trace[5].real -30*trace[6].real - 30*trace[7].real - 30*trace[8].real - 60*trace[9].real
        - 60*trace[10].real - 60*trace[11].real - 60*trace[12].real - 60*trace[13].real - 60*trace[14].real
        -90*trace[15].real +120*(trace[16].real + trace[17].real +trace[18].real  + trace[19].real) 
        +180*(trace[20].real + trace[21].real +trace[22].real  + trace[23].real +trace[24].real +trace[25].real)
        -360*(trace[26].real + trace[27].real +trace[28].real  + trace[29].real +trace[30].real)
        +720*trace[31].real;

  Trace.imag = -trace[0].imag + 6*trace[1].imag +6*trace[2].imag +15*trace[3].imag + 15*trace[4].imag
        +20*trace[5].imag -30*trace[6].imag - 30*trace[7].imag - 30*trace[8].imag - 60*trace[9].imag
        - 60*trace[10].imag - 60*trace[11].imag - 60*trace[12].imag - 60*trace[13].imag - 60*trace[14].imag
        -90*trace[15].imag +120*(trace[16].imag + trace[17].imag +trace[18].imag  + trace[19].imag)
        +180*(trace[20].imag + trace[21].imag +trace[22].imag  + trace[23].imag +trace[24].imag +trace[25].imag)
        -360*(trace[26].imag + trace[27].imag +trace[28].imag  + trace[29].imag +trace[30].imag)
        +720*trace[31].imag; 

  node0_printf("d6trM_inv_dmu6: mass %e,  R: %e  Im: %e ( %d of %d )\n", mass,
               Trace.real, Trace.imag, jpbp_reps+1, npbp_reps);

return(Trace);
}

double_complex trace(int n) {
  site *st;
  int i;
  complex cc;
  double_complex trace;

  trace.real=0.0; trace.imag = 0.0;

  FORALLSITES(i,st){
    cc = su3_dot( &(st->g_rand), &(st->dMdu_x));
    trace.real += cc.real;
    trace.imag += cc.imag;
  }
  g_doublesum(&trace.real);
  g_doublesum(&trace.imag);

  trace.real *= (1.0/(double)volume);
  trace.imag *= (1.0/(double)volume);

  node0_printf("TR_M_inv_d%d_dMdu_dmu%d  %e %e\n",n,n,trace.real,trace.imag);
return(trace);
}

void TR_M_inv_dMdu_deriv(field_offset g_rand, field_offset phi_off, field_offset xxx_off,
                         field_offset xxx1_off, Real mass, int jpbp_reps, int npbp_reps)
{
 //mat_invert_uml( F_OFFSET(g_rand), xxx_off, phi_off, mass, PRECISION,  fn); //xxx_off =M_inv*R
 double_complex tmp, temp[7][7];

 dn_dMdu_dmun (0, xxx_off, fn, fn_dmdu0 );
 temp[0][0]=trace(0);
 temp[1][0]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);
 temp[2][0]=d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);
 temp[3][0]=d3trM_inv_dmu3(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);
 temp[4][0]=d4trM_inv_dmu4(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);
 temp[5][0]=d5trM_inv_dmu5(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);
 temp[6][0]=d6trM_inv_dmu6(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps, 0);


 dn_dMdu_dmun (1, xxx_off, fn, fn_dmdu0 ); // result is in dM_M_inv
 temp[0][1]=trace(1);

 temp[1][1]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,1);
 temp[2][1]=d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,1);
 temp[3][1]=d3trM_inv_dmu3(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,1);
 temp[4][1]=d4trM_inv_dmu4(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,1);
 temp[5][1]=d5trM_inv_dmu5(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,1);


 dn_dMdu_dmun (2, xxx_off, fn, fn_dmdu0 ); // result is in dM_M_inv
 temp[0][2]=trace(2);

 temp[1][2]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,2);
 temp[2][2]=d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,2);
 temp[3][2]=d3trM_inv_dmu3(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,2);
 temp[4][2]=d4trM_inv_dmu4(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,2);


 dn_dMdu_dmun (3, xxx_off, fn, fn_dmdu0 ); //. result is in dM_M_inv
 temp[0][3]=trace(3);

 temp[1][3]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,3);
 temp[2][3]=d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,3);
 temp[3][3]=d3trM_inv_dmu3(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,3);


 dn_dMdu_dmun (4, xxx_off, fn, fn_dmdu0 ); //. result is in dM_M_inv
 temp[0][4]=trace(4);

 temp[1][4]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,4);
 temp[2][4]=d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,4);

 dn_dMdu_dmun (5, xxx_off, fn, fn_dmdu0 ); //
 temp[0][5]=trace(5);

 temp[1][5]=d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off         , xxx1_off, mass, jpbp_reps, npbp_reps,5);


 dn_dMdu_dmun (6, xxx_off, fn, fn_dmdu0 ); 
 temp[0][6]=trace(6);


  node0_printf ("d0_dMduM_i_dmu0: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    temp[0][0].real, temp[0][0].imag, jpbp_reps+1, npbp_reps);
 
  tmp.real = temp[1][0].real + temp[0][1].real;
  node0_printf ("\nRESULT: %e = %e + (%e)\n", tmp.real, temp[1][0].real, temp[0][1]);
  tmp.imag = temp[1][0].imag + temp[0][1].imag;
  node0_printf ("d1_dMduM_i_dmu1: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[2][0].real + 2*temp[1][1].real + temp[0][2].real;
 tmp.imag = temp[2][0].imag + 2*temp[1][1].imag + temp[0][2].imag;
 node0_printf ("d2_dMduM_i_dmu2: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[3][0].real + 3*temp[2][1].real + 3*temp[1][2].real +temp[0][3].real;
 tmp.imag = temp[3][0].imag + 3*temp[2][1].imag + 3*temp[1][2].imag +temp[0][3].imag;
 node0_printf ("d3_dMduM_i_dmu3: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);


 tmp.real = temp[4][0].real + 4*temp[3][1].real + 6*temp[2][2].real + 4*temp[1][3].real + temp[0][4].real;
 tmp.imag = temp[4][0].imag + 4*temp[3][1].imag + 6*temp[2][2].imag + 4*temp[1][3].imag + temp[0][4].imag;
 node0_printf ("d4_dMduM_i_dmu4: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[5][0].real + 5*temp[4][1].real + 10*temp[3][2].real + 5*temp[1][4].real + 10*temp[2][3].real
           +temp[0][5].real;
 tmp.imag = temp[5][0].imag + 5*temp[4][1].imag + 10*temp[3][2].imag + 5*temp[1][4].imag + 10*temp[2][3].imag
           +temp[0][5].imag;
 node0_printf ("d5_dMduM_i_dmu5: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

 tmp.real = temp[6][0].real + 6*temp[5][1].real + 6*temp[1][5].real + 15*temp[2][4].real + 15*temp[4][2].real
           +20*temp[3][3].real + temp[0][6].real;
 tmp.imag = temp[6][0].imag + 6*temp[5][1].imag + 6*temp[1][5].imag + 15*temp[2][4].imag + 15*temp[4][2].imag
           +20*temp[3][3].imag + temp[0][6].imag;
 node0_printf ("d6_dMduM_i_dmu6: mass %e,  R: %.9e  Im: %.9e ( %d of %d )\n",mass,
                    tmp.real, tmp.imag, jpbp_reps+1, npbp_reps);

}


void M_derivatives(field_offset phi_off, field_offset xxx_off, 
		   field_offset xxx1_off, Real mass,
		   ferm_links_t *fn, ferm_links_t *fn_dmdu0)
{
#ifndef FN              /* FN is assumed for quark number susc. */
  node0_printf("Problem with FN definition\n");
  terminate(1);
#endif
#ifndef NPBP_REPS       /* Need multiple repetitions for susceptibilities! */
  node0_printf("Problem with  NPBP_REPS definition\n");
  terminate(1);  
#endif

 
  int npbp_reps = npbp_reps_in;  /* Number of repetitions of stochastic estimate */

  int jpbp_reps;

  //  node0_printf("fatlink %e %e\n", lattice[0].fatlink[3].e[1][0].real, lattice[0].fatlink[3].e[1][0].imag);
  //node0_printf("longlink %e %e\n", lattice[0].longlink[3].e[1][0].real, lattice[0].longlink[3].e[1][0].imag); 


  for(jpbp_reps = 0; jpbp_reps < npbp_reps; jpbp_reps++){
   
    /* Make random source, and do inversion */
    grsource_imp( phi_off, mass, EVENANDODD );

    //node0_printf("Before FIRST INVERSION g_rand =%e\n", lattice[1].g_rand.c[0].imag);
    initialize(xxx_off);
    mat_invert_uml( F_OFFSET(g_rand), xxx_off, phi_off, mass, PRECISION, fn);

    //node0_printf("FINISHED FIRST INVERSION xxx1 =%e\n", lattice[1].xxx1.c[0].imag);
    //node0_printf("FINISHED FIRST INVERSION g_rand =%e\n", lattice[1].g_rand.c[0].imag);
 
    pbp(xxx_off, mass);
    
    d1trlnM_dmu1 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);
    //node0_printf("FINISHED FIRST TERM\n\n");
    d2trlnM_dmu2 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);
    //node0_printf("FINISHED SECOND TERM\n\n");
    d3trlnM_dmu3 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);
    d4trlnM_dmu4 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);
    d5trlnM_dmu5 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);
    d6trlnM_dmu6 (F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps);

    TR_M_inv_dMdu_deriv(F_OFFSET(g_rand), phi_off, xxx_off,  xxx1_off,  mass, jpbp_reps, npbp_reps); 

   //mat_invert_uml( xxx_off, xxx1_off, phi_off, mass, PRECISION, fn);
    d1trM_inv_dmu1(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
    d2trM_inv_dmu2(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
    d3trM_inv_dmu3(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
    d4trM_inv_dmu4(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
    d5trM_inv_dmu5(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
    d6trM_inv_dmu6(F_OFFSET(g_rand), phi_off, xxx_off, xxx1_off, mass, jpbp_reps, npbp_reps, -1);
  }

}

