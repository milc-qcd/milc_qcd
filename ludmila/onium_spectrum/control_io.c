#define CONTROL
#include <string.h>
#include "ks_includes.h"

/* Comment these out if you want to suppress detailed timing */
/*#define IOTIME*/
/*#define PRTIME*/

int main(int argc,char *argv[])
{
  void All_KS_hl_prop(field_offset snk, field_offset src, complex **propagator);
  int calculate_stag_prop();

  int prompt ,cgn,i;
  double starttime,endtime,dtime;
  site *s;
  float space_vol;

  int status, num_prop,t,color,spin, color1, spin1;
  int c,sp,c1,s1,pr;
  int key[4];

#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];

  complex pr_tmp; 
  complex *prop_rot[35];
  complex *prop_smear[35];
  ks_prop_file *spf;
  w_prop_file *fp_in_w;
  w_prop_file *fp_out_w;
  spin_wilson_vector localmat;
  wilson_propagator qtemp;

  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;
  space_vol = volume/nt;

  printf("beginning of program\n"); fflush(stdout);
  //setup_restrict_fourier(key, restrict);
  printf("start initialize_machine\n"); fflush(stdout);



  initialize_machine(argc,argv);
  printf("after initialize machine \n"); fflush(stdout);

  g_sync();
  /* set up */

  prompt = setup(); 
 
  printf("after set_up \n"); fflush(stdout);

  while( readin(prompt) == 0){
    
    starttime=dclock();
    pr=1;

     FORALLSITES(i,s){
       printf("(%d %d %d %d)\n)
	for(color=0;color<3;color++) for(color1=0;color1<3;color1++){




    if (pr==0){
    /**************************************************************/
    /* staggered test read*/
      spf =  r_serial_ks_fm_i(start_ks_prop_file);
      status = r_serial_ks_fm(spf, F_OFFSET(stag_propagator.e[0]), F_OFFSET(stag_propagator.e[1]),
			      F_OFFSET(stag_propagator.e[2]));   
      r_serial_ks_fm_f(spf);
			    /*write*/
      spf = r_serial_ks_i("/home/levkova/test_data/stag_prop/milc/qf_milc_stag_d_d_00001_m0.1000_nt0");

      //r_serial_ks_fm(spf, F_OFFSET(stag_propagator.e[0]), F_OFFSET(stag_propagator.e[1]),
      //	     F_OFFSET(stag_propagator.e[2]));
      
      for(color=0;color<3;color++){
	status = r_serial_ks(spf, color, F_OFFSET(stag_propagator_copy.e[color]));
      }
      FORALLSITES(i,s){
	for(color=0;color<3;color++) for(color1=0;color1<3;color1++){
	  
	  if((fabs(s->stag_propagator.e[color][color1].real - s->stag_propagator_copy.e[color][color1].real)>1e-8)||
	     fabs(s->stag_propagator.e[color][color1].imag- s->stag_propagator_copy.e[color][color1].imag)>1e-8){
	    printf("DIFFERENT VALUES IN STAG PROP! color %d, color1 %d\n",color,color1); terminate(1);
	  }
	  else{
	    printf("fm: %.9e %.9e\nml: %.9e %.9e\n\n",s->stag_propagator.e[color][color1].real,
		   s->stag_propagator.e[color][color1].imag,s->stag_propagator_copy.e[color][color1].real,
		   s->stag_propagator_copy.e[color][color1].imag);
	  }
	}
      }
      printf("COMPARISON COMPLETED!\n");
    }
    /**************************************************************/
    else{    
      /* Load the wilson heavy quark from a file*/ 
      //      r_prop_w_fm("/home/levkova/test_data/wilson_prop/arch_fm/qf_d_d_0.124_00001", F_OFFSET(quark_propagator_copy));

     //fp_out_w = w_open_prop(SAVE_SERIAL,"w_serial.out"); //in lattice.h

     /*for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
      {
#ifdef IOTIME
        save_propagator( SAVE_SERIAL, fp_out_w, 
                                    spin, color, F_OFFSET(quark_propagator.c[color].d[spin]),1);
#else
        save_propagator( SAVE_SERIAL,fp_out_w, 
                                    spin, color, F_OFFSET(quark_propagator.c[color].d[spin]),0);
#endif
      }

      w_close_prop(SAVE_SERIAL,fp_out_w);*/
    
      fp_in_w = r_open_prop(RELOAD_SERIAL,
			  "/home/levkova/test_data/wilson_prop/milc_serial/qf_milc_0.124_00001.e-17");
    printf("opened w_serial.out\n");

    for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
      {
#ifdef IOTIME
        status = reload_propagator( RELOAD_SERIAL,fp_in_w , 
                                    spin, color, F_OFFSET(quark_propagator.c[color].d[spin]),1);
#else
        status = reload_propagator( RELOAD_SERIAL,fp_in_w, 
                                    spin, color, F_OFFSET(quark_propagator.c[color].d[spin]),0);
#endif
      }
   printf("reloaded w_serial.out\n");
   r_close_prop(RELOAD_SERIAL,fp_in_w);
      
   //weyl2canopy_w_rot(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator_copy));
   rotate_wq(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator_copy));
      //canopy2weyl_w_rot(F_OFFSET(quark_propagator_copy),F_OFFSET(quark_propagator_copy));
      //weyl2canopy_w_rot(F_OFFSET(quark_propagator_copy),F_OFFSET(quark_propagator_copy));
    /*fp_in_w = r_open_prop(RELOAD_SERIAL,"/home/levkova/test_data/wilson_prop/milc_serial/qf_milc_0.124_00001.e-18");
    printf("opened w_serial.out\n");

    for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
      {
#ifdef IOTIME
        status = reload_propagator( RELOAD_SERIAL,fp_in_w , 
                                    spin, color, F_OFFSET(quark_propagator_copy.c[color].d[spin]),1);
#else
        status = reload_propagator( RELOAD_SERIAL,fp_in_w, 
                                    spin, color, F_OFFSET(quark_propagator_copy.c[color].d[spin]),0);
#endif
      }
   printf("reloaded w_serial.out\n");
   r_close_prop(RELOAD_SERIAL,fp_in_w);*/


    
      FORALLSITES(i,s) {
	/*	if(s->x==0&&s->y==0&&s->z==0&&s->t==0){  */
	  //hermitian conjugate
	/*for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	    for(color1=0;color1<3;color1++)for(spin1=0;spin1<4;spin1++){
	           qtemp.c[color].d[spin].d[spin1].c[color1].real =
		   s->quark_propagator_copy.c[color1].d[spin1].d[spin].c[color].real;
	           qtemp.c[color].d[spin].d[spin1].c[color1].imag =
		   s->quark_propagator_copy.c[color1].d[spin1].d[spin].c[color].imag;
	    }
	    s->quark_propagator_copy = qtemp;*/
	
      //gamma_5 multiplication
	/*r(color=0;color<3;color++){
	    mult_swv_by_gamma_l( &(s->quark_propagator.c[color]), &localmat, GAMMAFIVE);
	    mult_swv_by_gamma_r( &localmat, &(s->quark_propagator.c[color]), GAMMAFIVE);
	    }*/
	  //} 
	  /*else {*/
	  for(color=0;color<3;color++){
	    //mult_swv_by_gamma_l( &(s->quark_propagator.c[color]), &localmat, GAMMAFIVE);
	    mult_swv_by_gamma_l( &(s->quark_propagator.c[color]), &localmat, TUP);
	    //mult_swv_by_gamma_l( &localmat, &(s->quark_propagator.c[color]), TUP); 
	    //mult_swv_by_gamma_r( &(s->quark_propagator.c[color]), &localmat, GAMMAFIVE);  
  	    mult_swv_by_gamma_r( &localmat, &(s->quark_propagator.c[color]), TUP);        
	  }
	  //*/
      }
      
      //canopy2weyl_w_rot(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator));

     FORALLSITES(i,s) {   
       for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	for(color1=0;color1<3;color1++)for(spin1=0;spin1<4;spin1++)
	  if(fabs(s->quark_propagator.c[color].d[spin].d[spin1].c[color1].real -
		  s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real)>1e-6 ||
		 fabs(s->quark_propagator.c[color].d[spin].d[spin1].c[color1].imag -
		 s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag)>1e-6)
		{
		  printf("%d, %d, %d, %d (%d,%d,%d,%d) ",color,spin,spin1,color1,s->x,s->y,s->z,s->t);
		  printf("%.9e %.9e\n", s->quark_propagator.c[color].d[spin].d[spin1].c[color1].real,
			 s->quark_propagator.c[color].d[spin].d[spin1].c[color1].imag);
		  printf("                     %.9e %.9e\n\n", 
		   s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real,
		   s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag);
		  
		  terminate(1);
		}
       else{	     
	 printf("milc rot       %.9e %.9e\n", s->quark_propagator.c[color].d[spin].d[spin1].c[color1].real,
		s->quark_propagator.c[color].d[spin].d[spin1].c[color1].imag);
	 printf("canopy orig    %.9e %.9e\n\n", s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real,
		s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag);}
       
     }
	printf("COMPARISON COMPLETED!\n");
    }
  }
}
