#define CONTROL
#include <string.h>
#include "onium_includes.h"
#include "onium_generic.h"

/* Comment these out if you want to suppress detailed timing */
/*#define IOTIME*/
/*#define PRTIME*/

int main(int argc,char *argv[])
{
  void All_heavy_prop(field_offset snk, field_offset src, complex **propagator);


  int prompt ,k,ns,i;
  double starttime,endtime,dtime;
  site *s;
  float space_vol;

  int status, num_prop,t,color,spin, color1, spin1;

  int key[4];

#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];

  complex pr_tmp; 
  complex *prop_rot[35];
  complex *prop_smear[35];
  ks_prop_file *spf;
  w_prop_file *fp_in_w;
  wilson_propagator qtemp1;
 wilson_propagator *qdest;

  char *trace_kind[35]=
   {"P5_P5 p000","P5_P5 p100", "P5_P5 p110","P5_P5 p111", "P5_P5 p200",
    "P5_P5 p210","P5_P5 p211", "P5_P5 p220","P5_P5 p300", "P5_P5 p221",
    "P5_P5 p400","A4_P5 p000", "A4_P5 p100","A4_P5 p110", "A4_P5 p111",
    "A4_P5 p200","A1_P5 p100", "A1_P5 p110","A1_P5 p111", "A1_P5 p200",
    "P5_A1 p100","V1_V1 p000", "V1_V1 p100","V1_V1 p110", "V1_V1 p111",
    "V1_V1 p200","V1_V1 p210", "V1_V1 p211","V1_V1 p220", "V1_V1 p300",
    "V1_V1 p221","V1_V1 p400", "S_S p000"  ,"A1_A1 p000", "T23_T23 p000"}; 

  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;


  initialize_machine(argc,argv);
  g_sync();
  /* set up */
  prompt = setup(); 
  setup_restrict_fourier(key, restrict);
  space_vol = (float)(nx*ny*nz);

  while( readin(prompt) == 0){
    
    starttime=dclock();
  
    /**************************************************************/
    /*Allocate storage space for propagators*/
    for(num_prop=0;num_prop<35;num_prop++)
      {
	prop_rot[num_prop] = (complex *)malloc(nt*sizeof(complex));
	prop_smear[num_prop] = (complex *)malloc(nt*sizeof(complex));
	for(t=0;t<nt;t++){
	  prop_rot[num_prop][t].real = 0.0; 
	  prop_smear[num_prop][t].real = 0.0;
	  prop_rot[num_prop][t].imag = 0.0; 
	  prop_smear[num_prop][t].imag = 0.0; 
	}
      }
 
    /**************************************************************/
    // Load antiquark 
    if(a_format) { 
	  r_prop_w_fm(startfile_w[k], F_OFFSET(antiquark_propagator)); //loads in quark_propagator*/
	  /* Rotate to Milc format  V g0 G g0 V^+ */

	  FORALLSITES(i,s)
	  {	 
	    for(color=0;color<3;color++){
	      mult_swv_by_gamma_l( &(s->antiquark_propagator.c[color]), &(s->quark_propagator_copy.c[color]), TUP);
	      mult_swv_by_gamma_r( &(s->quark_propagator_copy.c[color]), &(s->antiquark_propagator.c[color]), TUP);
	    }

	  }
	  weyl2canopy_w_rot( F_OFFSET(antiquark_propagator),  F_OFFSET(antiquark_propagator));
    }
    else {
      fp_in_w = r_open_prop(RELOAD_SERIAL, a_startfile_w);
      
      for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	{
#ifdef IOTIME
	  status = reload_propagator( RELOAD_SERIAL,fp_in_w , 
				      spin, color, F_OFFSET(antiquark_propagator.c[color].d[spin]),1);
#else
	  status = reload_propagator( RELOAD_SERIAL,fp_in_w, 
				      spin, color, F_OFFSET(antiquark_propagator.c[color].d[spin]),0);
#endif
	}
      r_close_prop(RELOAD_SERIAL,fp_in_w); 
    }

    for(k=0; k<num_kap; k++){
      kappa = kap[k];
      if (format[k])
	{
	  /**************************************************************/
	  /* Load the wilson heavy quark from a file*/ 
	  r_prop_w_fm(startfile_w[k], F_OFFSET(quark_propagator)); //loads in quark_propagator*/
	  FORALLSITES(i,s)
	  {	 
	    for(color=0;color<3;color++){
	      mult_swv_by_gamma_l( &(s->quark_propagator.c[color]), &(s->quark_propagator_copy.c[color]), TUP);
	      mult_swv_by_gamma_r( &(s->quark_propagator_copy.c[color]), &(s->quark_propagator.c[color]), TUP);
	    }

	  }
	  weyl2canopy_w_rot( F_OFFSET(quark_propagator),  F_OFFSET(quark_propagator));
 
	}
      else{	
	/**************************************************************/
	/* Load (and Rotate) the heavy quarks. Result in quark_propagator_copy */    
	
	fp_in_w = r_open_prop(RELOAD_SERIAL, startfile_w[k] );
	
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
	r_close_prop(RELOAD_SERIAL,fp_in_w); 
      }

    rotate_w_quark(F_OFFSET(quark_propagator),F_OFFSET(quark_propagator_copy), d1[k]);  
    // result in quark_propagator
    rotate_w_quark(F_OFFSET(antiquark_propagator), F_OFFSET(antiquark_propagator_copy), a_d1);
 
    FORALLSITES(i,s) 
      {
	for(color=0;color<3;color++){
	  mult_swv_by_gamma_l( &(s->antiquark_propagator_copy.c[color]), &(qtemp1.c[color]), GAMMAFIVE);
	  mult_swv_by_gamma_r( &(qtemp1.c[color]), &(s->antiquark_propagator_copy.c[color]), GAMMAFIVE);
	}
      }     
    /**************************************************************/
    /*Calculate and print out the spectrum with the rotated heavy quark propagators*/
    
    All_heavy_prop(F_OFFSET(antiquark_propagator_copy), F_OFFSET(quark_propagator_copy), prop_rot);
    for(i=0;i<35;i++)
      {
	if(this_node==0) printf("\n\nTr %d, %s\n_________________________________\n",i, trace_kind[i]);
	for(t=0;t<nt;t++)
	  {
	    g_floatsum( &(prop_rot[i][t].real) );
	    g_floatsum( &(prop_rot[i][t].imag) );
	    if(this_node==0)
	      printf("%d %e %e\n", t,
		     prop_rot[i][t].real, prop_rot[i][t].imag);
	  }
      }
   
   for(i=0;i<35;i++)
	for(t=0;t<nt;t++){  prop_rot[i][t].real=0.0; prop_rot[i][t].imag = 0.0; } 

    /********************************************************************************************/
    /*Smear quarks, calculate and print out the spectrum with the smeared heavy quark propagators*/
   FORALLSITES(i,s) 
      {
	for(color=0;color<3;color++){
	  mult_swv_by_gamma_l( &(s->antiquark_propagator.c[color]), &(qtemp1.c[color]), GAMMAFIVE);
	  mult_swv_by_gamma_r( &(qtemp1.c[color]), &(s->antiquark_propagator.c[color]), GAMMAFIVE);
	}
      } 

   for(color=0;color<3;color++)for(spin=0;spin<4;spin++){
     restrict_fourier(F_OFFSET(quark_propagator.c[color].d[spin]),
		      F_OFFSET(mp), F_OFFSET(tmp),
		      sizeof(wilson_vector), FORWARDS);
      }

   for(ns=0; ns<num_smear;ns++){
	
     if(strcmp(smearfile[ns],"none")==0) continue;
     get_smearings_bi_serial(smearfile[ns]);
    
     restrict_fourier(F_OFFSET(w),
		      F_OFFSET(w1), F_OFFSET(w2),
		      sizeof(complex),FORWARDS);
     FORALLSITES(i,s){
       for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	 for(color1=0;color1<3;color1++)for(spin1=0;spin1<4;spin1++){
	      
	   pr_tmp = s->quark_propagator.c[color].d[spin].d[spin1].c[color1];
	      
	   s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real =
	     pr_tmp.real * s->w.real - pr_tmp.imag * s->w.imag;
	      
	   s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag =
	     pr_tmp.real * s->w.imag + pr_tmp.imag * s->w.real;
	 }
     }
	
     for(color=0;color<3;color++)for(spin=0;spin<4;spin++){
       restrict_fourier(F_OFFSET(quark_propagator_copy.c[color].d[spin]),
			F_OFFSET(mp), F_OFFSET(tmp),
			sizeof(wilson_vector), BACKWARDS);
     }	


     All_heavy_prop(F_OFFSET(antiquark_propagator), F_OFFSET(quark_propagator_copy), prop_smear);
     
     for(i=0;i<35;i++){
       if(this_node==0) printf("\n\nSMEAR_#%d, %s_k_%f\n_________________________________\n", 
			       ns, trace_kind[i],kap[k]);
       for(t=0;t<nt;t++)
	 {
	   g_floatsum( &prop_smear[i][t].real );
	   g_floatsum( &prop_smear[i][t].imag );
	   if(this_node==0)
	     printf("%d %e %e\n", t,
		    prop_smear[i][t].real/space_vol, prop_smear[i][t].imag/space_vol);
	 }
     }
     for(i=0;i<35;i++)
       for(t=0;t<nt;t++){ prop_smear[i][t].real=0.0; prop_smear[i][t].imag = 0.0; }    
     
   }/* loop ns*/
   
    }/*loop kappa*/
    
    
  }
}

