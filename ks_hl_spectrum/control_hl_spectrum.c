/***************** control_hl_spectrum.c ***********************************/

/* Main procedure for heavy-light mesons */
/* MIMD version 7 */

#define CONTROL
#include "ks_hl_spectrum_includes.h"
#include <string.h>

/*--------------------------------------------------------------------*/

int main(int argc,char *argv[])
{
  int prompt , k, ns, i;
  site *s;
  double inv_space_vol;
  
  int color,spin, color1, spin1;
  
  int key[4];
  int dummy[4];
  FILE *corr_fp;
  
  complex pr_tmp; 
  wilson_propagator *qdest;
  wilson_propagator qtemp1;

  wilson_vector *psi = NULL;
  w_prop_file *wpf;
  wilson_quark_source wqs;
  
  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;

  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  prompt = setup(); 
  setup_restrict_fourier(key, dummy);

  psi = create_wv_field();

  /* Initialize the source type */
  init_wqs(&wqs);

  while( readin(prompt) == 0){
    
    
    /**************************************************************/
    /*load staggered propagator*/
    
    reload_ksprop_to_site3(ks_prop_startflag, 
			   start_ks_prop_file, &ksqs, F_OFFSET(prop), 1);
    
    FORALLSITES(i,s){
      for(color = 0; color < 3; color++)for(k = 0; k < 3; k++)
	s->stag_propagator.e[color][k] = s->prop[color].c[k];
    }
    
    /* Initialize FNAL correlator file */
    
    corr_fp = open_fnal_meson_file(savefile_c);

    /* Load Wilson propagator for each kappa */

    for(k=0; k<num_kap; k++){
      kappa = kap[k];
      wpf = r_open_wprop(startflag_w[k], startfile_w[k]);
      for(spin=0;spin<4;spin++)
	for(color=0;color<3;color++){
	  if(reload_wprop_sc_to_field(startflag_w[k], wpf,
				      &wqs, spin, color, psi, 1) != 0)
	    terminate(1);
	  FORALLSITES(i,s){
	    copy_wvec(&psi[i],&lattice[i].quark_propagator.c[color].d[spin]);
	  }
	}
      r_close_wprop(startflag_w[k],wpf);

      /*******************************************************************/
      /* Rotate the heavy quark */
      
      rotate_w_quark(F_OFFSET(quark_propagator), 
		     F_OFFSET(quark_propagator_copy), d1[k]);  
      // result in quark_propagator_copy


      /**************************************************************/
      /*Calculate and print out the spectrum with the rotated heavy
        quark propagators*/

      spectrum_hl_rot(corr_fp, F_OFFSET(stag_propagator), 
		      F_OFFSET(quark_propagator_copy), k);
      
      
      /**************************************************************/
      /*Smear quarks, calculate and print out the spectrum with the
        smeared heavy quark propagators*/
      
      for(color=0;color<3;color++)for(spin=0;spin<4;spin++){
	restrict_fourier_site(F_OFFSET(quark_propagator.c[color].d[spin]),
			      sizeof(wilson_vector), FORWARDS);
      }
      
      for(ns=0; ns<num_smear;ns++){
	if(strcmp(smearfile[ns],"none")==0) continue;

	inv_space_vol = 1./((double)nx*ny*nz);

	/* Either read a smearing file, or take it to be a point sink */
	if(strlen(smearfile[ns]) != 0){

	   get_smearings_bi_serial(smearfile[ns]);
	
	   restrict_fourier_site(F_OFFSET(w),
				 sizeof(complex), FORWARDS);
	
	   FORALLSITES(i,s){
	     for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	      for(color1=0;color1<3;color1++)for(spin1=0;spin1<4;spin1++){
		  pr_tmp = 
		    s->quark_propagator.c[color].d[spin].d[spin1].c[color1];
	  
		  s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real =
		    pr_tmp.real * s->w.real - pr_tmp.imag * s->w.imag;
	      
		  s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag =
		    pr_tmp.real * s->w.imag + pr_tmp.imag * s->w.real;
		}
	   }
	  } else { /* Point sink */
	   FORALLSITES(i,s){
	     for(color=0;color<3;color++)for(spin=0;spin<4;spin++)
	      for(color1=0;color1<3;color1++)for(spin1=0;spin1<4;spin1++){
		  pr_tmp = 
		    s->quark_propagator.c[color].d[spin].d[spin1].c[color1];
		  
		  s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].real =
		    pr_tmp.real;
	      
		  s->quark_propagator_copy.c[color].d[spin].d[spin1].c[color1].imag =
		    pr_tmp.imag;
		}
	   }
	  }
	
	  for(color=0;color<3;color++)for(spin=0;spin<4;spin++){
	      restrict_fourier_site(F_OFFSET(quark_propagator_copy.c[color].d[spin]),
				    sizeof(wilson_vector), BACKWARDS);
	    }	
	  
	  FORALLSITES(i,s)
	  {
	    qdest = &(s->quark_propagator_copy);
	    qtemp1 = s->quark_propagator_copy;
	    for(spin=0;spin<4;spin++)for(color=0;color<3;color++)
 	      for(spin1=0;spin1<4;spin1++)for(color1=0;color1<3;color1++)
	      {
		qdest->c[color].d[spin1].d[spin].c[color1].real = 
		  qtemp1.c[color].d[spin].d[spin1].c[color1].real;
		qdest->c[color].d[spin1].d[spin].c[color1].imag = 
		  qtemp1.c[color].d[spin].d[spin1].c[color1].imag;
	      }
	  }
	  
	
      /**************************************************************/
      /*Calculate and print out the spectrum with the smeared sink */

	spectrum_hl_smear(corr_fp, F_OFFSET(stag_propagator), 
			  F_OFFSET(quark_propagator_copy), k, ns);
	
      }/* loop ns*/
      
    }/*loop kappa*/

    close_fnal_meson_file(corr_fp);
  } /* readin(prompt) == 0 */

  destroy_wv_field(psi);
  node0_printf("\nRUNNING COMPLETED\n"); fflush(stdout);

#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  
  normal_exit(0);
  return 0;
}
