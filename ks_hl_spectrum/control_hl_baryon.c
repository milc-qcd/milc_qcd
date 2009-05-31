/***************** control_hl_baryon.c **************************************/

/* Main procedure for heavy-light-light baryons */
/* MIMD version 7 */

/* Baryon two point function using Staggered light quark propagators
 * from Fermilab and Wilson heavy quark propagators, based on
 * Ludmila's ks_hl_spectrum code. */
/* 2.25.2006 Heechang Na */

/* This code is modified to print every dirac components of two-point
   functions.  But, it uses only one clover quark for heavy.
   11.14.2006 Heechang. */
/* This code take care of the correct spin indices for source and
   sink.  * * But, it still have the -gamma_mu convention as a
   convention.  * * 
   12.7.2007 Heechang */

/* This is a code for multi light masses and multi operators.  The
   output file has all two point functions for multiful cases, so that
   it needs another program to analyze.  
   2.21.2008 Heechang */

#define CONTROL
#include <string.h>
#include "ks_hl_spectrum_includes.h"

/*--------------------------------------------------------------------*/

int main(int argc,char *argv[])
{
  int prompt , k, i, j, strange, light, mu, nu;
  double starttime;
  site *s;
  double space_vol;
  
  int t, color, spin, spin1;

  /*file for output*/
  FILE *fp_out;
  char fname[100];
 
  /* propagators for O5, and Omu type operators */
  double_complex *(prop_5[4][4]), *(prop_mu[4][4][3][3]);
  
  ks_quark_source ksqs, ksqs2;
  wilson_quark_source wqs;
  ks_prop_field ksp;
  
  initialize_machine(&argc,&argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();
  prompt = setup(); 
  space_vol = (double)(nx*ny*nz);

  /* Initialize the source structures */
  init_wqs(&wqs);
  init_ksqs(&ksqs);
  init_ksqs(&ksqs2);

  while( readin(prompt) == 0){
    
    starttime=dclock();
    
    
    /**************************************************************/
    /*Allocate storage space for propagators*/
    for(spin=0;spin<4;spin++) 
      for(spin1=0;spin1<4;spin1++)
	{
	  prop_5[spin][spin1] = 
	    (double_complex *)malloc(nt*sizeof(double_complex));
	  for(t=0;t<nt;t++){
	    prop_5[spin][spin1][t].real = 0.0; 
	    prop_5[spin][spin1][t].imag = 0.0; 
	    
	    for(i=0;i<3;i++)for(j=0;j<3;j++){
		prop_mu[spin][spin1][i][j] = 
		  (double_complex *)malloc(nt*sizeof(double_complex));
		for(t=0;t<nt;t++){
		  prop_mu[spin][spin1][i][j][t].real = 0.0;
		  prop_mu[spin][spin1][i][j][t].imag = 0.0;
		}
	      }
	  }
	}
    
    node0_printf("BEGIN\n");

    /**************************************************************/
    /* Load Wilson propagator for each kappa (quark_propagator) */

    for(k=0; k<num_kap; k++){
      kappa = kap[k];

      /**************************************************************/
      /* Load the fm-style wilson heavy quark from a file*/ 
      
      reload_wprop_to_site(startflag_w[k], startfile_w[k], &wqs,
			   F_OFFSET(quark_propagator), 1);
      /******************************************************************/
      /*Baryon needs two staggered propagators. -H */
      /*This is for multiple light qaurks and strange quarks,
	so it takes loops for num_strange and num_light  -H 2.22.2008*/
      
      for(strange=0; strange<num_strange; strange++){
	
	/* Staggered strange quark (stag_strange_propagator) */
	ksp = create_ksp_field();
	reload_ksprop_to_ksp_field(start_ks_strange_flag[strange], 
				   start_ks_strange_file[strange], &ksqs, 
				   ksp, 1);
	FORALLSITES(i,s){
	  for(color = 0; color < 3; color++)
	    for(k = 0; k < 3; k++)
	      s->stag_strange_propagator.e[color][k] = ksp[color][i].c[k];
	}
	
	destroy_ksp_field(ksp);

	for(light=0; light<num_light; light++){

	  /* Staggered light quark (stag_light_propagator) */
	  ksp = create_ksp_field();
	  reload_ksprop_to_ksp_field(start_ks_light_flag[light], 
				     start_ks_light_file[light], &ksqs, 
				     ksp, 1);
	  FORALLSITES(i,s){
	    for(color = 0; color < 3; color++)
	      for(k = 0; k < 3; k++)
		s->stag_light_propagator.e[color][k] = ksp[color][i].c[k];
	  }
	  
	  destroy_ksp_field(ksp);
	  

	  /**************************************************************/
	  /*Calculate only light valence quarks: O5_ll Omu_ll Omu_lHH   */
	  /*It's only with first strange quark, because there is        */
	  /*no strange qaurks. -Heechang		                  					*/
	  /**************************************************************/
	  
	  if(strange==0){
	    
	    /* O5_ll ************************************************/
	    
	    /*get 2pt*/
	  ks_baryon_2point_O5(F_OFFSET(stag_light_propagator), 
			      F_OFFSET(stag_light_propagator), 
			      F_OFFSET(quark_propagator), prop_5);
	  
	  /*print out the results*/
	  if(this_node==0) {
	    printf("\n\n2Pt O5_ll m%f \n_________________________________\n", 
		   m_light[light]);
	    /*make name of file and open it*/
	    sprintf(fname,"O5_ll.m%f.out",m_light[light]);
	    fp_out = fopen(fname,"w");
	  }

	  for(t=0;t<nt;t++)
	    {
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  g_doublesum( &(prop_5[spin][spin1][t].real) );
		  g_doublesum( &(prop_5[spin][spin1][t].imag) );
		}
	    }
	  g_sync();
	  if(this_node==0){
	    for(t=0;t<nt;t++){
	      fprintf(fp_out," %d ", t);
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  fprintf(fp_out,"%e %e ", prop_5[spin][spin1][t].real, 
			 prop_5[spin][spin1][t].imag);
		}
	      fprintf(fp_out,"\n");
	    }
	  }
	  
	  for(t=0;t<nt;t++) 
	    for(spin1=0;spin1<4;spin1++)
	      for(spin=0;spin<4;spin++) 
		{
		  prop_5[spin][spin1][t].real=0.0; 
		  prop_5[spin][spin1][t].imag=0.0; 
		}
		if(this_node==0)fclose(fp_out);
	  g_sync();
	  
	  /* Omu_ll **************************************************/
	  
	  ks_baryon_2point_Omu(F_OFFSET(stag_light_propagator), 
			       F_OFFSET(stag_light_propagator), 
			       F_OFFSET(quark_propagator), prop_mu);
	  
	  if(this_node==0) {
	    printf("\n\n2Pt Omu_ll m%f \n_________________________________\n", 
		   m_light[light]);
	    /*make name of file and open it*/
	    sprintf(fname,"Omu_ll.m%f.out",m_light[light]);
	    fp_out = fopen(fname,"w");
	  }
	  for(t=0;t<nt;t++)
	    {
	      for(i=0;i<3;i++){
		for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		    g_doublesum( &(prop_mu[spin][spin1][i][i][t].real) );
		    g_doublesum( &(prop_mu[spin][spin1][i][i][t].imag) );
		  }
	      }
	    }
	  g_sync();
	  if(this_node==0){
	    for(i=0;i<3;i++){
	      for(t=0;t<nt;t++){
		fprintf(fp_out," %d %d %d ", i,i,t);
		for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		    fprintf(fp_out,"%e %e ", prop_mu[spin][spin1][i][i][t].real, 
			   prop_mu[spin][spin1][i][i][t].imag);
		  }
		fprintf(fp_out,"\n");
	      }
	    }
	  }
	  
	  for(t=0;t<nt;t++) 
	    for(i=0;i<3;i++) 
	      for(j=0;j<3;j++) 
		for(spin1=0;spin1<4;spin1++)
		  for(spin=0;spin<4;spin++) 
		    {
		      prop_mu[spin][spin1][i][j][t].real=0.0; 
		      prop_mu[spin][spin1][i][j][t].imag=0.0;
		    }
	  if(this_node==0)fclose(fp_out);
	  g_sync();
	  
	  /*  Omu_lHH ***************************************************/
	  
	  ks_baryon_2point_Omu_HH(F_OFFSET(stag_light_propagator), 
				  F_OFFSET(quark_propagator), prop_mu);
	  
	  if(this_node==0) {
	    printf("\n\n2Pt Omu_lHH m%f \n_________________________________\n", 
		   m_light[light]);
	    /*make name of file and open it*/
	    sprintf(fname,"Omu_lHH.m%f.out",m_light[light]);
	    fp_out = fopen(fname,"w");
	  }
	  for(t=0;t<nt;t++)
	    {
	      for(mu=0;mu<3;mu++)for(nu=0;nu<3;nu++){
		  for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		      g_doublesum(&(prop_mu[spin][spin1][mu][nu][t].real));
		      g_doublesum(&(prop_mu[spin][spin1][mu][nu][t].imag));
		    }
		}
	    }
	  g_sync();
	  if(this_node==0){
	    for(mu=0;mu<3;mu++)for(nu=0;nu<3;nu++){
		for(t=0;t<nt;t++){
		  fprintf(fp_out," %d %d %d ", mu ,nu ,t);
		  /* Note: Here the mu and nu are transposed. */
		  for(spin=0;spin<4;spin++)for(spin1=0;spin1<4;spin1++) {
		      fprintf(fp_out,"%e %e ", 
			     prop_mu[spin][spin1][nu][mu][t].real, 
			     prop_mu[spin][spin1][nu][mu][t].imag);
		    }
		  fprintf(fp_out,"\n");
		}
	      }
	  }
	  
	  for(t=0;t<nt;t++) 
	    for(mu=0;mu<3;mu++)
	      for(nu=0;nu<3;nu++)
		for(spin1=0;spin1<4;spin1++)
		  for(spin=0;spin<4;spin++) 
		    {
		      prop_mu[spin][spin1][mu][nu][t].real=0.0; 
		      prop_mu[spin][spin1][mu][nu][t].imag=0.0;
		    }
	  if(this_node==0)fclose(fp_out);
	  g_sync();
	  
	  }/*for the light only block */ 	
	  
	  
	  /**************************************************************/
	  /*Calculate light and strange valence quarks: O5_ls Omu_ls    */
	  /**************************************************************/
	  
	  /* O5_ls ********************************************************/
	  
	  ks_baryon_2point_O5(F_OFFSET(stag_light_propagator), 
			      F_OFFSET(stag_strange_propagator), 
			      F_OFFSET(quark_propagator), prop_5);
	  
	  /*print out the results*/
	  if(this_node==0) {
	    printf("\n\n2Pt O5_ls m%f  s%f \n_________________________________\n", 
		   m_light[light], m_strange[strange]);
	    /*make name of file and open it*/
	    sprintf(fname,"O5_ls.m%f.s%f.out",m_light[light],m_strange[strange]);
	    fp_out = fopen(fname,"w");
	  }
	  for(t=0;t<nt;t++)
	    {
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  g_doublesum( &(prop_5[spin][spin1][t].real) );
		  g_doublesum( &(prop_5[spin][spin1][t].imag) );
		}
	    }
	  g_sync();
	  if(this_node==0){
	    for(t=0;t<nt;t++){
	      fprintf(fp_out," %d ", t);
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  fprintf(fp_out,"%e %e ", prop_5[spin][spin1][t].real, 
			 prop_5[spin][spin1][t].imag);
		}
	      fprintf(fp_out,"\n");
	    }
	  }
	  
	  for(t=0;t<nt;t++) 
	    for(spin1=0;spin1<4;spin1++)
	      for(spin=0;spin<4;spin++) 
		{
		  prop_5[spin][spin1][t].real=0.0; 
		  prop_5[spin][spin1][t].imag=0.0;
		}
	  if(this_node==0)fclose(fp_out);
	  g_sync();
	  
	  /* Omu_ls *******************************************************/
	  
	  ks_baryon_2point_Omu(F_OFFSET(stag_light_propagator), 
			       F_OFFSET(stag_strange_propagator), 
			       F_OFFSET(quark_propagator), prop_mu);
	  
	  if(this_node==0) {
	    printf("\n\n2Pt Omu_ls m%f  s%f \n_________________________________\n", 
		   m_light[light], m_strange[strange]);
	    /*make name of file and open it*/
	    sprintf(fname,"Omu_ls.m%f.s%f.out",m_light[light],m_strange[strange]);
	    fp_out = fopen(fname,"w");
	  }
	  for(t=0;t<nt;t++)
	    {
	      for(i=0;i<3;i++){
		for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		    g_doublesum( &(prop_mu[spin][spin1][i][i][t].real) );
		    g_doublesum( &(prop_mu[spin][spin1][i][i][t].imag) );
		  }
	      }
	    }
	  g_sync();
	  if(this_node==0){
	    for(i=0;i<3;i++){
	      for(t=0;t<nt;t++){
		fprintf(fp_out," %d %d %d ", i,i,t);
		for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		    fprintf(fp_out,"%e %e ", prop_mu[spin][spin1][i][i][t].real, 
			   prop_mu[spin][spin1][i][i][t].imag);
		  }
		fprintf(fp_out,"\n");
	      }
	    }
	  }
	  
	  for(t=0;t<nt;t++) 
	    for(i=0;i<3;i++) 
	      for(j=0;j<3;j++) 
		for(spin1=0;spin1<4;spin1++)
		  for(spin=0;spin<4;spin++) 
		    {
		      prop_mu[spin][spin1][i][j][t].real=0.0; 
		      prop_mu[spin][spin1][i][j][t].imag=0.0;
		    }
	  if(this_node==0)fclose(fp_out);
	  g_sync();
	  
	  
	}/*loop for num_light*/
	
	
	
	/**********************************************/
	/*only strange valence quarks. Omu_ss Omu_sHH */
	/**********************************************/
	
	/* Omu_ss *************************************************************/
	
	ks_baryon_2point_Omu(F_OFFSET(stag_strange_propagator), 
			     F_OFFSET(stag_strange_propagator), 
			     F_OFFSET(quark_propagator), prop_mu);
	
	if(this_node==0){
	  printf("\n\n2Pt Omu_ss s%f \n_________________________________\n", 
		 m_strange[strange]);
	  /*make name of file and open it*/
	  sprintf(fname,"Omu_ss.s%f.out",m_strange[strange]);
	  fp_out = fopen(fname,"w");
	}
	for(t=0;t<nt;t++)
	  {
	    for(i=0;i<3;i++){
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  g_doublesum( &(prop_mu[spin][spin1][i][i][t].real) );
		  g_doublesum( &(prop_mu[spin][spin1][i][i][t].imag) );
		}
	    }
	  }
	g_sync();
	if(this_node==0){
	  for(i=0;i<3;i++){
	    for(t=0;t<nt;t++){
	      fprintf(fp_out," %d %d %d ", i,i,t);
	      for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		  fprintf(fp_out,"%e %e ", prop_mu[spin][spin1][i][i][t].real, 
			 prop_mu[spin][spin1][i][i][t].imag);
		}
	      fprintf(fp_out,"\n");
	    }
	  }
	}
	
	for(t=0;t<nt;t++) 
	  for(i=0;i<3;i++) 
	    for(j=0;j<3;j++) 
	      for(spin1=0;spin1<4;spin1++)
		for(spin=0;spin<4;spin++) 
		  {
		    prop_mu[spin][spin1][i][j][t].real=0.0; 
		    prop_mu[spin][spin1][i][j][t].imag=0.0;
		  }
	if(this_node==0)fclose(fp_out);
	g_sync();
	
	/* Omu_sHH  ******************************************************/
	
	
	ks_baryon_2point_Omu_HH(F_OFFSET(stag_strange_propagator), 
				F_OFFSET(quark_propagator), prop_mu);
	
	if(this_node==0) {
	  printf("\n\n2Pt Omu_sHH s%f \n_________________________________\n", 
		 m_strange[strange]);
	  /*make name of file and open it*/
	  sprintf(fname,"Omu_sHH.s%f.out",m_strange[strange]);
	  fp_out = fopen(fname,"w");
	}
	for(t=0;t<nt;t++)
	  {
	    for(mu=0;mu<3;mu++)for(nu=0;nu<3;nu++){
		for(spin1=0;spin1<4;spin1++)for(spin=0;spin<4;spin++) {
		    g_doublesum( &(prop_mu[spin][spin1][mu][nu][t].real) );
		    g_doublesum( &(prop_mu[spin][spin1][mu][nu][t].imag) );
		  }
	      }
	  }
	g_sync();
	if(this_node==0){
	  for(mu=0;mu<3;mu++)for(nu=0;nu<3;nu++){
	      for(t=0;t<nt;t++){
		fprintf(fp_out," %d %d %d ", mu ,nu ,t);
		/* Note: Here the mu and nu are transposed. */
		for(spin=0;spin<4;spin++)for(spin1=0;spin1<4;spin1++) {
		    fprintf(fp_out,"%e %e ", prop_mu[spin][spin1][nu][mu][t].real, 
			   prop_mu[spin][spin1][nu][mu][t].imag);
		  }
		fprintf(fp_out,"\n");
	      }
	    }
	}
	
	for(t=0;t<nt;t++) 
	  for(mu=0;mu<3;mu++)
	    for(nu=0;nu<3;nu++)
	      for(spin1=0;spin1<4;spin1++)
		for(spin=0;spin<4;spin++) 
		  {
		    prop_mu[spin][spin1][mu][nu][t].real=0.0; 
		    prop_mu[spin][spin1][mu][nu][t].imag=0.0;
		  }
	if(this_node==0)fclose(fp_out);
	g_sync();
	
      }/*loop for strange quarks*/
    } /*loop kappa*/
  } /*while(readin) */

  node0_printf("\nRUNNING COMPLETED\n"); fflush(stdout);

#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  
  normal_exit(0);
  return 0;
}

