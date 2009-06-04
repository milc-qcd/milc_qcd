/***************** control_mr.c *****************************************/
/* MIMD version 7 */
/*  set tabstop=2   for easy reading of this file */

/* Main procedure for quenched heavy-light
  SU3 Wilson fermions --- hopping param expansion
  for heavies 	

  This original code has been modified to calculate the a smearing
  matrix of static-light correlators. The smeared meson operators used
  in the calculation of the static B_B parameter are also calculated.

  The static-light variational code will only be active, when one
  light quark inversion is done (nkap=1).

  WARNING :: the code that Claude originally sent me had an internal loop over
  the kappa values inside the code. After the modification that I have made to the
  code, it shoule only be used for one kappa value at one time.

*/
/* MIMD version 7 */

/*
 * Background gauge field is supplied, assumed to be in coulomb gauge 
 */

#define CONTROL

#include "w_static_includes.h"
#include "lattice_qdp.h"



int main(int argc, char **argv)
{
  int meascount[MAX_NKAP];
  int prompt, count1, count2;
  Real avm_iters[MAX_NKAP];
  double starttime, endtime;

  int MaxMR, restart_flag;
  Real RsdMR;		/******/

  int spin, color, nk;		/******/
  int max_prop;

  int cl_cg = CL_CG;
  double ssplaq, stplaq;


  FILE *fp_m_out = NULL;  /*** meson IO stuff **/
  int fb_m_out = 0;	   /*** meson IO stuff **/

  w_prop_file *fp_in_w[MAX_NKAP];        /* For propagator files */
  w_prop_file *fp_out_w[MAX_NKAP];       /* For propagator files */
  
  double g_time ; 

  int i ;
  int MinMR;
  
/*** variables required for the static variational code ***/
  int nodata = 0 ;
  complex *meson = NULL;

/****** start of the execution of the code ************/


  initialize_machine(&argc, &argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  g_sync();

  /* set up */
  prompt = setup_h();

  /**DEBUG***/  
#ifdef DEBUGDEF
  light_quark_pion(0) ;
#endif

    /* loop over input sets */
  while( readin(prompt) == 0)
  {

    if( fixflag == COULOMB_GAUGE_FIX)
    {
      if(this_node == 0) 
	printf("Fixing to Coulomb gauge\n");
      g_time = -dclock();

      gaugefix(TUP,(Real)1.5,500,GAUGE_FIX_TOL);

      g_time += dclock();
      if(this_node==0)printf("Time to gauge fix = %e\n",g_time);
      invalidate_this_clov(gen_clov);
      
    }
    else
      if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");

    /* save lattice if requested */
    if( saveflag != FORGET )
    {
      save_lattice( saveflag, savefile, stringLFN );
    }


    /* call plaquette measuring process */
    d_plaquette(&ssplaq, &stplaq);
    if (this_node == 0)
      printf("START %e %e\n",(double) ssplaq, (double) stplaq);


    /******* set up code for the static variational calculation *****/
    if( nkap == 1 )
    {
      nodata = nt*nosmear*144 ;
      /** reserve memory for the smeared meson correlators on each node ****/
      if( ( meson = (complex *) calloc( (size_t) nodata, sizeof(complex) )  ) == NULL )
      {
	printf("ERROR: could not reserve buffer space for the meson smearing functions\n");
	terminate(1);
      }
      
      /** call a number of set up routines for the static variational code ***/
      setup_vary(meson, nodata);
      
    }  /*** end of set up section for the static-variational calculation ***/

    /***DEBUG check_calc_matrix() ;   ****/  


    starttime=dclock();



    MaxMR = niter;
    RsdMR = (Real) sqrt((double) rsqprop);

    if (this_node == 0)
      printf("Residue=%e\n",(double) RsdMR);
	     

    for (nk = 0; nk < nkap; nk++)
    {
      avm_iters[nk] = 0.0;
      meascount[nk] = 0;
    }
    max_prop = 12;
    count1 = 0;
    count2 = 0;


    for (spin = start_spin; spin < 4; spin++)
    {

      for (color = 0; color < 3; color++)
      {

	count1++;
	if (count1 == 1)
	  color += start_color;

	for (nk = 0; nk < nkap; nk++)
	{

	  count2++;
	  if (count2 == 1)
	    nk += start_kap;


	  kappa = cappa[nk];

	  meascount[nk]++;

	  /* open file for wilson propagators */
	  fp_in_w[nk]  = r_open_wprop(startflag_w[nk], startfile_w[nk]);

	  if ((spin + color) == 0)
	  {
	    /*** first pass of the code  **/
	    fp_out_w[nk] = w_open_wprop(saveflag_w[nk],  savefile_w[nk],
					wqs.type);

	    /* open file for meson output and write the header */
	    if (saveflag_m == SAVE_MESON_ASCII)
	    {
	      fp_m_out = w_ascii_m_i(savefile_m[nk], max_prop);
	      fb_m_out = -1;	/* i.e. file is NOT binary */
	    }
	    else if (saveflag_m == SAVE_MESON_BINARY)
	    {
	      fb_m_out = w_binary_m_i(savefile_m[nk], max_prop);
	      fp_m_out = NULL;	/* i.e. file is NOT ascii */
	    }
	    else
	    {
	      if( this_node == 0 ) 
		printf("ERROR in main saveflag_m = %d is out of range in initial opening\n",saveflag_m)  ;
	      terminate(1); 
	    }


	  } /*** end of spin =0 && color == 0 **/
	  else
	  {
	    fp_out_w[nk] = w_open_wprop(saveflag_w[nk],  savefile_w[nk],
					wqs.type);

	    /* open file for meson output for appending output*/
	    if (saveflag_m == SAVE_MESON_ASCII)
	    {
	      fp_m_out = a_ascii_m_i(savefile_m[nk], max_prop);
	      fb_m_out = -1;	/* i.e. file is NOT binary */
	    }
	    if (saveflag_m == SAVE_MESON_BINARY)
	    {
	      fb_m_out = a_binary_m_i(savefile_m[nk], max_prop);
	      fp_m_out = NULL;	/* i.e. file is NOT ascii */
	    }
	    else
	    {
	      if( this_node == 0 ) 
		printf("ERROR in main saveflag_m = %d is out of range in appending opening\n",saveflag_m)  ;
	      terminate(1); 
	    }



	  }  /*** end of spin && color not equal to zero ***/


	  if (this_node == 0)
	    printf("color=%d spin=%d kappa=%f nk=%d\n", color, spin, (double) kappa, nk);

	  /* load psi if requested */
	  init_wqs(&wqstmp2);
	  reload_wprop_sc_to_site(startflag_w[nk], fp_in_w[nk],&wqstmp2,
			    spin, color, F_OFFSET(psi),1);

	  if (nk == 0 || count2 == 1)
	    restart_flag = flag;
	  else
	    restart_flag = 1;

	  
	  /* Conjugate gradient inversion uses site structure
	     temporary"chi" */
	  
	  
	  /* Complete the source structure */
	  wqs.color = color;
	  wqs.spin = spin;
	  wqs.parity = EVENANDODD;
	  
	  /* For wilson_info */
	  wqstmp = wqs;
	  
	  /* If we are starting fresh, we want to set a mininum number of
	     iterations */
	  if(startflag_w[nk] == FRESH)MinMR = nt/2; else MinMR = 0;

	  /* Load inversion control structure */
	  qic.min = MinMR;
	  qic.max = MaxMR;
	  qic.nrestart = nrestart;
	  qic.resid = RsdMR;
	  qic.start_flag = restart_flag;
	    
	  /* Load Dirac matrix parameters */
	  dwp.Kappa = kappa;
	  
	  switch (cl_cg) {
	  case CG:
	    /* Load temporaries specific to inverter */
	    
	    /* compute the propagator.  Result in psi. */
	    avm_iters[nk] += 
	      (Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
					   w_source_h,&wqs,
					   cgilu_w_site,&qic,(void *)&dwp);
	    break;
	  case MR:
	    /* Load temporaries specific to inverter */
	    
	    /* compute the propagator.  Result in psi. */
	    avm_iters[nk] += 
	      (Real)wilson_invert_site_wqs(F_OFFSET(chi),F_OFFSET(psi),
					   w_source_h,&wqs,
					   mrilu_w_site,&qic,(void *)&dwp);
		break;
	      default:
		node0_printf("main(%d): Inverter choice %d not supported\n",
			     this_node,cl_cg);
	  }
	  
	  /* save psi if requested */
	  save_wprop_sc_from_site( saveflag_w[nk],fp_out_w[nk], &wqstmp2,
			  spin,color,F_OFFSET(psi),1);


	  light_meson(F_OFFSET(psi), color, spin, wqs.type, fp_m_out, fb_m_out);

	  if (this_node == 0)
	    printf("Light mesons found\n");

	/*** calculate the correlators required for the static variational code **/
	if( nkap == 1 )
	{

	  /** calculate the smeared meson correlators required for Bparam **/
	  calc_smeared_meson(meson, F_OFFSET(psi) ,  F_OFFSET(mp), color, spin);
			 
	  /** calculate the object required for the 2-pt variational calculation ***/
	  buildup_strip(F_OFFSET(psi)   ,  color,  spin); 

	} /** end of the partial calculations for the variationl project ***/



	  /*
	   * find source again since mrilu overwrites it; for hopping
	   * expansion 
	   */
	  /* source must be of definite parity */

	  wqs.parity = source_parity;
	  w_source_h(F_OFFSET(chi), &wqs);

	  hopping(F_OFFSET(chi), F_OFFSET(mp), F_OFFSET(psi), nhop,
		  kappa_c, wqs.parity, color, spin, wqs.type,
		  fp_m_out, fb_m_out);

	  /* close files */

	  r_close_wprop(startflag_w[nk], fp_in_w[nk]);
	  w_close_wprop(saveflag_w[nk],fp_out_w[nk]);

	  if (saveflag_m == SAVE_MESON_ASCII)
	    w_ascii_m_f(fp_m_out, savefile_m[nk]);
	  else if (saveflag_m == SAVE_MESON_BINARY)
	    w_binary_m_f(fb_m_out, savefile_m[nk]);


	  if (spin == end_spin && color == end_color && nk == end_kap)
	    goto end_of_loops;


	}
      }
    }				/* end of loop over spin, color, kappa */

end_of_loops:


    if (this_node == 0)
      printf("RUNNING COMPLETED\n");

    /**DEBUG***/  
#ifdef DEBUGDEF
    light_quark_pion(2) ;
#endif

    for (nk = 0; nk < nkap; nk++)
    {
      if (meascount[nk] > 0)
      {
	if (this_node == 0)
	  printf("total mr iters for measurement= %e\n",
		 (double) avm_iters[nk]);
	if (this_node == 0)
	  printf("average mr iters per spin-color= %e\n",
		 (double) avm_iters[nk] / (double) meascount[nk]);
      }
    }


    endtime=dclock();
    node0_printf("Time = %e seconds\n", (double) (endtime - starttime));

    fflush(stdout);

    /*** calculation section for the variational code *****/
    if( nkap == 1 )
    {
      
      /** sum up the smeared meson correlators over all the nodes ***/
      for(i=0 ; i < nodata ;++i)
      {
	g_complexsum(meson + i) ;
      }
	  

      /* write the smeared correlators to a single disk file ***/
      IF_MASTER
	write_smear_mesonx(meson);
      
      free(meson);  /*** free up the memory for the b-parameter correlators ***/
      
      calc_vary_matrix() ;  /** calculate the static variational matrix **/
      node0_printf(">> The end of the static variational code <<<<\n");
      
    }/** end of the final static variational code *****/


    node0_printf("Time = %e seconds\n",(double)(endtime-starttime));


    fflush(stdout);
    
  } /* end of while(prompt) */

  return 0;

}  /* end of main() */


