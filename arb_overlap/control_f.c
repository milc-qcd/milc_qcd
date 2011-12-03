/***************** control_f.c *****************************************/

/* MIMD version 7 */




#define CONTROL
#include "arb_ov_includes.h"
int read_gauge_info_i(FILE* fp, char *key, int* value);
int read_gauge_info_d(FILE* fp, char *key, double* value);

/* Comment these out if you want to suppress detailed timing */
#define IOTIME
#define PRTIME

int main(int argc, char *argv[])
{
  int prompt;
  Real avm_iters = 0.;
  double dssplaq,dstplaq;
  complex plp;
  int m_iters;

  double starttime,endtime;
#ifdef IOTIME
  double dtime;
#endif
  char tmpfile[MAXFILENAME];
  FILE* fp;


  register int i,j;
  register site *s;

  int mu,k,l;

  Real RsdCG;
  Real cdp,cdm;
  Real cdpp,cdmm;
  int trial_chirality, num_zero;



  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);


  g_sync();
    /* set up */
  prompt = setup_p();



    /* loop over input sets */

  while( readin(prompt) == 0)
    {



      starttime=dclock();

       /* gaugefix if requested:  */
       if( fixflag == COULOMB_GAUGE_FIX){
         if(this_node == 0)
            printf("Fixing to Coulomb gauge\n");
#ifdef IOTIME
          dtime = -dclock();
#endif
          gaugefix(TUP,(Real)1.5,500,(Real)GAUGE_FIX_TOL);


#ifdef IOTIME
          dtime += dclock();
          if(this_node==0)printf("Time to gauge fix = %e\n",dtime);
#endif
        }
       else if( fixflag == LANDAU_GAUGE_FIX){
         if(this_node == 0)
            printf("Fixing to landau gauge\n");
#ifdef IOTIME
          dtime = -dclock();
#endif
          gaugefix(8,(Real)1.5,500,(Real)GAUGE_FIX_TOL);


#ifdef IOTIME
          dtime += dclock();
          if(this_node==0)printf("Time to gauge fix = %e\n",dtime);
#endif
        }


      else
        if(this_node == 0)printf("GAUGE FIXING SKIPPED.\n");

#ifdef NHYP
/*  block_nhyp looks for the thin link in gauge_field_thin[mu],
    and puts the fat link in gauge_field[mu].
    For SF, fixed spatial links at t=0 handled in sf_make_boundary.c
*/
   for(mu=0;mu<4;mu++){
     FORALLSITES(i,s){
         su3mat_copy(&(s->link[mu]), gauge_field_thin[mu]+i);
      }
   }
   
   block_nhyp();
   
   for(mu=0;mu<4;mu++){
     FORALLSITES(i,s){
       su3mat_copy( gauge_field[mu]+i, &(s->link[mu]) );
     }
   }
#endif

/* APBC in time direction     */
   current_boundary=PLUS;
   current_boundary_x=PLUS;
   boundary_flip(MINUS);


   /* and correct the b.c idf you want periodic temporal b.c's */


#ifdef TPERIODIC
	  boundary_flip(PLUS);
	  node0_printf("PERIODIC TEMPORAL b.c.'s !\n\n");
#endif
           /* call plaquette measuring process */
	    plp = ploop();
            d_plaquette(&dssplaq,&dstplaq);
	    m_iters=0;

            if(this_node==0)printf("FATLINK %e %e %e %e %e\n",
                (double)plp.real,(double)plp.imag,(double)m_iters,
                dssplaq,dstplaq);
            /* Re(Polyakov) Im(Poyakov) cg_iters ss_plaq st_plaq */


	  /*
	  FORALLSITES(i,s){
	    for(j=0;j<4;j++){
	      printf("XX %d %d %d %d  : %d\n",s->x,s->y,s->z,s->t,j);
	      dumpmat(&(s->link[j])); 
	    }
	  }
	  */

        /* save lattice if requested */
        if( saveflag != FORGET ) save_lattice( saveflag, savefile, stringLFN  );

      setup_inner();

      if (current_topology<=-100)
      {
         sprintf(tmpfile,"%s.info",startfile);
         fp=fopen(tmpfile,"r");
         if (fp)
         {
             read_gauge_info_i(fp,"gauge.topology",&current_topology);
             read_gauge_info_i(fp,"md.time",&md_time);
         }
         else
             node0_printf("NO .info file : %s\n", tmpfile);
      }


#ifdef FIELD
  /* Allocate space for t_blocked_link 2
printf("off_max %d sites on node %d\n",off_max,sites_on_node);
*/
      t_blocked_link2 =
	(su3_matrix *)malloc(sites_on_node*off_max*sizeof(su3_matrix));

      if(t_blocked_link2==NULL){
	printf("NODE %d: no room for t_blocked_link2\n",this_node);
	terminate(1);
      }


      t_clov =
	(triangular *)malloc(sites_on_node*sizeof(triangular));
      if(t_clov==NULL){
	printf("NODE %d: no room for t_clov\n",this_node);
	terminate(1);
      }

      t_clov_diag =
	(diagonal *)malloc(sites_on_node*sizeof(diagonal));

      if(t_clov_diag==NULL){
	printf("NODE %d: no room for t_clov_diag\n",this_node);
	terminate(1);
      }
      delta0_tmp1=malloc(sites_on_node*sizeof(wilson_vector));
      delta0_htmp1=malloc(sites_on_node*sizeof(half_wilson_vector));
      delta0_htmp2=malloc(sites_on_node*sizeof(half_wilson_vector));
      delta0_htmp3=malloc(sites_on_node*sizeof(half_wilson_vector));
#endif
     

      setup_links(SIMPLE);
      make_clov1();
 
	/* Loop over mass */
      for(k=0;k<num_masses;k++){
	if(mass[k] <= -10.0) exit(1);
	RsdCG=resid[k];
	if(this_node==0)printf("mass= %g r0= %g residue= %g\n",
			       (double)mass[k], (double)wqs.r0, (double)RsdCG);
      }



	/* allocate space for eigenmodes and eigenvals of h(-r0) */

      eigVal0 = (double *)malloc(Nvecs_h0r0*sizeof(double));
      eigVec0 = (wilson_vector **)malloc(Nvecs_h0r0*sizeof(wilson_vector*));
      for(i=0;i<Nvecs_h0r0;i++)
	eigVec0[i]=
	  (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

      if(in_hr0_flag == FRESH){
        build_hr0(0, HIGHP);
      }
      else{
	read_eigen(eigVec0,Nvecs_h0r0,in_hr0_flag,in_hr0);
        build_hr0(1, HIGHP);
      }

      if( out_hr0_flag != FORGET ){
	write_eigen(eigVec0, Nvecs_h0r0, out_hr0_flag, out_hr0);
      }



      /* and eigenmodes and eigenvals of H(0)--if needed */

      /* reset Zolotarov to account for highest kept eigenmode of h(-r0) */
      re_setup_inner(0.9 * fabs(eigVal0[Nvecs_h0r0 - 1]), zolo_max_save);
      for(i=0;i<Norder;i++)
        node0_printf("shift[ %d ]=%e coeff= %e\n",i,shift[i],coeff[i]);
      node0_printf("scale dslash by %e\n",scalez);

      if(Nvecs_hov != 0){
	eigVal = (double *)malloc( (1+Nvecs_hov)*sizeof(double));
	eigVec = (wilson_vector **)malloc(  (1+Nvecs_hov)*sizeof(wilson_vector*));
	for(i=0;i<Nvecs_hov+1;i++)
	  eigVec[i]=
	    (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector)); 

	if(in_hov_flag == FRESH){
	  if(this_node==0) printf("now begin overlap eigenvalues\n");
	  build_h0();
#ifdef KNOWCHI
	  node0_printf("program thinks topology is %d\n",current_topology);
	  trial_chirality=current_topology;
#else
	  avm_iters += build_lowest_chi( &trial_chirality);
#endif
	  avm_iters += build_hov( &trial_chirality,&num_zero);
	}


	/* if you have the vectors in a file, don't make new ones.
	   Recall that the eigenmodes have been ``packed'' with the
	   upper two Dirac components being chirality +, the lower
	   components being the chirality - components. */

	else{
	  read_eigen(eigVec, Nvecs_hov, in_hov_flag, in_hov);
	  read_eigenval(eigVal, Nvecs_hov, in_hov);
	  for(j=0;j<Nvecs_hov;j++) broadcast_double(&(eigVal[j]));

	}



	/* we want to figure out the number of zeromodes if we have not done so already...) */
	  num_zero = 0;
	  for(j=0;j<Nvecs_hov;j++){
	    cdp=0.0;
	    cdm=0.0;
	    FORALLSITES(i,s){
	      cdpp=0.0; 
	      for(l=0;l<2;l++) for(k=0;k<3;k++){
		cdp += cabs_sq(&(eigVec[j][i].d[l].c[k]));
		cdpp += cabs_sq(&(eigVec[j][i].d[l].c[k]));
	      }
	      cdmm=0.0;
	      for(l=2;l<4;l++) for(k=0;k<3;k++){
		cdm += cabs_sq(&(eigVec[j][i].d[l].c[k]));
		cdmm += cabs_sq(&(eigVec[j][i].d[l].c[k]));
	      }
	    }
	    g_floatsum(&cdp);
	    g_floatsum(&cdm);
	    if(cdm< 1.e-6 && cdp >1.e-6){
	      num_zero += 1;
	      node0_printf("TESTNORM %e %e\n",(double)cdm,(double)cdp);}
	    else if (cdm >1.e-6 && cdp < 1.e-6){
	      num_zero -= 1;
	      node0_printf("TESTNORM %e %e\n",(double)cdm,(double)cdp);}
	    else if (cdm >1.e-6 && cdp > 1.e-6){
	      node0_printf("TESTNORM %e %e\n",(double)cdm,(double)cdp);}
	    else{
	      node0_printf("eigVec[%d] is a null vector!\n",j);
	      exit(1);
	    }
	  }
	  if(trial_chirality == 0 && num_zero != 0){
	    trial_chirality = num_zero/abs(num_zero);
	  }

	node0_printf("ZEROMODES %d\n", num_zero);
	/* record current_topology */
	current_topology=num_zero;
	node0_printf("reset current TOPOLOGY %d\n",current_topology);



      /*save overlap vectors if requested */

      if( out_hov_flag != FORGET ){
	write_eigen(eigVec, Nvecs_hov, out_hov_flag, out_hov);
	write_eigenval(eigVal, Nvecs_hov, out_hov);
	}


      if(num_zero != 0){
	/* Set the zero eigenvalues exactly to zero! */
	num_zero = abs(num_zero);
	for(j=0; j<num_zero; j++){
	  eigVal[j] = 0.0;
	}
      }


      } /*if(Nvecs_hov != 0) */



#ifdef INV
      build_params(-R0);
      make_clov1();
      avm_iters += congrad_xxx(F_OFFSET(chi),mass[0], trial_chirality);

#endif


      fflush(stdout);



      for(i=0;i<Nvecs_h0r0;i++) free(eigVec0[i]);
      free(eigVec0);
      free(eigVal0);

      if(Nvecs_hov != 0){
	for(i=0;i<Nvecs_hov;i++) free(eigVec[i]);
	free(eigVec);
	free(eigVal);
      }


	/* end of measurements */

      node0_printf("RUNNING COMPLETED\n");
      node0_printf("avm_iters %e\n",avm_iters);

      endtime=dclock();
      if(this_node==0){
	printf("Time = %e seconds\n",(double)(endtime-starttime));
      }
      fflush(stdout);

#ifdef FIELD
      free(t_blocked_link2);
      free(t_clov);
      free(t_clov_diag);
#endif
    }

  return 0;
}

