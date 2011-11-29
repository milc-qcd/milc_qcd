/* Routines intended for debugging purposes.
   May need MILC_GLOBAL_DEBUG flag to be on */



#include "debug.h"


void d_plaquette_field_minmax(su3_matrix **U_field,
       double *ss_plaq,double *st_plaq,
       double *ss_plaq_min, double *st_plaq_min,
       double *ss_plaq_max, double *st_plaq_max) {
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register site *s;
register int first_pass_s,first_pass_t;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ss_sum,st_sum;
double rtrace_s, rtrace_t, ss_min, st_min, ss_max, st_max;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
    first_pass_s=1;
    first_pass_t=1;

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_field( U_field[dir2], sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_field( U_field[dir1], sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(U_field[dir1][i]);
		m4 = &(U_field[dir2][i]);
		mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP ) {
                  rtrace_t = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  st_sum += rtrace_t;
                }
		else {
                  rtrace_s = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  ss_sum += rtrace_s;
                }
//                printf("Plaq i=%d, dir1=%d, dir2=%d: %f %f\n",
//                  i,dir1,dir2,rtrace_s,rtrace_t);
                /* set min and max values on the first pass */
                if( dir1==TUP ) {
                  if( 1==first_pass_t ) {
                    st_min = rtrace_t;
                    st_max = rtrace_t;
                    first_pass_t = 0;
                  }
                  else {
                    if( rtrace_t < st_min ) st_min = rtrace_t;
                    if( rtrace_t > st_max ) st_max = rtrace_t;
                  }
                }
                else {
                  if( 1==first_pass_s ) {
                    ss_min = rtrace_s;
                    ss_max = rtrace_s;
                    first_pass_s = 0;
                  }
                  else {
                    if( rtrace_s < ss_min ) ss_min = rtrace_s;
                    if( rtrace_s > ss_max ) ss_max = rtrace_s;
                  }
                }
            }

            cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }
    /* DEBUGGING: output min/max values on each node before global operation */
    //printf("MIN PLAQ (W links) on node %d:\t%f\t%f\n", mynode(), ss_min, st_min );
    //printf("MAX PLAQ (W links) on node %d:\t%f\t%f\n", mynode(), ss_max, st_max );

    /* to find min: flip sign, find maximum, flip sign */
    ss_min = -ss_min;
    st_min = -st_min;
    g_doublemax( &ss_min );
    g_doublemax( &st_min );
    g_doublemax( &ss_max );
    g_doublemax( &st_max );

    g_doublesum( &ss_sum );
    g_doublesum( &st_sum );
    *ss_plaq = ss_sum /((Real)(3*nx*ny*nz*nt));
    *st_plaq = st_sum /((double)(3*nx*ny*nz*nt));

    *ss_plaq_min = -ss_min;
    *st_plaq_min = -st_min;
    *ss_plaq_max = ss_max;
    *st_plaq_max = st_max;

    free(su3mat);
} /* d_plaq_field_minmax */



void d_plaquette_field_dump(su3_matrix **U_field, char *file_name_prefix) {
/* su3mat is scratch space of size su3_matrix */
su3_matrix *su3mat;
register int i,dir1,dir2;
register site *s;
register su3_matrix *m1,*m4;
su3_matrix mtmp;
double ss_sum,st_sum;
double rtrace;
msg_tag *mtag0,*mtag1;
    ss_sum = st_sum = 0.0;
FILE *fp;
char plaq_file_name[300];

    su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(su3mat == NULL)
      {
	printf("plaquette: can't malloc su3mat\n");
	fflush(stdout); terminate(1);
      }

    sprintf( plaq_file_name, "%s_node%04d.dat", file_name_prefix, this_node );
    fp = fopen( plaq_file_name, "wt" );


    for(dir1=YUP;dir1<=TUP;dir1++){
	for(dir2=XUP;dir2<dir1;dir2++){

	    mtag0 = start_gather_field( U_field[dir2], sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather_field( U_field[dir1], sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[1] );

	    FORALLSITES(i,s){
		m1 = &(U_field[dir1][i]);
		m4 = &(U_field[dir2][i]);
		mult_su3_an(m4,m1,&su3mat[i]);
	    }

	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_nn( &su3mat[i], (su3_matrix *)(gen_pt[0][i]),
		    &mtmp);

		if(dir1==TUP ) {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  st_sum += rtrace;
                }
		else {
                  rtrace = (double)
		    realtrace_su3((su3_matrix *)(gen_pt[1][i]),&mtmp);
                  ss_sum += rtrace;
                }

                /* output plaquette into file */
                fprintf( fp, "%18.12g\n", rtrace );
            }

            cleanup_gather(mtag0);
	    cleanup_gather(mtag1);
	}
    }


    fclose( fp );

    free(su3mat);
} /* d_plaquette4_field_dump */







#if FERM_ACTION == HISQ

/* Measure plaquettes and write into file */
void g_measure_plaq() {
  double ss_plaq, st_plaq, ss_plaq_min;
  double st_plaq_min, ss_plaq_max, st_plaq_max;
//  double **histogram,**hist_bounds;
//  int ipower, ihist, *Nhist, i;
//  FILE *fp;
//  char hist_file_name[300];
//  char hisq_level1_coeff_file[]="hisq_exp_level1_coeff.dat";
//  double hisq_level1_coeff[5];
//  char temp_str[300];
//  double hist_step, x;


  node0_printf( "Entering g_measure_plaq()\n" );

  /* Load fat and long links for fermion measurements if needed */
  invalidate_all_ferm_links(&fn_links);
#if FERM_ACTION == HISQ
  fn_links.hl.current_X_set = 0;
#endif
  load_ferm_links(&fn_links);

  /* phases out */
  custom_rephase( fn_links.hl.V_link, OFF, &(fn_links.hl.phases_in_V) );
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );

  d_plaquette_field_minmax( fn_links.hl.V_link, &ss_plaq, &st_plaq,
         &ss_plaq_min, &st_plaq_min, &ss_plaq_max, &st_plaq_max );
  if(this_node==0) {
    printf("PLAQ_V:\t%f\t%f\n", ss_plaq, st_plaq );
    printf("PLAQ_V_MIN:\t%f\t%f\n", ss_plaq_min, st_plaq_min );
    printf("PLAQ_V_MAX:\t%f\t%f\n", ss_plaq_max, st_plaq_max );
  }

  d_plaquette_field_minmax( fn_links.hl.W_unitlink, &ss_plaq, &st_plaq,
         &ss_plaq_min, &st_plaq_min, &ss_plaq_max, &st_plaq_max );
  if(this_node==0) {
    printf("PLAQ_W:\t%f\t%f\n", ss_plaq, st_plaq );
    printf("PLAQ_W_MIN:\t%f\t%f\n", ss_plaq_min, st_plaq_min );
    printf("PLAQ_W_MAX:\t%f\t%f\n", ss_plaq_max, st_plaq_max );
  }

  /* phases in */
  custom_rephase( fn_links.hl.V_link, ON, &(fn_links.hl.phases_in_V) );
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );

  node0_printf( "Exiting g_measure_plaq()\n" );
}





/* Measure plaquettes and write into file */
void g_measure_tune() {
//  double ss_plaq, st_plaq, ss_plaq_min;
//  double st_plaq_min, ss_plaq_max, st_plaq_max;
//  double **histogram,**hist_bounds;
//  int ipower, ihist, *Nhist, i;
//  FILE *fp;
//  char hist_file_name[300];
//  char hisq_level1_coeff_file[]="hisq_exp_level1_coeff.dat";
//  double hisq_level1_coeff[5];
//  char temp_str[300];
//  double hist_step, x;


  node0_printf( "Entering g_measure_tune()\n" );

  /* Load fat and long links for fermion measurements if needed */
  load_ferm_links(&fn_links);

  /* phases out U_link */
  custom_rephase( fn_links.hl.U_link, OFF, &(fn_links.hl.phases_in_U) );
  /* dump plaquettes into file */
  d_plaquette_field_dump(fn_links.hl.U_link, "plaq_U_00" );
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.U_link, ON, &(fn_links.hl.phases_in_U) );

  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );
  /* dump plaquettes into file */
  d_plaquette_field_dump(fn_links.hl.W_unitlink, "plaq_W_00" );
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );

  node0_printf( "Exiting g_measure_tune()\n" );
}





void g_measure_hisq_plaq() {
  int i;
  FILE *fp;
  char hisq_level1_coeff_file[]="hisq_exp_level1_coeff.dat";
  double hisq_level1_coeff[5];
  char temp_str[300];
  double fdf_eps=1e-4; /* epsilon for finite difference */


  node0_printf( "Entering g_measure_hisq_plaq()\n" );


  /* read coefficients of path table from file */
  fp = fopen( hisq_level1_coeff_file, "rt" );
  if( NULL==fp ) {
    node0_printf( "File with coefficients %s cannot be opened\n",
                  hisq_level1_coeff_file );
    terminate(1);
  }
  for( i=1; i<5; i++) {
    fgets( temp_str, 300, fp );
    hisq_level1_coeff[i] = atof( temp_str );
  }
  fclose(fp);
  hisq_level1_coeff[0] = 1.0
                        +6*hisq_level1_coeff[1] /* "+" here because of the phase convention */
                       -24*hisq_level1_coeff[2]
                       -48*hisq_level1_coeff[3]
                       -12*hisq_level1_coeff[4];
  if( 0==this_node ) {
    printf( "HISQ EXP coefficients read from file:\n" );
    for( i=0; i<5; i++) {
      printf( " hisq_level1_coeff[%d]=%lf\n", i, hisq_level1_coeff[i] );
    }
    printf( "Precompiled paths coefficients:\n" );
    for( i=0; i<5; i++) {
      printf( " act_path_coeff[%d]=%lf\n", i, (ks_act_paths.p1).act_path_coeff[i] );
    }
  }

  /* LONG COMMENT: we need function (which is 10%-quantile) and its gradient.
     For this purpose finite differences are used. We need variation (fdf_eps)
     with respect to 4 parameters (+/-fdf_eps) and function itself.
     Thus plaquettes on fat links W are evaluated at 9 values of paths coefficients.
     File naming convention: 00 - at present values of coefficients,
     1p - first coefficient is displaced by +fdf_eps,
     1m - first coefficient is displaced by -fdf_eps, etc. */

  /* fdf_block: function */
  /* change path coefficients */
  for( i=0; i<5; i++) {
    (ks_act_paths.p1).act_path_coeff[i] = hisq_level1_coeff[i];
  }
  node0_printf( "fdf_block: function\n" );
  /* Load fat and long links  */
  load_ferm_links(&fn_links);
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );
  /* dump plaquettes into file */
  d_plaquette_field_dump(fn_links.hl.W_unitlink, "plaq_W_00" );
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );

  /* fdf_block: function at 1p */
  /* change path coefficients */
  for( i=1; i<5; i++) {
    (ks_act_paths.p1).act_path_coeff[i] = hisq_level1_coeff[i];
  }
  (ks_act_paths.p1).act_path_coeff[1]+=fdf_eps;
  (ks_act_paths.p1).act_path_coeff[0] = 1.0
                        +6*(ks_act_paths.p1).act_path_coeff[1] /* "+" here because of the phase convention */
                       -24*(ks_act_paths.p1).act_path_coeff[2]
                       -48*(ks_act_paths.p1).act_path_coeff[3]
                       -12*(ks_act_paths.p1).act_path_coeff[4];
  node0_printf( "fdf_block: function at 1p\n" );
  invalidate_all_ferm_links(&fn_links);
  /* Load fat and long links  */
  load_ferm_links(&fn_links);
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );
  /* dump plaquettes into file */
  d_plaquette_field_dump(fn_links.hl.W_unitlink, "plaq_W_1p" );
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );

#ifdef AB_DEBUG
  /* ** measure min/max plaquette on fundamental links ** */
  d_plaquette_minmax( &ss_plaq, &st_plaq,
                      &ss_plaq_min, &st_plaq_min,
                      &ss_plaq_max, &st_plaq_max );
  if(this_node==0) {
    printf("PLAQ:\t%f\t%f\n", ss_plaq, st_plaq );
    printf("MIN PLAQ:\t%f\t%f\n", ss_plaq_min, st_plaq_min );
    printf("MAX PLAQ:\t%f\t%f\n", ss_plaq_max, st_plaq_max );
  }



  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );
  /* ** measure min/max plaquette on fat links ** */
  d_plaquette_field_minmax(fn_links.hl.W_unitlink, &ss_plaq, &st_plaq,
              &ss_plaq_min, &st_plaq_min, &ss_plaq_max, &st_plaq_max);
  if(this_node==0) {
    printf("PLAQ (W):\t%f\t%f\n", ss_plaq, st_plaq );
    printf("MIN PLAQ (W):\t%f\t%f\n", ss_plaq_min, st_plaq_min );
    printf("MAX PLAQ (W):\t%f\t%f\n", ss_plaq_max, st_plaq_max );
  }
  /* phases in W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );

  /* ** measure histograms on fundamental links ** */
  histogram = (double**) malloc(sizeof(double*)*HIST_Npowers);
  hist_bounds = (double**) malloc(sizeof(double*)*HIST_Npowers);
  Nhist=(int*)malloc(sizeof(int)*HIST_Npowers);
  x=6;
  for(ipower=0;ipower<HIST_Npowers;ipower++) {
    histogram[ipower] = (double*) malloc(sizeof(double)*HIST_Nhist);
    hist_bounds[ipower] = (double*) malloc(sizeof(double)*2);
    hist_bounds[ipower][0]=0.0;
    hist_bounds[ipower][1]=x;
    x*=6;
    Nhist[ipower]=50;
  }
  hist_bounds[1][1]=6;
  hist_bounds[2][1]=6;
  d_plaquette_hist( HIST_Npowers, Nhist, histogram, 
                    hist_bounds, &ss_plaq, &st_plaq );
  if( 0==this_node ) { /* output histogram files */
    for( ipower=0; ipower<HIST_Npowers; ipower++ ) {
      hist_step = 
        (hist_bounds[ipower][1]-hist_bounds[ipower][0])/Nhist[ipower];
      /* prepare the file name */
      sprintf( hist_file_name, "hist_plaq_power%02d.dat", ipower+1 );
      fp = fopen( hist_file_name, "wt" );
      for( ihist=0; ihist<Nhist[ipower]; ihist++ ) {
        x = hist_bounds[ipower][0] + (ihist+0.5) * hist_step;
        fprintf( fp, "%18.12g%18.12g\n", x, histogram[ipower][ihist] );
      }
      fclose( fp );
    }
    printf( "Histograms on fundamental links created\n" );
  }
  for(ipower=0;ipower<HIST_Npowers;ipower++) {
    free(histogram[ipower]);
    free(hist_bounds[ipower]);
  }
  free(Nhist);
  free(hist_bounds);
  free(histogram);



  /* ** measure histograms on fat links W ** */
  /* phases out W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, OFF, &(fn_links.hl.phases_in_W) );
  histogram = (double**) malloc(sizeof(double*)*HIST_Npowers);
  hist_bounds = (double**) malloc(sizeof(double*)*HIST_Npowers);
  Nhist=(int*)malloc(sizeof(int)*HIST_Npowers);
  x=6;
  for(ipower=0;ipower<HIST_Npowers;ipower++) {
    histogram[ipower] = (double*) malloc(sizeof(double)*HIST_Nhist);
    hist_bounds[ipower] = (double*) malloc(sizeof(double)*2);
    hist_bounds[ipower][0]=0.0;
    hist_bounds[ipower][1]=x;
    x*=6;
    Nhist[ipower]=50;
  }
  hist_bounds[1][1]=6;
  hist_bounds[2][1]=6;
  d_plaquette_field_hist( fn_links.hl.W_unitlink,
                    HIST_Npowers, Nhist, histogram,
                    hist_bounds, &ss_plaq, &st_plaq );
  if( 0==this_node ) { /* output histogram files */
    for( ipower=0; ipower<HIST_Npowers; ipower++ ) {
      hist_step =
        (hist_bounds[ipower][1]-hist_bounds[ipower][0])/Nhist[ipower];
      /* prepare the file name */
      sprintf( hist_file_name, "hist_plaqW_power%02d.dat", ipower+1 );
      fp = fopen( hist_file_name, "wt" );
      for( ihist=0; ihist<Nhist[ipower]; ihist++ ) {
        x = hist_bounds[ipower][0] + (ihist+0.5) * hist_step;
        fprintf( fp, "%18.12g%18.12g\n", x, histogram[ipower][ihist] );
      }
      fclose( fp );
    }
    printf( "Histograms on fat links created\n" );
  }
  for(ipower=0;ipower<HIST_Npowers;ipower++) {
    free(histogram[ipower]);
    free(hist_bounds[ipower]);
  }
  free(Nhist);
  free(hist_bounds);
  free(histogram);
  /* phases in W_unitlink */
  custom_rephase( fn_links.hl.W_unitlink, ON, &(fn_links.hl.phases_in_W) );
#endif /* AB_DEBUG */





  node0_printf( "Exiting g_measure_hisq_plaq()\n" );
}

#endif /* FERM_ACTION == HISQ */
