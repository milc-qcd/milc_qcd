/***************** make_ext_src.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "ext_src_includes.h"

/*--------------------------------------------------------------------*/
#if 0
static complex *create_w_sink_smearing(quark_source *wqs){
  double vol = (double)nx*ny*nz;
  complex *chi_cs;
  int i;
  site *s;
  double dtime;

  /* No smearing for delta function at the origin */
  if(wqs->type == POINT && wqs->x0 == 0 && wqs->y0 == 0 && wqs->z0 == 0)
    return NULL;

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(chi_cs == NULL){
    printf("create_w_sink_smearing(%d) No room for chi_cs\n",this_node);
    terminate(1);
  }

  dtime = start_timing();
  w_sink_field(chi_cs, wqs);
  print_timing(dtime,"w_sink_field");
  
  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  print_timing(dtime, "FFT of chi and multiply");

  /* Normalize, anticipating inverse FT */
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i], 1./vol, chi_cs[i]);
  }
  
  return chi_cs;
}

/*--------------------------------------------------------------------*/
/* Same as create_w_sink_smearing, but uses ks_sink_field with ksqs */
static complex *create_ks_sink_smearing(quark_source *ksqs){
  double vol = (double)nx*ny*nz;
  complex *chi_cs;
  int i;
  site *s;
  double dtime;

  /* No smearing for delta function at the origin */
  if(ksqs->type == POINT && ksqs->x0 == 0 && ksqs->y0 == 0 && ksqs->z0 == 0)
    return NULL;

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(chi_cs == NULL){
    printf("create_ks_sink_smearing(%d) No room for chi_cs\n",this_node);
    terminate(1);
  }

  dtime = start_timing();
  ks_sink_field(chi_cs, ksqs);
  print_timing(dtime,"ks_sink_field");
  
  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  print_timing(dtime, "FFT of chi and multiply");

  /* Normalize, anticipating inverse FT */
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i], 1./vol, chi_cs[i]);
  }

  return chi_cs;
}

/*--------------------------------------------------------------------*/
static void sink_smear_w_src(wilson_vector *wv, complex *chi_cs,
			     quark_source *wqs){
  int c,d;
  int i;
  site *s;
  complex z;
  double dtime = start_timing();

  /* No smearing for delta function at the origin */
  if(wqs->type == POINT && wqs->x0 == 0 && wqs->y0 == 0 && wqs->z0 == 0)
    return;

  /* Now convolute the quark source with the given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft the quark source on a all time slices (in place) */
  restrict_fourier_field((complex *)wv, sizeof(wilson_vector), FORWARDS);

  print_timing(dtime,"FFT");

  /* Now multiply the quark by the FT sink wave function */
  dtime = start_timing();
  FORALLSITES(i,s){
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	z = wv[i].d[d].c[c];
	CMUL(z, chi_cs[i], wv[i].d[d].c[c]);
      }
  }
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark source (in place) */
  restrict_fourier_field((complex *)wv, sizeof(wilson_vector), BACKWARDS);
  print_timing(dtime,"FFT");

}  

/*--------------------------------------------------------------------*/
static void sink_smear_ks_src(su3_vector *v, complex *chi_cs,
			      quark_source *ksqs){
  
  int c;
  int i;
  site *s;
  complex z;
  double dtime = start_timing();

  /* No smearing for delta function at the origin */
  if(ksqs->type == POINT && ksqs->x0 == 0 && ksqs->y0 == 0 && ksqs->z0 == 0)
    return;

  /* Now convolute the quark source with the given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark source (in place) */
  restrict_fourier_field((complex *)v, sizeof(su3_vector), FORWARDS);

  print_timing(dtime,"FFT");

  /* Now multiply the quark by the FT sink wave function */
  dtime = start_timing();
  FORALLSITES(i,s){
    for(c=0;c<3;c++){
      z = v[i].c[c];
      CMUL(z, chi_cs[i], v[i].c[c]);
    }
  }
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark source (in place) */
  restrict_fourier_field((complex *)v, sizeof(su3_vector), BACKWARDS);
  print_timing(dtime,"FFT");
}  
#endif

/*--------------------------------------------------------------------*/
/* Multiply the sink Wilson vectors by the sink gamma matrix */
static void mult_sink_gamma_wv(wilson_vector *wv, int snk_gam)
{
  site *s;
  int i;
  wilson_vector tmp;

  FORALLSITES(i,s){
    mult_w_by_gamma( wv+i, &tmp, snk_gam );
    *(wv+i) = tmp;
  }
}

/*--------------------------------------------------------------------*/

void extract_wprop_to_w_source(int startflag, char startfile[], int ncolor,
			       int nt0, quark_source *dst_qs,
			       quark_source_sink_op *snk_qs_op, int snk_gam,
			       int dst_type)
{
  int color, spin, ksource;
  int status;
  int j;
  int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
  int r0[4] = {0,0,0,0};            /* Dummy offset */
  wilson_vector *wv;
  //complex *chi_cs;
  w_prop_file *fp_in; 
  char *fileinfo;
  quark_source src_qs;
  double dtime = 0;
#ifdef IOTIME
  int timing = 1;
#else
  int timing = 0;
#endif
  char myname[] = "extract_wprop_to_w_source";

  wv = create_wv_field();

  /* Create sink smearing function */
  // chi_cs = create_w_sink_smearing(snk_qs_op);

  init_qs(&src_qs);
  src_qs.ncolor = ncolor;
  src_qs.nsource = ncolor*4;

  /* Open input Dirac propagator file */
  fp_in  = r_open_wprop(startflag, startfile);

  if(fp_in == NULL){
    node0_printf("%s: Failed to open %s\n", myname, startfile);
    terminate(1);
  }

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_qs);

  /* Open output extended Dirac source file */
  for(j = 0; j < nt0; j++)
    w_source_open_dirac(dst_qs+j, fileinfo);

  /* Loop over source colors and spins */
  for(ksource = 0; ksource < src_qs.nsource; ksource++){
    spin = convert_ksource_to_spin(ksource);
    color = convert_ksource_to_color(ksource);
    
    dtime = start_timing();
    status = reload_wprop_sc_to_field(startflag, fp_in, &src_qs, 
				      spin, color, wv, timing);
    if(status != 0)terminate(1);
    
    /* Smear as requested */
    dtime = start_timing();
    wv_field_op(wv, snk_qs_op, FULL, ALL_T_SLICES);
    
    /* Multiply by the sink gamma */
    dtime = start_timing();
    mult_sink_gamma_wv(wv, snk_gam);
    print_timing(dtime,"mult_sink_gamma_wv");
    
    /* Switch back to staggered basis for KS4 */
    /* This follows the spin convention for the naive extended
       propagator in the clover_invert2 code */
    if(dst_type == KS4_TYPE){
      convert_naive_to_staggered_wv(wv, ks_source_r, r0);
    }

    /* Write the extended source as a time slice of the propagator */
    
    for(j = 0; j < nt0; j++){
      dst_qs[j].ksource = ksource;
      dtime = start_timing();
      w_source_dirac( wv, dst_qs+j );
      print_timing(dtime,"w_source_write");
    }
    
  }
  
  /* Close files */
  r_close_wprop(startflag, fp_in);
  for(j = 0; j < nt0; j++)
    w_source_close(dst_qs+j);

  //  if(chi_cs != NULL)free(chi_cs);
  free_ws_XML(fileinfo);
  destroy_wv_field(wv);

}

/*--------------------------------------------------------------------*/

void extract_ksprop_to_ks_source(int startflag, char startfile[], int ncolor,
				 int nt0, quark_source *dst_qs,
				 quark_source_sink_op *snk_qs_op, 
				 int snk_spin_taste, int r0[])
{
  int color;
  int status;
  int j;
  su3_vector *v,*w;
  //complex *chi_cs;
  ks_prop_file *fp_in; 
  char *fileinfo;
  quark_source src_qs;
  double dtime;
#ifdef IOTIME
  int timing = 1;
#else
  int timing = 0;
#endif

  v = create_v_field();
  w = create_v_field();

  /* Create sink smearing function */
  //chi_cs = create_ks_sink_smearing(snk_qs);

  init_qs(&src_qs);
  src_qs.ncolor = ncolor;

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);

  /* Create metadata for staggered source file */
  fileinfo = create_kss_XML(startfile, dst_qs);

  /* Open output extended color vector source file */
  for(j = 0; j < nt0; j++)
    w_source_open_ks(dst_qs+j, fileinfo);

  /* Loop over source colors */

  for(color = 0;color < ncolor; color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, &src_qs, 
				      color, v, timing);
    if(status != 0)terminate(1);

    /* Smear */
    dtime = start_timing();
    v_field_op(v, snk_qs_op, FULL, ALL_T_SLICES);
    print_timing(dtime,"sink_smear_ks_src");

    /* Apply sink spin-taste operation */
    spin_taste_op_fn(NULL, snk_spin_taste, r0, w, v);
    /* The spin_taste_op phases were designed for tying together two
       forward propagators by first converting one of them to an
       antiquark propagator.  So they include the antiquark phase
       (-)^{x+y+z+t}.  For an extended source we don't want the
       antiquark phase.  The next step removes it.  That way the user
       input snk_spin_taste label describes the actual meson at the
       extended source. */
    spin_taste_op_fn(NULL, spin_taste_index("pion05"), r0, v, w);

    /* Write the extended source as a time slice of the propagator */
    for(j = 0; j < nt0; j++){
      dst_qs[j].color = color;
      dtime = start_timing();
      w_source_ks( v, dst_qs+j );
      print_timing(dtime,"ks_source_write");
    }
    
  }
  
  /* close files for staggered propagators */
  r_close_ksprop(startflag, fp_in);
  for(j = 0; j < nt0; j++)
    w_source_close(dst_qs+j);

  //if(chi_cs != NULL)free(chi_cs);
  free_kss_XML(fileinfo);
  destroy_v_field(w);
  destroy_v_field(v);
}

/*--------------------------------------------------------------------*/
/* Convert staggered propagator to naive.  We assume that the source
   of the staggered propagator has support on hypercube corners       */

void extract_ksprop_to_w_source(int startflag, char startfile[], int ncolor,
				int nt0, quark_source *dst_qs,
				quark_source_sink_op *snk_qs_op, int snk_gam,
				int dst_type)
{
  int color, spin, ksource;
  int status;
  int j;
  su3_vector *v[3];
  wilson_vector *wv;
  spin_wilson_vector *swv;
  //complex *chi_cs;
  ks_prop_file *fp_in; 
  char *fileinfo;
  int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
  int r0[4] = {0,0,0,0};            /* Dummy offset */
  quark_source src_qs;
  double dtime;
#ifdef IOTIME
  int timing = 1;
#else
  int timing = 0;
#endif
  char myname[] = "extract_ksprop_to_w_source";

  swv = create_swv_field();
  wv  = create_wv_field();

  for(color = 0; color < ncolor; color++)
    v[color] = create_v_field();

  /* Create sink smearing function */
  //chi_cs = create_w_sink_smearing(snk_qs_op);

  init_qs(&src_qs);
  src_qs.ncolor = ncolor;

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);
  if(fp_in == NULL){
    node0_printf("%s: Can't read file %s\n", myname, startfile);
    terminate(1);
  }
  
  /* Read the entire staggered propagator */
  for(color = 0; color < ncolor; color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, &src_qs, 
				      color, v[color], timing);
    if(status != 0)terminate(1);
    
  }
  
  r_close_ksprop(startflag, fp_in);

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_qs);

  /* Open output extended Dirac source files */
  for(j = 0; j < nt0; j++)
    w_source_open_dirac(dst_qs+j, fileinfo);
  
  /* Loop over source spins and colors.  This is no longer in the
     standard USQCD order! */
  
  for(ksource = 0; ksource < src_qs.nsource; ksource++){
    spin = convert_ksource_to_spin(ksource);
    color = convert_ksource_to_color(ksource);
    
    /* Convert KS prop to naive prop (su3_vector maps to
       spin_wilson_vector) for a given source color */
    
    dtime = start_timing();
    convert_ksprop_to_wprop_swv(swv, v[color], ks_source_r, r0);

    copy_wv_from_swv(wv, swv, spin);
    print_timing(dtime, "naive conversion and copy");
    
#if 0
    // DEBUG
    {
      wilson_vector wvtmp1, wvtmp2;
      int y[4] = {1,1,1,2};
      int k = node_index(y[0],y[1],y[2],y[3]);
      wvtmp1 = wv[k];
      if(y[3]%2 == 1){mult_w_by_gamma(&wvtmp1, &wvtmp2, GT); wvtmp1 = wvtmp2;}
      if(y[0]%2 == 1){mult_w_by_gamma(&wvtmp1, &wvtmp2, GX); wvtmp1 = wvtmp2;}
      if(y[1]%2 == 1){mult_w_by_gamma(&wvtmp1, &wvtmp2, GY); wvtmp1 = wvtmp2;}
      if(y[2]%2 == 1){mult_w_by_gamma(&wvtmp1, &wvtmp2, GZ); wvtmp1 = wvtmp2;}
      dump_wvec(&wvtmp1);
      dumpvec(&v[color][k]);
    }
#endif

    /* Smear */
    dtime = start_timing();
    wv_field_op(wv, snk_qs_op, FULL, ALL_T_SLICES);
    print_timing(dtime, "sink_smear_w_src");
    
    /* Multiply by the sink gamma */
    dtime = start_timing();
    mult_sink_gamma_wv(wv, snk_gam);
    print_timing(dtime, "mult_sink_gamma");

    /* Switch back to staggered basis for KS4 */
    /* This follows the spin convention for the naive extended
       propagator in the clover_invert2 code */
    if(dst_type == KS4_TYPE){
      convert_naive_to_staggered_wv(wv, ks_source_r, r0);
    }
    
    /* Write the extended source as a time slice of the propagator */
    for(j = 0; j < nt0; j++){
      dst_qs[j].ksource = ksource;
      dtime = start_timing();
      w_source_dirac( wv, dst_qs+j );
      print_timing(dtime, "w_source_write");
    }
  } /* ksource */
  
  for(j = 0; j < nt0; j++)
    w_source_close(dst_qs+j);
  
  //if(chi_cs != NULL)free(chi_cs);
  for(color = 0; color < ncolor; color++)
    destroy_v_field(v[color]); 
  destroy_wv_field(wv); 
  destroy_swv_field(swv);
}

