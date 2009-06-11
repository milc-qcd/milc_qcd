/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "ext_src_includes.h"

/*--------------------------------------------------------------------*/
static complex *create_w_sink_smearing(wilson_quark_source *wqs){
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
static complex *create_ks_sink_smearing(ks_quark_source *ksqs){
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
			     wilson_quark_source *wqs){
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
			      ks_quark_source *ksqs){
  
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

void extract_wprop_to_w_source(int startflag, char startfile[], 
			       int nt0, wilson_quark_source *dst_wqs,
			       wilson_quark_source *snk_wqs, int snk_gam)
{
  int color, spin;
  int status;
  int j;
  wilson_vector *wv;
  complex *chi_cs;
  w_prop_file *fp_in; 
  char *fileinfo;
  wilson_quark_source src_wqs;
  double dtime = 0;

  wv = create_wv_field();

  /* Create sink smearing function */
  chi_cs = create_w_sink_smearing(snk_wqs);

  init_wqs(&src_wqs);

  /* Open input Dirac propagator file */
  fp_in  = r_open_wprop(startflag, startfile);

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_wqs);

  /* Open output extended Dirac source file */
  for(j = 0; j < nt0; j++)
    w_open_w_source(dst_wqs+j, fileinfo);

  /* Loop over source colors and spins */
  for(spin=0;spin<4;spin++)
    for(color=0;color<3;color++){
      
      dtime = start_timing();
      status = reload_wprop_sc_to_field(startflag, fp_in, &src_wqs, 
					spin, color, wv, 1);
      if(status != 0)terminate(1);

      /* Smear as requested */
      dtime = start_timing();
      sink_smear_w_src(wv, chi_cs, snk_wqs);
      print_timing(dtime, "sink_smear_w_src");

      /* Multiply by the sink gamma */
      dtime = start_timing();
      mult_sink_gamma_wv(wv, snk_gam);
      print_timing(dtime,"mult_sink_gamma_wv");

      /* Write the extended source as a time slice of the propagator */

      for(j = 0; j < nt0; j++){
	dst_wqs[j].spin = spin;
	dst_wqs[j].color = color;
	dtime = start_timing();
	w_source_write( wv, dst_wqs+j );
	print_timing(dtime,"w_source_write");
      }
      
    }
  
  /* Close files */
  r_close_wprop(startflag, fp_in);
  for(j = 0; j < nt0; j++)
    w_close_w_source(dst_wqs+j);

  if(chi_cs != NULL)free(chi_cs);
  free_ws_XML(fileinfo);
  destroy_wv_field(wv);

}

/*--------------------------------------------------------------------*/

void extract_ksprop_to_ks_source(int startflag, char startfile[], 
				 int nt0, ks_quark_source *dst_ksqs,
				 ks_quark_source *snk_ksqs)
{
  int color;
  int status;
  int j;
  su3_vector *v;
  complex *chi_cs;
  ks_prop_file *fp_in; 
  char *fileinfo;
  ks_quark_source src_ksqs;
  double dtime;

  v = create_v_field();

  /* Create sink smearing function */
  chi_cs = create_ks_sink_smearing(snk_ksqs);

  init_ksqs(&src_ksqs);

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);

  /* Create metadata for Dirac source file */
  fileinfo = create_kss_XML(startfile, dst_ksqs);

  /* Open output extended Dirac source file */
  for(j = 0; j < nt0; j++)
    w_open_ks_source(dst_ksqs+j, fileinfo);

  /* Loop over source colors */

  for(color=0;color<3;color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, &src_ksqs, 
				      color, v, 1);
    if(status != 0)terminate(1);

    /* Smear */
    dtime = start_timing();
    sink_smear_ks_src(v, chi_cs, snk_ksqs);
    print_timing(dtime,"sink_smear_ks_src");

    /* Write the extended source as a time slice of the propagator */
    for(j = 0; j < nt0; j++){
      dst_ksqs[j].color = color;
      dtime = start_timing();
      ks_source_write( v, dst_ksqs+j );
      print_timing(dtime,"ks_source_write");
    }
    
  }
  
  /* close files for staggered propagators */
  r_close_ksprop(startflag, fp_in);
  for(j = 0; j < nt0; j++)
    w_close_ks_source(dst_ksqs+j);

  if(chi_cs != NULL)free(chi_cs);
  free_kss_XML(fileinfo);
  destroy_v_field(v);
}

/*--------------------------------------------------------------------*/
/* Convert staggered propagator to naive.  We assume that the source
   of the staggered propagator has support on hypercube corners       */

void extract_ksprop_to_w_source(int startflag, char startfile[], 
				int nt0, wilson_quark_source *dst_wqs,
				wilson_quark_source *snk_wqs, int snk_gam)
{
  int color, spin;
  int status;
  int j;
  su3_vector *v[3];
  wilson_vector *wv;
  spin_wilson_vector *swv;
  complex *chi_cs;
  ks_prop_file *fp_in; 
  char *fileinfo;
  int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
  ks_quark_source src_ksqs;
  double dtime;

  swv = create_swv_field();
  wv  = create_wv_field();

  for(color = 0; color < 3; color++)
    v[color]   = create_v_field();

  /* Create sink smearing function */
  chi_cs = create_w_sink_smearing(snk_wqs);

  init_ksqs(&src_ksqs);

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);
  
  /* Read the entire staggered propagator */
  for(color = 0; color < 3; color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, &src_ksqs, 
				      color, v[color], 1);
    if(status != 0)terminate(1);
    
  }
  
  r_close_ksprop(startflag, fp_in);

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_wqs);

  /* Open output extended Dirac source files */
  for(j = 0; j < nt0; j++)
    w_open_w_source(dst_wqs+j, fileinfo);
  
  /* Loop over source spins (Making spin the outer loop is awkward,
     but it is the standard order for USQCD propagators) */
  
  for(spin = 0; spin < 4; spin++){
    
    node0_printf("Spin %d\n", spin);
    
    for(color = 0; color < 3; color++){
      
      /* Convert KS prop to naive prop (su3_vector maps to
	 spin_wilson_vector) for a given source color */
      
      dtime = start_timing();
      convert_ksprop_to_wprop_swv(swv, v[color], ks_source_r);
      
      copy_wv_from_swv(wv, swv, spin);
      print_timing(dtime, "naive conversion and copy");
      
      /* Smear */
      dtime = start_timing();
      sink_smear_w_src(wv, chi_cs, snk_wqs);
      print_timing(dtime, "sink_smear_w_src");

      /* Multiply by the sink gamma */
      dtime = start_timing();
      mult_sink_gamma_wv(wv, snk_gam);
      print_timing(dtime, "mult_sink_gamma");

      /* Write the extended source as a time slice of the propagator */
      for(j = 0; j < nt0; j++){
	dst_wqs[j].color = color;
	dst_wqs[j].spin  = spin;
	dtime = start_timing();
	w_source_write( wv, dst_wqs+j );
	print_timing(dtime, "w_source_write");
      }
    }
  } /* spin */
  
  for(j = 0; j < nt0; j++)
    w_close_w_source(dst_wqs+j);
  
  if(chi_cs != NULL)free(chi_cs);
  for(color = 0; color < 3; color++)
    destroy_v_field(v[color]); 
  destroy_wv_field(wv); 
  destroy_swv_field(swv);
}

