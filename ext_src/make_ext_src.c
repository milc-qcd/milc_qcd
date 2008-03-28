/***************** make_prop.c ****************************************/
/* MIMD version 7 */
/* Read and/or generate a propagator */

#include "ext_src_includes.h"

/*--------------------------------------------------------------------*/
static void init_timeslice_fft(int t0){
  int slice[4] = {0,0,0,0};
  int key[4];

  /* Set up 3D Fourier transform for smearing.  Do it only on time slice t0 */
  slice[TUP] = t0;

  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 2;
  setup_restrict_fourier(key, slice);
}

/*--------------------------------------------------------------------*/
static void sink_smear_w_src(wilson_vector *wv, wilson_quark_source *wqs){
  
  int c,d;
  int i;
  site *s;
  complex *chi_cs, z;
  int t0 = wqs->t0;
  double vol = (double)nx*ny*nz;
  double dtime = start_timing();

  /* No smearing for delta function at the origin */
  if(wqs->type == POINT && wqs->x0 == 0 && wqs->y0 == 0 && wqs->z0 == 0)
    return;

  init_timeslice_fft(t0);

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark source with the given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft the quark source on a single time slice (in place) */
  restrict_fourier_field((complex *)wv, sizeof(wilson_vector), FORWARDS);

  print_timing(dtime,"FFT");

  w_sink_field(chi_cs, wqs);
  
  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply the quark by the sink wave function on the selected
     time slice */
  FORALLSITES(i,s)if(s->t == t0){
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

  /* Adjust normalization */
  FORALLSITES(i,s)if(s->t == t0){
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	CMULREAL(wv[i].d[d].c[c], 1./vol, wv[i].d[d].c[c]);
      }
  }

  print_timing(dtime,"FFT");
  free(chi_cs);
}  

/*--------------------------------------------------------------------*/
static void sink_smear_ks_src(su3_vector *v, ks_quark_source *ksqs){
  
  int c;
  int i;
  site *s;
  complex *chi_cs, z;
  int t0 = ksqs->t0;
  double vol = (double)nx*ny*nz;
  double dtime = start_timing();

  /* No smearing for delta function at the origin */
  if(ksqs->type == POINT && ksqs->x0 == 0 && ksqs->y0 == 0 && ksqs->z0 == 0)
    return;

  init_timeslice_fft(t0);

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark source with the given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark source (in place) */
  restrict_fourier_field((complex *)v, sizeof(su3_vector), FORWARDS);

  print_timing(dtime,"FFT");

  ks_sink_field(chi_cs, ksqs);
  
  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function on the selected time slice */
  FORALLSITES(i,s)if(s->t == t0){
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

  /* Adjust normalization */
  FORALLSITES(i,s)if(s->t == t0){
    for(c=0;c<3;c++)
      CMULREAL(v[i].c[c], 1./vol, v[i].c[c]);
  }

  print_timing(dtime,"FFT");
  free(chi_cs);
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
			       wilson_quark_source *dst_wqs,
			       wilson_quark_source *snk_wqs, int snk_gam)
{
  int color, spin;
  int status;
  wilson_vector *wv;
  w_prop_file *fp_in; 
  char *fileinfo;
  wilson_quark_source src_wqs;

  wv = create_wv_field();

  init_wqs(&src_wqs);

  /* Open input Dirac propagator file */
  fp_in  = r_open_wprop(startflag, startfile);

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_wqs);

  /* Open output extended Dirac source file */
  w_open_w_source(dst_wqs, fileinfo);

  /* Loop over source colors and spins */
  for(spin=0;spin<4;spin++)
    for(color=0;color<3;color++){
      
      status = reload_wprop_sc_to_field(startflag, fp_in, &src_wqs, 
					spin, color, wv, 1);
      if(status != 0)terminate(1);

      /* Smear the time slice */
      sink_smear_w_src(wv, snk_wqs);

      /* Multiply by the sink gamma */
      mult_sink_gamma_wv(wv, snk_gam);

      /* Write the extended source as a time slice of the propagator */

      dst_wqs->spin = spin;
      dst_wqs->color = color;
      w_source_write( wv, dst_wqs );
      
    }
  
  /* Close files */
  r_close_wprop(startflag, fp_in);
  w_close_w_source(dst_wqs);

  free_ws_XML(fileinfo);
  destroy_wv_field(wv);

}

/*--------------------------------------------------------------------*/

void extract_ksprop_to_ks_source(int startflag, char startfile[], 
				 ks_quark_source *dst_ksqs,
				 ks_quark_source *snk_ksqs)
{
  int color;
  int status;
  su3_vector *v;
  ks_prop_file *fp_in; 
  char *fileinfo;
  ks_quark_source src_ksqs;

  v = create_v_field();

  init_ksqs(&src_ksqs);

  /* Open files for KS propagators, if requested */
  fp_in  = r_open_ksprop(startflag, startfile);

  /* Create metadata for Dirac source file */
  fileinfo = create_kss_XML(startfile, dst_ksqs);

  /* Open output extended Dirac source file */
  w_open_ks_source(dst_ksqs, fileinfo);

  /* Loop over source colors */

  for(color=0;color<3;color++){
    
    /* Read color vector (and source as appropriate) from file */
    status = reload_ksprop_c_to_field(startflag, fp_in, &src_ksqs, 
				      color, v, 1);
    if(status != 0)terminate(1);

    /* Smear the time slice */
    sink_smear_ks_src(v, snk_ksqs);

    /* Write the extended source as a time slice of the propagator */
    
    dst_ksqs->color = color;
    ks_source_write( v, dst_ksqs );
    
  }
  
  /* close files for staggered propagators */
  r_close_ksprop(startflag, fp_in);
  w_close_ks_source(dst_ksqs);

  free_kss_XML(fileinfo);
  destroy_v_field(v);
}

/*--------------------------------------------------------------------*/
/* Convert staggered propagator to naive.  We assume that the source
   of the staggered propagator has support on hypercube corners       */

void extract_ksprop_to_w_source(int startflag, char startfile[], 
				wilson_quark_source *dst_wqs,
				wilson_quark_source *snk_wqs, int snk_gam)
{
  int color, spin;
  int status;
  su3_vector *v;
  wilson_vector *wv;
  spin_wilson_vector *swv;
  ks_prop_file *fp_in; 
  char *fileinfo;
  int ks_source_r[4] = {0,0,0,0};   /* Hypercube corners */
  ks_quark_source src_ksqs;
  double dtime;

  swv = create_swv_field();
  wv  = create_wv_field();
  v   = create_v_field();

  init_ksqs(&src_ksqs);

  /* Create metadata for Dirac source file */
  fileinfo = create_ws_XML(startfile, dst_wqs);

  /* Open output extended Dirac source file */
  w_open_w_source(dst_wqs, fileinfo);

  /* Loop over source spins (Making spin the outer loop is awkward,
     but it is the standard order for USQCD propagators) */

  for(spin = 0; spin < 4; spin++){

    node0_printf("Spin %d\n", spin);

    /* Open files for KS propagators, if requested */
    fp_in  = r_open_ksprop(startflag, startfile);

    for(color = 0; color < 3; color++){
    
      /* Read color vector (and source as appropriate) from file */
      status = reload_ksprop_c_to_field(startflag, fp_in, &src_ksqs, 
					color, v, 1);
      if(status != 0)terminate(1);
      
      /* Convert KS prop to naive prop (su3_vector maps to
	 spin_wilson_vector) for a given source color */
      
      dtime = start_timing();
      convert_ksprop_to_wprop_swv(swv, v, ks_source_r);
   
      copy_wv_from_swv(wv, swv, spin);
      print_timing(dtime, "naive conversion and copy");
      
      /* Smear the time slice */
      sink_smear_w_src(wv, snk_wqs);

      /* Multiply by the sink gamma */
      mult_sink_gamma_wv(wv, snk_gam);

      /* Write the extended source as a time slice of the propagator */
      dst_wqs->color = color;
      dst_wqs->spin  = spin;
      w_source_write( wv, dst_wqs );
    }

    r_close_ksprop(startflag, fp_in);
  } /* spin */
  
  w_close_w_source(dst_wqs);

  destroy_v_field(v); 
  destroy_wv_field(wv); 
  destroy_swv_field(swv);
}

