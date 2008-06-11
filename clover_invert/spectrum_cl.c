/***************** spectrum_cl.c *****************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

/* These procedures compute the zero momentum meson correlators for a
   "point" sink, rotated sink, and "smeared" sink and the zero
   momentum baryon correlators for a point sink.  Here point means
   both quark and antiquark (or three quarks) are tied together at the
   same sink point.  
   
   "Smeared" means the sink quark fields are convoluted with a
   specified function before tying them together and projecting to
   zero momentum.  NOTE: THIS IS NOT THE SAME AS SMEARING WITH A SINK
   WAVE FUNCTION IN THE RELATIVE COORDINATE, except when the smearing
   function is Gaussian and the width is adjusted appropriately.

*/

#include "cl_inv_includes.h"
#include <string.h>

#define MAX_MESON_PROP 10
#define MAX_BARYON_PROP 4

static char *mes_kind[MAX_MESON_PROP] = {"PION","PS505","PS055","PS0505",
			     "RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};
static char *bar_kind[MAX_BARYON_PROP] = {"PROTON","PROTON0","DELTA","DELTA0"};

static complex *pmes_prop[MAX_MESON_PROP];
static complex *rot_prop[MAX_MESON_PROP];
static complex *smes_prop[MAX_MESON_PROP];
static complex *bar_prop[MAX_BARYON_PROP];

/*--------------------------------------------------------------------*/
void spectrum_cl_init(){
  
  int num_prop;
  int t;
  int key[4];
#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];
  
  
  for(num_prop=0;num_prop<MAX_MESON_PROP;num_prop++){
    pmes_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
    rot_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
    smes_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
    for(t=0;t<nt;t++){
      pmes_prop[num_prop][t] = cmplx(0.0,0.0); 
      rot_prop[num_prop][t] = cmplx(0.0,0.0); 
      smes_prop[num_prop][t] = cmplx(0.0,0.0); 
    }
  }
  
  for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++){
    bar_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
    for(t=0;t<nt;t++)bar_prop[num_prop][t] = cmplx(0.0,0.0); 
  }
  
  /* Set up Fourier transform */
  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;
  setup_restrict_fourier(key, restrict);
  
}

/*--------------------------------------------------------------------*/
void spectrum_cl_baryon(wilson_prop_field qp, complex *bp[]){
  
  double dtime = start_timing();
  w_baryon(qp, qp, qp, bp);
  print_timing(dtime, "diagonal baryons");
}

/*--------------------------------------------------------------------*/
void rotate_prop(spin_wilson_vector *rp, wilson_prop_field qp, int color){
  
  int spin;
  int i;
  site *s;
  wilson_vector *psi, *mp, *tmp;
  
  psi = create_wv_field();
  mp  = create_wv_field();
  tmp = create_wv_field();
  
  /* Construct propagator for "rotated" fields,
     psi_rot = Dslash psi, with Dslash the naive operator. */
  for(spin=0;spin<4;spin++){
    copy_wv_from_wp(psi, qp, color, spin);
    
    /* Do Wilson Dslash on the psi field */
    dslash_w_field(psi, mp,  PLUS, EVENANDODD);
    dslash_w_field(psi, tmp, MINUS, EVENANDODD);
    
    /* From subtraction we get 2*Dslash */
    FORALLSITES(i,s){
      sub_wilson_vector(&(mp[i]), &(tmp[i]), &rp[i].d[spin]);
    }
  }
    
  cleanup_dslash_wtemps();
  destroy_wv_field(psi); 
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}

/*--------------------------------------------------------------------*/
void sink_smear_prop(wilson_prop_field qp){

  int color;
  int ci,si,sf,cf;
  int i;
  site *s;
  wilson_quark_source sink_wqs;
  complex *chi_cs;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark propagator with a Gaussian for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    /* qps points into qp, so qp is changed here */
    restrict_fourier_field((complex *)qps, sizeof(spin_wilson_vector), 
			   FORWARDS);
  }

  print_timing(dtime,"FFT");

  strcpy(sink_wqs.descrp, "gaussian");
  sink_wqs.type   = GAUSSIAN;
  sink_wqs.r0     = sink_r0;
  sink_wqs.x0     = 0;
  sink_wqs.y0     = 0;
  sink_wqs.z0     = 0;
  sink_wqs.t0     = 0;
  w_sink_field(chi_cs, &sink_wqs);
  
  /* We want chi(-k)* -- the complex conjugate of FFT of the
     complex conjugate of the quark sink. */

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function */
  for(ci=0;ci<3;ci++){
    FORALLSITES(i,s)
      for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
	    CMUL(qp[ci][i].d[si].d[sf].c[cf],
		 chi_cs[i],
		 qp[ci][i].d[si].d[sf].c[cf]);
	  }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Note: we do not have to FFT the quark propagator back,
     as long as we don't compute smeared baryons,
     since the meson construction works the same in momentum space!
     We do need to divide by an additional (space) volume
     factor, though! */
  
  free(chi_cs);
}  

/*--------------------------------------------------------------------*/
void spectrum_cl_diag_gen_meson(wilson_prop_field qp, complex *mp[]){

  int color;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  for(color=0;color<3;color++){
    qps = extract_swv_from_wp(qp, color);
    w_meson_field(qps, qps, mp);
  }
  print_timing(dtime, "diagonal mesons");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_diag_meson(wilson_prop_field qp){
  spectrum_cl_diag_gen_meson(qp, pmes_prop);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_diag_rot_meson(wilson_prop_field qp){
  spin_wilson_vector *rps, *qps;
  int color;
  double dtime = start_timing();

  /* Diagonal quark-mass mesons */
  rps = create_swv_field();
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    rotate_prop(rps, qp, color);
    w_meson_field(qps, rps, rot_prop);
  }
  destroy_swv_field(rps);

  print_timing(dtime,"diagonal mesons with rotns\n");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_diag_smeared_meson(wilson_prop_field qp){
    sink_smear_prop(qp);
    spectrum_cl_diag_gen_meson(qp, smes_prop);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_print(int t0, int k){

  Real space_vol;
  Real norm_fac[MAX_MESON_PROP];
  int t, t_sink;
  int num_prop;
  
  space_vol = (Real)(nx*ny*nz);
  for(num_prop=0;num_prop<MAX_MESON_PROP;num_prop++) 
    norm_fac[num_prop] = space_vol;
  norm_fac[4] *= 3.0;
  norm_fac[5] *= 3.0;
  norm_fac[8] *= 3.0;
  norm_fac[9] *= 3.0;
  
  /* print meson propagators */
  for(num_prop=0;num_prop<MAX_MESON_PROP;num_prop++) {
    for(t=0; t<nt; t++){
      t_sink = (t + t0) % nt;
      g_floatsum( &pmes_prop[num_prop][t_sink].real );
      pmes_prop[num_prop][t_sink].real  /= norm_fac[num_prop];
      g_floatsum( &pmes_prop[num_prop][t_sink].imag );
      pmes_prop[num_prop][t_sink].imag  /= norm_fac[num_prop];
      if(this_node == 0)
	printf("POINT%s %d %d  %e %e\n",mes_kind[num_prop],k,t_sink,
	       (double)pmes_prop[num_prop][t_sink].real,
	       (double)pmes_prop[num_prop][t_sink].imag);
    }
    for(t=0; t<nt; t++){
      t_sink = (t + t0) % nt;
      g_floatsum( &rot_prop[num_prop][t_sink].real );
      rot_prop[num_prop][t_sink].real  /= norm_fac[num_prop];
      g_floatsum( &rot_prop[num_prop][t_sink].imag );
      rot_prop[num_prop][t_sink].imag  /= norm_fac[num_prop];
      if(this_node == 0)
	printf("ROT_%s %d %d  %e %e\n",mes_kind[num_prop],k,t_sink,
	       (double)rot_prop[num_prop][t_sink].real,
	       (double)rot_prop[num_prop][t_sink].imag);
    }
    /* Print sink-smeared meson propagator if requested */
    if(strstr(spectrum_request,",sink_smear,") != NULL){
      for(t=0; t<nt; t++){
	t_sink = (t + t0) % nt;
	g_floatsum( &smes_prop[num_prop][t_sink].real );
	smes_prop[num_prop][t_sink].real  /=
	  (space_vol*norm_fac[num_prop]);
	g_floatsum( &smes_prop[num_prop][t_sink].imag );
	smes_prop[num_prop][t_sink].imag  /=
	  (space_vol*norm_fac[num_prop]);
	if(this_node == 0)
	  printf("SMEAR%s %d %d  %e %e\n",mes_kind[num_prop],k,t_sink,
		 (double)smes_prop[num_prop][t_sink].real,
		 (double)smes_prop[num_prop][t_sink].imag);
      }
    }
  }
  
  /* print baryon propagators */
  for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++){
    for(t=0; t<nt; t++){
      t_sink = (t + t0) % nt;
      g_floatsum( &bar_prop[num_prop][t_sink].real );
      bar_prop[num_prop][t_sink].real  /= space_vol;
      g_floatsum( &bar_prop[num_prop][t_sink].imag );
      bar_prop[num_prop][t_sink].imag  /= space_vol;
      if(this_node == 0)
	printf("POINT%s %d %d  %e %e\n",bar_kind[num_prop],k,t_sink,
	       (double)bar_prop[num_prop][t_sink].real,
	       (double)bar_prop[num_prop][t_sink].imag);
    }
  }
}

/*--------------------------------------------------------------------*/
void spectrum_cl_cleanup(){

  int num_prop;

  for(num_prop=0;num_prop<MAX_MESON_PROP;num_prop++){
    free(pmes_prop[num_prop]);
    free(rot_prop[num_prop]);
    free(smes_prop[num_prop]);
  }
  
  for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++)
    free(bar_prop[num_prop]);
}

/*--------------------------------------------------------------------*/
void spectrum_cl(wilson_prop_field qp, int t0, int k){

  spectrum_cl_init();
  spectrum_cl_baryon(qp, bar_prop);

  spectrum_cl_diag_meson(qp);
  spectrum_cl_diag_rot_meson(qp);

  /* Smeared meson propagators */
  if(strstr(spectrum_request,",sink_smear,") != NULL)
    spectrum_cl_diag_smeared_meson(qp);
  
  spectrum_cl_print(t0, k);
  spectrum_cl_cleanup();

}

