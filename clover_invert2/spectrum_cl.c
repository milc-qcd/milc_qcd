/***************** spectrum_cl.c *****************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

/* These procedures compute the zero momentum meson correlators for a
   "point" sink, rotated sink, and "smeared" sink and the zero
   momentum baryon correlators for a point sink.  

   Here point means both quark and antiquark (or three quarks) are
   tied together at the same sink point.
   
   "Smeared" means the the quark field is convoluted with a sink wave
   function before tying it to the antiquark.  Thus the smearing
   function is a true wave function in the relative coordinate.

*/

#include "cl_inv_includes.h"
#include <string.h>
#include <time.h>

#define MAX_BARYON_PROP 4
static char *bar_kind[MAX_BARYON_PROP] = 
  {"PROTON","PROTON0","DELTA","DELTA0"};
static complex *bar_prop[MAX_BARYON_PROP];
#define MAX_BARYON_PROP_OFFDIAG 6
static char *bar_kind_offdiag[MAX_BARYON_PROP_OFFDIAG] = 
  {"PROTON","PROTON0","DELTA","DELTA0", "LAMBDA","LAMBDA0"};
static complex *bar_prop100[MAX_BARYON_PROP];
static complex *bar_prop011[MAX_BARYON_PROP];

static complex ***pmes_prop = NULL;
static complex ***lrot_prop = NULL;
static complex ***rrot_prop = NULL;
static complex ***rot2_prop = NULL;
static complex ***smes_prop = NULL;

/*--------------------------------------------------------------------*/
/* indexing is prop[meson_type][momentum][time] */

static complex *** create_mes_prop(int nmes, int nmom, int ntime){
  complex ***prop;
  int m, p, t;

  prop = (complex ***)malloc(nmes*sizeof(complex **));
  if(prop == NULL)return prop;

  for(m = 0; m < nmes; m++){
    prop[m] = (complex **)malloc(nmom*sizeof(complex *));
    if(prop[m] == NULL)return NULL;
    
    for(p = 0; p < nmom; p++){
      prop[m][p] = (complex *)malloc(ntime*sizeof(complex));
      if(prop[m][p] == NULL)return NULL;

      for(t = 0; t < nt; t++){
	prop[m][p][t].real = 0.0; prop[m][p][t].imag = 0.0;
      }
    }
  }

  return prop;
}

/*--------------------------------------------------------------------*/

static void destroy_mes_prop(complex ***prop){
  int m, p;
  int nmes = param.num_meson_report;
  int nmom = param.num_mom;

  if(prop == NULL)return;

  for(m = 0; m < nmes; m++){
    if(prop[m] != NULL){
      for(p = 0; p < nmom; p++){
	if(prop[m][p] != NULL)free(prop[m][p]);
      }
      free(prop[m]);
    }
  }

  free(prop);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_init(int pair){
  
  int num_mes = param.num_meson_report; /* Number of meson correlators to store */
  int num_mom = param.num_mom;
  int t, num_prop;
  int *qkpair = param.qkpair[pair];
  int key[4];
#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];
  
  if(param.do_point_meson_spect[pair])
    pmes_prop = create_mes_prop(num_mes, num_mom, nt);

  if(param.do_smear_meson_spect[pair])
    smes_prop = create_mes_prop(num_mes, num_mom, nt);

  if(param.do_rot_meson_spect[pair])
    lrot_prop = create_mes_prop(num_mes, num_mom, nt);
  
  if(param.do_rot_meson_spect[pair])
    rrot_prop = create_mes_prop(num_mes, num_mom, nt);

  if(param.do_rot_meson_spect[pair])
    rot2_prop = create_mes_prop(num_mes, num_mom, nt);

  if(qkpair[0] == qkpair[1])
    {
      if(param.do_baryon_spect[pair])
	for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++){
	  bar_prop[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++)
	    bar_prop[num_prop][t] = cmplx(0.0,0.0); 
	}
    } 
  else 
    {

      /* Off-diagonal case */

      if(param.do_baryon_spect[pair])
	for(num_prop=0;num_prop<MAX_BARYON_PROP_OFFDIAG;num_prop++){
	  bar_prop100[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  bar_prop011[num_prop] = (complex *)malloc(nt*sizeof(complex));
	  for(t=0;t<nt;t++){
	    bar_prop100[num_prop][t] = cmplx(0.0,0.0); 
	    bar_prop011[num_prop][t] = cmplx(0.0,0.0); 
	  }
	}
    }
  
  /* Set up Fourier transform for smearing */
  if(param.do_smear_meson_spect[pair]){
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 0;
    setup_restrict_fourier(key, restrict);
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_baryon(wilson_prop_field qp){
  
  double dtime = start_timing();
  w_baryon(qp, qp, qp, bar_prop);
  print_timing(dtime, "diagonal baryons");
}

/*---------------------------------------------------------------------*/
static void spectrum_cl_offdiag_baryon(wilson_prop_field qp0, 
				       wilson_prop_field qp1)
{
  double dtime = start_timing();
  w_baryon_hl( qp1, qp0, qp0, bar_prop100);
  w_baryon_hl( qp0, qp1, qp1, bar_prop011);
  print_timing(dtime, "offdiag baryons");
}

/*--------------------------------------------------------------------*/
static void rotate_prop(spin_wilson_vector *rp, wilson_prop_field qp, 
			int color){
  
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
    //    dslash_w_field(psi, mp,  PLUS, EVENANDODD);
    //    dslash_w_field(psi, tmp, MINUS, EVENANDODD);
    dslash_w_3D_field(psi, mp,  PLUS, EVENANDODD);
    dslash_w_3D_field(psi, tmp, MINUS, EVENANDODD);

    /* From subtraction we get 2*Dslash */
    FORALLSITES(i,s){
      sub_wilson_vector(&(mp[i]), &(tmp[i]), &rp[i].d[spin]);
    }
  }

  cleanup_dslash_w_3D_temps();
  destroy_wv_field(psi); 
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
}

/*--------------------------------------------------------------------*/
static void sink_smear_prop(wilson_prop_field qp, wilson_quark_source *wqs){
  
  int color, spin;
  int ci,si,sf,cf;
  int i;
  site *s;
  complex *chi_cs;
  wilson_vector *chi;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  chi    = create_wv_field();
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

  /* Use chi, for spin=color=0 for the sink wave function */
  FORALLSITES(i,s) clear_wvec( &(chi[i]));

  spin = 0; color = 0;
  wqs->spin   = spin; 
  wqs->color  = color;
  w_sink_field(chi, wqs);
  
  /* Pack for FFT */
  FORALLSITES(i,s){
    chi_cs[i] = chi[i].d[spin].c[color];
  }

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  FORALLSITES(i,s){
    CONJG(chi_cs[i], chi_cs[i]);
  }
  
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

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark_propagator (in place) */
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    /* qps points into qp, so qp is changed here */
    restrict_fourier_field((complex *)qps, sizeof(spin_wilson_vector), 
			   BACKWARDS);
  }
  print_timing(dtime,"FFT");

  destroy_wv_field(chi); 
  free(chi_cs);
}  


/*--------------------------------------------------------------------*/

static void accum_gen_meson(complex ***mp, spin_wilson_vector *qp0, 
			    spin_wilson_vector *qp1){

  meson_cont_mom(mp, qp0, qp1, param.num_mom, param.meson_mom,
		 param.mom_index, param.num_meson, param.meson_index,
		 param.gam_snk, param.gam_src, param.meson_phase);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_gen_meson(wilson_prop_field qp, complex ***mp){
  
  int color;
  spin_wilson_vector *qps;
  double dtime = start_timing();

  for(color=0;color<3;color++){
    qps = extract_swv_from_wp(qp, color);
    accum_gen_meson(mp, qps, qps);
  }
  print_timing(dtime, "diagonal mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_meson(wilson_prop_field qp){
  spectrum_cl_diag_gen_meson(qp, pmes_prop);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_rot_meson(wilson_prop_field qp){
  spin_wilson_vector *rps, *qps;
  int color;
  double dtime = start_timing();

  /* Diagonal quark-mass mesons */
  rps = create_swv_field();
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    rotate_prop(rps, qp, color);
    accum_gen_meson(lrot_prop, rps, qps);
    accum_gen_meson(rrot_prop, qps, rps);
    accum_gen_meson(rot2_prop, rps, rps);
  }
  destroy_swv_field(rps);

  print_timing(dtime,"diagonal mesons with rotns");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_gen_meson(wilson_prop_field qp0, 
					  wilson_prop_field qp1, 
					  complex ***mp){

  spin_wilson_vector *qps0, *qps1;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  for(color = 0; color < 3; color++){
    qps0 = extract_swv_from_wp( qp0, color);
    qps1 = extract_swv_from_wp( qp1, color);

    /* Point meson propagator */
    accum_gen_meson(mp, qps0, qps1);
  }
  print_timing(dtime, "offdiag mesons");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_diag_smeared_meson(wilson_prop_field qp,
					   wilson_quark_source *wqs){
  wilson_prop_field qpcopy = create_wp_field_copy(qp);

  sink_smear_prop(qpcopy, wqs);
  spectrum_cl_offdiag_gen_meson(qp, qpcopy, smes_prop);
  destroy_wp_field(qpcopy);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_meson(wilson_prop_field qp0, 
				      wilson_prop_field qp1){

  /* Point meson */
  spectrum_cl_offdiag_gen_meson(qp0, qp1, pmes_prop);
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_rot_meson(wilson_prop_field qp0, 
					  wilson_prop_field qp1)
{
  spin_wilson_vector *rqs0, *rqs1, *qps0, *qps1;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  rqs0 = create_swv_field();
  rqs1 = create_swv_field();
  for(color = 0; color < 3; color++){
    qps0 = extract_swv_from_wp(qp0, color);
    qps1 = extract_swv_from_wp(qp1, color);

    /* Corresponding rotated propagators */
    rotate_prop(rqs0, qp0, color);
    accum_gen_meson(lrot_prop, rqs0, qps1);

    rotate_prop(rqs1, qp1, color);
    accum_gen_meson(rrot_prop, qps0, rqs1);

    accum_gen_meson(rot2_prop, rqs0, rqs1);
  }
  destroy_swv_field(rqs0);
  destroy_swv_field(rqs1);

  print_timing(dtime, "offdiag mesons with rotns: two kappa pairs");
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_offdiag_smeared_meson(wilson_prop_field qp0, 
					      wilson_prop_field qp1, 
					      wilson_quark_source *wqs){
  wilson_prop_field qpcopy = create_wp_field_copy(qp1);

  sink_smear_prop(qpcopy, wqs);
  spectrum_cl_offdiag_gen_meson(qp0, qpcopy, smes_prop);
  destroy_wp_field(qpcopy);
}

/*--------------------------------------------------------------------*/
static void print_start_meson_prop(int pair, int p, int m, char sink[]){
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("MOMENTUM: %s\n", param.mom_label[p]);
  printf("KAPPAS: %g %g\n",param.dcp[iq0].Kappa,param.dcp[iq1].Kappa);
  printf("SOURCES: %s %s\n",param.src_wqs[iq0].descrp,
	 param.src_wqs[iq1].descrp);
  printf("SINK: %s %s\n",sink, param.meson_label[m]);
}
		       
/*--------------------------------------------------------------------*/

static char *get_utc_datetime(void)
{
  time_t time_stamp;
  struct tm *gmtime_stamp;
  static char time_string[64];

  time(&time_stamp);
  gmtime_stamp = gmtime(&time_stamp);
  strncpy(time_string,asctime(gmtime_stamp),64);
  
  /* Remove trailing end-of-line character */
  if(time_string[strlen(time_string) - 1] == '\n')
    time_string[strlen(time_string) - 1] = '\0';
  return time_string;
}

/*--------------------------------------------------------------------*/
static FILE* open_fnal_meson_file(int pair, char sink[]){
  int iq0 = param.qkpair[pair][0];
  int iq1 = param.qkpair[pair][1];
  FILE *fp;

  /* Create the FNAL file only for rotated meson.  Only node 0
     writes. */
  if(this_node != 0 || param.saveflag_c[pair] == FORGET ||
     ! param.do_rot_meson_spect[pair])
    return NULL;

  fp = fopen(param.savefile_c[pair],"w");
  if(fp == NULL){
    printf("print_open_fnal_meson_file: ERROR. Can't open %s\n",
	   param.savefile_c[pair]);
    return NULL;
  }
  fprintf(fp,"# Job ID:               %s\n",param.job_id);
  fprintf(fp,"# date:                 %s UTC\n",get_utc_datetime());
  fprintf(fp,"# lattice size:         %d x %d x %d x %d\n", nx, ny, nz, nt);
  fprintf(fp,"# spatial volume:       %f\n",((float)nx)*ny*nz);
  fprintf(fp,"# sources: %s %s # sink: %s\n",param.src_wqs[iq0].descrp,
	  param.src_wqs[iq1].descrp,sink);
  fprintf(fp,"# rotated current: d1[%d]=%f  d1[%d]=%f\n",
	  iq0, param.d1[iq0], iq1, param.d1[iq1]);
  fprintf(fp,"# kappas: kappa[%d]=%g  kappa[%d]=%g\n",
	  iq0, param.dcp[iq0].Kappa, iq1, param.dcp[iq1].Kappa);
  return fp;
}
		       
/*--------------------------------------------------------------------*/
static void print_start_fnal_meson_prop(FILE *fp, int pair, int p, int m, 
					char sink[]){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;

  fprintf(fp,"# operator %d # element: %s  momentum: %s\n",
	  m,param.meson_label[m],param.mom_label[p]);
}
		       
/*--------------------------------------------------------------------*/
static void print_start_baryon_prop(int iq1, int iq2, int iq3, char sink[]){
  if(this_node != 0)return;
  printf("STARTPROP\n");
  printf("SOURCES: %s %s %s\n",param.src_wqs[iq1].descrp,
	 param.src_wqs[iq2].descrp, param.src_wqs[iq3].descrp );
  printf("KAPPAS: %g %g %g\n",param.dcp[iq1].Kappa,
	 param.dcp[iq2].Kappa, param.dcp[iq3].Kappa);
  printf("SINK: POINT %s\n",sink);
}
/*--------------------------------------------------------------------*/
static void print_fnal_meson_prop(FILE *fp, int pair, int t, complex c)
{
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fprintf(fp, "%d\t%e\t%e\n", t, (double)c.real, (double)c.imag);
}
/*--------------------------------------------------------------------*/
static void print_end_prop(void){
  node0_printf("ENDPROP\n");
}
/*--------------------------------------------------------------------*/
static void print_end_fnal_meson_prop(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fprintf(fp, "&\n");
}
/*--------------------------------------------------------------------*/
static void close_fnal_meson_file(FILE *fp, int pair){
  if(this_node != 0 || param.saveflag_c[pair] == FORGET)return;
  fclose(fp);
}
/*--------------------------------------------------------------------*/
static void spectrum_cl_print_diag(int pair){

  Real space_vol;
  Real norm_fac;
  FILE *corr_fp;
  int occur_meson[MAX_MESON];
  int occur_mom[MAX_MESON_MOMENTUM];
  int t;
  int i,m,b,p;
  int num_report = param.num_meson_report;
  int iqk;
  complex prop,tmp0, tmp1, tmp2;
  Real x;
  
  /* Count duplicates */
  for(m = 0; m < num_report; m++)
    occur_meson[m] = 0;
  for(i = 0; i < param.num_meson; i++)
    occur_meson[param.meson_index[i]]++;

  /* Count duplicates */
  for(m = 0; m < param.num_mom; m++)
    occur_mom[m] = 0;
  for(i = 0; i < param.num_mom; i++)
    occur_mom[param.mom_index[i]]++;

  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);
  
  /* Initialize FNAL correlator file */
  corr_fp = open_fnal_meson_file(pair, "local sink");

  /* print meson propagators */
  for(p = 0; p < param.num_mom_report; p++)
    for(m=0;m<num_report;m++) {

      norm_fac = occur_meson[m]*occur_mom[p];

      if(param.do_point_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "POINT");
	for(t=0; t<nt; t++){
	  prop = pmes_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	    node0_printf("%d %e %e\n",t,
			 (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }

      /* Rotation - operators only */
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "RIGHT_ROTOP");
	for(t=0; t<nt; t++){
	  prop = rrot_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "LEFT_ROTOP");
	for(t=0; t<nt; t++){
	  prop = lrot_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "BOTH_ROTOP");
	for(t=0; t<nt; t++){
	  prop = rot2_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      /* Apply rotation */
      if(param.do_point_meson_spect[pair] && 
         param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "ROTATED");
	print_start_fnal_meson_prop(corr_fp, pair, p, m, "point rotated");
	for(t=0; t<nt; t++){
	  iqk = param.qkpair[pair][0];
	  /* The MILC sign convention for gamma matrix in Dslash is
	     opposite FNAL, so we rotate with -d1 */
	  x = -param.d1[iqk]/4.;
	  CMULREAL( lrot_prop[m][p][t], x, tmp0);
	  CMULREAL( rrot_prop[m][p][t], x, tmp1);
	  CMULREAL( rot2_prop[m][p][t], x*x, tmp2);
	  CADD( tmp0, tmp1, prop );
	  CSUM( prop, tmp2 );
	  CSUM( prop, pmes_prop[m][p][t] );
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	  print_fnal_meson_prop(corr_fp, pair, t, prop);
	}
	print_end_prop();
	print_end_fnal_meson_prop(corr_fp, pair);
      }

      if(param.do_smear_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, param.snk_label[pair]);
	for(t=0; t<nt; t++){
	  prop = smes_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, space_vol*norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }


    } /* mesons and momenta */
  
  close_fnal_meson_file(corr_fp, pair);
  
  /* print baryon propagators */
  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop();
    }
  }
}

/*--------------------------------------------------------------------*/
static void spectrum_cl_print_offdiag(int pair){
  
  Real space_vol;
  Real norm_fac;
  int occur_meson[MAX_MESON];
  int occur_mom[MAX_MESON_MOMENTUM];
  int t;
  int i,m,p,b;
  int num_report = param.num_meson_report;
  complex prop, tmp0, tmp1, tmp2;
  int iqk0, iqk1;
  FILE *corr_fp;
  Real x0, x1;
  
  /* Count duplicates */
  for(m = 0; m < num_report; m++)
    occur_meson[m] = 0;
  for(i = 0; i < param.num_meson; i++)
    occur_meson[param.meson_index[i]]++;

  /* Count duplicates */
  for(m = 0; m < param.num_mom; m++)
    occur_mom[m] = 0;
  for(i = 0; i < param.num_mom; i++)
    occur_mom[param.mom_index[i]]++;

  /* Normalization factor */
  space_vol = (Real)(nx*ny*nz);

  /* Initialize FNAL correlator file */
  corr_fp = open_fnal_meson_file(pair, "local sink");

  /* print meson propagators */
  for(p = 0; p < param.num_mom_report; p++)
    for(m=0;m<num_report;m++) {
      
      norm_fac = occur_meson[m]*occur_mom[p];

      if(param.do_point_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "POINT");
	for(t=0; t<nt; t++){
	  prop = pmes_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      /* Rotation - operators only */
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "RIGHT_ROTOP");
	for(t=0; t<nt; t++){
	  prop = rrot_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "LEFT_ROTOP");
	for(t=0; t<nt; t++){
	  prop = lrot_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      if(param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "BOTH_ROTOP");
	for(t=0; t<nt; t++){
	  prop = rot2_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
      
      /* Apply rotation */
      if(param.do_point_meson_spect[pair] && 
         param.do_rot_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, "ROTATED");
	print_start_fnal_meson_prop(corr_fp, pair, p, m, "point rotated");
	for(t=0; t<nt; t++){
	  iqk0 = param.qkpair[pair][0];
	  iqk1 = param.qkpair[pair][1];
	  x0 = -param.d1[iqk0]/4.;
	  x1 = -param.d1[iqk1]/4.;
	  CMULREAL( lrot_prop[m][p][t], x0, tmp0);
	  CMULREAL( rrot_prop[m][p][t], x1, tmp1);
	  CMULREAL( rot2_prop[m][p][t], x0*x1, tmp2);
	  CADD( tmp0, tmp1, prop );
	  CSUM( prop, tmp2 );
	  CSUM( prop, pmes_prop[m][p][t] );
	  g_complexsum( &prop );
	  CDIVREAL(prop, norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	  print_fnal_meson_prop(corr_fp, pair, t, prop);
	}
	print_end_prop();
	print_end_fnal_meson_prop(corr_fp, pair);
      }

      if(param.do_smear_meson_spect[pair]){
	print_start_meson_prop(pair, p, m, param.snk_label[pair]);
	for(t=0; t<nt; t++){
	  prop = smes_prop[m][p][t];
	  g_complexsum( &prop );
	  CDIVREAL(prop, space_vol*norm_fac, prop);
	  node0_printf("%d %e %e\n",t,
		       (double)prop.real,(double)prop.imag);
	}
	print_end_prop();
      }
    } /* m < num_report */

  close_fnal_meson_file(corr_fp, pair);

  /* print baryon propagators */

  if(param.do_baryon_spect[pair]){
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][1], param.qkpair[pair][0], 
			      param.qkpair[pair][0], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop100[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop();
    }
    
    for(b=0;b<MAX_BARYON_PROP_OFFDIAG;b++){
      print_start_baryon_prop(param.qkpair[pair][0], param.qkpair[pair][1], 
			      param.qkpair[pair][1], bar_kind_offdiag[b]);
      for(t=0; t<nt; t++){
	prop = bar_prop011[b][t];
	g_complexsum( &prop );
	CDIVREAL(prop, space_vol, prop);
	node0_printf("%d %e %e\n",t,
		     (double)prop.real,(double)prop.imag);
      }
      print_end_prop();
    }
  }
}

/*--------------------------------------------------------------------*/
void spectrum_cl_cleanup(int pair){
  int num_prop;
  int *qkpair = param.qkpair[pair];

  if(param.do_point_meson_spect[pair]){
    destroy_mes_prop(pmes_prop);  pmes_prop = NULL;}
  if(param.do_smear_meson_spect[pair]){
    destroy_mes_prop(smes_prop);  smes_prop = NULL;}
  if(param.do_rot_meson_spect[pair]){
    destroy_mes_prop(rrot_prop); rrot_prop = NULL;
    destroy_mes_prop(lrot_prop); lrot_prop = NULL;
    destroy_mes_prop(rot2_prop); rot2_prop = NULL;
  }

  if(qkpair[0] == qkpair[1]){
    if(param.do_baryon_spect[pair])
      for(num_prop=0;num_prop<MAX_BARYON_PROP;num_prop++)
	free(bar_prop[num_prop]);
  } else {
    if(param.do_baryon_spect[pair])
      for(num_prop=0;num_prop<MAX_BARYON_PROP_OFFDIAG;num_prop++){
	free(bar_prop011[num_prop]);
	free(bar_prop100[num_prop]);
      }
  }
}

/*--------------------------------------------------------------------*/
void spectrum_cl(wilson_prop_field qp0, wilson_prop_field qp1, int pair)
{

  int *qkpair = param.qkpair[pair];

  spectrum_cl_init(pair);

  if(qkpair[0] == qkpair[1]){

    if(param.do_point_meson_spect[pair])
      spectrum_cl_diag_meson(qp0);

    if(param.do_rot_meson_spect[pair])
      spectrum_cl_diag_rot_meson(qp0);

    if(param.do_smear_meson_spect[pair])
      spectrum_cl_diag_smeared_meson(qp0, &param.snk_wqs[pair]);

    if(param.do_baryon_spect[pair])
      spectrum_cl_diag_baryon(qp0);

    spectrum_cl_print_diag(pair);

  } else {

    if(param.do_point_meson_spect[pair])
      spectrum_cl_offdiag_meson(qp0, qp1);

    if(param.do_rot_meson_spect[pair])
      spectrum_cl_offdiag_rot_meson(qp0, qp1);

    if(param.do_smear_meson_spect[pair])
      spectrum_cl_offdiag_smeared_meson(qp0, qp1, &param.snk_wqs[pair]);

    if(param.do_baryon_spect[pair])
      spectrum_cl_offdiag_baryon(qp0, qp1);
    
    spectrum_cl_print_offdiag(pair);
  }    
	    
  spectrum_cl_cleanup(pair);
}
