/***************** spectrum_cl_hl.c *************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

#include "cl_inv_includes.h"
#include <string.h>

#define MAX_MESON_PROP_HL 10
#define MAX_BARYON_PROP_HL 6
Real norm_fac[MAX_MESON_PROP_HL];

static char *mes_kind[MAX_MESON_PROP_HL] = {"PION","PS505","PS055","PS0505",
					    "RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};
static char *bar_kind[MAX_BARYON_PROP_HL] = {"PROTON","PROTON0","DELTA","DELTA0",
					     "LAMBDA","LAMBDA0"};

static complex *pmes_prop_hl[MAX_KAP][MAX_KAP][MAX_MESON_PROP_HL];
static complex *rrot_prop_hl[MAX_KAP][MAX_KAP][MAX_MESON_PROP_HL];
static complex *lrot_prop_hl[MAX_KAP][MAX_KAP][MAX_MESON_PROP_HL];
static complex *smes_prop_hl[MAX_KAP][MAX_KAP][MAX_MESON_PROP_HL];
static complex *bar_prop_hl[MAX_KAP][MAX_KAP][MAX_BARYON_PROP_HL];

/*---------------------------------------------------------------------*/
void spectrum_cl_hl_init(){

  int num_prop;
  int i,j,t;
  int key[4];
#define restrict rstrict /* C-90 T3D cludge */
  int restrict[4];


  for(num_prop=0;num_prop<MAX_MESON_PROP_HL;num_prop++)
    for(i=0;i<num_kap;i++)for(j=0;j<=i;j++){
	pmes_prop_hl[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	rrot_prop_hl[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	lrot_prop_hl[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	smes_prop_hl[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	for(t=0;t<nt;t++){
	  pmes_prop_hl[i][j][num_prop][t] = cmplx(0.0,0.0); 
	  rrot_prop_hl[i][j][num_prop][t] = cmplx(0.0,0.0); 
	  lrot_prop_hl[i][j][num_prop][t] = cmplx(0.0,0.0); 
	  smes_prop_hl[i][j][num_prop][t] = cmplx(0.0,0.0); 
	}
      }
  
  for(num_prop=0;num_prop<MAX_BARYON_PROP_HL;num_prop++)
    for(i=0;i<num_kap;i++)for(j=0;j<num_kap;j++){
	bar_prop_hl[i][j][num_prop] = (complex *)malloc(nt*sizeof(complex));
	for(t=0;t<nt;t++)bar_prop_hl[i][j][num_prop][t] = cmplx(0.0,0.0); 
      }
  

  /* Set up Fourier transform */
  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;
  setup_restrict_fourier(key, restrict);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_diag_baryon(wilson_prop_field qp, int k){
  double dtime = start_timing();
  w_baryon(qp, qp, qp, bar_prop_hl[k][k]);
  print_timing(dtime, "diagonal baryons");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_diag_meson(wilson_prop_field qp, int k){
  /* Point meson */
  spectrum_cl_diag_gen_meson(qp, pmes_prop_hl[k][k]);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_diag_rot_meson(wilson_prop_field qp, int k){
  spin_wilson_vector *rqs, *qps;
  int color;
  double dtime = start_timing();

  /* Diagonal quark-mass mesons */
  rqs = create_swv_field();
  for(color = 0; color < 3; color++){
    qps = extract_swv_from_wp(qp, color);
    /* Rotated propagators */
    rotate_prop(rqs, qp, color);
    /* For diagonal propagators, these two are identical */
    w_meson_field(qps, rqs, rrot_prop_hl[k][k]);
    w_meson_field(rqs, qps, lrot_prop_hl[k][k]);
  }
  destroy_swv_field(rqs);

  print_timing(dtime, "diagonal mesons with rotns");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_diag_smeared_meson(wilson_prop_field qp, int k){
  spectrum_cl_diag_gen_meson(qp, smes_prop_hl[k][k]);
}

/*---------------------------------------------------------------------*/
void spectrum_cl_hl_baryon(wilson_prop_field qp1, wilson_prop_field qp2,
			   complex *bpa[], complex *bpb[]){
  w_baryon_hl( qp2, qp1, qp1, bpa);
  w_baryon_hl( qp1, qp2, qp2, bpb);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_offdiag_baryon(wilson_prop_field qp1, 
				   wilson_prop_field qp2, 
				   int j, int k){
  double dtime = start_timing();
  spectrum_cl_hl_baryon(qp1, qp2, bar_prop_hl[j][k], bar_prop_hl[k][j]);
  print_timing(dtime, "hl baryons: one pair of kappas");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_offdiag_gen_meson(wilson_prop_field qp1, 
				      wilson_prop_field qp2, 
				      complex *mp[]){
  spin_wilson_vector *qps1, *qps2;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  for(color = 0; color < 3; color++){
    qps1 = extract_swv_from_wp( qp1, color);
    qps2 = extract_swv_from_wp( qp2, color);

    /* Point meson propagator */
    w_meson_field(qps1, qps2, mp);
  }
  print_timing(dtime, "hl mesons: two kappa pairs");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_offdiag_meson(wilson_prop_field qp1, 
				  wilson_prop_field qp2, 
				  int j, int k){
  /* Point meson */
  spectrum_cl_hl_offdiag_gen_meson(qp1, qp2, pmes_prop_hl[j][k]);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_offdiag_rot_meson(wilson_prop_field qp1, 
				      wilson_prop_field qp2, 
				      int j, int k){
  spin_wilson_vector *rqs, *qps1, *qps2;
  int color;
  double dtime = start_timing();

  /* Off diagonal quark mass mesons */
  rqs = create_swv_field();
  for(color = 0; color < 3; color++){
    qps1 = extract_swv_from_wp(qp1, color);
    qps2 = extract_swv_from_wp(qp2, color);

    /* Corresponding rotated propagators */
    rotate_prop(rqs, qp1, color);
    w_meson_field(rqs, qps2, lrot_prop_hl[j][k]);

    rotate_prop(rqs, qp2, color);
    w_meson_field(qps1, rqs, rrot_prop_hl[j][k]);
  }
  destroy_swv_field(rqs);

  print_timing(dtime, "hl mesons with rotns: two kappa pairs");
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_offdiag_smeared_meson(wilson_prop_field qp1, 
					  wilson_prop_field qp2, 
					  int j, int k){
  /* Both props are FT'd and qp2 has already been smeared */
  spectrum_cl_hl_offdiag_gen_meson(qp1, qp2, smes_prop_hl[j][k]);
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_print(int t0){
  
  Real space_vol;
  Real norm_fac[MAX_MESON_PROP_HL];
  int t, t_sink;
  int num_prop;
  int i, j;
  
  space_vol = (Real)(nx*ny*nz);
  for(num_prop=0;num_prop<MAX_MESON_PROP_HL;num_prop++) 
    norm_fac[num_prop] = space_vol;
  norm_fac[4] *= 3.0;
  norm_fac[5] *= 3.0;
  norm_fac[8] *= 3.0;
  norm_fac[9] *= 3.0;
  
  /* print meson propagators */
  for(num_prop=0;num_prop<MAX_MESON_PROP_HL;num_prop++)
    for(i=0;i<num_kap;i++)for(j=0;j<=i;j++){
      for(t=0; t<nt; t++){
	t_sink = (t + t0) % nt;
	g_floatsum( &pmes_prop_hl[i][j][num_prop][t_sink].real );
	pmes_prop_hl[i][j][num_prop][t_sink].real  /= norm_fac[num_prop];
	g_floatsum( &pmes_prop_hl[i][j][num_prop][t_sink].imag );
	pmes_prop_hl[i][j][num_prop][t_sink].imag  /= norm_fac[num_prop];
	if(this_node == 0)
	  printf("POINT%s %d %d %d  %e %e\n",
		 mes_kind[num_prop],i,j,t_sink,
		 (double)pmes_prop_hl[i][j][num_prop][t_sink].real,
		 (double)pmes_prop_hl[i][j][num_prop][t_sink].imag);
      }
      for(t=0; t<nt; t++){
	t_sink = (t + t0) % nt;
	g_floatsum( &rrot_prop_hl[i][j][num_prop][t_sink].real );
	rrot_prop_hl[i][j][num_prop][t_sink].real  /= norm_fac[num_prop];
	g_floatsum( &rrot_prop_hl[i][j][num_prop][t_sink].imag );
	rrot_prop_hl[i][j][num_prop][t_sink].imag  /= norm_fac[num_prop];
	if(this_node == 0)
	  printf("RROT_%s %d %d %d  %e %e\n",
		 mes_kind[num_prop],i,j,t_sink,
		 (double)rrot_prop_hl[i][j][num_prop][t_sink].real,
		 (double)rrot_prop_hl[i][j][num_prop][t_sink].imag);
      }
      for(t=0; t<nt; t++){
	t_sink = (t + t0) % nt;
	g_floatsum( &lrot_prop_hl[i][j][num_prop][t_sink].real );
	lrot_prop_hl[i][j][num_prop][t_sink].real  /= norm_fac[num_prop];
	g_floatsum( &lrot_prop_hl[i][j][num_prop][t_sink].imag );
	lrot_prop_hl[i][j][num_prop][t_sink].imag  /= norm_fac[num_prop];
	if(this_node == 0)
	  printf("LROT_%s %d %d %d  %e %e\n",
		 mes_kind[num_prop],i,j,t_sink,
		 (double)lrot_prop_hl[i][j][num_prop][t_sink].real,
		 (double)lrot_prop_hl[i][j][num_prop][t_sink].imag);
      }
      for(t=0; t<nt; t++){
	t_sink = (t + t0) % nt;
	g_floatsum( &smes_prop_hl[i][j][num_prop][t_sink].real );
	smes_prop_hl[i][j][num_prop][t_sink].real  /=
	  (space_vol*norm_fac[num_prop]);
	g_floatsum( &smes_prop_hl[i][j][num_prop][t_sink].imag );
	smes_prop_hl[i][j][num_prop][t_sink].imag  /=
	  (space_vol*norm_fac[num_prop]);
	if(this_node == 0)
	  printf("SMEAR%s %d %d %d  %e %e\n",
		 mes_kind[num_prop],i,j,t_sink,
		 (double)smes_prop_hl[i][j][num_prop][t_sink].real,
		 (double)smes_prop_hl[i][j][num_prop][t_sink].imag);
      }
    }
  
  /* print baryon propagators */
  for(num_prop=0;num_prop<6;num_prop++)
    for(i=0;i<num_kap;i++){
      for(j=0;j<i;j++){
	for(t=0; t<nt; t++){
	  t_sink = (t + t0) % nt;
	  g_floatsum( &bar_prop_hl[i][j][num_prop][t_sink].real );
	  bar_prop_hl[i][j][num_prop][t_sink].real  /= space_vol;
	  g_floatsum( &bar_prop_hl[i][j][num_prop][t_sink].imag );
	  bar_prop_hl[i][j][num_prop][t_sink].imag  /= space_vol;
	  if(this_node == 0)
	    printf("POINT%s %d %d %d %d  %e %e\n",
		   bar_kind[num_prop],i,j,j,t_sink,
		   (double)bar_prop_hl[i][j][num_prop][t_sink].real,
		   (double)bar_prop_hl[i][j][num_prop][t_sink].imag);
	}
      }
      if(num_prop<4){
	for(t=0; t<nt; t++){
	  t_sink = (t + t0) % nt;
	  g_floatsum( &bar_prop_hl[i][i][num_prop][t_sink].real );
	  bar_prop_hl[i][i][num_prop][t_sink].real  /= space_vol;
	  g_floatsum( &bar_prop_hl[i][i][num_prop][t_sink].imag );
	  bar_prop_hl[i][i][num_prop][t_sink].imag  /= space_vol;
	  if(this_node == 0)
	    printf("POINT%s %d %d %d %d  %e %e\n",
		   bar_kind[num_prop],i,i,i,t_sink,
		   (double)bar_prop_hl[i][i][num_prop][t_sink].real,
		   (double)bar_prop_hl[i][i][num_prop][t_sink].imag);
	}
      }
      for(j=i+1;j<num_kap;j++){
	for(t=0; t<nt; t++){
	  t_sink = (t + t0) % nt;
	  g_floatsum( &bar_prop_hl[i][j][num_prop][t_sink].real );
	  bar_prop_hl[i][j][num_prop][t_sink].real  /= space_vol;
	  g_floatsum( &bar_prop_hl[i][j][num_prop][t_sink].imag );
	  bar_prop_hl[i][j][num_prop][t_sink].imag  /= space_vol;
	  if(this_node == 0)
	    printf("POINT%s %d %d %d %d  %e %e\n",
		   bar_kind[num_prop],j,j,i,t_sink,
		   (double)bar_prop_hl[i][j][num_prop][t_sink].real,
		   (double)bar_prop_hl[i][j][num_prop][t_sink].imag);
	}
      }
    }
}

/*--------------------------------------------------------------------*/
void spectrum_cl_hl_cleanup(){

  int num_prop;
  int i,j;

  for(num_prop=0;num_prop<MAX_MESON_PROP_HL;num_prop++)
    for(i=0;i<num_kap;i++)for(j=0;j<=i;j++){
      free(pmes_prop_hl[i][j][num_prop]);
      free(rrot_prop_hl[i][j][num_prop]);
      free(lrot_prop_hl[i][j][num_prop]);
      free(smes_prop_hl[i][j][num_prop]);
    }
  
  for(num_prop=0;num_prop<MAX_BARYON_PROP_HL;num_prop++)
    for(i=0;i<num_kap;i++)for(j=0;j<num_kap;j++){
      free(bar_prop_hl[i][j][num_prop]);
    }
}

