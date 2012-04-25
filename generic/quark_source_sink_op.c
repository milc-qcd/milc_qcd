/********************** quark_source_sink_op.c **************************/
/* MIMD version 7 */

/* Operations on quark sources and sinks

   Do nothing:

   identity

   Complex field options leading to convolutions:

   complex_field_file              
   complex_field_fm_file 
   gaussian                        
   wavefunction_file               

   Operations exclusive to color vectors:

   funnywall1                 pion5, pioni5, pioni, pions, rhoi, rhos
   funnywall2                 pion05, pionij, pioni0, pion0, rhoi0, rho0

   Covariant modifications of the source:

   covariant_gaussian              Smear with exp(const * Laplacian)
   fat_covariant_gaussian          Same, but base covariant Laplacian on 
                                   APE-smeared links
   deriv1                          Apply covariant derivative
   deriv2_D                        Apply covariant D-type 2nd derivative (symm)
   deriv2_B                        Apply covariant B-type 2nd derivative 
                                   (antisymm)
   deriv3_A                        Apply covariant A0type 3rd derivative 
                                   (not supported)
   hop                             Multiply by hopping matrix
   rotate_3D                       Do 3D FNAL rotation

   General attributes:

   source values can be restricted to corners of elementary hypercubes.
   source values can be scaled by an overall constant.

 */

/* entry points 

   init_qss_op
   create_qss_op
   destroy_qss_op
   copy_qss_op
   copy_qss_op_list
   insert_qss_op

   v_field_op
   wv_field_op
   ksp_sink_op
   wp_sink_op
   get_wv_field_op
   get_v_field_op
   print_field_op_info
   print_field_op_info_list


 */

#include "generic_includes.h"
#include "../include/generic_quark_types.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#endif

/*-------------------------------------------------------------*/
/* Utilities for managing the quark_source_sink_op type */
/*-------------------------------------------------------------*/

void init_qss_op(quark_source_sink_op *qss_op){
  qss_op->type             = UNKNOWN;
  strcpy(qss_op->descrp, "none");
  qss_op->label[0]         = '\0';
  qss_op->a                = 0.;
  qss_op->d1               = 0.;
  qss_op->dir1             = -1;
  qss_op->dir2             = -1;
  qss_op->disp             = 0;
  qss_op->eps_naik         = 0.;
  qss_op->dhop             = 0;
  qss_op->iters            = 0;
  qss_op->r0               = 0.;
  qss_op->r_offset[0]      = 0;
  qss_op->r_offset[1]      = 0;
  qss_op->r_offset[2]      = 0;
  qss_op->r_offset[3]      = 0;
  qss_op->spin_taste       = -1;
  qss_op->source_file[0]   = '\0';
  qss_op->op               = NULL;
} /* init_qss_op */

void set_qss_op_offset(quark_source_sink_op *qss_op, int r0[]){
  qss_op->r_offset[0]      = r0[0];
  qss_op->r_offset[1]      = r0[1];
  qss_op->r_offset[2]      = r0[2];
  qss_op->r_offset[3]      = r0[3];
}

quark_source_sink_op *create_qss_op(void){
  quark_source_sink_op *qss_op;
  
  qss_op = (quark_source_sink_op *)malloc(sizeof(quark_source_sink_op));
  if(qss_op == NULL){
    printf("create_qss_op(%d): No room\n",this_node);
    terminate(1);
  }

  init_qss_op(qss_op);
  return qss_op;
} /* create_qss_op */

/* Destroy the linked list */

void destroy_qss_op(quark_source_sink_op *qss_op){
  quark_source_sink_op *op, *next_op;

  op = qss_op;
  while(op != NULL){
    next_op = op->op;
    free(op);
    op = next_op;
  }
} /* destroy_qss_op */

/* Create a (shallow) copy of the qss_op */

static quark_source_sink_op *copy_qss_op(quark_source_sink_op *src_qss_op){
  quark_source_sink_op *dst_qss_op;

  dst_qss_op = create_qss_op();
  *dst_qss_op = *src_qss_op;
  return dst_qss_op;
} /* copy_qss_op */

/* Copy the qss_op list, starting from the top */

quark_source_sink_op *copy_qss_op_list(quark_source_sink_op *src_qss_op){
  quark_source_sink_op *dst_qss_op, *op;
  
  if(src_qss_op == NULL)return NULL;

  op = dst_qss_op = copy_qss_op(src_qss_op);
  while(op->op != NULL){
    op = op->op = copy_qss_op(op->op);
  }
  return dst_qss_op;
} /* copy_qss_op_list */


/* Broadcast initial data from node 0 to other nodes */
/* We assume that the top-level pointer *qss_op has been
   broadcast from node 0 already */

void broadcast_quark_source_sink_op_recursive(quark_source_sink_op **qss_op){

  if(*qss_op == NULL)return;

  /* All nodes but node 0 make space for their operator data*/
  if(mynode() != 0)
    *qss_op = create_qss_op();

  /* Node 0 broadcasts its data */
  broadcast_bytes((char *)(*qss_op), sizeof(quark_source_sink_op));

  /* Descend the chain of operators */
  broadcast_quark_source_sink_op_recursive(&(*qss_op)->op);
}

/* Add operator to the end of the linked list */

void insert_qss_op(quark_source *qs, quark_source_sink_op *qss_op){
  quark_source_sink_op *op;

  if(qs->op == NULL){
    qs->op = qss_op;
    return;
  }

  op = qs->op;
  while(op->op != NULL){
    op = op->op;
  }
  op->op = qss_op;
} /* insert_qss_op */

/*--------------------------------------------------------------------*/
/* Smearing (convolution) operations                                  */
/*--------------------------------------------------------------------*/

#ifdef HAVE_KS

/* Smear a color vector field */

static void smear_v_field(su3_vector *v, complex *chi_cs){
  
  int cf;
  int i;
  site *s;
  complex z;
  Real x;
  double dtime = start_timing();
  int key[4] = {1,1,1,0};  /* 3D Fourier transform */
  
  /* Set up Fourier transform for smearing */
  setup_restrict_fourier(key, NULL);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			 FORWARDS);
  print_timing(dtime,"FFT");

  dtime = start_timing();

  /* Normalize (for FFT) */
  x = 1./(((Real)nx)*((Real)ny)*(Real)nz);
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i],x,chi_cs[i]);
  }
  
  print_timing(dtime,"source/sink convolution");

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply color vector by sink wave function */
  FORALLSITES(i,s){
    for(cf=0;cf<3;cf++){
      z = v[i].c[cf];
      CMUL(z, chi_cs[i], v[i].c[cf]);
    }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();

  /* fft color vector (in place) */

  restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			 BACKWARDS);

  print_timing(dtime,"FFT");
  cleanup_restrict_fourier();
}   /* smear_v_field */

#endif

#ifdef HAVE_DIRAC

/* Smear a Dirac vector field */

static void smear_wv_field(wilson_vector *wv, complex *chi_cs){
  
  int sf,cf;
  int i;
  site *s;
  complex z;
  Real x;
  double dtime = start_timing();
  int key[4] = {1,1,1,0};  /* 3D Fourier transform */
  
  /* Set up Fourier transform for smearing */
  setup_restrict_fourier(key, NULL);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  restrict_fourier_field((complex *)wv, sizeof(wilson_vector), 
			 FORWARDS);
  print_timing(dtime,"FFT");

  dtime = start_timing();

  /* Normalize (for FFT) */
  x = 1./(((Real)nx)*((Real)ny)*(Real)nz);
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i],x,chi_cs[i]);
  }
  
  print_timing(dtime,"source/sink convolution");

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply wilson vector by sink wave function */
  FORALLSITES(i,s){
    for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
	z = wv[i].d[sf].c[cf];
	CMUL(z, chi_cs[i], wv[i].d[sf].c[cf]);
      }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();

  /* fft wilson vector (in place) */

  restrict_fourier_field((complex *)wv, sizeof(wilson_vector), 
			 BACKWARDS);

  print_timing(dtime,"FFT");
  cleanup_restrict_fourier();
}   /* smear_wv_field */

#endif /* ifdef HAVE_DIRAC */

#if 0
/* Smear a ksprop field */

static void sink_smear_ksprop(ks_prop_field *ksprop, quark_source *qs){
  
  int color;
  int ci,cf;
  int i;
  site *s;
  complex *chi_cs, z;
  Real x;
  su3_vector *v;
  double dtime = start_timing();
  int key[4] = {1,1,1,0};  /* 3D Fourier transform */
  
  /* Set up Fourier transform for smearing */
  setup_restrict_fourier(key, NULL);

  chi_cs = (complex *)malloc(sizeof(complex)*sites_on_node);

  /* Now convolute the quark propagator with a given wave function for
     the smeared mesons. This is done with FFT's */
  
  /* fft quark_propagator (in place) */
  for(color = 0; color < ksprop->nc; color++){
    v = ksprop->v[color];
    /* v points into ksprop, so ksprop is changed here */
    restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			   FORWARDS);
  }

  print_timing(dtime,"FFT");

  dtime = start_timing();

  /* Build sink smearing wave function as a complex field repeated on
     each time slice */
  ks_sink_field(chi_cs, qs);

  /* Normalize (for FFT) */
  x = 1./(((Real)nx)*((Real)ny)*(Real)nz);
  FORALLSITES(i,s){
    CMULREAL(chi_cs[i],x,chi_cs[i]);
  }
  
  print_timing(dtime,"build sink field");

  dtime = start_timing();
  restrict_fourier_field(chi_cs, sizeof(complex), FORWARDS);
  
  /* Now multiply quark by sink wave function */
  for(ci=0;ci<3;ci++){
    FORALLSITES(i,s)
      for(cf=0;cf<ksprop->nc;cf++){
	z = ksprop->v[ci][i].c[cf];
	CMUL(z, chi_cs[i], ksprop->v[ci][i].c[cf]);
      }
  }
  
  print_timing(dtime, "FFT of chi and multiply");

  /* Inverse FFT */
  dtime = start_timing();
  /* fft quark_propagator (in place) */
  for(color = 0; color < ksprop->nc; color++){
    v = ksprop->v[color];
    /* v points into ksprop, so ksprop is changed here */
    restrict_fourier_field((complex *)v, sizeof(su3_vector), 
			   BACKWARDS);
  }
  print_timing(dtime,"FFT");
  cleanup_restrict_fourier();
  free(chi_cs);
}  /* sink_smear_ksprop */

#endif


#ifndef NO_GAUGE_FIELD

/*--------------------------------------------------------------------*/
/* Covariant derivative                                               */
/*--------------------------------------------------------------------*/

#ifdef HAVE_KS

/* On color vector field  */

static void cov_deriv_v(su3_vector *v_dst, su3_vector *v_src,
			su3_matrix *t_links, int dir, int disp,
			Real weights[]){
  int i, j;
  msg_tag *tagm, *tagp;
  site *s;
  su3_vector *vm,*vp,*v;
  char myname[] = "cov_deriv_v";

  /* Check for valid dir */
  if(dir != XUP && dir != YUP && dir != ZUP){
    printf("cov_deriv_v(%d): Bad direction %d\n",this_node,dir);
    terminate(1);
  }

  vp  = create_v_field();
  vm  = create_v_field();
  v   = create_v_field();

  if(t_links == NULL){
    printf("%s(%d): NULL t_links\n",myname,this_node);
    terminate(1);
  }

  /* Zero the destination field */
  clear_v_field(v_dst);
  
  /* vp <- v_src ;  vm <- v_src */
  copy_v_field(vp, v_src);
  copy_v_field(vm, v_src);
  
  /* Parallel transport over distance disp in direction dir */
  for(j = 0; j < disp; j++){
    /* Start gather from positive direction */
    tagp = start_gather_field( vp, sizeof(su3_vector), dir, 
			       EVENANDODD, gen_pt[dir] );
    
    /* Prepare gather from negative direction */
    /* v <- U vm */
    FORALLSITES(i,s){
      mult_adj_su3_mat_vec( &t_links[4*i+dir], vm + i, v + i);
    }
    
    /* Gather from negative direction */
    tagm = start_gather_field(v, sizeof(su3_vector), OPP_DIR(dir), 
			      EVENANDODD, gen_pt[OPP_DIR(dir)] );
    
    /* Wait for gather from negative direction */
    wait_gather(tagm);
    
    /* Copy gather result to vm */
    /* vm <- shift v = shift U vm */
    FORALLSITES(i,s){
      vm[i] = *((su3_vector *)gen_pt[OPP_DIR(dir)][i]);
    }
    
    cleanup_gather(tagm);
    
    /* Complete gather from positive direction and multiply by link mat */
    wait_gather(tagp);
    
    /* v <- U shift vp */
    FORALLSITES(i,s){
      mult_su3_mat_vec( &t_links[4*i+dir], 
			(su3_vector * )(gen_pt[dir][i]), 
			v + i); 
    }
    
    cleanup_gather(tagp);
    
    /* Copy result back to vp */
    /* vp <- U shift vp */
    FORALLSITES(i,s){
      v[i] = vp[i];
    }
    
    /* Take the difference and accumulate with weight to complete the
       derivative */
    /* v_dst <- v_dst + weights[j](vp - vm) */
    
    if(weights[j] != 0){
      FORALLSITES(i,s){
	sub_su3_vector(vp+i, vm+i, v+i);
	scalar_mult_add_su3_vector(v_dst+i, v+i, weights[j], v_dst+i);
      }
    }
  }
  
  destroy_v_field(vp);
  destroy_v_field(vm);
  destroy_v_field(v);
} /* cov_deriv_v */

#endif

#ifdef HAVE_DIRAC

/* On Dirac vector field */

static void cov_deriv_wv(wilson_vector *wv_dst, wilson_vector *wv_src,
			 su3_matrix *t_links, int dir, int disp,
			 Real weights[]){
  int i, j;
  msg_tag *tagm, *tagp;
  site *s;
  wilson_vector *wvm,*wvp,*wv;
  char myname[] = "cov_deriv_wv";

  /* Check for valid dir */
  if(dir != XUP && dir != YUP && dir != ZUP){
    printf("cov_deriv_wv(%d): Bad direction %d\n",this_node,dir);
    terminate(1);
  }

  wvp  = create_wv_field();
  wvm  = create_wv_field();
  wv   = create_wv_field();

  if(t_links == NULL){
    printf("%s(%d): NULL t_links\n",myname,this_node);
    terminate(1);
  }

  /* Zero the destination field */
  clear_wv_field(wv_dst);
  
  /* wvp <- wv_src ;  wvm <- wv_src */
  copy_wv_field(wvp, wv_src);
  copy_wv_field(wvm, wv_src);
  
  /* Parallel transport over distance disp in direction dir */
  for(j = 0; j < disp; j++){
    /* Start gather from positive direction */
    tagp = start_gather_field( wvp, sizeof(wilson_vector), dir, 
			       EVENANDODD, gen_pt[dir] );
    
    /* Prepare gather from negative direction */
    /* wv <- U wvm */
    FORALLSITES(i,s){
      mult_adj_mat_wilson_vec( &t_links[4*i+dir], wvm + i, wv + i);
    }
    
    /* Gather from negative direction */
    tagm = start_gather_field(wv, sizeof(wilson_vector), OPP_DIR(dir), 
			      EVENANDODD, gen_pt[OPP_DIR(dir)] );
    
    /* Wait for gather from negative direction */
    wait_gather(tagm);
    
    /* Copy gather result to wvm */
    /* wvm <- shift wv = shift U wvm */
    FORALLSITES(i,s){
      copy_wvec((wilson_vector *)gen_pt[OPP_DIR(dir)][i], wvm + i);
    }
    
    cleanup_gather(tagm);
    
    /* Complete gather from positive direction and multiply by link mat */
    wait_gather(tagp);
    
    /* wv <- U shift wvp */
    FORALLSITES(i,s){
      mult_mat_wilson_vec( &t_links[4*i+dir], 
			   (wilson_vector * )(gen_pt[dir][i]), 
			   wv + i); 
    }
    
    cleanup_gather(tagp);
    
    /* Copy result back to wvp */
    /* wvp <- U shift wvp */
    FORALLSITES(i,s){
      copy_wvec(wv + i, wvp + i);
    }
    
    /* Take the difference and accumulate with weight to complete the
       derivative */
    /* wv_dst <- wv_dst + weights[j](wvp - wvm) */
    
    if(weights[j] != 0){
      FORALLSITES(i,s){
	sub_wilson_vector(wvp+i, wvm+i, wv+i);
	scalar_mult_add_wvec(wv_dst+i, wv+i, weights[j], wv_dst+i);
      }
    }
  }
  
  destroy_wv_field(wvp);
  destroy_wv_field(wvm);
  destroy_wv_field(wv);
} /* cov_deriv_wv */

#endif /* ifdef HAVE_DIRAC */

#endif /* ifndef NO_GAUGE_FIELD */

#if 0
/* Covariant derivative of Wilson propagator field */
static void cov_deriv_wp(wilson_prop_field *wp_dst, wilson_prop_field *wp_src,
			 su3_matrix *t_links, int dir, int disp, 
			 Real weights[]){
  wilson_vector *wv_src, *wv_dst;
  int spin, color;

  if(wp_dst->nc != wp_src->nc){
    node0_printf("cov_deriv_wp: inconsistent number of colors\n");
    terminate(1);
  }
  wv_src = create_wv_field();
  wv_dst = create_wv_field();
  for(color = 0; color < wp_src->nc; color++)
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(wv_src, wp_src, color, spin);
      cov_deriv_wv(wv_dst, wv_src, t_links, dir, disp, weights);
      copy_wp_from_wv(wp_dst, wv_dst, color, spin);
    }
  
  destroy_wv_field(wv_src);
  destroy_wv_field(wv_dst);
} /* cov_deriv_wp */
#endif

/*--------------------------------------------------------------------*/
/* Covariant Gauss smearing wrapper                                   */
/*--------------------------------------------------------------------*/

#if 0

/* On wilson_prop_field type */

static void gauss_smear_wprop_field(wilson_prop_field *wp, 
				    su3_matrix *t_links, 
				    Real r0, int iters, int t0){
  wilson_vector *wv;
  int spin, color;

  wv = create_wv_field();
  for(color = 0; color < wp->nc; color++)
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(wv, wp, color, spin);
      gauss_smear_wv_field(wv, t_links, r0, iters, t0);
      copy_wp_from_wv(wp, wv, color, spin);
    }
  
  destroy_wv_field(wv);
} /* gauss_smear_wprop_field */

/* On ks_prop_field type */

static void gauss_smear_ksprop_field(ksprop_field *ksp, 
				     su3_matrix *t_links, 
				     Real r0, int iters, int t0){
  su3_vector *v;
  int color;

  v = create_v_field();
  for(color = 0; color < wp->nc; color++){
      copy_v_from_ksp(v, ksp, color);
      gauss_smear_v_field(v, t_links, r0, iters, t0);
      copy_ksp_from_v(ksp, wv, color);
    }
  
  destroy_v_field(v);

} /* guass_smear_ksprop_field */

#endif /* #if 0 */


#ifdef HAVE_DIRAC

#ifndef NO_GAUGE_FIELD

/*--------------------------------------------------------------------*/
/* 3D FNAL rotation                                                   */
/*--------------------------------------------------------------------*/

/* Do the 3D Fermilab rotation on a Wilson vector field.  That is,
   compute src <- (1 + a d1 * \gamma * D/2) src where D is a 3D
   covariant difference normalized to give twice the covariant
   derivative in the continuum limit. */

static void rotate_3D_wvec(wilson_vector *src, Real d1)
{
  int i;
  wilson_vector *mp, *tmp;
  char myname[] = "rotate_3D_wvec";
  
  if(src == NULL){
    node0_printf("%s: Error: called with NULL arg\n", myname);
    terminate(1);
  }

  mp  = create_wv_field();
  tmp = create_wv_field();

  /* Do Wilson Dslash on the source field */
  dslash_w_3D_field(src, mp,  PLUS, EVENANDODD);
  dslash_w_3D_field(src, tmp, MINUS, EVENANDODD);
  
  FORALLFIELDSITES(i){
    /* tmp <- mp - tmp = 2*Dslash*src */
    sub_wilson_vector(mp + i, tmp + i, tmp + i);
    /* src <- d1/4 * tmp + src */
    scalar_mult_add_wvec(src + i, tmp + i, d1/4., src + i);
  }

  cleanup_dslash_w_3D_temps();
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
} /* rotate_3D_wvector */

/*--------------------------------------------------------------------*/
/* Wilson hop                                                         */
/*--------------------------------------------------------------------*/

/* Apply the Wilson hopping matrix or its "derivatives" for fixed mu.
   That is,

   for dhop even, multiply by 

   (1 + gamma_mu) U_x,mu \delta_x,x+mu + (1 - gamma_mu) U^\dagger_(x-mu,mu) 

   and for dhop odd, multiply by

   (1 + gamma_mu) U_x,mu \delta_x,x+mu - (1 - gamma_mu) U^\dagger_(x-mu,mu) 
*/

static void hop_wvec(wilson_vector *src, int dhop, int mu)
{
  int i, sign = 1;
  wilson_vector *mp;
  char myname[] = "hop_wvec";
  
  if(src == NULL){
    node0_printf("%s: Error: called with NULL arg\n", myname);
    terminate(1);
  }

  if(dhop % 2 == 1)sign = -1;

  mp  = create_wv_field();

  /* Apply hopping matrix to the source field */
  hop_w_field(src, mp, PLUS, sign, EVENANDODD, mu);
  
  FORALLFIELDSITES(i){
    copy_wvec(mp + i, src + i);
  }

  destroy_wv_field(mp); 

} /* hop_wvec */

#endif

#endif

#if 0
/*--------------------------------------------------------------------*/
static void rotate_prop_field(wilson_prop_field *dst, wilson_prop_field *src, 
			      quark_source *qs){
  
  int spin, color;
  int i;
  site *s;
  wilson_vector *psi, *mp, *tmp;
  spin_wilson_vector *rp;
  /* The MILC sign convention for gamma matrix in Dslash is
     opposite FNAL, so we rotate with -d1 */
  Real d1 = -qs->d1;
  
  psi = create_wv_field();
  mp  = create_wv_field();
  tmp = create_wv_field();

  for(color = 0; color < src->nc; color++){
    rp = dst->swv[color];

    /* Construct propagator for "rotated" fields,
       psi_rot = Dslash psi, with Dslash the naive operator. */
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(psi, src, color, spin);
      
      /* Do Wilson Dslash on the psi field */
      dslash_w_3D_field(psi, mp,  PLUS, EVENANDODD);
      dslash_w_3D_field(psi, tmp, MINUS, EVENANDODD);
      
      FORALLSITES(i,s){
	/* From subtraction we get 2*Dslash */
	sub_wilson_vector(mp + i, tmp + i, &rp[i].d[spin]);
	/* Apply rotation */
	scalar_mult_add_wvec(psi + i, &rp[i].d[spin], d1/4., &rp[i].d[spin]);
      }
    }
  }
    
  cleanup_dslash_w_3D_temps();
  destroy_wv_field(psi); 
  destroy_wv_field(mp); 
  destroy_wv_field(tmp);
} /* rotate_wprop_field */
#endif

/*--------------------------------------------------------------------*/
/* Time-slice operations on color vector fields                       */
/*--------------------------------------------------------------------*/

#ifdef NO_GAUGE_FIELD
static int requires_gauge_field(int op_type){

  return
    op_type == COVARIANT_GAUSSIAN ||
    op_type == DERIV1 ||
    op_type == DERIV2_D ||
    op_type == DERIV2_B ||
    op_type == DERIV3_A ||
    op_type == FAT_COVARIANT_GAUSSIAN ||
    op_type == FUNNYWALL1 ||
    op_type == FUNNYWALL2 ||
    op_type == ROTATE_3D  ||
    op_type == HOPPING        ||
    op_type == SPIN_TASTE;
}
#endif

static int is_convolution(int op_type){

  return
    op_type == COMPLEX_FIELD_FILE ||
    op_type == COMPLEX_FIELD_FM_FILE ||
    //    op_type == EVENANDODD_WALL ||
    op_type == GAUSSIAN ||
    op_type == WAVEFUNCTION_FILE ;
}

static complex *get_convolving_field(quark_source_sink_op *qss_op, int t0){

  /* Unpack structure. */
  int op_type       = qss_op->type;
  Real a            = qss_op->a;
  Real r0           = qss_op->r0;
  char *sink_file   = qss_op->source_file;
#ifndef HAVE_QIO
  char myname[] = "get_convolving_field";
#endif

  complex *chi_cs   = create_c_field();

  if(op_type == COMPLEX_FIELD_FILE){
    /* Here we read a single complex field */
#ifdef HAVE_QIO
    restore_complex_scidac_to_field(sink_file, QIO_SERIAL, chi_cs, 1);
#else
    node0_printf("%s QIO compilation required for this source\n",myname);
    terminate(1);
#endif
  }
  else if(op_type == COMPLEX_FIELD_FM_FILE){
    r_source_cmplx_fm_to_field(sink_file, chi_cs, 1, 0, 0, 0, t0);
  }
  //  else if(op_type == EVENANDODD_WALL)
  //    even_and_odd_wall(chi_cs, t0);
  
  else if(op_type == GAUSSIAN) {
    gaussian_source(chi_cs, r0, 0, 0, 0, t0);
  }
  else if(op_type == WAVEFUNCTION_FILE){
    fnal_wavefunction(chi_cs, 0, 0, 0, t0, a, sink_file);
  }
  else {
    free(chi_cs);
    return NULL;
  }
  return chi_cs;
}

#ifdef HAVE_KS

static int apply_convolution_v(su3_vector *src, quark_source_sink_op *qss_op,
			     int t0){

  complex *chi_cs   = get_convolving_field(qss_op, t0);

  if(chi_cs == NULL){
    free(chi_cs);
    return 0;
  }

  /* Do the convolution */
  smear_v_field(src, chi_cs);
  free(chi_cs);

  return 1;
}

#endif

#ifdef HAVE_DIRAC

static int apply_convolution_wv(wilson_vector *src, 
			      quark_source_sink_op *qss_op,
			      int t0){

  complex *chi_cs   = get_convolving_field(qss_op, t0);

  if(chi_cs == NULL){
    free(chi_cs);
    return 0;
  }

  /* Do the convolution */

  smear_wv_field(src, chi_cs);
  free(chi_cs);

  return 1;

}

#endif

#ifdef HAVE_KS

static void evenandodd_wall_v_op(su3_vector *src, int t0){

  su3_vector *tslice_sum = (su3_vector *)malloc(sizeof(su3_vector)*nt);
  int i,t;

  for(t = 0; t < nt; t++){
    if(t == t0 || t0 == ALL_T_SLICES)
      clearvec(tslice_sum + t);
  }

  /* Sum local values on requested time slice(s) */
  FORALLFIELDSITES(i){
    t = lattice[i].t;
    if(t == t0 || t0 == ALL_T_SLICES)
      add_su3_vector(tslice_sum + t, src + i, tslice_sum + t);
  }

  /* Global sum */
  g_veccomplexsum( (complex *)tslice_sum, 3*nt);

  /* Copy sums back to field */
  FORALLFIELDSITES(i){
    t = lattice[i].t;
    if(t == t0 || t0 == ALL_T_SLICES)
      src[i] = tslice_sum[t];
  }

  free(tslice_sum);

}

#endif

#ifdef HAVE_DIRAC

static void evenandodd_wall_wv_op(wilson_vector *src, int t0){

  wilson_vector *tslice_sum = (wilson_vector *)malloc(sizeof(wilson_vector)*nt);
  int i,t;

  for(t = 0; t < nt; t++){
    if(t == t0 || t0 == ALL_T_SLICES)
      clear_wvec(tslice_sum + t);
  }

  /* Sum local values on requested time slice(s) */
  FORALLFIELDSITES(i){
    t = lattice[i].t;
    if(t == t0 || t0 == ALL_T_SLICES)
      add_wilson_vector(tslice_sum + t, src + i, tslice_sum + t);
  }

  /* Global sum */
  g_veccomplexsum( (complex *)tslice_sum, 12*nt);

  /* Copy sums back to field */
  FORALLFIELDSITES(i){
    t = lattice[i].t;
    if(t == t0 || t0 == ALL_T_SLICES)
      src[i] = tslice_sum[t];
  }

  free(tslice_sum);
}

#endif

#ifndef NO_GAUGE_FIELD

static int is_cov_deriv(int op_type){

  return
    op_type == DERIV1 ||
    op_type == DERIV2_D ||
    op_type == DERIV2_B ||
    op_type == DERIV2_D ;
}

#ifdef HAVE_KS

static int apply_cov_deriv_v(su3_vector *src, quark_source_sink_op *qss_op){
  int op_type       = qss_op->type;
  int i;
  site *s;

  if(op_type == DERIV1){
    su3_vector *v_dst = create_v_field();
    /* wv_dst <- D_dir1 src */
    cov_deriv_v(v_dst, src, ape_links, qss_op->dir1, qss_op->disp,
		qss_op->weights);
    /* src <- wv_dst */
    copy_v_field(src, v_dst);
    destroy_v_field(v_dst);
  }
  /* Do second derivative if requested */
  else if(op_type == DERIV2_D ||
	  op_type == DERIV2_B ){
    su3_vector *v_dst   = create_v_field();
    su3_vector *v_dst12 = create_v_field();
    su3_vector *v_dst21 = create_v_field();
    
    /* Take 21 derivative */
    /* v_dst <- D_dir1 src */
    cov_deriv_v(v_dst, src, ape_links, qss_op->dir1, qss_op->disp,
		qss_op->weights);
    /* v_dst21 <- D_dir2 D_dir1 src */
    cov_deriv_v(v_dst21, v_dst, ape_links, qss_op->dir2, qss_op->disp,
		qss_op->weights);
    /* Take 12 derivative */
    /* v_dst <- D_dir2 src */
    cov_deriv_v(v_dst, src, ape_links, qss_op->dir2, qss_op->disp,
		qss_op->weights);
    /* v_dst12 <- D_dir1 D_dir2 src */
    cov_deriv_v(v_dst12, v_dst, ape_links, qss_op->dir1, qss_op->disp,
		qss_op->weights);
    /* Make symmetric or antisymmetric tensor combinations --
       result in src */
    if(op_type == DERIV2_D){
      /* Symmetric tensor component */
      /* src <- v_dst12 + v_dst21 */
      FORALLSITES(i,s){
	add_su3_vector(v_dst12+i, v_dst21+i, src+i);
      }
    } else { /* DERIV2_B */
      /* Antisymmetric tensor component */
      /* src <- v_dst12 - v_dst21 */
      FORALLSITES(i,s){
	sub_su3_vector(v_dst12+i, v_dst21+i, src+i);
      }
    }
    
    destroy_v_field(v_dst);
    destroy_v_field(v_dst12);
    destroy_v_field(v_dst21);
  }
  else {
    return 0;
  }

  return 1;
}

#endif /* HAVE_KS */

#ifdef HAVE_DIRAC

static int apply_cov_deriv_wv(wilson_vector *src, quark_source_sink_op *qss_op){
  int i;
  site *s;

  int op_type       = qss_op->type;

  if(op_type == DERIV1){
    wilson_vector *wv_dst = create_wv_field();
    /* wv_dst <- D_dir1 src */
    cov_deriv_wv(wv_dst, src, ape_links, qss_op->dir1, qss_op->disp,
		 qss_op->weights);
    /* src <- wv_dst */
    copy_wv_field(src, wv_dst);
    destroy_wv_field(wv_dst);
  }
  /* Do second derivative if requested */
  else if(op_type == DERIV2_D ||
	  op_type == DERIV2_B ){
    wilson_vector *wv_dst   = create_wv_field();
    wilson_vector *wv_dst12 = create_wv_field();
    wilson_vector *wv_dst21 = create_wv_field();
    
    /* Take 21 derivative */
    /* wv_dst <- D_dir1 src */
    cov_deriv_wv(wv_dst, src, ape_links, qss_op->dir1, qss_op->disp,
		 qss_op->weights);
    /* wv_dst21 <- D_dir2 D_dir1 src */
    cov_deriv_wv(wv_dst21, wv_dst, ape_links, qss_op->dir2, qss_op->disp,
		 qss_op->weights);
    /* Take 12 derivative */
    /* wv_dst <- D_dir2 src */
    cov_deriv_wv(wv_dst, src, ape_links, qss_op->dir2, qss_op->disp,
		 qss_op->weights);
    /* wv_dst12 <- D_dir1 D_dir2 src */
    cov_deriv_wv(wv_dst12, wv_dst, ape_links, qss_op->dir1, qss_op->disp,
		 qss_op->weights);
    /* Make symmetric or antisymmetric tensor combinations --
       result in src */
    if(op_type == DERIV2_D){
      /* Symmetric tensor component */
      /* src <- wv_dst12 + wv_dst21 */
      FORALLSITES(i,s){
	add_wilson_vector(wv_dst12+i, wv_dst21+i, src+i);
      }
    } else { /* DERIV2_B */
      /* Antisymmetric tensor component */
      /* src <- wv_dst12 - wv_dst21 */
      FORALLSITES(i,s){
	sub_wilson_vector(wv_dst12+i, wv_dst21+i, src+i);
      }
    }
    
    destroy_wv_field(wv_dst);
    destroy_wv_field(wv_dst12);
    destroy_wv_field(wv_dst21);
  }

  else {
    return 0;
  }

  return 1;
}

#endif  /* ifdef HAVE_DIRAC */


static int is_cov_smear(int op_type){
  return
    op_type == COVARIANT_GAUSSIAN ||
    op_type == FAT_COVARIANT_GAUSSIAN;
}

#ifdef HAVE_KS


static int apply_cov_smear_v(su3_vector *src, quark_source_sink_op *qss_op,
			     int t0){

  /* Smearing is done with coordinate stride 2 to preserve taste */

  int op_type       = qss_op->type;
  int iters         = qss_op->iters;
  Real r0           = qss_op->r0;

  if(op_type == COVARIANT_GAUSSIAN){
    int iters = qss_op->iters;
    static su3_matrix *t_links;

    t_links = create_G_from_site();
    gauss_smear_v_field(src, t_links, r0, iters, t0);
    destroy_G(t_links);
  }

  else if(op_type == FAT_COVARIANT_GAUSSIAN )
    gauss_smear_v_field(src, ape_links, r0, iters, t0);

  else
    return 0;

  return 1;
}

#endif /* HAVE_KS */

#ifdef HAVE_DIRAC

static int apply_cov_smear_wv(wilson_vector *src, 
			    quark_source_sink_op *qss_op, int t0){

  /* Unpack structure. */
  int op_type       = qss_op->type;
  int iters         = qss_op->iters;
  Real r0           = qss_op->r0;

  if(op_type == COVARIANT_GAUSSIAN){
    int iters = qss_op->iters;

    static su3_matrix *t_links;

    t_links = create_G_from_site();
    gauss_smear_wv_field(src, t_links, r0, iters, t0);
    destroy_G(t_links);
  }

  else if(op_type == FAT_COVARIANT_GAUSSIAN )
    gauss_smear_wv_field(src, ape_links, r0, iters, t0);

  else
    return 0;

  return 1;
}

#endif /* ifdef HAVE_DIRAC */

#ifdef HAVE_KS

static int is_funnywall(int op_type){
  return
    op_type == FUNNYWALL1 ||
    op_type == FUNNYWALL2;
}

static int apply_funnywall(su3_vector *src, quark_source_sink_op *qss_op){
  int op_type = qss_op->type;

  /* (couples to pion5, pioni5, pioni, pions, rhoi, rhos) */
  if(op_type == FUNNYWALL1) {
    su3_vector *tvec1 = create_v_field();
    su3_vector *tvec2 = create_v_field();

    mult_pion5_field( qss_op->r_offset, src, tvec2 );
    mult_pioni_field( ZUP, qss_op->r_offset, src, tvec1, ape_links );
    add_v_fields( tvec2, tvec1, tvec2 );
    mult_rhoi_field( ZUP, qss_op->r_offset, src, tvec1 );
    add_v_fields( src, tvec1, tvec2 );
    
    destroy_v_field(tvec2);
    destroy_v_field(tvec1);
  }

  /* (couples to pion05, pionij, pioni0, pion0, rhoi0, rho0) */
  else if(op_type == FUNNYWALL2) {
    su3_vector *tvec1 = create_v_field();
    su3_vector *tvec2 = create_v_field();

    mult_pion05_field( qss_op->r_offset, src, tvec2 );
    mult_pioni0_field( ZUP, qss_op->r_offset, src, tvec1, ape_links );
    add_v_fields( tvec2, tvec1, tvec2 );
    mult_rhoi0_field( ZUP, qss_op->r_offset, src, tvec1 );
    add_v_fields( src, tvec1, tvec2 );
    
    destroy_v_field(tvec2);
    destroy_v_field(tvec1);
  }
  else {
    return 0;
  }

  return 1;
}

static void apply_spin_taste(su3_vector *src, quark_source_sink_op *qss_op){
  su3_vector *dst = create_v_field();
  
  /* At present we do not support "fn" type spin-taste operators here */
  spin_taste_op(qss_op->spin_taste, qss_op->r_offset, dst, src);
  
  copy_v_field(src, dst);
  destroy_v_field(dst);
  
}

/*--------------------------------------------------------------------*/
/* KS hop                                                             */
/*--------------------------------------------------------------------*/

/* Apply the staggered hopping matrix or its "derivatives" for fixed
   direction mu

   That is, for dhop = 0, multiply by 
   \alpha_x,mu (D_x,mu \delta_x,x+mu - D^\dagger_(x-mu,mu)) +
   \alpha_x,mu (D3_x,mu \delta_x,x+3mu - D3^\dagger_(x-3mu,mu)) 

   For dhop = 1, multiply by 
   \alpha_x,mu (D_x,mu \delta_x,x+mu + D^\dagger_(x-mu,mu)) +
   3*\alpha_x,mu (D3_x,mu \delta_x,x+3mu + D3^\dagger_(x-3mu,mu)) 

   For dhop = 2, multiply by 
   \alpha_x,mu (D_x,mu \delta_x,x+mu - D^\dagger_(x-mu,mu)) +
   9*\alpha_x,mu (D3_x,mu \delta_x,x+3mu - D3^\dagger_(x-3mu,mu)) 

   etc.

   where D are the fat links, D3 are the long links, and alpha_x,mu are
   the staggered phases.
*/

static void hop_vec(su3_vector *src, Real eps, int dhop, int mu)
{
  int sign;
  su3_vector *v;
  char myname[] = "hop_vec";
  Real wtfatf, wtfatb, wtlongf, wtlongb;
#if FERM_ACTION == HISQ
  int n_naiks = get_n_naiks_hisq(fn_links);
  double *eps_naik = get_eps_naik_hisq(fn_links);
  int inaik = index_eps_naik(eps_naik, n_naiks, eps);
#else
  int inaik = 0;
#endif
  imp_ferm_links_t *fn = get_fm_links(fn_links)[inaik];
  
  if(src == NULL){
    node0_printf("%s: Error: called with NULL arg\n", myname);
    terminate(1);
  }

  if(dhop % 2 == 0)
    sign = +1;
  else
    sign = -1;

  wtfatf = 1.;
  wtlongf = pow(3.,dhop);
  wtfatb = sign*wtfatf;
  wtlongb = sign*wtlongf;

  v = create_v_field();  /* Create and clear */

  /* v += hop*src */
  dslash_fn_dir(src, v, EVENANDODD, fn, mu, +1, wtfatf, wtlongf); /* forward */
  dslash_fn_dir(src, v, EVENANDODD, fn, mu, -1, wtfatb, wtlongb); /* backward */
  
  /* result in src */
  copy_v_field(src, v);

  destroy_v_field(v); 

} /* hop_vec */

#endif /* HAVE_KS */

#endif /* ifndef NO_GAUGE_FIELD */

#ifdef HAVE_KS

void v_field_op(su3_vector *src, quark_source_sink_op *qss_op, 
		int subset, int t0)
{
  char myname[] = "v_field_op";
  
  /* Unpack structure. */
  int op_type       = qss_op->type;

  if(op_type == IDENTITY)
    return;

  if(op_type == EVENANDODD_WALL)
    evenandodd_wall_v_op(src, t0);

  /* Convolution operations */

  else if(is_convolution(op_type))
    apply_convolution_v(src, qss_op, t0);

#ifndef NO_GAUGE_FIELD

  else if(is_cov_deriv(op_type))
    apply_cov_deriv_v(src, qss_op);

  else if(is_cov_smear(op_type))
    apply_cov_smear_v(src, qss_op, t0);

  else if(is_funnywall(op_type))
    apply_funnywall(src, qss_op);

  else if(op_type == SPIN_TASTE)
    apply_spin_taste(src, qss_op);

  else if(op_type == HOPPING)
    hop_vec(src, qss_op->eps_naik, qss_op->dhop, qss_op->dir1);

#endif

  else {
    node0_printf("%s: Unrecognized operator type %d\n", myname, op_type);
    terminate(1);
  }

  /* Apply subset mask */

  subset_mask_v(src, subset, t0);

} /* v_field_op */

#endif

#ifdef HAVE_DIRAC

/*--------------------------------------------------------------------*/
/* Time-slice operations on Dirac vector fields                       */
/*--------------------------------------------------------------------*/

void wv_field_op(wilson_vector *src, quark_source_sink_op *qss_op, 
		 int subset, int t0)
{
  char myname[] = "wv_field_op";
  
  /* Unpack structure. */
  int op_type       = qss_op->type;

  if(op_type == IDENTITY)
    return;

  if(op_type == EVENANDODD_WALL)
    evenandodd_wall_wv_op(src, t0);

  /* Convolution operations */

  else if(is_convolution(op_type))
    apply_convolution_wv(src, qss_op, t0);

#ifndef NO_GAUGE_FIELD
  else if(is_cov_smear(op_type))
    apply_cov_smear_wv(src, qss_op, t0);

  else if(is_cov_deriv(op_type))
    apply_cov_deriv_wv(src, qss_op);

  else if(op_type == HOPPING)
    hop_wvec(src, qss_op->dhop, qss_op->dir1);

  else if(op_type == ROTATE_3D)
    /* The MILC sign convention for gamma matrix in Dslash is
       opposite FNAL, so we rotate with -d1 */
    rotate_3D_wvec(src, -qss_op->d1);

#endif

  else {
    node0_printf("%s: Unrecognized operator type %d\n", myname, op_type);
    terminate(1);
  }

  /* Apply subset mask */

  subset_mask_wv(src, subset, t0);

} /* wv_field_op */

#endif

#ifdef HAVE_KS
/*--------------------------------------------------------------------*/
/* Sink operations on the ks_prop_field type                          */
/*--------------------------------------------------------------------*/

void ksp_sink_op(quark_source_sink_op *qss_op, ks_prop_field *ksp )
{
  int color;
  su3_vector *v = create_v_field();

  for(color = 0; color < ksp->nc; color++){
      copy_v_from_ksp(v, ksp, color);
      v_field_op( v, qss_op, FULL, ALL_T_SLICES);
      insert_ksp_from_v(ksp, v, color);
  }
  
  destroy_v_field(v);
} /* ksp_sink_op */

#endif

#ifdef HAVE_DIRAC

/*--------------------------------------------------------------------*/
/* Sink operations on the wilson_prop_field type                      */
/*--------------------------------------------------------------------*/

void wp_sink_op(quark_source_sink_op *qss_op, wilson_prop_field *wp )
{
  int color, spin;
  wilson_vector *wv = create_wv_field();

  for(color = 0; color < wp->nc; color++)
    for(spin = 0; spin < 4; spin++){
      copy_wv_from_wp(wv, wp, color, spin);
      wv_field_op( wv, qss_op, FULL, ALL_T_SLICES);
      copy_wp_from_wv(wp, wv, color, spin);
    }

  destroy_wv_field(wv);
} /* wp_sink_op */

#endif

/*--------------------------------------------------------------------*/
/* Collect operator parameters                                        */
/*--------------------------------------------------------------------*/

/* Parameters for a color vector field */

static int ask_field_op( FILE *fp, int prompt, int *source_type, char *descrp)
{
  char *savebuf;
  char myname[] = "ask_field_op";

  if (prompt==1){
    printf("enter ");
    printf("'covariant_gaussian', ");
    printf("'complex_field', ");
    printf("'complex_field_fm', ");
    printf("'evenandodd_wall', ");
    printf("'gaussian', ");
    printf("'identity', ");
    printf("'wavefunction', ");
    printf("\n     ");
    printf("'deriv1', ");
    printf("'deriv2_D', ");
    printf("'deriv2_B', ");
    printf("'deriv3_A', ");
    printf("\n     ");
    printf("'fat_covariant_gaussian', ");
    printf("'funnywall1', ");
    printf("'funnywall2', ");
    printf("'hop', ");
    printf("'rotate_3D', ");
    printf("'spin_taste',");
    printf("\n     ");
    printf(", for source type\n");
  }

  savebuf = get_next_tag(fp, "quark source op command", myname);
  if (savebuf == NULL)return 1;

  /* Convolutions */
  if(strcmp("complex_field",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FILE;
    strcpy(descrp,"complex_field");
  }
  else if(strcmp("complex_field_fm",savebuf) == 0 ) {
    *source_type = COMPLEX_FIELD_FM_FILE;
    strcpy(descrp,"complex_field_fm");
  }
  else if(strcmp("evenandodd_wall",savebuf) == 0 ) {
    *source_type = EVENANDODD_WALL;
    strcpy(descrp,"even_and_odd_wall");
  }
  else if(strcmp("gaussian",savebuf) == 0 ) {
    *source_type = GAUSSIAN;
    strcpy(descrp,"gaussian");
  }
  else if(strcmp("identity",savebuf) == 0 ){
    *source_type = IDENTITY;
    strcpy(descrp,"identity");
  }
  else if(strcmp("wavefunction",savebuf) == 0 ){
    *source_type = WAVEFUNCTION_FILE;
    strcpy(descrp,"wavefunction");
  }

  /* Local operators */
  else if(strcmp("covariant_gaussian",savebuf) == 0 ) {
    *source_type = COVARIANT_GAUSSIAN;
    strcpy(descrp,"covariant_gaussian");
  }
  else if(strcmp("deriv1",savebuf) == 0 ) {
    *source_type = DERIV1;
    strcpy(descrp,"deriv1");
  }
  else if(strcmp("deriv2_D",savebuf) == 0 ) {
    *source_type = DERIV2_D;
    strcpy(descrp,"deriv2_D");
  }
  else if(strcmp("deriv2_B",savebuf) == 0 ) {
    *source_type = DERIV2_B;
    strcpy(descrp,"deriv2_B");
  }
  else if(strcmp("deriv3_A",savebuf) == 0 ) {
    *source_type = DERIV3_A;
    strcpy(descrp,"deriv3_A");
  }
  else if(strcmp("fat_covariant_gaussian",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN;
    strcpy(descrp,"fat_covariant_gaussian");
  }
  /* Funnywall1 operator (couples to pion5, pioni5, pioni, pions, rhoi, rhos) */
  else if(strcmp("funnywall1",savebuf) == 0 ){
    *source_type = FUNNYWALL1;
    strcpy(descrp,"FUNNYWALL1");
  }
  /* Funnywall2 operator (couples to pion05, pionij, pioni0, pion0, rhoi0, rho0) */
  else if(strcmp("funnywall2",savebuf) == 0 ){
    *source_type = FUNNYWALL2;
    strcpy(descrp,"FUNNYWALL2");
  }
  else if(strcmp("hop",savebuf) == 0 ){
    *source_type = HOPPING;
    strcpy(descrp,"hop");
  }
  else if(strcmp("rotate_3D",savebuf) == 0 ){
    *source_type = ROTATE_3D;
    strcpy(descrp,"rotate_3D");
  }
  else if(strcmp("spin_taste",savebuf) == 0 ){
    *source_type = SPIN_TASTE;
    strcpy(descrp,"SPIN_TASTE");
  }
  else{
    printf("%s: ERROR IN INPUT: field operation command %s is invalid\n",myname,
	   savebuf); 
    return 1;
  }

  printf("%s\n",savebuf);

#ifdef NO_GAUGE_FIELD
  if(requires_gauge_field(*source_type)){
    printf("%s: ERROR IN INPUT: field operation not supported for this application\n"
	   ,myname);
    return 1;
  }
#endif
  
  return 0;

} /* ask_field_op */

#define IF_OK if(status==0)

static char *encode_dir(int dir){
  if(dir == XUP)return "x";
  else if(dir == YUP)return "y";
  else if(dir == ZUP)return "z";
  else return "?";
}


/* For parsing the derivative direction */
static int decode_dir(int *dir, char c_dir[]){
  int status = 0;
  if(strcmp(c_dir,"x")==0)
    *dir = XUP;
  else if(strcmp(c_dir,"y")==0)
    *dir = YUP;
  else if(strcmp(c_dir,"z")==0)
    *dir = ZUP;
  else {
    node0_printf("Expected 'x' or 'y' or 'z' but got %s\n",c_dir);
    status++;
  }
  return status;
} /* decode_dir */

static int get_field_op(int *status_p, FILE *fp, 
			int prompt, quark_source_sink_op *qss_op){

  Real a = 0;
  char c_dir0[] = " ";
  char c_dir1[2] = " ";
  char *c_dir[2] = {c_dir0, c_dir1};
  int  source_iters = 0;
  Real source_r0 = 0;
  int  op_type = qss_op->type;
  char source_file[MAXFILENAME] = "";
  int  status = *status_p;

  /*********************************************************************************/
  /* Operators common to Dirac and KS                                              */
  /*********************************************************************************/

  /* Convolutions */
  if ( op_type == COVARIANT_GAUSSIAN ){
    IF_OK status += get_f(fp, prompt, "r0", &source_r0);
    IF_OK status += get_i(fp, prompt, "source_iters", &source_iters);
  }
  else if ( op_type == COMPLEX_FIELD_FILE ||
	    op_type == COMPLEX_FIELD_FM_FILE ){
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
  }
  else if ( op_type == EVENANDODD_WALL ){
    /* No additional parameters needed */
  }
  else if ( op_type == GAUSSIAN ){
    /* width: psi=exp(-(r/r0)^2) */
    IF_OK status += get_f(fp, prompt,"r0", &source_r0 );
  }
  else if ( op_type == IDENTITY ){
    /* No additional parameters needed */
  }
  else if ( op_type == WAVEFUNCTION_FILE ){
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
    IF_OK status += get_f(fp, prompt, "a", &a);
  }
  
  /* Local operators */
  else if( op_type == DERIV1){
    /* Parameters for derivative */
    IF_OK status += get_vs(fp, prompt, "dir", c_dir, 1);
    IF_OK status += decode_dir(&qss_op->dir1, c_dir[0]);
    IF_OK status += get_i(fp, prompt, "disp", &qss_op->disp);
    IF_OK status += get_vf(fp, prompt, "weights", qss_op->weights, qss_op->disp);
  }
  else if( op_type == DERIV2_D ||
	   op_type == DERIV2_B ){
    /* Parameters for derivatives */
    IF_OK status += get_vs(fp, prompt, "dir", c_dir, 2);
    IF_OK status += decode_dir(&qss_op->dir1, c_dir[0]);
    IF_OK status += decode_dir(&qss_op->dir2, c_dir[1]);
    IF_OK status += get_i(fp, prompt, "disp", &qss_op->disp);
    IF_OK status += get_vf(fp, prompt, "weights", qss_op->weights, qss_op->disp);
  }
  else if( op_type == DERIV3_A){
    /* Parameters for derivative */
    IF_OK status += get_i(fp, prompt, "disp", &qss_op->disp);
    IF_OK status += get_vf(fp, prompt, "weights", qss_op->weights, qss_op->disp);
  }
  else if( op_type == FAT_COVARIANT_GAUSSIAN){
    /* Parameters for covariant Gaussian */
    IF_OK status += get_f(fp, prompt, "r0", &source_r0);
    IF_OK status += get_i(fp, prompt, "source_iters", &source_iters);
  }
  else if( op_type == HOPPING){
    /* Parameters for hopping matrix */
    IF_OK status += get_i(fp, prompt, "derivs",   &qss_op->dhop);
    IF_OK status += get_vs(fp, prompt, "dir", c_dir, 1);
    IF_OK status += decode_dir(&qss_op->dir1, c_dir[0]);
#if FERM_ACTION == HISQ
    IF_OK status += get_f(fp, prompt, "eps_naik", &qss_op->eps_naik);
#endif
  }

  else {
    return 0;
  }

  qss_op->a             = a;
  qss_op->iters         = source_iters;
  qss_op->r0            = source_r0;
  strcpy(qss_op->source_file,source_file);

  *status_p = status;
  return 1;
}/* get_field_op */
  
#ifdef HAVE_DIRAC

/*--------------------------------------------------------------------*/
/* Get parameters for operations on Dirac vector fields                     */
/*--------------------------------------------------------------------*/

int get_wv_field_op(FILE *fp, int prompt, quark_source_sink_op *qss_op){
  
  char myname[] = "get_wv_field_op";
  Real d1 = 0;
  int  op_type;
  char op_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_field_op(fp, prompt, &op_type, qss_op->descrp);
  IF_OK qss_op->type  = op_type;

  /* Get source parameters */
  IF_OK {

    /* Operators common to Dirac and KS */
    if( get_field_op(&status, fp, prompt, qss_op) );

    /* Other operators exclusive to Dirac */
    else if ( op_type == ROTATE_3D ){
      IF_OK status += get_f(fp, prompt, "d1", &d1);
    }

    else {
      printf("(%s)Dirac operator type %d not supported in this application\n",
	     myname, op_type);
      status++;
    }
  }

  qss_op->d1            = d1;

  IF_OK status += get_s(fp, prompt, "op_label", op_label);
  /* Source label '(null)' suppresses the label */
  if(strcmp(op_label,"(null)")==0)op_label[0]='\0';
  strncpy(qss_op->label, op_label,MAXSRCLABEL);
  
  return status;
} /* get_wv_field_op */

#endif

#ifdef HAVE_KS

/*--------------------------------------------------------------------*/
/* Get parameters for operations on color vector fields               */
/*--------------------------------------------------------------------*/

int get_v_field_op(FILE *fp, int prompt, quark_source_sink_op *qss_op){
  
  int  op_type;
  char op_label[MAXSRCLABEL];
  int  status = 0;
  
  /* Get quark source type */
  IF_OK status += ask_field_op(fp, prompt, &op_type, qss_op->descrp);
  IF_OK qss_op->type  = op_type;

  /* Get source parameters */
  IF_OK {

    /* Operators common to Dirac and KS */
    if( get_field_op(&status, fp, prompt, qss_op) );

    /* Other operators exclusive to KS */
    else if ( op_type == FUNNYWALL1 ||
	      op_type == FUNNYWALL2 );   /* No additional parameters needed for these */
    
    else if( op_type == SPIN_TASTE){
      char spin_taste_label[8];
      /* Parameters for spin-taste */
      IF_OK status += get_s(fp, prompt, "gamma", spin_taste_label);
      IF_OK {
	qss_op->spin_taste = spin_taste_index(spin_taste_label);
	if(qss_op->spin_taste < 0){
	  printf("\n: Unrecognized spin-taste label %s.\n", spin_taste_label);
	  status++;
	}
      }
    }

    else {
      printf("KS operator type %d not supported in this application\n", op_type);
      status++;
    }
  }

  IF_OK status += get_s(stdin, prompt, "op_label", op_label);
  /* Source label '(null)' suppresses the label */
  if(strcmp(op_label,"(null)")==0)op_label[0]='\0';
  strncpy(qss_op->label,op_label,MAXSRCLABEL);
  
  return status;
} /* get_v_field_op */

#endif


/*--------------------------------------------------------------------*/
/* Print the parameters of the operators used.
   This is intended to be part of the metadata
   in the FNAL-style correlator file*/
/*--------------------------------------------------------------------*/

#define NTAG 31  /* Print line space available for a tag */
/* Create a fixed-width tag for tidy output */
static char *make_tag(char prefix[], char tag[]){
  static char full_tag[NTAG];
  full_tag[0] = '\0';
  strncat(full_tag, prefix, NTAG-1);
  /* Insert a dash if the prefix is nonblank */
  if(strspn(prefix," ") != strlen(prefix))
    strncat(full_tag, "_", NTAG-1-strlen(full_tag));
  strncat(full_tag, tag, NTAG-1-strlen(full_tag));
  strncat(full_tag, ":", NTAG-1-strlen(full_tag));
  memset(full_tag+strlen(full_tag), ' ', NTAG-1-strlen(full_tag));
  full_tag[NTAG-1] = '\0';
  return full_tag;
}

/*--------------------------------------------------------------------*/
static int print_single_op_info(FILE *fp, char prefix[], 
				quark_source_sink_op *qss_op){
  int op_type = qss_op->type;
  int status = 1;

  if ( op_type == COVARIANT_GAUSSIAN ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g,\n", make_tag(prefix, "r0"), qss_op->r0);
    fprintf(fp,"%s%d\n", make_tag(prefix, "iters"), qss_op->iters);
  }
  else if ( op_type == COMPLEX_FIELD_FILE ||
	    op_type == COMPLEX_FIELD_FM_FILE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s\n", make_tag(prefix, "file"), qss_op->source_file);
  }
  else if( op_type == DERIV1){
    int k;
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "dir"), encode_dir(qss_op->dir1));
    fprintf(fp,"%s%d,\n", make_tag(prefix, "disp"), qss_op->disp);
    fprintf(fp,"%s%g", make_tag(prefix, "weights"), qss_op->weights[0]);
    for(k = 1; k < qss_op->disp; k++)fprintf(fp," %g", qss_op->weights[k]);
    fprintf(fp,"\n");
  }
  else if( op_type == DERIV2_D ||
	   op_type == DERIV2_B ){
    int k;
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "dir1"), encode_dir(qss_op->dir1));
    fprintf(fp,"%s%s,\n", make_tag(prefix, "dir2"), encode_dir(qss_op->dir2));
    fprintf(fp,"%s%d,\n", make_tag(prefix, "disp"), qss_op->disp);
    fprintf(fp,"%s%g", make_tag(prefix, "weights"), qss_op->weights[0]);
    for(k = 1; k < qss_op->disp; k++)fprintf(fp," %g", qss_op->weights[k]);
    fprintf(fp,"\n");
  }
  else if( op_type == DERIV3_A){
    int k;
    fprintf(fp,",\n");
    fprintf(fp,"%s%d,\n", make_tag(prefix, "disp"), qss_op->disp);
    fprintf(fp,"%s%g", make_tag(prefix, "weights"), qss_op->weights[0]);
    for(k = 1; k < qss_op->disp; k++)fprintf(fp," %g", qss_op->weights[k]);
    fprintf(fp,"\n");
  }
  else if( op_type == FAT_COVARIANT_GAUSSIAN){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g,\n", make_tag(prefix, "r0"), qss_op->r0);
    fprintf(fp,"%s%d\n", make_tag(prefix, "iters"), qss_op->iters);
  }
  else if ( op_type == GAUSSIAN ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g\n", make_tag(prefix, "r0"), qss_op->r0);
  }
  else if ( op_type == HOPPING ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d\n", make_tag(prefix, "derivs"), qss_op->dhop);
    fprintf(fp,"%s%s\n", make_tag(prefix, "dir"), encode_dir(qss_op->dir1));
#if FERM_ACTION == HISQ
    fprintf(fp,"%s%g\n", make_tag(prefix, "eps_naik"), qss_op->eps_naik);
#endif
  }
  else if ( op_type == ROTATE_3D ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g\n", make_tag(prefix, "d1"), qss_op->d1);
  }
#ifdef HAVE_KS
  else if ( op_type == SPIN_TASTE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s\n", make_tag(prefix, "gamma"), 
	    spin_taste_label(qss_op->spin_taste));
  }
#endif
  else if ( op_type == WAVEFUNCTION_FILE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "file"), qss_op->source_file);
    fprintf(fp,"%s%g\n", make_tag(prefix, "a"), qss_op->a);
  }
  else {
    fprintf(fp,"\n");
    status = 0;
  }

  return status;
}

/*--------------------------------------------------------------------*/
void print_field_op_info(FILE *fp, char prefix[], 
			 quark_source_sink_op *qss_op){

  quark_source_sink_op* qs_op = qss_op;

  /* Suppress printing of trivial operations */
  if(qs_op == NULL)return;
  if(( qs_op->type == IDENTITY ||
       qs_op->type == UNKNOWN    ) && qs_op->op == NULL)return;

  fprintf(fp,"%s[\n", make_tag(prefix, "ops"));
  
  while(qs_op != NULL){
    fprintf(fp,"{ operation:                  %s",qs_op->descrp);
    print_single_op_info(fp, "  ", qs_op);
    qs_op = qs_op->op;
    if(qs_op != NULL)fprintf(fp, "},\n");
    else fprintf(fp, "}\n");
  }

  fprintf(fp, "]\n");
}

/*--------------------------------------------------------------------*/
/* Same as above, but take data from a flat array of ops */
void print_field_op_info_list(FILE *fp, char prefix[], 
			      quark_source_sink_op *qss_op[], int n){

  int i;

  /* Suppress printing of trivial operations */
  if(qss_op[0] == NULL)return;
  if(n == 0)return;
  if(n == 1 && ( qss_op[0]->type == IDENTITY ||
		 qss_op[0]->type == UNKNOWN) )
     return;

  fprintf(fp,"%s[\n", make_tag(prefix, "ops"));
  
  for(i = 0; i < n; i++){
    fprintf(fp,"{ operation:                  %s",qss_op[i]->descrp);
    print_single_op_info(fp, "  ", qss_op[i]);
    if(i < n-1)fprintf(fp, "},\n");
    else fprintf(fp, "}\n");
  }

  fprintf(fp, "]\n");
}


/* quark_source_sink_op.c */
