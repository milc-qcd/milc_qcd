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

   Operations exclusive to Dirac vectors:

   dirac_inverse              Multiply by a Dirac propagator
   
   ext_src_dirac              Project to a time slice, apply momentum and gamma
   gamma                      Multiply by a gamma factor
   hop                        Multiply by the Dirac hopping matrix (or its derivative)
   ks_gamma                   Multiply by the staggered gamma factor
   ks_gamma_inv               Multiply by the inverse of the staggered gamma factor

   Operations exclusive to color vectors:

   aslash_ks                  Insert A_mu gamma_mu (currently only staggered)
   ext_src_ks                 Project to a time slice, apply momentum and gamma
   funnywall1                 pion5 + pioni5 + pioni + pions + rhoi + rhos
   funnywall2                 pion05 + pionij + pioni0 + pion0 + rhoi0 + rho0
   ks_inverse                 Multiply by a staggered propagator
   spin_taste
   spin_taste_extend

   Noncovariant modifications of the source:

   modulation                 Multiply by a complex field
   momentum                   Multiply by a Fourier phase
   project_t_slice            Project to a time slice
   
   Covariant modifications of the source:

   covariant_gaussian              Smear with exp(const * Laplacian)
   fat_covariant_gaussian          Same, but base covariant Laplacian on 
   fat_covariant_laplacian         Covariant Laplacian itself
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
   insert_qss_eps_naik_index

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
#include "../include/gammatypes.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#endif
#ifdef HAVE_DIRAC
#include "../include/generic_clover.h"
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
  qss_op->stride           = 1;
  qss_op->r_offset[0]      = 0;
  qss_op->r_offset[1]      = 0;
  qss_op->r_offset[2]      = 0;
  qss_op->r_offset[3]      = 0;
  qss_op->spin_taste       = -1;
  qss_op->gamma            = 0;
  qss_op->mom[0]           = 0.;
  qss_op->mom[1]           = 0.;
  qss_op->mom[2]           = 0.;
  qss_op->source_file[0]   = '\0';
  qss_op->dcp.Kappa        = 0;
  qss_op->dcp.Clov_c       = 1.;
  qss_op->dcp.U0           = 1.;
  qss_op->co[0]            = 0;
  qss_op->co[1]            = 0;
  qss_op->co[2]            = 0;
  qss_op->co[3]            = 0;
  qss_op->bp[0]            = 0.;
  qss_op->bp[1]            = 0.;
  qss_op->bp[2]            = 0.;
  qss_op->bp[3]            = 0.;
  qss_op->t0               = 0;
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

/* Accessor for Naik epsilon parameter in embedded KS inverse and hopping operator */
/* Returns 0 if a Naik epsilon is not used for this operator */
int get_qss_eps_naik(Real *eps_naik, quark_source_sink_op *qss_op){
  if(qss_op->type == KS_INVERSE || qss_op->type == HOPPING){
    *eps_naik = qss_op->eps_naik;
    return 1;
  }
  return 0;
}

void insert_qss_eps_naik_index(int index, quark_source_sink_op *qss_op){
  qss_op->ksp.naik_term_epsilon_index = index;
}

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

  chi_cs = create_c_field();

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
  destroy_c_field(chi_cs);
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
      vp[i] = v[i];
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
				    su3_matrix *t_links, int subset,
				    Real r0, int iters, int t0){
  wilson_vector *wv;
  int spin, color;
  int subset, stride;

  if(subset == FULL)stride = 1;
  else stride = 2;

  wv = create_wv_field();
  for(color = 0; color < wp->nc; color++)
    for(spin=0;spin<4;spin++){
      copy_wv_from_wp(wv, wp, color, spin);
      gauss_smear_wv_field(wv, stride, t_links, r0, iters, t0);
      copy_wp_from_wv(wp, wv, color, spin);
    }
  
  destroy_wv_field(wv);
} /* gauss_smear_wprop_field */

/* On ks_prop_field type */

static void gauss_smear_ksprop_field(ksprop_field *ksp, 
				     su3_matrix *t_links,
				     Real r0, int iters, int t0){
  su3_vector *v;
  int color, stride;

  v = create_v_field();
  for(color = 0; color < ksp->nc; color++){
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

static void hop_wvec(wilson_vector *src, int dhop, int mu, int fb)
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
  hop_w_field(src, mp, PLUS, sign, EVENANDODD, mu, fb);
  
  FORALLFIELDSITES(i){
    copy_wvec(mp + i, src + i);
  }

  destroy_wv_field(mp); 

} /* hop_wvec */

#endif

/*--------------------------------------------------------------------*/
/* Multiply by KS Gamma                                               */
/*--------------------------------------------------------------------*/

/* Apply the inverse of the Kawamoto Smit operator to the Dirac field */

static void ks_gamma(wilson_vector *src, int r0[])
{
  char myname[] = "ks_gamma";
  
  if(src == NULL){
    node0_printf("%s: Error: called with NULL arg\n", myname);
    terminate(1);
  }

  /* Apply KS Gamma to the src field */
  mult_by_ks_gamma_wv(src, r0);

} /* ks_gamma */

/*--------------------------------------------------------------------*/
/* Multiply by KS Gamma inverse                                               */
/*--------------------------------------------------------------------*/

/* Apply the inverse of the Kawamoto Smit operator to the Dirac field */

static void ks_gamma_inv(wilson_vector *src, int r0[])
{
  char myname[] = "ks_gamma_inv";
  
  if(src == NULL){
    node0_printf("%s: Error: called with NULL arg\n", myname);
    terminate(1);
  }

  /* Apply KS Gamma inverse to the source field */
  mult_by_ks_gamma_inv_wv(src, r0);

} /* ks_gamma_inv */

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
    op_type == FAT_COVARIANT_LAPLACIAN ||
    op_type == FUNNYWALL1 ||
    op_type == FUNNYWALL2 ||
    op_type == HOPPING    ||
    op_type == KS_INVERSE ||
    op_type == ROTATE_3D  ||
    op_type == SPIN_TASTE ||
    op_type == SPIN_TASTE_EXTEND;
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
  int stride        = qss_op->stride;
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
    fnal_wavefunction(chi_cs, 0, 0, 0, t0, stride, a, sink_file);
  }
  else {
    destroy_c_field(chi_cs);
    return NULL;
  }
  return chi_cs;
}

#ifdef HAVE_KS

static int apply_convolution_v(su3_vector *src, quark_source_sink_op *qss_op,
			     int t0){

  complex *chi_cs   = get_convolving_field(qss_op, t0);

  if(chi_cs == NULL){
    destroy_c_field(chi_cs);
    return 0;
  }

  /* Do the convolution */
  smear_v_field(src, chi_cs);
  destroy_c_field(chi_cs);

  return 1;
}

static void apply_modulation_v(su3_vector *src, 
			       quark_source_sink_op *qss_op,
			       int t0){
  
#ifndef HAVE_QIO
  char myname[] = "apply_modulation_wv";
#endif
  char *modulation_file  = qss_op->source_file;
  complex *chi_cs  = create_c_field();
  int i,cf;
  site *s;
  complex z;

  /* Read complex modulating field */

#ifdef HAVE_QIO
  restore_complex_scidac_to_field(modulation_file, QIO_SERIAL, chi_cs, 1);
#else
  node0_printf("%s QIO compilation required for this source\n", myname);
  terminate(1);
#endif
  
  /* src  <- src * chi_cs */
  FORALLSITES(i,s){
    if(s->t == t0 || t0 == ALL_T_SLICES){
      for(cf=0;cf<3;cf++){
	  z = src[i].c[cf];
	  CMUL(z, chi_cs[i], src[i].c[cf]);
	}
    }
  }
}

static void apply_momentum_v(su3_vector *src, 
			     quark_source_sink_op *qss_op,
			     int t0){
  int i; site *s;
  int c;
  Real pi = 3.141592653589793;
  Real th;
  Real px, py, pz;
  int x0, y0, z0;
  complex z,y;
  
  if(qss_op->mom[0] == 0 && qss_op->mom[1] == 0 && qss_op->mom[2] == 0)return;

  px = 2*pi*qss_op->mom[0]/nx;
  py = 2*pi*qss_op->mom[1]/ny;
  pz = 2*pi*qss_op->mom[2]/nz;

  x0 = qss_op->r_offset[0];
  y0 = qss_op->r_offset[1];
  z0 = qss_op->r_offset[2];

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++){
      z = src[i].c[c];
      CMUL(z,y,src[i].c[c]);
    }
  }
}
  
static void apply_tslice_projection_v(su3_vector *src, 
				      quark_source_sink_op *qss_op){
  int i;
  site *s;
  FORALLSITES(i,s){
    if(s->t != qss_op->t0)clearvec(src+i);
  }
}

static void apply_spin_taste(su3_vector *src, quark_source_sink_op *qss_op){
  su3_vector *dst = create_v_field();
  
  /* At present we do not support "fn" type spin-taste operators here */
  spin_taste_op(qss_op->spin_taste, qss_op->r_offset, dst, src);

  copy_v_field(src, dst);
  destroy_v_field(dst);
  
}

static void apply_spin_taste_extend(su3_vector *src, quark_source_sink_op *qss_op){
  su3_vector *dst = create_v_field();
  
  /* At present we do not support "fn" type spin-taste operators here */
  spin_taste_op(qss_op->spin_taste, qss_op->r_offset, dst, src);

  /* The spin_taste_op phases were designed for tying together two
       forward propagators by first converting one of them to an
       antiquark propagator.  So they include the antiquark phase
       (-)^{x+y+z+t}.  For an extended source we don't want the
       antiquark phase.  The next step removes it.  That way the user
       input snk_spin_taste label describes the actual meson at the
       extended source. */
  spin_taste_op(spin_taste_index("pion05"), qss_op->r_offset, src, dst);

  
  //  copy_v_field(src, dst);
  destroy_v_field(dst);
  
}

#define NMU 4

/* Create spin-taste indices for current */
static int *
get_spin_taste(void){
  
  /* Current spin-taste list */
  char *spin_taste_label[NMU] = {"GX-GX", "GY-GY", "GZ-GZ", "GT-GT"};
  static int spin_taste[NMU];
  int mu;
  
  /* Decode spin-taste label */
  for(mu = 0; mu < NMU; mu++){
    char dummy[6];
    strncpy(dummy, spin_taste_label[mu], 6);
    spin_taste[mu] = spin_taste_index(dummy);
  }
  
  return spin_taste;
}

static void apply_aslash_v(su3_vector *src, 
			   quark_source_sink_op *qss_op,
			   int t0){
  
#ifndef HAVE_QIO
  char myname[] = "apply_aslash";
#endif
  char *modulation_file  = qss_op->source_file;
  complex *chi_cs  = create_c_array_field(4);
  su3_vector *dst = create_v_field();
  su3_vector *tmp = create_v_field();
  int i,cf;
  site *s;
  int *spin_taste = get_spin_taste();
  int mu;

  /* Read complex modulating field */

#ifdef HAVE_QIO
  restore_complex_scidac_to_field(modulation_file, QIO_SERIAL, chi_cs, NMU);
#else
  node0_printf("%s QIO compilation required for this source\n", myname);
  terminate(1);
#endif

  for(mu = 0; mu < NMU; mu++){
    qss_op->spin_taste = spin_taste[mu];
    copy_v_field(tmp, src);
    /* tmp <- src * gamma_mu */
    apply_spin_taste(tmp, qss_op);

    /* dst += tmp * A_mu */
    FORALLSITES(i,s){
      if(s->t == t0 || t0 == ALL_T_SLICES){
	for(cf=0;cf<3;cf++){
	  complex z = src[i].c[cf];
	  CMUL(z, chi_cs[NMU*i + mu], tmp[i].c[cf]);
	  CSUM(dst[i].c[cf],z);
	}
      }
    }
  }
  destroy_v_field(tmp);
  spin_taste_op(spin_taste_index("pion05"), qss_op->r_offset, src, dst);
  destroy_v_field(dst);
  destroy_c_array_field(chi_cs, NMU);
}

static void apply_ext_src_v(su3_vector *src, 
			    quark_source_sink_op *qss_op){
  apply_tslice_projection_v(src, qss_op);
  apply_spin_taste_extend(src, qss_op);
  apply_momentum_v(src, qss_op, qss_op->t0);
}

#ifndef NO_GAUGE_FIELD

static int apply_ks_inverse(su3_vector *v, quark_source_sink_op *qss_op,
			    int t0){
  char myname[] = "apply_ks_inverse";
  int i; site *s;
  int tot_iters = 0;
  su3_vector *src = create_v_field();
  quark_invert_control *my_qic = &qss_op->qic;
  ks_param *my_ksp             = &qss_op->ksp;
  int inaik                    = my_ksp->naik_term_epsilon_index;
  Real *bdry_phase             = qss_op->bp;
  int *r0                      = qss_op->co;
  Real mybdry_phase[4];
  imp_ferm_links_t *fn;

  if(src == NULL){
    printf("%s: No room\n", myname);
    terminate(1);
  }

  /* Copy field to source */
  copy_v_field(src, v);
  clear_v_field(v);

  /* Clear all values except on specified time slice */
  if(t0 != ALL_T_SLICES)
    FORALLSITES(i,s){
      if(s->t != t0)
	clearvec(src+i);
    }

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];
  
  /* Get fn links appropraite to this Naik term epsilon */
  
  restore_fermion_links_from_site(fn_links, my_qic->prec);
  fn = get_fm_links(fn_links)[inaik];

  /* Apply twist to the boundary links of fn and reset origin of KS
     phases if requested */
  set_boundary_twist_fn(fn, mybdry_phase, r0);
  boundary_twist_fn(fn, ON);

  /* If we are twisting, apply the momentum twist to the source */
  /* See ks_spectrum/make_prop.c for an explanation */
  mybdry_phase[3] = 0; 
  rephase_v_field(src, mybdry_phase, r0, 1);
  mybdry_phase[3] = bdry_phase[3]; 
  
  /* Solve for propagator */
  tot_iters += mat_invert_uml_field(src, v, my_qic, my_ksp->mass, fn);
  
  /* Undo momentum twist on sink */
  mybdry_phase[3] = 0; 
  rephase_v_field(v, mybdry_phase, r0, -1);
  mybdry_phase[3] = bdry_phase[3]; 
  
  /* Undo twist on fn links */
  boundary_twist_fn(fn, OFF);

  destroy_v_field(src);
  return tot_iters;
}

#endif /* ifndef NO_GAUGE_FIELD */
#endif /* HAVE_KS */

#ifdef HAVE_DIRAC

static int apply_convolution_wv(wilson_vector *src, 
			      quark_source_sink_op *qss_op,
			      int t0){

  complex *chi_cs   = get_convolving_field(qss_op, t0);

  if(chi_cs == NULL){
    return 0;
  }

  /* Do the convolution */

  smear_wv_field(src, chi_cs);
  free(chi_cs);

  return 1;

}

static void apply_modulation_wv(wilson_vector *src, 
				quark_source_sink_op *qss_op,
				int t0){
  
#ifndef HAVE_QIO
  char myname[] = "apply_modulation_wv";
#endif
  char *modulation_file  = qss_op->source_file;
  complex *chi_cs  = create_c_field();
  int i,sf,cf;
  site *s;
  complex z;

  /* Read complex modulating field */

#ifdef HAVE_QIO
  restore_complex_scidac_to_field(modulation_file, QIO_SERIAL, chi_cs, 1);
#else
  node0_printf("%s QIO compilation required for this source\n", myname);
  terminate(1);
#endif
  
  /* src  <- src * chi_cs */
  FORALLSITES(i,s){
    if(s->t == t0 || t0 == ALL_T_SLICES){
      for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
	  z = src[i].d[sf].c[cf];
	  CMUL(z, chi_cs[i], src[i].d[sf].c[cf]);
	}
    }
  }
}

static void apply_momentum_wv(wilson_vector *src, 
			      quark_source_sink_op *qss_op,
			      int t0){
  int i; site *s;
  int c,d;
  Real pi = 3.14159265358979324;
  Real th;
  Real px, py, pz;
  int x0, y0, z0;
  complex z,y;
  
  if(qss_op->mom[0] == 0 && qss_op->mom[1] == 0 && qss_op->mom[2] == 0)return;

  px = 2*pi*qss_op->mom[0]/nx;
  py = 2*pi*qss_op->mom[1]/ny;
  pz = 2*pi*qss_op->mom[2]/nz;

  x0 = qss_op->r_offset[0];
  y0 = qss_op->r_offset[1];
  z0 = qss_op->r_offset[2];

  FORALLSITES(i,s)if(t0 == ALL_T_SLICES || s->t == t0){
    th = px*(s->x-x0) + py*(s->y-y0) + pz*(s->z-z0);
    y.real = cos(th); y.imag = sin(th);
    for(c=0;c<3;c++)
      for(d=0;d<4;d++){
	z = src[i].d[d].c[c];
	CMUL(z,y,src[i].d[d].c[c]);
      }
  }
}
  
static void apply_tslice_projection_wv(wilson_vector *src, 
				       quark_source_sink_op *qss_op){
  int i;
  site *s;
  FORALLSITES(i,s){
    if(s->t != qss_op->t0)clear_wvec(src+i);
  }
}

static void apply_gamma(wilson_vector *src, 
			quark_source_sink_op *qss_op){
  int i;
  site *s;
  int gam = qss_op->gamma;
  wilson_vector tmp;

  FORALLSITES(i,s){
    mult_w_by_gamma( src+i, &tmp, gam );
    src[i] = tmp;
  }
}

static void apply_ext_src_wv(wilson_vector *src, 
			     quark_source_sink_op *qss_op){
  apply_tslice_projection_wv(src, qss_op);
  apply_gamma(src, qss_op);
  apply_momentum_wv(src, qss_op, qss_op->t0);
}

#ifndef NO_GAUGE_FIELD 

static int apply_dirac_inverse(wilson_vector *wv,
			       quark_source_sink_op *qss_op,
			       int t0){
  char myname[] = "apply_dirac_inverse";
  int i; site *s;
  int tot_iters = 0;
  wilson_vector *src = create_wv_field();
  quark_invert_control *my_qic = &qss_op->qic;
  dirac_clover_param *my_dcp   = &qss_op->dcp;
  Real *bdry_phase             = qss_op->bp;
  int *r0                      = qss_op->co;
  Real mybdry_phase[4];

  if(src == NULL){
    printf("%s: No room\n", myname);
    terminate(1);
  }

  /* Copy field to source */
  copy_wv_field(src, wv);
  clear_wv_field(wv);

  /* Clear all values except on specified time slice */
  if(t0 != ALL_T_SLICES)
    FORALLSITES(i,s){
      if(s->t != t0)
	clear_wvec(src+i);
    }

  /* Local copy of bdry_phase */
  for(i = 0; i < 4; i++)
    mybdry_phase[i] = bdry_phase[i];

  /* Apply twist to boundary links of the gauge field in site structure */
  boundary_twist_site(mybdry_phase, r0, +1);

  /* Do momentum twist if requested */
  mybdry_phase[3] = 0; 
  rephase_wv_field(src, mybdry_phase, r0, 1);
  mybdry_phase[3] = bdry_phase[3]; 
  
  /* Solve for propagator */
  tot_iters += bicgilu_cl_field(src, wv, my_qic, my_dcp);
  report_status(my_qic);
  
  /* Undo momentum twist */
  mybdry_phase[3] = 0; 
  rephase_wv_field(wv, mybdry_phase, r0, -1);
  mybdry_phase[3] = bdry_phase[3]; 
  
  /* Undo twist to boundary links in gauge field */
  boundary_twist_site(mybdry_phase, r0, -1);
  
  destroy_wv_field(src);
  return tot_iters;
}

#endif /* not NO_GAUGE_FIELD */

#endif /* HAVE_DIRAC */

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
    op_type == FAT_COVARIANT_LAPLACIAN ||
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
    static su3_matrix *t_links;

    t_links = create_G_from_site();
    gauss_smear_v_field(src, t_links, r0, iters, t0);
    destroy_G(t_links);
  }

  else if(op_type == FAT_COVARIANT_GAUSSIAN ){
    gauss_smear_v_field(src, ape_links, r0, iters, t0);
  }

  else if(op_type == FAT_COVARIANT_LAPLACIAN ){
    laplacian_v_field(src, ape_links, t0);
  }

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
  int stride        = qss_op->stride;

  if(op_type == COVARIANT_GAUSSIAN){
    int iters = qss_op->iters;

    static su3_matrix *t_links;

    t_links = create_G_from_site();
    gauss_smear_wv_field(src, t_links, stride, r0, iters, t0);
    destroy_G(t_links);
  }

  else if(op_type == FAT_COVARIANT_GAUSSIAN )
    gauss_smear_wv_field(src, ape_links, stride, r0, iters, t0);

  else if(op_type == FAT_COVARIANT_LAPLACIAN )
    laplacian_wv_field(src, ape_links, stride, t0);

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
    rephase_field_offset( ape_links, ON, NULL, qss_op->r_offset );
    mult_pioni_field( ZUP, qss_op->r_offset, src, tvec1, ape_links );
    rephase_field_offset( ape_links, OFF, NULL, qss_op->r_offset );
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
    rephase_field_offset( ape_links, ON, NULL, qss_op->r_offset );
    mult_pioni0_field( ZUP, qss_op->r_offset, src, tvec1, ape_links );
    rephase_field_offset( ape_links, OFF, NULL, qss_op->r_offset );
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

static void hop_vec(su3_vector *src, ks_param *ksp, int dhop, int mu)
{
  int sign;
  su3_vector *v;
  char myname[] = "hop_vec";
  Real wtfatf, wtfatb, wtlongf, wtlongb;
#if FERM_ACTION == HISQ
  int inaik = ksp->naik_term_epsilon_index;
#else
  int inaik = 0;
#endif
  /* Note: we are not restoring the links here and don't set the precision */
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

  /* Modulation and projection */

  else if(op_type == MODULATION_FILE)
    apply_modulation_v(src, qss_op, t0);

  else if(op_type == MOMENTUM)
    apply_momentum_v(src, qss_op, t0);

  else if(op_type == EXT_SRC_KS)
    apply_ext_src_v(src, qss_op);

  else if(op_type == ASLASH_KS_FILE)
    apply_aslash_v(src, qss_op, t0);

  else if(op_type == PROJECT_T_SLICE)
    apply_tslice_projection_v(src, qss_op);

#ifndef NO_GAUGE_FIELD
  else if(op_type == KS_INVERSE)
    apply_ks_inverse(src, qss_op, t0);

  else if(is_cov_deriv(op_type))
    apply_cov_deriv_v(src, qss_op);

  else if(is_cov_smear(op_type))
    apply_cov_smear_v(src, qss_op, t0);

  else if(is_funnywall(op_type))
    apply_funnywall(src, qss_op);

  else if(op_type == SPIN_TASTE)
    apply_spin_taste(src, qss_op);

  else if(op_type == SPIN_TASTE_EXTEND)
    apply_spin_taste_extend(src, qss_op);

  else if(op_type == HOPPING)
    hop_vec(src, &qss_op->ksp, qss_op->dhop, qss_op->dir1);

#endif

  else {
    node0_printf("%s: Unrecognized operator type %d\n", myname, op_type);
    terminate(1);
  }

  /* Apply subset mask */

  // Not good. Cancels some of these ops
  //  subset_mask_v(src, subset, t0);

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

  /* Modulation and projection */

  else if(op_type == MODULATION_FILE)
    apply_modulation_wv(src, qss_op, t0);

  else if(op_type == MOMENTUM)
    apply_momentum_wv(src, qss_op, t0);

  else if(op_type == EXT_SRC_DIRAC)
    apply_ext_src_wv(src, qss_op);

  else if(op_type == PROJECT_T_SLICE)
    apply_tslice_projection_wv(src, qss_op);

  else if(op_type == GAMMA)
    apply_gamma(src, qss_op);

#ifndef NO_GAUGE_FIELD
  else if(op_type == DIRAC_INVERSE)
    apply_dirac_inverse(src, qss_op, t0);

  else if(is_cov_smear(op_type))
    apply_cov_smear_wv(src, qss_op, t0);

  else if(is_cov_deriv(op_type))
    apply_cov_deriv_wv(src, qss_op);

  else if(op_type == HOPPING)
    hop_wvec(src, qss_op->dhop, qss_op->dir1, qss_op->fb);

  else if(op_type == ROTATE_3D)
    /* The MILC sign convention for gamma matrix in Dslash is
       opposite FNAL, so we rotate with -d1 */
    rotate_3D_wvec(src, -qss_op->d1);

#endif

  else if(op_type == KS_GAMMA)
    ks_gamma(src, qss_op->r_offset);

  else if(op_type == KS_GAMMA_INV)
    ks_gamma_inv(src, qss_op->r_offset);

  else {
    node0_printf("%s: Unrecognized operator type %d\n", myname, op_type);
    terminate(1);
  }

  /* Apply subset mask */

  /* Do we ever need this? */
  //  subset_mask_wv(src, subset, t0);

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
    printf("'fat_covariant_laplacian', ");
    printf("'funnywall1', ");
    printf("'funnywall2', ");
    printf("'dirac_inverse', ");
    printf("'hop', ");
    printf("'ks_gamma', ");
    printf("'ks_gamma_inv', ");
    printf("'ks_inverse', ");
    printf("'rotate_3D', ");
    printf("'gamma', ");
    printf("'spin_taste', ");
    printf("'spin_taste_extend', ");
    printf("\n     ");
    printf("'momentum', ");
    printf("'modulation', ");
    printf("'project_t_slice', ");
    printf("'ext_src_ks', ");
    printf("'ext_src_dirac', ");
    printf("\n     ");
    printf(", for field op\n");
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
  else if(strcmp("dirac_inverse",savebuf) == 0 ) {
    *source_type = DIRAC_INVERSE;
    strcpy(descrp,"dirac_inverse");
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
  else if(strcmp("ks_inverse",savebuf) == 0 ) {
    *source_type = KS_INVERSE;
    strcpy(descrp,"ks_inverse");
  }
  else if(strcmp("wavefunction",savebuf) == 0 ){
    *source_type = WAVEFUNCTION_FILE;
    strcpy(descrp,"wavefunction");
  }
  /* Modulation and projection */
  else if(strcmp("modulation",savebuf) == 0 ){
    *source_type = MODULATION_FILE;
    strcpy(descrp,"modulation");
  }
  else if(strcmp("momentum",savebuf) == 0 ){
    *source_type = MOMENTUM;
    strcpy(descrp,"momentum");
  }
  else if(strcmp("ext_src_ks",savebuf) == 0 ){
    *source_type = EXT_SRC_KS;
    strcpy(descrp,"ext_src_ks");
  }
  else if(strcmp("aslash_ks",savebuf) == 0 ){
    *source_type = ASLASH_KS_FILE;
    strcpy(descrp,"aslash_ks");
  }
  else if(strcmp("ext_src_dirac",savebuf) == 0 ){
    *source_type = EXT_SRC_DIRAC;
    strcpy(descrp,"ext_src_dirac");
  }
  else if(strcmp("project_t_slice",savebuf) == 0 ){
    *source_type = PROJECT_T_SLICE;
    strcpy(descrp,"project_t_slice");
  }
  else if(strcmp("gamma",savebuf) == 0 ){
    *source_type = GAMMA;
    strcpy(descrp,"gamma");
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
  else if(strcmp("ks_gamma",savebuf) == 0 ) {
    *source_type = KS_GAMMA;
    strcpy(descrp,"ks_gamma");
  }
  else if(strcmp("ks_gamma_inv",savebuf) == 0 ) {
    *source_type = KS_GAMMA_INV;
    strcpy(descrp,"ks_gamma_inv");
  }
  else if(strcmp("fat_covariant_gaussian",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_GAUSSIAN;
    strcpy(descrp,"fat_covariant_gaussian");
  }
  else if(strcmp("fat_covariant_laplacian",savebuf) == 0 ) {
    *source_type = FAT_COVARIANT_LAPLACIAN;
    strcpy(descrp,"fat_covariant_laplacian");
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
    strcpy(descrp,"spin_taste");
  }
  else if(strcmp("spin_taste_extend",savebuf) == 0 ){
    *source_type = SPIN_TASTE_EXTEND;
    strcpy(descrp,"spin_taste_extend");
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
  else if(dir == TUP)return "t";
  else return "?";
}

#if FERM_ACTION == HISQ

static char *encode_sign_dir(int fb, int dir){
  static char sign_dir[3] = "  ";
  char *d = encode_dir(dir);

  if(fb == 0)return d;
  else if(fb == +1)sign_dir[0] = '+';
  else if(fb == -1)sign_dir[0] = '-';
  else return "?";
  sign_dir[1] = d[0];
  return sign_dir;
}

#endif


/* For parsing the derivative direction */
static int decode_dir(int *dir, char c_dir[]){
  int status = 0;
  if(strcmp(c_dir,"x")==0)
    *dir = XUP;
  else if(strcmp(c_dir,"y")==0)
    *dir = YUP;
  else if(strcmp(c_dir,"z")==0)
    *dir = ZUP;
  else if(strcmp(c_dir,"t")==0)
    *dir = TUP;
  else {
    node0_printf("Expected 'x' or 'y' or 'z' or 't' but got %s\n",c_dir);
    status++;
  }
  return status;
} /* decode_dir */

static int get_field_op(int *status_p, FILE *fp, 
			int prompt, quark_source_sink_op *qss_op){

  Real a = 0;
  char c_dir0[3] = " ";
  char c_dir1[3] = " ";
  char *c_dir[2] = {c_dir0, c_dir1};
  int  stride = 1;
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
    IF_OK status += get_i(fp, prompt, "stride", &stride);
    IF_OK status += get_f(fp, prompt, "r0", &source_r0);
    IF_OK status += get_i(fp, prompt, "source_iters", &source_iters);
    if(stride != 1 && this_node==0)printf("NOTE: new convention for r0!\n");
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
    IF_OK status += get_i(fp, prompt, "stride", &stride);
    IF_OK status += get_f(fp, prompt, "a", &a);
  }

  /* Modulation and projections */
  else if ( op_type == MODULATION_FILE ){
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
  }
  else if ( op_type == MOMENTUM){
    IF_OK status += get_vi(fp, prompt, "momentum", qss_op->mom, 3);
  }
  else if ( op_type == PROJECT_T_SLICE ){
    IF_OK status += get_i(fp, prompt, "t0", &qss_op->t0);
  }
#ifdef HAVE_DIRAC
  else if ( op_type == GAMMA ){
    char gam_op_lab[MAXGAMMA];
    IF_OK status += get_s(fp, prompt, "gamma", gam_op_lab);
    IF_OK {
      qss_op->gamma = gamma_index(gam_op_lab);
      if(qss_op->gamma < 0){
	printf("\n%s is not a valid gamma matrix label\n",gam_op_lab);
	status++;
      }
    }
  }
#endif
#ifdef HAVE_KS
  else if ( op_type == EXT_SRC_KS ){
    char gam_op_lab[MAXGAMMA];
    IF_OK status += get_s(fp, prompt, "spin_taste_extend", gam_op_lab);
    IF_OK {
      qss_op->spin_taste = spin_taste_index(gam_op_lab);
      if(qss_op->spin_taste < 0){
	printf("\n%s is not a valid spin-taste label\n",gam_op_lab);
	status++;
      }
    }
    IF_OK status += get_vi(fp, prompt, "momentum", qss_op->mom, 3);
    IF_OK status += get_i(fp, prompt, "t0", &qss_op->t0);
  }
  else if ( op_type == ASLASH_KS_FILE ){
    IF_OK status += get_s(fp, prompt, "load_source", source_file);
    IF_OK status += get_i(fp, prompt, "t0", &qss_op->t0);
  }
#endif
#ifdef HAVE_DIRAC
  else if ( op_type == EXT_SRC_DIRAC ){
    char gam_op_lab[MAXGAMMA];
    IF_OK status += get_s(fp, prompt, "gamma", gam_op_lab);
    IF_OK {
      qss_op->gamma = gamma_index(gam_op_lab);
      if(qss_op->gamma < 0){
	printf("\n%s is not a valid gamma matrix label\n",gam_op_lab);
	status++;
      }
    }
    IF_OK status += get_vi(fp, prompt, "momentum", qss_op->mom, 3);
    IF_OK status += get_i(fp, prompt, "t0", &qss_op->t0);
  }
#endif
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
  else if( op_type == DIRAC_INVERSE){
    char savebuf[128];
    /* Parameters for Dirac inverse */
    IF_OK status += get_s(stdin, prompt,"kappa", qss_op->kappa_label);
    IF_OK qss_op->dcp.Kappa = atof(qss_op->kappa_label);
    IF_OK status += get_f(stdin, prompt,"clov_c", &qss_op->dcp.Clov_c );
    IF_OK status += get_f(stdin, prompt,"u0", &qss_op->dcp.U0 );
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &qss_op->qic.max );
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &qss_op->qic.nrestart );
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			  &qss_op->qic.resid );
    IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			  &qss_op->qic.relresid );
    IF_OK status += get_i(stdin, prompt,"precision", &qss_op->qic.prec );
    IF_OK qss_op->qic.parity = EVENANDODD;
    IF_OK qss_op->qic.min = 0;
    IF_OK qss_op->qic.start_flag = 0;
    IF_OK qss_op->qic.nsrc = 1;
    IF_OK status += get_vi(stdin, prompt, "coordinate_origin", qss_op->co, 4);
    IF_OK status += get_vf(stdin, prompt, "momentum_twist", qss_op->bp, 3);
    IF_OK status += get_s(stdin, prompt,"time_bc", savebuf);
    IF_OK {
      /* NOTE: The Dirac built-in bc is periodic. */
      if(strcmp(savebuf,"antiperiodic") == 0)qss_op->bp[3] = 1;
      else if(strcmp(savebuf,"periodic") == 0)qss_op->bp[3] = 0;
    }
  }
  else if( op_type == FAT_COVARIANT_GAUSSIAN){
    /* Parameters for covariant Gaussian */
    IF_OK status += get_i(fp, prompt, "stride", &stride);
    IF_OK status += get_f(fp, prompt, "r0", &source_r0);
    IF_OK status += get_i(fp, prompt, "source_iters", &source_iters);
    if(stride != 1 && this_node==0)printf("NOTE: new convention for r0!\n");
  }
  else if( op_type == FAT_COVARIANT_LAPLACIAN){
    /* Parameters for covariant Gaussian */
    IF_OK status += get_i(fp, prompt, "stride", &stride);
  }
#ifdef HAVE_KS
  else if( op_type == SPIN_TASTE){
    char spin_taste_label[8];
    /* Parameters for spin-taste */
    IF_OK status += get_s(fp, prompt, "spin_taste", spin_taste_label);
    IF_OK {
      qss_op->spin_taste = spin_taste_index(spin_taste_label);
      if(qss_op->spin_taste < 0){
	printf("\n: Unrecognized spin-taste label %s.\n", spin_taste_label);
	status++;
      }
    }
  }
  else if( op_type == SPIN_TASTE_EXTEND){
    char spin_taste_label[8];
    /* Parameters for spin-taste */
    IF_OK status += get_s(fp, prompt, "spin_taste_extend", spin_taste_label);
    IF_OK {
      qss_op->spin_taste = spin_taste_index(spin_taste_label);
      if(qss_op->spin_taste < 0){
	printf("\n: Unrecognized spin-taste label %s.\n", spin_taste_label);
	status++;
      }
    }
  }
#endif
  else if( op_type == HOPPING){
    /* Parameters for hopping matrix */
    IF_OK status += get_i(fp, prompt, "derivs",   &qss_op->dhop);
    IF_OK status += get_vs(fp, prompt, "dir", c_dir, 1);
    /* Allow a + or - sign to specify a one-sided hop */
    IF_OK {
      if(c_dir[0][0] == '+'){
	qss_op->fb = +1;
	strncpy(c_dir[0], c_dir[0]+1, 2);
      } else if(c_dir[0][0] == '-'){
	qss_op->fb = -1;
	strncpy(c_dir[0], c_dir[0]+1, 2);
      } else {
	qss_op->fb = 0;
      }
    }
    IF_OK status += decode_dir(&qss_op->dir1, c_dir[0]);
#if FERM_ACTION == HISQ
    IF_OK status += get_f(fp, prompt, "eps_naik", &qss_op->eps_naik);
#endif
  }
  else if( op_type == KS_GAMMA || op_type == KS_GAMMA_INV ){
    /* No additional parameters needed */
  }

  else if( op_type == KS_INVERSE){
    char savebuf[128];
    /* Parameters for Dirac inverse */
    IF_OK status += get_s(stdin, prompt,"mass", qss_op->mass_label);
    IF_OK qss_op->ksp.mass = atof(qss_op->mass_label);
    IF_OK status += get_f(stdin, prompt,"naik_term_epsilon", &qss_op->eps_naik );
    IF_OK status += get_f(stdin, prompt,"u0", &qss_op->dcp.U0 );
    IF_OK status += get_i(stdin,prompt,"max_cg_iterations", 
			  &qss_op->qic.max );
    IF_OK status += get_i(stdin,prompt,"max_cg_restarts", 
			  &qss_op->qic.nrestart );
    IF_OK status += get_f(stdin, prompt,"error_for_propagator", 
			  &qss_op->qic.resid );
    IF_OK status += get_f(stdin, prompt,"rel_error_for_propagator", 
			  &qss_op->qic.relresid );
    IF_OK status += get_i(stdin, prompt,"precision", &qss_op->qic.prec );
    IF_OK qss_op->qic.parity = EVENANDODD;
    IF_OK qss_op->qic.min = 0;
    IF_OK qss_op->qic.start_flag = 0;
    IF_OK qss_op->qic.nsrc = 1;
    IF_OK status += get_vi(stdin, prompt, "coordinate_origin", qss_op->co, 4);
    IF_OK status += get_vf(stdin, prompt, "momentum_twist", qss_op->bp, 3);
    IF_OK status += get_s(stdin, prompt,"time_bc", savebuf);
    IF_OK {
      /* NOTE: The KS built-in bc is antiperiodic. */
      /* This is the reverse of the Dirac clover convention */
      if(strcmp(savebuf,"antiperiodic") == 0)qss_op->bp[3] = 0;
      else if(strcmp(savebuf,"periodic") == 0)qss_op->bp[3] = 1;
    }
  }
  else {
    return 0;
  }

  qss_op->a             = a;
  qss_op->iters         = source_iters;
  qss_op->r0            = source_r0;
  qss_op->stride        = stride;
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
    else if ( op_type == KS_GAMMA || op_type == KS_GAMMA_INV );
    /* No additional parameters needed for these */
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

#include "../include/flavor_ops.h"

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
    fprintf(fp,"%s%d,\n", make_tag(prefix, "stride"), qss_op->stride);
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
  else if( op_type == DIRAC_INVERSE){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "kappa"), qss_op->kappa_label);
    fprintf(fp,"%s%g,\n", make_tag(prefix, "clov_c"), qss_op->dcp.Clov_c);
    fprintf(fp,"%s%g,\n", make_tag(prefix, "u0"), qss_op->dcp.U0);
    fprintf(fp,"%s[%d, %d, %d, %d],\n", make_tag(prefix, "coordinate_origin"), 
	    qss_op->co[0], qss_op->co[1], qss_op->co[2], qss_op->co[3]);
    fprintf(fp,"%s[%f, %f, %f],\n", make_tag(prefix, "momentum_twist"), 
	    qss_op->bp[0], qss_op->bp[1], qss_op->bp[2]);
    if(qss_op->bp[3] == 1)
      fprintf(fp,"%s%s\n", make_tag(prefix, "time_bc"), "antiperiodic");
    else
      fprintf(fp,"%s%s\n", make_tag(prefix, "time_bc"), "periodic");
      
  }
  else if( op_type == FAT_COVARIANT_GAUSSIAN){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d,\n", make_tag(prefix, "stride"), qss_op->stride);
    fprintf(fp,"%s%g,\n", make_tag(prefix, "r0"), qss_op->r0);
    fprintf(fp,"%s%d\n", make_tag(prefix, "iters"), qss_op->iters);
  }
  else if( op_type == FAT_COVARIANT_LAPLACIAN){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d\n", make_tag(prefix, "stride"), qss_op->stride);
  }
  else if ( op_type == GAUSSIAN ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g\n", make_tag(prefix, "r0"), qss_op->r0);
  }
  else if ( op_type == HOPPING ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d,\n", make_tag(prefix, "derivs"), qss_op->dhop);
#if FERM_ACTION == HISQ
    fprintf(fp,"%s%s,\n", make_tag(prefix, "dir"), 
	    encode_sign_dir(qss_op->fb, qss_op->dir1));
    fprintf(fp,"%s%g\n", make_tag(prefix, "eps_naik"), qss_op->eps_naik);
#else
    fprintf(fp,"%s%s\n", make_tag(prefix, "dir"), encode_dir(qss_op->dir1));
#endif
  }
  else if( op_type == KS_INVERSE){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "mass"), qss_op->mass_label);
    fprintf(fp,"%s%g,\n", make_tag(prefix, "eps_naik"), qss_op->eps_naik);
    fprintf(fp,"%s[%d, %d, %d, %d],\n", make_tag(prefix, "coordinate_origin"), 
	    qss_op->co[0], qss_op->co[1], qss_op->co[2], qss_op->co[3]);
    fprintf(fp,"%s[%f, %f, %f],\n", make_tag(prefix, "momentum_twist"), 
	    qss_op->bp[0], qss_op->bp[1], qss_op->bp[2]);
    if(qss_op->bp[3] == 1)
      fprintf(fp,"%s%s\n", make_tag(prefix, "time_bc"), "antiperiodic");
    else
      fprintf(fp,"%s%s\n", make_tag(prefix, "time_bc"), "periodic");
      
  }
  else if ( op_type == ROTATE_3D ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%g\n", make_tag(prefix, "d1"), qss_op->d1);
  }
#ifdef HAVE_KS
  else if ( op_type == SPIN_TASTE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s\n", make_tag(prefix, "spin_taste"), 
	    spin_taste_label(qss_op->spin_taste));
  }
  else if ( op_type == SPIN_TASTE_EXTEND ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s\n", make_tag(prefix, "spin_taste_extend"), 
	    spin_taste_label(qss_op->spin_taste));
  }
  else if ( op_type == EXT_SRC_KS ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "spin_taste_extend"), 
	    spin_taste_label(qss_op->spin_taste));
    fprintf(fp,"%s[%d, %d, %d],\n", make_tag(prefix, "mom"), qss_op->mom[0],
	    qss_op->mom[1], qss_op->mom[2]);
    fprintf(fp,"%s%d\n", make_tag(prefix, "t0"), qss_op->t0);
  }
  else if ( op_type == ASLASH_KS_FILE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "file"), qss_op->source_file);
    fprintf(fp,"%s%d\n", make_tag(prefix, "t0"), qss_op->t0);
  }
#endif
  else if ( op_type == WAVEFUNCTION_FILE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "file"), qss_op->source_file);
    fprintf(fp,"%s%d,\n", make_tag(prefix, "stride"), qss_op->stride);
    fprintf(fp,"%s%g\n", make_tag(prefix, "a"), qss_op->a);
  }
  else if ( op_type == MODULATION_FILE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "file"), qss_op->source_file);
  }
  else if ( op_type == MOMENTUM ){
    fprintf(fp,",\n");
    fprintf(fp,"%s[%d, %d, %d]\n", make_tag(prefix, "mom"), qss_op->mom[0],
            qss_op->mom[1], qss_op->mom[2]);
  }
  else if ( op_type == PROJECT_T_SLICE ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d\n", make_tag(prefix, "t0"), qss_op->t0);
  }
  else if ( op_type == GAMMA ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s\n", make_tag(prefix, "gamma"), 
	    gamma_label(qss_op->gamma));
  }
  else if ( op_type == MOMENTUM ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%d %d %d\n", make_tag(prefix, "momentum"), 
	    qss_op->mom[0], qss_op->mom[1], qss_op->mom[2]);
  }
  else if ( op_type == EXT_SRC_DIRAC ){
    fprintf(fp,",\n");
    fprintf(fp,"%s%s,\n", make_tag(prefix, "gamma"), 
	    gamma_label(qss_op->gamma));
    fprintf(fp,"%s[%d, %d, %d],\n", make_tag(prefix, "mom"), qss_op->mom[0],
	    qss_op->mom[1], qss_op->mom[2]);
    fprintf(fp,"%s%d\n", make_tag(prefix, "t0"), qss_op->t0);
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
    if(qss_op[i]->type == IDENTITY)
      continue;
    fprintf(fp,"{ operation:                  %s",qss_op[i]->descrp);
    print_single_op_info(fp, "  ", qss_op[i]);
    if(i < n-1)fprintf(fp, "},\n");
    else fprintf(fp, "}\n");
  }

  fprintf(fp, "]\n");
}


/* quark_source_sink_op.c */
