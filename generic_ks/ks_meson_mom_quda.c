/********************** ks_meson_mom.c *********************************/
/* MIMD version 7 */

/* Reorganized from ks_meson_mom.c to support a call to QUDA */

/*  Tie together staggered propagators to make meson correlators
 *
 */

/*****************************************

This is a general function that can contract two quark
propagators together to form a meson correlator.

\sum_{x} exp( i p .x } 
     [ O_st src_1 ]^\dagger src_2  )

   where src_1 and src_2 are quark propagators and
   O_st is a spin-taste operator at the sink.

   Note that, unlike the analogous Dirac propagator routine, the
   source spin-taste can not be applied independently of the
   propagator calculation, so it is assumed already to be included in
   the propagators src_1 and/or src_2. Ordinarily, when converting a
   quark propagator running backwards in time to an antiquark running
   forwards in time, we take the adjoint and apply (-)^(x+y+z+t)
   factors at source and sink.  This sink sign factor is applied in
   the call to spin_taste_op_ape_fn, so we don't do it explicitly here.
   This means, for example, that when we specify O_st = pion5 we get
   simply src_1^\dagger src_2.


  Function arguments

     On input 
        src1 :: su3_vector 
        src2 :: su3_vector 
        no_q_momenta :: number of momentum values p in table to use
	q_momstore[p] :: table of momentum values p to use
        q_parity[p] :: reflection parity for each momentum component
        no_spin_taste_corr :: Number of unique spin-taste assignments (gt g)
        num_corr_mom[g] :: number of momentum/parity sets for each g
        corr_table[g] :: list of correlator indices c for gamma pair g
        p_index[c] :: p = p_index[c] is the momentum index for correlator c 
        fn_src1 :: Fat-Naik links for src1 needed for some spin-taste operators 
        fn_src2 :: Fat-Naik links for src2 Needed for some spin-taste operators
        meson_phase[c] :: phase factor to multiply the correlator before
                          accumulating.  Encoded as in gammatypes.h
        meson_factor[c] :: normalization factor for each correlator c
        corr_index[c] :: correlator index m where like propagators are summed
        r0[c] :: origin for defining FT and KS phases for correlator c

     On output
         prop :: complex vector to the data correlators with indexing
                 prop[m][time] where m is the correlator index.


*******************************************/

#include "generic_ks_includes.h"
#include <string.h>
#include "../include/openmp_defs.h"
#ifdef OMP
#include <omp.h>
#endif
#include "../include/static_cast.h"


/*******************************************/
// Uncomment the following include, then definitions below can be made into comments or removed
//#include <quda_milc_interface.h>

#include <limits.h>
#define QUDA_INVALID_ENUM INT_MIN

// enum_quda.h  describes corr_parity
typedef enum QudaFFTSymmType_t {
  QUDA_FFT_SYMM_ODD  = 1,  // sin(phase)
  QUDA_FFT_SYMM_EVEN = 2,  // cos(phase)
  QUDA_FFT_SYMM_EO   = 3,  // exp(-i phase)
  QUDA_FFT_SYMM_INVALID = QUDA_INVALID_ENUM
} QudaFFTSymmType;

// quda_milc_interface.h  parameters for propagator contractions with FT
typedef struct {
  int n_mom;  /* Number of sink momenta */
  int *mom_modes;  /* List of 4-component momenta as integers. Dimension 4*n_mom */
  QudaFFTSymmType *fft_type; /* The "parity" of the FT component */
  int *source_position;  /* The coordinate origin for the Fourier phases */
  double flops; /* Return value */
  double dtime; /* Return value */
} QudaContractArgs_t;

/** quda_milc_interface.h
 * @brief Tie together two staggered propagators including spatial Fourier phases.
 * The result is summed separately over each time slice and across all MPI ranks.
 * The FT is defined by a list of momentum indices (three-component integer vectors)
 * Included with the FT is a parity (symmetry) parameter for each momentum
 * component that selects an exp, cos, or sin factor for each direction
 *
 * @param[in] external_precision Precision of host fields passed to QUDA (2 - double, 1 - single)
 * @param[in,out] parameters for the contraction, including FT specification
 * @param[in] local storage of color spinor field.  three complex values * number of sites on node
 * @param[in] local storage of color spinor field.  three complex values * number of sites on node
 * @param[out] hadron correlator  Flattened double array as though [n_mom][nt][2] for 2 = re,im. 
 */
void qudaContractFT(int external_precision, // milc_precision
		    QudaContractArgs_t *cont_args,
		    void *const quark1, // su3_vector* antiquark
		    void *const quark2, // su3_vector* quark
		    double *corr // double_complex meson_q[]
		    );


/*******************************************/

void dump_QudaContractArgs(const QudaContractArgs_t * args)
{
  printf("QudaContractArgs:\n");
  printf("source_position: %3d %3d %3d %3d\n",
	 args->source_position[0],args->source_position[1],
	 args->source_position[2],args->source_position[3]);
  printf("n_mom: %d\n",args->n_mom);
  for(int k=0; k<args->n_mom; ++k) {
    printf("mom[%2d]: %2d %2d %2d %2d\n",k,
	   args->mom_modes[4*k+0],args->mom_modes[4*k+1],
	   args->mom_modes[4*k+2],args->mom_modes[4*k+3]);
    printf("sym[%2d]: %2d %2d %2d %2d\n",k,
	   args->fft_type[4*k+0],args->fft_type[4*k+1],
	   args->fft_type[4*k+2],args->fft_type[4*k+3]);
  }
}

/* Normalize the correlator contributions */

static double norm_v(double_complex tr[], double_complex meson_q[], 
		     int phase[], Real factor[],
		     int ct[], int nk, int nt)
{
  double_complex z = {0.,0.};
  double flops = 0;
  
  /* For each momentum in list, normalize, and phase */

  for(int k=0; k<nk; k++){
    int c = ct[k];
    int ph = phase[c];
    Real fact = factor[c];
    for(int t=0; t<nt; ++t)
      {
	int idx = k*nt + t;
	tr[idx] = meson_q[idx];
	switch(ph){
	case 0:
	  z =            tr[idx];
	  break;
	case 1:
	  TIMESPLUSI(    tr[idx], z);
	  break;
	case 2:
	  TIMESMINUSONE( tr[idx], z);
	  break;
	case 3:
	  TIMESMINUSI(   tr[idx], z);
	}
	CMULREAL(z,fact,tr[idx]);
      }
  }

  flops = 2*nk*nt;
  
  return flops;
  
} /* norm_v */

/*******************************************/
static double_complex *
create_meson_q(int nt, int num_corr_mom){
  char myname[] = "create_meson_q";

  /* Unlike the CPU version, meson_q here is indexed by the 
     actual momenta, rather than the hashed momentum. */
  double_complex* meson_q = static_cast(double_complex*,malloc(num_corr_mom*nt*sizeof(double_complex))); // index as meson_q[k*nt+t]
  if(meson_q == NULL){
    printf("%s(%d): No room for meson_q\n",myname,this_node);
    terminate(1);
  }
  
  for(int j=0; j<nt*num_corr_mom; j++)
    {   
      meson_q[j].real = 0.;
      meson_q[j].imag = 0.;
    }
  return meson_q;
}

/*******************************************/
static void
destroy_meson_q(double_complex *meson_q){
  if(meson_q == NULL)return;
  free(meson_q);
}
/*******************************************/
static Real
update_props(complex **prop, double_complex *meson_q, int nt, int num_corr_mom,
	     int meson_phase[], Real meson_factor[],
	     int *corr_table, int corr_index[]){

  Real flops = 0.;

  double_complex tr[num_corr_mom*nt];
  /* Normalize for all sink momenta q */
  flops += norm_v(tr, meson_q, meson_phase, meson_factor,
		  corr_table, num_corr_mom, nt);
  /* Accumulate in corr_index location */
  for(int k=0; k<num_corr_mom; k++)
    for(int t=0; t<nt; ++t)
      {
	int c = corr_table[k];
	int m = corr_index[c];
	int idx = k*nt + t;
	prop[m][t].real += tr[idx].real;
	prop[m][t].imag += tr[idx].imag;
      }

  flops += 2. * num_corr_mom * nt;
  
  return flops;
}

/*******************************************/

static QudaFFTSymmType
map_to_fft_symm(char p)
{
  QudaFFTSymmType s = QUDA_INVALID_ENUM;
  switch(p)
    {
    case EVEN:
      s = QUDA_FFT_SYMM_EVEN; break;
    case ODD:
      s = QUDA_FFT_SYMM_ODD; break;
    case EVENANDODD:
      s = QUDA_FFT_SYMM_EO; break;
    }
  return s;
}

static void
map_corr_mom_parity(int corr_mom[], QudaFFTSymmType corr_parity[], int num_corr_mom, int corr_table[],
		    int p_index[], int **q_momstore, char **q_parity){

  for(int k=0; k<num_corr_mom; k++) {
    int c = corr_table[k];
    int p = p_index[c];
    //corr_mom[k] = q_momstore[p];
    for(int j=0; j<3; ++j) { corr_mom[4*k+j] = q_momstore[p][j]; }
    corr_mom[4*k+3] = 0; // t-component: always zero
    //corr_parity[k] = q_parity[p];
    for(int j=0; j<3; ++j) { corr_parity[4*k+j] = map_to_fft_symm(q_parity[p][j]); }
    corr_parity[4*k+3] = QUDA_FFT_SYMM_EO; // t-component: pt = 0, but quda contractions allows pt
  }
}

/*******************************************/
void ks_meson_cont_mom(
  complex **prop,           /* prop[m][t] is where result is accumulated */
  su3_vector *src1,         /* quark propagator (to become antiquark) */
  su3_vector *src2,         /* quark propagator */
  int no_q_momenta,         /* no of unique mom/parity values (gt p) */
  int **q_momstore,         /* q_momstore[p] are the momentum components */
  char **q_parity,          /* q_parity[p] the parity of each mom component */
  int no_spin_taste_corr,   /* Number of unique spin-taste assignments (gt g) */
  int num_corr_mom[],       /* number of momentum/parity values (gt k)*/
  int **corr_table,         /* c = corr_table[k] correlator index */
  int p_index[],            /* p = p_index[c] is the momentum index */
  imp_ferm_links_t *fn_src1, /* Needed for some spin-taste operators */
  imp_ferm_links_t *fn_src2, /* Needed for some spin-taste operators */
  int spin_taste_snk[],     /* spin_taste_snk[c] gives the s/t assignment */
  int meson_phase[],        /* meson_phase[c] is the correlator phase */
  Real meson_factor[],      /* meson_factor[c] scales the correlator */
  int num_corr,             /* number of corrs - first index of prop */
  int corr_index[],         /* m = corr_index[c] is the correlator index */
  int r0[]                  /* origin for defining FT and KS phases */
		    )
{
  char myname[] = "meson_cont_mom";
  int i,k,c,m;
  site *s; 
  
  int spin_taste;
  int g,p,t;
  int max_threads,mythread;
  complex tmp;
  
  /* performance */
  double dtime;
  double flops;
#ifdef WMTIME
  double mflops;
#endif
  dtime = -dclock();
  flops = 0;

  /* Run through the sink spin-taste combinations */
  for(g = 0; g < no_spin_taste_corr; g++)
    {
      double_complex *meson_q = create_meson_q(nt, num_corr_mom[g]);
      double *dmeson_q = static_cast(double*,meson_q);

      /* Transfer momenta and parity from q_momstore to corr_mom array */
      int corr_mom[num_corr_mom[g]*4]; // all four components of the momentum
      QudaFFTSymmType corr_parity[num_corr_mom[g]*4];
      map_corr_mom_parity(corr_mom, corr_parity, num_corr_mom[g],
			  corr_table[g], p_index, q_momstore, q_parity);
      /* Load the contraction args */
      QudaContractArgs_t cont_args;
      
      cont_args.n_mom = num_corr_mom[g];
      cont_args.mom_modes = corr_mom;
      cont_args.fft_type = corr_parity;
      cont_args.source_position = r0;
      cont_args.flops = 0;
      cont_args.dtime = 0;
      //dump_QudaContractArgs(&cont_args);

      /* Apply spin tastes and call the GPU contraction routine */

      /* All spin-taste assignments with the same index g must be the same */
      c = corr_table[g][0];  
      spin_taste = spin_taste_snk[c];

      /* Handle the various parallel transport option */
      if(is_rhosfn_index(spin_taste) || is_rhosape_index(spin_taste)){
	/* Special treatment for vector-current operators */
	/* Apply backward sink spin-taste operator to src1 */

	su3_vector *q = create_v_field();
	spin_taste_op_ape_fn(fn_src1, backward_index(spin_taste), r0, q, src1);
	qudaContractFT(MILC_PRECISION, &cont_args, q, src2, dmeson_q);
	/* Apply forward sink spin-taste operator to src2 */
	spin_taste_op_ape_fn(fn_src2, forward_index(spin_taste), r0, q, src2);
	qudaContractFT(MILC_PRECISION, &cont_args, src1, q, dmeson_q);
	destroy_v_field(q);
	for(int j = 0; j < nt*num_corr_mom[g]; ++j){
	  CMULREAL(meson_q[j], 0.5, meson_q[j]);
	}
      } else if(is_rhosffn_index(spin_taste) || is_rhosfape_index(spin_taste)){
	/* Apply forward sink spin-taste operator to src2 */
	su3_vector *q = create_v_field();
	spin_taste_op_ape_fn(fn_src2, forward_index(spin_taste), r0, q, src2);
	qudaContractFT(MILC_PRECISION, &cont_args, src1, q, dmeson_q );
	destroy_v_field(q);
      } else if(is_rhosbfn_index(spin_taste) || is_rhosbape_index(spin_taste)){
	/* Apply backward sink spin-taste operator to src1 */
	su3_vector *q = create_v_field();
	spin_taste_op_ape_fn(fn_src1, backward_index(spin_taste), r0, q, src1);
	qudaContractFT(MILC_PRECISION, &cont_args, q, src2, dmeson_q);
	destroy_v_field(q);
      } else {
	/* Apply sink spin-taste operator to src1 */
	su3_vector *q = create_v_field();
	spin_taste_op_ape_fn(fn_src1, spin_taste, r0, q, src1);
	qudaContractFT(MILC_PRECISION, &cont_args, q, src2, dmeson_q);
	destroy_v_field(q);
       }
      flops += update_props(prop, meson_q, nt, num_corr_mom[g], meson_phase,
			    meson_factor, corr_table[g], corr_index);

      destroy_meson_q(meson_q);
    }  /**** end of the loop over the spin-taste table ******/
  
  dtime += dclock();

/**  flops = sites_on_node*(nmeson_evals*1536 + nmeson_q_evals*128) +
     68*no_q_momenta) + 14*nprop_evals*ntslices*no_q_momenta; **/

#ifdef WMTIME
  if(dtime > 0)mflops = flops/(dtime*1e6);
  else mflops = 0;

  node0_printf("WMTIME: time %.1e sec %g flops %.1f MF\n",
	       dtime,flops,mflops);fflush(stdout);
#endif
  
} /* end of ks_meson_mom_gpu function  */

