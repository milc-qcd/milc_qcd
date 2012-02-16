/******************************** flavor_ops2.c *************************/
/* MIMD version 7                                                       */

/* Implementation of the spin-taste ("flavor") (\Xi_\mu) operators.
   See Golterman & Smit Nulc. Phys. B245 (1984) 61.
   They are used for constructing the non-local pion sources
   and in the measurement of all the componenets of \bar{\psi}\psi

   CD 07/27/09 This version works with site-major fields
   Ths version also includes a set of single time-slice meson operators.
*/

#include "generic_ks_includes.h"
#include <string.h>

#ifndef NO_GAUGE_FIELD

/* This is the definition of the epsilon symbol */
static struct {
  int d[4];
  Real sign ;
} eps[24]={
  {{0,1,2,3},+1.0},
  {{0,3,1,2},+1.0},
  {{0,2,3,1},+1.0},
  {{0,3,2,1},-1.0},
  {{0,1,3,2},-1.0},
  {{0,2,1,3},-1.0},

  {{1,0,2,3},-1.0},
  {{1,3,0,2},-1.0},
  {{1,2,3,0},-1.0},
  {{1,3,2,0},+1.0},
  {{1,0,3,2},+1.0},
  {{1,2,0,3},+1.0},

  {{2,1,0,3},-1.0},
  {{2,3,1,0},-1.0},
  {{2,0,3,1},-1.0},
  {{2,3,0,1},+1.0},
  {{2,1,3,0},+1.0},
  {{2,0,1,3},+1.0},

  {{3,1,2,0},-1.0},
  {{3,0,1,2},-1.0},
  {{3,2,0,1},-1.0},
  {{3,0,2,1},+1.0},
  {{3,1,0,2},+1.0},
  {{3,2,1,0},+1.0}
};

#endif

/* Compute the hypercube coordinate relative to an offset.  We assume
   that all lattice dimensions are even, as they should be for
   staggered fermions! */
static short *
hyp_coord(site *s, int r0[]){
  static short h[4];
  h[XUP] = (s->x - r0[XUP]) & 0x1;
  h[YUP] = (s->y - r0[YUP]) & 0x1;
  h[ZUP] = (s->z - r0[ZUP]) & 0x1;
  h[TUP] = (s->t - r0[TUP]) & 0x1;
  return h;
}

/* Compute the parity of the site relative to an offset. */
static short 
hyp_parity(site *s, int r0[]){
  short p;
  if((r0[XUP] + r0[YUP] + r0[ZUP] + r0[TUP]) % 2 == 0)
    p = EVEN;
  else p = ODD;
  if(p == s->parity)
    return EVEN;
  else
    return ODD;
}

#ifndef NO_GAUGE_FIELD

/* Apply the symmetric shift operator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              *
 * The KS phases MUST BE in the links                     */
static void 
sym_shift_field(int dir, su3_vector *src, su3_vector *dest)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec = create_v_field();

  tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
  /* With ONE_SIDED_SHIFT defined, the shift is asymmetric */
#ifndef ONE_SIDED_SHIFT
  FORALLSITES(i,s)
    {
      mult_adj_su3_mat_vec( &(s->link[dir]), src+i, tvec+i ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
			      EVENANDODD ,gen_pt[1] );
#endif
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
      mult_su3_mat_vec( &(s->link[dir]), (su3_vector *)gen_pt[0][i], dest+i );
    }
#ifndef ONE_SIDED_SHIFT
  wait_gather(tag[1]);
  FORALLSITES(i,s)
    {
      add_su3_vector( dest+i, (su3_vector *)gen_pt[1][i], dest+i ) ;    
    }
  /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
  FORALLSITES(i,s)
    {
      scalar_mult_su3_vector( dest+i, .5, dest+i );
    }
#endif
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
  
  destroy_v_field(tvec);
}

/* Apply the symmetric shift with directions                            *
 * stored in the array d. Each shift is multiplied by \zeta_k           *
 * n is the number of shifts                                            *
 * This is the E_\mu(x,y)=\Xi_\mu operator defined by Golterman.        *
 * Nucl. Phys. B245  eq.3.5 and eq. 4.2b                                */
static void 
zeta_shift_field(int n, int *d, int r0[], su3_vector *src, 
		 su3_vector *dest )
{
  register int i,c ;
  register site *s;
  register Real sign ;
  short *h;
  su3_vector *tvec = create_v_field();
  
  
  for(c=0;c<n;c++)
    {  
      /* Do the shift in d[c] */ 
      if(c==0)
	/* first time from source */
	sym_shift_field(d[c], src, tvec);
      else
	/* second time from dest */
	sym_shift_field(d[c], dest, tvec);
      /* Multiply by \zeta_d[c]. Because the phases are               *
       * on we multiply by \zeta * \eta = \epsilon * (-1)^coord[d[c]] */
      FORALLSITES(i,s){
	h = hyp_coord(s, r0);
	/* The \epsilon */
	if(hyp_parity(s, r0)==EVEN) 
	  sign =  1.0 ;
	else
	  sign = -1.0 ;
	/* And the (-1)^coord[d[c]] */
	if(h[d[c]]==1) sign=-sign ;
	scalar_mult_su3_vector(tvec+i, sign, dest+i );
      }
    }
  destroy_v_field(tvec);
}

/* Apply the symmetric shift with directions                            *
 * stored in the array d. Each shift is multiplied by \eta_k            *
 * In fact since \eta_k are already absorbed into the U matrices        *
 * no action is needed. Just do the symmetric shift.                    *
 * n is the number of shifts                                            */
static void 
eta_shift_field(int n, int *d, su3_vector *src, su3_vector *dest )
{
  int c ;
  su3_vector *temp0 = create_v_field();
  su3_vector *temp1 = create_v_field();
  su3_vector *tmp;
  
  if(n==1)
    {
      sym_shift_field(d[0], src, dest);
    }
  else
    {
      sym_shift_field(d[0], src, temp0 );
      for(c=1;c<n-1;c++)
	{  
	  sym_shift_field(d[c], temp0, temp1);
	  /* switch the pointers */
	  tmp = temp0   ;
	  temp0 = temp1 ;
	  temp1 = tmp   ;
	}
      /* do the last shift */
      sym_shift_field(d[n-1], temp0, dest );
    }
  destroy_v_field(temp0);
  destroy_v_field(temp1);
}

/* Multiply by the Xi_mu taste operator */
void 
mult_flavor_vector_field(int mu, int r0[],
			 su3_vector *src, su3_vector *dest )
{
  int d[1] ;
  
  d[0] = mu ;
  zeta_shift_field(1,d,r0,src,dest) ;
}

/* Multiply by the 1/2(Xi_mu Xi_nu - Xi_nu Xi_mu)taste operator */
void 
mult_flavor_tensor_field(int mu, int nu, int r0[],
			 su3_vector *src, su3_vector *dest )
{ 
  register int i;
  register site *s;
  int d[2];
  su3_vector *tvec = create_v_field();

  d[0] = mu ; d[1]=nu ;
  zeta_shift_field(2, d, r0, src, dest) ;
 
  d[0] = nu ; d[1]=mu ;
  zeta_shift_field(2, d, r0, src, tvec) ;

  
  FORALLSITES(i,s){
    sub_su3_vector(dest+i, tvec+i, dest+i);
    scalar_mult_su3_vector(dest+i, 0.5, dest+i );
  }
  destroy_v_field(tvec);
}

/* Multiply by the Xi_mu Xi_5 taste operator */
void 
mult_flavor_pseudovector_field(int mu, int r0[],
			       su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;
  int p ; 
  su3_vector *tvec = create_v_field();

  clear_v_field(dest);
  
  for(p=0;p<24;p++)
    if(eps[p].d[0]==mu)
      {
	zeta_shift_field(3,&eps[p].d[1],r0,src,tvec) ;
	/* Multiply the extra 1/6 needed by the definition    *
	 * of the operator (number of permutations)           */
	FORALLSITES(i,s){
	  scalar_mult_sum_su3_vector(dest+i, tvec+i, eps[p].sign/6.0);
	}
      }
  destroy_v_field(tvec);
}

/* Multiply by the Xi_5 taste operator */
void 
mult_flavor_pseudoscalar_field(int r0[],
			       su3_vector *src, su3_vector *dest )
{ 
  register int i;
  register site *s;
  int p ; 
  su3_vector *tvec = create_v_field();

  clear_v_field(dest);

  for(p=0;p<24;p++)
     {
       zeta_shift_field(4,eps[p].d,r0,src,tvec) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector(dest+i, tvec+i, eps[p].sign/24.0);
       }
     } 
  destroy_v_field(tvec);
}

/* Multiply by the Gamma_mu spin operator */
void 
mult_spin_vector_field(int mu, int r0[], 
		       su3_vector *src, su3_vector *dest )
{
  int d[1] ;
  
  d[0] = mu ;
  eta_shift_field(1,d,src,dest) ;
}

/* Multiply by the 1/2(Gamma_mu Gamma_nu - Gamma_nu gamma_mu) spin operator */
void 
mult_spin_tensor_field(int mu, int nu, int r0[],
		       su3_vector *src, su3_vector *dest )
{ 
  register int i;
  register site *s;
  int d[2];
  su3_vector *tvec = create_v_field();
  
  d[0] = mu ; d[1]=nu ;
  eta_shift_field(2, d, src, dest) ;
 
  d[0] = nu ; d[1]=mu ;
  eta_shift_field(2, d, src, tvec);

  FORALLSITES(i,s){
    sub_su3_vector(dest+i, tvec+i, dest+i);
    scalar_mult_su3_vector(dest+i, 0.5, dest+i);
  }
  destroy_v_field(tvec);
}

/* Multiply by the Gamma_mu Gamma_5 spin operator */
void 
mult_spin_pseudovector_field(int mu, int r0[],
			     su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;
  int p ; 
  su3_vector *tvec = create_v_field();

  /*clean up dest */
  clear_v_field(dest);

  for(p=0;p<24;p++)
    if(eps[p].d[0]==mu)
      {
	eta_shift_field(3,&eps[p].d[1],src,tvec);
	/* Multiply the extra 1/6 needed by the definition    *
	   of the operator (number of permutations)           */
	FORALLSITES(i,s){
	  scalar_mult_sum_su3_vector(dest+i, tvec+i, eps[p].sign/6.0);
	}
      }
  destroy_v_field(tvec);
}

/* Multiply by the Gamma_5 spin operator */
void 
mult_spin_pseudoscalar_field(int r0[], 
			     su3_vector *src, su3_vector *dest )
{ 
  register int i;
  register site *s;
  int p ; 
  su3_vector *tvec = create_v_field();
  
  /*clean up dest */
  clear_v_field(dest);

  for(p=0;p<24;p++)
     {
       eta_shift_field(4,eps[p].d,src,tvec) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector(dest+i, tvec+i, eps[p].sign/24.0);
       }
     } 
  destroy_v_field(tvec);
}

#endif

/*------------------------------------------------------------------*/

/* "Multiply by" the quark-antiquark local pion operator */
/* Here the operator is gamma_5 x gamma_5                */
void 
mult_pion5_field( int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is (-1)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark,
     so nothing to do */
  
  copy_v_field(dest, src);
}

/* "Multiply by" the second quark-antiquark local pion operator */
/* Here the operator is 1 x 1 times (-1)^(x+y+z+t) (a scalar).
   It has the other parity state gamma_0 gamma_5 x gamma_0 gamma_5
   (a pion) -CD */
void 
mult_pion05_field( int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is (1), another (-1)^(x+y+z+t) for antiquark, 
     so multiply by (-1)^(x+y+z+t) */
  register int i;
  register site *s;
  FORALLSITES(i,s){
    if( hyp_parity(s, r0) == EVEN ){
      dest[i] = src[i];
    } else {
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
    }
  }
}

#ifndef NO_GAUGE_FIELD

/* "Multiply by" the one link pion operator.
   "pi_i5" or gamma_5 x gamma_i gamma_5 */
/* CD Are pionij and pioni5 switched here? */
void 
mult_pioni5_field( int fdir, int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is (-)^fdir, another (-1)^(x+y+z+t) for antiquark, 
     so multiply by (-1)^(x+y+z+t) (-)^fdir */
  register int i;
  register site *s;
  Real sign;
  
  /* apply the symmetric shift operator */
  sym_shift_field(fdir, src, dest);
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    sign = +1.0;
    if( h[fdir]==1 ) sign = -sign;
    if( hyp_parity(s, r0)==ODD ) sign = -sign; /* the extra gamma_5 */
    scalar_mult_su3_vector( dest+i, sign, dest+i );
  }
}

/* "Multiply by" the one link pion operator.
   fdir is 0 for ij = 12, 1 for ij = 20, 2 for ij = 01.
   "pi_ij" or gamma_0 gamma_5 x gamma_i gamma_j */
/* CD Are pionij and pioni5 switched here? */
void 
mult_pionij_field( int fdir, int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is (-)^fdir (-)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark, 
     so multiply by (-1)^fdir */
  register int i;
  register site *s;
  Real sign;
  
  /* apply the symmetric shift operator */
  sym_shift_field(fdir, src, dest);
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    sign = +1.0;
    if( h[fdir]==1) sign = -sign;
    /* two gamma_5's = nothing */
    scalar_mult_su3_vector( dest+i, sign, dest+i );
  }
}

/* "Multiply by" the two link pion operator. *
 * "pi_i or gamma_0 gamma_5 x gammma_i"                                    */
void 
mult_pioni_field( int fdir, int r0[], su3_vector *src, su3_vector *dest )
{
  register int i,d1,d2;
  register site *s;
  Real sign, epsilon_sign;
  short dir[2] ;
  su3_vector *tvec0 = create_v_field();
  su3_vector *tvec1 = create_v_field();
  
  switch(fdir)
    {
    case 0: dir[0]=1 ; dir[1] = 2 ; break ;
    case 1: dir[0]=2 ; dir[1] = 0 ; break ;
    case 2: dir[0]=0 ; dir[1] = 1 ; break ;
    default:   node0_printf("ERROR! invalid direction %i\n",fdir);
    }
  
  /*clean up dest */
  /**/
  clear_v_field(dest);
  
  /**/
  epsilon_sign=.5 ; /* divide by two for averaging (1,2), (2,1) paths */
  for(d1=0,d2=1;d1<2;d1++,d2--)
    {
      /*      printf("In mult_pioni: (d1,d2): (%i,%i)\n",dir[d1],dir[d2]); */
      /* apply the symmetric shift operator in dir2 */
      sym_shift_field(dir[d2], src, tvec0);
      /* multiply by \zeta_dir2 */
      FORALLSITES(i,s){
	short *h = hyp_coord(s, r0);
	/*   because the phases are on we multiply by  * 
	 *  \zeta * \eta = \epsilon * (-1)^coord[dir2] */
	if( hyp_parity(s, r0)==EVEN ) 
	  sign =  1.0 ;
	else
	  sign = -1.0 ;
	if( h[dir[d2]]==1 ) sign = -sign ;
	scalar_mult_su3_vector( tvec0+i, sign, tvec0+i );
      }
      /* apply the symmetric shift operator in dir1 */
      sym_shift_field(dir[d1], tvec0, tvec1);
      /* multiply by \zeta_dir1 */
      FORALLSITES(i,s){
	short *h = hyp_coord(s, r0);
	/*   because the phases are on we multiply by  * 
	 *  \zeta * \eta = \epsilon * (-1)^coord[dir1] *
	 *  Here we have to multiply with the extra    *
	 *  \epsilon for the anti-quark propagator.    *
	 * so we are left with only  (-1)^coord[dir1]  */
	sign=1.0 ;
	if(h[dir[d1]]==1) sign = -sign ;
	sign *= epsilon_sign ;
	scalar_mult_sum_su3_vector( dest+i, tvec1+i, sign );
      }
      epsilon_sign = -epsilon_sign ;
    }
  destroy_v_field(tvec0);
  destroy_v_field(tvec1);
}


/* "Multiply by" the two link pion operator. */
/* pi_i0 or gamma_5 x gamma_0 gamma_i */
void 
mult_pioni0_field( int fdir, int r0[], su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;
  
  /* We only need to multiply the pi_4 by \zeta_4\eta_4=(-1)^(x+y+z) */
  mult_pioni_field( fdir, r0, src, dest ) ;
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    if( (h[XUP] + h[YUP] + h[ZUP])%2 == 1 )
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
}


/* "Multiply by" the three link pion operator.  */
/* pion_s or gamma_0 gamma_5 x 1 */
void 
mult_pions_field(int r0[], su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;
  int c ;
  su3_vector *tvec = create_v_field();
  
  /* All permutation with appropriate sign */
  struct {
    int d[3];
    Real sign ;
  } p[6]={{{0,1,2},+1.0/6.0},
	  {{1,2,0},+1.0/6.0},
	  {{2,0,1},+1.0/6.0},
	  {{0,2,1},-1.0/6.0},
	  {{1,0,2},-1.0/6.0},
	  {{2,1,0},-1.0/6.0}}; /* The factor of 6 accounts for the *
				* multiplicity of the permutations */
  
  /*clean up dest */
  clear_v_field(dest);
  
  for(c=0;c<6;c++)
    {
      zeta_shift_field(3,p[c].d,r0,src,tvec);
      FORALLSITES(i,s){
	scalar_mult_sum_su3_vector(dest+i, tvec+i, p[c].sign );
      }
    }
  /* multiply by \epsilon for the anti-quark */
  FORALLSITES(i,s){
    if( hyp_parity(s,r0) == ODD )
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
  
  destroy_v_field(tvec);
}

/* "Multiply by" the three link pion operator.  */
/* gamma_5 x gamma_0 */
void 
mult_pion0_field(int r0[], su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;
  
  /* We only need to multiply the pi_5 by \zeta_4\eta_4=(-1)^(x+y+z) */
  mult_pions_field( r0, src, dest ) ;
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    if( (h[XUP] + h[YUP] + h[ZUP]) % 2 == 1 )
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
}

#endif

/* "Multiply by" the quark-antiquark local rho operator */
/* gamma_i x gamma_i */
void 
mult_rhoi_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
  register int i;
  register site *s;
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    //    if( ((((short *)&(s->x))[pdir]) & 0x1) == 0 ){
    if(h[pdir] == 0){
      dest[i] = src[i];
    } else {
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
    }
  }
}

/* "Multiply by" the second quark-antiquark local rho operator */
void 
mult_rhoi0_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest ){
  /* gamma_0 gamma_i x gamma_0 gamma_i */
  /* operator is gamma_0 gamma_pdir, another (-1)^(x+y+z+t) for antiquark */
  register int i;
  register site *s;
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    // if( ((((short *)&(s->x))[pdir] + s->x+s->y+s->z) & 0x1) == 0 ){
    if( (h[pdir] + h[XUP] + h[YUP] + h[ZUP]) % 2 == 0){
      dest[i] = src[i];
    } else {
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
    }
  }
}

#ifndef NO_GAUGE_FIELD

/* "Multiply by" the quark-antiquark one link rho operator */
/* This is the one-link conserved vector current for the naive action -CD */
/* It is gamma_i x 1 times (-1)^(x+y+z+t) for antiquark */
void 
mult_rhois_field( int fdir,  int r0[], su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;  
  
  /* apply the symmetric shift operator */
  sym_shift_field(fdir, src, dest);
  FORALLSITES(i,s){
    /* \eta_k already in the phases                   * 
     * only the \epsilon for the anti-quark is needed */
    if( hyp_parity(s, r0)==ODD )
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
}

/* "Multiply by" the second quark-antiquark one link rho operator */
/* gamma_i gamma_0 x gamma_0  times (-1)^(x+y+z+t) for antiquark ?? */
void 
mult_rho0_field( int fdir,  int r0[], su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;  
  
  /* apply the symmetric shift operator */
  sym_shift_field(fdir, src, dest);
  FORALLSITES(i,s){
    short *h = hyp_coord(s, r0);
    /* \eta_k already in the phases                                     * 
     * only the \epsilon for the anti-quark is needed times (-1)^(x+y+z)*
     * this is equal to (-1)^t                                          */
    if( h[TUP] == 1)
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
}

#endif

#if 0
/* "Multiply by" the quark-antiquark local a1 operator */
void 
mult_a1_field( int pdir, int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is gamma_pdir, (-1)^(x+y+z+t), another (-1)^(x+y+z+t)
     for antiquark */
  register int i;
  register site *s;
  if(this_node==0)printf("OOPS, mult_a1 NOT WRITTEN\n");
  terminate(0);
  FORALLSITES(i,s){
  }
}

/* "Multiply by" the quark-antiquark local b1 operator */
void 
mult_b1_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is gamma_pdir, gamma_0, (-1)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark */
  register int i;
  register site *s;
  if(this_node==0)printf("OOPS, mult_b1 NOT WRITTEN\n");
  terminate(0);
  FORALLSITES(i,s){
  }
}
#endif

/*------------------------------------------------------------------*/
/* Symbolic names for propagators.  prop_SOURCE_SINK */
/* Name by KS taste content e.g. Goldstone pion = pion5 */
/* operators:
   pion5:	local 0-+:  (taste)gamma_5     partner=0+-  phase=(1)
   pion05:	local 0-+:  gamma_0 gamma_5     partner=0++  phase=(-1)^(x+y+z+t)
   pioni5:	one-link 0-+:  gamma_i gamma_5	partner=0+-
   pionij:	one-link 0-+:  gamma_i gamma_j  partner=0++
   pioni:	two-link 0-+:  gamma_i          partner=0++
   pioni0:	two-link 0-+:  gamma_i gamma_0	partner=0+-
   pions:	three-link 0-+:  1 ("singlet")  partner=0++
   pion0:	three-link 0-+:  gamma_0	partner=0+-
   
   rhoi:	local 1--: gamma_i              partner=1+-  phase=(-1)^(dir) (VT)
   rhoi0:	local 1--: gamma_i gamma_0      partner=1++  phase=(-1)^(x+y+z+t+dir) (PV)
   rhois:	one-link 1--: 1 ("singlet")     partner=1+-
   rho0:	one-link 1--: gamma_0           partner=1++
*/

enum spin_taste_type {
  pion5,
  pion05,
  pioni5,
  pionij,
  pioni,
  pioni0,
  pions,
  pion0,
  rhoi,
  rhox,
  rhoy,
  rhoz,
  rhoi0,
  rhox0,
  rhoy0,
  rhoz0,
  rhoxs,
  rhoys,
  rhozs,
  rhots,
  rhois,
  rho0,
  rhoxsfn,
  rhoysfn,
  rhozsfn,
  rhotsfn,
  rhoisfn,
  MAX_SPIN_TASTE
};

static char *spin_taste_string[MAX_SPIN_TASTE]  = { 
  "pion5",
  "pion05",
  "pioni5",
  "pionij",
  "pioni",
  "pioni0",
  "pions",
  "pion0",
  "rhoi",
  "rhox",
  "rhoy",
  "rhoz",
  "rhoi0",
  "rhox0",
  "rhoy0",
  "rhoz0",
  "rhoxs",
  "rhoys",
  "rhozs",
  "rhots",
  "rhois",
  "rho0",
  "rhoxsfn",
  "rhoysfn",
  "rhozsfn",
  "rhotsfn",
  "rhoisfn",
};

/*------------------------------------------------------------------*/
/* Choices requiring a gauge field */

#ifdef NO_GAUGE_FIELD

/* All of the operators with covariant fn shifts obviously require the gauge field */
static int 
local_operator(int index){

  return
    index == pion5 ||
    index == pion05 ||
    index == rhox ||
    index == rhoy ||
    index == rhoz ||
    index == rhoi ||
    index == rhox0 ||
    index == rhoy0 ||
    index == rhoz0 ||
    index == rhoi0;
}
#endif

/*------------------------------------------------------------------*/
/* Map a label to the spin-taste index */

int 
spin_taste_index(char *label){
  int i;
  for(i = 0; i < MAX_SPIN_TASTE; i++){
    if(strcmp(label,spin_taste_string[i]) == 0)break;
  }
  if(i == MAX_SPIN_TASTE)
    return -1;  /* Error condition */

#ifdef NO_GAUGE_FIELD
  if(! local_operator(i)){
    printf("spin_taste_index: ERROR IN INPUT: field operation not supported for this application\n");
    return -1;
  }
#endif
  return i;
}

/* Map an index to the label */

char *
spin_taste_label(int index){
  return spin_taste_string[index];
}

/*------------------------------------------------------------------*/
/* Generic spin-taste operator                                      */

void 
spin_taste_op(int index, int r0[],
	      su3_vector *dest, su3_vector *src){
  switch(index){
  case pion5:
    mult_pion5_field(r0, src, dest);
    break;
  case pion05:
    mult_pion05_field(r0, src, dest);
    break;
#ifndef NO_GAUGE_FIELD
  case pioni5:
    mult_pioni5_field(ZUP, r0, src, dest);
    break;
  case pionij:
    mult_pionij_field(ZUP, r0, src, dest);
    break;
  case pioni:
    mult_pioni_field(ZUP, r0, src, dest);
    break;
  case pioni0:
    mult_pioni0_field(ZUP, r0, src, dest);
    break;
  case pions:
    mult_pions_field(r0, src, dest);
    break;
  case pion0:
    mult_pion0_field(r0, src, dest);
    break;
#endif
  case rhox:
    mult_rhoi_field(XUP, r0, src, dest);
    break;
  case rhoy:
    mult_rhoi_field(YUP, r0, src, dest);
    break;
  case rhoz:
  case rhoi:
    mult_rhoi_field(ZUP, r0, src, dest);
    break;
  case rhox0:
    mult_rhoi0_field(XUP, r0, src, dest);
    break;
  case rhoy0:
    mult_rhoi0_field(YUP, r0, src, dest);
    break;
  case rhoz0:
  case rhoi0:
    mult_rhoi0_field(ZUP, r0, src, dest);
    break;
#ifndef NO_GAUGE_FIELD
  case rhoxs:
    mult_rhois_field(XUP, r0, src, dest);
    break;
  case rhoys:
    mult_rhois_field(YUP, r0, src, dest);
    break;
  case rhozs:
  case rhois:
    mult_rhois_field(ZUP, r0, src, dest);
    break;
  case rhots:
    mult_rhois_field(TUP, r0, src, dest);
    break;
  case rho0:
    mult_rho0_field(ZUP, r0, src, dest);
    break;
#endif
  default:
    printf("spin_taste_op(%d): Bad spin-taste index %d\n",this_node, index);
    terminate(1);
  }
}

/********************************************************************/
/* Fat-Naik variants                                                */
/* 10/3/09 C DeTar                                                  */
/********************************************************************/

#ifndef NO_GAUGE_FIELD

/* Apply the symmetric shift operator in direction "dir" *
 * Fat and long links are used instead of unfattened links          */

static void 
sym_shift_fn_field(imp_ferm_links_t *fn, int dir, 
		   su3_vector *src, su3_vector *dest)
{
  char myname[] = "sym_shift_fn_field";
  if(fn == NULL){
    node0_printf("%s: Called with NULL FN links\n", myname);
    terminate(1);
  }
  
  clear_v_field(dest);

  /* Apply Fat-Naik shift operation to src in forward dir */
  //dslash_fn_dir(src, dest, EVENANDODD, fn, dir, +1, 1., 1.);
  dslash_fn_dir(src, dest, EVENANDODD, fn, dir, +1, 1., 0.);

  /* Add to it the Fat-Naik shift operation in backward dir */
  //dslash_fn_dir(src, dest, EVENANDODD, fn, dir, -1, 1., 1.);
}

/* "Multiply by" the quark-antiquark Fat-Naik rho operator */
/* This is the conserved vector current for the asqtad and HISQ actions -CD */

static void 
mult_rhois_fn_field( imp_ferm_links_t *fn, int fdir, int r0[],
		     su3_vector *src, su3_vector *dest )
{
  register int i;
  register site *s;  
  
  /* apply the symmetric shift FN operator (uses fat and long links) */
  sym_shift_fn_field(fn, fdir, src, dest);
  FORALLSITES(i,s){
    /* \eta_k already in the phases                   * 
     * only the \epsilon for the anti-quark is needed */
    if(s->parity==ODD)
      scalar_mult_su3_vector( dest+i, -1.0, dest+i );
  }
}

#endif

/*------------------------------------------------------------------*/
/* Generic spin-taste operator                                      */

void 
spin_taste_op_fn( imp_ferm_links_t *fn, int index, int r0[],
		  su3_vector *dest, su3_vector *src){
  switch(index){
#ifndef NO_GAUGE_FIELD
  case rhoxsfn:
    mult_rhois_fn_field(fn, XUP, r0, src, dest);
    break;
  case rhoysfn:
    mult_rhois_fn_field(fn, YUP, r0, src, dest);
    break;
  case rhozsfn:
  case rhoisfn:
    mult_rhois_fn_field(fn, ZUP, r0, src, dest);
    break;
  case rhotsfn:
    mult_rhois_fn_field(fn, TUP, r0, src, dest);
    break;
#endif
  default:
    /* For all non-FN operators */
    spin_taste_op(index, r0, dest, src);
  }
}




