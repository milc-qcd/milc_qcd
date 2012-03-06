/******************************** flavor_ops2.c *************************/
/* MIMD version 7                                                       */

/* Formerly flavor_ops2.c */

/* Implementation of the spin-taste ("flavor") (\Xi_\mu) operators.
   See Golterman & Smit Nulc. Phys. B245 (1984) 61.
   They are used for constructing the non-local pion sources
   and in the measurement of all the componenets of \bar{\psi}\psi

   CD 07/27/09 This version works with site-major fields
   Ths version also includes a set of single time-slice meson operators.

   CD 12/30/11 Generalized to any spin-taste combination
*/

/* Entry points

   general_spin_taste_op  Any spin/taste and any gauge links but not FN
   spin_taste_op          Limited spin/taste set and APE links
   spin_taste_op_fn       Any spin/taste and APE links except FN links with FN operators

   spin_taste_index       Convert character label to index
   spin_taste_label       Convert index to label

*/


#include "generic_ks_includes.h"
#include <string.h>

#include "../include/gammatypes.h"

/*------------------------------------------------------------------*/
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

/* Compute the parity of the site s relative to offset r0. */
/* Return 0x0 for even and 0x1 for odd */
static short 
hyp_parity_bit(site *s, int r0[]){
  short p;
  if((r0[XUP] + r0[YUP] + r0[ZUP] + r0[TUP]) % 2 == 0)
    p = EVEN;
  else p = ODD;
  if(p == s->parity)
    return 0x0;
  else
    return 0x1;
}

#ifndef NO_GAUGE_FIELD

/*------------------------------------------------------------------*/
/* Apply the symmetric shift operator in direction "dir" *
 * This is the explicit version                           *
 * Covar
 * The KS phases MUST BE in the links                     */
static void 
sym_shift_field(int dir, su3_vector *dest, su3_vector *src, su3_matrix *links)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec = create_v_field();

  tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD, gen_pt[0] );
  /* With ONE_SIDED_SHIFT defined, the shift is asymmetric */
#ifndef ONE_SIDED_SHIFT
  FORALLFIELDSITES(i)
    {
      mult_adj_su3_mat_vec( links+4*i+dir, src+i, tvec+i ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
			      EVENANDODD ,gen_pt[1] );
#endif
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
      mult_su3_mat_vec( links+4*i+dir, (su3_vector *)gen_pt[0][i], dest+i );
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

/*------------------------------------------------------------------*/
/* Apply the symmetric shift with directions                            *
 * stored in the array d. Each shift is multiplied by \zeta_k           *
 * n is the number of shifts                                            *
 * This is the E_\mu(x,y)=\Xi_\mu operator defined by Golterman.        *
 * Nucl. Phys. B245  eq.3.5 and eq. 4.2b                                */

static void 
zeta_shift_field(int n, int *d, int r0[], su3_vector *dest, 
		 su3_vector *src, su3_matrix *links )
{
  int i,c ;
  site *s;
  short *h;
  su3_vector *tvec = create_v_field();
  
  for(c=0;c<n;c++)
    {  
      /* Do the shift in d[c] */ 
      if(c==0)
	/* first time from source */
	sym_shift_field(d[c], tvec, src, links);
      else
	/* other times from dest */
	sym_shift_field(d[c], tvec, dest, links);
      /* Multiply by \zeta_d[c]. Because the phases are               *
       * on we multiply by \zeta * \eta = \epsilon * (-1)^coord[d[c]] */
      FORALLSITES(i,s){
	h = hyp_coord(s, r0);
	/* epsilon times (-1)^coord[d[c]] for others */
	if( hyp_parity_bit(s, r0) ^ h[d[c]] )
	  scalar_mult_su3_vector(tvec+i, -1., dest+i );
	else
	  dest[i] = tvec[i];
      }
    }
  destroy_v_field(tvec);
}

#endif

/*------------------------------------------------------------------*/
static Real
spin_sign(int spin, int r0[], site *s){
  /* Compute (-)^[spin dot (x-r0)] epsilon(x-r0)^spin */
  /* Same as prod_\mu [eta_\mu(x-r0) zeta_\mu(x-r0)]^(s_\mu) */
  int j, mask;
  Real sign = 1.;
  short *h = hyp_coord(s, r0);
  short hp = hyp_parity_bit(s, r0);
  
  /* For each nonzero gamma_mu bit in "spin",
     a factor of (-)^(x[mu]-r0[mu]) epsilon(x) */
  mask = 1;
  FORALLUPDIR(j){
    if( (spin & mask) && (hp ^ h[j]) ) sign = -sign;
    mask <<= 1;
  }

  return sign;
}

#ifndef NO_GAUGE_FIELD

static void
spin_sign_field(int spin, int r0[], su3_vector *dest, su3_vector *src){
  /* Apply the spin_sign operation to an entire field */
  int i;
  site *s;

  FORALLSITES(i,s){
    if(spin_sign(spin, r0, s) > 0)
      dest[i] = src[i];
    else 
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
  }
}

#endif

/*------------------------------------------------------------------*/
static short
antiquark_sign_flip(int r0[], site *s){
  return hyp_parity_bit(s, r0) == 0x1;
}

#ifndef NO_GAUGE_FIELD

/* Extra (-)^(x+y+z+t) to allow for antiquark */
static void
antiquark_sign_flip_field(int r0[], su3_vector *dest, su3_vector *src){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    if( antiquark_sign_flip(r0, s) )
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
    else 
      dest[i] = src[i];
  }
}

/*------------------------------------------------------------------*/
static void
sign_flip_field(su3_vector *dest, su3_vector *src){
  int i;
  site *s;
  FORALLSITES(i,s){
    scalar_mult_su3_vector(src+i, -1.0, dest+i );
  }
}

#endif

/*------------------------------------------------------------------*/
static void 
local(int spin, int r0[], su3_vector *dest, su3_vector *src){

  int i;
  site *s;
  Real sign;

  FORALLSITES(i,s){
    sign = spin_sign(spin, r0, s);
    if(antiquark_sign_flip(r0, s))sign = -sign;
    if( sign > 0 )
      dest[i] = src[i];
    else 
      scalar_mult_su3_vector( src+i, -1.0, dest+i );
  }
}

#ifndef NO_GAUGE_FIELD

/*------------------------------------------------------------------*/
static void 
one_link(int spin, int dir, int r0[], su3_vector *dest, 
	 su3_vector *src, su3_matrix *links){

  int n = 1;
  int c[1] = {dir};
  su3_vector *tvec = create_v_field();

  spin_sign_field(spin, r0, tvec, src);
  zeta_shift_field(n, c, r0, dest, tvec, links);
  antiquark_sign_flip_field(r0, dest, dest);

  destroy_v_field(tvec);
}

/*------------------------------------------------------------------*/
static void 
two_link(int spin, int dir1, int dir2, int r0[], su3_vector *dest, 
	 su3_vector *src, su3_matrix *links){
  int i;
  site *s;
  int n = 2;
  int c[2];
  su3_vector *tvec0 = create_v_field();
  su3_vector *tvec1 = create_v_field();

  spin_sign_field(spin, r0, dest, src);

  c[0] = dir2; c[1] = dir1;
  zeta_shift_field(n, c, r0, tvec0, dest, links);

  c[0] = dir1; c[1] = dir2;
  zeta_shift_field(n, c, r0, tvec1, dest, links);

  FORALLSITES(i,s){
    sub_su3_vector( tvec0+i, tvec1+i, dest+i );
    if( antiquark_sign_flip(r0,s) )
      scalar_mult_su3_vector( dest+i, -0.5, dest+i );
    else
      scalar_mult_su3_vector( dest+i,  0.5, dest+i );
  }

  destroy_v_field(tvec1);
  destroy_v_field(tvec0);
}

/*------------------------------------------------------------------*/
static void 
three_link(int spin, int r0[], su3_vector *dest, su3_vector *src, su3_matrix *links){
  /* For each gamma_mu in the "spin" gamma matrix, the 
     operator gets a factor (-)^x[mu] times a factor (-)^(x+y+z+t).
     Then an overall factor (-)^(x+y+z+t) for the antiquark. */

  int i, c;
  site *s;
  su3_vector *tvec0 = create_v_field();
  su3_vector *tvec1 = create_v_field();
  
  /* All permutation with appropriate sign */
  struct {
    int d[3];
    Real sign ;
  } p[6]={{{XUP,YUP,ZUP},+1.0/6.0},
	  {{YUP,ZUP,XUP},+1.0/6.0},
	  {{ZUP,XUP,YUP},+1.0/6.0},
	  {{XUP,ZUP,YUP},-1.0/6.0},
	  {{YUP,XUP,ZUP},-1.0/6.0},
	  {{ZUP,YUP,XUP},-1.0/6.0}}; /* The factor of 6 accounts for the *
				* multiplicity of the permutations */
  
  clear_v_field(dest);
  spin_sign_field(spin, r0, tvec0, src);
  for(c=0;c<6;c++)
    {
      zeta_shift_field(3,p[c].d,r0,tvec1,tvec0,links);
      FORALLSITES(i,s){
	scalar_mult_sum_su3_vector(dest+i, tvec1+i, p[c].sign );
      }
    }
  /* multiply by \epsilon for the anti-quark */
  antiquark_sign_flip_field(r0, dest, dest);
  
  destroy_v_field(tvec1);
  destroy_v_field(tvec0);
}

#endif

/*------------------------------------------------------------------*/
/* Apply a general spin-taste operator to a field */
void
general_spin_taste_op(enum gammatype spin_index, enum gammatype taste_index, int r0[],
		      su3_vector *dest, su3_vector *src, su3_matrix *links){

  /* Convert gamma label to hexadecimal */
  short spin = gamma_hex(spin_index);
  short taste = gamma_hex(taste_index);

  /* Compute offset */

  int offset = spin ^ taste;

  switch(offset)
    {
    case 0: local(spin, r0, dest, src); break;
#ifndef NO_GAUGE_FIELD
    case 1: one_link(spin, XUP, r0, dest, src, links); break;
    case 2: one_link(spin, YUP, r0, dest, src, links); break;
    case 3: two_link(spin, XUP, YUP, r0, dest, src, links); break;
    case 4: one_link(spin, ZUP, r0, dest, src, links); break;
    case 5: two_link(spin, ZUP, XUP, r0, dest, src, links); break;
    case 6: two_link(spin, YUP, ZUP, r0, dest, src, links); break;
    case 7: three_link(spin, r0, dest, src, links); break;
#endif
    default: printf("Time split operator not supported\n");
    }
}

/*------------------------------------------------------------------*/
/* Procedures for backward compatibility with former flavor_ops2.c  */
/*------------------------------------------------------------------*/

/* operators:
   pion5:	local 0-+:  (flavor)gamma_5     partner=0+-  phase=(1)
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

/* Convert dir to corresponding gamma index */

static enum gammatype 
dir2gi(int dir, char *myname){

  enum gammatype gi = 0;

  switch(dir){
  case XUP: gi = GX; break;
  case YUP: gi = GY; break;
  case ZUP: gi = GZ; break;
  default: 
    node0_printf("%s: bad dir = %d\n",myname,dir);
    terminate(1);
  }
  return gi;
}

static enum gammatype 
dir2gi0(int dir, char *myname){

  enum gammatype gi0 = 0;

  switch(dir){
  case XUP: gi0 = GXT; break;
  case YUP: gi0 = GYT; break;
  case ZUP: gi0 = GZT; break;
  default: 
    node0_printf("%s: bad dir = %d\n",myname,dir);
    terminate(1);
  }
  return gi0;
}

#ifndef NO_GAUGE_FIELD

static enum gammatype 
dir2g5i(int dir, char *myname){

  enum gammatype g5i = 0;
  
  switch(dir){
  case XUP: g5i = G5X; break;
  case YUP: g5i = G5Y; break;
  case ZUP: g5i = G5Z; break;
  default: 
    node0_printf("%s: bad dir = %d\n",myname,dir);
    terminate(1);
  }
  
  return g5i;
}

static enum gammatype 
dir2gij(int dir, char *myname){

  enum gammatype gij = 0;
  
  switch(dir){
  case XUP: gij = GYZ; break;
  case YUP: gij = GZX; break;
  case ZUP: gij = GXY; break;
  default: 
    node0_printf("%s: bad dir = %d\n",myname,dir);
    terminate(1);
  }

  return gij;
}

#endif

/* "Multiply by" the quark-antiquark local pion operator */
/* Here the operator is gamma_5 x gamma_5 times (-1)^(x+y+z+t) */
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

  general_spin_taste_op( G1, G1, r0, dest, src, NULL);
}

#ifndef NO_GAUGE_FIELD

/* "Multiply by" the one link pion operator.
   "pi_i5" or gamma_5 x gamma_i gamma_5 */
/* Need a minus sign for backward compatibility */
void 
mult_pioni5_field( int fdir, int r0[], su3_vector *src, 
		   su3_vector *dest, su3_matrix *links ){

  /* operator is (-)^fdir, another (-1)^(x+y+z+t) for antiquark, 
     so multiply by (-1)^(x+y+z+t) (-)^fdir */
  char myname[] = "mult_pioni5_field";
  enum gammatype g5i = dir2g5i(fdir, myname);

  general_spin_taste_op( G5, g5i, r0, dest, src, links);
  sign_flip_field(dest, dest);
}

/* "Multiply by" the one link pion operator.
   fdir is 0 for ij = 12, 1 for ij = 20, 2 for ij = 01.
   "pi_ij" or gamma_0 gamma_5 x gamma_i gamma_j 
   CD: Actually, by convention, this one is the parity-switched operator
   1 x gamma_k
*/
void 
mult_pionij_field( int fdir, int r0[], su3_vector *src, 
		   su3_vector *dest, su3_matrix *links ){
  /* operator is (-)^fdir (-)^(x+y+z+t), another (-1)^(x+y+z+t) for antiquark, 
     so multiply by (-1)^fdir */
  char myname[] = "mult_pionij_field";
  enum gammatype gi = dir2gi(fdir, myname);
  
  general_spin_taste_op( G1, gi, r0, dest, src, links);
}

/* "Multiply by" the two link pion operator. *
 * "pi_i or gamma_0 gamma_5 x gamma_i"                                    
   CD: Actually, by convention, this one is the parity-switched operator
   1 x gamma_j gamma_k
*/

void 
mult_pioni_field( int fdir, int r0[], su3_vector *src, 
		  su3_vector *dest, su3_matrix *links )
{

  char myname[] = "mult_pioni_field";
  enum gammatype gij = dir2gij(fdir, myname);

  general_spin_taste_op( G1, gij, r0, dest, src, links);

}


/* "Multiply by" the two link pion operator. */
/* pi_i0 or gamma_5 x gamma_0 gamma_i 
   CD: Actually, by convention, this one is the parity-switched operator
   gamma_0 x gamma_5 gamma_i
*/
void 
mult_pioni0_field( int fdir, int r0[], su3_vector *src, 
		   su3_vector *dest, su3_matrix *links )
{
  char myname[] = "mult_pioni0_field";
  enum gammatype g5i = dir2g5i(fdir, myname);
  
  general_spin_taste_op( GT, g5i, r0, dest, src, links);
}


/* "Multiply by" the three link pion operator.  */
/* pion_s or gamma_0 gamma_5 x 1 
   CD: Actually, by convention, this one is the parity-switched operator
   1 x gamma_0 gamma_5 */
void 
mult_pions_field(int r0[], su3_vector *src, su3_vector *dest,
		 su3_matrix *links)
{
  general_spin_taste_op( G1, G5T, r0, dest, src, links);
}

/* "Multiply by" the three link pion operator.  */
/* gamma_5 x gamma_0 
   CD: Actually, by convention, this one is the parity-switched operator
   gamma_0 x gamma_5 and then it needs an overall sign flip
*/
void 
mult_pion0_field(int r0[], su3_vector *src, su3_vector *dest, 
		 su3_matrix *links )
{
  general_spin_taste_op( GT, G5, r0, dest, src, links);
  sign_flip_field(dest, dest);
}
#endif

/* "Multiply by" the quark-antiquark local rho operator */
/* gamma_i x gamma_i */
void 
mult_rhoi_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest ){
  /* operator is gamma_pdir, another (-1)^(x+y+z+t) for antiquark */

  char myname[] = "mult_rhoi_field";
  enum gammatype gi = dir2gi(pdir, myname);

  general_spin_taste_op( gi, gi, r0, dest, src, NULL);

}

/* "Multiply by" the second quark-antiquark local rho operator */
void 
mult_rhoi0_field( int pdir,  int r0[], su3_vector *src, su3_vector *dest ){
  /* gamma_0 gamma_i x gamma_0 gamma_i */
  /* operator is gamma_0 gamma_pdir, another (-1)^(x+y+z+t) for antiquark */

  char myname[] = "mult_rhoi0_field";
  enum gammatype gi0 = dir2gi0(pdir, myname);

  general_spin_taste_op( gi0, gi0, r0, dest, src, NULL);
}

#ifndef NO_GAUGE_FIELD

/* "Multiply by" the quark-antiquark one link rho operator */
/* It is gamma_i x 1 times (-1)^(x+y+z+t) for antiquark */
void 
mult_rhois_field( int fdir,  int r0[], su3_vector *src, 
		  su3_vector *dest, su3_matrix *links )
{
  char myname[] = "mult_rhois_field";
  enum gammatype gi = dir2gi(fdir, myname);

  general_spin_taste_op( gi, G1, r0, dest, src, links);
}

/* "Multiply by" the second quark-antiquark one link rho operator */
/* gamma_i gamma_0 x gamma_0  times (-1)^(x+y+z+t) for antiquark ?? */
/* Need a minus sign for backward compatibility */
void 
mult_rho0_field( int fdir,  int r0[], su3_vector *src, 
		 su3_vector *dest, su3_matrix *links )
{
  char myname[] = "mult_rhoi0_field";
  enum gammatype gi0 = dir2gi0(fdir, myname);

  general_spin_taste_op( gi0, GT, r0, dest, src, links);
  sign_flip_field(dest, dest);
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
/* Map a bilinear spin-taste label to the spin-taste index */

/* Two label styles are in use:

   The more general label specifies Gamma(spin)-Gamma(taste) where
   Gamma is one of the 16 gamma labels in generic_wilson/gammas.c.
   For example "G5T-GXY" would be a single-time-slice pseudoscalar
   with taste gamma_x * gamma_y.

   For backward compatibility we also support the old labeling so that
   the same bilinear can be specified with "pionij".
   
*/

/*------------------------------------------------------------------*/
/* An encoding scheme to index both operator styles */

/* The compatibility operators are indexed below 128 */
/* The gamma-gamma operators are indexed from 128 up */

static int
encode_gamma_gamma_index(int s, int t){
  return s*16+t+128;
}

static int
decode_gamma_spin_index(int index){
  return (index-128)/16;
}

static int
decode_gamma_taste_index(int index){
  return (index-128)%16;
}

/* True if the index is gamma-gamma type */
static int
is_gamma_gamma_index(int index){
  return (index >= 128);
}

/*------------------------------------------------------------------*/
static char *
gamma_gamma_string(int index){

  static char label[32];
  int gamma_spin_index, gamma_taste_index;

  if( ! is_gamma_gamma_index(index) )
    return "\0";

  gamma_spin_index = decode_gamma_spin_index(index);
  gamma_taste_index = decode_gamma_taste_index(index);

  strcpy(label, gamma_label(gamma_spin_index));
  strcat(label, "-");
  strcpy(label, gamma_label(gamma_taste_index));

  return label;
}

/*------------------------------------------------------------------*/
/* Look up the operator label in the table */

static int
compatibility_style(char *label){
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

/*------------------------------------------------------------------*/
/* Parse the label using the gamma matrix table                     */

static int
gamma_gamma_style(char *label){

  char *hyphen;
  int gamma_spin_index, gamma_taste_index;

  hyphen = strstr(label, "-");
 
  if(hyphen == NULL){
    printf("Can't parse the spin_taste label\n");
    return -1;
  }

  /* Temporarily set hyphen to null */
  *hyphen = '\0';

  /* Look up gammas */
  gamma_spin_index = gamma_index(label);
  gamma_taste_index = gamma_index(hyphen+1);

  /* Restore hyphen */
  *hyphen = '-';

  if(gamma_spin_index == -1 || gamma_taste_index == -1){
    printf("Can't parse the spin_taste label\n");
    return -1;
  }

#ifdef NO_GAUGE_FIELD
  /* If the operator requires a shift, we must have a gauge field */
  if( ( gamma_hex(gamma_spin_index) ^ gamma_hex(gamma_taste_index) ) != 0){
    printf("spin_taste_index: ERROR IN INPUT: field operation not supported for this application\n");
    return -1;
  }
#endif

  /* Encode the gamma-gamma style operator index */

  return encode_gamma_gamma_index(gamma_spin_index, gamma_taste_index);

}

/*------------------------------------------------------------------*/
int 
spin_taste_index(char *label){

  if(label[0] != 'G')
    return compatibility_style(label);
  else
    return gamma_gamma_style(label);
}

/*------------------------------------------------------------------*/
/* Map an index to the corresponding label */

char *
spin_taste_label(int index){
  if(is_gamma_gamma_index(index))
    return gamma_gamma_string(index);
  else
    return spin_taste_string[index];
}

/*------------------------------------------------------------------*/
/* Compatiblity spin-taste operator                                 */

static void 
spin_taste_op_links(int index, int r0[], su3_vector *dest, 
		    su3_vector *src, su3_matrix *links){
  switch(index){
  case pion5:
  copy_v_field(dest, src);
    mult_pion5_field(r0, src, dest);
    break;
  case pion05:
    mult_pion05_field(r0, src, dest);
    break;
#ifndef NO_GAUGE_FIELD
  case pioni5:
    mult_pioni5_field(ZUP, r0, src, dest, links);
    break;
  case pionij:
    mult_pionij_field(ZUP, r0, src, dest, links);
    break;
  case pioni:
    mult_pioni_field(ZUP, r0, src, dest, links);
    break;
  case pioni0:
    mult_pioni0_field(ZUP, r0, src, dest, links);
    break;
  case pions:
    mult_pions_field(r0, src, dest, links);
    break;
  case pion0:
    mult_pion0_field(r0, src, dest, links);
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
    mult_rhois_field(XUP, r0, src, dest, links);
    break;
  case rhoys:
    mult_rhois_field(YUP, r0, src, dest, links);
    break;
  case rhozs:
  case rhois:
    mult_rhois_field(ZUP, r0, src, dest, links);
    break;
  case rhots:
    mult_rhois_field(TUP, r0, src, dest, links);
    break;
  case rho0:
    mult_rho0_field(ZUP, r0, src, dest, links);
    break;
#endif
  default:
    printf("spin_taste_op(%d): Bad spin-taste index %d\n",this_node, index);
    terminate(1);
  }
}


/*------------------------------------------------------------------*/
/* Spin-taste operator without the FN option                        */
void
spin_taste_op(int index, int r0[], su3_vector *dest, su3_vector *src){
#ifdef NO_GAUGE_FIELD
    spin_taste_op_links(index, r0, dest, src, NULL);
#else
    /* Use APE links for shifts */
    spin_taste_op_links(index, r0, dest, src, ape_links);
#endif
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
  //char myname[] = "sym_shift_fn_field";
  if(fn == NULL){
    node0_printf("sym_shift_fn_field: Called with NULL FN links\n");
    terminate(1);
  }
  
  clear_v_field(dest);

  /* Apply Fat-Naik shift operation to src in forward dir */
  //dslash_fn_dir(src, dest, EVENANDODD, fn, dir, +1, 1., 1.);
  dslash_fn_dir(src, dest, EVENANDODD, fn, dir, +1, 1., 0.);

  /* Add to it the Fat-Naik shift operation in backward dir */
  //dslash_fn_dir(src, dest, EVENANDODD, fn, dir, -1, -1., -1.);
}

/* "Multiply by" the quark-antiquark Fat-Naik rho operator */
/* This is part of the conserved vector current for the asqtad and
   HISQ actions -CD */

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
/* gamma_gamma spin-taste operator                                  */
static void
gamma_gamma_spin_taste_op(int index, int r0[], 
			  su3_vector *dest, su3_vector *src){

  int spin_index = decode_gamma_spin_index(index);
  int taste_index = decode_gamma_taste_index(index);

#ifdef NO_GAUGE_FIELD
  general_spin_taste_op(spin_index, taste_index, r0, dest, src, NULL);
#else
  /* Use APE links for shifts */
  general_spin_taste_op(spin_index, taste_index, r0, dest, src,
			ape_links);
#endif
}

/*------------------------------------------------------------------*/
/* Spin-taste operator including the FN option                      */

void 
spin_taste_op_fn( imp_ferm_links_t *fn, int index, int r0[],
		  su3_vector *dest, su3_vector *src){

  if(is_gamma_gamma_index(index))
    gamma_gamma_spin_taste_op(index, r0, dest, src);

  else {
    
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
}

