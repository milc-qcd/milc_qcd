/******** flavor_ops.c *************/
/* MIMD version 7*/

/* THIS VERSION IS DEPRECATED.  USE flavor_ops2.c INSTEAD */ 
/* Implementation of the flavor (\Xi_\mu) operators. 
   See Golterman & Smit Nulc. Phys. B245 (1984) 61.
   They are used for constructing the non-local pion sources
   and in the measurement of all the componenets of \bar{\psi}\psi
*/

#include "generic_ks_includes.h"

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

/* Apply the symmetric shift opperator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              */
void sym_shift(int dir, field_offset src,field_offset dest)
{
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec;
  
  tvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

  tag[0] = start_gather_site( src, sizeof(su3_vector), dir, EVENANDODD ,gen_pt[0] );
  FORALLSITES(i,s)
    {
      mult_adj_su3_mat_vec( &(s->link[dir]), (su3_vector *)F_PT(s,src), 
			    &(tvec[i]) ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
				  EVENANDODD ,gen_pt[1] );
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
    mult_su3_mat_vec( &(s->link[dir]), (su3_vector *)gen_pt[0][i], 
		      (su3_vector *)F_PT(s,dest) ) ;    
    }
  wait_gather(tag[1]);
  FORALLSITES(i,s)
    {
      add_su3_vector( (su3_vector *)F_PT(s,dest), (su3_vector *)gen_pt[1][i], 
		      (su3_vector *)F_PT(s,dest) ) ;    
    }
  /* Now devide by 2 eq. (4.2b) of Golderman's Meson paper*/
 FORALLSITES(i,s)
   {
     scalar_mult_su3_vector( (su3_vector *)F_PT(s,dest), .5,
			     (su3_vector *)F_PT(s,dest) );
   }
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
  free(tvec);
}

/* It applies the symmetric shift with directions                       *
 * stored in the array d. Each shift is multiplied by \zeta_k           *
 * n is the number of shifts                                            *
 * This is the E_\mu(x,y)=\Xi_\mu operator defined by Golterman.        *
 * Nucl. Phys. B245  eq.3.5 and eq. 4.2b                                */
void zeta_shift(int n, int *d, field_offset src, field_offset dest )
{
  register int i,c ;
  short coords[4];
  register site *s;
  register Real sign ;

  if(dest==F_OFFSET(tempvec[0])||src==F_OFFSET(tempvec[0]))
    { 
      node0_printf("zeta_shift(): ERROR! dest or src should not be ");
      node0_printf("tempvec[0]\n" );
      terminate(101) ;
    }
  for(c=0;c<n;c++)
    {  
      /* Do the shift in d[c] */ 
      if(c==0)
	/* first time from source */
	sym_shift(d[c], src, F_OFFSET(tempvec[0]));
      else
	/* second time from dest */
	sym_shift(d[c], dest, F_OFFSET(tempvec[0]));
      /* Multiply by \zeta_d[c]. Because the phases are               *
       * on we multiply by \zeta * \eta = \epsilon * (-1)^coord[d[c]] */
      FORALLSITES(i,s){
	coords[XUP] = s->x; coords[YUP] = s->y;
	coords[ZUP] = s->z; coords[TUP] = s->t;
	/* The \epsilon */
	if (s->parity==EVEN) 
	  sign =  1.0 ;
	else
	  sign = -1.0 ;
	/* And the (-1)^coord[d[c]] */
	if(coords[d[c]]%2==1) sign=-sign ;
	scalar_mult_su3_vector(&(s->tempvec[0]), sign, 
			       (su3_vector *)F_PT(s,dest) );
      }
    }
}

/* It applies the symmetric shift with directions                       *
 * stored in the array d. Each shift is multiplied by \eta_k            *
 * In fact since \eta_k are already absorbed into the U matrices do     *
 * no action is needed. Just do the symmetric shift.                    *
 * n is the number of shifts                                            */
void eta_shift(int n, int *d, field_offset src, field_offset dest )
{
  int c ;
  field_offset temp0,temp1,tmp ;
  
  if(n==1)
    {
      sym_shift(d[0], src, dest);
    }
  else
    {
      temp0 = F_OFFSET(tempvec[0]) ;
      temp1 = F_OFFSET(tempvec[1]) ;
      if(dest==temp0||dest==temp1)
	{ 
	  node0_printf("eta_shift(): ERROR! dest should not be ");
	  node0_printf("tempvec[0] or tempvec[1]\n" );
	  terminate(102) ;
	}
      if(src==temp0||src==temp1)
	{ 
	  node0_printf("eta_shift(): ERROR! src should not be ");
	  node0_printf("tempvec[0] or tempvec[1]\n" );
	  terminate(101) ;
	}
      sym_shift(d[0], src, temp0 );
      for(c=1;c<n-1;c++)
	{  
	  sym_shift(d[c], temp0, temp1);
	  /* switch the pointers */
	  tmp = temp0   ;
	  temp0 = temp1 ;
	  temp1 = tmp   ;
	}
      /* do the last shift */
      sym_shift(d[n-1], temp0, dest );
    }
}

/* Multiply by the Xi_mu flavor operator */
void mult_flavor_vector(int mu, field_offset src, field_offset dest )
{
  int d[1] ;
  
  d[0] = mu ;
  zeta_shift(1,d,src,dest) ;
}

/* Multiply by the 1/2(Xi_mu Xi_nu - Xi_nu Xi_mu)flavor operator */
void mult_flavor_tensor(int mu, int nu, field_offset src, field_offset dest )
{ 
  register int i;
  register site *s;
  int d[2];

  if(dest==F_OFFSET(tempvec[1])||src==F_OFFSET(tempvec[1]))
    { 
      node0_printf("mult_flavor_pseudovector(): ERROR! dest or src should ");
      node0_printf("not be tempvec[1]\n" );
      terminate(101) ;
    }

  d[0] = mu ; d[1]=nu ;
  zeta_shift(2, d,src,dest) ;
 
  d[0] = nu ; d[1]=mu ;
  zeta_shift(2,d,src,F_OFFSET(tempvec[1])) ;

  
  FORALLSITES(i,s){
    scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest),
			       &(s->tempvec[1]), -1.0);
    scalar_mult_su3_vector((su3_vector *)F_PT(s,dest), 0.5, 
			   (su3_vector *)F_PT(s,dest) );
  }
}

/* Multiply by the Xi_mu Xi_5 flavor operator */
void mult_flavor_pseudovector(int mu, field_offset src, field_offset dest )
{
  register int i;
  register site *s;
  int p ; 

  if(dest==F_OFFSET(tempvec[1])||src==F_OFFSET(tempvec[1]))
    { 
      node0_printf("mult_flavor_pseudovector(): ERROR! dest or src should ");
      node0_printf("not be tempvec[1]\n" );
      terminate(101) ;
    }
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }   
   for(p=0;p<24;p++)
      if(eps[p].d[0]==mu)
	{
	  zeta_shift(3,&eps[p].d[1],src,F_OFFSET(tempvec[1])) ;
	  /* Multiply the extra 1/6 needed by the definition    *
	   * of the operator (number of permutations)           */
	  FORALLSITES(i,s){
	    scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
				       &(s->tempvec[1]), eps[p].sign/6.0);
	  }
	}
}

/* Multiply by the Xi_5 flavor operator */
void mult_flavor_pseudoscalar(field_offset src, field_offset dest )
{ 
  register int i;
  register site *s;
  int p ; 
  
  if(dest==F_OFFSET(tempvec[1])||src==F_OFFSET(tempvec[1]))
    { 
      node0_printf("mult_flavor_pseudoscalar(): ERROR! dest or src should ");
      node0_printf("not be tempvec[1]\n" );
      terminate(101) ;
    }
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }   
  for(p=0;p<24;p++)
     {
       zeta_shift(4,eps[p].d,src,F_OFFSET(tempvec[1])) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
				    &(s->tempvec[1]), eps[p].sign/24.0);
       }
     } 
}

/* Multiply by the Gamma_mu spin operator */
void mult_spin_vector(int mu, field_offset src, field_offset dest )
{
  int d[1] ;
  
  d[0] = mu ;
  eta_shift(1,d,src,dest) ;
}

/* Multiply by the 1/2(Gamma_mu Gamma_nu - Gamma_nu gamma_mu) spin operator */
void mult_spin_tensor(int mu, int nu, field_offset src, field_offset dest )
{ 
  register int i;
  register site *s;
  int d[2];

  if(dest==F_OFFSET(tempvec[2])||src==F_OFFSET(tempvec[2]))
    { 
      node0_printf("mult_flavor_pseudovector(): ERROR! dest or src should ");
      node0_printf("not be tempvec[2]\n" );
      terminate(101) ;
    }

  d[0] = mu ; d[1]=nu ;
  eta_shift(2, d,src,dest) ;
 
  d[0] = nu ; d[1]=mu ;
  eta_shift(2,d,src,F_OFFSET(tempvec[2])) ;

  FORALLSITES(i,s){
    scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest),
			       &(s->tempvec[2]), -1.0);
    scalar_mult_su3_vector((su3_vector *)F_PT(s,dest), 0.5, 
			   (su3_vector *)F_PT(s,dest) );
  }
}

/* Multiply by the Gamma_mu Gamma_5 spin operator */
void mult_spin_pseudovector(int mu, field_offset src, field_offset dest )
{
  register int i;
  register site *s;
  int p ; 

  if(dest==F_OFFSET(tempvec[2])||src==F_OFFSET(tempvec[2]))
    { 
      node0_printf("mult_spin_pseudovector(): ERROR! dest or src should ");
      node0_printf("not be tempvec[2]\n" );
      terminate(101) ;
    }
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }   
   for(p=0;p<24;p++)
      if(eps[p].d[0]==mu)
	{
	  eta_shift(3,&eps[p].d[1],src,F_OFFSET(tempvec[2])) ;
	  /* Multiply the extra 1/6 needed by the definition    *
	     of the operator (number of permutations)           */
	  FORALLSITES(i,s){
	    scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
				       &(s->tempvec[2]), eps[p].sign/6.0);
	  }
	}
}

/* Multiply by the Gamma_5 spin operator */
void mult_spin_pseudoscalar(field_offset src, field_offset dest )
{ 
  register int i;
  register site *s;
  int p ; 
  
  if(dest==F_OFFSET(tempvec[2])||src==F_OFFSET(tempvec[2]))
    { 
      node0_printf("mult_spin_pseudoscalar(): ERROR! dest or src should ");
      node0_printf("not be tempvec[2]\n" );
      terminate(101) ;
    }
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec((su3_vector *)F_PT(s,dest));
  }   
  for(p=0;p<24;p++)
     {
       eta_shift(4,eps[p].d,src,F_OFFSET(tempvec[2])) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector((su3_vector *)F_PT(s,dest), 
				    &(s->tempvec[2]), eps[p].sign/24.0);
       }
     } 
}
