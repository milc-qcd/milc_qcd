/************************* check_invert2.c *******************************/
/* MIMD version 7 */
#include "ks_imp_includes.h"	/* definitions files and prototypes */


/* FOR TESTING: multiply src by Mdagger M and check against dest */
void check_invert2( su3_vector *src, su3_vector *dest, 
		    Real mass, Real tol, int parity,
		    imp_ferm_links_t *fn){
  register int i,k,flag;
  Real r_diff, i_diff;
  double sum,sum2,dflag,dmaxerr,derr;
  int otherparity = 0;
  su3_vector *temp = create_v_field();
  su3_vector *cg_p = create_v_field();

  switch(parity){
  case EVEN: otherparity=ODD; break;
  case ODD: otherparity=EVEN; break;
  case EVENANDODD: otherparity=EVENANDODD; break;
  }
					     
  dslash_field( src, temp, otherparity, fn);
  dslash_field( temp, cg_p, parity, fn);

  FORSOMEFIELDPARITY(i,parity){
    scalar_mult_su3_vector( cg_p+i, -1.0, cg_p+i);
    scalar_mult_add_su3_vector( cg_p+i, src+i, +4.0*mass*mass, cg_p+i );
  }
  sum2=sum=0.0;
  dmaxerr=0;
  flag = 0;
  FORSOMEFIELDPARITY(i,parity){
    for(flag=0,k=0;k<3;k++){
      r_diff = dest[i].c[k].real - cg_p[i].c[k].real;
      i_diff = dest[i].c[k].imag - cg_p[i].c[k].imag;
      if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
	       i,k, dest[i].c[k].real,dest[i].c[k].imag,
	       cg_p[i].c[k].real, cg_p[i].c[k].imag);
	flag++;
      }
      derr = r_diff*r_diff + i_diff*i_diff;
      if(derr>dmaxerr)dmaxerr=derr;
      sum += derr;
    }
    sum2 += magsq_su3vec( dest+i );
  }
  g_doublesum( &sum );
  g_doublesum( &sum2 );
  dflag=flag;
  g_doublesum( &dflag );
  g_doublemax( &dmaxerr );
  if(this_node==0){
    printf("Inversion checked, frac. error = %e\n",sqrt(sum/sum2));
    printf("Flagged comparisons = %d\n",(int)dflag);
    printf("Max err. = %e frac. = %e\n",sqrt(dmaxerr),
	   sqrt(dmaxerr*volume/sum2));
    fflush(stdout);
  }
  destroy_v_field(temp);
  destroy_v_field(cg_p);
}
