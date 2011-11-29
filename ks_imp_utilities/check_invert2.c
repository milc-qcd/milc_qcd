/************************* check_invert2.c *******************************/
/* MIMD version 7 */
#include "ks_imp_includes.h"	/* definitions files and prototypes */


/* FOR TESTING: multiply src by Mdagger M and check against dest */
void check_invert2( field_offset src, field_offset dest, 
		    field_offset temp, Real mass,
		    Real tol, int parity,
		    imp_ferm_links_t *fn){
  register int i,k,flag;
  register site *s;
  Real r_diff, i_diff;
  double sum,sum2,dflag,dmaxerr,derr;
  int otherparity = 0;

  switch(parity){
  case EVEN: otherparity=ODD; break;
  case ODD: otherparity=EVEN; break;
  case EVENANDODD: otherparity=EVENANDODD; break;
  }
					     
  dslash_site( src, temp, otherparity, fn);
  dslash_site( temp, F_OFFSET(cg_p), parity, fn);

  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector( &(s->cg_p), -1.0, &(s->cg_p));
    scalar_mult_add_su3_vector( &(s->cg_p), (su3_vector *)F_PT(s,src), 
				+4.0*mass*mass, &(s->cg_p) );
  }
  sum2=sum=0.0;
  dmaxerr=0;
  flag = 0;
  FORSOMEPARITY(i,s,parity){
    for(flag=0,k=0;k<3;k++){
      r_diff = ((su3_vector *)F_PT(s,dest))->c[k].real
	- s->cg_p.c[k].real;
      i_diff = ((su3_vector *)F_PT(s,dest))->c[k].imag
	- s->cg_p.c[k].imag;
      if( fabs(r_diff) > tol || fabs(i_diff) > tol ){
	printf("site %d color %d  expected ( %.4e , %.4e ) got ( %.4e , %.4e )\n",
	       i,k,
	       ((su3_vector *)F_PT(s,dest))->c[k].real,
	       ((su3_vector *)F_PT(s,dest))->c[k].imag,
	       s->cg_p.c[k].real,s->cg_p.c[k].imag);
	flag++;
      }
      derr = r_diff*r_diff + i_diff*i_diff;
      if(derr>dmaxerr)dmaxerr=derr;
      sum += derr;
    }
    sum2 += magsq_su3vec( (su3_vector *)F_PT(s,dest) );
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
}
