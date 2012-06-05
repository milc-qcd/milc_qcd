/***************** u1plaq.c *************************************/

/* Calculate u1 plaquette                                       */

/* MIMD version 7 */
/* ************************************************************	*/
/*								*/
/*								*/
/* Code for average space-space and space-time plaquettes 	*/
/*								*/
/*				b				*/
/*      		     o ->-- o				*/
/*       		     :	    :				*/
/*        ^	   	   c ^	    v a				*/
/*     	  |		     :	    :				*/
/*     	 dir1 		     o ==<= o				*/
/*             		       	d    				*/
/* 								*/
/*            dir2->						*/
/* 								*/
/* Author: S. Basak					      	*/
/* Last updated on 07.24.07					*/
/* CD modified 5/24/12						*/
/* 								*/
/* ************************************************************	*/

#include "generic_u1_includes.h"

void u1plaq(Real *splaq,Real *tplaq)
{

  int i,dir1,dir2;
  double ssplaq,stplaq;
  double s2plaq,t2plaq;
  msg_tag *mtag0,*mtag1;
  Real pre;
  Real *Atmp = create_r_field();
  Real *Astpl = create_r_field();

  ssplaq=stplaq=0.0;
  s2plaq=t2plaq=0.0;
  
  /* plaq. in dir1-dir2 plane */
  for(dir1=YUP;dir1<=TUP;dir1++){
    for(dir2=XUP;dir2<dir1;dir2++){
      
      mtag0 = declare_strided_gather(u1_A+dir2, 4*sizeof(Real),
				     sizeof(Real), dir1,
				     EVENANDODD,gen_pt[0]);
      prepare_gather(mtag0);
      do_gather(mtag0);
      /* = b */
      mtag1 = declare_strided_gather(u1_A+dir1, 4*sizeof(Real),
				     sizeof(Real), dir2,
				     EVENANDODD,gen_pt[1]);
      prepare_gather(mtag1);
      do_gather(mtag1);
      /* = a */
      FORALLFIELDSITES(i){
	Atmp[i] = u1_A[4*i+dir1] - u1_A[4*i+dir2];;	/* = (d^)c */
      }
      
      wait_gather(mtag0);
      FORALLFIELDSITES(i){
	Astpl[i] = Atmp[i] + *((Real *)gen_pt[0][i]);
      }					/* = (d^)cb */
      
      wait_gather(mtag1);
      FORALLFIELDSITES(i){
	Atmp[i] = Astpl[i] - *(Real *)gen_pt[1][i];  /* = (d^)cb(a^) */
      }					
      
      FORALLFIELDSITES(i){
	pre = cos(Atmp[i]);
	if(dir1==TUP) stplaq += pre;
	else	      ssplaq += pre;
	if(dir1==TUP) t2plaq += pre*pre;
	else	      s2plaq += pre*pre;
      }
      
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
      
    } /* dir2-loop ends */
  } /* dir1-loop ends */

  destroy_r_field(Atmp);
  destroy_r_field(Astpl);
  
  g_doublesum(&ssplaq);
  g_doublesum(&stplaq);
  *splaq=ssplaq/((Real)(3.*nx*ny*nz*nt));
  *tplaq=stplaq/((Real)(3.*nx*ny*nz*nt));

} /* end of plaq() */

/* ************************************************************	*/

