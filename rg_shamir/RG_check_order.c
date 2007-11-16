/************************* RG_check_order.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */

#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

#define TOL 10e-6

void RG_gauge_paths(QDP_ColorMatrix *wlink[NRG][RG_Ncn], 
		    QDP_ColorMatrix *rg_link[NRG][RG_Nd],
		    QDP_Sub_Block QDP_block[NRG+1])
{
int i,j,len;
QDP_ColorMatrix *unit,*link_qdp[RG_Nd];
QLA_Real trace, c = 1.0;
QLA_Complex one;

  unit = QDP_create_M();
  QLA_C_eq_R(&one,&c);  
  QDP_M_eq_c(unit, &one, QDP_all);

  for(i=0; i< RG_Nd; ++i)
   {
    link_qdp[i] = QDP_create_M();
    set_M_from_site(link_qdp[i],F_OFFSET(link[i]),EVENANDODD);
   }

  for(j=0; j<nrg; ++j)
  for(i=0; i<RG_Nd; ++i)
   {
    rg_link[j][i] = QDP_create_M();
    SQDP_M_eq_zero(rg_link[j][i],QDP_block[j+1]);
   }


  RG_gauge(rg_link,link_qdp,QDP_block);
  
  for(i=0; i<RG_Nd; ++i)
    QDP_destroy_M(link_qdp[i]);

/* Create paths W(2x,2x+r) that connect the 16 corners of the hypercube */
  for (i=0; i < nrg; i++)
  {
  len = intpow(2,nrg-i-1);

  printf("node %d:paths from lattice %d to lattice %d\n",this_node,2*len,len);
  fflush(stdout);

   for (j = 0; j < RG_Ncn; j++)
    wlink[i][j] = QDP_create_M();
   
  RG_create_path(wlink[i],rg_link[i],QDP_block[i],len);
/*
  for (j = 0; j < RG_Ncn; j++)
   {
   SQDP_r_eq_re_M_dot_M(&trace,unit,wlink[i][j],QDP_block[i]); 
   printf("trace for wlink %d %d is %lf\n",i,j,trace);
   }
*/

  }


  
return;

}

QLA_Real RG_check_order(QDP_ColorMatrix *wlink,QDP_ColorMatrix *link_qdp[RG_Nd],QDP_Sub_Block s, int len,int cn)
{
int j,k,rn,sum,pow2,status=0;
int v[RG_Nd][RG_Nd],v_s[RG_Nd];
double vol;
QLA_Real trace,c =1.0;
QLA_Complex one;
QDP_ColorMatrix *unit,*temp0,*temp1,*temp;
QDP_Shift offset;


  vol = (double)(nx*ny*nz*nt);
  pow2 = intpow(2.0,(double)(len-1));
  vol /= pow(2.0*(double)len,4.0);
  //printf("vol %lf\n",vol);

  temp = QDP_create_M();
  temp1 = QDP_create_M();
  temp0 = QDP_create_M();

   
  unit = QDP_create_M();
  QLA_C_eq_R(&one,&c);  
  SQDP_M_eq_c(unit, &one, s);
    
   for (j = 0; j < RG_Nd; j++)
   {
    rn = len * rvs[cn][j];
    for (k = 0; k < RG_Nd; k++) v[j][k] = 0;
    v[j][j] = rn;
    }
/*  
  printf("corner %d\n",cn);
  printf("(%d %d %d %d)\n",rvs[cn][0],rvs[cn][1],rvs[cn][2],rvs[cn][3]);
  printf("\n");
  for (j = 0; j < RG_Nd; j++)
   printf("(%d %d %d %d)\n",v[j][0],v[j][1],v[j][2],v[j][3]);
  printf("\n");
*/
   SQDP_M_eq_M(temp,unit,s); 
   

   if (cn != 0 )
   {
   for (k = 0; k < RG_Nd; k++) v_s[k] = 0;
   len*rvs[cn][k];

   for (j = 0; j < RG_Nd; j++)
   {

    if( v[j][j] != 0 )
    {
    
 /*   
    printf("Multiply \n",j);
    printf("U_%d(%d %d %d %d)\n",j,v_s[0],v_s[1],v_s[2],v_s[3]);
 */

    offset = QDP_create_shift(v_s); 
    SQDP_M_eq_M_times_sM(temp1,temp,link_qdp[j],offset,QDP_forward,s);
    QDP_destroy_shift(offset);
    SQDP_M_eq_M(temp,temp1,s); 
    
    for (k = 0; k < RG_Nd; k++) v_s[k] +=  v[j][k];
    }

   }
  // printf("\n");
  
     
    sum = 0;
    for (k = 0; k < RG_Nd; k++) sum += (v_s[k]-len*rvs[cn][k]);

    if (sum != 0 )
    {
    printf("Error: Loop not completed\n");
    exit(-1);
    } 
   
   } 
    SQDP_r_eq_re_M_dot_M(&trace,temp,wlink,s); 
//    SQDP_r_eq_re_M_dot_M(&trace,temp,wlink,s); 
  
  QDP_destroy_M(unit);
  QDP_destroy_M(temp);
  QDP_destroy_M(temp0);
  QDP_destroy_M(temp1);

 
return (trace/vol);

}


int RG_check_path()
{
  int i,j,status=0;
  int len;
  QLA_Real check[NRG][RG_Ncn],checkg[NRG][RG_Ncn];
  QDP_Sub_Block QDP_block[NRG+1];
  QDP_ColorMatrix *wlink[NRG][RG_Ncn], *rg_link[NRG][RG_Nd];
  
  for (i = 0; i < nrg+1; i++)
    RG_create_block(&QDP_block[i],i);
  
  //  test(QDP_block);
  
  printf("START CHECK:-----------------order paths-----------------\n");
  for (i = 0; i < nrg; i++)
    {
      
      for (j = 0; j < RG_Ncn; j++)
	{
	  wlink[i][j] = QDP_create_M();
	  QDP_M_eq_zero(wlink[i][j],QDP_all);
	}
      
      for (j = 0; j < RG_Nd; j++)
	{
	  rg_link[i][j] = QDP_create_M();
	  QDP_M_eq_zero(rg_link[i][j],QDP_all);
	}
      
    }
  
  
  rephase(OFF);
  RG_gauge_paths(wlink,rg_link,QDP_block);
  
  printf("paths created\n");
  for (i = 0; i < nrg; i++)
    {
      len = intpow(2,nrg-i-1);
      for (j = 0; j < RG_Ncn; j++)
	check[i][j] = RG_check_order(wlink[i][j],
				     rg_link[i],
				     QDP_block[i],len,j);
    }
  printf("Traces evaluated\n");
  
  for (i = 0; i < nrg; i++)
    {
      
      for (j = 0; j < RG_Ncn; j++)
	QDP_M_eq_zero(wlink[i][j],QDP_all);
      
      for (j = 0; j < RG_Nd; j++)
	QDP_M_eq_zero(rg_link[i][j],QDP_all);
      
    }
  
  rand_gauge(F_OFFSET(rgt));
  RG_gauge_paths(wlink,rg_link,QDP_block);
  rephase(ON);
  
  for (i = 0; i < nrg; i++)
    {
      len = intpow(2,nrg-i-1);
      for (j = 0; j < RG_Ncn; j++)
	checkg[i][j] = RG_check_order(wlink[i][j],
				      rg_link[i],
				      QDP_block[i],len,j);
    }
  
  for (i = 0; i < nrg; i++)
    {
      len = intpow(2,nrg-i-1);
      for (j = 0; j < RG_Ncn; j++)
	{
	  printf("For len = %d and paths %d diff:  %.9e = %.9e\n",len,j,checkg[i][j],check[i][j]);
	  if ( fabs(checkg[i][j] - check[i][j] ) > TOL )
	    {
	      //printf("For nrg = %d and paths %d diff: %e=%e\n",i,j,checkg[i][j],check[i][j]);
	      status = 1;
	    }
	}
    }
  
  printf("END CHECK:-----------------order paths-----------------\n");
  
  return status;
} 
