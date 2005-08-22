/************************* RG_path.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
/*****************************************************************************/
/*  Create the following paths:                                              
  paths of length 1 {d} d=x,y,z,t  (4!/(4-1)! possible)                    
  paths of length 2 {d1,d2} d1,d2=x,y,z,t  (4!/(4-2)! possible)             
  paths of length 3 {d1,d2,d3} d1,d2,d3=x,y,z,t  (4!/(4-3)! possible)      
  paths of length 4 {d1,d2,d3,d4} d1,d2,d3,d4=x,y,z,t 4! possible)
  Different permutations of {d_i} reach the some corner:                    
  paths of lenght 1 ( 1 permutation)          
  paths of lenght 2 ( 2! permutation)          
  paths of lenght 3 ( 3! permutation)          
  paths of lenght 4 ( 4! permutation)                                      */
/*****************************************************************************/
#include <stdio.h>
#include <qdp.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"


int f_corner(int rn[RG_Nd],int len)
{
 int c,i,k,value,r[RG_Nd];

//printf("len %d vec=%d %d %d %d\n",len,rn[0],rn[1],rn[2],rn[3]);

 for (i=0;i<RG_Ncn;i++) 
 {
  for (k=0;k<RG_Nd;k++) r[k] = len*rvs[i][k];
  if ((rn[0] == r[0])&&(rn[1] == r[1])&&(rn[2] == r[2])&&(rn[3] == r[3]))
   {
    value = i;
 // printf("i = %d vec=%d %d %d %d\n",i,rvs[i][0],rvs[i][1],rvs[i][2],rvs[i][3]);
    break;
   }

 } 

return value;

}

void RG_path2(QDP_ColorMatrix *w[RG_Ncn],QDP_ColorMatrix *link[RG_Nd],QDP_Sub_Block s,int x[RG_Nd],int cn,int len)
{
int k,rn[RG_Nd];
QDP_Shift offset;
int dir1,dir2;
QLA_Real trace;
QDP_ColorMatrix *temp;

 dir1 = x[0];
 dir2 = x[1];

 temp = QDP_create_M();
 if (temp == NULL) terminate(-1);
 
 QDP_M_eq_zero(w[cn],QDP_all);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
  rn[dir1] = len;

 offset = QDP_create_shift(rn);
 SQDP_M_eq_sM(temp,link[dir2],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],link[dir1],temp,s);

 /* Does not work properly !!!??? */
 //SQDP_M_peq_M_times_sM(w[cn],w[dir1+1],link[dir2],offset,QDP_forward,s);

 QDP_destroy_shift(offset);
 
 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[dir2] = len;
 
 offset = QDP_create_shift(rn);
 SQDP_M_eq_sM(temp,link[dir1],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],link[dir2],temp,s);
 QDP_destroy_shift(offset);


 QDP_destroy_M(temp);


return;
}

void RG_path3(QDP_ColorMatrix *w[RG_Ncn],QDP_ColorMatrix *link[RG_Nd],QDP_Sub_Block s,int x[RG_Nd],int cn,int len)
{
int k,cn2,rn[RG_Nd];
QDP_Shift offset;
int dir1,dir2,dir3;
QDP_ColorMatrix *temp;

 dir1 = x[0];
 dir2 = x[1];
 dir3 = x[2];

 temp = QDP_create_M();
 if (temp == NULL) terminate(-1);
 
 SQDP_M_eq_zero(w[cn],s);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[dir1] = len;
 rn[dir2] = len;

 cn2 = f_corner(rn,len);

 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn2],link[dir3],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[dir3],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn2],temp,s);
 QDP_destroy_shift(offset);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[dir1] = len;
 rn[dir3] = len;

 cn2 = f_corner(rn,len);
 
 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn2],w[dir2+1],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[dir2],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn2],temp,s);
 QDP_destroy_shift(offset);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[dir2] = len;
 rn[dir3] = len;

 cn2 = f_corner(rn,len);
 
 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn2],w[dir1+1],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[dir1],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn2],temp,s);
 QDP_destroy_shift(offset);



 QDP_destroy_M(temp);

return;

}
 
void RG_path4(QDP_ColorMatrix *w[RG_Ncn],QDP_ColorMatrix *link[RG_Nd],QDP_Sub_Block s,int cn,int len)
{
int k, cn3,rn[RG_Nd];
QDP_Shift offset;
int dir1,dir2,dir3,dir4;
QDP_ColorMatrix *temp;

 
 temp = QDP_create_M();
 if (temp == NULL) terminate(-1);
 SQDP_M_eq_zero(w[cn],s);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[0] = len;
 rn[1] = len;
 rn[2] = len;

 cn3 = f_corner(rn,len);

 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn3],w[4],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[3],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn3],temp,s);
 QDP_destroy_shift(offset);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[0] = len;
 rn[1] = len;
 rn[3] = len;

 cn3 = f_corner(rn,len);
 
 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn3],w[3],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[2],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn3],temp,s);
 QDP_destroy_shift(offset);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[0] = len;
 rn[2] = len;
 rn[3] = len;


 cn3 = f_corner(rn,len);
 
 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn3],w[2],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[1],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn3],temp,s);
 QDP_destroy_shift(offset);

 for (k = 0; k < RG_Nd; k++) rn[k] = 0;
 rn[1] = len;
 rn[2] = len;
 rn[3] = len;


 cn3 = f_corner(rn,len);
 
 offset = QDP_create_shift(rn);
 //SQDP_M_peq_M_times_sM(w[cn],w[cn3],w[1],offset,QDP_forward,s);
 SQDP_M_eq_sM(temp,link[0],offset,QDP_forward,s);
 SQDP_M_peq_M_times_M(w[cn],w[cn3],temp,s);
 QDP_destroy_shift(offset);


 QDP_destroy_M(temp);

return;

}

void RG_create_path(QDP_ColorMatrix *pr_wlink[RG_Ncn], QDP_ColorMatrix *link_qdp[RG_Nd], QDP_Sub_Block s,int len)
{
  int i,k,t,path_len,x[RG_Nd];
  int space_only;
  QDP_ColorMatrix *wlink[RG_Ncn];
  QLA_Real trace,norm,c = 1.0;
  QLA_Real fact2 = 1.0/2.0;
  QLA_Real fact3 = 1.0/6.0;
  QLA_Real fact4 = 1.0/24.0;
  QLA_Complex unit;

 
//printf("Start create paths\n"); 
   for (i = 0; i < RG_Ncn; i++)
   {
    wlink[i] = QDP_create_M();
    if (wlink[i] == NULL) terminate(-1);
   }
  
  QLA_C_eq_R(&unit,&c);
  SQDP_M_eq_c(wlink[0],&unit,s);

  for (i=1;i<5;i++)
   SQDP_M_eq_M(wlink[i],link_qdp[i-1],s);
   
  for (i = 5 ; i < RG_Ncn ; i++)
  {
   path_len = 0;
   for (k = 0 ; k < RG_Nd ; k++)
    if (rvs[i][k] != 0) 
     {
     x[path_len] = k;
     path_len += rvs[i][k]; 
     }

   if (path_len == 2)
    RG_path2(wlink,link_qdp,s,x,i,len);

   if (path_len == 3)
    RG_path3(wlink,link_qdp,s,x,i,len);

   if (path_len == 4)
    RG_path4(wlink,link_qdp,s,i,len);

  }
  
  for (i = 0 ; i < RG_Ncn ; i++)
  {
   path_len = 0;
   for (k = 0 ; k < RG_Nd ; k++)
    if (rvs[i][k] != 0) path_len += rvs[i][k];
 
   if ( path_len == 0 ) norm = 1.0; 
   if ( path_len == 1 ) norm = 1.0; 
   if ( path_len == 2 ) norm = fact2; 
   if ( path_len == 3 ) norm = fact3; 
   if ( path_len == 4 ) norm = fact4; 
   
   SQDP_M_eq_r_times_M(pr_wlink[i],&norm,wlink[i],s);
   
/* Comment if check the paths */
   QDP_M_eq_zero(wlink[i],QDP_all);
  }

/* Comment if check the paths........ */

   space_only = RG_Ncn;
   project_qdp(pr_wlink,wlink,&space_only);
 
   for (i = 0; i < RG_Ncn; i++)
    {
    QDP_M_eq_zero(pr_wlink[i],QDP_all);
    SQDP_M_eq_M(pr_wlink[i],wlink[i],s);
    }

/*......................until here */
   for (i = 0; i < RG_Ncn; i++)
    QDP_destroy_M(wlink[i]);

return;

}

