/************************* RG_path.old.c *******************************/
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

shift_v create_shift(int x[4], int n,int len)
{
shift_v d;
int ss,s,v,rn[16][4];

d.x = (int *) malloc(n*sizeof(int));

for (s=0; s<n; s++)
 d.x[s] = x[s];

for (s=0; s< RG_Nd; s++)
 d.s[s] = 0;

for (s=0; s<n; s++)
{
v = x[s];
d.s[v] = 1*len;
}

for (s=0; s<RG_Ncn; s++)
for (ss=0; ss<RG_Nd; ss++)
rn[s][ss] = len * rvs[s][ss];


for (s=0; s<RG_Ncn; s++)
if ((d.s[0] == rn[s][0])&&(d.s[1] == rn[s][1])&&(d.s[2] == rn[s][2])&&(d.s[3] == rn[s][3])) d.rv = s;

return d;
}

int find_count (shift_v *d,int x[4],int n)
{
int s,ss,sf,c;

if ( n == 1) sf = 4;
if ( n == 2) sf = 12;
if ( n == 3) sf = 24;

for (s=0 ; s < sf; s++)
{
c = 0;

for (ss=0 ; ss < n; ss++)
if ( d[s].x[ss] == x[ss]) c++;

if (c == n)  break;

}

return s;
}
 
void RG_create_path(QDP_ColorMatrix *pr_wlink[RG_Ncn], QDP_ColorMatrix *link_qdp[RG_Nd], QDP_Sub_Block s,int len)
{
  int i,j,k,t,x[4];
  int count,c2,space_only;
  QDP_ColorMatrix *path_1[4];
  QDP_ColorMatrix *path_2[12];
  QDP_ColorMatrix *path_3[24];
  QDP_ColorMatrix *path_4[24];
  QDP_ColorMatrix *wlink[RG_Ncn];
  QDP_Shift offset;
  shift_v *d1,*d2,*d3,*d4;
  QLA_Real c = 1.0;
  QLA_Real fact2 = 1.0/2.0;
  QLA_Real fact3 = 1.0/6.0;
  QLA_Real fact4 = 1.0/24.0;
  QLA_Complex unit;



   d1 = (shift_v *) malloc(4*sizeof(shift_v));
   d2 = (shift_v *) malloc(12*sizeof(shift_v));
   d3 = (shift_v *) malloc(24*sizeof(shift_v));
   d4 = (shift_v *) malloc(24*sizeof(shift_v));
  
   for (i = 0; i < RG_Ncn; i++)
    wlink[i] = QDP_create_M();

   for (i = 0; i < 4; i++)
    path_1[i] = QDP_create_M();
   for (i = 0; i < 4; i++)
    SQDP_M_eq_M(path_1[i],link_qdp[i],s);

   for (i = 0; i < 12; i++)
    path_2[i] = QDP_create_M();

   for (i = 0; i < 24; i++)
   {
   path_3[i] = QDP_create_M();
   path_4[i] = QDP_create_M();
   }


//   printf("Start building paths %d\n",this_node); fflush(stdout);
   for (i = 0; i < RG_Nd; i++)
   {
    x[0] = i; 
    d1[i] = create_shift(x,1,len); 
   }
//   printf("First shift %d\n",this_node);fflush(stdout);

   count = 0;
   for (i = 0; i < RG_Nd; i++)
   {
   x[0] = i;
   c2 = find_count(d1,x,1);
   offset = QDP_create_shift(d1[c2].s);
   for (j = 0; j < RG_Nd ; j++) if ( j != i)
   {
   x[1] = j;
   d2[count] = create_shift(x,2,len); 
   SQDP_M_eq_M_times_sM(path_2[count],path_1[c2],link_qdp[j],offset,QDP_forward,s);
   count ++;
   }
   QDP_destroy_shift(offset);
   }
  // printf("Second shift %d\n",this_node);fflush(stdout);

   

   count = 0;
   for (i = 0; i < RG_Nd; i++)
   for (j = 0; j < RG_Nd; j++) if (j != i)
   {
   x[0] = i; x[1] = j;
   c2 = find_count(d2,x,2);
   offset = QDP_create_shift(d2[c2].s);
   for (k = 0; k < RG_Nd; k++) if (k != i) if (k != j)
   {
   x[2] = k;
   d3[count] = create_shift(x,3,len);
   SQDP_M_eq_M_times_sM(path_3[count],path_2[c2],link_qdp[k],offset,QDP_forward,s);
   count++;
   }
   QDP_destroy_shift(offset);
   }
//   printf("Third shift %d\n",this_node);fflush(stdout);

   count = 0;
   for (i = 0; i < RG_Nd; i++)
   for (j = 0; j < RG_Nd; j++) if (j != i)
   for (k = 0; k < RG_Nd; k++) if (k != i) if (k != j)
   {
   x[0] = i; x[1] = j; x[2] = k;
   c2 = find_count(d3,x,3);
   offset = QDP_create_shift(d3[c2].s);
   for (t = 0; t < RG_Nd; t++) if (t != i) if (t != j) if (t != k)
   {
   x[3] = t;
   d4[count] = create_shift(x,4,len);
   SQDP_M_eq_M_times_sM(path_4[count],path_3[c2],link_qdp[t],offset,QDP_forward,s);
   count++;
   }
   QDP_destroy_shift(offset);
   }
//   printf("Fourth shift %d\n",this_node);fflush(stdout);

   QLA_C_eq_R(&unit,&c);
   SQDP_M_eq_c(wlink[0],&unit,s);

   for (i=1;i<5;i++)
   SQDP_M_eq_M(wlink[i],path_1[i-1],s);
   
   for (i=5;i<RG_Ncn;i++)
   {
   SQDP_M_eq_zero(wlink[i],s);

   for (j=0;j<12;j++) if(d2[j].rv == i)
    SQDP_M_peq_r_times_M(wlink[i],&fact2,path_2[j],s);

   for (j=0;j<24;j++) if(d3[j].rv == i)
    SQDP_M_peq_r_times_M(wlink[i],&fact3,path_3[j],s);

   for (j=0;j<24;j++) if(d4[j].rv == i)
    SQDP_M_peq_r_times_M(wlink[i],&fact4,path_3[j],s);

   }

   
   space_only = RG_Ncn;
//   printf("projection %d\n",this_node);fflush(stdout);
   project_qdp(wlink, pr_wlink,&space_only);
 
   for (i = 0; i < 4; i++)
    QDP_destroy_M(path_1[i]);
   
   for (i = 0; i < 12; i++)
    QDP_destroy_M(path_2[i]);

   for (i = 0; i < 24; i++)
   {
    QDP_destroy_M(path_3[i]);
    QDP_destroy_M(path_4[i]);
   }

   for (i = 0; i < RG_Ncn; i++)
    QDP_destroy_M(wlink[i]);

   free(d1);
   free(d2);
   free(d3);
   free(d4);


return;

}

