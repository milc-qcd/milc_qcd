/************************* RG_trace_new.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#include "RG_include.h"

/*
    spin_taste_table S T

    prints the eta/etap matrix for the specified gamma_S and gamma_T */

/* Encode gamma matrices as a single hex digit */

/* Convert hex encoding to four integers */
void gamma_to_int(RG_gamma g, int r[4])
{
  int i;


  for(i = 0; i < 4; i++){
    r[i] = (g & 8)/8;
    g <<= 1;
  }

}


RG_gamma int_to_gamma(int g)
{
  int i,j,c,v[4],vs[4];

/* the order I want is (r[3],r[0],r[1],r[2]) */
for(i = 0; i < 4; i++)
    vs[i] = rvs[g][i];
  
/*
for(i = 0; i < 3; i++)
    vs[i+1] = rvs[g][i];
   
  vs[0] = rvs[g][3];
*/

  for(i=0;i<16;i++)
   {
    gamma_to_int(i,v);
   
    if ((v[0]==vs[0])&&(v[1]==vs[1])&&(v[2]==vs[2])&&(v[3]==vs[3])) 
    {
    return i;
    break;
    }
  }
}


gamma_phase gamma_x_gamma(gamma_phase *g1, gamma_phase *g2)
{
  int r[4], s[4];
  gamma_phase g12;
  int phase12;

  g12.g = g1->g ^ g2->g;
  g12.sign = g1->sign * g2->sign;  /* To be corrected below */

  gamma_to_int(g1->g, r);
  gamma_to_int(g2->g, s);

  phase12 = (s[3]*(r[2]+r[1]+r[0]) + s[0]*(r[2] + r[1]) + s[1]*r[2]) % 2;
  //phase12 = (s[0]*(r[3]+r[2]+r[1]) + s[1]*(r[3] + r[2]) + s[2]*r[3]) % 2;
  if(phase12 == 1)g12.sign = -g12.sign;

  return g12;
}


gamma_phase adjoint(gamma_phase *g)
{
  gamma_phase gadj;
  int phaseadj;
  int r[4];

  gadj.g = g->g;
  gadj.sign = g->sign;

  gamma_to_int(g->g, r);
  
  phaseadj = (r[0]*(r[1] + r[2] + r[3]) + r[1]*(r[2] + r[3]) + r[2]*r[3]) % 2;

  if(phaseadj == 1)gadj.sign = -gadj.sign;

  return gadj;
}
void RG_trace(int *sign,int reta[4],int cmp[4])
{

  RG_gamma S, T, eta, etap;
  gamma_phase gS, gadjS, gT, geta, getap, gtmp1, gtmp2;
  int rS[4], rT[4],radjS[4],retap[4];

  S = cmp[0];
  T = cmp[1];

  etap =  int_to_gamma(cmp[2]);

  gamma_to_int(S, rS);
  gamma_to_int(T, rT);

  gS.g = S;
  gS.sign = 1;
  gT.g = T;
  gT.sign = 1;


  gadjS = adjoint(&gS);
  gamma_to_int(gadjS.g, radjS);
    getap.g = etap;
    getap.sign = 1;
    gtmp1 = gamma_x_gamma(&gadjS, &getap);
    gtmp2 = gamma_x_gamma(&gtmp1, &gT);

    gamma_to_int(etap, retap);

/*
fprintf(stderr,"(%d,%d,%d,%d) == (%d,%d,%d,%d) \n",
	   rvs[cmp[2]][0],rvs[cmp[2]][1],rvs[cmp[2]][2],rvs[cmp[2]][3],
	   retap[0],retap[1],retap[2],retap[3]); 
 */ 

  gamma_to_int(gtmp2.g, reta);

/*
fprintf(stderr,"Tr[ S^dag(%d,%d,%d,%d) etap(%d,%d,%d,%d) T(%d %d %d %d) eta^dag(%d %d %d %d) ] = %2d\n",
	   rS[0],rS[1],rS[2],rS[3],
	   retap[0],retap[1],retap[2],retap[3], 
	   rT[0],rT[1],rT[2],rT[3],
	   reta[0],reta[1],reta[2],reta[3],
	   gtmp2.sign
	   );
 */
  *sign = gtmp2.sign;

return; 
}

void RG_trace1(int *sign,int cmp[4])
{

  RG_gamma S, T, eta, etap;
  gamma_phase gS, gadjS, gT, geta, getap, gadjeta,gtmp1, gtmp2,gtmp3;
  int rS[4], rT[4],radjS[4],retap[4],radjeta[4],reta[4],r[4];

  S = cmp[0];
  T = cmp[1];

  etap =  int_to_gamma(cmp[2]);
  eta =  int_to_gamma(cmp[3]);

  gamma_to_int(S, rS);
  gamma_to_int(T, rT);
  gamma_to_int(eta, reta);

  gS.g = S;
  gS.sign = 1;
  gT.g = T;
  gT.sign = 1;
  getap.g = etap;
  getap.sign = 1;
  geta.g = eta;
  geta.sign = 1;

  gadjS = adjoint(&gS);
  gamma_to_int(gadjS.g, radjS);
  gtmp1 = gamma_x_gamma(&gadjS, &getap);
  gtmp2 = gamma_x_gamma(&gtmp1, &gT);

  gadjeta = adjoint(&geta);
  gamma_to_int(gadjeta.g, radjeta);
  gtmp3 = gamma_x_gamma(&gtmp2, &gadjeta);

/*
  fprintf(stderr,"%d (%d,%d,%d,%d) == (%d,%d,%d,%d) =  (%d,%d,%d,%d) \n",
	   gadjeta.sign,radjeta[0],radjeta[1],radjeta[2],radjeta[3],
	   reta[0],reta[1],reta[2],reta[3], 
	   rvs[cmp[3]][0],rvs[cmp[3]][1],rvs[cmp[3]][2],rvs[cmp[3]][3] );
*/
/*
fprintf(stderr,"(%d,%d,%d,%d) == (%d,%d,%d,%d) \n",
	   rvs[cmp[2]][0],rvs[cmp[2]][1],rvs[cmp[2]][2],rvs[cmp[2]][3],
	   retap[0],retap[1],retap[2],retap[3]); 
 */ 

  gamma_to_int(etap, retap);
  gamma_to_int(gtmp3.g, r);


 if (gtmp3.g == 0) 
 {
 fprintf(stderr,"Tr[ S^dag(%d,%d,%d,%d) etap(%d,%d,%d,%d) T(%d %d %d %d) etap^dag(%d %d %d %d) ] = %2d(%d,%d,%d,%d)\n",
	   radjS[0],radjS[1],radjS[2],radjS[3],
	   retap[0],retap[1],retap[2],retap[3], 
	   rT[0],rT[1],rT[2],rT[3],
	   radjeta[0],radjeta[1],radjeta[2],radjeta[3],
	   gtmp3.sign,r[0],r[1],r[2],r[3]);

  *sign = gtmp3.sign;
 }
  else
   *sign = 0;

return; 
}
