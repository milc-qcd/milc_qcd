/***************** twopt.c ********************************************/

/* Computes meson propagators from specified quark propagators */
/* MIMD version 7 */

#include "ks_hl_spectrum_includes.h"
#include <time.h>

#define GAMMAFIVE -1
#define eve  0
#define od   1
#define UNDEFINED  2

/* Computes 
   Tr_sc(Gamma_src * W^+(x,y,z,t) * Gamma_snk * q^+(x,0)* Q_HL(x,0),
   where q is a staggered propagator and Q_HL could be the heavy quark,
   at a given space-time point */

static complex  
KS_2pt_trace(su3_matrix * antiquark, wilson_propagator * quark, 
		      int * g_snk, int n_snk, int *g_src, int n_src, int *p, site *s)
{
  int t;
  int my_x;
  int my_y;
  int my_z;
  
  complex trace;
  int s0;
  int c0,c1,i;

  wilson_propagator temp,temp1;
  su3_matrix mat, mat1;
  

  t = s->t;
  my_x = s->x;
  my_y = s->y;
  my_z = s->z;
  
  temp = *quark;

  //multiply by gamma_snk

   for(i=0;i<n_snk;i++)
     for(c0=0;c0<3;c0++){
      mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), g_snk[i]);
      temp.c[c0] = temp1.c[c0]; 
    } 
   
   //multiply by Omega field
   if((t % 2) == 1)
     for(c0=0;c0<3;c0++){
       mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), TUP);
       temp.c[c0] = temp1.c[c0]; 
     }
    
   if((my_x % 2) == 1)
     for(c0=0;c0<3;c0++){
       mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), XUP); 
       temp.c[c0] = temp1.c[c0]; 
     }
   
   if((my_y % 2) == 1)
     for(c0=0;c0<3;c0++){
       mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), YUP); 
       temp.c[c0] = temp1.c[c0]; 
     }
   
   if((my_z % 2) == 1)
     for(c0=0;c0<3;c0++){
       mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), ZUP);  
       temp.c[c0] = temp1.c[c0]; 
     } 
 
   //mulptiply by gamma_src
   for(c0=0;c0<3;c0++)
     for(i=0;i<n_src;i++)
       {
	 mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), g_src[i]);   
	 temp.c[c0] = temp1.c[c0];
       }
   
   for(c0=0;c0<3;c0++) 
     for(c1=0;c1<3;c1++){
       trace.real = 0.0;
       trace.imag = 0.0;
       
       for(s0=0;s0<4;s0++){
	 trace.real += temp.c[c0].d[s0].d[s0].c[c1].real;
	 trace.imag += temp.c[c0].d[s0].d[s0].c[c1].imag;
       }
       
       mat.e[c0][c1].real = trace.real;
       mat.e[c0][c1].imag = trace.imag;
       
     }
   
   mult_su3_na(&mat, antiquark, &mat1); //antiquark is just the staggered prop su3 matrix
   
   trace = trace_su3(&mat1);

   return(trace);
}


/* Summation on spatial slice. Caution: elements of prop are not set initially to zero here! 
   In prop averages of propagators are accumulated!  */

static void 
KS_2pt_func(field_offset snk, field_offset src, int *g_snk, int n_snk,
	    int *g_src, int n_src, int *p, complex *prop, int parity)
{
  int i, t, my_x,my_y,my_z;
  site *s;
  complex trace, trace1, epx,epx1;
  Real pi, mom[3];

  su3_matrix *antiquark;
  wilson_propagator  *quark;


  /* printf("incoming propagator\n");
     for(i=0;i<nt;i++) printf("%.9e %.9e\n", prop[i].real,prop[i].imag);*/

  pi = 4.0 * atan( 1.);
  mom[0] = -2.*pi/(Real)nx;  
  mom[1] = -2.*pi/(Real)ny;  
  mom[2] = -2.*pi/(Real)nz; 
 
  FORALLSITES(i,s){
    t = s->t;
    my_x = s->x;
    my_y = s->y;
    my_z = s->z;
  
    epx1.real = 0.0;
    epx1.imag = (mom[0]*(Real)p[0]*(Real)my_x + mom[1]*(Real)p[1]*(Real)my_y 
		+ mom[2]*(Real)p[2]*(Real)my_z);
    epx = cexp(&epx1);

    antiquark = (su3_matrix *)F_PT(s, snk);
    quark = (wilson_propagator *)F_PT(s,src);
    trace = KS_2pt_trace(antiquark, quark, g_snk,n_snk, g_src,n_src, p, s);
    
    if((t+ my_x + my_y + my_z)%2 == 1) {
       trace.real *= -1.0;
       trace.imag *= -1.0;
    }
    if (parity == UNDEFINED){
      trace1.real = trace.real*epx.real - trace.imag*epx.imag;
      trace1.imag = trace.real*epx.imag + trace.imag*epx.real;
    }
    else{
      if(parity == eve) {
	trace1.real = trace.real*epx.real; 
	trace1.imag = trace.imag*epx.real;

      }
      else {
	if(parity == od) {
	trace1.real = - trace.imag*epx.imag;
	trace1.imag =   trace.real*epx.imag;
	}
	else {printf("Unknown parity\n");terminate(1);}
      }
    }
    prop[t].real += trace1.real;
    prop[t].imag += trace1.imag;
  }

  /* printf("outgoing propagator\n");
     for(i=0;i<nt;i++) printf("%.9e %.9e\n", prop[i].real,prop[i].imag);*/
}

/* Caution: elements of propagator are not set initially to zero here! */

static void 
All_KS_hl_prop(field_offset snk, field_offset src, complex **propagator)
{

  int i;
  int g_snk[4],g_src[4];
  int n_snk,n_src;
  Real temp;

  int p000[3] ={0,0,0}; 

  int p100[3][3] = {{1,0,0},  //0
		    {0,1,0},  //1
		    {0,0,1}}; //2


  int p110[9][3] = {{1,1,0},  //0
		    {1,0,1},  //1
		    {1,0,-1}, //2
		    {1,-1,0}, //3
		    {0,1,1},  //4
		    {0,1,-1}, //5
		    {-1,1,0}, //add 6
		    {0,-1,1}, //add 7
		    {-1,0,1}}; //add 8
  
  int p111[7][3] = {{1,1,1},  //0
		    {1,1,-1}, //1
		    {1,-1,1}, //2
		    {1,-1,-1}, //3
		    {-1,1,1}, //add 4
		    {-1,1,-1}, //add 5
		    {-1,-1,1}}; //add 6

  int p200[3][3] = {{2,0,0}, //0
		    {0,2,0}, //1
		    {0,0,2}}; //2
  
  int p210[12][3] = {{2,1,0}, //0
		     {2,0,1}, //1
		     {2,0,-1}, //2
		     {2,-1,0}, //3
		     {1,2,0},  //4
		     {1,0,2},  //5
		     {1,0,-2}, //6
		     {1,-2,0}, //7
		     {0,2,1},  //8
		     {0,2,-1}, //9 
		     {0,1,2}, //10
		     {0,1,-2}}; //11

  int p211[12][3] = {{2,1,1}, //0
		     {2,1,-1}, //1
		     {2,-1,1}, //2
		     {2,-1,-1}, //3
		     {1,2,1},  //4
		     {1,2,-1}, //5
		     {1,1,2},  //6
		     {1,1,-2}, //7
		     {1,-1,2}, //8
		     {1,-1,-2}, //9
		     {1,-2,1}, //10
		     {1,-2,-1}}; //11

  int p220[6][3] = {{2,2,0},  //0
		    {2,0,2},  //1
		    {2,0,-2}, //2
		    {2,-2,0}, //3
		    {0,2,2},  //4
		    {0,2,-2}}; //5
  
  int p300[3][3] = {{3,0,0},  //0
		    {0,3,0},  //1
		    {0,0,3}}; //2

  int p221[12][3] = {{2,2,1},  //0
		     {2,2,-1}, //1
		     {2,1,2},  //2
		     {2,1,-2}, //3
		     {2,-1,2}, //4
		     {2,-1,-2}, //5
		     {2,-2,1},  //6
		     {2,-2,-1}, //7
		     {1,2,2},   //8
		     {1,2,-2},  //9
		     {1,-2,2},  //10
		     {1,-2,-2}}; //11



  int p400[3][3] = {{4,0,0}, //0
		    {0,4,0}, //1
		    {0,0,4}}; //2



  /*======= (pi) P5_P5 p000 ===================*/
  n_snk=1; n_src=1;
  g_snk[0]=GAMMAFIVE;
  g_src[0]=GAMMAFIVE;

  KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p000, propagator[0], UNDEFINED);
  
  
  //*======= (pi) P5_P5 p100 ===================*/

  for(i=0;i<3;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p100[i], propagator[1], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[1][i].real = propagator[1][i].real/3.0;
      propagator[1][i].imag = propagator[1][i].imag/3.0;
    }
  
  /*======= (pi) P5_P5 p110 ===================*/

  for(i=0;i<6;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p110[i], propagator[2], eve);

  for(i=0;i<nt;i++) 
    {
      propagator[2][i].real = propagator[2][i].real/6.;
      propagator[2][i].imag = propagator[2][i].imag/6.;
    }

  /*======= (pi) P5_P5 p111 ===================*/

  for(i=0;i<4;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p111[i], propagator[3], eve);

  for(i=0;i<nt;i++) 
    {
      propagator[3][i].real = propagator[3][i].real/4.;
      propagator[3][i].imag = propagator[3][i].imag/4.;
    }

  /*======= (pi) P5_P5 p200 ===================*/

  for(i=0;i<3;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p200[i], propagator[4], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[4][i].real = propagator[4][i].real/3.;
      propagator[4][i].imag = propagator[4][i].imag/3.;
    }
  
  /*======= (pi) P5_P5 p210 ===================*/

  for(i=0;i<12 ;i++)
    KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[5], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[5][i].real = propagator[5][i].real/12.;
      propagator[5][i].imag = propagator[5][i].imag/12.;
    }
  
  /*======= (pi) P5_P5 p211 ===================*/

  for(i=0;i<12 ;i++)
    KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p211[i], propagator[6], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[6][i].real = propagator[6][i].real/12.;
      propagator[6][i].imag = propagator[6][i].imag/12.;
    }

  /*======= (pi) P5_P5 p220 ===================*/  
 
  for(i=0;i<6;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p220[i], propagator[7], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[7][i].real = propagator[7][i].real/6;
      propagator[7][i].imag = propagator[7][i].imag/6;
    }

  /*======= (pi) P5_P5 p300 ===================*/
  
  for(i=0;i<3;i++)
    KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p300[i], propagator[8], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[8][i].real = propagator[8][i].real/3;
      propagator[8][i].imag = propagator[8][i].imag/3;
    }
  
  /*======= (pi) P5_P5 p221 ===================*/

  for(i=0;i<12 ;i++)
    KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p221[i], propagator[9], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[9][i].real = propagator[9][i].real/12;
      propagator[9][i].imag = propagator[9][i].imag/12;
    }

  /*======= (pi) P5_P5 p400 ===================*/

  for(i=0;i<3;i++)
    KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[10], eve);

  
  for(i=0;i<nt;i++) 
    {
      propagator[10][i].real = propagator[10][i].real/3;
      propagator[10][i].imag = propagator[10][i].imag/3;
    }
  
  /*======= A4_P5 p000 ===================*/
  g_src[0] = GAMMAFIVE ;
  n_src = 1;
  g_snk[1] = GAMMAFIVE ;
  g_snk[0] = TUP;
  n_snk = 2;

   KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p000, propagator[11], UNDEFINED);

   //mulply by i
    for(i=0;i<nt;i++) 
     { 
       temp = propagator[11][i].real;
       propagator[11][i].real = -propagator[11][i].imag;
       propagator[11][i].imag = temp;
     }
  

    /*======= A4_P5 p100 ===================*/
   for(i=0;i<3;i++)
     KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p100[i], propagator[12], eve);
 
   //mulply by i and average over momenta
    for(i=0;i<nt;i++) 
     {
       temp = propagator[12][i].real;
       propagator[12][i].real = -propagator[12][i].imag/3.;
       propagator[12][i].imag = temp/3.;
     }
  
   /*======= A4_P5 p110 ===================*/

   for(i=0;i<6;i++)
     KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p110[i], propagator[13], eve);  
   for(i=0;i<nt;i++) 
     {
       temp = propagator[13][i].real;
       propagator[13][i].real = -propagator[13][i].imag/6;
       propagator[13][i].imag = temp/6;
     }
   
   /*======= A4_P5 p111 ===================*/

   for(i=0;i<4;i++)
     KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p111[i], propagator[14], eve);

   for(i=0;i<nt;i++) 
     {
       temp = propagator[14][i].real;
       propagator[14][i].real = -propagator[14][i].imag/4.;
       propagator[14][i].imag = temp/4.;
     }

   /*======= A4_P5 p200 ===================*/

   for(i=0;i<3;i++)
     KS_2pt_func(snk, src, g_snk, n_snk, g_src, n_src, p200[i], propagator[15], eve);
   
  
   for(i=0;i<nt;i++) 
     {
       temp = propagator[15][i].real;
       propagator[15][i].real = -propagator[15][i].imag/3.;
       propagator[15][i].imag = temp/3.;
     }
  

   /*======= A1_P5 p100 ===================*/
   g_src[0] = GAMMAFIVE ;
   n_src = 1;
   g_snk[1] = GAMMAFIVE;
   g_snk[0] = XUP;
   n_snk = 2;

   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[0], propagator[16], od);
  
   g_snk[0] = YUP;
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[1], propagator[16], od);
  
   g_snk[0] = ZUP;
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[2], propagator[16], od);
  
   for(i=0;i<nt;i++) 
     {
       propagator[16][i].real /= -3.;
       propagator[16][i].imag /= -3.;
     }
  
   /*======= A1_P5 p110 ===================*/

   g_snk[0] = XUP;
   for(i=0;i<4;i++)
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[i], propagator[17], od);

   g_snk[0] = YUP;
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[0], propagator[17], od); 
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[4], propagator[17], od);
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[5], propagator[17], od);
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[6], propagator[17], od);

   g_snk[0] = ZUP;
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[1], propagator[17], od);
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[4], propagator[17], od);
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[7], propagator[17], od);
   KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[8], propagator[17], od);

   for(i=0;i<nt;i++) 
     {
       propagator[17][i].real = -propagator[17][i].real/12;
       propagator[17][i].imag = -propagator[17][i].imag/12;
     }

   /*======= A1_P5 p111 ===================*/
   g_snk[0] = XUP;
   for(i=0;i<4;i++)
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[18], od);

   g_snk[0] = YUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[0], propagator[18], od); 
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[1], propagator[18], od);
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[4], propagator[18], od);
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[5], propagator[18], od);

   g_snk[0] = ZUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[0], propagator[18], od);
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[2], propagator[18], od);
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[4], propagator[18], od);
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[6], propagator[18], od);

     for(i=0;i<nt;i++) 
       {
	 propagator[18][i].real = -propagator[18][i].real/12;
	 propagator[18][i].imag = -propagator[18][i].imag/12;
       }

     /*======= A1_P5 p200 ===================*/
     g_snk[0] = XUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[0], propagator[19], od);

     g_snk[0] = YUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[1], propagator[19], od);

     g_snk[0] = ZUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[2], propagator[19], od);

     for(i=0;i<nt;i++) 
       {
	 propagator[19][i].real = -propagator[19][i].real/3;
	 propagator[19][i].imag = -propagator[19][i].imag/3;
       }

     /*======= P5_A1 p100 ===================*/

     g_snk[0] = GAMMAFIVE ;
     n_snk = 1;
     g_src[1] = XUP;
     g_src[0] = GAMMAFIVE;
     n_src = 2;

     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[0], propagator[20], od);

     g_snk[0] = YUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[1], propagator[20], od);

     g_snk[0] = ZUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[2], propagator[20], od);

     for(i=0;i<nt;i++) 
       {
	 propagator[20][i].real = propagator[20][i].real/3;
	 propagator[20][i].imag = propagator[20][i].imag/3;
       }

     /*======= (ro) V1_V1 p000 ===================*/

     g_snk[0] = XUP ;
     n_snk = 1;
     g_src[0] = XUP;
     n_src = 1;

     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);

     g_snk[0] = YUP;  g_src[0] = YUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);
 

     /*======= (ro) V1_V1 p100 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);

     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);
 
     g_snk[0] = ZUP;  g_src[0] = ZUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[22][i].real /= 3.;
	 propagator[22][i].imag /= 3.;
       }

     /*======= (ro) V1_V1 p110 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP;
    for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[23][i].real /= 6;
	propagator[23][i].imag /= 6;
      }

    /*======= (ro) V1_V1 p111 ===================*/
     g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<4;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);

   g_snk[0] = YUP;  g_src[0] = YUP; 
   for(i=0;i<4;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);

   g_snk[0] = ZUP;  g_src[0] = ZUP;
   for(i=0;i<4;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);


    for(i=0;i<nt;i++) 
      {
	propagator[24][i].real = propagator[24][i].real/4;
	propagator[24][i].imag = propagator[24][i].imag/4;
      }

    /*======= (ro) V1_V1 p200 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<3;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<3;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP; 
    for(i=0;i<3;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);
    
    for(i=0;i<nt;i++) 
      {
	propagator[25][i].real = propagator[25][i].real/3;
	propagator[25][i].imag = propagator[25][i].imag/3;
      }

    /*======= (ro) V1_V1 p210 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP; 
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP; 
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[26][i].real = propagator[26][i].real/12;
	propagator[26][i].imag = propagator[26][i].imag/12;
      }

    /*======= (ro) V1_V1 p211 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP; 
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p211[i], propagator[27], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP; 
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p211[i], propagator[27], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<12;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p211[i], propagator[27], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[27][i].real = propagator[27][i].real/12;
	propagator[27][i].imag = propagator[27][i].imag/12;
      }

    /*======= (ro) V1_V1 p220 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p220[i], propagator[28], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP; 
     for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p220[i], propagator[28], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP;
     for(i=0;i<6;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p220[i], propagator[28], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[28][i].real = propagator[28][i].real/6;
	 propagator[28][i].imag = propagator[28][i].imag/6;
       }

    /*======= (ro) V1_V1 p300 ===================*/
    
     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);
    
     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[29][i].real = propagator[29][i].real/3;
	 propagator[29][i].imag = propagator[29][i].imag/3;
       }
     
    /*======= (pi) V1_V1 p221 ===================*/
    
    g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p221[i], propagator[30], eve);
   
    g_snk[0] = YUP;  g_src[0] = YUP; 
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p221[i], propagator[30], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP;
    for(i=0;i<12;i++)
      KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p221[i], propagator[30], eve);
    
    for(i=0;i<nt;i++) 
      {
	propagator[30][i].real = propagator[30][i].real/12;
	propagator[30][i].imag = propagator[30][i].imag/12;
      }
    
   /*======= (ro) V1_V1 p400 ===================*/
    
     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[31], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src , p400[i], propagator[31], eve);
    
     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<3;i++)
       KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[31], eve);
     
     for(i=0;i<nt;i++) 
       {
	 propagator[31][i].real = propagator[31][i].real/3;
	 propagator[31][i].imag = propagator[31][i].imag/3;
       }

    /*======= (a0) S_S p000 =================*/

     n_snk=0; n_src=0;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[32], UNDEFINED);

    for(i=0;i<nt;i++) 
       {
	 propagator[32][i].real *= -1.;
	 propagator[32][i].imag *= -1.;
       }

     /*======= (a1) A1_A1 p000 =================*/
 
     g_snk[0] = GAMMAFIVE ;
     g_snk[1] = XUP;
     n_snk = 2;
     g_src[0] = XUP;
     g_src[1] = GAMMAFIVE;
     n_src = 2;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);

     g_snk[1] = YUP;    g_src[0] = YUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);

     g_snk[1] = ZUP;    g_src[0] = ZUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[33][i].real *= -1.;
	 propagator[33][i].imag *= -1.;
       }



     /*======= (b1) T23_T23 p000 =================*/

     g_snk[0] = YUP ;
     g_snk[1] = ZUP;
     n_snk = 2;
     g_src[0] = ZUP;
     g_src[1] = YUP;
     n_src = 2;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);

     g_snk[0] = XUP ;
     g_snk[1] = ZUP;
     g_src[0] = ZUP;
     g_src[1] = XUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);

     g_snk[0] = XUP ;
     g_snk[1] = YUP;
     g_src[0] = YUP;
     g_src[1] = XUP;
     KS_2pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);


}

/*--------------------------------------------------------------------*/

static char *get_utc_datetime(void)
{
  time_t time_stamp;
  struct tm *gmtime_stamp;
  static char time_string[64];

  time(&time_stamp);
  gmtime_stamp = gmtime(&time_stamp);
  strncpy(time_string,asctime(gmtime_stamp),64);
  
  /* Remove trailing end-of-line character */
  if(time_string[strlen(time_string) - 1] == '\n')
    time_string[strlen(time_string) - 1] = '\0';
  return time_string;
}

/*--------------------------------------------------------------------*/
FILE* 
open_fnal_meson_file(char filename[]){
  FILE *fp;

  /* Create the FNAL file only for rotated meson.  Only node 0
     writes. */
  if(this_node != 0 || saveflag_c == FORGET)
    return NULL;

  fp = fopen(filename,"w");
  if(fp == NULL){
    printf("print_open_fnal_meson_file: ERROR. Can't open %s\n",
	   filename);
    return NULL;
  }
  fprintf(fp,"# meta-data\n");
  fprintf(fp,"header.date      \"%s UTC\"\n",get_utc_datetime());
  fprintf(fp,"header.lattice   %d,%d,%d,%d\n", nx, ny, nz, nt);
  fprintf(fp,"header.job_id    %s\n", job_id);
  return fp;
}

/*--------------------------------------------------------------------*/
static void 
print_start_fnal_meson_group(FILE *fp, char sink[])
{
  if(this_node != 0 || saveflag_c == FORGET)return;
  fprintf(fp,"# %s\n",sink);
}

/*--------------------------------------------------------------------*/
static void 
print_start_fnal_meson_prop(FILE *fp, char kind[], 
			    char sink_label[], 
			    char src_label_w[],
			    char kap_label[], char mass_label[]){
  if(this_node != 0 || saveflag_c == FORGET)return;
  char gamma_kind[16];
  char mom_kind[16];

  /* Parse the "kind" string*/
  sscanf(kind,"%s %s",gamma_kind,mom_kind);

  fprintf(fp,"# %s_%s_%s_k%s_m%s_%s\n",
	  gamma_kind, sink_label, src_label_w, kap_label, 
	  mass_label, mom_kind);
}
		       
/*--------------------------------------------------------------------*/
static void 
print_fnal_meson_prop(FILE *fp, int t, complex c)
{
  if(this_node != 0 || saveflag_c == FORGET)return;
  fprintf(fp, "%d\t%e\t%e\n", t, (double)c.real, (double)c.imag);
}

/*--------------------------------------------------------------------*/
/* Does nothing for now */
static void 
print_end_fnal_meson_prop(FILE *fp){
  if(this_node != 0 || saveflag_c == FORGET)return;
}

/*--------------------------------------------------------------------*/
void 
close_fnal_meson_file(FILE *fp){
  if(this_node != 0 || saveflag_c == FORGET)return;
  if(fp != NULL)fclose(fp);
}

/*--------------------------------------------------------------------*/
static void
clear_mes_prop(complex **prop, int nmes, int ntime){
  int m, t;

  for(m = 0; m < nmes; m++)for(t = 0; t < nt; t++){
      prop[m][t].real = 0.0; prop[m][t].imag = 0.0;
    }
}

/*--------------------------------------------------------------------*/
static complex **
create_mes_prop(int nmes, int ntime){
  complex **prop;
  int m;

  prop = (complex **)malloc(nmes*sizeof(complex *));
  if(prop == NULL)return prop;

  for(m = 0; m < nmes; m++){
    prop[m] = (complex *)malloc(ntime*sizeof(complex));
    if(prop[m] == NULL)return NULL;
  }

  clear_mes_prop(prop, nmes, ntime);

  return prop;
}

/*--------------------------------------------------------------------*/

static void 
destroy_mes_prop(complex **prop, int nmes){
  int m;

  if(prop == NULL)return;

  for(m = 0; m < nmes; m++)
    if(prop[m] != NULL) free(prop[m]);

  free(prop);
}

/*--------------------------------------------------------------------*/
static void
print_hl_rot(FILE *corr_fp, complex **prop_rot, int k)
{
  int i, t;

  char *trace_kind_rot[35] = {
    "P5_P5 p000","P5_P5 p100", "P5_P5 p110","P5_P5 p111", "P5_P5 p200",
    "P5_P5 p210","P5_P5 p211", "P5_P5 p220","P5_P5 p300", "P5_P5 p221",
    "P5_P5 p400","A4_P5 p000", "A4_P5 p100","A4_P5 p110", "A4_P5 p111",
    "A4_P5 p200","A1_P5 p100", "A1_P5 p110","A1_P5 p111", "A1_P5 p200",
    "P5_A1 p100","V1_V1 p000", "V1_V1 p100","V1_V1 p110", "V1_V1 p111",
    "V1_V1 p200","V1_V1 p210", "V1_V1 p211","V1_V1 p220", "V1_V1 p300",
    "V1_V1 p221","V1_V1 p400", "S_S p000"  ,"A1_A1 p000", "T23_T23 p000"
  }; 

  int do_fnal_print[35] = {0,0,0,0,0,
			   0,0,0,0,0,
			   0,1,1,1,1,
			   1,1,1,1,1,
			   0,0,0,0,0,
			   0,0,0,0,0,
			   0,0,0,0,0};

  if(log_correlators)
    node0_printf("BEGIN\n");
      
  print_start_fnal_meson_group(corr_fp,"rotated-sinks");

  for(i=0;i<35;i++)
    {
      if(do_fnal_print[i])
	print_start_fnal_meson_prop(corr_fp, trace_kind_rot[i], 
				    "d",
				    src_label_w[k], kap_label[k], 
				    mass_label);
      
      if(log_correlators)
	node0_printf("\n\nTr %d, %s_k_%f\n_________________________________\n",
		     i, trace_kind_rot[i],kap[k]);
      for(t=0;t<nt;t++)
	{
	  g_complexsum( &prop_rot[i][t] );
	  if(this_node==0){
	    if(do_fnal_print[i])
	      print_fnal_meson_prop(corr_fp, t, prop_rot[i][t]);
	    if(log_correlators)
	      printf("%d %e %e\n", t,
		     prop_rot[i][t].real, prop_rot[i][t].imag);
	  }
	}
      if(do_fnal_print[i])
	print_end_fnal_meson_prop(corr_fp);
    }
  
  clear_mes_prop(prop_rot, 35, nt);

}

/*--------------------------------------------------------------------*/
static void
print_hl_smear(FILE *corr_fp, complex **prop_smear, int k, int ns)
{
  int i, t;
  double space_vol = (double)(nx*ny*nz);

  char *trace_kind_smear[35] = {
    "pi p000",   "pi p100",    "pi p110",   "pi p111",    "pi p200",
    "pi p210",   "pi p211",    "pi p220",   "pi p300",    "pi p221",
    "pi p400",   "A4_P5 p000", "A4_P5 p100","A4_P5 p110", "A4_P5 p111",
    "A4_P5 p200","A1_P5 p100", "A1_P5 p110","A1_P5 p111", "A1_P5 p200",
    "P5_A1 p100","ro_1 p000",  "ro_1 p100", "ro_1 p110",  "ro_1 p111",
    "ro_1 p200", "ro_1 p210",  "ro_1 p211", "ro_1 p220",  "ro_1 p300",
    "ro_1 p221", "ro_1 p400",  "S_S p000",  "A1_A1 p000", "T23_T23 p000"
  }; 

  int do_fnal_print[35] = {1,1,1,1,1,
			   1,1,1,1,1,
			   1,0,0,0,0,
			   0,0,0,0,0,
			   0,1,1,1,1,
			   1,1,1,1,1,
			   1,1,0,0,0};

  print_start_fnal_meson_group(corr_fp,"smeared-sinks");

  for(i=0;i<35;i++){
    if(do_fnal_print[i])
      print_start_fnal_meson_prop(corr_fp, trace_kind_smear[i], 
				  sink_label[ns],
				  src_label_w[k], kap_label[k], mass_label);
      
    if(log_correlators)
      node0_printf("\n\nSMEAR_#%d, %s_k_%f\n_________________________________\n", 
		   ns, trace_kind_smear[i],kap[k]);
    for(t=0;t<nt;t++)
      {
	g_complexsum( &prop_smear[i][t] );
	CDIVREAL(prop_smear[i][t],space_vol,prop_smear[i][t]);
		  
	if(this_node==0){
	  if(do_fnal_print[i])
	    print_fnal_meson_prop(corr_fp, t, prop_smear[i][t]);
	  if(log_correlators)
	    printf("%d %e %e\n", t,
		   prop_smear[i][t].real, 
		   prop_smear[i][t].imag);
	}
      }
    if(do_fnal_print[i])
      print_end_fnal_meson_prop(corr_fp);
  }

  clear_mes_prop(prop_smear, 35, nt);

}

/*--------------------------------------------------------------------*/
void
spectrum_hl_rot(FILE *corr_fp, field_offset snk, field_offset src, int k)
{
  complex **prop_rot;

  prop_rot = create_mes_prop(35, nt);

  All_KS_hl_prop(snk, src, prop_rot);
  print_hl_rot(corr_fp, prop_rot, k);

  destroy_mes_prop(prop_rot, 35);
}

/*--------------------------------------------------------------------*/
void
spectrum_hl_smear(FILE *corr_fp, field_offset snk, field_offset src, 
		  int k, int ns)
{
  complex **prop_smear;

  prop_smear = create_mes_prop(35, nt);
  
  All_KS_hl_prop(snk, src, prop_smear);
  print_hl_smear(corr_fp, prop_smear, k, ns);

  destroy_mes_prop(prop_smear, 35);
}
