#include <onium_generic.h>

#define GAMMAFIVE -1
#define eve  0
#define od   1
#define UNDEFINED  2

/* Computes 
   Tr_sc(Gamma_src * A_H(x,0) * Gamma_snk * Q_H(x,0),
   where A_H(x,0) and Q_H could be the heavy quarks,
   at a given space-time point */

complex  two_pt_trace(wilson_propagator * antiquark, wilson_propagator * quark, 
		      int * g_snk, int n_snk, int *g_src, int n_src, int *p, site *s)
{
  int t;
  int my_x;
  int my_y;
  int my_z;
  
  complex trace;
  int s0,s1,s2,s3;
  int c0,c1,i;

  wilson_propagator temp,temp1;
  //su3_matrix mat, mat1;
  

  t = s->t;
  my_x = s->x;
  my_y = s->y;
  my_z = s->z;
  
  temp1 = *quark;

 
  for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++) for(s0=0;s0<4;s0++)for(s1=0;s1<4;s1++)
    {
      temp.c[c0].d[s0].d[s1].c[c1].real = temp1.c[c0].d[s1].d[s0].c[c1].real;
      temp.c[c0].d[s0].d[s1].c[c1].imag = temp1.c[c0].d[s1].d[s0].c[c1].imag;
    }
  

  //multiply by gamma_snk

   for(i=0;i<n_snk;i++)
     for(c0=0;c0<3;c0++){
      mult_swv_by_gamma_l( &(temp.c[c0]), &(temp1.c[c0]), g_snk[i]);
      temp.c[c0] = temp1.c[c0]; 
    } 
   
  //multiply by gamma_src 
   for(c0=0;c0<3;c0++)
     for(i=0;i<n_src;i++)
       {
	 mult_swv_by_gamma_r( &(temp.c[c0]), &(temp1.c[c0]), g_src[i]);   
	 temp.c[c0] = temp1.c[c0];
       }
   trace.real = 0.0;
   trace.imag = 0.0;  
 
   for(c0=0;c0<3;c0++) 
     for(c1=0;c1<3;c1++)       
       for(s0=0;s0<4;s0++)
	 for(s1=0;s1<4;s1++){
	   trace.real +=  temp.c[c0].d[s0].d[s1].c[c1].real * antiquark->c[c0].d[s1].d[s0].c[c1].real+
	                  temp.c[c0].d[s0].d[s1].c[c1].imag * antiquark->c[c0].d[s1].d[s0].c[c1].imag ;
	   trace.imag += -temp.c[c0].d[s0].d[s1].c[c1].real * antiquark->c[c0].d[s1].d[s0].c[c1].imag +
	                  temp.c[c0].d[s0].d[s1].c[c1].imag * antiquark->c[c0].d[s1].d[s0].c[c1].real ;
	 }

            
   return(trace);
}


/* Summation on spatial slice. Caution: elements of prop are not set initially to zero here! 
   In prop averages of propagators are accumulated!  */

void two_pt_func(field_offset snk, field_offset src, int *g_snk, int n_snk,
		 int *g_src, int n_src, int *p, complex *prop, int parity)
{
  int i, t, my_x,my_y,my_z;
  site *s;
  complex trace, trace1, epx,epx1;
  float pi, mom[3];

  wilson_propagator *antiquark;
  wilson_propagator  *quark;


  /* printf("incoming propagator\n");
      for(i=0;i<nt;i++) printf("%.9e %.9e\n", prop[i].real,prop[i].imag);*/

  pi = 4.0 * atan( 1.);
  mom[0] = -2.*pi/(float)nx;  
  mom[1] = -2.*pi/(float)ny;  
  mom[2] = -2.*pi/(float)nz; 
 
  FORALLSITES(i,s){
    t = s->t;
    my_x = s->x;
    my_y = s->y;
    my_z = s->z;
  
    epx1.real = 0.0;
    epx1.imag = (mom[0]*(float)p[0]*(float)my_x + mom[1]*(float)p[1]*(float)my_y 
		+ mom[2]*(float)p[2]*(float)my_z);
    epx = cexp(&epx1);

    antiquark = (wilson_propagator *)F_PT(s, snk);
    quark = (wilson_propagator *)F_PT(s,src);
    trace = two_pt_trace(antiquark, quark, g_snk,n_snk, g_src,n_src, p, s);
    
    
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

}

/* Caution: elements of propagator are not set initially to zero here! */

void All_heavy_prop(field_offset snk, field_offset src, complex **propagator)
{

  int i, j;
  int g_snk[4],g_src[4];
  int n_snk,n_src;
  float temp;

  int p000[3] ={0,0,0}; 

  int p100[3][3] = {{1,0,0},  //0
		    {0,1,0},  //1
		    {0,0,1}}; //2



 int p_100[3][3] = {{-1,0,0},  //0
		    {0,-1,0},  //1
		    {0,0,-1}}; //2

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

  two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p000, propagator[0], UNDEFINED);
  
  
  //*======= (pi) P5_P5 p100 ===================*/

  for(i=0;i<3;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p100[i], propagator[1], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[1][i].real = propagator[1][i].real/3.0;
      propagator[1][i].imag = propagator[1][i].imag/3.0;
    }
  
  /*======= (pi) P5_P5 p110 ===================*/

  for(i=0;i<6;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p110[i], propagator[2], eve);

  for(i=0;i<nt;i++) 
    {
      propagator[2][i].real = propagator[2][i].real/6.;
      propagator[2][i].imag = propagator[2][i].imag/6.;
    }

  /*======= (pi) P5_P5 p111 ===================*/

  for(i=0;i<4;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p111[i], propagator[3], eve);

  for(i=0;i<nt;i++) 
    {
      propagator[3][i].real = propagator[3][i].real/4.;
      propagator[3][i].imag = propagator[3][i].imag/4.;
    }

  /*======= (pi) P5_P5 p200 ===================*/

  for(i=0;i<3;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p200[i], propagator[4], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[4][i].real = propagator[4][i].real/3.;
      propagator[4][i].imag = propagator[4][i].imag/3.;
    }
  
  /*======= (pi) P5_P5 p210 ===================*/

  for(i=0;i<12 ;i++)
    two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[5], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[5][i].real = propagator[5][i].real/12.;
      propagator[5][i].imag = propagator[5][i].imag/12.;
    }
  
  /*======= (pi) P5_P5 p211 ===================*/

  for(i=0;i<12 ;i++)
    two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p211[i], propagator[6], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[6][i].real = propagator[6][i].real/12.;
      propagator[6][i].imag = propagator[6][i].imag/12.;
    }

  /*======= (pi) P5_P5 p220 ===================*/  
 
  for(i=0;i<6;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p220[i], propagator[7], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[7][i].real = propagator[7][i].real/6;
      propagator[7][i].imag = propagator[7][i].imag/6;
    }

  /*======= (pi) P5_P5 p300 ===================*/
  
  for(i=0;i<3;i++)
    two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p300[i], propagator[8], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[8][i].real = propagator[8][i].real/3;
      propagator[8][i].imag = propagator[8][i].imag/3;
    }
  
  /*======= (pi) P5_P5 p221 ===================*/

  for(i=0;i<12 ;i++)
    two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p221[i], propagator[9], eve);
  
  for(i=0;i<nt;i++) 
    {
      propagator[9][i].real = propagator[9][i].real/12;
      propagator[9][i].imag = propagator[9][i].imag/12;
    }

  /*======= (pi) P5_P5 p400 ===================*/

  for(i=0;i<3;i++)
    two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[10], eve);

  
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

   two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p000, propagator[11], UNDEFINED);

   //mulply by i
    for(i=0;i<nt;i++) 
     { 
       temp = propagator[11][i].real;
       propagator[11][i].real = -propagator[11][i].imag;
       propagator[11][i].imag = temp;
     }
  

    /*======= A4_P5 p100 ===================*/
   for(i=0;i<3;i++)
     two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p100[i], propagator[12], eve);
 
   //mulply by i and average over momenta
    for(i=0;i<nt;i++) 
     {
       temp = propagator[12][i].real;
       propagator[12][i].real = -propagator[12][i].imag/3.;
       propagator[12][i].imag = temp/3.;
     }
  
   /*======= A4_P5 p110 ===================*/

   for(i=0;i<6;i++)
     two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p110[i], propagator[13], eve);  
   for(i=0;i<nt;i++) 
     {
       temp = propagator[13][i].real;
       propagator[13][i].real = -propagator[13][i].imag/6;
       propagator[13][i].imag = temp/6;
     }
   
   /*======= A4_P5 p111 ===================*/

   for(i=0;i<4;i++)
     two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p111[i], propagator[14], eve);

   for(i=0;i<nt;i++) 
     {
       temp = propagator[14][i].real;
       propagator[14][i].real = -propagator[14][i].imag/4.;
       propagator[14][i].imag = temp/4.;
     }

   /*======= A4_P5 p200 ===================*/

   for(i=0;i<3;i++)
     two_pt_func(snk, src, g_snk, n_snk, g_src, n_src, p200[i], propagator[15], eve);
   
  
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

   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[0], propagator[16], od);
  
   g_snk[0] = YUP;
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[1], propagator[16], od);
  
   g_snk[0] = ZUP;
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[2], propagator[16], od);
  
   for(i=0;i<nt;i++) 
     {
       propagator[16][i].real /= -3.;
       propagator[16][i].imag /= -3.;
     }
  
   /*======= A1_P5 p110 ===================*/

   g_snk[0] = XUP;
   for(i=0;i<4;i++)
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[i], propagator[17], od);

   g_snk[0] = YUP;
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[0], propagator[17], od); 
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[4], propagator[17], od);
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[5], propagator[17], od);
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p110[6], propagator[17], od);

   g_snk[0] = ZUP;
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[1], propagator[17], od);
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[4], propagator[17], od);
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[7], propagator[17], od);
   two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[8], propagator[17], od);

   for(i=0;i<nt;i++) 
     {
       propagator[17][i].real = -propagator[17][i].real/12;
       propagator[17][i].imag = -propagator[17][i].imag/12;
     }

   /*======= A1_P5 p111 ===================*/
   g_snk[0] = XUP;
   for(i=0;i<4;i++)
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[18], od);

   g_snk[0] = YUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[0], propagator[18], od); 
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[1], propagator[18], od);
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[4], propagator[18], od);
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[5], propagator[18], od);

   g_snk[0] = ZUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[0], propagator[18], od);
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[2], propagator[18], od);
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[4], propagator[18], od);
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[6], propagator[18], od);

     for(i=0;i<nt;i++) 
       {
	 propagator[18][i].real = -propagator[18][i].real/12;
	 propagator[18][i].imag = -propagator[18][i].imag/12;
       }

     /*======= A1_P5 p200 ===================*/
     g_snk[0] = XUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[0], propagator[19], od);

     g_snk[0] = YUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[1], propagator[19], od);

     g_snk[0] = ZUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[2], propagator[19], od);

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

     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[0], propagator[20], od);

     g_snk[0] = YUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[1], propagator[20], od);

     g_snk[0] = ZUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[2], propagator[20], od);

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

     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);

     g_snk[0] = YUP;  g_src[0] = YUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[21], eve);
 
  
     /*======= (ro) V1_V1 p100 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);

     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);
 
     g_snk[0] = ZUP;  g_src[0] = ZUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p100[i], propagator[22], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[22][i].real /= 3.;
	 propagator[22][i].imag /= 3.;
       }

     /*======= (ro) V1_V1 p110 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP;
    for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p110[i], propagator[23], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[23][i].real = propagator[23][i].real/6;
	propagator[23][i].imag = propagator[23][i].imag/6;
      }

    /*======= (ro) V1_V1 p111 ===================*/
     g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<4;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);

   g_snk[0] = YUP;  g_src[0] = YUP; 
   for(i=0;i<4;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);

   g_snk[0] = ZUP;  g_src[0] = ZUP;
   for(i=0;i<4;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p111[i], propagator[24], eve);


    for(i=0;i<nt;i++) 
      {
	propagator[24][i].real = propagator[24][i].real/4;
	propagator[24][i].imag = propagator[24][i].imag/4;
      }

    /*======= (ro) V1_V1 p200 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<3;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<3;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP; 
    for(i=0;i<3;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p200[i], propagator[25], eve);
    
    for(i=0;i<nt;i++) 
      {
	propagator[25][i].real = propagator[25][i].real/3;
	propagator[25][i].imag = propagator[25][i].imag/3;
      }

    /*======= (ro) V1_V1 p210 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP; 
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP;
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP; 
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p210[i], propagator[26], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[26][i].real = propagator[26][i].real/12;
	propagator[26][i].imag = propagator[26][i].imag/12;
      }

    /*======= (ro) V1_V1 p211 ===================*/

    g_snk[0] = XUP; g_src[0] = XUP; 
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p211[i], propagator[27], eve);
    
    g_snk[0] = YUP;  g_src[0] = YUP; 
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p211[i], propagator[27], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<12;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p211[i], propagator[27], eve);

    for(i=0;i<nt;i++) 
      {
	propagator[27][i].real = propagator[27][i].real/12;
	propagator[27][i].imag = propagator[27][i].imag/12;
      }

    /*======= (ro) V1_V1 p220 ===================*/

     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p220[i], propagator[28], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP; 
     for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p220[i], propagator[28], eve);

     g_snk[0] = ZUP;  g_src[0] = ZUP;
     for(i=0;i<6;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p220[i], propagator[28], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[28][i].real = propagator[28][i].real/6;
	 propagator[28][i].imag = propagator[28][i].imag/6;
       }

    /*======= (ro) V1_V1 p300 ===================*/
    
     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);
    
     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p300[i], propagator[29], eve);

     for(i=0;i<nt;i++) 
       {
	 propagator[29][i].real = propagator[29][i].real/3;
	 propagator[29][i].imag = propagator[29][i].imag/3;
       }
     
    /*======= (pi) V1_V1 p221 ===================*/
    
    g_snk[0] = XUP; g_src[0] = XUP;
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p221[i], propagator[30], eve);
   
    g_snk[0] = YUP;  g_src[0] = YUP; 
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src,  p221[i], propagator[30], eve);

    g_snk[0] = ZUP;  g_src[0] = ZUP;
    for(i=0;i<12;i++)
      two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p221[i], propagator[30], eve);
    
    for(i=0;i<nt;i++) 
      {
	propagator[30][i].real = propagator[30][i].real/12;
	propagator[30][i].imag = propagator[30][i].imag/12;
      }
    
   /*======= (ro) V1_V1 p400 ===================*/
    
     g_snk[0] = XUP; g_src[0] = XUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[31], eve);
    
     g_snk[0] = YUP;  g_src[0] = YUP;
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src , p400[i], propagator[31], eve);
    
     g_snk[0] = ZUP;  g_src[0] = ZUP; 
     for(i=0;i<3;i++)
       two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p400[i], propagator[31], eve);
     
     for(i=0;i<nt;i++) 
       {
	 propagator[31][i].real = propagator[31][i].real/3;
	 propagator[31][i].imag = propagator[31][i].imag/3;
       }

    /*======= (a0) S_S p000 =================*/

     n_snk=0; n_src=0;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[32], UNDEFINED);

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
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);

     g_snk[1] = YUP;    g_src[0] = YUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);

     g_snk[1] = ZUP;    g_src[0] = ZUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[33], eve);


     /*======= (b1) T23_T23 p000 =================*/

     g_snk[0] = YUP ;
     g_snk[1] = ZUP;
     n_snk = 2;
     g_src[0] = ZUP;
     g_src[1] = YUP;
     n_src = 2;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);

     g_snk[0] = XUP ;
     g_snk[1] = ZUP;
     g_src[0] = ZUP;
     g_src[1] = XUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);

     g_snk[0] = XUP ;
     g_snk[1] = YUP;
     g_src[0] = YUP;
     g_src[1] = XUP;
     two_pt_func(snk, src,g_snk, n_snk, g_src, n_src, p000, propagator[34], eve);

     //print to file ?


}

int calculate_stag_prop() /* return the C.G. iteration number */
{
  float mass_x2;
  register complex cc,cc2;
  float finalrsq, th;
  register int i,x,y,z,t,icol,cgn;
  site *st;


  /* Fix ZUP Coulomb gauge - gauge links only*/
  rephase( OFF );
  gaugefix(ZUP,(float)1.8,500,(float)GAUGE_FIX_TOL,
	   F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),0,NULL,NULL,0,NULL,NULL);
  rephase( ON );
  
  mass_x2 = 2.*mass;
  cgn=0;
  /* Phase increment for minimum Matsubara frequency */
  th = PI/nt;

  for(icol=0; icol<3; icol++)
    {
      
      /* initialize phi and xxx */
      clear_latvec( F_OFFSET(phi), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(t=0;t<nt;t+=2)
	/**for(x=0;x<1;x+=2)for(y=0;y<1;y+=2)for(t=0;t<1;t+=2)**/
	{
	  if( node_number(x,y,0,t) != mynode() )continue;
	  i=node_index(x,y,0,t);
	  /* Modulate source with Matsubara phase */
	  lattice[i].phi.c[icol].real = -cos((double)t*th);
	}

      /* do a C.G. (source in phi, result in xxx) */
      cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			niter,rsqprop,EVEN,&finalrsq);
      /* Multiply by -Madjoint */
      dslash_ks_( F_OFFSET(xxx), F_OFFSET(ttt), ODD);
      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), EVEN);
      
      /* fill the hadron matrix */
      copy_latvec( F_OFFSET(ttt), F_OFFSET(stag_propagator.e[icol]), EVENANDODD);
    } /* end loop on icol */

  return(cgn);
}
