#include "generic_wilson_includes.h"

//V G V+
void weyl2canopy_w_rot(field_offset src, field_offset dest)
{
  int s0,s1,c0,c1,i;
  site *s;
  complex z1,z2;
  wilson_propagator q, qtmp;

  FORALLSITES(i,s){
   q = *(wilson_propagator *)F_PT(s,src);

    for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
      {
	for(s0=0;s0<4;s0++){
	  qtmp.c[c0].d[0].d[s0].c[c1].real =
	     q.c[c0].d[1].d[s0].c[c1].real+
	     q.c[c0].d[3].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[0].d[s0].c[c1].imag =
	     q.c[c0].d[1].d[s0].c[c1].imag-
	     q.c[c0].d[3].d[s0].c[c1].real;
 
	  qtmp.c[c0].d[1].d[s0].c[c1].real =
	    -q.c[c0].d[0].d[s0].c[c1].real-
	     q.c[c0].d[2].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[1].d[s0].c[c1].imag =
	    -q.c[c0].d[0].d[s0].c[c1].imag+
	     q.c[c0].d[2].d[s0].c[c1].real;

 	  qtmp.c[c0].d[2].d[s0].c[c1].real =
	     q.c[c0].d[1].d[s0].c[c1].real-
	     q.c[c0].d[3].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[2].d[s0].c[c1].imag =
	     q.c[c0].d[1].d[s0].c[c1].imag+
	     q.c[c0].d[3].d[s0].c[c1].real;
 

	  qtmp.c[c0].d[3].d[s0].c[c1].real =
	    -q.c[c0].d[0].d[s0].c[c1].real+
	     q.c[c0].d[2].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[3].d[s0].c[c1].imag =
	    -q.c[c0].d[0].d[s0].c[c1].imag-
	     q.c[c0].d[2].d[s0].c[c1].real;
	}

	for(s0=0;s0<4;s0++){
	     q.c[c0].d[s0].d[0].c[c1].real =
     0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real-
	  qtmp.c[c0].d[s0].d[3].c[c1].imag);
	  
	     q.c[c0].d[s0].d[0].c[c1].imag =
     0.5*(qtmp.c[c0].d[s0].d[1].c[c1].imag+
	  qtmp.c[c0].d[s0].d[3].c[c1].real);
 
	     q.c[c0].d[s0].d[1].c[c1].real =
    0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	  qtmp.c[c0].d[s0].d[2].c[c1].imag);
	  
	     q.c[c0].d[s0].d[1].c[c1].imag =
    0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].imag-
	  qtmp.c[c0].d[s0].d[2].c[c1].real);

	     q.c[c0].d[s0].d[2].c[c1].real =
     0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real+
	  qtmp.c[c0].d[s0].d[3].c[c1].imag);
	  

	     q.c[c0].d[s0].d[2].c[c1].imag =
     0.5*(qtmp.c[c0].d[s0].d[1].c[c1].imag-
	  qtmp.c[c0].d[s0].d[3].c[c1].real);
 

	    q .c[c0].d[s0].d[3].c[c1].real =
    0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real-
	  qtmp.c[c0].d[s0].d[2].c[c1].imag);
	  
	     q.c[c0].d[s0].d[3].c[c1].imag =
    0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].imag+
	  qtmp.c[c0].d[s0].d[2].c[c1].real);
	}
      }
    *(wilson_propagator *)F_PT(s,dest)=q;
  }

}

//V+ G V

void bj_to_w_rot(field_offset src, field_offset dest)
{
  int s0,s1,c0,c1,i;
  site *s;
  complex z1,z2;
  wilson_propagator q, qtmp;

  FORALLSITES(i,s){
     q = *(wilson_propagator *)F_PT(s,src);

    for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
      {
	for(s0=0;s0<4;s0++){
	  qtmp.c[c0].d[0].d[s0].c[c1].real =
	    (q.c[c0].d[1].d[s0].c[c1].imag-
	     q.c[c0].d[3].d[s0].c[c1].imag);
	  
	  qtmp.c[c0].d[0].d[s0].c[c1].imag =
	   (-q.c[c0].d[1].d[s0].c[c1].real+
	     q.c[c0].d[3].d[s0].c[c1].real);
 
	  qtmp.c[c0].d[1].d[s0].c[c1].real =
	    -q.c[c0].d[0].d[s0].c[c1].imag+
	     q.c[c0].d[2].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[1].d[s0].c[c1].imag =
             q.c[c0].d[0].d[s0].c[c1].real-
	     q.c[c0].d[2].d[s0].c[c1].real;

 	  qtmp.c[c0].d[2].d[s0].c[c1].real =
	     q.c[c0].d[1].d[s0].c[c1].imag+
	     q.c[c0].d[3].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[2].d[s0].c[c1].imag =
	    -q.c[c0].d[1].d[s0].c[c1].real-
	     q.c[c0].d[3].d[s0].c[c1].real;
 

	  qtmp.c[c0].d[3].d[s0].c[c1].real =
	    -q.c[c0].d[0].d[s0].c[c1].imag-
	     q.c[c0].d[2].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[3].d[s0].c[c1].imag =
	     q.c[c0].d[0].d[s0].c[c1].real+
	     q.c[c0].d[2].d[s0].c[c1].real;
	}

	for(s0=0;s0<4;s0++){
	      q.c[c0].d[s0].d[0].c[c1].real =
     0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag+
	   qtmp.c[c0].d[s0].d[3].c[c1].imag);
	  
	      q.c[c0].d[s0].d[0].c[c1].imag =
     0.5*( qtmp.c[c0].d[s0].d[1].c[c1].real-
	   qtmp.c[c0].d[s0].d[3].c[c1].real);
 
	      q.c[c0].d[s0].d[1].c[c1].real =
      0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag-
	   qtmp.c[c0].d[s0].d[2].c[c1].imag);
	  
	      q.c[c0].d[s0].d[1].c[c1].imag =
     0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	   qtmp.c[c0].d[s0].d[2].c[c1].real);

 	      q.c[c0].d[s0].d[2].c[c1].real =
     0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag-
	   qtmp.c[c0].d[s0].d[3].c[c1].imag);
	  

	      q.c[c0].d[s0].d[2].c[c1].imag =
      0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real+
	   qtmp.c[c0].d[s0].d[3].c[c1].real);
 

	      q.c[c0].d[s0].d[3].c[c1].real =
      0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag+
	   qtmp.c[c0].d[s0].d[2].c[c1].imag);
	  
	      q.c[c0].d[s0].d[3].c[c1].imag =
     0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real-
	   qtmp.c[c0].d[s0].d[2].c[c1].real);
	}
      }
    *(wilson_propagator *)F_PT(s,dest)=q;
  }

}


void canopy2weyl_w_rot(field_offset src, field_offset dest)
{
  int s0,s1,c0,c1,i;
  site *s;
  complex z1,z2;
  wilson_propagator q, qtmp;

  FORALLSITES(i,s){
     q = *(wilson_propagator *)F_PT(s,src);

    for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
      {
	for(s0=0;s0<4;s0++){
	  qtmp.c[c0].d[0].d[s0].c[c1].real =
	-1.*(q.c[c0].d[1].d[s0].c[c1].real+
	     q.c[c0].d[3].d[s0].c[c1].real);
	  
	  qtmp.c[c0].d[0].d[s0].c[c1].imag =
	-1.*(q.c[c0].d[1].d[s0].c[c1].imag+
	     q.c[c0].d[3].d[s0].c[c1].imag);
 
	  qtmp.c[c0].d[1].d[s0].c[c1].real =
	     q.c[c0].d[0].d[s0].c[c1].real+
	     q.c[c0].d[2].d[s0].c[c1].real;
	  
	  qtmp.c[c0].d[1].d[s0].c[c1].imag =
             q.c[c0].d[0].d[s0].c[c1].imag+
	     q.c[c0].d[2].d[s0].c[c1].imag;

 	  qtmp.c[c0].d[2].d[s0].c[c1].real =
	     q.c[c0].d[1].d[s0].c[c1].imag-
	     q.c[c0].d[3].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[2].d[s0].c[c1].imag =
	    -q.c[c0].d[1].d[s0].c[c1].real+
	     q.c[c0].d[3].d[s0].c[c1].real;
 

	  qtmp.c[c0].d[3].d[s0].c[c1].real =
	    -q.c[c0].d[0].d[s0].c[c1].imag+
	     q.c[c0].d[2].d[s0].c[c1].imag;
	  
	  qtmp.c[c0].d[3].d[s0].c[c1].imag =
	     q.c[c0].d[0].d[s0].c[c1].real-
	     q.c[c0].d[2].d[s0].c[c1].real;
	}

	for(s0=0;s0<4;s0++){
	      q.c[c0].d[s0].d[0].c[c1].real =
     0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].real-
	   qtmp.c[c0].d[s0].d[3].c[c1].real);
	  
	      q.c[c0].d[s0].d[0].c[c1].imag =
     0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag-
	   qtmp.c[c0].d[s0].d[3].c[c1].imag);
 
	      q.c[c0].d[s0].d[1].c[c1].real =
      0.5*(qtmp.c[c0].d[s0].d[0].c[c1].real+
	   qtmp.c[c0].d[s0].d[2].c[c1].real);
	  
	      q.c[c0].d[s0].d[1].c[c1].imag =
      0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag+
	   qtmp.c[c0].d[s0].d[2].c[c1].imag);

 	      q.c[c0].d[s0].d[2].c[c1].real =
     0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag+
	   qtmp.c[c0].d[s0].d[3].c[c1].imag);
	  

	      q.c[c0].d[s0].d[2].c[c1].imag =
      0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real-
	   qtmp.c[c0].d[s0].d[3].c[c1].real);
 

	      q.c[c0].d[s0].d[3].c[c1].real =
      0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag-
	   qtmp.c[c0].d[s0].d[2].c[c1].imag);
	  
	      q.c[c0].d[s0].d[3].c[c1].imag =
     0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	   qtmp.c[c0].d[s0].d[2].c[c1].real);
	}
      }
    *(wilson_propagator *)F_PT(s,dest)=q;
  }

}
