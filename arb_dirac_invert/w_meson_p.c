/*********** w_meson_p.c *************/
/* MIMD version 6 */
/* UMH April 96 
 2/27/98 Simplify and avoid bug in Origin compiler UMH

TD modify to do many momenta 
 */

#include "arb_dirac_inv_includes.h"

void w_meson_p(field_offset src1,field_offset src2,
Real prop[10][MAX_P][MAX_NT])
/* src1 and src2 are type spin_wilson_vector */
{

int gamma_in[4],gamma_out[4];
int n_in,n_out;
void meson_cont_p(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		Real prop[MAX_P][MAX_NT]);

/* gamma_in = source Dirac matrix, gamma_out = sink Dirac matrix */

/* note (gamma_5 gamma_mu)^\dagger = gamma_5 gamma_mu */

    /* PION */
    n_in=1;n_out=1;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= GAMMAFIVE;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[0]);


    /* PS505 */
    n_in=1;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_out[0]= TUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[1]);


    /* PS055 */
    n_in=2;n_out=1;
    gamma_in[0]= TUP;
    gamma_in[1]= GAMMAFIVE;
    gamma_out[0]= GAMMAFIVE;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[2]);


    /* PS0505 */
    n_in=2;n_out=2;
    gamma_in[0]= TUP;
    gamma_in[1]= GAMMAFIVE;
    gamma_out[0]= TUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[3]);


    /* RHO33 (we use also the other polarizations) */
    n_in=1;n_out=1;
    gamma_in[0]= ZUP;
    gamma_out[0]= ZUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[0]= XUP;
    gamma_out[0]= XUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);

    gamma_in[0]= YUP;
    gamma_out[0]= YUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[4]);


    /* RHO0303 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= TUP;
    gamma_in[1]= ZUP;
    gamma_out[0]= TUP;
    gamma_out[1]= ZUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);

    gamma_in[1]= XUP;
    gamma_out[1]= XUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);

    gamma_in[1]= YUP;
    gamma_out[1]= YUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[5]);


    /* SCALAR */
    n_in=0;n_out=0;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[6]);


    /* SCALA0 */
    n_in=1;n_out=1;
    gamma_in[0]= TUP;
    gamma_out[0]= TUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[7]);


    /* PV35 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= GAMMAFIVE;
    gamma_in[1]= ZUP;
    gamma_out[0]= ZUP;
    gamma_out[1]= GAMMAFIVE;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);

    gamma_in[1]= XUP;
    gamma_out[0]= XUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);

    gamma_in[1]= YUP;
    gamma_out[0]= YUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[8]);


    /* B12 (we use also the other polarizations) */
    n_in=2;n_out=2;
    gamma_in[0]= XUP;
    gamma_in[1]= YUP;
    gamma_out[0]= YUP;
    gamma_out[1]= XUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

    gamma_in[0]= YUP;
    gamma_in[1]= ZUP;
    gamma_out[0]= ZUP;
    gamma_out[1]= YUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

    gamma_in[0]= ZUP;
    gamma_in[1]= XUP;
    gamma_out[0]= XUP;
    gamma_out[1]= ZUP;
    meson_cont_p(src1,src2,gamma_in,gamma_out,n_in,n_out,prop[9]);

} /* w_meson */

/******** meson_cont_p.c *************/
/* MIMD version 4 */
/* UMH April 96 */

void meson_cont_p(field_offset src1,field_offset src2,
		int *gamma_in,int *gamma_out,int n_in,int n_out,
		Real prop[MAX_P][MAX_NT]) 
/* src1 and src2 are type spin_wilson_vector */
{

register int i,t;
register site *s; 

int my_t,j;
int cf, sf, si;
int i_in,i_out;

complex g1,g2;

/* for nonzero momentum */
Real cx,cy,cz,cxy,cyz,cxz,c111;

spin_wilson_vector localmat;       /* temporary storage for quark */
spin_wilson_vector quark;          /* temporary storage for quark */
spin_wilson_vector antiquark;      /* temporary storage for antiquark */

void mult_by_gamma_l( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);
void mult_by_gamma_r( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);



    FORALLSITES(i,s){
	my_t = s->t;

	/*first, dirac multiplication by the source gamma matrices (on left) */

	/*  antiquark = c.c. of quark propagator */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++){
	    CONJG(((spin_wilson_vector *)F_PT(s,src1))->d[si].d[sf].c[cf],
		antiquark.d[si].d[sf].c[cf]);
         }

	/* left multiply antiquark by source gamma matrices,
	   beginning with gamma_5 for quark -> antiquark */
	mult_by_gamma_l( &antiquark, &localmat, GAMMAFIVE);

	/* right dirac multiplication by gamma-5 (finishing up antiquark) */
	mult_by_gamma_r( &localmat, &antiquark, GAMMAFIVE);


	/* left multiply by the particular source dirac matrices */
	for(i_in=0;i_in<n_in;i_in++){
	    mult_by_gamma_l( &antiquark, &localmat, gamma_in[i_in]);
	    antiquark = localmat;
	}

	/* right dirac multiplication by the sink gamma matrices */
	for(i_out=0;i_out<n_out;i_out++) {
	    mult_by_gamma_r( &antiquark, &localmat, gamma_out[i_out]);
	    antiquark = localmat;
	}

	/* copy into quark */
	/* 2/27/98 Simplify and avoid bug in Origin compiler UMH */
        quark = *(spin_wilson_vector *)F_PT(s,src2);

	/* trace over propagators */
	for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	for(cf=0;cf<3;cf++)
	{
	    g1 = antiquark.d[si].d[sf].c[cf];
	    g2 = quark.d[si].d[sf].c[cf];

		for(j=0;j<3;j++){
		cz=cos(2.0*PI/(Real)nz*(Real)(s->z)*(Real)j);
		cx=cos(2.0*PI/(Real)nx*(Real)(s->x)*(Real)j);
		cy=cos(2.0*PI/(Real)ny*(Real)(s->y)*(Real)j);

 
		prop[j][my_t] += (g1.real*g2.real - g1.imag*g2.imag)*
                                (cx+cy+cz)/3.0;
		}
		cxy=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y));
		cxz=cos(2.0*PI/(Real)nz*(Real)(s->x +s->z));
		cyz=cos(2.0*PI/(Real)nz*(Real)(s->y +s->z));
		c111=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y + s->z));
		prop[3][my_t] += (g1.real*g2.real - g1.imag*g2.imag)*
                                (cxy+cyz+cxz)/3.0;
		prop[4][my_t] += (g1.real*g2.real  -g1.imag*g2.imag)*
                                c111;

	}

    }

} /* meson_cont_p */


/************* mb_gamma_l.c  (in su3.a) **************************/
/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the row index
  (This is the first index, or equivalently, multiplication on the left)
  usage:  mb_gamma_l( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE

 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/

void mult_by_gamma_l( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir)
{
register int c2,s2;	/* column indices, color and spin */

  switch(dir){
    case XUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSI(  src->d[3].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[2].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSI( src->d[1].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case YUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESMINUSONE( src->d[3].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[2].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[0].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case ZUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSI(  src->d[2].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESMINUSI( src->d[3].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSI( src->d[0].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESPLUSI(  src->d[1].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case TUP:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSONE( src->d[2].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[3].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[0].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESPLUSONE( src->d[1].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    case GAMMAFIVE:
	for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	    TIMESPLUSONE(  src->d[0].d[s2].c[c2],
		dest->d[0].d[s2].c[c2] );
	    TIMESPLUSONE(  src->d[1].d[s2].c[c2],
		dest->d[1].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[2].d[s2].c[c2],
		dest->d[2].d[s2].c[c2] );
	    TIMESMINUSONE( src->d[3].d[s2].c[c2],
		dest->d[3].d[s2].c[c2] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_LEFT()\n");
  }
}

/************* mb_gamma_r.c  (in su3.a) **************************/
/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/
void mult_by_gamma_r( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir)
{
register int i; /*color*/
register int s1;	/* row  spin indices*/

  switch(dir){
    case XUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSI( src->d[s1].d[3].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESMINUSI( src->d[s1].d[2].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].d[1].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESPLUSI(  src->d[s1].d[0].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case YUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[2].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[1].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].d[0].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case ZUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESMINUSI( src->d[s1].d[2].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSI(  src->d[s1].d[3].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSI(  src->d[s1].d[0].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSI( src->d[s1].d[1].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case TUP:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESPLUSONE( src->d[s1].d[2].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESPLUSONE( src->d[s1].d[0].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESPLUSONE( src->d[s1].d[1].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    case GAMMAFIVE:
	for(i=0;i<3;i++)for(s1=0;s1<4;s1++){
	    TIMESPLUSONE(  src->d[s1].d[0].c[i],
		dest->d[s1].d[0].c[i] );
	    TIMESPLUSONE(  src->d[s1].d[1].c[i],
		dest->d[s1].d[1].c[i] );
	    TIMESMINUSONE( src->d[s1].d[2].c[i],
		dest->d[s1].d[2].c[i] );
	    TIMESMINUSONE( src->d[s1].d[3].c[i],
		dest->d[s1].d[3].c[i] );
	}
	break;
    default:
	printf("BAD CALL TO MULT_BY_GAMMA_RIGHT()\n");
  }
}

