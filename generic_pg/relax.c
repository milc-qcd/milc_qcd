/************************** relax.c *******************************/
/* Microcanonical overrelaxation by doing successive SU(2) gauge hits */
/* MIMD version 6 */
/* T. DeGrand March 1991 */
/* UMH: Combined with Schroedinger functional version, Jan 2000 */


#include "generic_pg_includes.h"

#define Nc 3

void relax(int NumStp)
{
  /* Do overrelaxation by SU(2) subgroups */
int NumTrj,Nhit, index1, ina, inb,ii;
int parity;
Real a0,a1,a2,a3,asq,r;
register int dir,i;
register site *st;
su3_matrix action;  su2_matrix u;

Nhit = 3;
	for( NumTrj = 0 ; NumTrj < NumStp; NumTrj++)
	for(parity=ODD;parity<=EVEN;parity++)
	{
	FORALLUPDIR(dir)
		{
              /* compute the gauge force */
		dsdu_qhb(dir,parity);
              /* now for the overrelaxed updating */
		for(index1=0;index1<Nhit;index1++)
		{
                      /*  pick out an SU(2) subgroup */
			ina=(index1+1) % Nc;
			inb=(index1+2) % Nc;
			if(ina > inb){ ii=ina; ina=inb; inb=ii;}


#ifdef SCHROED_FUN
		FORSOMEPARITY(i,st,parity) if(dir==TUP || st->t>0){
#else
		FORSOMEPARITY(i,st,parity){
#endif
			mult_su3_na( &(st->link[dir]), &(st->staple), &action );

/*decompose the action into SU(2) subgroups using Pauli matrix expansion */
/* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
			a0 =  action.e[ina][ina].real + action.e[inb][inb].real;
			a3 =  action.e[ina][ina].imag - action.e[inb][inb].imag;
			a1 =  action.e[ina][inb].imag + action.e[inb][ina].imag;
			a2 =  action.e[ina][inb].real - action.e[inb][ina].real;
      


      /* Normalize and complex conjugate u */

			asq = a0*a0 + a1*a1 + a2*a2 + a3*a3;
			r = sqrt((double)asq );
			a0 = a0/r; a1 = -a1/r; a2 = -a2/r; a3 = -a3/r;
      /* Elements of SU(2) matrix */

			u.e[0][0] = cmplx( a0, a3);
			u.e[0][1] = cmplx( a2, a1);
			u.e[1][0] = cmplx(-a2, a1);
			u.e[1][1] = cmplx( a0,-a3);
    
      /* Do SU(2) hit on all links twice (to overrelax)  */

			left_su2_hit_n(&u,ina,inb,&(st->link[dir]));
			left_su2_hit_n(&u,ina,inb,&(st->link[dir])); 

			} /*   st */
		} /*  hits */
		} /*  direction */
	} /* parity, NumTrj*/
 
} /* relax */

