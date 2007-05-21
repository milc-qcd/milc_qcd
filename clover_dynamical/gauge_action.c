/************************** gauge_action.c *******************************/
/* MIMD version 7 */
/*  uses staple.e[0][0] to accumulate, and ordinary gathers*/

#include "cl_dyn_includes.h"

void gauge_action(double *result)
{
 double sum_tr;
register int i,k;
int rep;
register site *s;
complex trace;
double g_action,plaq;
double action,act2,total_action;
int dirs[10],sign[10],length;


/* these are for loop_table  */
int ln,iloop;




	plaq=0.0;
	sum_tr=0.0;
	g_action=0.0;



/* gauge action */

   for(iloop=0;iloop<nloop;iloop++){
      length=loop_length[iloop];
   /* loop over rotations and reflections */
      for(ln=0;ln<loop_num[iloop];ln++){
 /* set up dirs and sign  */
       for(k=0;k<length;k++){
                  if(loop_table[iloop][ln][k] < 4 ){sign[k]=1;
                         dirs[k]=(loop_table[iloop][ln][k] )% 4;}
                   else {sign[k]=-1;
                         dirs[k]=(7-loop_table[iloop][ln][k] )% 4;}
                             }

        path(dirs,sign,length);


	FORALLSITES(i,s){
	trace=trace_su3( &s->tempmat1 );
	action=  3.0 - (double)trace.real;
        total_action= (double)loop_coeff[iloop][0]*action;
        act2=action;
        for(rep=1;rep<nreps;rep++){
                act2 *= action;
		total_action += (double)loop_coeff[iloop][rep]*act2;
	}


        g_action  += total_action;



	} /* sites */

        }} /* ln and iloop */

    g_doublesum( &g_action );

/* printf("g_action=%e\n",g_action); */
	*result   =g_action;
/* printf("result %e\n", *result); */


} /* gauge_action */


